#[macro_use]
extern crate log;
extern crate bam;
use bam::Record;
use bam::RecordReader;
use bam::{
    header::HeaderEntry,
    record::tags::{StringType, TagValue},
    Header, RecordWriter,
};
use bio::alphabets::dna::revcomp;
use byteorder::WriteBytesExt;
use io::BufReader;
use itertools::Itertools;
use regex::Regex;
use std::env;
use std::io::Write;
use std::{
    fs::File,
    io, iter,
    process::{Command, Stdio},
};

struct RecordIter<'a, I: Iterator<Item = &'a Record>>(I);

impl<'a, I> RecordReader for RecordIter<'a, I>
where
    I: Iterator<Item = &'a Record>,
{
    fn read_into(&mut self, record: &mut Record) -> std::io::Result<bool> {
        if let Some(next_record) = self.0.next() {
            // record = next_record.1;
            std::mem::replace(record, next_record.clone());
            Ok(true)
        } else {
            Ok(false)
        }
    }

    fn pause(&mut self) {}
}

/// Iterator over records.
impl<'a, I> Iterator for RecordIter<'a, I>
where
    I: Iterator<Item = &'a Record>,
{
    type Item = std::io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Record::new();
        match self.read_into(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub(crate) fn write_iterator<W, I>(writer: &mut W, mut iterator: I) -> io::Result<()>
where
    W: Write,
    I: Iterator<Item = u8>,
{
    const SIZE: usize = 1024;
    let mut buffer = [0_u8; SIZE];
    loop {
        for i in 0..SIZE {
            match iterator.next() {
                Some(value) => buffer[i] = value,
                None => {
                    return writer.write_all(&buffer[..i]);
                }
            }
        }
        writer.write_all(&buffer)?;
    }
}

fn select_base<'a>(
    read: &'a Vec<char>,
    seqs: Vec<(usize, &'a Vec<char>)>,
    alpha: f64, //Here, if 1.5 sigma, mean - 1.5 * stdev
                //    sigma: f64, // Here, if 1.5 sigma, 1.5 * stdev
) -> &'a Vec<char> {
    let d: usize = seqs.iter().map(|t| t.0).sum();
    let threshold = d as f64 * (1f64 - alpha);
    let minor = d - seqs[0].0;
    //eprintln!("{}, {}, {}, {:?},  {:?}", threshold, minor, d, read, seqs);
    eprintln!(
        "Agg:\t{}\t{}\t{}\t{}\t{}",
        d,
        seqs[0].0,
        minor,
        seqs[0].0 as f64 / d as f64,
        minor as f64 / d as f64
    );
    if minor as f64 <= threshold {
        return seqs[0].1;
    } else {
        return read;
    }
}

fn calculate_sam<'a, F>(primary: Vec<(Record, Option<&str>, Vec<u8>)>, lambda: F) -> Vec<Record>
where
    F: Fn(u32) -> Option<&'a str>,
{
    let mut primary = primary.clone();
    if primary.iter().all(|t| t.0.ref_id() == 1) {
        return primary.into_iter().map(|t| t.0).collect::<Vec<_>>();
    }
    primary.sort_by(|a, b| {
        (a.0.cigar().soft_clipping(!a.0.flag().is_reverse_strand())
            + a.0.cigar().hard_clipping(!a.0.flag().is_reverse_strand()))
        .cmp(
            &((b.0.cigar().soft_clipping(!b.0.flag().is_reverse_strand()))
                + (b.0.cigar().hard_clipping(!b.0.flag().is_reverse_strand()))),
        )
    });
    let cigars = primary
        .iter()
        .map(|t| {
            let mut readable: Vec<u8> = Vec::new();
            t.0.cigar().write_readable(&mut readable);

            let strand = if t.0.flag().is_reverse_strand() {
                "-"
            } else {
                "+"
            };
            let nm = match t.0.tags().get(b"NM") {
                Some(TagValue::String(array_view, StringType::String)) => array_view,
                _ => &[],
            };

            [
                lambda(t.0.ref_id() as u32).unwrap(),
                &t.0.start().to_string(),
                strand,
                &String::from_utf8_lossy(&readable),
                &t.0.mapq().to_string(),
                &String::from_utf8_lossy(&nm),
            ]
            .join(",")
        })
        .collect::<Vec<_>>();
    primary
        .into_iter()
        .enumerate()
        .map(|(i, mut t)| {
            t.0.tags_mut().push_string(
                b"SA",
                cigars
                    .iter()
                    .enumerate()
                    .filter(|t| t.0 != i)
                    .map(|t| t.1)
                    .join(";")
                    .as_bytes(),
            );
            t.0
        })
        .collect::<Vec<_>>()
}

fn calculate_primary<'a>(
    primary: Vec<(Record, Option<&str>, Vec<u8>)>,
    name_vec: Vec<u8>, //fasta_reader: bio::io::fasta::IndexedReader<File>,
    realigner: &str,
    alpha: f64,
) {
    let name = String::from_utf8_lossy(&name_vec);
    let mut process = match Command::new(realigner)
        .args(&[name.to_string()])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
    {
        Err(why) => panic!("couldn't spawn {}: {}", realigner, why),
        Ok(process) => process,
    };
    {
        // let output = io::BufWriter::new(io::stdout());
        let output = io::BufWriter::new(process.stdin.take().unwrap());
        let mut header = Header::new();

        header.push_entry(HeaderEntry::ref_sequence(
            name.to_string(),
            primary[0].0.query_len(),
        ));

        let mut writer = bam::SamWriter::build().from_stream(output, header).unwrap();

        // let primary = primary.iter().find(|t| t.seq());
        let mut read = vec![];
        let len = primary.len();

        for (record, ref_id, ref_seq) in primary {
            let record = &mut record.clone();
            let mut readable: Vec<u8> = Vec::new();
            let cigar = record.cigar();
            if record.flag().is_reverse_strand() {
                for (len, op) in cigar.iter().rev() {
                    write!(readable, "{}", len).unwrap();
                    readable.write_u8(op.to_byte()).unwrap();
                }
            } else {
                record.cigar().write_readable(&mut readable);
            }
            let readable_string = String::from_utf8_lossy(&readable);
            let readable_string2 = readable_string.replace("D", "Z");
            let readable_string3 = readable_string2.replace("I", "D");
            let readable_string4 = readable_string3.replace("Z", "I");
            let clipping = Regex::new(r"\d*[H|S]").unwrap();
            let clip_removed = clipping.replace_all(&readable_string4, "");

            //record.cigar().clear();
            let bytes = clip_removed.bytes().into_iter();
            //record.cigar().extend_from_text(bytes);
            let record_str = format!(
                "{}:{}-{}",
                ref_id.unwrap(),
                record.start(),
                record.calculate_end()
            );
            // record.set_start(record.cigar().soft_clipping(true) as i32);
            record.set_start(record.cigar().soft_clipping(true) as i32);
            record.set_cigar(bytes);

            record.tags_mut().remove(b"MD");

            record.set_ref_id(0);
            record.set_name(record_str.bytes());

            let ref_seq = if record.flag().is_reverse_strand() {
                revcomp(ref_seq)
            } else {
                ref_seq
            };
            // eprintln!("{} {}", ref_seq.len(), record.calculate_end() - record.start());

            // record.sequence().clear();
            // record.sequence().extend_from_text(ref_seq);
            //record.reset_seq();
            //record.set_seq(ref_seq);
            if record.sequence().to_vec().len() > 0 {
                if record.flag().is_reverse_strand() {
                    read = record.sequence().rev_compl(..).collect::<_>();
                } else {
                    read = record.sequence().to_vec();
                }
            }
            record.set_seq_qual(
                ref_seq,
                iter::empty(), // record.qualities().to_readable().into_iter().map(|q| q - 33),
            );

            writer.write(&record).unwrap();
        }

        let mut record = Record::new();
        record.set_name(name.bytes());
        record.set_ref_id(0);
        record.set_start(0);
        record.set_cigar(format!("{}M", read.len()).bytes());
        record.set_seq_qual(
            read,
            iter::empty(), // record.qualities().to_readable().into_iter().map(|q| q - 33),
        );

        writer.write(&record).unwrap();

        writer.flush().unwrap();
        writer.finish().unwrap();
        std::mem::drop(writer);
        println!(">{} {}", name.to_string(), len);
    }
    // std::mem::drop(process.stdin);
    //let output = process
    //.wait_with_output()
    //.expect("failed to wait on child");

    let stdout = BufReader::new(process.stdout.take().unwrap());
    let reader = bam::SamReader::from_stream(stdout).unwrap();

    let mut records = vec![];
    for i in reader {
        records.push(i.unwrap());
    }
    records.sort_by_key(|i| i.start());

    for column in bam::Pileup::with_filter(&mut RecordIter(records.iter()), |_record| true) {
        let column = column.unwrap();
        /*eprintln!("Column at {}:{}, {} records", column.ref_id(),
        column.ref_pos() + 1, column.entries().len());*/
        let mut seqs = Vec::with_capacity(column.entries().len()); //vec![];
        let mut read = vec![];
        for entry in column.entries().iter() {
            if entry.record().name().to_vec() == name_vec {
                read = entry
                    .sequence()
                    .unwrap()
                    .map(|nt| nt as char)
                    .collect::<Vec<_>>();
            }
            let seq: Vec<_> = entry.sequence().unwrap().map(|nt| nt as char).collect();
            //let qual: Vec<_> = entry.qualities().unwrap().iter()
            //    .map(|q| (q + 33) as char).collect();
            //eprintln!("    {:?}: {:?}", entry.record(), seq);
            seqs.push(seq);
        }
        let unique_elements = seqs.iter().cloned().unique().collect_vec();
        let mut unique_frequency = vec![];
        for unique_elem in unique_elements.iter() {
            unique_frequency.push((
                seqs.iter().filter(|&elem| elem == unique_elem).count(),
                unique_elem,
            ));
        }
        unique_frequency.sort_by_key(|t| t.0);
        unique_frequency.reverse();
        if unique_frequency.len() > 0 {
            print!(
                "{}",
                select_base(&read, unique_frequency, alpha)
                    .iter()
                    .collect::<String>()
            );
        }
    }
    println!("");

    //stdout.get_mut().wait_with_output().unwrap();
    process.wait_with_output().unwrap();

    // process.kill();
}

fn bam_stats(path: String) {
    let bam_stream = BufReader::with_capacity(1000000, File::open(path).unwrap());
    let reader = bam::BamReader::from_stream(bam_stream, 2).unwrap();
    let mut total_read_length: u64 = 0;
    let mut total_primary_aligned_read_length: u64 = 0;
    let mut unaligned_length: u64 = 0;
    let mut primary_alignment = 0;
    let mut unaligned_reads = 0;

    for record in reader {
        let record = record.unwrap();
        if !record.flag().is_secondary() && !record.flag().is_supplementary() && record.flag().is_mapped() {
            //let read_original_length = record.
            total_read_length += record.query_len() as u64 + record.cigar().hard_clipping(true) as u64 + record.cigar().hard_clipping(false) as u64;
            total_primary_aligned_read_length += (record.aligned_query_end() - record.aligned_query_start()) as u64;
            primary_alignment += 1;
        } else if !record.flag().is_mapped()
        {
            // total_read_length += record.query_len();
            unaligned_length += record.query_len() as u64;
            unaligned_reads += 1;
        }
    }
    println!("Total read length\t{}\nAligned length\t{}\nUnaligned_length\t{}\nPrimary alignment ratio\t{}\nPrimary alignment ratio with unaligned\t{}",
     total_read_length, total_primary_aligned_read_length, unaligned_length, total_primary_aligned_read_length as f64 / total_read_length as f64, total_primary_aligned_read_length as f64 / (total_read_length + unaligned_length) as f64);
    println!("# of primary alignment\t{}\n# of unaligned reads\t{}", primary_alignment, unaligned_reads);
    return
}

fn main() {
    let args: Vec<String> = env::args().collect();
    // We assume that input is sorted by "samtools sort -n".
    //let reader = bam::IndexedReader::from_path(args[1].clone()).unwrap();
    let command = &args[1];
    let sa_merge = command == "attachsa";
    let realigner = command == "realign";
    let stats = command == "stats" || command == "stat";
    if !sa_merge && !realigner && !stats {
        let string = "    bam-handler
        
    USAGE:
        bam-handler <SUBCOMMAND>

    SUBCOMMANDS:
        attatchsa Attach SA-tag in bam file from an output of LAST-split.
        realign   Realign aligned reads in bam file.
        stats     Collects statistics from a BAM file.
        ";
        println!("{}", string);
        return;
    }
    if stats {
        bam_stats(args[2].clone());
        return;
    }
    let bam_stream = BufReader::with_capacity(1000000, File::open(args[2].clone()).unwrap());
    //let reader = bam::BamReader::from_path(args[1].clone()).unwrap();
    let reader = bam::BamReader::from_stream(bam_stream, 2).unwrap();
    let x = args
        .get(5)
        .and_then(|a| a.parse::<usize>().ok())
        .unwrap_or(5usize);
    let alpha = args
        .get(6)
        .and_then(|a| a.parse::<f64>().ok())
        .unwrap_or(0.85);
    let sigma = args
        .get(7)
        .and_then(|a| a.parse::<f64>().ok())
        .unwrap_or(0.075);
    let dummy = "dummy.fa".to_string();
    let path = args.get(3).unwrap_or(&dummy);
    let realigner = args.get(4).unwrap_or(&dummy);

    //let sa_merge = args.len() <= 3;
    let mut writer = bam::BamWriter::build()
        .write_header(true)
        .from_stream(std::io::stdout(), reader.header().clone())
        .unwrap();

    let mut previous_name = vec![];
    let mut primary = vec![];
    let header = reader.header().clone();
    let closure = |x: u32| header.reference_name(x);

    for record in reader {
        let record = record.unwrap();
        // let closure = |x: u32| reader.header().reference_name(x);
        let ref_name = closure(record.ref_id() as u32);
        // eprintln!("{:?} {}", ref_name, record.ref_id());
        //let ref_name = &reader.header().reference_name(record.ref_id() as u32);
        let mut ref_seq = vec![];
        if !sa_merge{
        if let Some(ref_name) = ref_name {
            let mut fasta_reader = bio::io::fasta::IndexedReader::from_file(path).unwrap();
            fasta_reader.fetch(
                ref_name,
                record.start() as u64,
                record.calculate_end() as u64,
            );
            fasta_reader.read(&mut ref_seq);
        }
    }
        // eprintln!("{:?} {}", record, record.flag().is_supplementary());
        if previous_name == record.name() {
            if !record.flag().is_supplementary() {
                primary.push((record, ref_name, ref_seq));
            }
        } else {
            if primary.len() > x {
                if sa_merge {
                    for i in calculate_sam(primary, closure) {
                        writer.write(&i).unwrap();
                    }
                } else {
                    calculate_primary(primary, previous_name, &args[3], alpha - sigma);
                }
            } else if primary.len() > 0 {
                if sa_merge {
                    for i in calculate_sam(primary, closure) {
                        writer.write(&i).unwrap();
                    }
                } else {
                    println!(">{}", String::from_utf8_lossy(&previous_name));
                    for (record, _, _) in primary {
                        if record.sequence().len() > 0 {
                            let seq = record.sequence();
                            if record.flag().is_reverse_strand() {
                                write_iterator(
                                    &mut io::stdout(),
                                    revcomp((0..seq.len()).map(|i| seq.at(i))).into_iter(),
                                )
                                .unwrap();
                            } else {
                                write_iterator(
                                    &mut io::stdout(),
                                    (0..seq.len()).map(|i| seq.at(i)),
                                )
                                .unwrap();
                            }
                        }
                    }
                    println!("");
                }
            }
            let previous = record.clone();
            previous_name = previous.name().to_vec().clone();
            primary = vec![];
            if !record.flag().is_supplementary() {
                //read_tree.insert(record);
                primary.push((record, ref_name, ref_seq));
            }
        }
    }
    if sa_merge {
        for i in calculate_sam(primary, closure) {
            writer.write(&i).unwrap();
        }
    } else {
        if primary.len() > x {
            calculate_primary(primary, previous_name, &realigner, alpha - sigma);
        } else {
            println!(">{}", String::from_utf8_lossy(&previous_name));
            for (record, _, _) in primary {
                if record.sequence().len() > 0 {
                    //record.sequence().write_readable(&mut io::stdout());
                    let seq = record.sequence();
                    if record.flag().is_reverse_strand() {
                        write_iterator(
                            &mut io::stdout(),
                            revcomp((0..seq.len()).map(|i| seq.at(i))).into_iter(),
                        )
                        .unwrap();
                    } else {
                        write_iterator(&mut io::stdout(), (0..seq.len()).map(|i| seq.at(i)))
                            .unwrap();
                    }
                }
            }
            println!("");
        }
    }
}
