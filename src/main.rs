extern crate bam;
use byteorder::WriteBytesExt;
use std::io::Write;
use bam::Record;
use bam::{header::HeaderEntry, Header, RecordWriter};
use bio::alphabets::dna::revcomp;
use io::BufReader;
use itertools::Itertools;
use regex::Regex;
use std::env;
use bam::RecordReader;
use std::{
    collections::BTreeMap,
    fs::File,
    io,
    process::{Command, Stdio}, iter,
};

trait Test {
    fn new() -> Self;
}

struct A {
    data: Vec<u32>,
}

struct B {
    data: Vec<i32>,
}

impl Test for A {
    fn new() -> A {
        return A { data: vec![] };
    }
}

impl Test for B {
    fn new() -> B {
        return B { data: vec![] };
    }
}

impl Test3 for B {
    fn new(&mut self) -> Option<bool> {
        Some(true)
    }
}

trait Test2 {
    fn new(&mut self) -> Option<bool>;
    fn to_stream<W: io::Write>(&self, stream: &mut W) -> io::Result<()>;
    fn from_stream<R: io::Read>(&mut self, stream: &mut R) -> io::Result<bool>;
}

trait Test3 {
    fn new(&mut self) -> Option<bool>;
}

struct D {
    data: Test3,
}
/*
impl D {
    pub fn new() -> D {
        D{data: Test3{}}
    }
}
*/
struct C<T: Test> {
    data: T,
}

impl<T: Test> C<T> {
    pub fn new() -> C<T> {
        C { data: T::new() }
    }
}

fn calculate_primary<'a>(
    primary: Vec<(Record, Option<&str>, Vec<u8>)>,
    name_vec: Vec<u8>, //fasta_reader: bio::io::fasta::IndexedReader<File>,
    realigner: &str,
) {
    let name = String::from_utf8_lossy(&name_vec);
    let process = match Command::new(realigner)
        .args(&[name.to_string()])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
    {
        Err(why) => panic!("couldn't spawn {}: {}", realigner, why),
        Ok(process) => process,
    };

    // let output = io::BufWriter::new(io::stdout());
    let output = io::BufWriter::new(process.stdin.unwrap());
    let mut header = Header::new();

    header.push_entry(HeaderEntry::ref_sequence(
        name.to_string(),
        primary[0].0.query_len(),
    ));

    let mut writer = bam::SamWriter::build().from_stream(output, header).unwrap();

    // let primary = primary.iter().find(|t| t.seq());
    let mut read = vec![];


    for (record, ref_id, ref_seq) in primary {
        let record = &mut record.clone();
        let mut readable: Vec<u8> = Vec::new();
        let mut cigar = record.cigar();
        if record.flag().is_reverse_strand() {
        for (len, op) in cigar.iter().rev() {
            write!(readable, "{}", len).unwrap();
            readable.write_u8(op.to_byte()).unwrap();
        } } else {
        record.cigar().write_readable(&mut readable);
        }
        let readable_string = String::from_utf8_lossy(&readable);
        let readable_string2 = readable_string.replace("D", "Z");
        let readable_string3 = readable_string2.replace("I", "D");
        let readable_string4 =readable_string3.replace("Z", "I");
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

    let stdout = BufReader::new(process.stdout.unwrap());
    let mut reader = bam::SamReader::from_stream(stdout).unwrap();
    
    println!(">{}", name.to_string());
    let mut record = bam::Record::new();
    loop {
    // reader: impl RecordReader
    // New record is saved into record.
    match reader.read_into(&mut record) {
        // No more records to read.
        Ok(false) => break,
        Ok(true) => {
            if record.name() == name_vec.as_slice() {
                record.sequence().write_readable(&mut io::stdout());
                println!("");
            }
        },
        Err(e) => panic!("{}", e),
    }
    // Do somethind with the record.
    }
    /*
    
    for column in bam::Pileup::with_filter(&mut reader, |record| record.flag().no_bits(1796)) {
        let column = column.unwrap();
        //println!("Column at {}:{}, {} records", column.ref_id(),
        //    column.ref_pos() + 1, column.entries().len());

        for entry in column.entries().iter() {
            let seq: Vec<_> = entry.sequence().unwrap().map(|nt| nt as char).collect();
            //let qual: Vec<_> = entry.qualities().unwrap().iter()
            //    .map(|q| (q + 33) as char).collect();
            //println!("    {:?}: {:?}, {:?}", entry.record(), seq, qual);
            let unique_elements = seq.iter().cloned().unique().collect_vec();
            let mut unique_frequency = vec![];
            for unique_elem in unique_elements.iter() {
                unique_frequency.push((
                    seq.iter().filter(|&elem| elem == unique_elem).count(),
                    unique_elem,
                ));
            }
            unique_frequency.sort_by_key(|t| t.0);
            print!("{}", unique_frequency[0].1);
        }
        println!("");
    }
    */
}

fn main() {
    let args: Vec<String> = env::args().collect();
    // We assume that input is sorted by "samtools sort -n".
    //let reader = bam::IndexedReader::from_path(args[1].clone()).unwrap();
    let bam_stream = BufReader::with_capacity(1000000, File::open(args[1].clone()).unwrap());
    //let reader = bam::BamReader::from_path(args[1].clone()).unwrap();
    let reader = bam::BamReader::from_stream(bam_stream, 2).unwrap();
    let x = args.get(4)
        .and_then(|a| a.parse::<usize>().ok())
        .unwrap_or(5usize);
    /*for bin in reader.index().references()[0].bins().values() {
        println!("{}\t{}", bin.bin_id(), bin.chunks().len());
    }
    

    println!("next");
    println!("{}", reader.index());*/
    let mut fasta_reader = bio::io::fasta::IndexedReader::from_file(&args[2]).unwrap();
    // let mut read_tree = BTreeMap::new();
    // let mut previous_name: &[u8] = &[];
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
        if let Some(ref_name) = ref_name {
            fasta_reader.fetch(
                ref_name,
                record.start() as u64,
                record.calculate_end() as u64,
            );
            fasta_reader.read(&mut ref_seq);
        }
        // eprintln!("{:?} {}", record, record.flag().is_supplementary());
        if previous_name == record.name() {
            if !record.flag().is_supplementary() {
                //read_tree.insert(record);
                primary.push((record, ref_name, ref_seq));
            }
        } else {
            // let closure = |x: u32| reader.header().reference_name(x);
            if primary.len() > x {
                calculate_primary(primary, previous_name, &args[3]);
            } else {
                println!(">{}", String::from_utf8_lossy(&previous_name));
                for (record,_,_) in primary {
                    if record.sequence().len() > 0 {
                        record.sequence().write_readable(&mut io::stdout());
                    }
                }
                println!("");
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
    if primary.len() > 5 {
        calculate_primary(primary, previous_name, &args[3]);
    }

    /*
    if env::var_os("VIRTUAL_ENV").is_some() {
        println!("a")
    } else {
        println!("b")
    }
    let k: C<B> = C::new();
    */
    //let k2 = D{data: B};
}
