extern crate bam;

use bam::Record;
use bam::{header::HeaderEntry, Header, RecordWriter};
use io::BufReader;
use itertools::Itertools;
use std::env;
use std::{
    collections::BTreeMap,
    fs::File,
    io,
    process::{Command, Stdio},
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

fn calculate_primary<'a, F>(primary: Vec<Record>, fasta_reader: bio::io::fasta::IndexedReader<File>, ref_name: &str, lambda: F) where
    F: Fn(usize) -> Option<&'a str> {
    let output = io::BufWriter::new(io::stdout());
    let mut header = Header::new();
    let name = String::from_utf8_lossy(primary[0].name());
    header.push_entry(HeaderEntry::ref_sequence(
        name.to_string(),
        primary[0].query_len(),
    ));

    let mut writer = bam::SamWriter::build().from_stream(output, header).unwrap();
    for record in primary {
        let mut readable: Vec<u8> = Vec::new();
        record.cigar().write_readable(&mut readable);
        let readable_string = String::from_utf8_lossy(&readable);
        readable_string.replace("D", "Z");
        readable_string.replace("I", "D");
        readable_string.replace("Z", "I");
        record.cigar().clear();
        record.cigar().extend_from_text(readable_string.to_string());

        let mut ref_seq = vec![];
        fasta_reader.fetch(lambda(record.ref_id() as usize).unwrap(), record.start() as u64, record.calculate_end() as u64);
        fasta_reader.read(&mut ref_seq);
        record.sequence().clear();
        record.sequence().extend_from_text(ref_seq);

        writer.write(&record).unwrap();
    }

    let process = match Command::new("realigner")
        .args(&[name.to_string()])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
    {
        Err(why) => panic!("couldn't spawn realigner: {}", why),
        Ok(process) => process,
    };

    let stdout = BufReader::new(process.stdout.unwrap());
    let mut reader = bam::SamReader::from_stream(stdout).unwrap();
    println!(">{}", name.to_string());
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
                unique_frequency.push((seq.iter().filter(|&elem| elem == unique_elem).count(), unique_elem));
            }
            unique_frequency.sort_by_key(|t | t.0);
            print!("{}", unique_frequency[0].1);
        }
        println!("");
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    // We assume that input is sorted by "samtools sort -n".
    //let reader = bam::IndexedReader::from_path(args[1].clone()).unwrap();
    let bam_stream = BufReader::with_capacity(1000000, File::open(args[1].clone()).unwrap());
    //let reader = bam::BamReader::from_path(args[1].clone()).unwrap();
    let reader = bam::BamReader::from_stream(bam_stream, 2).unwrap();
    /*for bin in reader.index().references()[0].bins().values() {
        println!("{}\t{}", bin.bin_id(), bin.chunks().len());
    }

    println!("next");
    println!("{}", reader.index());*/
    let fasta_reader = bio::io::fasta::IndexedReader(args[2].clone());
    // let mut read_tree = BTreeMap::new();
    let mut previous_name = [];
    let mut primary = vec![];
    for record in reader {
        let record = record.unwrap();
        if previous_name != record.name() {
            if !record.flag().is_supplementary() {
                //read_tree.insert(record);
                primary.push(record);
            } else {
                let closure = |x: &str| reader.header().reference_id(x);

                calculate_primary(primary, fasta_reader, closure);
                primary = vec![record];
            }
        }
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
