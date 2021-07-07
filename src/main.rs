#[macro_use]
extern crate log;
extern crate bam;

extern crate protobuf;

mod protos;
use protos::vg::Alignment;

use bam::Record;
use bam::RecordReader;
use bam::{
    header::HeaderEntry,
    record::tags::{StringType, TagValue},
    Header, RecordWriter,
};
use bio::alphabets::dna::revcomp;
use byteorder::WriteBytesExt;
use io::{BufReader, Read};
use itertools::Itertools;
use protobuf::Message;
use regex::Regex;
use std::env;
use std::io::Write;
use std::{
    fs::File,
    io, iter,
    process::{Command, Stdio},
};
//use vg::Alignment;
//mod vg;
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

fn main() {
    let args: Vec<String> = env::args().collect();
    // We assume that input is sorted by "samtools sort -n".
    //let reader = bam::IndexedReader::from_path(args[1].clone()).unwrap();
    let mut pb_stream = BufReader::with_capacity(1000000, File::open(args[1].clone()).unwrap());
    //let reader = bam::BamReader::from_path(args[1].clone()).unwrap();
    //let reader = bam::BamReader::from_stream(bam_stream, 2).unwrap();
    let x = args
        .get(4)
        .and_then(|a| a.parse::<usize>().ok())
        .unwrap_or(5usize);
    let alpha = args
        .get(5)
        .and_then(|a| a.parse::<f64>().ok())
        .unwrap_or(0.85);
    let sigma = args
        .get(6)
        .and_then(|a| a.parse::<f64>().ok())
        .unwrap_or(0.075);
    let dummy = "dummy.fa".to_string();
    let path = args.get(2).unwrap_or(&dummy);

    let mut writer = bam::SamWriter::build()
        .write_header(false)
        .from_stream(std::io::stdout(), Header::new())
        .unwrap();
    // we can build a bytes reader directly out of the bytes

    let mut buffer = vec![];
    pb_stream.read_to_end(&mut buffer).unwrap();

    let mut person = Alignment::new();
    person.merge_from_bytes(&buffer[..]).unwrap();
    println!("{:?}", person)

    // now using the generated module decoding is as easy as:
    //let foobar = FooBar::from_reader(&mut reader, &bytes).expect("Cannot read FooBar");

    // if instead the buffer contains a length delimited stream of message we could use:
    //while !reader.is_eof() {
    //    let foobar: Alignment = reader.read_message(&buffer).expect("Cannot read Alignment");
    //    let record = Record::new();
    //    record.set_name(foobar.);

    //    writer.write(&record).unwrap();
    //}
}
