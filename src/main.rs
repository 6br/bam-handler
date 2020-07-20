#![feature(rustc_private)]
#![feature(plugin)]
extern crate bam;

use std::env;
use std::io;
use std::collections::BTreeMap;

trait Test {
    fn new() -> Self;
}

struct A {
    data: Vec<u32>
}

struct B {
    data: Vec<i32>
}

impl Test for A {
    fn new() -> A {
        return A{data: vec![]}
    }
}

impl Test for B {
    fn new() -> B {
        return B{data: vec![]}
    }
}

impl Test3 for B {
    fn new(&mut self) -> Option<bool>{
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
    data: Test3
}
/*
impl D {
    pub fn new() -> D {
        D{data: Test3{}}
    }
}
*/
struct C<T: Test> {
    data: T
}

impl<T:Test> C<T> {
    pub fn new() -> C<T> {
        C{data: T::new()}
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut reader = bam::IndexedReader::from_path(args[1].clone()).unwrap();
    for bin in reader.index().references()[0].bins().values() {
    //    println!("{}\t{}", bin.bin_id(), bin.chunks().len());
    }
    //println!("next");
    //println!("{}", reader.index());

    let region = [("chr1", 121_700_000, 123_400_000), ("chr5", 46100000, 51400000), ("chr19", 24200000, 28100000)];
    let margin = 1_000_000;
    let result = vec![args[1].clone()];

    let number: i32 = match &args[2].parse() {
                Ok(n) => {
                    *n
                },
                Err(_) => {
                    1
                },
            };
    
    for i in region.iter() {
    let chr = i.0; // "chr1";
    let id = reader.header().reference_id(chr).unwrap();
    let mut btree = BTreeMap::<u32, (i32,i32,i32,u32,u32,u32,String)>::new();
    let mut btree2 = BTreeMap::<i32, (i32,i32,i32,u32,u32,u32,String)>::new();
    for record in reader.fetch(&bam::Region::new(id, i.1 - margin, i.2 + margin)).unwrap() {
        let record = record.unwrap();
                    // REF_START, REF_END, REF_LEN, QRY_START, QRY_END, QRY_LEN
        let tuple = (record.start(), record.calculate_end(), record.calculate_end() - record.start(), record.aligned_query_start(), record.aligned_query_end(), record.query_len(), std::str::from_utf8(&record.name()).unwrap_or("*NAME NOT UTF-8*").to_string());
        btree.insert(tuple.5, tuple.clone());
        btree2.insert(tuple.2, tuple);
    }
    let mut index = 0;
    // result.push(write!("READ_LEN: {}, REF {}:{}-{} (len:{})", key, chr, value.0, value.1, value.2, ))
    for (key, value) in btree.iter().rev() {
        index += 1;
        println!("{}:, READ_LEN: {}, REF {}:{}-{} (len:{}) {}", args[1], key, chr, value.0, value.1, value.2, value.6);
        if (index > number) {
            break;
        }
    }  
    index = 0;
    for (key, value) in btree2.iter().rev() {
        index += 1;
        println!("{}:, READ_LEN: {}, REF {}:{}-{} (len:{}) {}", args[1], key, chr, value.0, value.1, value.2, value.6);
        if (index > number) {
            break;
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

