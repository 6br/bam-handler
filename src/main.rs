#![feature(rustc_private)]
#![feature(plugin)]
extern crate bam;

use std::env;
use std::io;

use bam::RecordWriter;

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
    let reader = bam::IndexedReader::from_path(args[1].clone()).unwrap();
    for bin in reader.index().references()[0].bins().values() {
        println!("{}\t{}", bin.bin_id(), bin.chunks().len());
    }
    println!("next");
    println!("{}", reader.index());


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

