#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate log;

extern crate serde_json;
extern crate libbigwig;
extern crate bitpacking;
extern crate regex;

pub mod bigbed;
pub mod range;

use bitpacking::{BitPacker4x, BitPacker};
use flate2::Compression;
use flate2::write::ZlibEncoder;
use bigbed::libbigbed;
use range::{Region};
use bigbed::{Feature};
use std::io::prelude::*;
/*
  TEST CASES
0. 普通に呼べるか? -> hello_rust()
1. Rust の構造体を他から呼び出しできるか？ -> decrement_start()
2. flate2 を外部から呼び出しできるか？ -> compress_bytes()
3. bigbed をパースできるか？ -> load_bigbed()
4. Pfor を外部呼び出しで動作させられるか？ ->  bit_packing()

*/

#[repr(u32)]
pub enum Foo {
    A = 1,
    B,
    C,
}

#[no_mangle]
pub extern fn hello_rust() -> *const u8 {
    "Hello, world!\0".as_ptr()
}


#[no_mangle]
pub unsafe extern "C" fn print_foo(foo: *const Foo) {
    println!(
        "{}",
        match *foo {
            Foo::A => "a",
            Foo::B => "b",
            Foo::C => "c",
        }
    );
}

#[no_mangle]
pub extern fn bit_packing(data: &[u32]) -> Vec<u8> {
    let mut original = vec![0u32, data.len() as u32];
    original.copy_from_slice(data);

    let bitpacker = BitPacker4x::new();
    let num_bits: u8 = bitpacker.num_bits(&original);

    let mut compressed = vec![0u8; 4 * BitPacker4x::BLOCK_LEN];
    let _compressed_len = bitpacker.compress(&original, &mut compressed[..], num_bits);

    return compressed
}

#[no_mangle]
pub extern fn decrement_start(mut region: Region) -> Region {
    region.start_minus();
    return region
}

#[no_mangle]
pub extern fn load_bigbed(path: String, region: Region) -> Vec<Feature> {
    return libbigbed(
        path,
        &region, 
        "".to_owned(),
    )
}

#[no_mangle]
pub extern fn compress_bytes(words: &[u8]) -> Result<Vec<u8>, std::io::Error> {
    let mut e = ZlibEncoder::new(Vec::new(), Compression::default());
    e.write_all(words)?;
    let compressed_bytes = e.finish();
    return compressed_bytes
}

#[cfg(test)]
mod tests {
    use super::libbigwig;
    use crate::range::{Region};
    use crate::bigbed::{Feature};
    use crate::bigbed::libbigbed;

    #[test]
    fn it_doesnot_work() {
        let vec: Vec<Feature> = vec![];
        assert_eq!(
            vec,
            libbigbed(
                "test/ensGene.bb".to_string(),
                &Region {
                    path: "Y".to_owned(),
                    start: 2712790,
                    stop: 2712894,
                },
                "".to_owned(),
            )
        );
    }

    #[test]
    fn it_works() {
        let raw_attr: String = "ENST00000387529\t0\t+\t2712894\t2712894\t0\t1\t104,\t0,\tENSG00000210264\tnull"
            .to_owned();
        let attr: Vec<String> = raw_attr.split("\t").map(|s| s.to_string()).collect();
        let feat: Feature = Feature {
            start_offset: 0,
            stop_offset: 0,
            id: 0,
            name: "test/ensGene.bb".to_owned(),
            is_reverse: None,
            value: None,
            attributes: attr,
        };
        let vec: Vec<Feature> = vec![feat];
        assert_eq!(
            vec,
            libbigbed(
                "test/ensGene.bb".to_string(),
                &Region {
                    path: "Y".to_owned(),
                    start: 2712790,
                    stop: 2712894,
                },
                "chr".to_owned(),
            )
        );
    }

}