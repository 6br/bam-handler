#![feature(rustc_private)]
#![feature(plugin)]

#![crate_type="staticlib"]

#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate log;
extern crate serde_json;
extern crate libbigwig;
extern crate bitpacking;
extern crate regex;

mod lib;

use bitpacking::{BitPacker4x, BitPacker};
use flate2::Compression;
use flate2::write::ZlibEncoder;
use lib::bigbed::libbigbed;
use lib::range::{Region};
use lib::bigbed::{Feature};
use std::io::prelude::*;

fn main() {
    println!("Hello, world!");
}

/*
  TEST CASES
0. 普通に呼べるか? -> hello_rust()
1. Rust の構造体を他から呼び出しできるか？ -> start_decrement()
2. flate2 を外部から呼び出しできるか？ -> flate()
3. bigbed をパースできるか？ -> load_bigbed()
4. Pfor を外部呼び出しで動作させられるか？ ->  bit_packings()

*/

#[no_mangle]
pub extern fn hello_rust() -> *const u8 {
    "Hello, world!\0".as_ptr()
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
pub extern fn start_decrement(mut region: Region) -> Region {
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
pub extern fn flate(words: &[u8]) -> Result<Vec<u8>, std::io::Error> {
    let mut e = ZlibEncoder::new(Vec::new(), Compression::default());
    e.write_all(words)?;
    let compressed_bytes = e.finish();
    return compressed_bytes
}

#[cfg(test)]
mod tests {
    use super::libbigwig;
    use crate::lib::range::{Region};
    use crate::lib::bigbed::{Feature};
    use crate::lib::bigbed::libbigbed;

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