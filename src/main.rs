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