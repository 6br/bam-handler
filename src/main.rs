#![crate_type="staticlib"]

#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate serde_json;

extern crate bitpacking;

use bitpacking::{BitPacker4x, BitPacker, UnsafeBitPacker};
//use lib::range::Region;
mod lib;

fn main() {
    println!("Hello, world!");
}

#[no_mangle]
pub extern fn hello_rust() -> *const u8 {
    "Hello, world!\0".as_ptr()
}

#[no_mangle]
pub extern fn bit_packing(data: &Vec<u32>) -> Vec<u8> {
    let mut original = vec![0u32, data.len()];
    original.copy_from_slice(data);

    let mut compressed = vec![0u8; (UnsafeBitPacker::BLOCK_LEN as usize) * 4];
    let numbits = UnsafeBitPacker::num_bits(&original[..]);

    UnsafeBitPacker::compress(&original[..], &mut compressed[..], numbits);

    return compressed;
}

#[no_mangle]
pub extern fn flip_region(region: Region) -> Region {

}

/*
  TEST CASES

1. Rust の構造体を他から呼び出しできるか？
2. flate2 を外部から呼び出しできるか？
3. bigwig をパースできるか？
4. Pfor を外部呼び出しで動作させられるか？ ->  bit_packing
5. 

*/

#[cfg(test)]
mod tests {
    use super::libbigbed;
    use lib::{Region, ConfigFeature};
    use features::{Feature};

    #[test]
    fn it_doesnot_work() {
        let vec: Vec<Feature> = vec![];
        assert_eq!(
            vec,
            libbigbed(
                &ConfigFeature{name: "feature".to_owned(), url: "test/ensGene.bb".to_owned(), chr_prefix: None, viz: None},
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
                &ConfigFeature{name: "test/ensGene.bb".to_owned(), url: "test/ensGene.bb".to_owned(), chr_prefix: None, viz: None},
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