pub mod range;
pub mod bigbed;


/*
  TEST CASES
0. 普通に呼べるか? -> hello_rust()
1. Rust の構造体を他から呼び出しできるか？ -> decrement_start()
2. flate2 を外部から呼び出しできるか？ -> compress_bytes()
3. bigbed をパースできるか？ -> load_bigbed()
4. Pfor を外部呼び出しで動作させられるか？ ->  bit_packing()

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