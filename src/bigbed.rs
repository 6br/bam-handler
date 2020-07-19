use super::range::Region;
use libbigwig::*;
use std::ffi::{CStr, CString};

#[derive(Debug, PartialEq, Serialize, Deserialize)]
pub struct Feature {
    pub start_offset: u64,
    pub stop_offset: u64,
    pub id: u64,
    pub name: String,
    pub is_reverse: Option<bool>,
    pub attributes: Vec<String>,
    pub value: Option<f32>,
}

pub fn libbigwig_simple(path: String, coord: &Region, prefix: String) -> Vec<Feature> {
    let path_loc = CString::new(path.clone()).unwrap();
    let read_only = CString::new("r").unwrap();
    let path_str = CString::new(prefix + coord.path.as_ref()).unwrap();
    let mut vec: Vec<Feature> = vec![];
    unsafe {
        if bwInit(1 << 17) != 0 {
            return vec;
        }
        let fp = bwOpen(path_loc.into_raw(), None, read_only.into_raw());
        if fp.is_null() {
            return vec;
        }
        //let intervals = bwGetValues(
        let intervals = bwGetOverlappingIntervals(
            fp,
            path_str.into_raw(),
            coord.start as u32,
            coord.stop as u32,
        );
        if !intervals.is_null() {
            for i in 0..(*intervals).l {
                let start_offset = *(*intervals).start.offset(i as isize) as u64;
                let stop_offset = *(*intervals).end.offset(i as isize) as u64;
                let value = *(*intervals).value.offset(i as isize);
                vec.push(Feature {
                    start_offset: start_offset,
                    stop_offset: stop_offset,
                    id: i as u64,
                    name: path.clone(),
                    attributes: vec![],
                    is_reverse: None,
                    value: Some(value),
                });
            }
        }

        if intervals.is_null() {
            bwDestroyOverlappingIntervals(intervals);
        }
        bwClose(fp);
        bwCleanup();
    }
    return vec;
}

pub fn libbigbed(path: String, coord: &Region, prefix: String) -> Vec<Feature> {
    debug!("{:?} {:?}", path, coord);
    let path_loc = CString::new(path.clone()).unwrap();
    let path_str = CString::new(prefix + coord.path.as_ref()).unwrap();
    let mut vec: Vec<Feature> = vec![];
    unsafe {
        if bwInit(1 << 17) != 0 {
            return vec;
        }
        let fp = bbOpen(path_loc.into_raw(), None);
        if fp.is_null() {
            return vec;
        }
        let intervals = bbGetOverlappingEntries(
            fp,
            path_str.into_raw(),
            coord.start as u32,
            coord.stop as u32,
            1,
        );
        if !intervals.is_null() {
            for i in 0..(*intervals).l {
                let start_offset = if *(*intervals).start.offset(i as isize) <= coord.start as u32 {
                    0
                } else {
                    *(*intervals).start.offset(i as isize) as u64 - coord.start
                };
                let stop_offset = if coord.stop <= *(*intervals).end.offset(i as isize) as u64 {
                    0
                } else {
                    coord.stop - *(*intervals).end.offset(i as isize) as u64
                };
                let name = CStr::from_ptr(*(*intervals).str.offset(i as isize))
                    .to_str()
                    .unwrap()
                    .to_owned();
                let splitted_attr = name.split("\t").map(|s| s.to_string()).collect();
                vec.push(Feature {
                    start_offset: start_offset,
                    stop_offset: stop_offset,
                    id: i as u64,
                    name: path.clone(),
                    attributes: splitted_attr,
                    is_reverse: None,
                    value: None,
                });
            }
        }
        if intervals.is_null() {
            bbDestroyOverlappingEntries(intervals);
        }
        bwClose(fp);
        bwCleanup();
    }
    return vec;
}
