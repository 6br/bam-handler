use libbigwig::*;
use lib::range::{Region}
use std::ffi::{CString, CStr};

fn libbigwig_simple(feature: &ConfigFeature, coord: &Region, prefix: String) -> Vec<Feature> {
    let path = &feature.url;
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
            coord.stop as u32
        );
        if !intervals.is_null() {
            for i in 0..(*intervals).l {
                let start_offset =
                    *(*intervals).start.offset(i as isize) as u64;
                let stop_offset =
                    *(*intervals).end.offset(i as isize) as u64;
                let value = *(*intervals).value.offset(i as isize);
                vec.push(Feature {
                    start_offset: start_offset,
                    stop_offset: stop_offset,
                    id: i as u64,
                    name: feature.name.clone(),
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

fn libbigbed(feature: &ConfigFeature, coord: &Region, prefix: String) -> Vec<Feature> {
    let path = &feature.url;
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
                let start_offset =
                    if *(*intervals).start.offset(i as isize) <= coord.start as u32 {
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
                    name: feature.name.clone(),
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