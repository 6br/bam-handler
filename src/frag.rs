use crate::bam::RecordWriter;
use bam::{
    index::bin_to_region,
    record::{cigar::Operation, Cigar},
};
use genomic_range::StringRegion;
use std::{
    fs::File,
    io::Error,
    io::{self, ErrorKind},
    time::Instant,
};

fn ranges(beg: i32, end: i32) -> Vec<(i32, i32)> {
    let mut iter = region_to_bins(beg, end);
    let mut vec = vec![];
    while let Some(item) = iter.next() {
        if iter.t >= 3 {
            vec.push(item);
        }
    }
    let mut array: Vec<(i32, i32)> = vec.into_iter().map(|t| bin_to_region(t)).collect();
    let array_length = array.len();
    array[0].0 = beg;
    array[array_length - 1].1 = end;
    return array;
}

/// Returns a BAI bin for the record with alignment `[beg-end)`.
pub fn region_to_bin(beg: i32, end: i32) -> u32 {
    let end = end - 1;
    let mut res = 0_i32;
    for i in (14..27).step_by(3) {
        if beg >> i == end >> i {
            res = ((1 << 29 - i) - 1) / 7 + (beg >> i);
            break;
        }
    }
    res as u32
}

/// Returns all possible BAI bins for the region `[beg-end)`.
pub fn region_to_bins(start: i32, end: i32) -> BinsIter {
    BinsIter {
        i: -1,
        t: 0,
        start,
        end,
        curr_bin: 0,
        bins_end: 0,
    }
}

/// Iterator over bins.
pub struct BinsIter {
    i: i32,
    pub t: i32,
    start: i32,
    end: i32,
    curr_bin: u32,
    bins_end: u32,
}

impl Iterator for BinsIter {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_bin == self.bins_end {
            if self.i >= 4 {
                return None;
            }
            self.i += 1;
            self.t += 1 << (self.i * 3);
            self.curr_bin = (self.t + (self.start >> 26 - 3 * self.i)) as u32 - 1;
            self.bins_end = (self.t + (self.end >> 26 - 3 * self.i)) as u32;

            if self.i == 0 {
                return Some(0);
            }
        }
        self.curr_bin += 1;
        Some(self.curr_bin)
    }
}

pub fn bench(input: String, range: String) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = bam::IndexedReader::from_path(input).unwrap();
    let lambda = |range: String| {
        eprintln!("{:?}", range);
        let ref_id = reader
            .header()
            .reference_id(&range)
            .ok_or_else(|| Error::new(ErrorKind::Other, "Invalid reference id."))?;
        let end = reader
            .header()
            .reference_len(ref_id)
            .ok_or_else(|| Error::new(ErrorKind::Other, "Invalid reference id."))?;
        Ok(StringRegion {
            path: range,
            start: 1,
            end: end as u64,
        })
    };
    let prefetch_range = StringRegion::new(&range)
        .or_else(
            |_| -> std::result::Result<genomic_range::StringRegion, Box<dyn std::error::Error>> {
                lambda(range)
            },
        )
        .unwrap();
    let ref_id = reader
        .header()
        .reference_id(prefetch_range.path.as_ref())
        .ok_or_else(|| Error::new(ErrorKind::Other, "Invalid reference id."))?;
    let ref_len = reader
        .header()
        .reference_len(ref_id)
        .ok_or_else(|| Error::new(ErrorKind::Other, "Invalid reference length."))?;
    let start = if prefetch_range.start == 0 {
        warn!("Invalid reference length: start == 0",);
        1
    } else {
        prefetch_range.start as u32
    };
    let end = if prefetch_range.end as u32 > ref_len {
        warn!(
            "Invalid reference length: end > reference length ({} > {})",
            prefetch_range.end as u32, ref_len
        );
        ref_len
    } else {
        prefetch_range.end as u32
    };
    let viewer = reader.fetch(&bam::bam_reader::Region::new(ref_id, start, end))?;

    let start = Instant::now();
    for record in viewer {
        record?.calculate_end();
    }
    let duration = start.elapsed();

    println!("Time elapsed in calculate_end() is: {:?}", duration);
    return Ok(());
}

impl std::iter::FusedIterator for BinsIter {}

pub fn frag(input: String, output_path: String) {
    let mut reader = bam::IndexedReader::from_path(input).unwrap();
    let out_writer = File::create(&output_path).unwrap();
    let output = io::BufWriter::with_capacity(1048576, out_writer);

    let mut writer = bam::BamWriter::build()
        .write_header(false)
        .from_stream(output, reader.header().clone())
        .unwrap();

    for record in reader.full() {
        let mut record = record.unwrap();
        let true_start = record.start();
        let true_end = record.calculate_end();
        let ranges = ranges(true_start, true_end);
        record.tags_mut().push_num(b"TE", true_end as i32);
        record.tags_mut().push_num(b"TS", true_start as i32);

        let mut cigars: Vec<(u32, Operation)> = record.cigar().iter().collect();
        let mut current_pos = true_start;
        let mut cigar_index = 0;

        for (beg, end) in ranges {
            let mut record = record.clone();
            record.set_start(beg);
            let mut cigar = Cigar::new();
            loop {
                if cigars.len() >= cigar_index {
                    break;
                }
                let cigar_len = match cigars[0].1 {
                    Operation::AlnMatch => cigars[0].0,
                    Operation::Insertion => 0,
                    Operation::Deletion => cigars[0].0,
                    Operation::Skip => 0,
                    Operation::Soft => 0,
                    Operation::Hard => 0,
                    Operation::Padding => 0,
                    Operation::SeqMatch => cigars[0].0,
                    Operation::SeqMismatch => cigars[0].0,
                } as i32;

                if current_pos + cigar_len > end {
                    cigar.push(
                        cigars[0].0 - (end - cigar_len - current_pos) as u32,
                        cigars[0].1,
                    );
                    cigars[0].0 -= (end - cigar_len - current_pos) as u32;
                    current_pos += end - cigar_len - current_pos;
                } else {
                    cigar.push(cigars[0].0, cigars[0].1);
                    current_pos += cigar_len;
                    cigar_index += 1;
                }
            }
            let mut readable: Vec<u8> = Vec::new();
            record.cigar().write_readable(&mut readable);
            record.set_cigar(readable);
            assert!(end == record.calculate_end());
            writer.write(&record).unwrap();
        }
    }
}
