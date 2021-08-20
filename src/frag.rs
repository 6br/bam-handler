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
        if iter.i >= 4 {
            vec.push(item);
        }
    }
    let mut array: Vec<(i32, i32)> = vec.into_iter().map(bin_to_region).collect();
    let array_length = array.len();
    array[0].0 = beg;
    array[array_length - 1].1 = end;
    array
}

/// Returns a BAI bin for the record with alignment `[beg-end)`.
/*
pub fn region_to_bin(beg: i32, end: i32) -> u32 {
    let end = end - 1;
    let mut res = 0_i32;
    for i in (14..27).step_by(3) {
        if beg >> i == end >> i {
            res = ((1 << (29 - i)) - 1) / 7 + (beg >> i);
            break;
        }
    }
    res as u32
}
*/
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
    pub i: i32,
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
            self.curr_bin = (self.t + (self.start >> (26 - 3 * self.i))) as u32 - 1;
            self.bins_end = (self.t + (self.end >> (26 - 3 * self.i))) as u32;

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
    Ok(())
}

impl std::iter::FusedIterator for BinsIter {}

pub fn frag(input_path: String, output_path: String, sequence_squash: bool) {
    let reader = bam::BamReader::from_path(input_path, 4).unwrap();
    let out_writer = File::create(&output_path).unwrap();
    let output = io::BufWriter::with_capacity(1048576, out_writer);
    //     let sequence_squash = true;
    let left_right_padding = if sequence_squash {
        Operation::Hard
    } else {
        Operation::Soft
    };

    let mut writer = bam::BamWriter::build()
        .write_header(true)
        .from_stream(output, reader.header().clone())
        .unwrap();

    for record in reader {
        let mut record = record.unwrap();
        let true_start = record.start();
        let true_end = record.calculate_end();
        if !record.flag().is_mapped() {
            writer.write(&record).unwrap();
            continue;
        }
        let ranges = ranges(true_start, true_end);
        record.tags_mut().push_num(b"TE", true_end as i32);
        record.tags_mut().push_num(b"TS", true_start as i32);
        record.tags_mut().remove(b"MD");

        let left_hard_clip = record.cigar().hard_clipping(true);
        let right_hard_clip = record.cigar().hard_clipping(false);
        let mut cigars: Vec<(u32, Operation)> = record
            .cigar()
            .iter()
            .filter(|t| t.1 != Operation::Hard)
            .collect();

        let mut current_pos = true_start;
        let sequence_len = record.sequence().len();
        let sequence = record.sequence();
        // let qualities = record.qualities();
        // assert!(sequence_len > 0, "{:?}", record.name());
        let mut consumed_query_len = 0; //left_hard_clip;
        let mut cigar_index = 0;
        //eprintln!("{:?}", cigars);
        if ranges.len() == 1 {
            writer.write(&record).unwrap();
        } else {
            for (beg, end) in ranges {
                let mut record = record.clone();
                record.set_start(beg);
                let mut cigar = Cigar::new();
                cigar.push(left_hard_clip, Operation::Hard);
                if consumed_query_len > 0 {
                    cigar.push(consumed_query_len, left_right_padding);
                }
                loop {
                    if cigars.len() <= cigar_index || current_pos == end {
                        break;
                    }
                    let cigar_len = match cigars[cigar_index].1 {
                        Operation::AlnMatch => cigars[cigar_index].0,
                        Operation::Insertion => 0,
                        Operation::Deletion => cigars[cigar_index].0,
                        Operation::Skip => 0,
                        Operation::Soft => 0,
                        Operation::Hard => 0,
                        Operation::Padding => 0,
                        Operation::SeqMatch => cigars[cigar_index].0,
                        Operation::SeqMismatch => cigars[cigar_index].0,
                    } as i32;
                    // eprintln!("{} {} {}", current_pos, cigar_len, end);

                    if current_pos + cigar_len > end {
                        cigar.push(
                            cigars[cigar_index].0 - (cigar_len - (end - current_pos)) as u32,
                            cigars[cigar_index].1,
                        );
                        cigars[cigar_index].0 = (cigar_len - (end - current_pos)) as u32;
                        current_pos += end - current_pos;
                    } else {
                        cigar.push(cigars[cigar_index].0, cigars[cigar_index].1);
                        current_pos += cigar_len;
                        cigar_index += 1;
                    }
                }
                let mut readable: Vec<u8> = Vec::new();
                cigar.write_readable(&mut readable).unwrap();
                //let readable_string = String::from_utf8_lossy(&readable);
                //eprintln!("{}", readable_string);
                //let bytes = readable_string.bytes().into_iter();
                record.set_cigar(readable).unwrap();
                let previous_consumed_query_len = consumed_query_len;
                consumed_query_len = record.cigar().calculate_query_len()
                    + if sequence_squash {
                        previous_consumed_query_len
                    } else {
                        0
                    };
                debug!(
                    "{} {} {} {}",
                    consumed_query_len,
                    sequence_len,
                    previous_consumed_query_len,
                    record.cigar().calculate_query_len()
                );
                assert!(previous_consumed_query_len < consumed_query_len);
                if sequence_len > 0 && sequence_len as u32 - consumed_query_len > 0 {
                    let remaining_query_len = (sequence_len as u32) - consumed_query_len;
                    cigar.push(remaining_query_len, left_right_padding);
                    debug!("{} {}", consumed_query_len, remaining_query_len);
                }
                cigar.push(right_hard_clip, Operation::Hard);
                let mut readable2: Vec<u8> = Vec::new();
                cigar.write_readable(&mut readable2).unwrap();
                let readable_string = String::from_utf8_lossy(&readable2);
                debug!("{}", readable_string);
                record.set_cigar(readable2).unwrap();

                if sequence_squash && sequence_len > 0 {
                    let subseq = sequence
                        .subseq(previous_consumed_query_len as usize..consumed_query_len as usize);
                    debug!(
                        "{} {} {} {}",
                        previous_consumed_query_len,
                        consumed_query_len,
                        sequence_len,
                        record.cigar().calculate_query_len(),
                    );
                    //let qualities = &qualities.raw()
                    //    [previous_consumed_query_len as usize..consumed_query_len as usize]; // //.collect();
                    record
                        .set_seq_qual(
                            subseq,
                            std::iter::empty(),
                            // (consumed_query_len as usize - previous_consumed_query_len as usize),
                            //                                * 2,
                        )
                        .unwrap();
                }
                assert!(
                    record.cigar().calculate_query_len() == record.sequence().len() as u32
                        || record.sequence().len() == 0,
                    "CIGAR len {}, record len {}",
                    record.cigar().calculate_query_len(),
                    record.sequence().len()
                );
                assert!(
                    end == record.start() + record.cigar().calculate_ref_len() as i32,
                    "Expected end {}, calculated end {}, true end {}",
                    end,
                    record.start() + record.cigar().calculate_ref_len() as i32,
                    true_end
                );
                writer.write(&record).unwrap();
            }
        }
        writer.finish().unwrap();
    }
}
