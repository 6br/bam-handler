// Copyright 2014-2016 Johannes Köster, Christopher Schröder, Henning Timm.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! FASTA format reading and writing.
//!
//! # Example
//!
//! ```
//! use std::io;
//! use bio::io::fasta;
//! let reader = fasta::Reader::new(io::stdin());
//! ```

use std::borrow::Borrow;
use std::cmp::min;
use std::collections;
use std::convert::AsRef;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::path::Path;

use csv;

/// Type alias for an owned text, i.e. ``Vec<u8>``.
pub type Text = Vec<u8>;
/// Type alias for a text slice, i.e. ``&[u8]``.
pub type TextSlice<'a> = &'a [u8];

/// Type alias for an iterator over a sequence, i.e. ``Iterator<Item=&u8>``.
pub trait TextIterator<'a>: Iterator<Item = &'a u8> {}
impl<'a, I: Iterator<Item = &'a u8>> TextIterator<'a> for I {}

/// Type alias for a type that can be coerced into a `TextIterator`.
/// This includes ``&Vec<u8>``, ``&[u8]``, ``Iterator<Item=&u8>``.
pub trait IntoTextIterator<'a>: IntoIterator<Item = &'a u8> {}
impl<'a, T: IntoIterator<Item = &'a u8>> IntoTextIterator<'a> for T {}

/// Remove a trailing newline from the given string in place.
pub fn trim_newline(s: &mut String) {
    if s.ends_with('\n') {
        s.pop();
    }
}

/// Maximum size of temporary buffer used for reading indexed FASTA files.
const MAX_FASTA_BUFFER_SIZE: usize = 512;

/// A FASTA reader.
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    reader: io::BufReader<R>,
    line: String,
}

impl Reader<fs::File> {
    /// Read FASTA from given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::new)
    }
}

impl<R: io::Read> Reader<R> {
    /// Create a new Fasta reader given an instance of `io::Read`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// let reader = Reader::new(fasta_file);
    /// # }
    /// ```
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line: String::new(),
        }
    }

    /// Read next FASTA record into the given `Record`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # use bio::io::fasta::Record;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// # let mut reader = Reader::new(fasta_file);
    /// let mut record = Record::new();
    /// reader.read(&mut record);
    ///
    /// assert_eq!(record.id(), "id");
    /// assert_eq!(record.desc().unwrap(), "desc");
    /// assert_eq!(record.seq().to_vec(), b"AAAA");
    /// # }
    /// ```
    pub fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        if self.line.is_empty() {
            r#try!(self.reader.read_line(&mut self.line));
            if self.line.is_empty() {
                return Ok(());
            }
        }

        if !self.line.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at record start.",
            ));
        }
        record.id = self.line[1..]
            .trim_right()
            .splitn(2, ' ')
            .nth(0)
            .map(|s| s.to_owned())
            .unwrap();
        record.desc = self.line[1..]
            .trim_right()
            .splitn(2, ' ')
            .nth(1)
            .map(|s| s.to_owned());
        loop {
            self.line.clear();
            r#try!(self.reader.read_line(&mut self.line));
            if self.line.is_empty() || self.line.starts_with('>') {
                break;
            }
            record.seq.push_str(self.line.trim_right());
        }

        Ok(())
    }

    /// Return an iterator over the records of this Fasta file.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # use bio::io::fasta::Record;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// # let reader = Reader::new(fasta_file);
    /// for record in reader.records() {
    ///     let record = record.unwrap();
    ///     assert_eq!(record.id(), "id");
    ///     assert_eq!(record.desc().unwrap(), "desc");
    ///     assert_eq!(record.seq().to_vec(), b"AAAA");
    /// }
    /// # }
    /// ```
    pub fn records(self) -> Records<R> {
        Records {
            reader: self,
            error_has_occured: false,
        }
    }
}

/// A FASTA index as created by SAMtools (.fai).
#[derive(Debug, Clone)]
pub struct Index {
    inner: Vec<IndexRecord>,
    name_to_rid: collections::HashMap<String, usize>,
}

impl Index {
    /// Open a FASTA index from a given `io::Read` instance.

    pub fn new<R: io::Read>(fai: R) -> csv::Result<Self> {
        let mut inner = vec![];
        let mut name_to_rid = collections::HashMap::new();

        let mut fai_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(fai);
        for (rid, row) in fai_reader.records().enumerate() {
            let string_record = row?;
            let record: IndexRecord = IndexRecord {
                name: string_record[0].to_string(),
                len: string_record[1].parse().unwrap(),
                offset: string_record[2].parse().unwrap(),
                line_bases: string_record[3].parse().unwrap(),
                line_bytes: string_record[4].parse().unwrap(),
            };
            name_to_rid.insert(record.name.clone(), rid);
            inner.push(record);
        }
        Ok(Index { inner, name_to_rid })
    }

    /// Open a FASTA index from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: &P) -> csv::Result<Self> {
        fs::File::open(path)
            .map_err(csv::Error::from)
            .and_then(Self::new)
    }

    /// Open a FASTA index given the corresponding FASTA file path.
    /// That is, for ref.fasta we expect ref.fasta.fai.
    pub fn with_fasta_file<P: AsRef<Path>>(fasta_path: &P) -> csv::Result<Self> {
        let mut fai_path = fasta_path.as_ref().as_os_str().to_owned();
        fai_path.push(".fai");

        Self::from_file(&fai_path)
    }

    /// Return a vector of sequences described in the index.
    pub fn sequences(&self) -> Vec<Sequence> {
        // sort kv pairs by rid to preserve order
        self.inner
            .iter()
            .map(|record| Sequence {
                name: record.name.clone(),
                len: record.len,
            })
            .collect()
    }
}

/// A FASTA reader with an index as created by SAMtools (.fai).
#[derive(Debug)]
pub struct IndexedReader<R: io::Read + io::Seek> {
    reader: io::BufReader<R>,
    pub index: Index,
    fetched_idx: Option<IndexRecord>,
    start: Option<u64>,
    stop: Option<u64>,
}

impl IndexedReader<fs::File> {
    /// Read from a given file path. This assumes the index ref.fasta.fai to be
    /// present for FASTA ref.fasta.
    pub fn from_file<P: AsRef<Path>>(path: &P) -> csv::Result<Self> {
        let index = Index::with_fasta_file(path)?;
        fs::File::open(path)
            .map(|f| Self::with_index(f, index))
            .map_err(csv::Error::from)
    }
}

impl<R: io::Read + io::Seek> IndexedReader<R> {
    /// Read from a FASTA and its index, both given as `io::Read`. FASTA has to
    /// be `io::Seek` in addition.
    pub fn new<I: io::Read>(fasta: R, fai: I) -> csv::Result<Self> {
        let index = r#try!(Index::new(fai));
        Ok(IndexedReader {
            reader: io::BufReader::new(fasta),
            index,
            fetched_idx: None,
            start: None,
            stop: None,
        })
    }

    /// Read from a FASTA and its index, the first given as `io::Read`, the
    /// second given as index object.
    pub fn with_index(fasta: R, index: Index) -> Self {
        IndexedReader {
            reader: io::BufReader::new(fasta),
            index,
            fetched_idx: None,
            start: None,
            stop: None,
        }
    }

    /// Fetch an interval from the sequence with the given name for reading.
    /// (stop position is exclusive).
    pub fn fetch(&mut self, seq_name: &str, start: u64, stop: u64) -> io::Result<()> {
        let idx = self.idx(seq_name)?;
        self.start = Some(start);
        self.stop = Some(stop);
        self.fetched_idx = Some(idx);
        Ok(())
    }

    /// Fetch an interval from the sequence with the given record index for reading.
    /// (stop position is exclusive).
    pub fn fetch_by_rid(&mut self, rid: usize, start: u64, stop: u64) -> io::Result<()> {
        let idx = self.idx_by_rid(rid)?;
        self.start = Some(start);
        self.stop = Some(stop);
        self.fetched_idx = Some(idx);
        Ok(())
    }

    /// Fetch the whole sequence with the given name for reading.
    pub fn fetch_all(&mut self, seq_name: &str) -> io::Result<()> {
        let idx = self.idx(seq_name)?;
        self.start = Some(0);
        self.stop = Some(idx.len);
        self.fetched_idx = Some(idx);
        Ok(())
    }

    /// Fetch the whole sequence with the given record index for reading.
    pub fn fetch_all_by_rid(&mut self, rid: usize) -> io::Result<()> {
        let idx = self.idx_by_rid(rid)?;
        self.start = Some(0);
        self.stop = Some(idx.len);
        self.fetched_idx = Some(idx);
        Ok(())
    }

    /// Read the fetched sequence into the given vector.
    pub fn read(&mut self, seq: &mut Text) -> io::Result<()> {
        let idx = self.fetched_idx.clone();
        match (idx, self.start, self.stop) {
            (Some(idx), Some(start), Some(stop)) => self.read_into_buffer(idx, start, stop, seq),
            _ => Err(io::Error::new(
                io::ErrorKind::Other,
                "No sequence fetched for reading.",
            )),
        }
    }

    /// Return an iterator yielding the fetched sequence.
    pub fn read_iter(&mut self) -> io::Result<IndexedReaderIterator<R>> {
        let idx = self.fetched_idx.clone();
        match (idx, self.start, self.stop) {
            (Some(idx), Some(start), Some(stop)) => self.read_into_iter(idx, start, stop),
            _ => Err(io::Error::new(
                io::ErrorKind::Other,
                "No sequence fetched for reading.",
            )),
        }
    }

    fn read_into_buffer(
        &mut self,
        idx: IndexRecord,
        start: u64,
        stop: u64,
        seq: &mut Text,
    ) -> io::Result<()> {
        if stop > idx.len {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "FASTA read interval was out of bounds",
            ));
        } else if start > stop {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Invalid query interval",
            ));
        }

        let mut bases_left = stop - start;
        let mut line_offset = self.seek_to(&idx, start)?;

        seq.clear();
        while bases_left > 0 {
            bases_left -= self.read_line(&idx, &mut line_offset, bases_left, seq)?;
        }

        Ok(())
    }

    fn read_into_iter(
        &mut self,
        idx: IndexRecord,
        start: u64,
        stop: u64,
    ) -> io::Result<IndexedReaderIterator<R>> {
        if stop > idx.len {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "FASTA read interval was out of bounds",
            ));
        } else if start > stop {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Invalid query interval",
            ));
        }

        let bases_left = stop - start;
        let line_offset = self.seek_to(&idx, start)?;
        let capacity = min(
            MAX_FASTA_BUFFER_SIZE,
            min(bases_left, idx.line_bases) as usize,
        );

        Ok(IndexedReaderIterator {
            reader: self,
            record: idx,
            bases_left,
            line_offset,
            buf: Vec::with_capacity(capacity),
            buf_idx: 0,
        })
    }

    /// Return the IndexRecord for the given sequence name or io::Result::Err
    fn idx(&self, seqname: &str) -> io::Result<IndexRecord> {
        match self.index.name_to_rid.get(seqname) {
            Some(rid) => self.idx_by_rid(*rid),
            None => Err(io::Error::new(
                io::ErrorKind::Other,
                "Unknown sequence name.",
            )),
        }
    }

    /// Return the IndexRecord for the given record index or io::Result::Err
    fn idx_by_rid(&self, rid: usize) -> io::Result<IndexRecord> {
        match self.index.inner.get(rid) {
            Some(record) => Ok(record.clone()),
            None => Err(io::Error::new(
                io::ErrorKind::Other,
                "Invalid record index in fasta file.",
            )),
        }
    }

    /// Seek to the given position in the specified FASTA record. The position
    /// of the cursor on the line that the seek ended on is returned.
    fn seek_to(&mut self, idx: &IndexRecord, start: u64) -> io::Result<u64> {
        assert!(start <= idx.len);

        let line_offset = start % idx.line_bases;
        let line_start = start / idx.line_bases * idx.line_bytes;
        let offset = idx.offset + line_start + line_offset;
        r#try!(self.reader.seek(io::SeekFrom::Start(offset)));

        Ok(line_offset)
    }

    /// Tries to read up to `bases_left` bases from the current line into `buf`,
    /// returning the actual number of bases read. Depending on the amount of
    /// whitespace per line, the current `line_offset`, and the amount of bytes
    /// returned from `BufReader::fill_buf`, this function may return Ok(0)
    /// multiple times in a row.
    fn read_line(
        &mut self,
        idx: &IndexRecord,
        line_offset: &mut u64,
        bases_left: u64,
        buf: &mut Vec<u8>,
    ) -> io::Result<u64> {
        let (bytes_to_read, bytes_to_keep) = {
            let src = self.reader.fill_buf()?;
            if src.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "FASTA file is truncated.",
                ));
            }

            let bases_on_line = idx.line_bases - min(idx.line_bases, *line_offset);
            let bases_in_buffer = min(src.len() as u64, bases_on_line);

            let (bytes_to_read, bytes_to_keep) = if bases_in_buffer <= bases_left {
                let bytes_to_read = min(src.len() as u64, idx.line_bytes - *line_offset);

                (bytes_to_read, bases_in_buffer)
            } else {
                (bases_left, bases_left)
            };

            buf.extend_from_slice(&src[..bytes_to_keep as usize]);
            (bytes_to_read, bytes_to_keep)
        };

        self.reader.consume(bytes_to_read as usize);

        assert!(bytes_to_read > 0);
        *line_offset += bytes_to_read;
        if *line_offset >= idx.line_bytes {
            *line_offset = 0;
        }

        Ok(bytes_to_keep)
    }
}

/// Record of a FASTA index.
#[derive(Debug, Clone)]
struct IndexRecord {
    name: String,
    len: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

/// A sequence record returned by the FASTA index.
#[derive(Debug, PartialEq)]
pub struct Sequence {
    pub name: String,
    pub len: u64,
}

pub struct IndexedReaderIterator<'a, R: io::Read + io::Seek + 'a> {
    reader: &'a mut IndexedReader<R>,
    record: IndexRecord,
    bases_left: u64,
    line_offset: u64,
    buf: Vec<u8>,
    buf_idx: usize,
}

impl<'a, R: io::Read + io::Seek + 'a> IndexedReaderIterator<'a, R> {
    fn fill_buffer(&mut self) -> io::Result<()> {
        assert!(self.bases_left > 0);

        self.buf.clear();
        let bases_to_read = min(self.buf.capacity() as u64, self.bases_left);

        // May loop one or more times; see IndexedReader::read_line.
        while self.buf.is_empty() {
            self.bases_left -= self.reader.read_line(
                &self.record,
                &mut self.line_offset,
                bases_to_read,
                &mut self.buf,
            )?;
        }

        self.buf_idx = 0;
        Ok(())
    }
}

impl<'a, R: io::Read + io::Seek + 'a> Iterator for IndexedReaderIterator<'a, R> {
    type Item = io::Result<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buf_idx < self.buf.len() {
            let item = Some(Ok(self.buf[self.buf_idx]));
            self.buf_idx += 1;
            item
        } else if self.bases_left > 0 {
            if let Err(e) = self.fill_buffer() {
                self.bases_left = 0;
                self.buf_idx = self.buf.len();

                return Some(Err(e));
            }

            self.buf_idx = 1;
            Some(Ok(self.buf[0]))
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let hint = self.bases_left as usize + (self.buf.len() - self.buf_idx);

        (hint, Some(hint))
    }
}

/// A Fasta writer.
#[derive(Debug)]
pub struct Writer<W: io::Write> {
    writer: io::BufWriter<W>,
}

impl Writer<fs::File> {
    /// Write to the given file path.
    pub fn to_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::create(path).map(Writer::new)
    }
}

impl<W: io::Write> Writer<W> {
    /// Create a new Fasta writer.
    pub fn new(writer: W) -> Self {
        Writer {
            writer: io::BufWriter::new(writer),
        }
    }

    /// Directly write a Fasta record.
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.write(record.id(), record.desc(), record.seq())
    }

    /// Write a Fasta record with given id, optional description and sequence.
    pub fn write(&mut self, id: &str, desc: Option<&str>, seq: TextSlice) -> io::Result<()> {
        r#try!(self.writer.write_all(b">"));
        r#try!(self.writer.write_all(id.as_bytes()));
        if desc.is_some() {
            r#try!(self.writer.write_all(b" "));
            r#try!(self.writer.write_all(desc.unwrap().as_bytes()));
        }
        r#try!(self.writer.write_all(b"\n"));
        r#try!(self.writer.write_all(seq));
        r#try!(self.writer.write_all(b"\n"));

        Ok(())
    }

    /// Flush the writer, ensuring that everything is written.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

/// A FASTA record.
#[derive(Default, Clone, Debug)]
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: String,
}

impl Record {
    /// Create a new instance.
    pub fn new() -> Self {
        Record {
            id: String::new(),
            desc: None,
            seq: String::new(),
        }
    }

    /// Create a new Fasta record from given attributes.
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: TextSlice) -> Self {
        let desc = match desc {
            Some(desc) => Some(desc.to_owned()),
            _ => None,
        };
        Record {
            id: id.to_owned(),
            desc,
            seq: String::from_utf8(seq.to_vec()).unwrap(),
        }
    }

    /// Check if record is empty.
    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_none() && self.seq.is_empty()
    }

    /// Check validity of Fasta record.
    pub fn check(&self) -> Result<(), &str> {
        if self.id().is_empty() {
            return Err("Expecting id for Fasta record.");
        }
        if !self.seq.is_ascii() {
            return Err("Non-ascii character found in sequence.");
        }

        Ok(())
    }

    /// Return the id of the record.
    pub fn id(&self) -> &str {
        self.id.as_ref()
    }

    /// Return descriptions if present.
    pub fn desc(&self) -> Option<&str> {
        match self.desc.as_ref() {
            Some(desc) => Some(&desc),
            None => None,
        }
    }

    /// Return the sequence of the record.
    pub fn seq(&self) -> TextSlice {
        self.seq.as_bytes()
    }

    /// Clear the record.
    fn clear(&mut self) {
        self.id.clear();
        self.desc = None;
        self.seq.clear();
    }
}

/// An iterator over the records of a Fasta file.
pub struct Records<R: io::Read> {
    reader: Reader<R>,
    error_has_occured: bool,
}

impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        if self.error_has_occured {
            None
        } else {
            let mut record = Record::new();
            match self.reader.read(&mut record) {
                Ok(()) if record.is_empty() => None,
                Ok(()) => Some(Ok(record)),
                Err(err) => {
                    self.error_has_occured = true;
                    Some(Err(err))
                }
            }
        }
    }
}

lazy_static! {
    static ref COMPLEMENT: Vec<u8> = {
        let mut comp = Vec::new();
        comp.resize(256, 0);
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in b"AGCUYRWSKMDVHBNZ".iter().zip(b"UCGARYWSMKHBDVNZ".iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}

/// Return complement of given RNA alphabet character (IUPAC alphabet supported).
pub fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

/// Calculate reverse complement of given text (IUPAC alphabet supported).
pub fn revcomp<C, T>(text: T) -> Vec<u8>
where
    C: Borrow<u8>,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator,
{
    text.into_iter()
        .rev()
        .map(|a| complement(*a.borrow()))
        .collect()
}
