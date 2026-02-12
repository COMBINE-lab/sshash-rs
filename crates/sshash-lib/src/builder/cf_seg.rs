//! Parser for cuttlefish `.cf_seg` segment files
//!
//! The `.cf_seg` format encodes the unitigs (segments) of a compacted de Bruijn
//! graph. Each line contains a numeric segment identifier followed by the unitig
//! DNA sequence, separated by whitespace:
//!
//! ```text
//! 47966863    CGCACATCCGTATCATGAAG...
//! 19947264    TATGAGGGTGGGAAGGTTGC...
//! ```
//!
//! Segment IDs are unique within a file but are neither ordered nor contiguous.
//! The sequences are the spectrum-preserving string set (SPSS) over which
//! sshash builds its dictionary.

use anyhow::{Context, Result, bail};
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Parsed contents of a `.cf_seg` file.
///
/// The `i`-th entry in `segment_ids` corresponds to the `i`-th entry in
/// `sequences`. When these sequences are passed to
/// [`DictionaryBuilder::build_from_sequences`](super::DictionaryBuilder::build_from_sequences),
/// the resulting sshash string ID for each sequence equals its index in
/// `sequences`, establishing the segment_id → string_id mapping.
#[derive(Debug, Clone)]
pub struct CfSegData {
    /// Cuttlefish segment identifiers, one per unitig, in file order.
    pub segment_ids: Vec<u64>,
    /// DNA sequences (unitigs), one per line, in file order.
    pub sequences: Vec<String>,
}

impl CfSegData {
    /// Number of segments.
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// Returns `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
}

/// Parse a `.cf_seg` file into segment IDs and sequences.
///
/// # Arguments
/// * `path` - Path to the `.cf_seg` file
///
/// # Errors
/// Returns an error if the file cannot be opened, a line is malformed
/// (missing ID or sequence), or a segment ID cannot be parsed as `u64`.
pub fn parse_cf_seg<P: AsRef<Path>>(path: P) -> Result<CfSegData> {
    let path = path.as_ref();
    let file = std::fs::File::open(path)
        .with_context(|| format!("Failed to open cf_seg file: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut segment_ids = Vec::new();
    let mut sequences = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line
            .with_context(|| format!("Failed to read line {} of {}", line_num + 1, path.display()))?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let (id_str, seq) = line.split_once(|c: char| c.is_ascii_whitespace())
            .with_context(|| {
                format!(
                    "Malformed cf_seg line {} in {}: expected '<id> <sequence>'",
                    line_num + 1,
                    path.display()
                )
            })?;

        let seg_id: u64 = id_str.parse()
            .with_context(|| {
                format!(
                    "Invalid segment ID '{}' on line {} of {}",
                    id_str,
                    line_num + 1,
                    path.display()
                )
            })?;

        let seq = seq.trim();
        if seq.is_empty() {
            bail!(
                "Empty sequence on line {} of {}",
                line_num + 1,
                path.display()
            );
        }

        super::parse::validate_dna_sequence(seq.as_bytes())
            .with_context(|| {
                format!(
                    "Invalid DNA in segment {} on line {} of {}",
                    seg_id,
                    line_num + 1,
                    path.display()
                )
            })?;

        segment_ids.push(seg_id);
        sequences.push(seq.to_uppercase());
    }

    Ok(CfSegData {
        segment_ids,
        sequences,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parse_cf_seg_basic() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "100\tACGTACGT")?;
        writeln!(f, "42\tTGCATGCA")?;
        writeln!(f, "999\tAAAACCCC")?;
        f.flush()?;

        let data = parse_cf_seg(f.path())?;
        assert_eq!(data.len(), 3);
        assert_eq!(data.segment_ids, vec![100, 42, 999]);
        assert_eq!(data.sequences, vec!["ACGTACGT", "TGCATGCA", "AAAACCCC"]);
        Ok(())
    }

    #[test]
    fn test_parse_cf_seg_spaces() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "1 ACGT")?;
        writeln!(f, "2\tTGCA")?;
        f.flush()?;

        let data = parse_cf_seg(f.path())?;
        assert_eq!(data.len(), 2);
        assert_eq!(data.segment_ids, vec![1, 2]);
        Ok(())
    }

    #[test]
    fn test_parse_cf_seg_skips_blank_lines() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "1\tACGT")?;
        writeln!(f)?;
        writeln!(f, "2\tTGCA")?;
        f.flush()?;

        let data = parse_cf_seg(f.path())?;
        assert_eq!(data.len(), 2);
        Ok(())
    }

    #[test]
    fn test_parse_cf_seg_lowercase() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "1\tacgtacgt")?;
        f.flush()?;

        let data = parse_cf_seg(f.path())?;
        assert_eq!(data.sequences[0], "ACGTACGT");
        Ok(())
    }

    #[test]
    fn test_parse_cf_seg_malformed_no_seq() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "42").unwrap();
        f.flush().unwrap();

        assert!(parse_cf_seg(f.path()).is_err());
    }

    #[test]
    fn test_parse_cf_seg_bad_id() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "notanumber\tACGT").unwrap();
        f.flush().unwrap();

        assert!(parse_cf_seg(f.path()).is_err());
    }

    #[test]
    fn test_parse_cf_seg_invalid_dna() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "1\tACGTNACGT").unwrap();
        f.flush().unwrap();

        assert!(parse_cf_seg(f.path()).is_err());
    }
}
