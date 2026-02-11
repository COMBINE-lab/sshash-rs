//! FASTA/FASTQ parsing with automatic decompression
//!
//! Reads DNA sequences from FASTA or FASTQ files, with transparent
//! gzip decompression. Validates DNA alphabet (A, C, G, T only).

use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use std::path::Path;

/// Parse a FASTA/FASTQ file and call a function for each valid DNA sequence
///
/// # Arguments
/// * `path` - Path to input file (may be gzipped)
/// * `callback` - Function called for each sequence, receives (name, sequence)
///
/// # Errors
/// Returns error if:
/// - File cannot be opened
/// - File format is invalid
/// - Sequence contains non-DNA characters
pub fn parse_sequences<P, F>(path: P, mut callback: F) -> Result<()>
where
    P: AsRef<Path>,
    F: FnMut(&[u8], &[u8]) -> Result<()>,
{
    let path = path.as_ref();
    
    // needletail automatically handles gzip decompression
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("Failed to open sequence file: {}", path.display()))?;
    
    while let Some(record) = reader.next() {
        let record = record
            .with_context(|| format!("Failed to parse sequence record in {}", path.display()))?;
        
        // Validate DNA alphabet (seq is Cow<[u8]>, convert to &[u8])
        let seq = record.seq();
        validate_dna_sequence(&seq)
            .with_context(|| format!("Invalid DNA sequence in {}", path.display()))?;
        
        // Call user callback with borrowed slices
        callback(record.id(), &seq)?;
    }
    
    Ok(())
}

/// Validate that a sequence contains only valid DNA bases (A, C, G, T)
///
/// # Arguments
/// * `seq` - Sequence bytes to validate
///
/// # Errors
/// Returns error if sequence contains non-ACGT characters
pub fn validate_dna_sequence(seq: &[u8]) -> Result<()> {
    for (i, &base) in seq.iter().enumerate() {
        match base {
            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' => {}
            _ => {
                return Err(anyhow::anyhow!(
                    "Invalid DNA base '{}' at position {}. Only A, C, G, T are allowed.",
                    base as char,
                    i
                ));
            }
        }
    }
    Ok(())
}

/// Count sequences and total bases in a file
///
/// # Arguments
/// * `path` - Path to input file
///
/// # Returns
/// `(num_sequences, total_bases)`
pub fn count_sequences<P: AsRef<Path>>(path: P) -> Result<(usize, usize)> {
    let mut num_sequences = 0;
    let mut total_bases = 0;
    
    parse_sequences(path, |_name, seq| {
        num_sequences += 1;
        total_bases += seq.len();
        Ok(())
    })?;
    
    Ok((num_sequences, total_bases))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;
    
    #[test]
    fn test_validate_dna_sequence_valid() {
        assert!(validate_dna_sequence(b"ACGT").is_ok());
        assert!(validate_dna_sequence(b"acgt").is_ok());
        assert!(validate_dna_sequence(b"ACGTacgt").is_ok());
    }
    
    #[test]
    fn test_validate_dna_sequence_invalid() {
        assert!(validate_dna_sequence(b"ACGTN").is_err());  // N
        assert!(validate_dna_sequence(b"ACGT ").is_err());   // Space
        assert!(validate_dna_sequence(b"ACG-T").is_err());   // Dash
    }
    
    #[test]
    fn test_parse_fasta_file() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, ">seq1")?;
        writeln!(temp_file, "ACGT")?;
        writeln!(temp_file, ">seq2")?;
        writeln!(temp_file, "TGCA")?;
        temp_file.flush()?;
        
        let mut sequences = Vec::new();
        parse_sequences(temp_file.path(), |name, seq| {
            sequences.push((name.to_vec(), seq.to_vec()));
            Ok(())
        })?;
        
        assert_eq!(sequences.len(), 2);
        assert_eq!(sequences[0].0, b"seq1");
        assert_eq!(sequences[0].1, b"ACGT");
        assert_eq!(sequences[1].0, b"seq2");
        assert_eq!(sequences[1].1, b"TGCA");
        
        Ok(())
    }
    
    #[test]
    fn test_count_sequences() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, ">seq1")?;
        writeln!(temp_file, "ACGT")?;
        writeln!(temp_file, ">seq2")?;
        writeln!(temp_file, "TGCATGCA")?;
        temp_file.flush()?;
        
        let (num_seqs, total_bases) = count_sequences(temp_file.path())?;
        assert_eq!(num_seqs, 2);
        assert_eq!(total_bases, 12);  // 4 + 8
        
        Ok(())
    }
}
