//! String encoding module for building the Spectrum-Preserving String Set
//!
//! Encodes DNA sequences into the 2-bit representation, extracts k-mers,
//! and builds the offsets structure for string boundaries.

use crate::encoding;
use crate::kmer::{Kmer, KmerBits};
use crate::offsets::OffsetsVector;
use crate::spectrum_preserving_string_set::SpectrumPreservingStringSet;
use anyhow::Result;

/// Encoder for building SPSS from DNA sequences
///
/// Accumulates sequences, encodes them to 2-bit format, tracks offsets,
/// and builds the final SPSS structure.
///
/// IMPORTANT: Bases are packed contiguously across sequences without
/// byte-boundary padding. This means `base_idx / 4` always gives the correct
/// byte and `(base_idx % 4) * 2` gives the correct bit offset within that byte.
pub struct Encoder<const K: usize>
where
    Kmer<K>: KmerBits,
{
    /// Encoded strings (2-bit packed, contiguous across all sequences)
    strings: Vec<u8>,
    
    /// Offset to start of each string (in bases from the beginning of all strings)
    offsets: OffsetsVector,
    
    /// Total number of k-mers
    num_kmers: u64,
    
    /// Total number of strings added
    num_strings: u64,
    
    /// Total number of bases added so far (tracks bit position for contiguous packing)
    total_bases: u64,
}

impl<const K: usize> Encoder<K>
where
    Kmer<K>: KmerBits,
{
    /// Create a new encoder
    pub fn new() -> Self {
        Self {
            strings: Vec::new(),
            offsets: OffsetsVector::new(),  // Already starts with [0]
            num_kmers: 0,
            num_strings: 0,
            total_bases: 0,
        }
    }
    
    /// Add a DNA sequence to the encoder
    ///
    /// Bases are packed contiguously into the byte buffer without padding
    /// between sequences, so decode_kmer can use simple base_idx/4 arithmetic.
    ///
    /// # Arguments
    /// * `sequence` - DNA sequence (A, C, G, T only)
    ///
    /// # Errors
    /// Returns error if sequence contains invalid characters or is too short
    pub fn add_sequence(&mut self, sequence: &[u8]) -> Result<()> {
        let seq_len = sequence.len();
        
        // Skip sequences too short to contain a k-mer
        if seq_len < K {
            return Ok(());
        }
        
        // Pack each base contiguously into the byte buffer
        for (i, &base) in sequence.iter().enumerate() {
            let encoded = encoding::encode_base(base).map_err(|_| {
                anyhow::anyhow!("Invalid base at position {}: {:?}", i, base as char)
            })?;
            
            let base_idx = self.total_bases as usize;
            let byte_idx = base_idx / 4;
            let bit_offset = (base_idx % 4) * 2;
            
            // Extend buffer if needed
            if byte_idx >= self.strings.len() {
                self.strings.push(0);
            }
            
            self.strings[byte_idx] |= encoded << bit_offset;
            self.total_bases += 1;
        }
        
        // Update offsets (offset is in bases, not bytes)
        self.offsets.push(self.total_bases);
        
        // Count k-mers in this string
        let kmers_in_string = if seq_len >= K {
            (seq_len - K + 1) as u64
        } else {
            0
        };
        
        self.num_kmers += kmers_in_string;
        self.num_strings += 1;
        
        Ok(())
    }
    
    /// Get the current number of k-mers
    pub fn num_kmers(&self) -> u64 {
        self.num_kmers
    }
    
    /// Get the current number of strings
    pub fn num_strings(&self) -> u64 {
        self.num_strings
    }
    
    /// Build the final SpectrumPreservingStringSet
    ///
    /// Consumes the encoder and returns the SPSS.
    pub fn build(self, m: usize) -> SpectrumPreservingStringSet {
        SpectrumPreservingStringSet::from_parts(
            self.strings,
            self.offsets,
            K,
            m,
        )
    }
}

impl<const K: usize> Default for Encoder<K>
where
    Kmer<K>: KmerBits,
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_encoder_creation() {
        let encoder = Encoder::<31>::new();
        assert_eq!(encoder.num_kmers(), 0);
        assert_eq!(encoder.num_strings(), 0);
    }
    
    #[test]
    fn test_encoder_add_sequence() {
        let mut encoder = Encoder::<7>::new();
        
        // Add sequence "ACGTACGT" (length 8, contains 2 k=7-mers)
        encoder.add_sequence(b"ACGTACGT").unwrap();
        
        assert_eq!(encoder.num_strings(), 1);
        assert_eq!(encoder.num_kmers(), 2);  // 8 - 7 + 1 = 2
    }
    
    #[test]
    fn test_encoder_skip_short_sequence() {
        let mut encoder = Encoder::<31>::new();
        
        // Add sequence shorter than k
        encoder.add_sequence(b"ACGT").unwrap();  // Length 4 < 31
        
        assert_eq!(encoder.num_strings(), 0);  // Not counted
        assert_eq!(encoder.num_kmers(), 0);
    }
    
    #[test]
    fn test_encoder_multiple_sequences() {
        let mut encoder = Encoder::<5>::new();
        
        encoder.add_sequence(b"ACGTACGT").unwrap();  // 8 bases, 4 k=5-mers
        encoder.add_sequence(b"TGCA").unwrap();      // 4 bases < 5, skipped
        encoder.add_sequence(b"AAAAAAA").unwrap();   // 7 bases, 3 k=5-mers
        
        assert_eq!(encoder.num_strings(), 2);  // Only 2 sequences >= k
        assert_eq!(encoder.num_kmers(), 7);    // 4 + 3 = 7
    }
    
    #[test]
    fn test_encoder_build_spss() {
        let mut encoder = Encoder::<7>::new();
        encoder.add_sequence(b"ACGTACGT").unwrap();
        encoder.add_sequence(b"TGCATGCA").unwrap();
        
        let spss = encoder.build(5);  // m=5
        
        assert_eq!(spss.num_strings(), 2);
        assert_eq!(spss.total_bases(), 16);  // 8 + 8
    }
}
