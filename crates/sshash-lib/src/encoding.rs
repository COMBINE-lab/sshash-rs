//! DNA nucleotide encoding
//!
//! This module implements the 2-bit encoding scheme for DNA nucleotides,
//! matching the C++ SSHash implementation exactly.
//!
//! Default encoding (SSHash custom):
//! - A (65/97)  -> 00
//! - C (67/99)  -> 01
//! - G (71/103) -> 11
//! - T (84/116) -> 10

use thiserror::Error;

/// Error type for encoding operations
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum EncodingError {
    /// The input byte is not a valid DNA base (A/C/G/T)
    #[error("Invalid DNA base: {0:?}")]
    InvalidBase(u8),
    /// The input string is not a valid k-mer
    #[error("Invalid k-mer string: {0}")]
    InvalidKmer(String),
    /// The input string length does not match the expected k-mer length
    #[error("K-mer length mismatch: expected {expected}, got {actual}")]
    LengthMismatch {
        /// Expected k-mer length
        expected: usize,
        /// Actual string length
        actual: usize,
    },
}

/// Encode a single DNA nucleotide to 2 bits
///
/// Uses the SSHash custom encoding by default:
/// A -> 00, C -> 01, G -> 11, T -> 10
#[inline]
pub const fn encode_base(base: u8) -> Result<u8, EncodingError> {
    match base {
        b'A' | b'a' => Ok(0b00),
        b'C' | b'c' => Ok(0b01),
        b'G' | b'g' => Ok(0b11),
        b'T' | b't' => Ok(0b10),
        _ => Err(EncodingError::InvalidBase(base)),
    }
}

/// Decode a 2-bit value to DNA nucleotide (uppercase)
#[inline]
pub const fn decode_base(bits: u8) -> u8 {
    match bits & 0b11 {
        0b00 => b'A',
        0b01 => b'C',
        0b11 => b'G',
        0b10 => b'T',
        _ => unreachable!(),
    }
}

/// Get the complement of a DNA base (encoded)
#[inline]
pub const fn complement_base(bits: u8) -> u8 {
    // For our encoding: A(00) <-> T(10), C(01) <-> G(11)
    // XOR with 0b10 gives the complement
    bits ^ 0b10
}

/// Encode a DNA string to a bit-packed representation
///
/// # Errors
/// Returns an error if the string contains invalid bases
pub fn encode_string(s: &str) -> Result<Vec<u8>, EncodingError> {
    encode_sequence(s.as_bytes())
}

/// Encode a DNA sequence (byte slice) to a bit-packed representation
///
/// # Errors
/// Returns an error if the sequence contains invalid bases
pub fn encode_sequence(sequence: &[u8]) -> Result<Vec<u8>, EncodingError> {
    let mut result = Vec::with_capacity(sequence.len().div_ceil(4));
    let mut current = 0u8;
    let mut bit_pos = 0;

    for (i, &base) in sequence.iter().enumerate() {
        let encoded = encode_base(base).map_err(|_| {
            EncodingError::InvalidKmer(format!("Invalid base at position {}: {:?}", i, base as char))
        })?;

        current |= encoded << bit_pos;
        bit_pos += 2;

        if bit_pos == 8 {
            result.push(current);
            current = 0;
            bit_pos = 0;
        }
    }

    if bit_pos > 0 {
        result.push(current);
    }

    Ok(result)
}

/// Decode a bit-packed representation back to a DNA string
pub fn decode_string(data: &[u8], length: usize) -> String {
    let mut result = String::with_capacity(length);
    let mut bit_pos = 0;
    let mut byte_idx = 0;

    for _ in 0..length {
        if byte_idx >= data.len() {
            break;
        }

        let bits = (data[byte_idx] >> bit_pos) & 0b11;
        result.push(decode_base(bits) as char);

        bit_pos += 2;
        if bit_pos == 8 {
            bit_pos = 0;
            byte_idx += 1;
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A').unwrap(), 0b00);
        assert_eq!(encode_base(b'a').unwrap(), 0b00);
        assert_eq!(encode_base(b'C').unwrap(), 0b01);
        assert_eq!(encode_base(b'c').unwrap(), 0b01);
        assert_eq!(encode_base(b'G').unwrap(), 0b11);
        assert_eq!(encode_base(b'g').unwrap(), 0b11);
        assert_eq!(encode_base(b'T').unwrap(), 0b10);
        assert_eq!(encode_base(b't').unwrap(), 0b10);

        // Invalid bases
        assert!(encode_base(b'N').is_err());
        assert!(encode_base(b'X').is_err());
        assert!(encode_base(b'0').is_err());
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(0b00), b'A');
        assert_eq!(decode_base(0b01), b'C');
        assert_eq!(decode_base(0b11), b'G');
        assert_eq!(decode_base(0b10), b'T');
    }

    #[test]
    fn test_complement_base() {
        assert_eq!(complement_base(0b00), 0b10); // A -> T
        assert_eq!(complement_base(0b10), 0b00); // T -> A
        assert_eq!(complement_base(0b01), 0b11); // C -> G
        assert_eq!(complement_base(0b11), 0b01); // G -> C
    }

    #[test]
    fn test_encode_decode_roundtrip() {
        let sequences = vec!["ACGT", "AAAA", "TTTT", "ACGTACGT", "GATTACA"];

        for seq in sequences {
            let encoded = encode_string(seq).unwrap();
            let decoded = decode_string(&encoded, seq.len());
            assert_eq!(decoded, seq.to_uppercase());
        }
    }

    #[test]
    fn test_encode_mixed_case() {
        let lower = encode_string("acgt").unwrap();
        let upper = encode_string("ACGT").unwrap();
        assert_eq!(lower, upper);
    }

    #[test]
    fn test_encode_invalid() {
        assert!(encode_string("ACGTN").is_err());
        assert!(encode_string("ACGT X").is_err());
    }
}
