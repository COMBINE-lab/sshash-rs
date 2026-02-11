//! K-mer representation with const generics and optimal storage
//!
//! This module implements k-mer types using const generics to support
//! all odd k-mer sizes from 3 to 63. Storage is automatically selected
//! (u64 for K ≤ 31, u128 for K > 31) for optimal memory usage.

use crate::encoding::{decode_base, encode_base, EncodingError};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::ops::{BitAnd, BitOr, Not, Shl, Shr};

/// Trait defining optimal storage type for a given K
///
/// This trait is implemented for all valid odd K values from 3 to 63.
/// - K ≤ 31: uses u64 (8 bytes)
/// - K > 31: uses u128 (16 bytes)
pub trait KmerBits: Sized {
    /// The underlying storage type (u64 or u128)
    type Storage: Copy
        + Ord
        + Hash
        + fmt::Debug
        + fmt::Display
        + fmt::Binary
        + From<u8>
        + From<u64>
        + BitAnd<Output = Self::Storage>
        + BitOr<Output = Self::Storage>
        + Not<Output = Self::Storage>
        + Shl<usize, Output = Self::Storage>
        + Shr<usize, Output = Self::Storage>
        + Shr<i32, Output = Self::Storage>;

    /// Number of bits in the storage type
    const BITS: usize;

    /// Maximum k-mer length this storage can hold
    const MAX_K: usize = Self::BITS / 2 - 1;
    
    /// Convert storage to u8 (truncates)
    fn to_u8(val: Self::Storage) -> u8;
    
    /// Convert storage to u64
    fn to_u64(val: Self::Storage) -> u64;
    
    /// Convert storage to u128
    fn to_u128(val: Self::Storage) -> u128;
    
    /// Convert u8 to storage
    fn from_u8(val: u8) -> Self::Storage;
    
    /// Convert u64 to storage
    fn from_u64(val: u64) -> Self::Storage;
    
    /// Convert u128 to storage
    fn from_u128(val: u128) -> Self::Storage;
    
    /// Shift left on storage (works directly with native type)
    fn shl(val: Self::Storage, bits: usize) -> Self::Storage;
    
    /// Shift right on storage (works directly with native type)
    fn shr(val: Self::Storage, bits: usize) -> Self::Storage;
    
    /// Bitwise AND on storage
    fn bitand(a: Self::Storage, b: Self::Storage) -> Self::Storage;
    
    /// Bitwise OR on storage
    fn bitor(a: Self::Storage, b: Self::Storage) -> Self::Storage;
    
    /// Bitwise NOT on storage
    fn bitnot(a: Self::Storage) -> Self::Storage;
}

/// Implement KmerBits for all K ≤ 31 (use u64)
macro_rules! impl_kmer_bits_u64 {
    ($($k:literal),* $(,)?) => {
        $(
            impl KmerBits for Kmer<$k> {
                type Storage = u64;
                const BITS: usize = 64;
                
                #[inline]
                fn to_u8(val: Self::Storage) -> u8 {
                    val as u8
                }
                
                #[inline]
                fn to_u64(val: Self::Storage) -> u64 {
                    val
                }
                
                #[inline]
                fn to_u128(val: Self::Storage) -> u128 {
                    val as u128
                }
                
                #[inline]
                fn from_u8(val: u8) -> Self::Storage {
                    val as u64
                }
                
                #[inline]
                fn from_u64(val: u64) -> Self::Storage {
                    val
                }
                
                #[inline]
                fn from_u128(val: u128) -> Self::Storage {
                    val as u64
                }
                
                #[inline]
                fn shl(val: Self::Storage, bits: usize) -> Self::Storage {
                    val << bits
                }
                
                #[inline]
                fn shr(val: Self::Storage, bits: usize) -> Self::Storage {
                    val >> bits
                }
                
                #[inline]
                fn bitand(a: Self::Storage, b: Self::Storage) -> Self::Storage {
                    a & b
                }
                
                #[inline]
                fn bitor(a: Self::Storage, b: Self::Storage) -> Self::Storage {
                    a | b
                }
                
                #[inline]
                fn bitnot(a: Self::Storage) -> Self::Storage {
                    !a
                }
            }
        )*
    };
}

impl_kmer_bits_u64!(3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31);

/// Implement KmerBits for all K > 31 (use u128)
macro_rules! impl_kmer_bits_u128 {
    ($($k:literal),* $(,)?) => {
        $(
            impl KmerBits for Kmer<$k> {
                type Storage = u128;
                const BITS: usize = 128;
                
                #[inline]
                fn to_u8(val: Self::Storage) -> u8 {
                    val as u8
                }
                
                #[inline]
                fn to_u64(val: Self::Storage) -> u64 {
                    val as u64
                }
                
                #[inline]
                fn to_u128(val: Self::Storage) -> u128 {
                    val
                }
                
                #[inline]
                fn from_u8(val: u8) -> Self::Storage {
                    val as u128
                }
                
                #[inline]
                fn from_u64(val: u64) -> Self::Storage {
                    val as u128
                }
                
                #[inline]
                fn from_u128(val: u128) -> Self::Storage {
                    val
                }
                
                #[inline]
                fn shl(val: Self::Storage, bits: usize) -> Self::Storage {
                    val << bits
                }
                
                #[inline]
                fn shr(val: Self::Storage, bits: usize) -> Self::Storage {
                    val >> bits
                }
                
                #[inline]
                fn bitand(a: Self::Storage, b: Self::Storage) -> Self::Storage {
                    a & b
                }
                
                #[inline]
                fn bitor(a: Self::Storage, b: Self::Storage) -> Self::Storage {
                    a | b
                }
                
                #[inline]
                fn bitnot(a: Self::Storage) -> Self::Storage {
                    !a
                }
            }
        )*
    };
}

impl_kmer_bits_u128!(33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63);

/// K-mer representation with compile-time size and optimal storage
///
/// The storage type (u64 or u128) is automatically selected based on K:
/// - K ≤ 31: u64 (more memory efficient)
/// - K > 31: u128 (necessary for larger k-mers)
///
/// # Example
/// ```
/// use sshash_lib::kmer::Kmer;
///
/// // K=31 uses u64 internally
/// let kmer31: Kmer<31> = Kmer::from_str("ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
///
/// // K=63 uses u128 internally
/// let kmer63: Kmer<63> = Kmer::from_str(
///     "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"
/// ).unwrap();
/// ```
#[derive(Clone, Copy)]
pub struct Kmer<const K: usize>
where
    Kmer<K>: KmerBits,
{
    bits: <Kmer<K> as KmerBits>::Storage,
}

impl<const K: usize> Kmer<K>
where
    Kmer<K>: KmerBits,
{
    /// Create a new k-mer from raw bits
    #[inline]
    pub fn new(bits: <Kmer<K> as KmerBits>::Storage) -> Self {
        Self { bits }
    }

    /// Create a k-mer from a u128 value (truncated to storage size)
    #[inline]
    pub fn from_bits(bits: u128) -> Self {
        Self { bits: <Kmer<K> as KmerBits>::from_u128(bits) }
    }

    /// Get the raw bits
    #[inline]
    pub fn bits(&self) -> <Kmer<K> as KmerBits>::Storage {
        self.bits
    }

    /// Convert to u64 (panics if K > 31)
    #[inline]
    pub fn as_u64(&self) -> u64 {
        if K > 31 {
            panic!("Cannot convert Kmer<{}> to u64, use as_u128()", K);
        }
        unsafe { *(&self.bits as *const _ as *const u64) }
    }

    /// Convert to u128
    #[inline]
    pub fn as_u128(&self) -> u128 {
        if K <= 31 {
            unsafe { *(&self.bits as *const _ as *const u64) as u128 }
        } else {
            unsafe { *(&self.bits as *const _ as *const u128) }
        }
    }

    /// Create a k-mer from a DNA string
    ///
    /// This is an inherent method so callers don't need to import [`std::str::FromStr`].
    ///
    /// # Errors
    /// Returns an error if the string length doesn't match K or contains invalid bases.
    #[inline]
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Result<Self, EncodingError> {
        <Self as std::str::FromStr>::from_str(s)
    }

    /// Create a k-mer from a DNA string (alias for [`from_str`](Self::from_str))
    #[inline]
    pub fn from_string(s: &str) -> Result<Self, EncodingError> {
        Self::from_str(s)
    }

    /// Get the reverse complement of this k-mer
    ///
    /// Uses bit-parallel operations: complement via XOR, then reverse 2-bit pairs.
    #[inline]
    pub fn reverse_complement(&self) -> Self {
        if K <= 31 {
            // u64 fast path
            let mut x = <Kmer<K> as KmerBits>::to_u64(self.bits);
            // Complement: XOR with 0xAAAA... (10 repeating) flips A<->T, C<->G
            x ^= 0xAAAA_AAAA_AAAA_AAAAu64;
            // Reverse 2-bit pairs within u64 using byte-level operations
            // Step 1: Swap adjacent 2-bit pairs
            x = ((x >> 2) & 0x3333_3333_3333_3333u64) | ((x & 0x3333_3333_3333_3333u64) << 2);
            // Step 2: Swap adjacent 4-bit nibbles
            x = ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0Fu64) | ((x & 0x0F0F_0F0F_0F0F_0F0Fu64) << 4);
            // Step 3: Reverse bytes
            x = x.swap_bytes();
            // Shift right to align K bases (remove padding from MSB side)
            x >>= 64 - K * 2;
            Self { bits: <Kmer<K> as KmerBits>::from_u64(x) }
        } else {
            // u128 path for K > 31
            let mut x = <Kmer<K> as KmerBits>::to_u128(self.bits);
            x ^= 0xAAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAAu128;
            x = ((x >> 2) & 0x3333_3333_3333_3333_3333_3333_3333_3333u128)
              | ((x & 0x3333_3333_3333_3333_3333_3333_3333_3333u128) << 2);
            x = ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0Fu128)
              | ((x & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0Fu128) << 4);
            x = x.swap_bytes();
            x >>= 128 - K * 2;
            Self { bits: <Kmer<K> as KmerBits>::from_u128(x) }
        }
    }

    /// Get the canonical representation (minimum of forward and reverse complement)
    pub fn canonical(&self) -> Self {
        let rc = self.reverse_complement();
        if self.bits < rc.bits {
            *self
        } else {
            rc
        }
    }

    /// Extract a base at a specific position (0-indexed)
    pub fn get_base(&self, pos: usize) -> u8 {
        assert!(pos < K, "Position {} out of bounds for k-mer of length {}", pos, K);
        let shift = pos * 2;
        <Kmer<K> as KmerBits>::to_u8(
            <Kmer<K> as KmerBits>::bitand(
                <Kmer<K> as KmerBits>::shr(self.bits, shift),
                <Kmer<K> as KmerBits>::from_u8(0b11u8)
            )
        )
    }

    /// Set a base at a specific position (0-indexed)
    pub fn set_base(&mut self, pos: usize, base: u8) {
        assert!(pos < K, "Position {} out of bounds for k-mer of length {}", pos, K);
        assert!(base <= 0b11, "Base value must be 0-3");
        
        let shift = pos * 2;
        // Clear the bits at position, then set new value
        let mask = <Kmer<K> as KmerBits>::bitnot(
            <Kmer<K> as KmerBits>::shl(
                <Kmer<K> as KmerBits>::from_u8(0x3u8),
                shift
            )
        );
        self.bits = <Kmer<K> as KmerBits>::bitand(self.bits, mask);
        
        let new_bits = <Kmer<K> as KmerBits>::shl(
            <Kmer<K> as KmerBits>::from_u8(base),
            shift
        );
        self.bits = <Kmer<K> as KmerBits>::bitor(self.bits, new_bits);
    }
    
    /// Create an empty k-mer (all zeros)
    #[inline]
    pub fn empty() -> Self {
        Self {
            bits: <Kmer<K> as KmerBits>::from_u8(0),
        }
    }
    
    /// Append a base to the k-mer (shift left and add new base at position 0)
    ///
    /// This shifts the entire k-mer left by 2 bits and adds the new base
    /// at the lowest position. The highest base is lost.
    ///
    /// # Arguments
    /// * `base` - 2-bit encoded base (0=A, 1=C, 3=G, 2=T)
    #[inline]
    pub fn append_base(self, base: u8) -> Self {
        assert!(base <= 0b11, "Base value must be 0-3");
        
        // Shift left by 2 bits
        let shifted = <Kmer<K> as KmerBits>::shl(self.bits, 2);
        
        // Mask to keep only K bases (2K bits)
        let mask_bits = 2 * K;
        let mask = if mask_bits >= <Kmer<K> as KmerBits>::BITS {
            <Kmer<K> as KmerBits>::from_u8(0xFF) // All ones (won't happen for valid K)
        } else {
            <Kmer<K> as KmerBits>::from_u64((1u64 << mask_bits) - 1)
        };
        
        let masked = <Kmer<K> as KmerBits>::bitand(shifted, mask);
        
        // OR in the new base at position 0
        let new_bits = <Kmer<K> as KmerBits>::bitor(
            masked,
            <Kmer<K> as KmerBits>::from_u8(base)
        );
        
        Self { bits: new_bits }
    }
}

impl<const K: usize> PartialEq for Kmer<K>
where
    Kmer<K>: KmerBits,
{
    fn eq(&self, other: &Self) -> bool {
        self.bits == other.bits
    }
}

impl<const K: usize> Eq for Kmer<K> where Kmer<K>: KmerBits {}

impl<const K: usize> PartialOrd for Kmer<K>
where
    Kmer<K>: KmerBits,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const K: usize> Ord for Kmer<K>
where
    Kmer<K>: KmerBits,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.bits.cmp(&other.bits)
    }
}

impl<const K: usize> Hash for Kmer<K>
where
    Kmer<K>: KmerBits,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bits.hash(state);
    }
}

impl<const K: usize> fmt::Debug for Kmer<K>
where
    Kmer<K>: KmerBits,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Kmer<{}>(\"{}\")", K, self)
    }
}

impl<const K: usize> fmt::Display for Kmer<K>
where
    Kmer<K>: KmerBits,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut bits = self.bits;
        for _ in 0..K {
            let base_bits = <Kmer<K> as KmerBits>::to_u8(
                <Kmer<K> as KmerBits>::bitand(
                    bits,
                    <Kmer<K> as KmerBits>::from_u8(0b11u8)
                )
            );
            write!(f, "{}", decode_base(base_bits) as char)?;
            bits = <Kmer<K> as KmerBits>::shr(bits, 2);
        }
        Ok(())
    }
}

impl<const K: usize> Default for Kmer<K>
where
    Kmer<K>: KmerBits,
{
    fn default() -> Self {
        Self {
            bits: <Kmer<K> as KmerBits>::Storage::from(0u8),
        }
    }
}

impl<const K: usize> std::str::FromStr for Kmer<K>
where
    Kmer<K>: KmerBits,
{
    type Err = EncodingError;

    /// Create a k-mer from a DNA string
    ///
    /// # Errors
    /// Returns an error if:
    /// - The string length doesn't match K
    /// - The string contains invalid bases
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() != K {
            return Err(EncodingError::LengthMismatch {
                expected: K,
                actual: s.len(),
            });
        }

        let mut bits = <Kmer<K> as KmerBits>::from_u8(0);

        for (i, &base) in s.as_bytes().iter().enumerate() {
            let encoded = encode_base(base)?;
            let shifted = <Kmer<K> as KmerBits>::shl(
                <Kmer<K> as KmerBits>::from_u8(encoded),
                i * 2
            );
            bits = <Kmer<K> as KmerBits>::bitor(bits, shifted);
        }

        Ok(Self { bits })
    }
}

/// Type alias for common k-mer sizes
pub type Kmer31 = Kmer<31>;
/// Type alias for 21-mers
pub type Kmer21 = Kmer<21>;
/// Type alias for 63-mers
pub type Kmer63 = Kmer<63>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_storage_types() {
        // Verify K ≤ 31 uses u64
        assert_eq!(<Kmer<3> as KmerBits>::BITS, 64);
        assert_eq!(<Kmer<31> as KmerBits>::BITS, 64);
        assert_eq!(std::mem::size_of::<Kmer<31>>(), 8);

        // Verify K > 31 uses u128
        assert_eq!(<Kmer<33> as KmerBits>::BITS, 128);
        assert_eq!(<Kmer<63> as KmerBits>::BITS, 128);
        assert_eq!(std::mem::size_of::<Kmer<63>>(), 16);
    }

    #[test]
    fn test_kmer_from_str() {
        let kmer: Kmer<5> = Kmer::from_str("ACGTG").unwrap();
        assert_eq!(kmer.to_string(), "ACGTG");

        let kmer: Kmer<31> = Kmer::from_str("ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        assert_eq!(kmer.to_string(), "ACGTACGTACGTACGTACGTACGTACGTACG");
    }

    #[test]
    fn test_kmer_reverse_complement() {
        let kmer: Kmer<5> = Kmer::from_str("ACGTG").unwrap();
        let rc = kmer.reverse_complement();
        assert_eq!(rc.to_string(), "CACGT");

        // Another test with K=7
        let kmer: Kmer<7> = Kmer::from_str("ACGTACG").unwrap();
        let rc = kmer.reverse_complement();
        assert_eq!(rc.to_string(), "CGTACGT");
    }

    #[test]
    fn test_kmer_canonical() {
        let kmer: Kmer<5> = Kmer::from_str("ACGTG").unwrap();
        let canon = kmer.canonical();
        
        // The canonical should be the minimum of ACGTG and its RC (CACGT)
        let rc = kmer.reverse_complement();
        assert!(canon == kmer || canon == rc);
        assert!(canon.bits <= kmer.bits && canon.bits <= rc.bits);
    }

    #[test]
    fn test_kmer_case_insensitive() {
        let lower: Kmer<5> = Kmer::from_str("acgtg").unwrap();
        let upper: Kmer<5> = Kmer::from_str("ACGTG").unwrap();
        assert_eq!(lower, upper);
    }

    #[test]
    fn test_kmer_length_mismatch() {
        let result: Result<Kmer<5>, _> = Kmer::from_str("ACGT");
        assert!(result.is_err());
        
        let result: Result<Kmer<5>, _> = Kmer::from_str("ACGTGG");
        assert!(result.is_err());
    }

    #[test]
    fn test_kmer_get_set_base() {
        let mut kmer: Kmer<5> = Kmer::from_str("AAAAA").unwrap();
        assert_eq!(kmer.get_base(0), 0b00); // A
        
        kmer.set_base(2, 0b10); // Set middle to T
        assert_eq!(kmer.to_string(), "AATAA");
        
        kmer.set_base(4, 0b11); // Set last to G
        assert_eq!(kmer.to_string(), "AATAG");
    }

    #[test]
    fn test_kmer_ordering() {
        let kmer1: Kmer<5> = Kmer::from_str("AAAAA").unwrap();
        let kmer2: Kmer<5> = Kmer::from_str("AAAAC").unwrap();
        let kmer3: Kmer<5> = Kmer::from_str("TTTTT").unwrap();
        
        assert!(kmer1 < kmer2);
        assert!(kmer2 < kmer3);
        assert!(kmer1 < kmer3);
    }
}
