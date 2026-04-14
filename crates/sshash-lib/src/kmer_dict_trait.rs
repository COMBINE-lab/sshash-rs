//! Traits abstracting the k-mer dictionary interface.
//!
//! These exist so downstream callers (notably piscem-rs) can be generic over
//! the dictionary implementation. The existing [`Dictionary`] is the reference
//! implementation; a future small-reference variant will provide a second
//! implementation with a different internal representation but the same
//! observable behavior.

use crate::dictionary::Dictionary;
use crate::kmer::{Kmer, KmerBits};
use crate::streaming_query::{LookupResult, StreamingQueryEngine};

/// The streaming-query half of a k-mer dictionary.
///
/// Consecutive k-mer lookups share state via `reset` + `lookup`; concrete
/// implementations are free to exploit that (e.g. minimizer reuse, unitig
/// extension, rolling k-mer updates).
pub trait KmerStreamingQuery {
    /// Drop any streaming state so the next lookup parses its k-mer from scratch.
    fn reset(&mut self);

    /// Look up a single k-mer (ASCII bytes, length = `k`).
    fn lookup(&mut self, kmer_bytes: &[u8]) -> LookupResult;

    /// Number of lookups that required a full dictionary search (slow path).
    fn num_searches(&self) -> u64;

    /// Number of lookups resolved by extending along the current unitig (fast path).
    fn num_extensions(&self) -> u64;
}

/// A k-mer dictionary: maps canonical k-mers to their position in an SPSS.
///
/// This trait captures the call surface that the piscem-rs mapping pipeline
/// uses today. It is deliberately narrow — everything outside the mapping
/// hot path remains on the concrete types.
pub trait KmerDictionary {
    /// Streaming-query engine type. Parameterized by the same `K` the caller
    /// dispatches on; the lifetime borrows from `&self`.
    type Query<'a, const K: usize>: KmerStreamingQuery
    where
        Self: 'a,
        Kmer<K>: KmerBits;

    /// The k-mer length this dictionary was built for.
    fn k(&self) -> usize;

    /// The minimizer length this dictionary was built for.
    fn m(&self) -> usize;

    /// Number of SPSS strings (unitigs) in the dictionary.
    fn num_strings(&self) -> u64;

    /// Whether the dictionary was built in canonical mode.
    fn canonical(&self) -> bool;

    /// Decode the k-mer at the given absolute base position within the SPSS.
    ///
    /// Hot path in piscem's `HitSearcher`.
    fn kmer_at_pos<const K: usize>(&self, absolute_base_pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits;

    /// Construct a new streaming-query engine borrowing from this dictionary.
    fn create_streaming_query<const K: usize>(&self) -> Self::Query<'_, K>
    where
        Kmer<K>: KmerBits;
}

// ---------------------------------------------------------------------------
// Blanket impls for the existing sshash types.
// ---------------------------------------------------------------------------

impl<'a, const K: usize> KmerStreamingQuery for StreamingQueryEngine<'a, K>
where
    Kmer<K>: KmerBits,
{
    #[inline]
    fn reset(&mut self) {
        StreamingQueryEngine::reset(self)
    }

    #[inline]
    fn lookup(&mut self, kmer_bytes: &[u8]) -> LookupResult {
        StreamingQueryEngine::lookup(self, kmer_bytes)
    }

    #[inline]
    fn num_searches(&self) -> u64 {
        StreamingQueryEngine::num_searches(self)
    }

    #[inline]
    fn num_extensions(&self) -> u64 {
        StreamingQueryEngine::num_extensions(self)
    }
}

impl KmerDictionary for Dictionary {
    type Query<'a, const K: usize>
        = StreamingQueryEngine<'a, K>
    where
        Self: 'a,
        Kmer<K>: KmerBits;

    #[inline]
    fn k(&self) -> usize {
        Dictionary::k(self)
    }

    #[inline]
    fn m(&self) -> usize {
        Dictionary::m(self)
    }

    #[inline]
    fn num_strings(&self) -> u64 {
        Dictionary::num_strings(self)
    }

    #[inline]
    fn canonical(&self) -> bool {
        Dictionary::canonical(self)
    }

    #[inline]
    fn kmer_at_pos<const K: usize>(&self, absolute_base_pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        Dictionary::kmer_at_pos::<K>(self, absolute_base_pos)
    }

    #[inline]
    fn create_streaming_query<const K: usize>(&self) -> Self::Query<'_, K>
    where
        Kmer<K>: KmerBits,
    {
        Dictionary::create_streaming_query::<K>(self)
    }
}
