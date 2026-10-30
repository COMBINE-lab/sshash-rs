//! DArray: constant-time Select on bitvectors.
//!
//! Port of the C++ `darray` from Okanohara & Sadakane (2007). Splits the
//! bitvector into blocks of L=1024 set bits. Dense blocks (span < 2^16) use
//! a base position + 16-bit subblock offsets + popcount scan. Sparse blocks
//! store all positions directly.
//!
//! Space: ≤1.5625 bits per indexed bit (dense); 64 bits per bit (sparse).

const BLOCK_SIZE: usize = 1024;
const SUBBLOCK_SIZE: usize = 32;
const MAX_IN_BLOCK_DISTANCE: u64 = 1 << 16;

pub(crate) struct DArray {
    num_positions: u64,
    block_inventory: Vec<i64>,
    subblock_inventory: Vec<u16>,
    overflow_positions: Vec<u64>,
}

impl DArray {
    pub fn build(bitvec: &[u64], num_bits: usize, negate: bool) -> Self {
        let mut num_positions: u64 = 0;
        let mut cur_block: Vec<u64> = Vec::with_capacity(BLOCK_SIZE);
        let mut block_inventory: Vec<i64> = Vec::new();
        let mut subblock_inventory: Vec<u16> = Vec::new();
        let mut overflow_positions: Vec<u64> = Vec::new();

        for (word_idx, &raw) in bitvec.iter().enumerate() {
            let mut word = if negate { !raw } else { raw };
            let base = (word_idx as u64) << 6;

            while word != 0 {
                let lsb = word.trailing_zeros() as u64;
                let pos = base + lsb;
                if pos >= num_bits as u64 {
                    break;
                }

                cur_block.push(pos);
                num_positions += 1;

                if cur_block.len() == BLOCK_SIZE {
                    Self::flush_block(
                        &mut cur_block,
                        &mut block_inventory,
                        &mut subblock_inventory,
                        &mut overflow_positions,
                    );
                }

                word &= word - 1;
            }
        }

        if !cur_block.is_empty() {
            Self::flush_block(
                &mut cur_block,
                &mut block_inventory,
                &mut subblock_inventory,
                &mut overflow_positions,
            );
        }

        Self {
            num_positions,
            block_inventory,
            subblock_inventory,
            overflow_positions,
        }
    }

    fn flush_block(
        cur_block: &mut Vec<u64>,
        block_inventory: &mut Vec<i64>,
        subblock_inventory: &mut Vec<u16>,
        overflow_positions: &mut Vec<u64>,
    ) {
        let span = cur_block.last().unwrap() - cur_block.first().unwrap();
        if span < MAX_IN_BLOCK_DISTANCE {
            let base = cur_block[0];
            block_inventory.push(base as i64);
            for i in (0..cur_block.len()).step_by(SUBBLOCK_SIZE) {
                subblock_inventory.push((cur_block[i] - base) as u16);
            }
        } else {
            block_inventory.push(-(overflow_positions.len() as i64) - 1);
            overflow_positions.extend_from_slice(cur_block);
            for _ in (0..cur_block.len()).step_by(SUBBLOCK_SIZE) {
                subblock_inventory.push(u16::MAX);
            }
        }
        cur_block.clear();
    }

    #[inline]
    pub fn select(&self, bitvec: &[u64], i: u64, negate: bool) -> u64 {
        debug_assert!(i < self.num_positions);
        let block = (i / BLOCK_SIZE as u64) as usize;
        let block_pos = self.block_inventory[block];

        if block_pos < 0 {
            let overflow_pos = (-block_pos - 1) as usize;
            return self.overflow_positions[overflow_pos + (i & (BLOCK_SIZE as u64 - 1)) as usize];
        }

        let subblock = (i / SUBBLOCK_SIZE as u64) as usize;
        let start_pos = block_pos as u64 + self.subblock_inventory[subblock] as u64;
        let mut remainder = i & (SUBBLOCK_SIZE as u64 - 1);
        if remainder == 0 {
            return start_pos;
        }

        let mut word_idx = (start_pos >> 6) as usize;
        let word_shift = start_pos & 63;
        let raw = bitvec[word_idx];
        let mut word = if negate { !raw } else { raw };
        word &= u64::MAX << word_shift;

        loop {
            let popcnt = word.count_ones() as u64;
            if remainder < popcnt {
                break;
            }
            remainder -= popcnt;
            word_idx += 1;
            let raw = bitvec[word_idx];
            word = if negate { !raw } else { raw };
        }

        (word_idx as u64) << 6 | select_in_word(word, remainder)
    }

    #[inline]
    #[allow(dead_code)]
    pub fn num_positions(&self) -> u64 {
        self.num_positions
    }

    pub fn heap_bytes(&self) -> usize {
        self.block_inventory.len() * std::mem::size_of::<i64>()
            + self.subblock_inventory.len() * std::mem::size_of::<u16>()
            + self.overflow_positions.len() * std::mem::size_of::<u64>()
    }

    pub fn write_to<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        write_u64(writer, self.num_positions)?;
        write_vec_i64(writer, &self.block_inventory)?;
        write_vec_u16(writer, &self.subblock_inventory)?;
        write_vec_u64(writer, &self.overflow_positions)?;
        Ok(())
    }

    pub fn read_from<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        let num_positions = read_u64(reader)?;
        let block_inventory = read_vec_i64(reader)?;
        let subblock_inventory = read_vec_u16(reader)?;
        let overflow_positions = read_vec_u64(reader)?;
        Ok(Self { num_positions, block_inventory, subblock_inventory, overflow_positions })
    }
}

// ---------------------------------------------------------------------------
// Little-endian serialization helpers
// ---------------------------------------------------------------------------

pub(crate) fn write_u64<W: std::io::Write>(w: &mut W, v: u64) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

pub(crate) fn write_u32<W: std::io::Write>(w: &mut W, v: u32) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

pub(crate) fn read_u64<R: std::io::Read>(r: &mut R) -> std::io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

pub(crate) fn read_u32<R: std::io::Read>(r: &mut R) -> std::io::Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

pub(crate) fn write_vec_u64<W: std::io::Write>(w: &mut W, v: &[u64]) -> std::io::Result<()> {
    write_u64(w, v.len() as u64)?;
    for &x in v { w.write_all(&x.to_le_bytes())?; }
    Ok(())
}

pub(crate) fn read_vec_u64<R: std::io::Read>(r: &mut R) -> std::io::Result<Vec<u64>> {
    let len = read_u64(r)? as usize;
    let mut v = Vec::with_capacity(len);
    let mut buf = [0u8; 8];
    for _ in 0..len {
        r.read_exact(&mut buf)?;
        v.push(u64::from_le_bytes(buf));
    }
    Ok(v)
}

fn write_vec_i64<W: std::io::Write>(w: &mut W, v: &[i64]) -> std::io::Result<()> {
    write_u64(w, v.len() as u64)?;
    for &x in v { w.write_all(&x.to_le_bytes())?; }
    Ok(())
}

fn read_vec_i64<R: std::io::Read>(r: &mut R) -> std::io::Result<Vec<i64>> {
    let len = read_u64(r)? as usize;
    let mut v = Vec::with_capacity(len);
    let mut buf = [0u8; 8];
    for _ in 0..len {
        r.read_exact(&mut buf)?;
        v.push(i64::from_le_bytes(buf));
    }
    Ok(v)
}

fn write_vec_u16<W: std::io::Write>(w: &mut W, v: &[u16]) -> std::io::Result<()> {
    write_u64(w, v.len() as u64)?;
    for &x in v { w.write_all(&x.to_le_bytes())?; }
    Ok(())
}

fn read_vec_u16<R: std::io::Read>(r: &mut R) -> std::io::Result<Vec<u16>> {
    let len = read_u64(r)? as usize;
    let mut v = Vec::with_capacity(len);
    let mut buf = [0u8; 2];
    for _ in 0..len {
        r.read_exact(&mut buf)?;
        v.push(u16::from_le_bytes(buf));
    }
    Ok(v)
}

/// Find the position of the k-th set bit (0-indexed) in a 64-bit word.
///
/// Dispatches to hardware-accelerated implementations when available:
/// - x86_64 with BMI2: uses PDEP (single instruction)
/// - aarch64 and x86_64 without BMI2: byte-wise broadword selection
/// - Fallback: bit-clearing loop
#[inline(always)]
fn select_in_word(word: u64, k: u64) -> u64 {
    // Compile-time path: when target-cpu=native enables BMI2 at compile time,
    // this path is taken without any runtime check.
    #[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))]
    {
        unsafe {
            use core::arch::x86_64::_pdep_u64;
            let deposited = _pdep_u64(1u64 << k, word);
            return deposited.trailing_zeros() as u64;
        }
    }

    // Runtime-detected path: when compiled without -C target-cpu=native but
    // running on BMI2 hardware.
    #[cfg(all(target_arch = "x86_64", not(target_feature = "bmi2")))]
    {
        if std::is_x86_feature_detected!("bmi2") {
            return unsafe { select_in_word_pdep(word, k) };
        }
        select_in_word_broadword(word, k)
    }

    // aarch64: no PDEP equivalent; use the byte-wise broadword selection
    #[cfg(target_arch = "aarch64")]
    {
        return select_in_word_broadword(word, k);
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        select_in_word_broadword(word, k)
    }
}

#[cfg(all(target_arch = "x86_64", not(target_feature = "bmi2")))]
#[target_feature(enable = "bmi2")]
#[inline]
unsafe fn select_in_word_pdep(word: u64, k: u64) -> u64 {
    use core::arch::x86_64::_pdep_u64;
    let deposited = _pdep_u64(1u64 << k, word);
    deposited.trailing_zeros() as u64
}

/// Byte-wise in-word select (Vigna's broadword select with a byte table).
///
/// Per-byte popcounts (SWAR) are turned into per-byte prefix sums by one
/// multiplication; the byte holding the `k`-th one is found with a single
/// subtract-and-mask over all eight bytes, and a 2 KiB table finishes inside
/// that byte. Unlike a halving search this has no chain of dependent
/// popcounts: measured ~43% lower latency than the previous six-step version
/// on Apple M5 and ~39% on x86-64 without BMI2.
#[inline(always)]
#[allow(dead_code)]
fn select_in_word_broadword(word: u64, k: u64) -> u64 {
    const L8: u64 = 0x0101_0101_0101_0101;
    const H8: u64 = 0x8080_8080_8080_8080;
    debug_assert!(k < word.count_ones() as u64);
    let mut s = word - ((word >> 1) & 0x5555_5555_5555_5555);
    s = (s & 0x3333_3333_3333_3333) + ((s >> 2) & 0x3333_3333_3333_3333);
    s = (s + (s >> 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    // Byte i of `prefix` = number of ones in bytes 0..=i (at most 64, so no carry).
    let prefix = s.wrapping_mul(L8);
    // Byte i has its high bit set iff prefix_i <= k; counting them gives the
    // index of the byte that contains the k-th one.
    let le = ((k.wrapping_mul(L8) | H8) - prefix) & H8;
    let bit = le.count_ones() as u64 * 8;
    let ones_before = ((prefix << 8) >> bit) & 0xFF;
    let byte = (word >> bit) & 0xFF;
    bit + SELECT_IN_BYTE[(byte | ((k - ones_before) << 8)) as usize] as u64
}

/// `SELECT_IN_BYTE[b | r << 8]` = position of the `r`-th one of byte `b`.
static SELECT_IN_BYTE: [u8; 2048] = {
    let mut t = [0u8; 2048];
    let mut b = 0;
    while b < 256 {
        let (mut r, mut pos) = (0, 0);
        while pos < 8 {
            if (b >> pos) & 1 == 1 {
                t[b | (r << 8)] = pos as u8;
                r += 1;
            }
            pos += 1;
        }
        b += 1;
    }
    t
};

#[cfg(test)]
mod tests {
    use super::*;

    fn make_bitvec(positions: &[u64]) -> (Vec<u64>, usize) {
        if positions.is_empty() {
            return (vec![0], 0);
        }
        let max_pos = *positions.last().unwrap();
        let num_words = (max_pos / 64 + 2) as usize;
        let mut bitvec = vec![0u64; num_words];
        for &p in positions {
            bitvec[(p / 64) as usize] |= 1u64 << (p % 64);
        }
        (bitvec, (max_pos + 1) as usize)
    }

    #[test]
    fn select_in_word_matches_naive() {
        let naive = |w: u64, k: u64| {
            let mut seen = 0;
            (0..64).find(|&p| (w >> p) & 1 == 1 && { seen += 1; seen == k + 1 }).unwrap()
        };
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        let mut words = vec![u64::MAX, 1, 1 << 63, 0x8000_0000_0000_0001, 0xFF00_0000_0000_00FF];
        for _ in 0..20_000 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            words.extend([x, x & (x >> 3) & (x >> 7), x | (x << 5) | (x >> 11)]);
        }
        for w in words.into_iter().filter(|&w| w != 0) {
            for k in 0..w.count_ones() as u64 {
                assert_eq!(select_in_word_broadword(w, k), naive(w, k), "w={w:#x} k={k}");
                assert_eq!(select_in_word(w, k), naive(w, k), "w={w:#x} k={k}");
            }
        }
    }

    #[test]
    fn test_select1_small() {
        let ones = vec![0u64, 2, 4, 6, 9];
        let (bitvec, num_bits) = make_bitvec(&ones);
        let da = DArray::build(&bitvec, num_bits, false);
        assert_eq!(da.num_positions(), 5);
        for (i, &expected) in ones.iter().enumerate() {
            assert_eq!(da.select(&bitvec, i as u64, false), expected);
        }
    }

    #[test]
    fn test_select0_small() {
        let ones = vec![0u64, 2, 4, 6, 9];
        let num_bits = 10;
        let num_words = 2;
        let mut bitvec = vec![0u64; num_words];
        for &p in &ones {
            bitvec[(p / 64) as usize] |= 1u64 << (p % 64);
        }
        let zeros: Vec<u64> = (0..10u64).filter(|p| !ones.contains(p)).collect();
        let da = DArray::build(&bitvec, num_bits, true);
        assert_eq!(da.num_positions(), zeros.len() as u64);
        for (i, &expected) in zeros.iter().enumerate() {
            assert_eq!(
                da.select(&bitvec, i as u64, true),
                expected,
                "select0({i})"
            );
        }
    }

    #[test]
    fn test_select1_multi_block() {
        let mut positions: Vec<u64> = (0..2048).map(|i| i * 3 + (i % 7)).collect();
        positions.sort_unstable();
        positions.dedup();
        let (bitvec, num_bits) = make_bitvec(&positions);
        let da = DArray::build(&bitvec, num_bits, false);
        assert_eq!(da.num_positions() as usize, positions.len());
        for (i, &expected) in positions.iter().enumerate() {
            assert_eq!(
                da.select(&bitvec, i as u64, false),
                expected,
                "select1({i})"
            );
        }
    }

    #[test]
    fn test_select_sparse_block() {
        // Create positions with huge gaps to trigger sparse blocks
        let mut positions: Vec<u64> = Vec::new();
        for i in 0..1024u64 {
            positions.push(i * 100);
        }
        let (bitvec, num_bits) = make_bitvec(&positions);
        let da = DArray::build(&bitvec, num_bits, false);
        assert_eq!(da.num_positions(), 1024);
        for (i, &expected) in positions.iter().enumerate() {
            assert_eq!(
                da.select(&bitvec, i as u64, false),
                expected,
                "select1({i})"
            );
        }
    }
}
