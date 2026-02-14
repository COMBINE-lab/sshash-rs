use sshash_lib::Dictionary;
use sshash_lib::serialization;
use std::io::{self, Cursor, Write, Seek, SeekFrom, BufReader, BufWriter, Read};
use epserde::prelude::*;
use tracing::info;

fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let path = std::env::args().nth(1).unwrap_or_else(|| "ecoli1.ssi".to_string());
    info!("Loading dictionary from: {}", path);

    let dict = Dictionary::load(&path)?;
    info!("Loaded. k={}, m={}, canonical={}", dict.k(), dict.m(), dict.canonical());
    info!("Num minimizers: {}", dict.num_minimizers());

    // Access internals via public fields/methods
    // We need to measure sizes by serializing each component separately

    // 1. Serialize just the SPSS
    let mut spss_buf = Vec::new();
    {
        let spss = dict.spss();
        spss.serialize(&mut spss_buf)
            .map_err(|e| anyhow::anyhow!("SPSS serialize error: {}", e))?;
    }

    // 2. Serialize control map without MPHF
    let mut cm_buf = Vec::new();
    {
        let cm = dict.control_map_ref();
        cm.serialize_without_mphf(&mut cm_buf)?;
    }

    // 3. Serialize index without MPHF
    let mut idx_buf = Vec::new();
    {
        let idx = dict.index_ref();
        idx.serialize_without_mphf(&mut idx_buf)?;
    }

    // 4. Measure MPHF
    let mphf_path = serialization::mphf_container_path(&path);
    let mphf_size = std::fs::metadata(&mphf_path)?.len();

    // 5. Measure sub-components of index
    let idx = dict.index_ref();

    // EF sequence
    let ef_bytes = idx.ef_bytes();

    // offsets BitFieldVec
    let offsets_bytes = idx.offsets_bytes();

    // skew_index (without MPHF)
    let mut skew_buf = Vec::new();
    idx.skew_index.serialize_without_mphf(&mut skew_buf)?;

    // 6. Measure MPHF partitions individually
    let mut mphf_reader = BufReader::new(std::fs::File::open(&mphf_path)?);
    let mphf_container = serialization::read_mphf_container(&mut mphf_reader)?;

    // Re-read container to get partition entries
    let mut mphf_file = std::fs::File::open(&mphf_path)?;
    mphf_file.read(&mut [0u8; 8])?; // skip magic
    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];
    mphf_file.read_exact(&mut buf4)?; // ver major
    mphf_file.read_exact(&mut buf4)?; // ver minor
    mphf_file.read_exact(&mut buf4)?; // num partitions
    let num_parts = u32::from_le_bytes(buf4);

    let mut partition_sizes = Vec::new();
    for _ in 0..num_parts {
        mphf_file.read_exact(&mut buf4)?; // partition_id
        let pid = u32::from_le_bytes(buf4);
        mphf_file.read_exact(&mut buf8)?; // byte_offset
        mphf_file.read_exact(&mut buf8)?; // byte_size
        let bsz = u64::from_le_bytes(buf8);
        partition_sizes.push((pid, bsz));
    }

    // 7. SPSS sub-components
    let spss = dict.spss();
    let spss_strings_bytes = spss.strings_len();
    let spss_num_strings = spss.num_strings();
    let spss_total_bases = spss.total_bases();
    let spss_offsets_entries = spss_num_strings + 1;

    // Main index file size
    let index_path = serialization::index_file_path(&path);
    let index_file_size = std::fs::metadata(&index_path)?.len();

    let num_kmers = 4_877_400u64;
    let num_minimizers = dict.num_minimizers();

    // Print results
    println!("================================================================================");
    println!("RUST SSHash INDEX SIZE BREAKDOWN - ecoli1 (k={}, m={})", dict.k(), dict.m());
    println!("================================================================================");
    println!();
    println!("{:<45} {:>12} {:>8} {:>10}", "Component", "Bytes", "MB", "bits/kmer");
    println!("{:<45} {:>12} {:>8} {:>10}", "-".repeat(45), "-".repeat(12), "-".repeat(8), "-".repeat(10));

    let bpk = |bytes: u64| -> f64 { bytes as f64 * 8.0 / num_kmers as f64 };
    let mb = |bytes: u64| -> f64 { bytes as f64 / 1024.0 / 1024.0 };

    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "Header", 37, mb(37), bpk(37));
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "SPSS (total, epserde)", spss_buf.len(), mb(spss_buf.len() as u64), bpk(spss_buf.len() as u64));
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  strings Vec<u8>", spss_strings_bytes, mb(spss_strings_bytes), bpk(spss_strings_bytes));
    println!("  -> {} strings, {} total bases", spss_num_strings, spss_total_bases);
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  offsets Vec<u64> (est)", spss_offsets_entries * 8, mb(spss_offsets_entries * 8), bpk(spss_offsets_entries * 8));
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  k,m + epserde overhead", spss_buf.len() as u64 - spss_strings_bytes - spss_offsets_entries * 8, mb(spss_buf.len() as u64 - spss_strings_bytes - spss_offsets_entries * 8), bpk(spss_buf.len() as u64 - spss_strings_bytes - spss_offsets_entries * 8));
    println!();

    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "ControlMap (without MPHF)", cm_buf.len(), mb(cm_buf.len() as u64), bpk(cm_buf.len() as u64));
    println!("  -> {} controls x 16 bytes each", num_minimizers);
    let controls_array_bytes = num_minimizers * 16;
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  controls[] array only", controls_array_bytes, mb(controls_array_bytes), bpk(controls_array_bytes));
    println!();

    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "SparseAndSkewIndex (without MPHF)", idx_buf.len(), mb(idx_buf.len() as u64), bpk(idx_buf.len() as u64));
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  num_super_kmers_before_bucket (EF)", ef_bytes, mb(ef_bytes as u64), bpk(ef_bytes as u64));
    println!("  -> {} buckets", idx.num_buckets());
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  offsets (BitFieldVec)", offsets_bytes, mb(offsets_bytes as u64), bpk(offsets_bytes as u64));
    println!("  -> {} entries x {} bits", idx.num_offsets(), idx.offsets.bit_width());
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  skew_index (without MPHF)", skew_buf.len(), mb(skew_buf.len() as u64), bpk(skew_buf.len() as u64));
    println!("  -> {} partitions", idx.skew_index.num_partitions());
    println!();

    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "MPHF container (total)", mphf_size, mb(mphf_size), bpk(mphf_size));
    for (pid, bsz) in &partition_sizes {
        let label = if *pid == 0 { "control map MPHF" } else { "skew index MPHF" };
        println!("{:<45} {:>12} {:>8.3} {:>10.2}", format!("  partition {} ({})", pid, label), *bsz, mb(*bsz), bpk(*bsz));
        if *bsz > 0 {
            let keys = if *pid == 0 { num_minimizers } else { 0 }; // estimate
            if keys > 0 {
                println!("  -> {:.2} bits/key for {} keys", *bsz as f64 * 8.0 / keys as f64, keys);
            }
        }
    }
    println!();

    let total = index_file_size + mphf_size;
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "TOTAL (main + MPHF)", total, mb(total), bpk(total));

    Ok(())
}
