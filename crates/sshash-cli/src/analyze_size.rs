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
    
    // control_codewords
    let mut cc_buf = Vec::new();
    idx.control_codewords.serialize(&mut cc_buf)
        .map_err(|e| anyhow::anyhow!("CC serialize error: {}", e))?;
    
    // mid_load_buckets
    let mut mlb_buf = Vec::new();
    idx.mid_load_buckets.serialize(&mut mlb_buf)
        .map_err(|e| anyhow::anyhow!("MLB serialize error: {}", e))?;
    
    // begin_buckets_of_size
    let bbs_size = 4 + idx.begin_buckets_of_size.len() * 4;
    
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
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  control_codewords (epserde)", cc_buf.len(), mb(cc_buf.len() as u64), bpk(cc_buf.len() as u64));
    println!("  -> {} entries x {} bits", idx.control_codewords.len(), idx.control_codewords.bit_width());
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  begin_buckets_of_size", bbs_size, mb(bbs_size as u64), bpk(bbs_size as u64));
    println!("{:<45} {:>12} {:>8.3} {:>10.2}", "  mid_load_buckets (epserde)", mlb_buf.len(), mb(mlb_buf.len() as u64), bpk(mlb_buf.len() as u64));
    println!("  -> {} entries x {} bits", idx.mid_load_buckets.len(), idx.mid_load_buckets.bit_width());
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
    println!();
    
    // C++ comparison
    println!("================================================================================");
    println!("COMPARISON WITH C++ (3.6 MB total, 5.9 bits/kmer)");
    println!("================================================================================");
    println!();
    
    let cpp_total: u64 = 3_600_000;
    let cpp_mphf = (0.40 * num_kmers as f64 / 8.0) as u64;
    let cpp_strings_offsets = (0.06 * num_kmers as f64 / 8.0) as u64;
    let cpp_control = (3.42 * num_kmers as f64 / 8.0) as u64;
    let cpp_mid_load = (0.009 * num_kmers as f64 / 8.0) as u64;
    let cpp_strings = (2.01 * num_kmers as f64 / 8.0) as u64;
    
    println!("{:<35} {:>12} {:>12} {:>12} {:>8}", "Component", "C++ bytes", "Rust bytes", "Excess", "Factor");
    println!("{:<35} {:>12} {:>12} {:>12} {:>8}", "-".repeat(35), "-".repeat(12), "-".repeat(12), "-".repeat(12), "-".repeat(8));
    
    let cmp = |name: &str, cpp: u64, rust: u64| {
        let excess = rust as i64 - cpp as i64;
        let factor = if cpp > 0 { rust as f64 / cpp as f64 } else { f64::INFINITY };
        println!("{:<35} {:>12} {:>12} {:>+12} {:>7.1}x", name, cpp, rust, excess, factor);
    };
    
    cmp("MPHF", cpp_mphf, mphf_size);
    cmp("strings_offsets", cpp_strings_offsets, spss_offsets_entries * 8);
    cmp("control_codewords", cpp_control, cc_buf.len() as u64);
    cmp("strings (2-bit packed)", cpp_strings, spss_strings_bytes);
    cmp("mid_load_buckets", cpp_mid_load, mlb_buf.len() as u64);
    println!("{:<35} {:>12} {:>12} {:>+12} {:>8}", "controls[] (NOT IN C++)", "0", format!("{}", controls_array_bytes), format!("+{}", controls_array_bytes), "N/A");
    println!();
    cmp("TOTAL", cpp_total, total);
    
    // Root cause analysis
    println!();
    println!("================================================================================");
    println!("ROOT CAUSE ANALYSIS");
    println!("================================================================================");
    println!();
    
    let excess = total as i64 - cpp_total as i64;
    println!("Total excess: {} bytes ({:.2} MB)", excess, excess as f64 / 1024.0 / 1024.0);
    println!();
    
    let c1 = controls_array_bytes as i64;
    let c2 = mphf_size as i64 - cpp_mphf as i64;
    let c3 = (spss_offsets_entries * 8) as i64 - cpp_strings_offsets as i64;
    let c4 = spss_strings_bytes as i64 - cpp_strings as i64;
    
    let mut causes: Vec<(&str, i64, &str)> = vec![
        ("controls[] Vec<MinimizerControl>", c1, 
         "C++ stores no per-minimizer metadata. bucket_id IS the MPHF output."),
        ("MPHF (fmph vs PTHash)", c2,
         "ph::fmph uses ~7 bits/key vs PTHash's ~2.8 bits/key."),
        ("Offsets (Vec<u64> vs Elias-Fano)", c3,
         "Rust uses 64 bits per offset. C++ uses Elias-Fano (~0.06 bits/kmer)."),
        ("Strings data difference", c4,
         "Encoding difference between implementations."),
    ];
    
    causes.sort_by(|a, b| b.1.cmp(&a.1));
    
    let mut accounted = 0i64;
    for (name, size, explanation) in &causes {
        let pct = *size as f64 / excess as f64 * 100.0;
        println!("  {}: {} bytes ({:.1}%)", name, size, pct);
        println!("    {}", explanation);
        println!();
        accounted += size;
    }
    
    let remaining = excess - accounted;
    println!("  Other (epserde headers, metadata, padding): {} bytes ({:.1}%)", 
             remaining, remaining as f64 / excess as f64 * 100.0);
    
    Ok(())
}
