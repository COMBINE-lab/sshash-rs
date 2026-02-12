use clap::{Parser, Subcommand};
use sshash_lib::{BuildConfiguration, Dictionary, DictionaryBuilder, Kmer, KmerBits};
use sshash_lib::streaming_query::{LookupResult, StreamingQueryEngine};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use tracing::{info, debug, warn};

#[derive(Parser)]
#[command(name = "sshash")]
#[command(version = "0.1.0")]
#[command(about = "SSHash: Sparse and Skew Hashing of k-mers", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build a dictionary from input file
    Build {
        /// Input FASTA/FASTQ file
        #[arg(short, long)]
        input: String,

        /// K-mer length
        #[arg(short, long)]
        k: usize,

        /// Minimizer length  
        #[arg(short, long)]
        m: usize,

        /// Output file
        #[arg(short, long)]
        output: Option<String>,

        /// Use canonical k-mers (k-mer or reverse complement, whichever is smaller)
        #[arg(long, default_value = "false")]
        canonical: bool,

        /// Number of threads (0 = all available cores)
        #[arg(short = 't', long, default_value = "0")]
        threads: usize,

        /// RAM limit in GiB for external sorting (0 = unlimited/in-memory)
        #[arg(short = 'r', long, default_value = "8")]
        ram_limit: usize,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },

    /// Query a dictionary
    Query {
        /// Index file
        #[arg(short, long)]
        index: String,

        /// Query file (FASTA/FASTQ or plain text with one k-mer per line)
        #[arg(short, long)]
        query: String,

        /// Use streaming query (FASTA/FASTQ only)
        #[arg(long, default_value = "false")]
        streaming: bool,
    },

    /// Check correctness of a dictionary
    Check {
        /// Index file
        #[arg(short, long)]
        index: String,

        /// Input file used to build the dictionary
        #[arg(short = 'f', long)]
        input: String,

        /// Use streaming query for validation (FASTA/FASTQ only)
        #[arg(long, default_value = "false")]
        streaming: bool,
    },

    /// Run performance benchmarks
    Bench {
        /// Index file
        #[arg(short, long)]
        index: String,
    },
}

fn main() -> anyhow::Result<()> {
    // Initialize tracing: use RUST_LOG if set, otherwise default to info
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Build { input, k, m, output, canonical, threads, ram_limit, verbose } => {
            build_command(input, k, m, output, canonical, threads, ram_limit, verbose)?;
        }
        Commands::Query { index, query, streaming } => {
            query_command(index, query, streaming)?;
        }
        Commands::Check { index, input, streaming } => {
            check_command(index, input, streaming)?;
        }
        Commands::Bench { index } => {
            bench_command(index)?;
        }
    }

    Ok(())
}

/// Build a dictionary from FASTA/FASTQ input
#[allow(clippy::too_many_arguments)]
fn build_command(
    input: String,
    k: usize,
    m: usize,
    output: Option<String>,
    canonical: bool,
    threads: usize,
    ram_limit: usize,
    verbose: bool,
) -> anyhow::Result<()> {
    info!("Building SSHash dictionary...");
    info!("  Input: {}", input);
    info!("  k: {}", k);
    info!("  m: {}", m);
    info!("  Canonical: {}", canonical);
    info!("  RAM limit: {} GiB", ram_limit);
    
    // Read sequences - detect format by trying FASTA/FASTQ first, then plain text
    let sequences = parse_sequences_file(&input)?;
    info!("  Loaded {} sequences", sequences.len());
    
    // Create builder configuration
    let mut config = BuildConfiguration::new(k, m)
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    config.canonical = canonical;
    config.verbose = verbose;
    config.ram_limit_gib = ram_limit;
    config.num_threads = threads;
    
    // Build dictionary
    let builder = DictionaryBuilder::new(config)
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    let dict = builder.build_from_sequences(sequences)
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    
    // Determine output filename
    let output_path = output.unwrap_or_else(|| {
        let input_path = Path::new(&input);
        let stem = input_path.file_stem().unwrap().to_str().unwrap();
        format!("{}.sshash", stem)
    });
    
    // Save dictionary
    info!("Saving dictionary to {}...", output_path);
    dict.save(&output_path)?;
    
    info!("Dictionary built successfully!");
    
    // Print space breakdown
    dict.print_space_breakdown();
    
    Ok(())
}

/// Run performance benchmarks
fn bench_command(index: String) -> anyhow::Result<()> {
    let index = normalize_index_path(&index);
    info!("Loading dictionary from {}...", index);
    
    let dict = Dictionary::load(&index)?;
    let k = dict.k();
    info!("Dictionary loaded (k={}, m={}, canonical={})", k, dict.m(), dict.canonical());
    
    dict.print_space_breakdown();

    sshash_lib::dispatch_on_k!(k, K => {
        bench_with_k::<K>(&dict)
    })
}

fn bench_with_k<const K: usize>(dict: &Dictionary) -> anyhow::Result<()>
where
    Kmer<K>: KmerBits,
{
    use std::time::Instant;
    use std::hint::black_box;
    
    let k = dict.k();
    let num_queries: usize = 1_000_000;
    let runs: usize = 5;
    
    // Generate positive lookup queries by walking all strings in SPSS
    info!("Generating {} positive lookup queries...", num_queries);
    let mut positive_queries: Vec<Kmer<K>> = Vec::with_capacity(num_queries);
    {
        let num_strings = dict.num_strings();
        let mut rng_state: u64 = 42;
        for _ in 0..num_queries {
            // Simple LCG random
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let string_id = (rng_state >> 33) % num_strings;
            let string_len = dict.string_length(string_id);
            if string_len < k {
                continue;
            }
            let pos = ((rng_state >> 16) as usize) % (string_len - k + 1);
            let kmer = dict.access_kmer::<K>(string_id, pos);
            positive_queries.push(kmer);
        }
    }
    info!("  Generated {} queries", positive_queries.len());
    
    // Benchmark positive lookup
    {
        let start = Instant::now();
        for _ in 0..runs {
            for kmer in &positive_queries {
                let res = dict.lookup(kmer);
                black_box(res);
            }
        }
        let elapsed = start.elapsed();
        let ns_per_kmer = elapsed.as_nanos() as f64 / (runs * positive_queries.len()) as f64;
        println!("positive lookup (avg_nanosec_per_kmer) = {:.3}", ns_per_kmer);
    }
    
    // Generate negative lookup queries (random k-mers)
    info!("Generating {} negative lookup queries...", num_queries);
    let mut negative_queries: Vec<Kmer<K>> = Vec::with_capacity(num_queries);
    {
        let mut rng_state: u64 = 12345;
        for _ in 0..num_queries {
            let mut bits: u128 = 0;
            for _ in 0..k {
                rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let base = (rng_state >> 33) & 3;
                bits = (bits << 2) | base as u128;
            }
            negative_queries.push(Kmer::<K>::from_bits(bits));
        }
    }
    
    // Benchmark negative lookup
    {
        let start = Instant::now();
        for _ in 0..runs {
            for kmer in &negative_queries {
                let res = dict.lookup(kmer);
                black_box(res);
            }
        }
        let elapsed = start.elapsed();
        let ns_per_kmer = elapsed.as_nanos() as f64 / (runs * negative_queries.len()) as f64;
        println!("negative lookup (avg_nanosec_per_kmer) = {:.3}", ns_per_kmer);
    }
    
    // Benchmark EF (offsets locate)
    // This isolates the partition_point search to see if S+ tree would help
    {
        info!("Benchmarking offsets locate (Elias-Fano)...");
        let num_strings = dict.num_strings();
        // Total positions is approximately total bases in SPSS
        // We use num_strings * avg_string_length as an approximation
        let mut total_positions: u64 = 0;
        for sid in 0..num_strings {
            total_positions += dict.string_length(sid) as u64;
        }
        if total_positions == 0 {
            total_positions = 1; // Avoid div-by-zero
        }
        
        // Generate random positions to look up
        let num_locate_queries = 1_000_000usize;
        let mut positions: Vec<u64> = Vec::with_capacity(num_locate_queries);
        let mut rng_state: u64 = 999;
        for _ in 0..num_locate_queries {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let pos = (rng_state >> 16) % total_positions;
            positions.push(pos);
        }
        
        let start = Instant::now();
        for _ in 0..runs {
            for &pos in &positions {
                let res = dict.locate_string(pos);
                black_box(res);
            }
        }
        let elapsed = start.elapsed();
        let ns_per_locate = elapsed.as_nanos() as f64 / (runs * positions.len()) as f64;
        println!("locate (Elias-Fano) (avg_nanosec_per_query) = {:.3}", ns_per_locate);

        println!("  offsets array size = {} strings", num_strings);
        println!("  log2(size) = {:.1} comparisons expected", (num_strings as f64).log2());
    }
    
    Ok(())
}

/// Query k-mers from a dictionary
fn query_command(index: String, query: String, streaming: bool) -> anyhow::Result<()> {
    let index = normalize_index_path(&index);
    info!("Loading dictionary from {}...", index);
    
    let dict = Dictionary::load(&index)?;
    let k = dict.k();
    info!("Dictionary loaded (k={}, m={}, canonical={})", k, dict.m(), dict.canonical());

    sshash_lib::dispatch_on_k!(k, K => {
        query_with_k::<K>(&dict, &query, streaming)
    })
}

fn query_with_k<const K: usize>(dict: &Dictionary, query: &str, streaming: bool) -> anyhow::Result<()>
where
    Kmer<K>: KmerBits,
{
    let k = dict.k();

    if streaming {
        return query_with_k_streaming::<K>(dict, query, k);
    }
    
    // Read query k-mers
    let queries = if query.ends_with(".fa") || query.ends_with(".fasta") || query.ends_with(".fq") || query.ends_with(".fastq")
        || query.ends_with(".fa.gz") || query.ends_with(".fasta.gz") || query.ends_with(".fq.gz") || query.ends_with(".fastq.gz")
        || query.ends_with(".spss") {
        // Parse sequence file (FASTA/FASTQ or plain text) and extract k-mers
        extract_kmers_from_sequences(query, k)?
    } else {
        // Parse as plain text, one k-mer per line
        parse_kmer_file(query)?
    };
    
    info!("Querying {} k-mers...", queries.len());
    
    let mut found = 0;
    let mut not_found = 0;
    
    for (i, kmer_str) in queries.iter().enumerate() {
        if kmer_str.len() != k {
            continue;
        }
        
        // Parse k-mer
        if let Ok(kmer) = Kmer::<K>::from_string(kmer_str) {
            let pos = dict.lookup(&kmer);
            if pos != sshash_lib::constants::INVALID_UINT64 {
                found += 1;
                if i < 10 {
                    debug!("  {} -> position {}", kmer_str, pos);
                }
            } else {
                not_found += 1;
                debug!("  {} -> NOT FOUND (index {})", kmer_str, i);
            }
        }
    }
    
    if queries.len() > 10 {
        debug!("  ... (showing first 10 results)");
    }
    
    println!("\nResults:");
    println!("  Found: {}", found);
    println!("  Not found: {}", not_found);
    println!("  Total: {}", queries.len());
    println!("  Hit rate: {:.2}%", (found as f64 / queries.len() as f64) * 100.0);
    
    Ok(())
}

fn query_with_k_streaming<const K: usize>(
    dict: &Dictionary,
    query: &str,
    k: usize,
) -> anyhow::Result<()>
where
    Kmer<K>: KmerBits,
{
    // Strip .gz suffix for extension check
    let name = query.strip_suffix(".gz").unwrap_or(query);
    if !(name.ends_with(".fa")
        || name.ends_with(".fasta")
        || name.ends_with(".fq")
        || name.ends_with(".fastq"))
    {
        return Err(anyhow::anyhow!("Streaming query requires FASTA/FASTQ input"));
    }

    let sequences = parse_fasta_file(query)?;
    let mut engine: StreamingQueryEngine<K> = StreamingQueryEngine::new(dict);

    let mut found = 0u64;
    let mut not_found = 0u64;
    let mut shown = 0usize;

    for seq in sequences {
        if seq.len() < k {
            continue;
        }
        engine.reset();
        for i in 0..=(seq.len() - k) {
            let kmer_str = &seq[i..i + k];
            let result = engine.lookup(kmer_str);
            if result.is_found() {
                found += 1;
                if shown < 10 {
                    debug!("  {} -> position {}", kmer_str, result.kmer_offset);
                    shown += 1;
                }
            } else {
                not_found += 1;
                if shown < 10 {
                    debug!("  {} -> NOT FOUND", kmer_str);
                    shown += 1;
                }
            }
        }
    }

    let total = found + not_found;
    if total > 10 {
        debug!("  ... (showing first 10 results)");
    }

    println!("\nResults:");
    println!("  Found: {}", found);
    println!("  Not found: {}", not_found);
    println!("  Total: {}", total);
    if total > 0 {
        println!("  Hit rate: {:.2}%", (found as f64 / total as f64) * 100.0);
    }

    let stats = engine.stats();
    println!("\nStreaming stats:");
    println!("  Searches: {}", stats.num_searches);
    println!("  Extensions: {}", stats.num_extensions);
    println!("  Invalid: {}", stats.num_invalid);
    println!("  Negative: {}", stats.num_negative);

    Ok(())
}

/// Check dictionary correctness by re-querying all k-mers from input
fn check_command(index: String, input: String, streaming: bool) -> anyhow::Result<()> {
    let index = normalize_index_path(&index);
    info!("Checking dictionary correctness...");
    info!("  Index: {}", index);
    info!("  Input: {}", input);

    let dict = Dictionary::load(&index)?;
    let k = dict.k();
    info!("Dictionary loaded (k={}, m={}, canonical={})", k, dict.m(), dict.canonical());
    
    sshash_lib::dispatch_on_k!(k, K => {
        check_with_k::<K>(&dict, &input, streaming)
    })
}

fn check_with_k<const K: usize>(dict: &Dictionary, input: &str, streaming: bool) -> anyhow::Result<()>
where
    Kmer<K>: KmerBits,
{
    let k = dict.k();

    if streaming {
        return check_with_k_streaming::<K>(dict, input, k);
    }
    
    // Extract all k-mers from input
    let kmers = extract_kmers_from_sequences(input, k)?;
    info!("Extracted {} k-mers from input", kmers.len());
    
    // Query all k-mers
    let mut found = 0;
    let mut not_found = 0;
    let mut errors = Vec::new();
    
    for (i, kmer_str) in kmers.iter().enumerate() {
        if let Ok(kmer) = Kmer::<K>::from_string(kmer_str) {
            let pos = dict.lookup(&kmer);
            if pos != sshash_lib::constants::INVALID_UINT64 {
                found += 1;
            } else {
                not_found += 1;
                if errors.is_empty() {
                    let mini = dict.debug_extract_minimizer(&kmer);
                    let has_mphf = dict.debug_has_control_map_mphf();
                    let in_map = dict.debug_control_map_lookup(mini.value);
                    let bucket_info = dict.debug_bucket_info(&kmer);
                    debug!("Debug (first miss):");
                    debug!("  minimizer value: {}", mini.value);
                    debug!("  minimizer pos_in_kmer: {}", mini.pos_in_kmer);
                    debug!("  control map MPHF present: {}", has_mphf);
                    debug!("  control map contains minimizer: {}", in_map);
                    if let Some((minimizer, bucket_id, control_code)) = bucket_info {
                        debug!("  bucket_id: {}", bucket_id);
                        debug!("  control_code: {}", control_code);
                        debug!("  control_code LSB: {}", control_code & 0b11);
                        debug!("  minimizer (from bucket): {}", minimizer);
                    }
                }
                if errors.len() < 10 {
                    errors.push(kmer_str.clone());
                }
            }
        }
        
        if (i + 1) % 100000 == 0 {
            info!("  Checked {} k-mers...", i + 1);
        }
    }
    
    println!("\n=== Check Results ===");
    println!("  Total k-mers: {}", kmers.len());
    println!("  Found: {}", found);
    println!("  Not found: {}", not_found);
    
    if not_found > 0 {
        warn!("CORRECTNESS CHECK FAILED! {} k-mers were not found", not_found);
        println!("\n✗ CORRECTNESS CHECK FAILED!");
        println!("  {} k-mers were not found in the dictionary", not_found);
        if !errors.is_empty() {
            println!("\nFirst {} missing k-mers:", errors.len());
            for kmer in &errors {
                println!("  {}", kmer);
            }
        }
        std::process::exit(1);
    } else {
        println!("\n✓ CORRECTNESS CHECK PASSED!");
        println!("  All k-mers from input were found in the dictionary");
    }
    
    Ok(())
}

fn normalize_index_path(index: &str) -> String {
    if let Some(stripped) = index.strip_suffix(".ssi.mphf") {
        stripped.to_string()
    } else if let Some(stripped) = index.strip_suffix(".ssi") {
        stripped.to_string()
    } else {
        index.to_string()
    }
}

fn check_with_k_streaming<const K: usize>(
    dict: &Dictionary,
    input: &str,
    k: usize,
) -> anyhow::Result<()>
where
    Kmer<K>: KmerBits,
{
    if !(input.ends_with(".fa")
        || input.ends_with(".fasta")
        || input.ends_with(".fq")
        || input.ends_with(".fastq")
        || input.ends_with(".fa.gz")
        || input.ends_with(".fasta.gz")
        || input.ends_with(".fq.gz")
        || input.ends_with(".fastq.gz"))
    {
        return Err(anyhow::anyhow!("Streaming check requires FASTA/FASTQ input"));
    }

    let sequences = parse_fasta_file(input)?;
    let mut engine: StreamingQueryEngine<K> = StreamingQueryEngine::new(dict);

    let mut total = 0u64;
    let mut mismatches = 0u64;
    let mut errors: Vec<String> = Vec::new();

    for seq in sequences {
        if seq.len() < k {
            continue;
        }
        engine.reset();
        for i in 0..=(seq.len() - k) {
            let kmer_str = &seq[i..i + k];
            total += 1;

            // Direct lookup: use query() to get a full LookupResult
            let direct_result = if let Ok(kmer) = Kmer::<K>::from_string(kmer_str) {
                dict.query(&kmer)
            } else {
                LookupResult::not_found()
            };

            let streaming_result = engine.lookup(kmer_str);

            // Compare fields matching C++ equal_lookup_result:
            // kmer_id, kmer_id_in_string, kmer_orientation, string_id, string_begin, string_end
            // (NOT kmer_offset)
            let equal = direct_result.kmer_id == streaming_result.kmer_id
                && direct_result.kmer_id_in_string == streaming_result.kmer_id_in_string
                && direct_result.kmer_orientation == streaming_result.kmer_orientation
                && direct_result.string_id == streaming_result.string_id
                && direct_result.string_begin == streaming_result.string_begin
                && direct_result.string_end == streaming_result.string_end;

            if !equal {
                mismatches += 1;
                if errors.len() < 10 {
                    errors.push(format!(
                        "{} (direct: id={},sid={},in_str={} | streaming: id={},sid={},in_str={})",
                        kmer_str,
                        direct_result.kmer_id, direct_result.string_id,
                        direct_result.kmer_id_in_string,
                        streaming_result.kmer_id, streaming_result.string_id,
                        streaming_result.kmer_id_in_string,
                    ));
                }
            }

            if total % 100000 == 0 {
                info!("  Checked {} k-mers...", total);
            }
        }
    }

    println!("\n=== Streaming Check Results ===");
    println!("  Total k-mers: {}", total);
    println!("  Mismatches: {}", mismatches);

    if mismatches > 0 {
        warn!("STREAMING CHECK FAILED! {} mismatches found", mismatches);
        println!("\n✗ STREAMING CHECK FAILED!");
        if !errors.is_empty() {
            println!("\nFirst {} mismatches:", errors.len());
            for err in &errors {
                println!("  {}", err);
            }
        }
        std::process::exit(1);
    } else {
        println!("\n✓ STREAMING CHECK PASSED!");
        println!("  Streaming results match direct lookup");
    }

    Ok(())
}

/// Parse FASTA/FASTQ file and return sequences
fn parse_fasta_file(path: &str) -> anyhow::Result<Vec<String>> {
    use needletail::{parse_fastx_file};
    
    let mut sequences = Vec::new();
    let mut reader = parse_fastx_file(path)?;
    
    while let Some(record) = reader.next() {
        let record = record?;
        let seq_bytes = record.seq();
        let seq = std::str::from_utf8(&seq_bytes)?;
        sequences.push(seq.to_uppercase());
    }
    
    Ok(sequences)
}

/// Parse plain text file with one sequence per line (SPSS format)
fn parse_plain_text_sequences(path: &str) -> anyhow::Result<Vec<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        let seq = line.trim().to_uppercase();
        if !seq.is_empty() {
            sequences.push(seq);
        }
    }
    
    Ok(sequences)
}

/// Parse sequences file - auto-detects FASTA/FASTQ vs plain text format
fn parse_sequences_file(path: &str) -> anyhow::Result<Vec<String>> {
    // Try FASTA/FASTQ first
    match parse_fasta_file(path) {
        Ok(seqs) => Ok(seqs),
        Err(_) => {
            // Fall back to plain text (one sequence per line)
            info!("  File is not FASTA/FASTQ, trying plain text format...");
            parse_plain_text_sequences(path)
        }
    }
}

/// Extract all k-mers from a sequence file (FASTA/FASTQ or plain text)
fn extract_kmers_from_sequences(path: &str, k: usize) -> anyhow::Result<Vec<String>> {
    let sequences = parse_sequences_file(path)?;
    let mut kmers = Vec::new();
    
    for seq in sequences {
        if seq.len() >= k {
            for i in 0..=(seq.len() - k) {
                kmers.push(seq[i..i+k].to_string());
            }
        }
    }
    
    Ok(kmers)
}

/// Parse a plain text file with one k-mer per line
fn parse_kmer_file(path: &str) -> anyhow::Result<Vec<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut kmers = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        let kmer = line.trim().to_uppercase();
        if !kmer.is_empty() {
            kmers.push(kmer);
        }
    }
    
    Ok(kmers)
}
