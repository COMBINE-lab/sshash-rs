//! Build configuration for SSHash dictionary construction
//!
//! Mirrors the C++ `build_configuration` struct with parameters for
//! dictionary building, minimizer computation, and resource limits.

use crate::constants::DEFAULT_SEED;
use std::path::PathBuf;

/// Configuration parameters for building an SSHash dictionary
#[derive(Debug, Clone)]
pub struct BuildConfiguration {
    /// K-mer length (must be odd, between 3 and 63)
    pub k: usize,
    
    /// Minimizer length (must be odd, m < k)
    pub m: usize,
    
    /// Seed for hash functions
    pub seed: u64,
    
    /// Number of threads for parallel operations (0 = all available cores)
    pub num_threads: usize,
    
    /// RAM limit in GiB for external sorting
    pub ram_limit_gib: usize,
    
    /// PTHash lambda parameter (trade-off for MPHF construction)
    /// Typically 3.5-4.0 for minimal size, higher for faster queries
    pub lambda: f64,
    
    /// Ignored since 0.7.0: the index is always canonical (the unified
    /// C++-v6-style minimizer scheme answers both strands at forward-mode
    /// cost, so there is no modality to choose).
    #[deprecated(since = "0.7.0", note = "the index is always canonical; this field is ignored")]
    pub canonical: bool,
    
    /// Use partitioned MPHF for parallel construction (default: true)
    pub partitioned_mphf: bool,

    /// Build weighted dictionary (with k-mer abundance/weights)
    pub weighted: bool,
    
    /// Verbose output during construction
    pub verbose: bool,

    /// Force the external-sort build path regardless of estimated size.
    /// For tests and diagnostics (the two paths must produce identical
    /// dictionaries); not exposed on the CLI.
    #[doc(hidden)]
    pub force_external_sort: bool,
    
    /// Directory for temporary files during construction
    pub tmp_dirname: PathBuf,
}

impl Default for BuildConfiguration {
    #[allow(deprecated)]
    fn default() -> Self {
        Self {
            k: 31,
            m: 19,  // Must be odd, less than k
            seed: DEFAULT_SEED,
            num_threads: 0, // 0 = use all available cores
            ram_limit_gib: 8,
            lambda: 6.0,  // C++ default
            canonical: true,
            partitioned_mphf: true,
            weighted: false,
            verbose: true,
            force_external_sort: false,
            tmp_dirname: PathBuf::from("sshash_tmp"),
        }
    }
}

impl BuildConfiguration {
    /// Create a new build configuration with the specified k-mer and minimizer lengths
    pub fn new(k: usize, m: usize) -> Result<Self, String> {
        let config = Self {
            k,
            m,
            ..Self::default()
        };
        config.validate()?;
        Ok(config)
    }
    
    /// Validate the configuration parameters
    pub fn validate(&self) -> Result<(), String> {
        // Check k is odd and in valid range
        if self.k % 2 == 0 {
            return Err(format!("k must be odd, got k={}", self.k));
        }
        if self.k < 3 || self.k > 63 {
            return Err(format!("k must be in range [3, 63], got k={}", self.k));
        }
        
        // Check m is less than k
        if self.m >= self.k {
            return Err(format!("m must be less than k, got m={}, k={}", self.m, self.k));
        }

        // Minimizers are u64 values (and reverse_complement_mmer needs a
        // sub-word shift), so m must fit 31 bases.
        if self.m == 0 || self.m > 31 {
            return Err(format!("m must be in range [1, 31], got m={}", self.m));
        }
        
        // Check lambda is reasonable
        if self.lambda < 1.0 || self.lambda > 100.0 {
            return Err(format!("lambda should be in range [1.0, 100.0], got {}", self.lambda));
        }
        
        Ok(())
    }
    
    /// Log configuration parameters via tracing
    pub fn print(&self) {
        tracing::info!("Build Configuration:");
        tracing::info!("  k = {}", self.k);
        tracing::info!("  m = {}", self.m);
        tracing::debug!("  seed = {}", self.seed);
        if self.num_threads == 0 {
            tracing::info!("  num_threads = all available cores");
        } else {
            tracing::info!("  num_threads = {}", self.num_threads);
        }
        tracing::debug!("  ram_limit_gib = {}", self.ram_limit_gib);
        tracing::debug!("  lambda = {}", self.lambda);
        tracing::debug!("  weighted = {}", self.weighted);
        tracing::debug!("  verbose = {}", self.verbose);
        tracing::debug!("  tmp_dirname = {:?}", self.tmp_dirname);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_default_config() {
        let config = BuildConfiguration::default();
        assert_eq!(config.k, 31);
        assert_eq!(config.m, 19);
        assert!(config.validate().is_ok());
    }
    
    #[test]
    fn test_new_config() {
        let config = BuildConfiguration::new(21, 11).unwrap();
        assert_eq!(config.k, 21);
        assert_eq!(config.m, 11);
    }
    
    #[test]
    fn test_validate_even_k() {
        let config = BuildConfiguration { k: 30, ..BuildConfiguration::default() };
        assert!(config.validate().is_err());
    }
    
    #[test]
    fn test_validate_even_m() {
        let config = BuildConfiguration { m: 20, ..BuildConfiguration::default() };
        assert!(config.validate().is_ok());
    }
    
    #[test]
    fn test_validate_m_ge_k() {
        let config = BuildConfiguration { k: 21, m: 21, ..BuildConfiguration::default() };
        assert!(config.validate().is_err());
    }
    
    #[test]
    fn test_validate_k_out_of_range() {
        let config = BuildConfiguration { k: 65, ..BuildConfiguration::default() };
        assert!(config.validate().is_err());
        
        let config = BuildConfiguration { k: 1, ..BuildConfiguration::default() };
        assert!(config.validate().is_err());
    }
}
