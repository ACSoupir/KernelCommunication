# Benchmark: soupir_moran_multi() vs multiScaleAC()

This directory contains benchmark comparisons between the legacy `soupir_moran_multi()` function and the new `multiScaleAC()` function with RcppEigen backend.

## Files

- `benchmark_comparison.Rmd` - R Markdown source file for the benchmark report
- `benchmark_comparison.html` - Generated HTML report (open in browser to view)
- `run_benchmark.R` - Runner script to execute the benchmark

## Running the Benchmark

To regenerate the benchmark report, run from the package root:

```bash
Rscript docs/benchmark/run_benchmark.R
```

Or from R:

```r
source("docs/benchmark/run_benchmark.R")
```

## Customizing Parameters

Edit `run_benchmark.R` to modify benchmark parameters:

- `h_max` - Maximum bandwidth value (default: 250)
- `h_by` - Bandwidth increment (default: 25 for speed, use 5 for full resolution)
- `perms` - Number of permutations (default: 50)
- `bench_iterations` - Number of benchmark iterations (default: 3)
- `use_subset` - Use data subset for faster execution (default: TRUE)
- `subset_size` - Number of cells to use when subsetting (default: 2000)

## Benchmark Configuration

The default benchmark uses:
- **Dataset**: ccRCCExample (subset of 2000 cells from 7122 total)
- **Bandwidth range**: 25 to 250 in increments of 25 (10 bandwidths)
- **Permutations**: 50
- **Benchmark iterations**: 3
- **Kernel**: Gaussian
- **Distance**: Euclidean

## Key Metrics

The benchmark report compares:

1. **Execution Time** - Median and mean elapsed time
2. **Memory Allocation** - Memory usage in MB
3. **Accuracy** - Absolute and relative differences in observed Moran's I values
4. **Statistical Equivalence** - Correlation, MAE, RMSE between methods

## Expected Results

Based on validation testing, the two functions should produce:
- Identical observed Moran's I values within machine precision (~1e-14)
- `multiScaleAC()` should be faster due to RcppEigen optimizations
- Memory usage may vary depending on implementation details

## Troubleshooting

If the benchmark fails:
1. Ensure all required packages are installed (devtools, bench, rmarkdown, ggplot2, dplyr, tidyr)
2. Verify the KernelCommunication package loads successfully
3. Reduce `subset_size` or increase `h_by` for faster execution
4. Check R version compatibility (requires R >= 3.5)
