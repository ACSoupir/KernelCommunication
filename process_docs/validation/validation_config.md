# Validation Framework for multiScaleAC() vs soupir_moran_multi()

## Test Configuration

### Datasets
- Synthetic data (N = 100, 200, 500, 1000, 2000)
- ccRCCExample (full: 7,122 cells)

### Bandwidth Configurations
- Small: c(5, 10, 20)
- Medium: c(5, 10, 20, 40, 80)
- Large: c(10, 20, 40, 80, 160)

### Validation Threshold
- **Pass:** max_abs_diff < 1e-14
- **Warning:** 1e-14 <= max_abs_diff < 1e-10
- **Fail:** max_abs_diff >= 1e-10
