# SMFnorm
Normalization of accessibility methylation data from genome-wide, single molecule, sequence data such as MAPit, NOMe-seq, Fiber-seq, etc..

# SMFnorm: Single Molecule Footprinting Data Normalization Package

## Overview
SMFnorm is a specialized R package designed for normalizing methylation data from genome-wide footprinting, such as MAPit-CpGiant or MAPit-RRBS, NOMe-seq or Fiber-seq. It implements a quantile-bin normalization that removes differences in enzyme labeling efficiency between samples while preserving local biological signals such as differentially accessible regions.

## Core Functions

### Data Loading and Preparation
load_data(): Efficiently loads and processes methylation call files with both C++ and R implementations
- C++ implementation provides up to 85% faster loading times
- Graceful fallback to pure R when C++ compilation isn't available
- Optional parallel processing with multiple cores
- Flexible sample selection through a sample sheet format
- Handles both GCH and HCG methylation contexts
- Supports filtering by group IDs for targeted analysis
- Performs initial data validation and organization
- Returns a structured list of samples with attached metadata

### Performance
The package implements a high-performance C++ reader that can significantly reduce loading times for large methylation files:

| Implementation | Time per file | Time for 10 files |
|----------------|---------------|-------------------|
| C++ loading    | ~19 seconds   | ~3.5 minutes      |
| R loading      | ~2.5 minutes  | ~25 minutes       |

For users without C++ compilation tools, the package automatically falls back to a pure R implementation, ensuring compatibility across all environments.

### Example Usage to Load
```r
# Load all samples using C++ implementation (default)
all_samples <- load_data(
  dir_path = "path/to/methylation_files/",
  sample_sheet = "sample_metadata.csv"
)

# Load only specific groups
group_samples <- load_data(
  dir_path = "path/to/methylation_files/",
  sample_sheet = "sample_metadata.csv",
  groups = c("M1", "M2")
)

# Filter by type and use pure R implementation
hcg_samples <- load_data(
  dir_path = "path/to/methylation_files/",
  sample_sheet = "sample_metadata.csv",
  type = "HCG",
  use_cpp = FALSE
)

# Load a single file
single_data <- load_single_file(
  single_file = "path/to/methylation_file.gz"
)

```
### Normalization Process

#### 1. Coverage Filter
- `normalize_methylation_data()` keeps a site only if every sample has at least
  `min_coverage` reads there, so all samples keep the same sites
- Read counts and coverage are not rescaled: `cov` stays the observed depth
- For depth-scaled coverage tracks in a genome browser, use
  `create_bigwig(type = "cov", normalize_coverage = TRUE)`

#### 2. Methylation Rate Normalization (`normalize_methylation_rates`)
- Quantile-bin normalization of labeling-efficiency differences
- Two steps:
  - Within-group: replicates of each group are normalized to one another
  - Between-group (a choice; see below): all samples are normalized together,
    removing global differences between conditions
- Sites are binned by mean rate (dynamic number of equal-count bins). Within a
  bin, sample *s*'s rate at site *i* becomes

  x′ = m̄ + α (x − mₛ)

  where mₛ is sample *s*'s own mean over the bin and m̄ is the average of
  those means across samples.
  - Replacing mₛ by m̄ is the normalization. It shifts each sample, bin by bin,
    to the common level, and is applied in full at every α. Local differences
    such as DMRs survive because their sites are a tiny fraction of any bin.
  - α only scales each site's deviation from its sample's bin mean. At α = 1
    (default) deviations are kept exactly; this is *not* the identity, since
    every sample is still shifted to the common level. At α = 0 sites collapse
    to the bin's common level.
- Within a bin the deviations are mostly sampling noise, so α < 1 shrinks noise
  and makes downstream tests overconfident (more false positives). Keep α = 1
  for differential analysis.
- `rate_between_groups` has no default and must be chosen. The between-group
  step cannot tell a global efficiency difference from a genuine broad
  biological difference: `TRUE` removes both (local differences such as DMRs
  are kept). `FALSE` keeps both, and a technical difference between groups
  then gives more false positives than unnormalized data, because the
  replicates have been made to agree. Use `TRUE` when a difference in group
  levels is more likely technical, `FALSE` when a broad biological difference
  is expected. `plot_efficiency(data, reference = "all")` shows the difference
  between groups.

### Visualization and Diagnostics

#### Within-Group Report (`plot_within_group_report`)
- Focused on technical replicate normalization
- Includes:
  - QQ plots comparing replicates
  - Density distributions
  - Multiple visualizations per experimental group

#### Between-Group Report (`plot_between_group_report`)
- Evaluates group-level normalization effects
- Features:
  - Group mean comparisons
  - Distribution plots
  - Box plots showing group-level changes

### Key Features
1. **Efficient Processing**
   - Optimized for large datasets
   - Uses data.table for fast operations
   - Memory-efficient visualization using hexbin plots

2. **Flexible Normalization**
   - Adjustable parameters for different experimental needs
   - Structure preservation to maintain biological signals
   - Independent control of within and between-group normalization

3. **Comprehensive Visualization**
   - Dense data visualization using hexbin plots
   - Separate reports for within and between-group effects
   - Customizable plotting parameters

4. **Quality Control**
   - Diagnostic information during normalization
   - Visual assessment of normalization effects
   - Multiple metrics for evaluating normalization success

## Technical Details

### Normalization Parameters
- `sites_per_quantile`: Controls quantile bin size (default = 1000)
- `within_alpha`: Shrinkage of site-level deviations in the within-group step (default = 1, none)
- `between_alpha`: Shrinkage of site-level deviations in the between-group step (default = 1, none)
- Dynamic quantile calculation based on dataset size

### Visualization Options
- Customizable color schemes (viridis color palettes)
- Adjustable bin sizes for density plots
- Optional data sampling for quick visualization
- PDF export capability for reports

## Future Developments
- Additional normalization methods
- Enhanced visualization options
- Performance optimizations for very large datasets
- Extended diagnostic capabilities

## Usage
The package is designed for methylation data analysis workflows, particularly suited for:
- MAPit experimental data
- Multi-group experimental designs
- Datasets with technical replicates
- Large-scale methylation studies

This implementation provides a robust framework for normalizing methylation data while maintaining biological significance and offering comprehensive diagnostic tools for quality assessment.
