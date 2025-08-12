# FASTA Cleaner - Modular Pipeline

A comprehensive Python toolkit for cleaning and analysing FASTA sequence alignments using multiple filtering approaches. This modular pipeline is designed to help systematically identify and remove possible contaminant reads from their alignments.

## Overview

FASTA Cleaner now consists of **6 specialised modules** that can be used independently or sequentially as part of a complete pipeline:

1. **Human COX1 Filter** - Removes human contamination sequences
2. **AT Content Filter** - Filters sequences with divergent nucleotide composition  
3. **Statistical Outlier Filter** - Removes sequences that are statistical outliers
4. **Reference Filter** (Optional) - Compares sequences against reference standards
5. **Consensus Generator** - Creates final consensus sequences and cleaned reads
6. **Metrics Aggregator** - Combines statistics from all filtering steps

## Pipeline Architecture

```
Input FASTA Files
       ↓
Human COX1 Filter → human_filtered.fasta + metrics
       ↓
AT Content Filter → at_filtered.fasta + metrics  
       ↓
Statistical Outlier Filter → outlier_filtered.fasta + metrics
       ↓
Reference Filter (Optional) → reference_filtered.fasta + metrics
       ↓
Consensus Generator → consensus.fasta + cleaned_reads.fasta + metrics
       ↓
Metrics Aggregator → combined_statistics.csv
```

## Features

### Core Capabilities
- **Modular Design**: Each filter can be used independently or as part of the complete pipeline
- **Memory-Aware Processing**: Automatic handling of large files with intelligent chunking and sequential processing
- **Parallel Processing**: Multi-threaded processing with SLURM integration
- **Comprehensive Logging**: Detailed logs and metrics for each processing step
- **Flexible Configuration**: Customizable thresholds and parameters for each filter


### Filtering Methods

#### 1. Human COX1 Contamination Detection
**How it works:**
- Compares input sequences against a hard-coded human COX1 sequences using local alignment
- Calculates similarity scores based on sequence identity over aligned regions
- Removes sequences that exceed the similarity threshold, indicating potential human contamination
**Details:**
- Calculates percentage identity: `(identical_bases / aligned_length) × 100`
- Only considers high-quality alignments above minimum length thresholds
- Handles gaps and ambiguous bases appropriately in similarity calculations
**Key Parameters:**
- `--human-threshold` (default: 0.95): Similarity cutoff (0.0-1.0)
  - 0.95 = Remove sequences with ≥95% similarity to human COX1
- `--threads`: Parallel processing for large datasets

---

#### 2. AT Content Analysis
**How it works:**
- Generates consensus sequence from input alignment using specified threshold
- Calculates AT content for both consensus and individual sequences
- Compares AT content only in overlapping (non-gap) regions between sequences
- Applies filtering based on selected mode and threshold
**Details:**
```
1. Generate consensus: For each position, use most frequent base if ≥ consensus_threshold
2. Calculate AT content:
   - Consensus AT = (A + T bases) / (total non-gap bases in consensus)
   - Sequence AT = (A + T bases) / (total non-gap bases in overlap with consensus)
3. Calculate difference: |sequence_AT - consensus_AT|
4. Apply filtering based on mode
```
**Filtering Modes:**
- **`absolute`** (default): Remove sequences if AT content differs from consensus by more than threshold in either direction
  - Removes sequences with: `|seq_AT - consensus_AT| > threshold`
  - Best for removing sequences with unusual base composition
- **`higher`**: Remove only sequences with AT content above consensus + threshold
  - Removes sequences with: `seq_AT > (consensus_AT + threshold)`
  - Useful for removing AT-rich contamination (e.g., some bacterial sequences)
- **`lower`**: Remove only sequences with AT content below consensus - threshold
  - Removes sequences with: `seq_AT < (consensus_AT - threshold)`
  - Useful for removing GC-rich contamination or poor-quality sequences

**Key Parameters:**
- `--at-threshold` (default: 0.1): Maximum allowed AT content difference (0.0-1.0)
  - 0.1 = Moderate (remove sequences differing by >10%)
- `--consensus-threshold` (default: 0.5): Threshold for consensus generation (0.0-1.0)
  - 0.5 = Use bases present in ≥50% of sequences

---

#### 3. Statistical Outlier Detection
**How it works:**
- Generates consensus sequence and calculates position-specific residue frequencies
- Computes deviation scores for each sequence against the consensus
- Uses percentile-based thresholds to identify and remove statistical outliers
**Details:**
**Step 1: Position Frequency Analysis**
```
For each alignment position i:
  freq[i][base] = count(base_at_position_i) / total_sequences_with_base_at_i
```
**Step 2: Deviation Score Calculation**
- **Unweighted Deviation**: Simple mismatch proportion
  ```
  unweighted_score = mismatches / total_comparable_positions
  ```
- **Weighted Deviation**: Conservation-weighted mismatch score
  ```
  For each position i where seq[i] != consensus[i]:
    weight = frequency_of_consensus_base_at_position_i
    weighted_score += weight
  weighted_score = weighted_score / sum_of_all_weights
  ```
**Step 3: Outlier Detection**
- Calculate percentile threshold for both weighted and unweighted scores
- Remove sequences exceeding either threshold
**Key Parameters:**
- `--outlier-percentile` (default: 90.0): Percentile threshold for outlier detection (0-100)
  - 90.0 = Remove sequences in top 10% of deviation scores
- `--consensus-threshold` (default: 0.5): Consensus generation threshold

---

#### 4. Reference Sequence Comparison (Optional)
**How it works:**
- Compares input sequences against high-quality reference sequences using the same statistical framework as outlier detection
- Supports two filtering modes - keeping similar sequences OR removing similar sequences
- Requires external reference sequence files matching sample names

**Details:**
**Step 1: Reference Preparation**
- Loads reference sequence from `{sample_name}_reference.fasta`
- Pads reference to match alignment length
- Reference provides the comparison standard
**Step 2: Frequency Calculation**
- Uses input alignment to calculate position-specific frequencies
- Provides conservation weighting based on the actual dataset
**Step 3: Deviation Scoring**
- **Unweighted**: `mismatches_with_reference / comparable_positions`
- **Weighted**: Conservation-weighted deviations from reference
  ```
  For each position i where seq[i] != reference[i]:
    weight = frequency_of_reference_base_in_alignment_at_position_i
    weighted_score += weight
  ```
**Step 4: Filtering Based on Mode**
- **`keep_similar` mode**: Remove sequences that deviate excessively from reference (outlier removal)
  - Calculate high percentile threshold (e.g., 90th percentile)
  - Remove sequences with deviation > threshold
- **`remove_similar` mode**: Remove sequences that closely match reference (contamination removal)
  - Calculate low percentile threshold (e.g., 10th percentile when threshold-percentile=90)
  - Remove sequences with deviation ≤ threshold
**Key Parameters:**
- `--filter-mode` (required): Choose filtering strategy
  - `keep_similar`: Keep sequences similar to reference, remove outliers
  - `remove_similar`: Remove sequences similar to reference (contamination filtering)
- `--reference-dir`: Directory containing reference sequences
  - Must contain files named `{sample}_reference.fasta`
  - Each file should contain single reference sequence
- `--threshold-percentile` (default: 90.0): Percentile threshold behavior depends on mode:
  - `keep_similar`: Remove sequences above 90th percentile (remove most different 10%)
  - `remove_similar`: Remove sequences below 10th percentile (remove most similar 10%)
- `--consensus-threshold` (default: 0.5): For frequency calculations from input data
**Reference File Requirements:**
```
reference_dir/
├── sample1_reference.fasta    # Single reference sequence
├── sample2_reference.fasta    # One file per sample
└── sample3_reference.fasta    # Filename must match sample base name
```

#### 5. Consensus Generation
- Creates high-quality consensus sequences
- Calculates comprehensive coverage statistics
- Generates cleaned sequence collections
- No line wrapping in FASTA output for compatibility

#### 6. Metrics Aggregation
- Combines statistics from all filtering steps
- Provides comprehensive sample-level summaries
- Tracks sequences through the entire pipeline

## Installation

### Prerequisites
- Python 3.6 or higher
- BioPython
- NumPy  
- Pandas

```bash
pip install biopython numpy pandas
```

### Installation Steps
1. Clone this repository:
```bash
git clone https://github.com/bge-barcoding/fasta-cleaner.git
cd fasta-cleaner
```

2. Install dependencies:
```bash
pip install biopython numpy pandas
```

## Usage

### Input Format
All scripts expect a **file list** as input - a text file containing full paths to FASTA alignment files (output by the BGEE pipeline), one per line:

```
/path/to/sample1_align.fasta
/path/to/sample2_align.fasta
/path/to/sample3_align.fasta
```

### Fasta-cleaner (sequential) usage
```bash
# Step 1: Human COX1 filtering
python human_cox1_filter.py \
    --input-files-list raw_alignments.txt \
    --output-dir step1_human_filtered/ \
    --filtered-files-list human_filtered.log \
    --metrics-csv human_metrics.csv \
    --threads 4

# Step 2: AT content filtering  
python at_content_filter.py \
    --input-files human_filtered.log \
    --output-dir step2_at_filtered/ \
    --at-threshold 0.1 \
    --at-mode absolute \
    --threads 8

# Step 3: Statistical outlier filtering
python statistical_outlier_filter.py \
    --input-files-list step2_at_filtered/at_filtered_paths.log \
    --output-dir step3_outlier_filtered/ \
    --filtered-files-list outlier_filtered.log \
    --summary-csv outlier_summary.csv \
    --metrics-csv outlier_metrics.csv \
    --threads 8

# Optional: Reference filtering (choose one mode)
# For quality control (keep sequences similar to high-quality references):
# python reference_filter.py \
#     --input-files-list outlier_filtered.log \
#     --output-dir reference_filtered/ \
#     --filtered-files-list reference_filtered.log \
#     --metrics-csv reference_metrics.csv \
#     --reference-dir reference_sequences/ \
#     --filter-mode keep_similar \
#     --threads 8

# For contamination removal (remove sequences similar to contaminants):
# python reference_filter.py \
#     --input-files-list outlier_filtered.log \
#     --output-dir contamination_filtered/ \
#     --filtered-files-list contamination_filtered.log \
#     --metrics-csv contamination_metrics.csv \
#     --reference-dir contaminant_sequences/ \
#     --filter-mode remove_similar \
#     --threads 8

# Step 4: Generate final consensus sequences
python consensus_generator.py \
    --input-files-list outlier_filtered.log \
    --output-dir final_consensus/ \
    --consensus-fasta all_consensus_sequences.fasta \
    --consensus-metrics consensus_metrics.csv \
    --threads 8

# Step 5: Aggregate all metrics
python aggregate_filter_metrics.py \
    --human-metrics human_metrics.csv \
    --at-metrics step2_at_filtered/at_filter_summary.csv \
    --outlier-metrics outlier_summary.csv \
    --consensus-metrics consensus_metrics.csv \
    --output-dir final_results/ \
    --combined-statistics pipeline_statistics.csv
```

## Output Files
- **Consensus sequences**: Individual FASTAs and combined multi-FASTA
- **Cleaned reads**: Individual cleaned sequence files
- **Combined statistics**: CSV metrics


## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

### Development Guidelines
- Each module should remain independently functional
- Maintain consistent command-line interfaces
- Include comprehensive logging and error handling
- Add tests for new functionality

## License

This project is licensed under the MIT License - see the LICENSE file for details.

**Authors**: Ben Price & Daniel Parsons @ NHMUK  
**Version**: 2.0.0 (Modular Pipeline)  
**License**: MIT
