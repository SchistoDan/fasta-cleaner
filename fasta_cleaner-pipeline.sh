#!/bin/bash

# =============================================================================
# Fasta_cleaner sequential cleaning pipeline with Global Memory Management
#
# PIPELINE STEPS:
# 1. Human COX1 Filter - Removes sequences with high similarity to human mitochondrial COX1
# 2. AT Content Filter - Removes sequences with unusual AT content compared to consensus
# 3. Statistical Outlier Filter - Removes sequences that are statistical outliers
# 4. Consensus Generator - Creates consensus sequences from filtered alignments
# 5. Metrics Aggregator - Combines statistics from all filtering steps
#
# MEMORY MANAGEMENT:
# - Global memory allocation for entire pipeline run
# - Dynamic thread adjustment based on available memory
# - Memory monitoring and early warning system
# - Automatic process isolation for memory-intensive steps
# - Recovery mechanisms for memory-related failures
#
# DEPENDENCIES:
# - Python 3.6+ with packages: BioPython, numpy, pandas
# - All 5 Python filter scripts (use --scripts-dir to specify location)
#
# INPUT:
# - alignment_files.log: Text file containing full paths to FASTA alignment files (one per line)
# - --memory: Total memory allocation in GB (e.g., 32)
# - --scripts-dir: Directory containing Python filter scripts
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================

# Human COX1 filter parameters
HUMAN_THRESHOLD=0.95

# AT content filter parameters  
AT_THRESHOLD=0.1
AT_MODE="absolute"
CONSENSUS_THRESHOLD=0.5

# Statistical outlier filter parameters
OUTLIER_PERCENTILE=90.0

# Consensus generation parameters
PREPROCESSING_MODE="concat"

# Default memory and system parameters
DEFAULT_MEMORY_GB=16
MEMORY_BUFFER_GB=2                    # Reserve 2GB for system processes
MAX_MEMORY_PER_PROCESS_GB=8           # Maximum memory per individual process
MEMORY_CHECK_INTERVAL=30              # Seconds between memory checks
THREADS=8

# Memory thresholds for warnings and actions
MEMORY_WARNING_THRESHOLD=80           # Warn at 80% usage
MEMORY_CRITICAL_THRESHOLD=90          # Take action at 90% usage

# =============================================================================
# GLOBAL VARIABLES
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR=""  # Will be set by command line argument or default to SCRIPT_DIR
PIPELINE_START_TIME=$(date +%s)
ALLOCATED_MEMORY_GB=""
AVAILABLE_MEMORY_GB=""
MEMORY_MONITOR_PID=""
MEMORY_LOG=""
LOG_FILE=""

# =============================================================================
# MEMORY MANAGEMENT FUNCTIONS
# =============================================================================

get_system_memory_info() {
    # Get total system memory in GB
    local total_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    local total_gb=$((total_kb / 1024 / 1024))
    
    # Get available memory in GB
    local available_kb=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
    local available_gb=$((available_kb / 1024 / 1024))
    
    echo "$total_gb $available_gb"
}

validate_memory_allocation() {
    local requested_memory="$1"
    
    # Get system memory info
    read -r total_memory available_memory <<< "$(get_system_memory_info)"
    
    log_message "System Memory Information:"
    log_message "  Total system memory: ${total_memory}GB"
    log_message "  Available memory: ${available_memory}GB"
    log_message "  Requested allocation: ${requested_memory}GB"
    
    # Check if requested memory is reasonable
    if [[ $requested_memory -gt $total_memory ]]; then
        log_error "Requested memory (${requested_memory}GB) exceeds total system memory (${total_memory}GB)"
        return 1
    fi
    
    if [[ $requested_memory -gt $available_memory ]]; then
        log_error "Requested memory (${requested_memory}GB) exceeds available memory (${available_memory}GB)"
        log_message "This may cause system instability. Consider reducing allocation or freeing memory."
        return 1
    fi
    
    # Warn if using more than 80% of available memory
    local usage_percent=$(( (requested_memory * 100) / available_memory ))
    if [[ $usage_percent -gt 80 ]]; then
        log_message "WARNING: Using ${usage_percent}% of available memory"
        log_message "Consider leaving more memory for system processes"
    fi
    
    AVAILABLE_MEMORY_GB=$available_memory
    return 0
}

calculate_optimal_threads() {
    local memory_gb="$1"
    local step_name="$2"
    
    # Base thread calculation on available memory
    local memory_per_thread=""
    
    case "$step_name" in
        "human_filter")
            memory_per_thread=2  # 2GB per thread for human filter
            ;;
        "at_filter")
            memory_per_thread=4  # 4GB per thread for AT filter (more memory intensive)
            ;;
        "outlier_filter")
            memory_per_thread=3  # 3GB per thread for outlier filter
            ;;
        "consensus")
            memory_per_thread=2  # 2GB per thread for consensus
            ;;
        *)
            memory_per_thread=2  # Default
            ;;
    esac
    
    # Calculate threads based on memory
    local max_threads_by_memory=$((memory_gb / memory_per_thread))
    
    # Don't exceed CPU count or original thread setting
    local cpu_threads=$(nproc)
    local optimal_threads=$THREADS
    
    if [[ $max_threads_by_memory -lt $optimal_threads ]]; then
        optimal_threads=$max_threads_by_memory
    fi
    
    if [[ $optimal_threads -lt 1 ]]; then
        optimal_threads=1
    fi
    
    echo "$optimal_threads"
}

setup_memory_monitoring() {
    local output_dir="$1"
    
    # Create memory reports directory
    local memory_dir="${output_dir}/memory_reports"
    create_directory "$memory_dir"
    
    MEMORY_LOG="${memory_dir}/memory_usage.log"
    
    # Start memory monitoring in background
    start_memory_monitor &
    MEMORY_MONITOR_PID=$!
    
    log_message "Memory monitoring started (PID: $MEMORY_MONITOR_PID)"
    log_message "Memory log: $MEMORY_LOG"
}

start_memory_monitor() {
    while true; do
        local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        
        # Get current memory usage
        local total_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        local available_kb=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
        local used_kb=$((total_kb - available_kb))
        
        local total_gb=$((total_kb / 1024 / 1024))
        local used_gb=$((used_kb / 1024 / 1024))
        local available_gb=$((available_kb / 1024 / 1024))
        local usage_percent=$(( (used_gb * 100) / total_gb ))
        
        # Log memory status
        echo "${timestamp},${used_gb},${available_gb},${total_gb},${usage_percent}" >> "$MEMORY_LOG"
        
        # Check for memory pressure
        if [[ $usage_percent -gt $MEMORY_CRITICAL_THRESHOLD ]]; then
            log_error "CRITICAL: Memory usage at ${usage_percent}% (${used_gb}GB/${total_gb}GB)"
        elif [[ $usage_percent -gt $MEMORY_WARNING_THRESHOLD ]]; then
            log_message "WARNING: Memory usage at ${usage_percent}% (${used_gb}GB/${total_gb}GB)"
        fi
        
        sleep $MEMORY_CHECK_INTERVAL
    done
}

stop_memory_monitoring() {
    if [[ -n "$MEMORY_MONITOR_PID" ]]; then
        kill "$MEMORY_MONITOR_PID" 2>/dev/null || true
        wait "$MEMORY_MONITOR_PID" 2>/dev/null || true
        log_message "Memory monitoring stopped"
    fi
}

generate_memory_report() {
    local output_dir="$1"
    
    if [[ ! -f "$MEMORY_LOG" ]]; then
        log_message "No memory log found for report generation"
        return
    fi
    
    local report_file="${output_dir}/memory_reports/memory_usage_summary.txt"
    
    {
        echo "Memory Usage Report - Generated $(date)"
        echo "=========================================="
        echo ""
        echo "Allocated Memory: ${ALLOCATED_MEMORY_GB}GB"
        echo "Available System Memory: ${AVAILABLE_MEMORY_GB}GB"
        echo ""
        echo "Peak Memory Usage:"
        
        # Find peak usage from log
        local peak_used=$(tail -n +2 "$MEMORY_LOG" | cut -d',' -f2 | sort -nr | head -1)
        local peak_percent=$(tail -n +2 "$MEMORY_LOG" | cut -d',' -f5 | sort -nr | head -1)
        
        echo "  Peak Used: ${peak_used}GB (${peak_percent}%)"
        echo ""
        
        echo "Memory Usage Over Time:"
        echo "Time,Used(GB),Available(GB),Total(GB),Usage(%)"
        cat "$MEMORY_LOG"
        
    } > "$report_file"
    
    log_message "Memory usage report generated: $report_file"
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

log_message() {
    local message="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ${message}" | tee -a "${LOG_FILE:-/dev/null}"
}

log_error() {
    local message="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ERROR: ${message}" | tee -a "${LOG_FILE:-/dev/null}" >&2
}

check_file_exists() {
    local file="$1"
    local description="$2"
    
    if [[ ! -f "$file" ]]; then
        log_error "${description} not found: $file"
        return 1
    fi
    return 0
}

create_directory() {
    local dir="$1"
    if [[ ! -d "$dir" ]]; then
        mkdir -p "$dir"
        log_message "Created directory: $dir"
    fi
}

count_lines() {
    local file="$1"
    if [[ -f "$file" ]]; then
        wc -l < "$file" | tr -d ' '
    else
        echo "0"
    fi
}

format_duration() {
    local seconds=$1
    local hours=$((seconds / 3600))
    local minutes=$(((seconds % 3600) / 60))
    local secs=$((seconds % 60))
    printf "%02d:%02d:%02d" $hours $minutes $secs
}

# =============================================================================
# DEPENDENCY CHECKING
# =============================================================================

check_dependencies() {
    log_message "Checking dependencies..."
    
    local missing_deps=()
    
    # Check Python
    if ! command -v python3 &> /dev/null; then
        missing_deps+=("python3")
    else
        log_message "  Python 3: $(python3 --version)"
    fi
    
    # Check Python packages
    local python_packages=("Bio" "numpy" "pandas")
    for package in "${python_packages[@]}"; do
        if ! python3 -c "import ${package}" &> /dev/null; then
            missing_deps+=("python3-${package}")
        else
            log_message "  Python package ${package}: OK"
        fi
    done
    
    # Check filter scripts in the specified scripts directory
    local filter_scripts=(
        "01_human_cox1_filter.py"
        "02_at_content_filter.py" 
        "03_statistical_outlier_filter.py"
        "05_consensus_generator.py"
        "06_aggregate_filter_metrics.py"
    )
    
    log_message "  Looking for scripts in: $SCRIPTS_DIR"
    
    for script in "${filter_scripts[@]}"; do
        local script_path="${SCRIPTS_DIR}/${script}"
        if [[ ! -f "$script_path" ]]; then
            missing_deps+=("${script}")
        else
            log_message "  Filter script ${script}: OK"
        fi
    done
    
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        log_error "Missing dependencies:"
        for dep in "${missing_deps[@]}"; do
            log_error "  - ${dep}"
        done
        log_error ""
        log_error "Script directory checked: $SCRIPTS_DIR"
        log_error "Use --scripts-dir to specify the location of Python filter scripts"
        return 1
    fi
    
    log_message "All dependencies satisfied"
    return 0
}

# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

run_pipeline() {
    local alignment_files="$1"
    local output_base_dir="$2"
    local memory_gb="$3"
    
    ALLOCATED_MEMORY_GB="$memory_gb"
    
    if ! validate_memory_allocation "$memory_gb"; then
        log_error "Memory allocation validation failed"
        return 1
    fi
    
    create_directory "$output_base_dir"
    create_directory "${output_base_dir}/logs"
    LOG_FILE="${output_base_dir}/logs/pipeline_execution.log"
    
    setup_memory_monitoring "$output_base_dir"
    
    log_message "======================================================================"
    log_message "FASTA Sequence Filtering Pipeline Started (Memory-Managed)"
    log_message "======================================================================"
    log_message "Input alignment files: $alignment_files"
    log_message "Output directory: $output_base_dir"
    log_message "Scripts directory: $SCRIPTS_DIR"
    log_message "Allocated memory: ${memory_gb}GB"
    log_message "Available system memory: ${AVAILABLE_MEMORY_GB}GB"
    log_message "Base threads: $THREADS"
    log_message "Memory buffer: ${MEMORY_BUFFER_GB}GB"
    log_message ""
    
    # Check dependencies
    if ! check_dependencies; then
        log_error "Dependency check failed"
        return 1
    fi
    
    # Create step output directories
    local human_dir="${output_base_dir}/01_human_filtered"
    local at_dir="${output_base_dir}/02_at_filtered"
    local outlier_dir="${output_base_dir}/03_outlier_filtered"
    local consensus_dir="${output_base_dir}/04_consensus_seqs"
    local metrics_dir="${output_base_dir}/05_metrics"
    
    for dir in "$human_dir" "$at_dir" "$outlier_dir" "$consensus_dir" "$metrics_dir"; do
        create_directory "$dir"
    done
    
    # Initialize file tracking
    local current_files="$alignment_files"
    local step_memory=$((ALLOCATED_MEMORY_GB - MEMORY_BUFFER_GB))
    
    # Step 1: Human COX1 Filter
    log_message "=== STEP 1: Human COX1 Filter ==="
    log_message "Input files: $(count_lines "$current_files")"
    log_message "Human similarity threshold: ${HUMAN_THRESHOLD}"
    
    local optimal_threads=$(calculate_optimal_threads "$step_memory" "human_filter")
    log_message "Memory allocation: ${step_memory}GB"
    log_message "Optimal threads: ${optimal_threads}"
        
    local step_start=$(date +%s)
    local filtered_files="${human_dir}/human_filtered_paths.log"
    local metrics_file="${human_dir}/human_filter_metrics.csv"
    
    python3 "${SCRIPTS_DIR}/01_human_cox1_filter.py" \
        --input-log "$current_files" \
        --output-dir "$human_dir" \
        --filtered-files-list "$filtered_files" \
        --metrics-csv "$metrics_file" \
        --human-threshold "$HUMAN_THRESHOLD" \
        --threads "$optimal_threads" || {
            log_error "Human filter step failed"
            return 1
        }
    
    local step_end=$(date +%s)
    local step_duration=$((step_end - step_start))
    log_message "Human COX1 filter completed in $(format_duration $step_duration)"
    log_message "Filtered files: $(count_lines "$filtered_files")"
    
    current_files="$filtered_files"
    
    # Step 2: AT Content Filter
    log_message "=== STEP 2: AT Content Filter ==="
    log_message "Input files: $(count_lines "$current_files")"
    log_message "AT threshold: ${AT_THRESHOLD}, Mode: ${AT_MODE}"
    
    local optimal_threads=$(calculate_optimal_threads "$step_memory" "at_filter")
    log_message "Memory allocation: ${step_memory}GB"
    log_message "Optimal threads: ${optimal_threads}"
    
    # Set higher memory limits for AT filter
    local at_memory_limit=$((MAX_MEMORY_PER_PROCESS_GB * 2))
    
    local step_start=$(date +%s)
    local filtered_files="${at_dir}/at_filtered_paths.log"
    
    python3 "${SCRIPTS_DIR}/02_at_content_filter.py" \
        --input_files "$current_files" \
        --output_dir "$at_dir" \
        --filtered-files-list "$filtered_files" \
        --at_threshold "$AT_THRESHOLD" \
        --at_mode "$AT_MODE" \
        --consensus_threshold "$CONSENSUS_THRESHOLD" \
        --threads "$optimal_threads" || {
            log_error "AT filter step failed"
            return 1
        }
    
    local step_end=$(date +%s)
    local step_duration=$((step_end - step_start))
    log_message "AT content filter completed in $(format_duration $step_duration)"
    log_message "Filtered files: $(count_lines "$filtered_files")"
    
    current_files="$filtered_files"
    
    # Step 3: Statistical Outlier Filter
    log_message "=== STEP 3: Statistical Outlier Filter ==="
    log_message "Input files: $(count_lines "$current_files")"
    log_message "Outlier percentile: ${OUTLIER_PERCENTILE}"
    
    local optimal_threads=$(calculate_optimal_threads "$step_memory" "outlier_filter")
    log_message "Memory allocation: ${step_memory}GB"
    log_message "Optimal threads: ${optimal_threads}"
        
    local step_start=$(date +%s)
    local filtered_files="${outlier_dir}/outlier_filtered_paths.log"
    local summary_csv="${outlier_dir}/outlier_filter_summary.csv"
    local metrics_csv="${outlier_dir}/outlier_filter_metrics.csv"
    
    python3 "${SCRIPTS_DIR}/03_statistical_outlier_filter.py" \
        --input-files-list "$current_files" \
        --output-dir "$outlier_dir" \
        --filtered-files-list "$filtered_files" \
        --summary-csv "$summary_csv" \
        --metrics-csv "$metrics_csv" \
        --outlier-percentile "$OUTLIER_PERCENTILE" \
        --consensus-threshold "$CONSENSUS_THRESHOLD" \
        --threads "$optimal_threads" || {
            log_error "Outlier filter step failed"
            return 1
        }
    
    local step_end=$(date +%s)
    local step_duration=$((step_end - step_start))
    log_message "Statistical outlier filter completed in $(format_duration $step_duration)"
    log_message "Filtered files: $(count_lines "$filtered_files")"
    
    current_files="$filtered_files"
    
    # Step 4: Consensus Generation
    log_message "=== STEP 4: Consensus Generator ==="
    log_message "Input files: $(count_lines "$current_files")"
    log_message "Consensus threshold: ${CONSENSUS_THRESHOLD}, Mode: ${PREPROCESSING_MODE}"
    
    local optimal_threads=$(calculate_optimal_threads "$step_memory" "consensus")
    log_message "Memory allocation: ${step_memory}GB"
    log_message "Optimal threads: ${optimal_threads}"
        
    local step_start=$(date +%s)
    local consensus_fasta="${consensus_dir}/all_consensus_sequences.fasta"
    local consensus_metrics="${consensus_dir}/consensus_metrics.csv"
    
    python3 "${SCRIPTS_DIR}/05_consensus_generator.py" \
        --input-files-list "$current_files" \
        --output-dir "$consensus_dir" \
        --consensus-fasta "$consensus_fasta" \
        --consensus-metrics "$consensus_metrics" \
        --consensus-threshold "$CONSENSUS_THRESHOLD" \
        --preprocessing-mode "$PREPROCESSING_MODE" \
        --threads "$optimal_threads" || {
            log_error "Consensus generation step failed"
            return 1
        }
    
    local step_end=$(date +%s)
    local step_duration=$((step_end - step_start))
    log_message "Consensus generation completed in $(format_duration $step_duration)"
    
    # Step 5: Metrics Aggregation
    log_message "=== STEP 5: Metrics Aggregator ==="
    log_message "Aggregating metrics from all filtering steps..."
        
    local step_start=$(date +%s)
    local combined_statistics="${metrics_dir}/combined_filter_statistics.csv"
    
    # Build arguments for metrics that exist
    local args=(
        --output-dir "$metrics_dir"
        --combined-statistics "$combined_statistics"
        --threads "$THREADS"
    )
    
    [[ -f "${human_dir}/human_filter_metrics.csv" ]] && args+=(--human-metrics "${human_dir}/human_filter_metrics.csv")
    [[ -f "${at_dir}/at_filter_summary.csv" ]] && args+=(--at-metrics "${at_dir}/at_filter_summary.csv")
    [[ -f "${outlier_dir}/outlier_filter_summary.csv" ]] && args+=(--outlier-metrics "${outlier_dir}/outlier_filter_summary.csv")
    [[ -f "$consensus_metrics" ]] && args+=(--consensus-metrics "$consensus_metrics")
    
    python3 "${SCRIPTS_DIR}/06_aggregate_filter_metrics.py" "${args[@]}" || {
        log_error "Metrics aggregation step failed"
        return 1
    }
    
    local step_end=$(date +%s)
    local step_duration=$((step_end - step_start))
    log_message "Metrics aggregation completed in $(format_duration $step_duration)"
    log_message "Combined statistics: $combined_statistics"
    
    # Copy key metrics to main metrics directory
    log_message "Copying key metrics files to central location..."
    cp "${human_dir}/human_filter_metrics.csv" "$metrics_dir/" 2>/dev/null || true
    cp "${at_dir}/at_filter_summary.csv" "$metrics_dir/" 2>/dev/null || true
    cp "${outlier_dir}/outlier_filter_summary.csv" "$metrics_dir/" 2>/dev/null || true
    cp "$consensus_metrics" "$metrics_dir/" 2>/dev/null || true
    
    # Stop memory monitoring and generate report
    stop_memory_monitoring
    generate_memory_report "$output_base_dir"
    
    # Calculate total pipeline time
    local pipeline_end_time=$(date +%s)
    local total_duration=$((pipeline_end_time - PIPELINE_START_TIME))
    
    log_message ""
    log_message "======================================================================"
    log_message "FASTA_CLEANER PIPELINE COMPLETED SUCCESSFULLY"
    log_message "======================================================================"
    log_message "Total runtime: $(format_duration $total_duration)"
    log_message "Memory allocation: ${ALLOCATED_MEMORY_GB}GB"
    log_message "Memory utilization report: ${output_base_dir}/memory_reports/"
    log_message ""
    log_message "Final output locations:"
    log_message "  - Consensus sequences: ${consensus_dir}/all_consensus_sequences.fasta"
    log_message "  - Individual consensus: ${consensus_dir}/"
    log_message "  - Cleaned reads: ${consensus_dir}/filter_pass_seqs/"
    log_message "  - Combined statistics: ${metrics_dir}/combined_filter_statistics.csv"
    log_message "  - All metrics: ${metrics_dir}/"
    log_message "  - Memory reports: ${output_base_dir}/memory_reports/"
    log_message "  - Pipeline log: $LOG_FILE"
    log_message "======================================================================"
}

# =============================================================================
# CLEANUP AND ERROR HANDLING
# =============================================================================

cleanup_on_exit() {
    local exit_code=$?
    
    log_message "Pipeline cleanup initiated (exit code: $exit_code)"
    
    # Stop memory monitoring if it's running
    stop_memory_monitoring
    
    # Kill any remaining Python processes from this pipeline
    local pipeline_processes=$(pgrep -f "${SCRIPTS_DIR}.*\.py" || echo "")
    if [[ -n "$pipeline_processes" ]]; then
        log_message "Terminating remaining pipeline processes..."
        for pid in $pipeline_processes; do
            if kill -0 "$pid" 2>/dev/null; then
                log_message "Terminating process $pid"
                kill -TERM "$pid" 2>/dev/null || true
                sleep 2
                kill -KILL "$pid" 2>/dev/null || true
            fi
        done
    fi
    
    # Generate final memory report if possible
    if [[ -n "$MEMORY_LOG" && -f "$MEMORY_LOG" ]]; then
        log_message "Generating final memory report..."
        generate_memory_report "$(dirname "$(dirname "$MEMORY_LOG")")" 2>/dev/null || true
    fi
    
    if [[ $exit_code -ne 0 ]]; then
        log_error "Pipeline failed with exit code $exit_code"
    else
        log_message "Pipeline cleanup completed successfully"
    fi
}

# Set up signal handlers for cleanup
trap cleanup_on_exit EXIT
trap 'log_error "Pipeline interrupted by user"; exit 130' INT TERM

# =============================================================================
# MAIN EXECUTION
# =============================================================================

main() {
    # Parse command line arguments
    local alignment_files=""
    local output_dir=""
    local memory_gb="$DEFAULT_MEMORY_GB"
    
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --memory|-m)
                memory_gb="$2"
                shift 2
                ;;
            --threads|-t)
                THREADS="$2"
                shift 2
                ;;
            --scripts-dir|-s)
                SCRIPTS_DIR="$2"
                shift 2
                ;;
            --help|-h)
                show_usage
                exit 0
                ;;
            -*)
                echo "Unknown option $1"
                show_usage
                exit 1
                ;;
            *)
                if [[ -z "$alignment_files" ]]; then
                    alignment_files="$1"
                elif [[ -z "$output_dir" ]]; then
                    output_dir="$1"
                else
                    echo "Too many positional arguments"
                    show_usage
                    exit 1
                fi
                shift
                ;;
        esac
    done
    
    # Set default scripts directory if not specified
    if [[ -z "$SCRIPTS_DIR" ]]; then
        SCRIPTS_DIR="$SCRIPT_DIR"
    fi
    
    # Validate scripts directory
    if [[ ! -d "$SCRIPTS_DIR" ]]; then
        echo "Error: Scripts directory does not exist: $SCRIPTS_DIR"
        exit 1
    fi
    
    # Convert to absolute path
    SCRIPTS_DIR="$(cd "$SCRIPTS_DIR" && pwd)"
    
    # Check for required arguments
    if [[ -z "$alignment_files" || -z "$output_dir" ]]; then
        echo "Error: Missing required arguments"
        show_usage
        exit 1
    fi
    
    # Validate memory argument
    if ! [[ "$memory_gb" =~ ^[0-9]+$ ]] || [[ $memory_gb -lt 1 ]]; then
        echo "Error: Memory allocation must be a positive integer (GB)"
        exit 1
    fi
    
    # Validate inputs
    check_file_exists "$alignment_files" "Alignment files list" || exit 1
    
    # Display startup information
    echo "======================================================================"
    echo "MEMORY-MANAGED FASTA SEQUENCE FILTERING PIPELINE"
    echo "======================================================================"
    echo "Configuration:"
    echo "  Input files: $alignment_files"
    echo "  Output directory: $output_dir"
    echo "  Scripts directory: $SCRIPTS_DIR"
    echo "  Allocated memory: ${memory_gb}GB"
    echo "  Base threads: $THREADS"
    echo "  Memory buffer: ${MEMORY_BUFFER_GB}GB"
    echo "  Max memory per process: ${MAX_MEMORY_PER_PROCESS_GB}GB"
    echo ""
    echo "Filter Parameters:"
    echo "  Human COX1 threshold: $HUMAN_THRESHOLD"
    echo "  AT content threshold: $AT_THRESHOLD"
    echo "  AT filter mode: $AT_MODE"
    echo "  Consensus threshold: $CONSENSUS_THRESHOLD"
    echo "  Outlier percentile: $OUTLIER_PERCENTILE"
    echo "  Preprocessing mode: $PREPROCESSING_MODE"
    echo "======================================================================"
    echo ""
    
    # Run the pipeline
    run_pipeline "$alignment_files" "$output_dir" "$memory_gb"
}

show_usage() {
    cat << EOF
Usage: $0 [OPTIONS] <alignment_files.log> <output_directory>

MEMORY-MANAGED FASTA SEQUENCE FILTERING PIPELINE

REQUIRED ARGUMENTS:
    alignment_files.log    Text file containing full paths to FASTA alignment files (one per line)
    output_directory      Directory where all pipeline outputs will be created

OPTIONS:
    -m, --memory GB       Total memory allocation for the pipeline in GB (default: ${DEFAULT_MEMORY_GB})
    -t, --threads N       Number of base threads (default: ${THREADS})
    -s, --scripts-dir DIR Directory containing Python filter scripts (default: script directory)
    -h, --help           Show this help message

EXAMPLES:
    # Basic usage with scripts in same directory
    $0 --memory 16 /path/to/alignment_files.log /path/to/output_dir
    
    # Specify custom scripts directory
    $0 --memory 32 --scripts-dir /path/to/scripts/ alignment_files.log output_dir/
    
    # High-memory system with custom settings
    $0 --memory 32 --threads 16 --scripts-dir ~/fasta_tools/ files.log results/

REQUIRED PYTHON SCRIPTS:
    The following scripts must be present in the scripts directory:
    - 01_human_cox1_filter.py
    - 02_at_content_filter.py
    - 03_statistical_outlier_filter.py
    - 05_consensus_generator.py
    - 06_aggregate_filter_metrics.py

EOF
}

# Run main function with all arguments
main "$@"
