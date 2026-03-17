#!/bin/bash

# =============================================================================
# easy_net.sh — SparCC co-occurrence network pipeline with ASV filtering
# =============================================================================
# Usage:
#   bash easy_net.sh -d <folder> -s <sparcc_path> [--auto] [-n <min_seqs>] [-ab <min_rel_abund_%>]
#
# Options:
#   -d  <folder>     Folder containing .txt ASV tables (required)
#   -s  <path>       Path to the SparCC directory (required)
#   --auto           Automatically find the lowest abundance cutoff
#                    that keeps all tables below 1000 ASVs
#   -n  <int>        Minimum number of reads per ASV per sample
#   -ab <float>      Minimum relative abundance in % per sample (e.g. 0.1 = 0.1%)
#   -h               Show this help message
#
# Examples:
#   bash easy_net.sh -d ./data -s /opt/SparCC --auto
#   bash easy_net.sh -d ./data -s /opt/SparCC -ab 0.1
#   bash easy_net.sh -d ./data -s /opt/SparCC -n 5 -ab 0.05
#   bash easy_net.sh -d ./data -s /opt/SparCC --auto -n 3
#
# Notes:
#   - Filters -n and -ab are applied together (AND) if both are provided.
#   - --auto uses absolute counts; -n is applied before counting if provided.
# =============================================================================

SPARCC_PATH=""

log_info()  { echo "[INFO]  $*"; }
log_ok()    { echo "[OK]    $*"; }
log_warn()  { echo "[WARN]  $*"; }
log_error() { echo "[ERROR] $*"; }
log_step()  { echo ""; echo ">>> $*"; }

# --- Defaults ---
folder=""
max_asv=1000
auto_mode=false
min_seqs=""
min_ab_pct=""

# =============================================================================
# Help
# =============================================================================
usage() {
cat <<EOF

easy_net.sh — SparCC co-occurrence network pipeline with ASV filtering

USAGE:
  bash easy_net.sh -d <folder> -s <sparcc_path> [--auto] [-n <min_seqs>] [-ab <min_ab_%>]

REQUIRED OPTIONS:
  -d  <folder>         Folder containing the .txt ASV tables
  -s  <path>           Path to the SparCC directory (e.g. /opt/SparCC)

OPTIONAL OPTIONS:
  --auto               Automatically find the lowest absolute abundance cutoff
                       so that all tables have < ${max_asv} ASVs
  -n  <int>            Minimum number of reads (absolute count) per ASV per sample
  -ab <float>          Minimum relative abundance in % per sample (e.g. 0.1 = 0.1%)
  -h                   Show this help message

EXAMPLES:
  bash easy_net.sh -d ./data -s /opt/SparCC --auto
  bash easy_net.sh -d ./data -s /opt/SparCC -ab 0.1
  bash easy_net.sh -d ./data -s /opt/SparCC -n 5 -ab 0.05
  bash easy_net.sh -d ./data -s /opt/SparCC --auto -n 3

NOTES:
  - If both -n and -ab are provided, both filters are applied (AND logic).
  - --auto searches for the lowest absolute cutoff; -n is applied before counting.

EOF
}

# =============================================================================
# Argument parsing (supports --auto and long flags)
# =============================================================================
ARGS=("$@")
i=0
while [ $i -lt ${#ARGS[@]} ]; do
    arg="${ARGS[$i]}"
    case "$arg" in
        -d)
            i=$((i+1)); folder="${ARGS[$i]}" ;;
        -s)
            i=$((i+1)); SPARCC_PATH="${ARGS[$i]}" ;;
        --auto)
            auto_mode=true ;;
        -n)
            i=$((i+1)); min_seqs="${ARGS[$i]}" ;;
        -ab)
            i=$((i+1)); min_ab_pct="${ARGS[$i]}" ;;
        -h|--help)
            usage; exit 0 ;;
        *)
            log_error "Unknown option: $arg"
            usage; exit 1 ;;
    esac
    i=$((i+1))
done

# =============================================================================
# Input validation
# =============================================================================
if [ -z "$folder" ]; then
    log_error "Input folder is required. Use -d <folder>."
    usage
    exit 1
fi

if [ -z "$SPARCC_PATH" ]; then
    log_error "SparCC path is required. Use -s <sparcc_path>."
    usage
    exit 1
fi

if [ ! -d "$folder" ]; then
    log_error "Folder not found: $folder"
    exit 1
fi

shopt -s nullglob
txt_files=("${folder}"/*.txt)
shopt -u nullglob
if [ ${#txt_files[@]} -eq 0 ]; then
    log_error "No .txt files found in: $folder"
    exit 1
fi

if [ ! -d "$SPARCC_PATH" ]; then
    log_error "SparCC not found at: $SPARCC_PATH"
    exit 1
fi

# =============================================================================
# Function: count ASVs by absolute abundance cutoff (+ optional -n)
# =============================================================================
count_asvs_by_min_abundance() {
    local file=$1
    local min_abund=$2
    local min_n=${3:-0}

    awk -F'\t' -v min_abund="$min_abund" -v min_n="$min_n" '
    NR==1 { next }
    {
        keep=0
        for (i=2; i<=NF; i++) {
            val=$i+0
            if (val >= min_abund && (min_n==0 || val >= min_n)) {
                keep=1; break
            }
        }
        if (keep) count++
    }
    END { print count+0 }
    ' "$file"
}

# =============================================================================
# Function: count ASVs by relative abundance (%)
# =============================================================================
count_asvs_by_rel_abundance() {
    local file=$1
    local min_pct=$2
    local min_n=${3:-0}

    awk -F'\t' -v min_pct="$min_pct" -v min_n="$min_n" '
    NR==1 { next }
    {
        sum=0
        for (i=2; i<=NF; i++) sum += $i+0

        keep=0
        if (sum > 0) {
            for (i=2; i<=NF; i++) {
                val=$i+0
                rel=(val/sum)*100
                if (rel >= min_pct && (min_n==0 || val >= min_n)) {
                    keep=1; break
                }
            }
        }
        if (keep) count++
    }
    END { print count+0 }
    ' "$file"
}

# =============================================================================
# Function: filter file and write to filtered/
# =============================================================================
filter_file() {
    local file=$1
    local out=$2
    local mode=$3       # "abs" or "rel"
    local cutoff=$4
    local min_n=${5:-0}

    if [ "$mode" = "rel" ]; then
        awk -F'\t' -v min_pct="$cutoff" -v min_n="$min_n" '
        NR==1 { print; next }
        {
            sum=0
            for (i=2; i<=NF; i++) sum += $i+0
            keep=0
            if (sum > 0) {
                for (i=2; i<=NF; i++) {
                    val=$i+0
                    rel=(val/sum)*100
                    if (rel >= min_pct && (min_n==0 || val >= min_n)) {
                        keep=1; break
                    }
                }
            }
            if (keep) print
        }
        ' "$file" > "$out"
    else
        awk -F'\t' -v min_abund="$cutoff" -v min_n="$min_n" '
        NR==1 { print; next }
        {
            keep=0
            for (i=2; i<=NF; i++) {
                val=$i+0
                if (val >= min_abund && (min_n==0 || val >= min_n)) {
                    keep=1; break
                }
            }
            if (keep) print
        }
        ' "$file" > "$out"
    fi
}

# =============================================================================
# Pipeline start
# =============================================================================
echo ""
echo "============================================"
echo "     easy_net -- SparCC Network Pipeline    "
echo "============================================"
log_info "Input folder  : $folder"
log_info "SparCC path   : $SPARCC_PATH"
log_info ".txt files    : ${#txt_files[@]}"
log_info "Auto mode     : $auto_mode"
[ -n "$min_seqs" ]   && log_info "Min reads     : $min_seqs"
[ -n "$min_ab_pct" ] && log_info "Min rel abund : ${min_ab_pct}%"
echo ""

# =============================================================================
# Filtering logic
# =============================================================================
filter_mode=""
cutoff_value=""

if $auto_mode; then
    log_step "Auto mode: searching for best absolute abundance cutoff"
    log_info "Target: all tables with < ${max_asv} ASVs"
    echo ""

    best_cutoff=""
    for cutoff in 2 5 10 20 50 100 150 200 300 400 500 1000; do
        ok=true
        echo -ne "  Testing cutoff = ${cutoff} ... "
        for file in "${txt_files[@]}"; do
            asvs=$(count_asvs_by_min_abundance "$file" "$cutoff" "${min_seqs:-0}")
            echo -ne "$(basename $file): ${asvs} ASVs  "
            if [ "$asvs" -gt "$max_asv" ]; then
                ok=false
                break
            fi
        done
        if $ok; then
            echo " OK"
            best_cutoff="$cutoff"
            break
        else
            echo " exceeds ${max_asv}"
        fi
    done

    if [ -z "$best_cutoff" ]; then
        log_error "No cutoff reduced all tables to < ${max_asv} ASVs."
        log_warn  "Consider increasing max_asv in the script or reviewing your dataset."
        exit 1
    fi

    log_ok "Auto cutoff selected: ${best_cutoff} (absolute abundance)"
    filter_mode="abs"
    cutoff_value="$best_cutoff"

elif [ -n "$min_ab_pct" ]; then
    log_step "Filtering by minimum relative abundance: ${min_ab_pct}%"
    filter_mode="rel"
    cutoff_value="$min_ab_pct"

elif [ -n "$min_seqs" ]; then
    log_step "Filtering by minimum read count: ${min_seqs}"
    filter_mode="abs"
    cutoff_value="$min_seqs"
    min_seqs=""

else
    log_warn "No filter specified (--auto, -ab, or -n). Proceeding without filtering."
    filter_mode="none"
fi

# =============================================================================
# Apply filtering
# =============================================================================
filtered_folder="${folder}/filtered"
mkdir -p "$filtered_folder"

log_step "Applying filters -- output: ${filtered_folder}/"
echo ""

for file in "${txt_files[@]}"; do
    base=$(basename "$file")
    out="${filtered_folder}/${base}"

    if [ "$filter_mode" = "none" ]; then
        cp "$file" "$out"
        log_info "Copied (no filter): $base"
    else
        before=$(awk 'NR>1' "$file" | wc -l)
        filter_file "$file" "$out" "$filter_mode" "$cutoff_value" "${min_seqs:-0}"
        after=$(awk 'NR>1' "$out" | wc -l)
        log_ok "$base : ${before} -> ${after} ASVs"
    fi
done

echo ""
log_ok "Filtering done. Files saved in: ${filtered_folder}/"

# =============================================================================
# SparCC
# =============================================================================
log_step "Starting SparCC analyses"
echo ""

cd "$filtered_folder" || exit 1

shopt -s nullglob
filtered_files=(*.txt)
shopt -u nullglob

if [ ${#filtered_files[@]} -eq 0 ]; then
    log_error "No .txt files found in filtered/."
    exit 1
fi

total=${#filtered_files[@]}
count=0

for file in "${filtered_files[@]}"; do
    count=$((count+1))
    base_name="${file%.txt}"

    echo ""
    echo "--- [$count/$total] Processing: $file ---"

    net_dir="${base_name}_net"
    mkdir -p "${net_dir}/perm"
    mkdir -p "${net_dir}/pvalues"

    # SparCC correlations
    log_info "Computing SparCC correlations..."
    python "${SPARCC_PATH}/SparCC.py" "${file}" \
        --cor_file="${net_dir}/cor_sparcc.out"
    log_ok "Correlations saved: ${net_dir}/cor_sparcc.out"

    # Bootstraps
    log_info "Generating 100 permutations..."
    python "${SPARCC_PATH}/MakeBootstraps.py" "${file}" \
        -n 100 \
        -t permutation_#.txt \
        -p "${net_dir}/perm/"
    log_ok "Permutations saved in: ${net_dir}/perm/"

    # SparCC on permutations
    log_info "Running SparCC on each permutation (0-99)..."
    for f in $(seq 0 99); do
        python "${SPARCC_PATH}/SparCC.py" \
            "${net_dir}/perm/permutation_${f}.txt" \
            -i 100 \
            --cor_file="${net_dir}/pvalues/perm_cor${f}.txt"
    done
    log_ok "Permutation SparCC done."

    # PseudoPvals
    log_info "Computing pseudo p-values (two-sided)..."
    python "${SPARCC_PATH}/PseudoPvals.py" \
        "${net_dir}/cor_sparcc.out" \
        "${net_dir}/pvalues/perm_cor#.txt" \
        100 \
        -o "${net_dir}/pvals_two_sided.txt" \
        -t two_sided
    log_ok "P-values saved: ${net_dir}/pvals_two_sided.txt"

    # Significant pairs
    log_info "Extracting significant pairs..."
    cd "${net_dir}" || exit 1
    python "${SPARCC_PATH}/get_significant_pairs.py"
    log_ok "Significant pairs extracted in: ${net_dir}/"
    cd ..

done

echo ""
echo "============================================"
log_ok "Pipeline completed successfully!"
echo "============================================"
echo ""
