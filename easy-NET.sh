#!/bin/bash

# =============================================================================
# easy_net.sh — SparCC Pipeline with sum-of-replicates filtering
# =============================================================================
# Usage:
#   bash easy_net.sh -d <folder> [--abs | --abr] [-n <sum>] [-ab <pct>]
#
# Automatic filtering modes (choose one):
#   --abs    Automatic sweep by absolute read sum (1 to 10000, step 1)
#            Selects the lowest cutoff that leaves <= 999 ASVs in all treatments
#   --abr    Automatic sweep by relative abundance (0.0001% to 1%, step 0.0001%)
#            Selects the lowest cutoff % that leaves <= 999 ASVs in all treatments
#
# Manual modes (no sweep):
#   -n  <int>      Minimum absolute read sum per ASV (e.g., -n 150)
#   -ab <float>    Minimum relative abundance in % per ASV (e.g., -ab 0.05)
#
# Other:
#   -d  <folder>   Folder containing .txt files (required)
#   -h             Display this help
#
# Filtering logic:
#   Each ASV is evaluated by the SUM of all its replicates (columns).
#   The same cutoff is applied to all treatments, ensuring comparability.
#   At the end, the script reports how many ASVs (%) were retained per treatment.
#
# Examples:
#   bash easy_net.sh -d ./data --abs
#   bash easy_net.sh -d ./data --abr
#   bash easy_net.sh -d ./data -n 150
#   bash easy_net.sh -d ./data -ab 0.02
#
# Notes:
#   - SparCC should be in: /Dados/bioinformatic_tools/SparCC
# =============================================================================

SPARCC_PATH="/home/thierry_bioinfo/works/bioinformatics_tools/SparCC/" ####CHANGE HERE YOUR SPARCC PATH

log_info()    { echo "[INFO]    $*"; }
log_ok()      { echo "[OK]      $*"; }
log_warn()    { echo "[WARNING] $*"; }
log_error()   { echo "[ERROR]   $*"; }
log_step()    { echo ""; echo ">>> $*"; }

# --- Defaults ---
folder=""
max_asv=999
mode_abs=false
mode_abr=false
min_seqs=""
min_ab_pct=""

# =============================================================================
# Help
# =============================================================================
usage() {
cat <<EOF

easy_net.sh — SparCC Pipeline with sum-of-replicates filtering

USAGE:
  bash easy_net.sh -d <folder> [--abs | --abr | -n <sum> | -ab <pct>]

REQUIRED OPTIONS:
  -d  <folder>    Folder containing .txt files

AUTOMATIC FILTERING MODES (sweep, choose one):
  --abs          Absolute sum sweep: 1 to 10000 (step 1)
                 Selects the lowest cutoff with <= ${max_asv} ASVs in all treatments
  --abr          Relative abundance sweep: 0.0001% to 1% (step 0.0001%)
                 Selects the lowest cutoff % with <= ${max_asv} ASVs in all treatments

MANUAL MODES (no sweep):
  -n  <int>      Minimum absolute read sum per ASV (e.g., -n 150)
  -ab <float>    Minimum relative abundance in % per ASV (e.g., -ab 0.02)

OTHER:
  -h             Display this help

EXAMPLES:
  bash easy_net.sh -d ./data --abs
  bash easy_net.sh -d ./data --abr
  bash easy_net.sh -d ./data -n 150
  bash easy_net.sh -d ./data -ab 0.02

NOTES:
  - The filter is applied over the SUM of all replicates of each ASV.
  - The same cutoff is used in all treatments (comparability guaranteed).
  - For --abr: relative abundance = ASV_sum / total_treatment_sum * 100
  - SparCC expected in: ${SPARCC_PATH}

EOF
}

# =============================================================================
# Argument parsing
# =============================================================================
ARGS=("$@")
i=0
while [ $i -lt ${#ARGS[@]} ]; do
    arg="${ARGS[$i]}"
    case "$arg" in
        -d)
            i=$((i+1)); folder="${ARGS[$i]}" ;;
        --abs)
            mode_abs=true ;;
        --abr)
            mode_abr=true ;;
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

# Validate: --abs and --abr are mutually exclusive
if $mode_abs && $mode_abr; then
    log_error "--abs and --abr are mutually exclusive. Choose only one."
    usage
    exit 1
fi

# =============================================================================
# Initial validations
# =============================================================================
if [ -z "$folder" ]; then
    log_error "Input folder is required. Use -d <folder>."
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
    log_error "SparCC not found in: $SPARCC_PATH"
    exit 1
fi

# =============================================================================
# Automatic patch in get_significant_pairs.py (adds delimiter='\t')
# Necessary because genfromtxt without delimiter fails with empty fields between tabs
# =============================================================================
patch_get_significant_pairs() {
    local script="${SPARCC_PATH}/get_significant_pairs.py"

    if grep -q "delimiter='\\\\t'" "$script" 2>/dev/null; then
        log_info "get_significant_pairs.py is already patched (delimiter='\\t'). Skipping."
        return 0
    fi

    # Make backup on first run
    if [ ! -f "${script}.bak" ]; then
        cp "$script" "${script}.bak"
        log_info "Backup created: ${script}.bak"
    fi

    # Replace the genfromtxt line by adding delimiter='\t'
    sed -i "s/genfromtxt(path, names=True/genfromtxt(path, delimiter='\\t', names=True/" "$script"

    if grep -q "delimiter='\\\\t'" "$script"; then
        log_ok "Patch applied to get_significant_pairs.py (delimiter='\\t' added)"
    else
        log_warn "Automatic patch failed — manually check line 13 of get_significant_pairs.py"
    fi
}

# =============================================================================
# Function: check if cor_sparcc.out has real data (not just tabs)
# =============================================================================
check_cor_file() {
    local cor_file=$1
    # If the file has no digits, it's empty
    if ! grep -qP '\d' "$cor_file" 2>/dev/null; then
        return 1
    fi
    return 0
}

# =============================================================================
# Function: count ASVs retained by minimum absolute sum filter
# An ASV is retained if the SUM of all its replicates >= min_sum
# =============================================================================
count_asvs_by_sum() {
    local file=$1
    local min_suma=$2

    awk -F'\t' -v min_suma="$min_suma" '
    NR==1 { next }
    {
        suma=0
        for (i=2; i<=NF; i++) suma += $i+0
        if (suma >= min_suma) count++
    }
    END { print count+0 }
    ' "$file"
}

# =============================================================================
# Function: count ASVs retained by relative sum filter (% of total sum)
# The ASV sum must represent >= min_pct % of the treatment's total sum
# =============================================================================
count_asvs_by_sum_rel() {
    local file=$1
    local min_pct=$2

    awk -F'\t' -v min_pct="$min_pct" '
    NR==1 { next }
    {
        suma=0
        for (i=2; i<=NF; i++) suma += $i+0
        total += suma
        sumas[NR] = suma
    }
    END {
        for (r in sumas) {
            if (total > 0 && (sumas[r]/total)*100 >= min_pct) count++
        }
        print count+0
    }
    ' "$file"
}

# =============================================================================
# Function: filter file by sum and save to filtered/ folder
# Modes: abs (absolute sum), rel (individual %), rel_acum (cumulative %)
# =============================================================================
filter_file() {
    local file=$1
    local out=$2
    local mode=$3
    local cutoff=$4

    if [ "$mode" = "rel" ]; then
        # Keep ASVs whose individual % >= cutoff
        awk -F'\t' -v min_pct="$cutoff" '
        NR==1 { header=$0; next }
        {
            suma=0
            for (i=2; i<=NF; i++) suma += $i+0
            total += suma
            sumas[NR] = suma
            lines[NR] = $0
        }
        END {
            print header
            for (r=2; r<=NR; r++) {
                if (r in sumas && total > 0 && (sumas[r]/total)*100 >= min_pct)
                    print lines[r]
            }
        }
        ' "$file" > "$out"

    elif [ "$mode" = "rel_acum" ]; then
        # Keep top N ASVs (by decreasing %) that cover >= cutoff% cumulative
        # Uses Python to sort by % and select the correct lines
        python3 - "$file" "$out" "$cutoff" << 'PYEOF'
import sys

fpath, outpath, cutoff = sys.argv[1], sys.argv[2], float(sys.argv[3])

rows = []
with open(fpath) as fh:
    header = fh.readline()
    for line in fh:
        parts = line.rstrip('\r\n').split('\t')
        s = sum(float(x) for x in parts[1:] if x.strip())
        rows.append((s, line))

suma_total = sum(r[0] for r in rows)

# Sort by decreasing sum, accumulate %, stop when reaching cutoff
rows_sorted = sorted(rows, key=lambda x: x[0], reverse=True)
keep = set()
acum = 0.0
for i, (s, line) in enumerate(rows_sorted):
    keep.add(i)
    acum += (s / suma_total * 100) if suma_total > 0 else 0
    if acum >= cutoff:
        break

# Write maintaining the original file order
with open(outpath, 'w') as fout:
    fout.write(header)
    for i, (s, line) in enumerate(rows_sorted):
        if i in keep:
            fout.write(line)
PYEOF

    else
        # abs mode: keep ASVs with sum >= cutoff
        awk -F'\t' -v min_suma="$cutoff" '
        NR==1 { print; next }
        {
            suma=0
            for (i=2; i<=NF; i++) suma += $i+0
            if (suma >= min_suma) print
        }
        ' "$file" > "$out"
    fi
}

# =============================================================================
# Function: statistical profile of a file (ASVs, total sum, min, max per ASV)
# =============================================================================
print_treatment_profile() {
    local file=$1
    local label=$2
    awk -F'\t' -v label="$label" '
    NR==1 { next }
    {
        suma=0
        for (i=2; i<=NF; i++) suma += $i+0
        total += suma
        if (NR==2 || suma < min_s) min_s = suma
        if (NR==2 || suma > max_s) max_s = suma
        count++
    }
    END {
        printf "  %-30s ASVs: %4d  total_sum: %10d  min_sum: %6d  max_sum: %8d\n",
               label, count, total, min_s, max_s
    }
    ' "$file"
}

# =============================================================================
# Pipeline start — with duplicated logging to file
# =============================================================================
log_file="${folder}/easy_net_$(date '+%Y%m%d_%H%M%S').log"

# Redirect stdout+stderr to terminal AND log simultaneously
exec > >(tee -a "$log_file") 2>&1

echo ""
echo "============================================"
echo "   SparCC Pipeline — Sum-based Filtering    "
echo "============================================"
echo "  Log saved in: $log_file"
echo "  Date/time   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "--------------------------------------------"
log_info "Input folder     : $folder"
log_info ".txt files       : ${#txt_files[@]}"
$mode_abs && log_info "Mode             : --abs (absolute sum sweep 1-10000)"
$mode_abr && log_info "Mode             : --abr (relative abundance sweep 0.0001%-1%)"
[ -n "$min_seqs" ]   && log_info "Mode             : -n manual  | Minimum sum: $min_seqs reads"
[ -n "$min_ab_pct" ] && log_info "Mode             : -ab manual | Min. rel. abundance: ${min_ab_pct}%"
log_info "Criterion        : sum of all replicates per ASV"
echo ""

# --- Initial treatment profiles ---
echo "--------------------------------------------"
echo "  Treatment profiles (raw data):"
echo "--------------------------------------------"
for file in "${txt_files[@]}"; do
    print_treatment_profile "$file" "$(basename "$file")"
done
echo "--------------------------------------------"
echo ""

# Apply patch to get_significant_pairs.py before anything else
log_step "Checking/patching get_significant_pairs.py"
patch_get_significant_pairs

# =============================================================================
# FILTERING LOGIC
# =============================================================================
filter_mode=""
cutoff_value=""

if $mode_abs; then
    log_step "Mode --abs: sweeping absolute sum from 1 to 10000 (step 1)"
    log_info "Goal: all treatments with <= ${max_asv} ASVs"
    echo ""

    best_cutoff=$(python3 - "${txt_files[@]}" << PYEOF
import sys, bisect

files = sys.argv[1:]
max_asv = ${max_asv}

# Pre-calculate absolute sums for each file (single read)
sumas = {}
for fpath in files:
    vals = []
    with open(fpath) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            parts = line.rstrip('\r\n').split('\t')
            s = sum(float(x) for x in parts[1:] if x.strip())
            vals.append(int(s))
    vals.sort()
    sumas[fpath] = vals
    suma_total = sum(vals)
    sys.stderr.write(
        f"  {fpath.split('/')[-1]}: {len(vals)} ASVs  "
        f"total_sum={suma_total:,}  min={vals[0] if vals else 0}  max={vals[-1] if vals else 0}\n"
    )

# Sweep 1-10000 with binary search
best = None
for cut in range(1, 10001):
    ok = True
    for vals in sumas.values():
        idx = bisect.bisect_left(vals, cut)
        if (len(vals) - idx) > max_asv:
            ok = False
            break
    if ok:
        best = cut
        break

print(-1 if best is None else best)
PYEOF
)

    if [ "$best_cutoff" = "-1" ] || [ -z "$best_cutoff" ]; then
        log_error "No cutoff between 1 and 10000 reduced all treatments to <= ${max_asv} ASVs."
        log_warn "Consider using -n with value > 10000 for higher cutoffs."
        exit 1
    fi

    log_ok "Cutoff --abs selected: sum >= ${best_cutoff} reads (same cutoff for all treatments)"
    filter_mode="abs"
    cutoff_value="$best_cutoff"

elif $mode_abr; then
    log_step "Mode --abr: sweep by cumulative relative abundance (0.01% to 100%, step 0.01%)"
    log_info "Logic: ASVs ordered from most to least abundant, summing % until reaching cutoff"
    log_info "Goal: smallest set of ASVs (all treatments <= ${max_asv}) covering X% of reads"
    echo ""

    best_cutoff=$(python3 - "${txt_files[@]}" << PYEOF
import sys, bisect

files = sys.argv[1:]
max_asv = ${max_asv}

# Pre-calculate cumulative abundance curve per treatment
# ASVs ordered from most to least abundant; acum[n] = cumulative % of top n+1 ASVs
curvas = {}
for fpath in files:
    vals = []
    with open(fpath) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            parts = line.rstrip('\r\n').split('\t')
            s = sum(float(x) for x in parts[1:] if x.strip())
            vals.append(s)
    suma_total = sum(vals)
    if suma_total == 0:
        sys.stderr.write(f"  WARNING: {fpath.split('/')[-1]} has total sum = 0\n")
        curvas[fpath] = []
        continue
    pcts = sorted((v / suma_total * 100 for v in vals), reverse=True)
    acum = []
    s = 0.0
    for p in pcts:
        s += p
        acum.append(round(s, 8))
    curvas[fpath] = acum
    sys.stderr.write(
        f"  {fpath.split('/')[-1]}: {len(acum)} ASVs  "
        f"total_sum={int(suma_total):,}  "
        f"min_abr={pcts[-1]:.6f}%  max_abr={pcts[0]:.4f}%\n"
    )

# Sweep from highest to lowest cutoff (0.01% to 100%, step 0.01%)
# Find the most restrictive cutoff (highest %) where N_retained <= max_asv in all
# N_retained = bisect_left(acum, cutoff) + 1  (how many ASVs cover >= cutoff%)
best = None
step = 0.01
cut = 100.0
while cut >= 0.0 - 1e-9:
    ok = True
    for acum in curvas.values():
        idx = bisect.bisect_left(acum, cut - 1e-9)
        n_mantidos = idx + 1
        if n_mantidos > max_asv:
            ok = False
            break
    if ok:
        best = round(cut, 4)
        break
    cut = round(cut - step, 8)

if best is None:
    print(-1)
else:
    # Print result per treatment to stderr (goes to log)
    sys.stderr.write(f"\n  Cutoff selected: {best}% cumulative\n")
    for fpath, acum in curvas.items():
        idx = bisect.bisect_left(acum, best - 1e-9)
        n = idx + 1
        sys.stderr.write(
            f"  {fpath.split('/')[-1]}: {n} ASVs retained  "
            f"({acum[idx]:.4f}% actual cumulative)\n"
        )
    print(best)
PYEOF
)

    if [ "$best_cutoff" = "-1" ] || [ -z "$best_cutoff" ]; then
        log_error "No cutoff between 0.01% and 100% left all treatments with <= ${max_asv} ASVs."
        log_warn "Consider increasing max_asv in the script."
        exit 1
    fi

    log_ok "Cutoff --abr selected: ${best_cutoff}% cumulative abundance (same cutoff for all)"
    filter_mode="rel_acum"
    cutoff_value="$best_cutoff"

elif [ -n "$min_ab_pct" ]; then
    log_step "Mode -ab manual: minimum relative abundance = ${min_ab_pct}% of treatment total sum"
    filter_mode="rel"
    cutoff_value="$min_ab_pct"

elif [ -n "$min_seqs" ]; then
    log_step "Mode -n manual: minimum read sum = ${min_seqs}"
    filter_mode="abs"
    cutoff_value="$min_seqs"

else
    log_warn "No filtering mode specified (--abs, --abr, -n or -ab). Proceeding without filtering."
    filter_mode="none"
fi

# =============================================================================
# APPLY FILTERING + REPORT % RETAINED
# =============================================================================
filtered_folder="${folder}/filtered"
mkdir -p "$filtered_folder"

log_step "Applying filters and saving to: ${filtered_folder}/"
echo ""

# Report header
echo "------------------------------------------------------------"
printf "  %-30s %6s %6s %8s\n" "File" "Before" "After" "Retained"
echo "------------------------------------------------------------"

for file in "${txt_files[@]}"; do
    base=$(basename "$file")
    out="${filtered_folder}/${base}"

    if [ "$filter_mode" = "none" ]; then
        cp "$file" "$out"
        before=$(awk 'NR>1' "$file" | wc -l)
        printf "  %-30s %6d %6d %7s\n" "$base" "$before" "$before" "100.0%"
    else
        before=$(awk 'NR>1' "$file" | wc -l)
        filter_file "$file" "$out" "$filter_mode" "$cutoff_value"
        after=$(awk 'NR>1' "$out" | wc -l)
        pct=$(awk -v a="$after" -v b="$before" 'BEGIN{printf "%.1f%%", (b>0)?(a/b)*100:0}')
        printf "  %-30s %6d %6d %8s\n" "$base" "$before" "$after" "$pct"
    fi
done

echo "------------------------------------------------------------"
echo ""
log_ok "Filtering completed. Files in: ${filtered_folder}/"

# =============================================================================
# SPARCC
# =============================================================================
log_step "Starting SparCC analyses"
echo ""

cd "$filtered_folder" || exit 1

shopt -s nullglob
filtered_files=(*.txt)
shopt -u nullglob

if [ ${#filtered_files[@]} -eq 0 ]; then
    log_error "No .txt files in filtered/ folder."
    exit 1
fi

total=${#filtered_files[@]}
count=0
failed_files=()

for file in "${filtered_files[@]}"; do
    count=$((count+1))
    base_name="${file%.txt}"

    echo ""
    echo "--- [$count/$total] Processing: $file ---"

    net_dir="${base_name}_net"
    mkdir -p "${net_dir}/perm"
    mkdir -p "${net_dir}/pvalues"

    # --- Main SparCC ---
    log_info "Calculating SparCC correlations..."
    python "${SPARCC_PATH}/SparCC.py" "${file}" \
        --cor_file="${net_dir}/cor_sparcc.out"
    if [ $? -ne 0 ]; then
        log_error "SparCC.py failed on: $file — skipping."
        failed_files+=("$file")
        continue
    fi

    # Check if output file has real data
    if ! check_cor_file "${net_dir}/cor_sparcc.out"; then
        log_error "cor_sparcc.out is empty (no numeric values) for: $file"
        log_warn "Possible cause: incorrect input table format (separator, header, etc.)"
        log_warn "Check: head -2 ${file} | cat -A"
        failed_files+=("$file")
        continue
    fi
    log_ok "Correlations calculated: ${net_dir}/cor_sparcc.out"

    # --- Bootstraps ---
    log_info "Generating 100 permutations..."
    python "${SPARCC_PATH}/MakeBootstraps.py" "${file}" \
        -n 100 \
        -t permutation_#.txt \
        -p "${net_dir}/perm/"
    if [ $? -ne 0 ]; then
        log_error "MakeBootstraps.py failed on: $file — skipping."
        failed_files+=("$file")
        continue
    fi
    log_ok "Permutations generated in: ${net_dir}/perm/"

    # --- SparCC on permutations ---
    log_info "Calculating SparCC for each permutation (0–99)..."
    perm_failed=0
    for f in $(seq 0 99); do
        python "${SPARCC_PATH}/SparCC.py" \
            "${net_dir}/perm/permutation_${f}.txt" \
            -i 100 \
            --cor_file="${net_dir}/pvalues/perm_cor${f}.txt"
        if [ $? -ne 0 ]; then
            log_warn "SparCC failed on permutation ${f} — continuing."
            perm_failed=$((perm_failed+1))
        fi
    done
    if [ $perm_failed -gt 0 ]; then
        log_warn "${perm_failed}/100 permutations failed for: $file"
    fi
    log_ok "Permutation SparCC completed."

    # --- PseudoPvals ---
    log_info "Calculating pseudo p-values (two-sided)..."
    python "${SPARCC_PATH}/PseudoPvals.py" \
        "${net_dir}/cor_sparcc.out" \
        "${net_dir}/pvalues/perm_cor#.txt" \
        100 \
        -o "${net_dir}/pvals_two_sided.txt" \
        -t two_sided
    if [ $? -ne 0 ]; then
        log_error "PseudoPvals.py failed on: $file — skipping."
        failed_files+=("$file")
        continue
    fi
    log_ok "P-values: ${net_dir}/pvals_two_sided.txt"

    # --- Significant pairs ---
    log_info "Extracting significant pairs..."
    cd "${net_dir}" || exit 1
    python "${SPARCC_PATH}/get_significant_pairs.py"
    if [ $? -ne 0 ]; then
        log_error "get_significant_pairs.py failed in: ${net_dir}/"
        log_warn "Check the format of cor_sparcc.out and pvals_two_sided.txt"
        cd ..
        failed_files+=("$file")
        continue
    fi
    log_ok "Significant pairs extracted in: ${net_dir}/"
    cd ..

done

# =============================================================================
# Final summary
# =============================================================================
echo ""
echo "============================================"
if [ ${#failed_files[@]} -eq 0 ]; then
    log_ok "Pipeline completed successfully! ($total/$total files)"
else
    log_warn "Pipeline completed with errors."
    log_warn "Files that failed (${#failed_files[@]}/${total}):"
    for f in "${failed_files[@]}"; do
        echo "    - $f"
    done
fi

# Treatment profiles after filtering
echo ""
echo "--------------------------------------------"
echo "  Treatment profiles (after filtering):"
echo "--------------------------------------------"
for file in "${filtered_folder}"/*.txt; do
    print_treatment_profile "$file" "$(basename "$file")"
done
echo "--------------------------------------------"
echo ""
echo "  Complete log saved in:"
echo "  $log_file"
echo "============================================"
echo ""
