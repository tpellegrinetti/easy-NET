#!/bin/bash

# =============================================================================
# easy_net.sh — Pipeline SparCC com filtragem por soma de repetições
# =============================================================================
# Uso:
#   bash easy_net.sh -d <pasta> [--abs | --abr] [-n <soma>] [-ab <pct>]
#
# Modos de filtragem automática (escolha um):
#   --abs    Varredura automática por soma absoluta de reads (1 a 10000, passo 1)
#            Seleciona o menor cutoff que deixa <= 999 ASVs em todos os tratamentos
#   --abr    Varredura automática por abundância relativa (0.0001% a 1%, passo 0.0001%)
#            Seleciona o menor cutoff % que deixa <= 999 ASVs em todos os tratamentos
#
# Modos manuais (sem varredura):
#   -n  <int>      Soma absoluta mínima de reads por ASV (ex: -n 150)
#   -ab <float>    Abundância relativa mínima em % por ASV (ex: -ab 0.05)
#
# Outros:
#   -d  <pasta>    Pasta contendo os arquivos .txt (obrigatório)
#   -h             Exibe esta ajuda
#
# Lógica de filtragem:
#   Cada ASV é avaliado pela SOMA de todas as suas repetições (colunas).
#   O mesmo cutoff é aplicado a todos os tratamentos, garantindo comparabilidade.
#   Ao final, o script informa quantos ASVs (%) foram mantidos por tratamento.
#
# Exemplos:
#   bash easy_net.sh -d ./dados --abs
#   bash easy_net.sh -d ./dados --abr
#   bash easy_net.sh -d ./dados -n 150
#   bash easy_net.sh -d ./dados -ab 0.02
#
# Notas:
#   - O SparCC deve estar em: /Dados/bioinformatic_tools/SparCC
# =============================================================================

SPARCC_PATH="/home/thierry_bioinfo/works/bioinformatics_tools/SparCC/"

log_info()    { echo "[INFO]  $*"; }
log_ok()      { echo "[OK]    $*"; }
log_warn()    { echo "[AVISO] $*"; }
log_error()   { echo "[ERRO]  $*"; }
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

easy_net.sh — Pipeline SparCC com filtragem por soma de repetições

USO:
  bash easy_net.sh -d <pasta> [--abs | --abr | -n <soma> | -ab <pct>]

OPCOES OBRIGATORIAS:
  -d  <pasta>    Pasta contendo os arquivos .txt

MODOS DE FILTRAGEM AUTOMATICA (varredura, escolha um):
  --abs          Varredura de soma absoluta: 1 a 10000 (passo 1)
                 Seleciona o menor cutoff com <= ${max_asv} ASVs em todos os tratamentos
  --abr          Varredura de abundancia relativa: 0.0001% a 1% (passo 0.0001%)
                 Seleciona o menor cutoff % com <= ${max_asv} ASVs em todos os tratamentos

MODOS MANUAIS (sem varredura):
  -n  <int>      Soma absoluta minima de reads por ASV (ex: -n 150)
  -ab <float>    Abundancia relativa minima em % por ASV (ex: -ab 0.02)

OUTROS:
  -h             Exibe esta ajuda

EXEMPLOS:
  bash easy_net.sh -d ./dados --abs
  bash easy_net.sh -d ./dados --abr
  bash easy_net.sh -d ./dados -n 150
  bash easy_net.sh -d ./dados -ab 0.02

NOTAS:
  - O filtro e aplicado sobre a SOMA de todas as repeticoes de cada ASV.
  - O mesmo cutoff e usado em todos os tratamentos (comparabilidade garantida).
  - Para --abr: abundancia relativa = soma_ASV / soma_total_tratamento * 100
  - SparCC esperado em: ${SPARCC_PATH}

EOF
}

# =============================================================================
# Parse de argumentos
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
            log_error "Opção desconhecida: $arg"
            usage; exit 1 ;;
    esac
    i=$((i+1))
done

# Validar: --abs e --abr são mutuamente exclusivos
if $mode_abs && $mode_abr; then
    log_error "--abs e --abr são mutuamente exclusivos. Escolha apenas um."
    usage
    exit 1
fi

# =============================================================================
# Validações iniciais
# =============================================================================
if [ -z "$folder" ]; then
    log_error "A pasta de entrada é obrigatória. Use -d <pasta>."
    usage
    exit 1
fi

if [ ! -d "$folder" ]; then
    log_error "Pasta não encontrada: $folder"
    exit 1
fi

shopt -s nullglob
txt_files=("${folder}"/*.txt)
shopt -u nullglob
if [ ${#txt_files[@]} -eq 0 ]; then
    log_error "Nenhum arquivo .txt encontrado em: $folder"
    exit 1
fi

if [ ! -d "$SPARCC_PATH" ]; then
    log_error "SparCC não encontrado em: $SPARCC_PATH"
    exit 1
fi

# =============================================================================
# Patch automático em get_significant_pairs.py (adiciona delimiter='\t')
# Necessário porque genfromtxt sem delimiter falha com campos vazios entre tabs
# =============================================================================
patch_get_significant_pairs() {
    local script="${SPARCC_PATH}/get_significant_pairs.py"

    if grep -q "delimiter='\\\\t'" "$script" 2>/dev/null; then
        log_info "get_significant_pairs.py já está patcheado (delimiter='\\t'). Pulando."
        return 0
    fi

    # Faz backup na primeira vez
    if [ ! -f "${script}.bak" ]; then
        cp "$script" "${script}.bak"
        log_info "Backup criado: ${script}.bak"
    fi

    # Substitui a linha do genfromtxt adicionando delimiter='\t'
    sed -i "s/genfromtxt(path, names=True/genfromtxt(path, delimiter='\\t', names=True/" "$script"

    if grep -q "delimiter='\\\\t'" "$script"; then
        log_ok "Patch aplicado em get_significant_pairs.py (delimiter='\\t' adicionado)"
    else
        log_warn "Patch automático falhou — verifique manualmente a linha 13 de get_significant_pairs.py"
    fi
}

# =============================================================================
# Função: checar se cor_sparcc.out tem dados reais (não só tabs)
# =============================================================================
check_cor_file() {
    local cor_file=$1
    # Se o arquivo não tem nenhum dígito, está vazio
    if ! grep -qP '\d' "$cor_file" 2>/dev/null; then
        return 1
    fi
    return 0
}

# =============================================================================
# Função: contar ASVs mantidos pelo filtro de soma absoluta mínima
# Um ASV é mantido se a SOMA de todas as suas repetições >= min_soma
# =============================================================================
count_asvs_by_sum() {
    local file=$1
    local min_soma=$2

    awk -F'\t' -v min_soma="$min_soma" '
    NR==1 { next }
    {
        soma=0
        for (i=2; i<=NF; i++) soma += $i+0
        if (soma >= min_soma) count++
    }
    END { print count+0 }
    ' "$file"
}

# =============================================================================
# Função: contar ASVs mantidos pelo filtro de soma relativa (% da soma total)
# A soma do ASV deve representar >= min_pct % da soma total do tratamento
# =============================================================================
count_asvs_by_sum_rel() {
    local file=$1
    local min_pct=$2

    awk -F'\t' -v min_pct="$min_pct" '
    NR==1 { next }
    {
        soma=0
        for (i=2; i<=NF; i++) soma += $i+0
        total += soma
        somas[NR] = soma
    }
    END {
        for (r in somas) {
            if (total > 0 && (somas[r]/total)*100 >= min_pct) count++
        }
        print count+0
    }
    ' "$file"
}

# =============================================================================
# Função: filtrar arquivo por soma e gravar na pasta filtered/
# Modos: abs (soma absoluta), rel (% individual), rel_acum (% acumulada)
# =============================================================================
filter_file() {
    local file=$1
    local out=$2
    local mode=$3
    local cutoff=$4

    if [ "$mode" = "rel" ]; then
        # Manter ASVs cuja % individual >= cutoff
        awk -F'\t' -v min_pct="$cutoff" '
        NR==1 { header=$0; next }
        {
            soma=0
            for (i=2; i<=NF; i++) soma += $i+0
            total += soma
            somas[NR] = soma
            lines[NR] = $0
        }
        END {
            print header
            for (r=2; r<=NR; r++) {
                if (r in somas && total > 0 && (somas[r]/total)*100 >= min_pct)
                    print lines[r]
            }
        }
        ' "$file" > "$out"

    elif [ "$mode" = "rel_acum" ]; then
        # Manter os top N ASVs (por % decrescente) que cobrem >= cutoff% acumulado
        # Usa Python para ordenar por % e selecionar as linhas corretas
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

soma_total = sum(r[0] for r in rows)

# Ordenar por soma decrescente, acumular %, parar ao atingir cutoff
rows_sorted = sorted(rows, key=lambda x: x[0], reverse=True)
keep = set()
acum = 0.0
for i, (s, line) in enumerate(rows_sorted):
    keep.add(i)
    acum += (s / soma_total * 100) if soma_total > 0 else 0
    if acum >= cutoff:
        break

# Gravar mantendo a ordem original do arquivo
with open(outpath, 'w') as fout:
    fout.write(header)
    for i, (s, line) in enumerate(rows_sorted):
        if i in keep:
            fout.write(line)
PYEOF

    else
        # Modo abs: manter ASVs com soma >= cutoff
        awk -F'\t' -v min_soma="$cutoff" '
        NR==1 { print; next }
        {
            soma=0
            for (i=2; i<=NF; i++) soma += $i+0
            if (soma >= min_soma) print
        }
        ' "$file" > "$out"
    fi
}

# =============================================================================
# Função: perfil estatístico de um arquivo (ASVs, soma total, min, max por ASV)
# =============================================================================
print_treatment_profile() {
    local file=$1
    local label=$2
    awk -F'\t' -v label="$label" '
    NR==1 { next }
    {
        soma=0
        for (i=2; i<=NF; i++) soma += $i+0
        total += soma
        if (NR==2 || soma < min_s) min_s = soma
        if (NR==2 || soma > max_s) max_s = soma
        count++
    }
    END {
        printf "  %-30s ASVs: %4d  soma_total: %10d  soma_min: %6d  soma_max: %8d\n",
               label, count, total, min_s, max_s
    }
    ' "$file"
}

# =============================================================================
# Início do pipeline — com log duplicado em arquivo
# =============================================================================
log_file="${folder}/easy_net_$(date '+%Y%m%d_%H%M%S').log"

# Redireciona stdout+stderr para terminal E para o log simultaneamente
exec > >(tee -a "$log_file") 2>&1

echo ""
echo "============================================"
echo "   Pipeline SparCC — Filtragem por Soma    "
echo "============================================"
echo "  Log salvo em: $log_file"
echo "  Data/hora   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "--------------------------------------------"
log_info "Pasta de entrada : $folder"
log_info "Arquivos .txt    : ${#txt_files[@]}"
$mode_abs && log_info "Modo             : --abs (varredura soma absoluta 1-10000)"
$mode_abr && log_info "Modo             : --abr (varredura abundância relativa 0.0001%-1%)"
[ -n "$min_seqs" ]   && log_info "Modo             : -n manual  | Soma mínima: $min_seqs reads"
[ -n "$min_ab_pct" ] && log_info "Modo             : -ab manual | Abund. rel. mínima: ${min_ab_pct}%"
log_info "Critério         : soma de todas as repetições por ASV"
echo ""

# --- Perfil inicial dos tratamentos ---
echo "--------------------------------------------"
echo "  Perfil dos tratamentos (dados brutos):"
echo "--------------------------------------------"
for file in "${txt_files[@]}"; do
    print_treatment_profile "$file" "$(basename "$file")"
done
echo "--------------------------------------------"
echo ""

# Aplica patch no get_significant_pairs.py antes de qualquer coisa
log_step "Verificando/patcheando get_significant_pairs.py"
patch_get_significant_pairs

# =============================================================================
# LÓGICA DE FILTRAGEM
# =============================================================================
filter_mode=""
cutoff_value=""

if $mode_abs; then
    log_step "Modo --abs: varrendo soma absoluta de 1 a 10000 (passo 1)"
    log_info "Meta: todos os tratamentos com <= ${max_asv} ASVs"
    echo ""

    best_cutoff=$(python3 - "${txt_files[@]}" << PYEOF
import sys, bisect

files = sys.argv[1:]
max_asv = ${max_asv}

# Pré-calcular somas absolutas de cada arquivo (leitura única)
somas = {}
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
    somas[fpath] = vals
    soma_total = sum(vals)
    sys.stderr.write(
        f"  {fpath.split('/')[-1]}: {len(vals)} ASVs  "
        f"soma_total={soma_total:,}  min={vals[0] if vals else 0}  max={vals[-1] if vals else 0}\n"
    )

# Varredura 1-10000 com busca binária
best = None
for cut in range(1, 10001):
    ok = True
    for vals in somas.values():
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
        log_error "Nenhum cutoff entre 1 e 10000 reduziu todos os tratamentos para <= ${max_asv} ASVs."
        log_warn "Considere usar -n com valor > 10000 para cutoffs maiores."
        exit 1
    fi

    log_ok "Cutoff --abs selecionado: soma >= ${best_cutoff} reads (mesmo corte para todos os tratamentos)"
    filter_mode="abs"
    cutoff_value="$best_cutoff"

elif $mode_abr; then
    log_step "Modo --abr: varredura por abundância relativa acumulada (0.01% a 100%, passo 0.01%)"
    log_info "Lógica: ASVs ordenados do mais ao menos abundante, somando % até atingir o cutoff"
    log_info "Meta: menor conjunto de ASVs (todos os tratamentos <= ${max_asv}) que cubra X% das reads"
    echo ""

    best_cutoff=$(python3 - "${txt_files[@]}" << PYEOF
import sys, bisect

files = sys.argv[1:]
max_asv = ${max_asv}

# Pré-calcular curva de abundância acumulada por tratamento
# ASVs ordenados do mais abundante ao menos; acum[n] = % acumulada dos top n+1 ASVs
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
    soma_total = sum(vals)
    if soma_total == 0:
        sys.stderr.write(f"  AVISO: {fpath.split('/')[-1]} tem soma total = 0\n")
        curvas[fpath] = []
        continue
    pcts = sorted((v / soma_total * 100 for v in vals), reverse=True)
    acum = []
    s = 0.0
    for p in pcts:
        s += p
        acum.append(round(s, 8))
    curvas[fpath] = acum
    sys.stderr.write(
        f"  {fpath.split('/')[-1]}: {len(acum)} ASVs  "
        f"soma_total={int(soma_total):,}  "
        f"abr_min={pcts[-1]:.6f}%  abr_max={pcts[0]:.4f}%\n"
    )

# Varredura do maior cutoff para o menor (0.01% a 100%, passo 0.01%)
# Busca o cutoff mais restritivo (maior %) onde N_mantidos <= max_asv em todos
# N_mantidos = bisect_left(acum, cutoff) + 1  (quantos ASVs cobrem >= cutoff%)
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
    # Imprimir resultado por tratamento no stderr (vai para o log)
    sys.stderr.write(f"\n  Cutoff selecionado: {best}% acumulado\n")
    for fpath, acum in curvas.items():
        idx = bisect.bisect_left(acum, best - 1e-9)
        n = idx + 1
        sys.stderr.write(
            f"  {fpath.split('/')[-1]}: {n} ASVs mantidos  "
            f"({acum[idx]:.4f}% acumulado real)\n"
        )
    print(best)
PYEOF
)

    if [ "$best_cutoff" = "-1" ] || [ -z "$best_cutoff" ]; then
        log_error "Nenhum cutoff entre 0.01% e 100% deixou todos os tratamentos com <= ${max_asv} ASVs."
        log_warn "Considere aumentar max_asv no script."
        exit 1
    fi

    log_ok "Cutoff --abr selecionado: ${best_cutoff}% de abundância acumulada (mesmo corte para todos)"
    filter_mode="rel_acum"
    cutoff_value="$best_cutoff"

elif [ -n "$min_ab_pct" ]; then
    log_step "Modo -ab manual: abundância relativa mínima = ${min_ab_pct}% da soma total do tratamento"
    filter_mode="rel"
    cutoff_value="$min_ab_pct"

elif [ -n "$min_seqs" ]; then
    log_step "Modo -n manual: soma mínima de reads = ${min_seqs}"
    filter_mode="abs"
    cutoff_value="$min_seqs"

else
    log_warn "Nenhum modo de filtragem especificado (--abs, --abr, -n ou -ab). Prosseguindo sem filtragem."
    filter_mode="none"
fi

# =============================================================================
# APLICAR FILTRAGEM + RELATÓRIO DE % MANTIDOS
# =============================================================================
filtered_folder="${folder}/filtered"
mkdir -p "$filtered_folder"

log_step "Aplicando filtros e gravando em: ${filtered_folder}/"
echo ""

# Cabeçalho do relatório
echo "------------------------------------------------------------"
printf "  %-30s %6s %6s %8s\n" "Arquivo" "Antes" "Depois" "Mantido"
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
log_ok "Filtragem concluída. Arquivos em: ${filtered_folder}/"

# =============================================================================
# SPARCC
# =============================================================================
log_step "Iniciando análises SparCC"
echo ""

cd "$filtered_folder" || exit 1

shopt -s nullglob
filtered_files=(*.txt)
shopt -u nullglob

if [ ${#filtered_files[@]} -eq 0 ]; then
    log_error "Nenhum arquivo .txt na pasta filtered/."
    exit 1
fi

total=${#filtered_files[@]}
count=0
failed_files=()

for file in "${filtered_files[@]}"; do
    count=$((count+1))
    base_name="${file%.txt}"

    echo ""
    echo "--- [$count/$total] Processando: $file ---"

    net_dir="${base_name}_net"
    mkdir -p "${net_dir}/perm"
    mkdir -p "${net_dir}/pvalues"

    # --- SparCC principal ---
    log_info "Calculando correlações SparCC..."
    python "${SPARCC_PATH}/SparCC.py" "${file}" \
        --cor_file="${net_dir}/cor_sparcc.out"
    if [ $? -ne 0 ]; then
        log_error "SparCC.py falhou em: $file — pulando."
        failed_files+=("$file")
        continue
    fi

    # Verifica se o arquivo de saída tem dados reais
    if ! check_cor_file "${net_dir}/cor_sparcc.out"; then
        log_error "cor_sparcc.out está vazio (sem valores numéricos) para: $file"
        log_warn "Possível causa: formato incorreto da tabela de entrada (separador, header, etc.)"
        log_warn "Verifique: head -2 ${file} | cat -A"
        failed_files+=("$file")
        continue
    fi
    log_ok "Correlações calculadas: ${net_dir}/cor_sparcc.out"

    # --- Bootstraps ---
    log_info "Gerando 100 permutações..."
    python "${SPARCC_PATH}/MakeBootstraps.py" "${file}" \
        -n 100 \
        -t permutation_#.txt \
        -p "${net_dir}/perm/"
    if [ $? -ne 0 ]; then
        log_error "MakeBootstraps.py falhou em: $file — pulando."
        failed_files+=("$file")
        continue
    fi
    log_ok "Permutações geradas em: ${net_dir}/perm/"

    # --- SparCC nas permutações ---
    log_info "Calculando SparCC para cada permutação (0–99)..."
    perm_failed=0
    for f in $(seq 0 99); do
        python "${SPARCC_PATH}/SparCC.py" \
            "${net_dir}/perm/permutation_${f}.txt" \
            -i 100 \
            --cor_file="${net_dir}/pvalues/perm_cor${f}.txt"
        if [ $? -ne 0 ]; then
            log_warn "SparCC falhou na permutação ${f} — continuando."
            perm_failed=$((perm_failed+1))
        fi
    done
    if [ $perm_failed -gt 0 ]; then
        log_warn "${perm_failed}/100 permutações falharam para: $file"
    fi
    log_ok "SparCC permutacional concluído."

    # --- PseudoPvals ---
    log_info "Calculando pseudo p-valores (two-sided)..."
    python "${SPARCC_PATH}/PseudoPvals.py" \
        "${net_dir}/cor_sparcc.out" \
        "${net_dir}/pvalues/perm_cor#.txt" \
        100 \
        -o "${net_dir}/pvals_two_sided.txt" \
        -t two_sided
    if [ $? -ne 0 ]; then
        log_error "PseudoPvals.py falhou em: $file — pulando."
        failed_files+=("$file")
        continue
    fi
    log_ok "P-valores: ${net_dir}/pvals_two_sided.txt"

    # --- Pares significativos ---
    log_info "Extraindo pares significativos..."
    cd "${net_dir}" || exit 1
    python "${SPARCC_PATH}/get_significant_pairs.py"
    if [ $? -ne 0 ]; then
        log_error "get_significant_pairs.py falhou em: ${net_dir}/"
        log_warn "Verifique o formato de cor_sparcc.out e pvals_two_sided.txt"
        cd ..
        failed_files+=("$file")
        continue
    fi
    log_ok "Pares significativos extraídos em: ${net_dir}/"
    cd ..

done

# =============================================================================
# Sumário final
# =============================================================================
echo ""
echo "============================================"
if [ ${#failed_files[@]} -eq 0 ]; then
    log_ok "Pipeline concluído com sucesso! ($total/$total arquivos)"
else
    log_warn "Pipeline concluído com erros."
    log_warn "Arquivos que falharam (${#failed_files[@]}/${total}):"
    for f in "${failed_files[@]}"; do
        echo "    - $f"
    done
fi

# Perfil dos tratamentos após filtragem
echo ""
echo "--------------------------------------------"
echo "  Perfil dos tratamentos (após filtragem):"
echo "--------------------------------------------"
for file in "${filtered_folder}"/*.txt; do
    print_treatment_profile "$file" "$(basename "$file")"
done
echo "--------------------------------------------"
echo ""
echo "  Log completo salvo em:"
echo "  $log_file"
echo "============================================"
echo ""
