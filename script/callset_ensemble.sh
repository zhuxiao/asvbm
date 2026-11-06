#!/bin/bash
set -euo pipefail
usage() {
    cat <<'EOF'
Usage:
  callset_ensemble.sh -B BAM -f REFERENCE -T TOOLS -x TYPE \
                      -i INPUT_VCF [INPUT_VCF ...] \
                      [-C CHROMOSOMES] [-o OUTPUT]

Required options:
  -B FILE        Input BAM file
  -f FILE        Reference FASTA file
  -T STRING      Tools separated by semicolons, e.g., "tool1;tool2;tool3"
  -x STRING      Sequence type, e.g., "ccs" or "ont"
  -i FILE ...    One or more input VCF files

Optional options:
  -C STRING      Chromosomes, e.g., "chr1-22;chrX;chrY"
  -o DIR         Output directory
  -h             Show this help
EOF
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

require_value() {
    [[ $# -ge 2 && -n "$2" ]] || die "option $1 requires a value"
}

_arg_bam=""
_arg_reference=""
_arg_chromosomes=""
_arg_tools=""
_arg_output=""
_arg_seq_type=""
_arg_input_vcfs=()

while (($# > 0)); do
    case "$1" in
        -B)
            require_value "$@"
            _arg_bam="$2"
            shift 2
            ;;
        -f)
            require_value "$@"
            _arg_reference="$2"
            shift 2
            ;;
        -C)
            require_value "$@"
            _arg_chromosomes="$2"
            shift 2
            ;;
        -T)
            require_value "$@"
            _arg_tools="$2"
            shift 2
            ;;
        -o)
            require_value "$@"
            _arg_output="$2"
            shift 2
            ;;
        -x)
            require_value "$@"
            _arg_seq_type="$2"
            shift 2
            ;;
        -i)
            shift
            _input_count=0

            # Collect all VCF paths after -i until the next option.
            while (($# > 0)); do
                case "$1" in
                    -B|-f|-C|-T|-o|-x|-i|-h|--help)
                        break
                        ;;
                    --)
                        shift
                        if (($# > 0)); then
                            _arg_input_vcfs+=("$@")
                            _input_count=$((_input_count + $#))
                        fi
                        set --
                        break
                        ;;
                    -*)
                        die "unknown option after -i: $1"
                        ;;
                    *)
                        _arg_input_vcfs+=("$1")
                        _input_count=$((_input_count + 1))
                        shift
                        ;;
                esac
            done

            ((_input_count > 0)) ||
                die "option -i requires at least one input VCF"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            (($# == 0)) ||
                die "unexpected positional arguments: $*"
            ;;
        -*)
            die "unknown option: $1"
            ;;
        *)
            die "unexpected positional argument: $1"
            ;;
    esac
done

[[ -n "$_arg_bam" ]]       || die "missing required option -B"
[[ -n "$_arg_reference" ]] || die "missing required option -f"
[[ -n "$_arg_tools" ]]     || die "missing required option -T"
[[ -n "$_arg_seq_type" ]]  || die "missing required option -x"
(( ${#_arg_input_vcfs[@]} > 0 )) ||
    die "missing required option -i"

# Restore the parameter order used by the original pipeline.
# Empty -C and -o values are retained so that their defaults can be applied.
set -- \
    "$_arg_bam" \
    "$_arg_reference" \
    "$_arg_chromosomes" \
    "$_arg_tools" \
    "$_arg_output" \
    "$_arg_seq_type" \
    "${_arg_input_vcfs[@]}"

unset \
    _arg_bam \
    _arg_reference \
    _arg_chromosomes \
    _arg_tools \
    _arg_output \
    _arg_seq_type \
    _arg_input_vcfs \
    _input_count

# ================= 1. Parameter Validation & Extraction =================
if [ $# -lt 7 ]; then
    echo "❌ Error: Insufficient arguments!"
    usage
    exit 1
fi

BAM_FILE="$1"
REF_FILE="$2"
CHROMOSOMES="$3"
TOOLS="$4"
OUTPUT_PREFIX="$5"
SEQ_TYPE="$6"
shift 6
ORIGINAL_INPUT_FILES=("$@")

# Check sequence type
if [[ "$SEQ_TYPE" != "ccs" && "$SEQ_TYPE" != "ont" ]]; then
    echo "❌ Error: Sequence type must be 'ccs' or 'ont', but got: $SEQ_TYPE"
    exit 1
fi

CURRENT_DIR=$(pwd)
IFS=';' read -r -a TOOL_ARRAY <<< "$TOOLS"

# If -C is omitted, do not pass -C to asvbm.
# asvbm will then use all chromosomes by default.
ASVBM_C_ARGS=()
if [[ -n "$CHROMOSOMES" ]]; then
    ASVBM_C_ARGS=(-C "$CHROMOSOMES")
fi

# ================= Create Unified Output Directories =================
# With -o PREFIX:
#   PREFIX_output/PREFIX
#
# Without -o:
#   OUTPUT_DIR is "output".
#   All generated files and directories are kept under output/.
# Later internal creates use "_fp_merge" and "_final" suffixes
#   so that their outputs do not overwrite each other.
if [[ -z "$OUTPUT_PREFIX" ]]; then
    OUTPUT_PREFIX="output"
    OUTPUT_DIR="output"
else
    OUTPUT_DIR="${OUTPUT_PREFIX}"
fi

ASVBM_OUTPUT_BASE="${OUTPUT_DIR%/}/${OUTPUT_PREFIX}_base"
FP_MERGE_OUTPUT="${OUTPUT_DIR%/}/${OUTPUT_PREFIX}_fp_merge"
FINAL_OUTPUT="${OUTPUT_DIR%/}/${OUTPUT_PREFIX}_final"
STAT_OUTPUT="${OUTPUT_DIR%/}/${OUTPUT_PREFIX}_stat"
CUTEFC_DIR="${OUTPUT_DIR%/}/output_cuteFC"

CUTEFC_DIR="$(realpath -m -- "$CUTEFC_DIR")"

STEP2_OUTPUT_ARGS=(-o "$ASVBM_OUTPUT_BASE")

mkdir -p "$OUTPUT_DIR" "$CUTEFC_DIR"
echo "📁 Output base: $(pwd)/${ASVBM_OUTPUT_BASE#./}"

# ================= Dynamically Set Parameters Based on Data Type =================
if [[ "$SEQ_TYPE" == "ccs" ]]; then
    INS_BIAS=1000; INS_RATIO=0.9; DEL_BIAS=1000; DEL_RATIO=0.5
else
    INS_BIAS=1000; INS_RATIO=0.5; DEL_BIAS=1000; DEL_RATIO=0.5
fi

# If data type is ont, append -I 0.7 parameter for asvbm create
ASVBM_I_PARAM=""
if [[ "$SEQ_TYPE" == "ont" ]]; then
    ASVBM_I_PARAM="-I 0.7"
fi

echo "🚀 [Pipeline Start] Sample: $OUTPUT_PREFIX | Data Type: $SEQ_TYPE | Number of Input Files: ${#ORIGINAL_INPUT_FILES[@]}"

# ================= 2. Step 1 & Step 2 Pre-processing =================
STEP2_INPUT_FILES=()
for i in "${!ORIGINAL_INPUT_FILES[@]}"; do
    INPUT_FILE="${ORIGINAL_INPUT_FILES[$i]}"
    INDEX=$((i + 1))

    echo "🧬 [Step 1.$INDEX] Running cuteFC on file: $(basename "$INPUT_FILE")"
    CUTEFC_OUTPUT="${CUTEFC_DIR}/${OUTPUT_PREFIX}_${SEQ_TYPE}_cuteFC_${INDEX}.vcf"

CUTEFC_INPUT="$INPUT_FILE"
CUTEFC_BAM="$BAM_FILE"
CUTEFC_REF="$REF_FILE"

CUTEFC_INPUT="$(realpath -m -- "$CUTEFC_INPUT")"
CUTEFC_BAM="$(realpath -m -- "$CUTEFC_BAM")"
CUTEFC_REF="$(realpath -m -- "$CUTEFC_REF")"

cuteFC -t 32 -l 20 -md 500 -mi 500 \
    --max_cluster_bias_INS "$INS_BIAS" \
    --diff_ratio_merging_INS "$INS_RATIO" \
    --max_cluster_bias_DEL "$DEL_BIAS" \
    --diff_ratio_merging_DEL "$DEL_RATIO" \
    --genotype -Ivcf "$CUTEFC_INPUT" \
    "$CUTEFC_BAM" "$CUTEFC_REF" \
    "$CUTEFC_OUTPUT" "$CUTEFC_DIR"

    FILTERED_VCF="${CUTEFC_DIR}/${OUTPUT_PREFIX}_${SEQ_TYPE}_cuteFC_${INDEX}_filtered.vcf"
    echo "📊 [Step 2.$INDEX-1] Filtering AF >= 0.25 or RE >= 10..."
    bcftools filter -i '(INFO/AF >= 0.25) || (INFO/RE >= 10)' "$CUTEFC_OUTPUT" \
        -o "$FILTERED_VCF" -O v

    FINAL_VCF="${CUTEFC_DIR}/${OUTPUT_PREFIX}_${SEQ_TYPE}_cuteFC_${INDEX}_filtered2.vcf"
    echo "🔬 [Step 2.$INDEX-2] Filtering IMPRECISE variants..."
    awk '
        /^##/ || /^#CHROM/ { print; next }
        {
            has_precise = 1
            n = split($8, info, /;/)
            for (j = 1; j <= n; j++) {
                if (info[j] == "IMPRECISE") {
                    has_precise = 0
                    break
                }
            }
            if (has_precise) print
        }
    ' "$FILTERED_VCF" > "$FINAL_VCF"

    STEP2_INPUT_FILES+=("$FINAL_VCF")
done
echo "✅ Step 1 & 2 pre-processing completed!"

# ================= 3. Step 2: Run asvbm create =================
DUMMY_VCF="${OUTPUT_DIR}/output.vcf"
echo "##fileformat=VCFv4.2" > "$DUMMY_VCF"
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample" >> "$DUMMY_VCF"

echo "🔄 [Step 2] Running asvbm create..."
asvbm create $ASVBM_I_PARAM \
    -B "$BAM_FILE" \
    "${ASVBM_C_ARGS[@]}" \
    -T "$TOOLS" \
    "${STEP2_INPUT_FILES[@]}" \
    "${STEP2_OUTPUT_ARGS[@]}" \
    "$DUMMY_VCF" \
    "$REF_FILE"
echo "✅ asvbm create completed!"

# ================= 4. Step 3: Pre-filtering + asvbm stat =================
BENCHMARK_VCF="${ASVBM_OUTPUT_BASE}/SharedFN_RefinedBenchmark/benchmark.vcf"
BENCHMARK_TOOLSNUM_2PLUS="${ASVBM_OUTPUT_BASE}/SharedFN_RefinedBenchmark/benchmark_toolsNum_2+.vcf"

if [ ! -f "$BENCHMARK_VCF" ]; then
    echo "❌ Error: Benchmark file not found: $BENCHMARK_VCF"
    exit 1
fi

echo "📊 [Step 3.0] Filtering benchmark_toolsNum_2+..."
bcftools filter -i 'INFO/TOOLSNUM >= 2' "$BENCHMARK_VCF" \
    -o "$BENCHMARK_TOOLSNUM_2PLUS" \
    -O v

echo "📊 [Step 3] Running asvbm stat..."
asvbm stat \
    -B "$BAM_FILE" \
    "${ASVBM_C_ARGS[@]}" \
    -T "$TOOLS" \
    "${ORIGINAL_INPUT_FILES[@]}" \
    -o "$STAT_OUTPUT" \
    "$BENCHMARK_TOOLSNUM_2PLUS" \
    "$REF_FILE"
echo "✅ asvbm stat completed!"

# ================= 5. Step 4: Rename files and run asvbm create =================

STEP4_INPUT_FILES=()
for i in "${!TOOL_ARRAY[@]}"; do
    TOOL_NAME="${TOOL_ARRAY[$i]}"
    INDEX=$((i + 1))

    OLD_FILE="${STAT_OUTPUT}/${TOOL_NAME}/2_num_stat/num_stat_FP.vcf"
    NEW_FILE="${STAT_OUTPUT}/${TOOL_NAME}/2_num_stat/num_stat_FP_${INDEX}.vcf"

    if [ -f "$OLD_FILE" ]; then
        mv "$OLD_FILE" "$NEW_FILE"
        echo "   📝 Renamed: $(basename "$NEW_FILE")"
        STEP4_INPUT_FILES+=("$NEW_FILE")
    else
        echo "⚠️  Warning: File not found $OLD_FILE, skipped."
    fi
done

echo "🔄 [Step 4] Running asvbm create (FP merge)..."
asvbm create $ASVBM_I_PARAM \
    -B "$BAM_FILE" \
    "${ASVBM_C_ARGS[@]}" \
    -T "$TOOLS" \
    "${STEP4_INPUT_FILES[@]}" \
    -o "$FP_MERGE_OUTPUT" \
    "$DUMMY_VCF" \
    "$REF_FILE"
echo "✅ asvbm create (FP merge) completed!"

# ================= 6. Step 5: Filter TOOLSNUM >= 4 using awk =================
echo "📊 [Step 5] Filtering TOOLSNUM >= 4..."
FP_BENCHMARK="${FP_MERGE_OUTPUT}/SharedFN_RefinedBenchmark/benchmark.vcf"
FP_BENCH_DIR="${FP_MERGE_OUTPUT}/SharedFN_RefinedBenchmark"

if [ ! -f "$FP_BENCHMARK" ]; then
    echo "❌ Error: FP Benchmark file not found: $FP_BENCHMARK"
    exit 1
fi

awk 'BEGIN { FS="\t"; OFS="\t" }
/^#/ { print; next }
{
    n = split($8, info, ";")
    for (i = 1; i <= n; i++) {
        if (info[i] ~ /^TOOLSNUM=/) {
            split(info[i], val, "=")
            if (val[2] + 0 >= 4) {
                print
                break
            }
        }
    }
}' "$FP_BENCHMARK" > "${FP_BENCH_DIR}/Benchmark_toolsNum_4+.vcf"

# ================= 7. Step 6: Filter TOOLSNUM == 2 and == 3 =================
echo "📊 [Step 6] Filtering TOOLSNUM == 2 and == 3..."

awk 'BEGIN { FS="\t"; OFS="\t" }
/^#/ { print; next }
{
    n = split($8, info, ";")
    for (i = 1; i <= n; i++) {
        if (info[i] ~ /^TOOLSNUM=/) {
            split(info[i], val, "=")
            if (val[2] + 0 == 2) {
                print
                break
            }
        }
    }
}' "$FP_BENCHMARK" > "${FP_BENCH_DIR}/Benchmark_toolsNum_2.vcf"

awk 'BEGIN { FS="\t"; OFS="\t" }
/^#/ { print; next }
{
    n = split($8, info, ";")
    for (i = 1; i <= n; i++) {
        if (info[i] ~ /^TOOLSNUM=/) {
            split(info[i], val, "=")
            if (val[2] + 0 == 3) {
                print
                break
            }
        }
    }
}' "$FP_BENCHMARK" > "${FP_BENCH_DIR}/Benchmark_toolsNum_3.vcf"

# ================= 8. Step 7: Run cuteFC on TOOLSNUM == 2 and == 3 =================
echo "🧬 [Step 7] Running cuteFC on TOOLSNUM == 2 and == 3..."

TOOLS2_VCF="${FP_BENCH_DIR}/Benchmark_toolsNum_2.vcf"
TOOLS3_VCF="${FP_BENCH_DIR}/Benchmark_toolsNum_3.vcf"

cuteFC -t 32 -l 20 -md 500 -mi 500 \
    --max_cluster_bias_INS "$INS_BIAS" \
    --diff_ratio_merging_INS "$INS_RATIO" \
    --max_cluster_bias_DEL "$DEL_BIAS" \
    --diff_ratio_merging_DEL "$DEL_RATIO" \
    --genotype -Ivcf "$TOOLS2_VCF" \
    "$BAM_FILE" "$REF_FILE" \
    "${FP_BENCH_DIR}/Benchmark_toolsNum_2_cuteFC.vcf" \
    "$FP_BENCH_DIR"

cuteFC -t 32 -l 20 -md 500 -mi 500 \
    --max_cluster_bias_INS "$INS_BIAS" \
    --diff_ratio_merging_INS "$INS_RATIO" \
    --max_cluster_bias_DEL "$DEL_BIAS" \
    --diff_ratio_merging_DEL "$DEL_RATIO" \
    --genotype -Ivcf "$TOOLS3_VCF" \
    "$BAM_FILE" "$REF_FILE" \
    "${FP_BENCH_DIR}/Benchmark_toolsNum_3_cuteFC.vcf" \
    "$FP_BENCH_DIR"

# ================= 9. Step 8: Filter PRECISE variants using awk =================
echo "🔬 [Step 8] Filtering PRECISE variants ..."

awk '
    /^##/ || /^#CHROM/ { print; next }
    {
        has_precise = 1
        n = split($8, info, /;/)
        for (i = 1; i <= n; i++) {
            if (info[i] == "IMPRECISE") {
                has_precise = 0
                break
            }
        }
        if (has_precise) print
    }
' "${FP_BENCH_DIR}/Benchmark_toolsNum_2_cuteFC.vcf" \
    > "${FP_BENCH_DIR}/Benchmark_toolsNum_2_cuteFC_filteredPRECISE.vcf"

awk '
    /^##/ || /^#CHROM/ { print; next }
    {
        has_precise = 1
        n = split($8, info, /;/)
        for (i = 1; i <= n; i++) {
            if (info[i] == "IMPRECISE") {
                has_precise = 0
                break
            }
        }
        if (has_precise) print
    }
' "${FP_BENCH_DIR}/Benchmark_toolsNum_3_cuteFC.vcf" \
    > "${FP_BENCH_DIR}/Benchmark_toolsNum_3_cuteFC_filteredPRECISE.vcf"

# ================= 10. Step 9: Final asvbm create merge =================
echo "🏁 [Step 9] Running final merge asvbm create..."

STEP9_INPUT_FILES=(
    "$BENCHMARK_TOOLSNUM_2PLUS"
    "${FP_BENCH_DIR}/Benchmark_toolsNum_4+.vcf"
    "${FP_BENCH_DIR}/Benchmark_toolsNum_3_cuteFC_filteredPRECISE.vcf"
    "${FP_BENCH_DIR}/Benchmark_toolsNum_2_cuteFC_filteredPRECISE.vcf"
)

asvbm create -e $ASVBM_I_PARAM \
    -B "$BAM_FILE" \
    "${ASVBM_C_ARGS[@]}" \
    "${STEP9_INPUT_FILES[@]}" \
    -o "$FINAL_OUTPUT" \
    "$DUMMY_VCF" \
    "$REF_FILE"

BENCHMARK_FINAL="${FINAL_OUTPUT}/SharedFN_RefinedBenchmark/benchmark.vcf"
ALLELE_FINAL="${FINAL_OUTPUT}/SharedFN_RefinedBenchmark/allele.vcf"

if [[ ! -f "$BENCHMARK_FINAL" ]]; then
    die "final benchmark file not found: $BENCHMARK_FINAL"
fi

if [[ ! -f "$ALLELE_FINAL" ]]; then
    die "final allele file not found: $ALLELE_FINAL"
fi

ln -sfn "$(realpath -- "$BENCHMARK_FINAL")" \
    "${OUTPUT_DIR}/benchmark.vcf"

ln -sfn "$(realpath -- "$ALLELE_FINAL")" \
    "${OUTPUT_DIR}/allele.vcf"

rm -f "$DUMMY_VCF"
echo ""
echo "🎉 All 9 steps of the pipeline completed successfully!"
echo "Final output benchmark set: $(pwd)/${FINAL_OUTPUT#./}/SharedFN_RefinedBenchmark/benchmark.vcf"
