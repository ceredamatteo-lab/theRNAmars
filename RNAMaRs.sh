NAME=$1
DIR_RNAMARS=$2
INPUT_EXONS_FILE=$3
CELL_LINE=$4
DESEQ_FILE=$5
DIR_ECLIP=$6
CPUS=$7
RNAMARS_RESULTS=${DIR_RNAMARS}/results/RNAmars/${NAME}_on_${CELL_LINE}/


echo ""
echo "██████╗  ███╗   ██╗  █████╗  ███╗   ███╗  ██████╗  ████████╗ ██╗ ███████╗ ███████╗"
echo "██╔══██╗ ████╗  ██║ ██╔══██╗ ████╗ ████║ ██╔═══██╗ ╚══██╔══╝ ██║ ██╔════╝ ██╔════╝"
echo "██████╔╝ ██╔██╗ ██║ ███████║ ██╔████╔██║ ██║   ██║    ██║    ██║ █████╗   ███████╗"
echo "██╔══██╗ ██║╚██╗██║ ██╔══██║ ██║╚██╔╝██║ ██║   ██║    ██║    ██║ ██╔══╝   ╚════██║"
echo "██║  ██║ ██║ ╚████║ ██║  ██║ ██║ ╚═╝ ██║ ╚██████╔╝    ██║    ██║ ██║      ███████║"
echo "╚═╝  ╚═╝ ╚═╝  ╚═══╝ ╚═╝  ╚═╝ ╚═╝     ╚═╝  ╚═════╝     ╚═╝    ╚═╝ ╚═╝      ╚══════╝"
echo ""


echo
echo "========================= ======================="
echo "Parameters                Values"
echo "========================= ======================="
echo "Run name                 = ${NAME}"
echo "RNAmotifs directory      = ${DIR_RNAMARS}"
echo "Input exons file         = ${INPUT_EXONS_FILE}"
echo "Cell line                = ${CELL_LINE}"
echo "RNAMaRs results          = ${RNAMARS_RESULTS}"
echo "DESeq file               = ${DESEQ_FILE}"
echo "eCLIP peaks directory    = ${DIR_ECLIP}"
echo "# cores                  = ${CPUS}"
echo "========================= ======================="
echo

# ---------------------------------------------------
# Run
# ---------------------------------------------------

cd "${DIR_RNAMARS}"
mkdir -p "${RNAMARS_RESULTS}"

if [ "${CELL_LINE}" == 'HepG2'  ]; then
proteins=("HNRNPC" "HNRNPK" "HNRNPU" "NCBP2"  "PRPF8"  "PTBP1"  "QKI"    "RBFOX2" "RBM22"  "SF3A3"  "SF3B4"  "SRSF1"  "U2AF1"  "U2AF2"  "UCHL5" )
elif [ "${CELL_LINE}" == 'K562'  ]; then
proteins=("AGGF1"  "EFTUD2" "FXR1"   "HNRNPU" "PRPF8"  "PTBP1"  "PUS1"   "RBM15"  "SF3B4"  "SRSF1"  "TARDBP" "U2AF1"  "U2AF2" )
fi

hw_list=(5 15 25 35)
ew_list=(30 50 100 200 300)

PARAMS_FILE="Tables/RNAmotifs_optimal_parameters.csv"
PROTEINS_PATTERN=$(IFS="|"; echo "${proteins[*]}")
OPT_PARAMS=$(awk -F',' -v pat="$PROTEINS_PATTERN" 'NR>1 && $2 ~ pat {print $1}' "$PARAMS_FILE" | sort -u)



echo
echo "
=========================
WAIT FOR RNAMOTIFS TO END
Number of runs: $(printf "%s\n" "${OPT_PARAMS[@]}" | sort -u | wc -l)
========================="
echo


JOB_CPUS=4
PARALLEL_JOBS=$(( CPUS / JOB_CPUS ))
(( PARALLEL_JOBS < 1 )) && PARALLEL_JOBS=1

for hw in "${hw_list[@]}"; do
  for ew in "${ew_list[@]}"; do
    key="hw_${hw}_ew_${ew}"
    if grep -qw "$key" <<< "$OPT_PARAMS"; then
      echo "${DIR_RNAMARS}/Scripts/RNAmotifs.sh \
        ${DIR_RNAMARS} ${NAME}_hw_${hw}_ew_${ew} ${INPUT_EXONS_FILE} \
        hg19 ${hw} 0.5 ${ew} 10000 ${JOB_CPUS} \
        0.1 0.0005 1000 1000 30 300"
    fi
  done
done | parallel -j "${PARALLEL_JOBS}"



echo
echo "[*] RNAmotifs runs completed."
echo


echo ""
echo "██████╗ ███╗   ██╗ █████╗ ███╗   ███╗ █████╗ ██████╗ ███████╗"
echo "██╔══██╗████╗  ██║██╔══██╗████╗ ████║██╔══██╗██╔══██╗██╔════╝"
echo "██████╔╝██╔██╗ ██║███████║██╔████╔██║███████║██████╔╝███████╗"
echo "██╔══██╗██║╚██╗██║██╔══██║██║╚██╔╝██║██╔══██║██╔══██╗╚════██║"
echo "██║  ██║██║ ╚████║██║  ██║██║ ╚═╝ ██║██║  ██║██║  ██║███████║"
echo "╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝"
echo ""


Rscript --vanilla ${DIR_RNAMARS}/Scripts/compute_association_scores.R --input ${INPUT_EXONS_FILE} \
  --name ${NAME} \
  --RNAmotifs_res ${DIR_RNAMARS}/results/RNAmotifs \
  --repository ${DIR_RNAMARS}/ \
  --cell_line ${CELL_LINE} \
  --cores ${CPUS} \
  --eclip_path "$DIR_ECLIP/$CELL_LINE/" \
  --output ${RNAMARS_RESULTS}


Rscript --vanilla ${DIR_RNAMARS}/Scripts/generate_heatmap.R --input ${INPUT_EXONS_FILE} \
  --name ${NAME} \
  --RNAmotifs_res ${DIR_RNAMARS}/results/RNAmotifs \
  --repository ${DIR_RNAMARS}/ \
  --cell_line ${CELL_LINE} \
  --deseq_file ${DESEQ_FILE} \
  --cores ${CPUS} \
  --output ${RNAMARS_RESULTS}

echo
echo "[*] RNAMaRs completed! Results are available at ${RNAMARS_RESULTS}"
