NAME=$1
DIR_RNAMARS=$2
INPUT_EXONS_FILE=$3
CELL_LINE=$4
DIR_ECLIP=$5
CPUS=$6
DESEQ_FILE=${7:-NA}
RNAMARS_RESULTS=${DIR_RNAMARS}/results/RNAmars/${NAME}_on_${CELL_LINE}/


echo ""
echo "██████╗ ███╗   ██╗ █████╗ ███╗   ███╗ █████╗ ██████╗ ███████╗"
echo "██╔══██╗████╗  ██║██╔══██╗████╗ ████║██╔══██╗██╔══██╗██╔════╝"
echo "██████╔╝██╔██╗ ██║███████║██╔████╔██║███████║██████╔╝███████╗"
echo "██╔══██╗██║╚██╗██║██╔══██║██║╚██╔╝██║██╔══██║██╔══██╗╚════██║"
echo "██║  ██║██║ ╚████║██║  ██║██║ ╚═╝ ██║██║  ██║██║  ██║███████║"
echo "╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝"
echo ""


cd "${DIR_RNAMARS}"
mkdir -p "${RNAMARS_RESULTS}"

Rscript --vanilla ${DIR_RNAMARS}/Scripts/compute_association_scores.R --input ${INPUT_EXONS_FILE} \
 --name ${NAME} \
 --RNAmotifs_res ${DIR_RNAMARS}/results/RNAmotifs/ \
 --repository ${DIR_RNAMARS}/ \
 --cell_line ${CELL_LINE} \
 --cores ${CPUS} \
 --eclip_path "$DIR_ECLIP/$CELL_LINE/" \
 --output ${RNAMARS_RESULTS}

 
Rscript --vanilla ${DIR_RNAMARS}/Scripts/generate_heatmap.R --input ${INPUT_EXONS_FILE} \
  --name ${NAME} \
  --RNAmotifs_res ${DIR_RNAMARS}/results/RNAmotifs/ \
  --repository ${DIR_RNAMARS}/ \
  --cell_line ${CELL_LINE} \
  --deseq_file "${DESEQ_FILE}" \
  --cores ${CPUS} \
  --output ${RNAMARS_RESULTS}

echo
echo "[*] RNAMaRs completed! Results are available at ${RNAMARS_RESULTS}"
