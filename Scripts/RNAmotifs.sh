# Configuration
# --------------------

echo ""
echo "██████╗  ███╗   ██╗  █████╗  ███╗   ███╗  ██████╗  ████████╗ ██╗ ███████╗ ███████╗"
echo "██╔══██╗ ████╗  ██║ ██╔══██╗ ████╗ ████║ ██╔═══██╗ ╚══██╔══╝ ██║ ██╔════╝ ██╔════╝"
echo "██████╔╝ ██╔██╗ ██║ ███████║ ██╔████╔██║ ██║   ██║    ██║    ██║ █████╗   ███████╗"
echo "██╔══██╗ ██║╚██╗██║ ██╔══██║ ██║╚██╔╝██║ ██║   ██║    ██║    ██║ ██╔══╝   ╚════██║"
echo "██║  ██║ ██║ ╚████║ ██║  ██║ ██║ ╚═╝ ██║ ╚██████╔╝    ██║    ██║ ██║      ███████║"
echo "╚═╝  ╚═╝ ╚═╝  ╚═══╝ ╚═╝  ╚═╝ ╚═╝     ╚═╝  ╚═════╝     ╚═╝    ╚═╝ ╚═╝      ╚══════╝"
echo ""

# Before running this script, please adjust folders (paths) in the r.config and n.config files.

RNAmotifs=$1 # SET THE ABSOLUTE PATH TO RNAMOTIFS FOLDER

cd $RNAmotifs

export PYTHONPATH=$RNAmotifs

result_dir=$2

splicing_change_infile=$3

# Search Configuration options
#
# ================ ============== ======
# variable         type           value
# ================ ============== ======
# folder_root      str            root folder of m3, e.g. /home/user/m3
# regions_name     str            name of the regions folder and file (regions/<regions_name>/<regions_name>.tab) from which regions will be read
# genome           str            which genome to use (e.g. mm9, hg19,...)
# cluster_hw       int            half window size used for clustering, default 15
# h_min            int            minimal h at the thresholding step, default 4
# pth              []             list of pth values to compute the thresholding for, default [0.1, 0.25, 0.5, 0.75, 1]
# bootstrap        int            number of bootstrapping steps, default 1000
# ================ ============== ======

species=$4

if [ $species='hg19' ]; then
   organism="Human"
else
   organism="Mouse"
fi

half_window=$5

min_height=4

percent_coverage=$6

enrichment_wd=$7

bootstraps=$8

ncores=$9

p_fisher=${10}

p_empirical=${11}

inIntron=${12}

inExon=${13}

inExon_MAP=${14}

inIntron_MAP=${15}

echo "========================= ======================="
echo "Parameters                Values"
echo "========================= ======================="
echo "Splicing change infile   = ${splicing_change_infile}"
echo "Genome ref               = ${species}"
echo "Organism                 = ${organism}"
echo "Cluster half window      = ${half_window}"
echo "Cluster min height       = ${min_height}"
echo "Percent cluster coverage = ${percent_coverage}"
echo "Enrichment window        = ${enrichment_wd}"
echo "Bootstraps               = ${bootstraps}"
echo "N cores                  = ${ncores}"
echo "P.Value (Fisher)         = ${p_fisher}"
echo "P.Value (Empirical)      = ${p_empirical}"
echo "N cores                  = ${ncores}"
echo "inIntron                 = ${inIntron}"
echo "inExon                   = ${inExon}"
echo "inIntron                 = ${inIntron_MAP}"
echo "inExon                   = ${inExon_MAP}"
echo "========================= ======================="

# ---------------------------------------------------
# Run
# ---------------------------------------------------
echo "[*] Config ..."

N_config="${RNAmotifs}/results/RNAmotifs/config/${result_dir}.n.config"
R_config="${RNAmotifs}/results/RNAmotifs/config/${result_dir}.r.config"

if [ ! -e results/RNAmotifs ] ; then
  mkdir results/RNAmotifs
fi

if [ ! -e results/RNAmotifs/$result_dir ] ; then
  mkdir results/RNAmotifs/$result_dir
fi

if [ ! -e results/RNAmotifs/config ] ; then
  mkdir results/RNAmotifs/config
fi

# writing gMotifs configuration files
if [ -f "$N_config" ] ; then
  rm $N_config $R_config
fi

if [ ! -f "$N_config" ] ; then
  # if not create the file
  touch "$N_config"
fi

echo "organism=$organism" >> $N_config
echo "splicing_change_file=$splicing_change_infile" >> $N_config
echo "tetramer_folder=$RNAmotifs/tetramers/$result_dir" >> $N_config
echo "results_folder=$RNAmotifs/results/RNAmotifs/$result_dir" >> $N_config

cp $N_config $R_config

echo "search=N" >> $N_config
echo "search=R" >> $R_config

echo "enrichment_window=$enrichment_wd" >> $N_config
echo "enrichment_window=$enrichment_wd" >> $R_config

echo "inIntron=$inIntron" >> $N_config
echo "inIntron=$inIntron" >> $R_config

echo "inExon=$inExon" >> $N_config
echo "inExon=$inExon" >> $R_config

# ---------------------------------------------------
# prepare non-overlapping regions from splicing file
# ---------------------------------------------------

echo "[*] Preparing regions from splicing file ..."
if [ -e m3_light/regions/$result_dir ] ; then
	echo "Dir already present"
else
	mkdir m3_light/regions/$result_dir
fi

python m3_light/regions/prepare_regions.py $splicing_change_infile m3_light/regions/$result_dir/overlapping.tab
python m3_light/regions/merge_overlapping_regions.py m3_light/regions/$result_dir/overlapping.tab m3_light/regions/$result_dir/$result_dir.tab
rm m3_light/regions/$result_dir/overlapping.tab


# ---------------------------------------------------
# find tetramers in the genome
# ---------------------------------------------------
echo "[*] Searching tetramers ..."
TMP_FILE="tmp_${result_dir}.py"
echo "import m3_light" > $TMP_FILE
echo "m3_light.find_tetramers('${result_dir}', '${species}', '${half_window}','${min_height}','${percent_coverage}')" >> $TMP_FILE
python $TMP_FILE
rm $TMP_FILE

# move results to tetramers folder
if [ ! -e tetramers ] ; then
  mkdir tetramers
fi

if [ ! -e tetramers/$result_dir ] ; then
  mkdir tetramers/$result_dir
  mkdir tetramers/$result_dir/r
  mkdir tetramers/$result_dir/nr
fi

m3_light_pwd="m3_light/results/${result_dir}/motifs.final/pth_${percent_coverage}"

cp $m3_light_pwd/A* tetramers/$result_dir/nr
cp $m3_light_pwd/C* tetramers/$result_dir/nr
cp $m3_light_pwd/T* tetramers/$result_dir/nr
cp $m3_light_pwd/G* tetramers/$result_dir/nr

cp $m3_light_pwd/W* tetramers/$result_dir/r
cp $m3_light_pwd/Y* tetramers/$result_dir/r
cp $m3_light_pwd/S* tetramers/$result_dir/r
cp $m3_light_pwd/R* tetramers/$result_dir/r


# # ---------------------------------------------------
# # counting tetramer occurences in the genome
# # ---------------------------------------------------

echo "[*] Counting tetramer occurrences in the genome ..."

# tetramer: calculate the positional specific tetramer occurrences along the RNA splicing map
./gMotifs/build/tetramer -c $N_config
./gMotifs/build/tetramer -c $R_config

# couting:  count the number of exons with tetramer occurrencies in the region of interest R1,R2,R3
./gMotifs/build/counting -c $N_config
./gMotifs/build/counting -c $R_config

# ---------------------------------------------------
# bootstrap and tetramer selection
# ---------------------------------------------------

rm -fr tetramers/$result_dir
rm -fr m3_light/results/${result_dir}
rm -fr m3_light/regions/${result_dir}

cd Scripts

echo "[*] Boostrapping ..."
Rscript bootstrap-FDR.R ../results/RNAmotifs $result_dir $bootstraps $ncores

echo "[*] RNAmaps ..."
Rscript selection_of_tetramers.R ../results/RNAmotifs $result_dir $bootstraps $p_fisher $p_empirical $inExon_MAP $inIntron_MAP
Rscript sign_reg_plot.R ../results/RNAmotifs $result_dir
Rscript selection_of_tetramers_old_subset.R ../results/RNAmotifs $result_dir $bootstraps $p_fisher $p_empirical $inExon_MAP $inIntron_MAP
