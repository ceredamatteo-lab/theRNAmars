# RNAMaRs
RNAMaRs is a computational framework that integrates multivalent RNA motifs (MRMs) discovery with in vivo RBP binding and functional splicing evidence to infer MRM-RBP associations.

<img width="2730" height="755" alt="RNAMaRs_overview" src="https://github.com/user-attachments/assets/464ac58d-81a4-4a5a-817f-9620b382e0f3" />


## Installation

To start using RNAMaRs, clone this repository on your machine:

```
git clone https://github.com/ceredamatteo-lab/theRNAmars.git
```
It is also recommended to install all required packages using the provided Conda `.yml` file:

```
conda env create -f RNAmars_env.yml
conda activate RNAmars_env
```

RNAMaRs builds on top of RNAmotifs tool ([Cereda et al., 2014](https://link.springer.com/article/10.1186/gb-2014-15-1-r20)) results. By cloning the present repository, you already have all necessary files to install and run RNAmotifs without the need to clone the original folder.

First, add RNAMaRs complete path to your Python path environment variable:

```
export PYTHONPATH=$PYTHONPATH:/path_to_RNAmars_root_folder
```

Then, from RNAMaRs root, move to `m3_light/genomes` folder and run `hg19_download.sh`: this script will automatically install hg19 genome under `m3_light/genomes/hg19`.

```
cd m3_light/genomes
bash hg19.download.sh
```


## Input

RNAMaRs requires a semicolon-delimited input file (in `.txt` format) where each row represents a cassette exon.
Specifically, columns must contain the following information:

+ `rowID`
+ `eventID`
+ `chr`
+ `strand`
+ `upstream_exon_end`
+ `exonStart_0base`
+ `exonEnd`
+ `downstream_exon_start`
+ `dIrank`: a user-defined column specifying if the exon is enhanced (1), silenced (-1) or constitutive (0).

This is an example of how the input file should look like:

```
1;33;chrX;-;153583440;153585618;153585642;153585801;1
2;104;chrX;+;149901202;149905066;149905252;149905713;1
3;107;chrX;+;149919315;149919501;149919712;149924160;-1
4;116;chrX;-;148564749;148568455;148568629;148571844;1
5;180;chrX;+;123022568;123025087;123025159;123026580;0
```



## Usage
RNAMaRs runs from command-line with `RNAMaRs.sh`. At first, the script runs RNAmotifs internally; then, the RNAMaRs module computes and plots the association score (AS) between RNAmotifs-based enriched tetramers and RBP binding evidence previously computed.

```
./RNAMaRs.sh $NAME $DIR_RNAMARS $INPUT_EXONS_FILE $CELL_LINE $DESEQ_FILE $DIR_ECLIP $CPUS
```

`RNAMaRs.sh` requires the following command-line parameters:

| Syntax | Description | Example
| ----------- | ----------- | ----------- |
NAME | Name of the folder used to save RNAmotifs and RNAMaRs results | "HNRNPK_silencing_72h"
DIR_RNAMARS | Path to RNAMaRs tool | "path/to/local/theRNAmars/"
INPUT_EXONS_FILE | Text file containing alternative and constitutive exons (organized as described in `Input` section) | "HNRNPK_silencing_72h_input_rnamotifs.txt"
CELL_LINE | Cell line to use as reference (it must be: `HepG2` or `K562`) | "HepG2"
DESEQ_FILE | Path to differential gene expression table (optional). File can be both in `.tsv` or `.rds` format. If not available, just leave `""` | "DESeq2_HNRNPK_silencing_72h_DEG.csv"
DIR_ECLIP | Path to processed iCount peaks. This must contain `HepG2` and `K562` subdirectories | "DIR_RNAMARS/Sources/eCLIP/"
CPUS | How many CPUs to use | 14

Additionally, user could run RNAmotifs and RNAMaRs separately using `01_RNAmotifs.sh` and `02_RNAmars.sh` from command-line as follows:

```
./01_RNAmotifs.sh $NAME $DIR_RNAMOTIFS $INPUT_EXONS_FILE $CELL_LINE $CPUS
```

```
./02_RNAMaRs.sh $NAME $DIR_RNAMARS $INPUT_EXONS_FILE $CELL_LINE $DESEQ_FILE $DIR_ECLIP $CPUS
```

Bear in mind that RNAMaRs module will not work if RNAmotifs has not been previously run.

`DESEQ_FILE` must be in the following format:

```
baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	gene_name	deg
ENSG00000000003.14	1207.23729680793	-0.594757925435271	0.13697550117035	-4.34207519120954	1.41143208787687e-05	3.9999563754732e-05	TSPAN6	TRUE
ENSG00000000419.12	1366.72558786783	0.0858651981950839	0.123396414060123	0.695848407339028	0.486523770255723	0.575036991227159	DPM1	FALSE
```
with gene IDs as rownames. Stratification of genes as DEGs (TRUE/FALSE) can be added in the `deg` column for visualization purposes. If `deg` column is not provided, RNAMaRs will automatically consider as DEG genes with `padj` ≤ 0.1 and `|log2FoldChange|` ≥ 0.1.

## Results
Results from the two modules are stored in `DIR_RNAMARS/results/RNAmotifs/` and `DIR_RNAMARS/results/RNAmars/`, respectively.

`DIR_RNAMARS/results/RNAmotifs/` is further divided in additional result folders, one per each optimal parameter combination.

`DIR_RNAMARS/results/RNAmars/NAME_on_CELL_LINE` contains the results of the RNAMaRs run, consisting in:

| Format | Description |
|--------|------------|
| `summary_results_NAME_hw_x_ew_y.rds` | Summary files of RNAmotifs results, one per each optimal parameter combination x,y |
| `SCORE1_{enh_sil}_signal_recovery_rate_{NAME}.rds` | Signal recovery rates for the specific RBP |
| `SCORE2_{enh_sil}_profile_similarity_{NAME}.rds` | Cosine similarity between RNAmotifs-derived tetramer enrichment scores and reference cell line binding profile |
| `final_mat_prot_score_{enh/sil}.rds` | Protein score matrix |
| `final_mat_tet_score_{enh/sil}.rds` | Tetramer score matrix |
| `final_mat_{enh/sil}.rds` | Association score matrix |
| `{CELL_LINE}_{NAME}_{enh/sil}.pdf` | Association score heatmaps |

## Output
The graphical output of RNAMaRs are two heatmaps - one for enhanced exons and one for silenced exons - that quantifies the association between MRMs (columns) and cognate RBPs (rows). In this example, these heatmaps were obtained using exons from HNRNPK depletion in PC3 cells, as described in the paper: you can find the file at `input/HNRNPK_silencing_72h_input_rnamotifs.txt`.

<p align="center"><img width="770" height="228" alt="RNAMaRs_overview" src="https://github.com/user-attachments/assets/d5c54331-7eb1-4cc0-b5b6-5a90f786a9f1" /></p>

From left to right, the plot reports:
+ Results of differential gene expression analysis by means of `|log2FoldChange|`. Direction of regulation is also shown.
+ Binding profile of the RBP in the reference cell line, restricted to regulation type (enhanced/silenced) in the three regions under analysis: upstream intron (R1), exon (R2) and downstream intron (R3).
+ A bar plot of the mean AS of the RBP with all enriched tetramers. This value was used to sort RBPs along the rows.
+ Main heatmap: dot plot reporting the RBP-tetramer AS. This value was obtained by the product of SCORE1 (signal recovery rate) and SCORE2 (cosine similarity). Dot color and size of dots are proportional to AS.
+ Top heatmap reports combined regional enrichments (CRE) of MRMs across RNAmotifs runs, aggregated using Fisher's method. Top bar plot reports summed CRE values across regions, namely cumulative CRE (CCRE).
+ RNAmotifs-based splicing maps.


## Contributors
RNAMaRs has been designed by Dr Matteo Cereda and Uberto Pozzoli and developed with Tommaso Becchi, Gabriele Boscagli, Mariachiara Grieco and Francesca Priante.

Main developer: Tommaso Becchi, Gabriele Boscagli and Mariachiara Grieco.

Contributing developers: Francesca Priante.

Contributions are always welcome!

## License
Please read the [License](https://github.com/ceredamatteo-lab/theRNAmars?tab=MIT-1-ov-file) first.

