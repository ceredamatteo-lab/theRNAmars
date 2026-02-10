## ___________________________
## Script name: generate_heatmap.R
## Purpose of script: Uses the Signal Recovery Rate (SRR) and Cosine Similarity (CS) score matrices derived from the previous step (compute_association_scores) to aggregate them and produce the final heatmap.
## ___________________________

library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="RNAmotifs input txt file", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="RNAmotifs and RNAMaRs run name", metavar="character"),
  make_option(c("-d", "--RNAmotifs_res"), type="character", default=NULL,
              help="RNAmotifs results folder", metavar="character"),
  make_option(c("-r", "--repository"), type="character", default=NULL,
              help="RNAmars repository", metavar="character"),
  make_option(c("-c", "--cell_line"), type="character", default=NULL,
              help="Cell line to compare input exons to (HepG2 or K562)", metavar="character"),
  make_option(c("-e", "--deseq_file"), type="character", default=NULL,
              help="DESEQ2 expression file", metavar="character"),
  make_option(c("-p", "--cores"), type="character", default=NULL,
              help="Number of cpus", metavar="numerics"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

READY_RNAMOTIFS_INPUT = opt$input
name                  = opt$name
DIR_RESULTS_RNAMOTIFS = opt$RNAmotifs_res
repo                  = opt$repository
cell_line_target      = opt$cell_line
deseq_file            = opt$deseq_file
cpus                  = as.numeric(opt$cores)
DIR_OUTPUT            = opt$output


# deseq_file = '' # either '' or the tsv with these columns:
## baseMean	log2FoldChange	lfcSE	stat	pvalue	padj

if (is.na(deseq_file) || deseq_file == "NA") {
  deseq_file <- ""
}

# Set up =========================================================
setwd(repo)
source('Scripts/conf/config_RNAmars.R')


# TRAINING RBPS DEFINITION ================
if (cell_line_target == "HepG2"){
  rbps = c("HNRNPC", "HNRNPK", "HNRNPU", "NCBP2", "PRPF8", "PTBP1", "QKI", "RBFOX2", "RBM22", "SF3A3", "SF3B4", "SRSF1", "U2AF1", "U2AF2", "UCHL5")
} else if (cell_line_target == "K562"){
  rbps = c("AGGF1", "EFTUD2", "FXR1", "HNRNPU", "PRPF8", "PTBP1", "PUS1", "RBM15", "SF3B4", "SRSF1", "TARDBP", "U2AF1", "U2AF2")
}


## Optimal set of params
OPT    = read.csv(file.path(repo, "Tables/RNAmotifs_optimal_parameters.csv"), stringsAsFactors = FALSE)
PARAMS = subset(OPT, true_rbp %in% rbps)
colnames(PARAMS)[colnames(PARAMS)=='true_rbp']<- 'pr'

## Peak
PEAK = readRDS(paste0("Rdata/",cell_line_target, "_binding_profile_PEAK_normalized.rds"))


# Splicing maps  ================

message(noquote("\n[*] Creating splicing maps using RBP-specific optimized parameters..."))
all_splicingMaps = do.call(rbind,lapply(unique(OPT$params), function(y){
  params = unlist(str_split(y, '_'))[c(2,4)]
 
  path_folder   = paste0(DIR_RESULTS_RNAMOTIFS, name, '_',y)
  
  if(length(grep("^All_group_enh-",list.files(path=path_folder),value=T))==0){return()}
  
  e_spl_map  = read.csv(paste0(path_folder,"/",grep("^All_group_enh-", list.files(path=path_folder),value=T)[1]))[,1]
  ES_spl_map = read.csv(paste0(path_folder,"/",grep("^All_group_both-",list.files(path=path_folder),value=T)[1]))[,1]
  
  if(length(grep("^Sil_group_sil-",list.files(path=path_folder),value=T))!=0){
    SIL_spl_map = read.csv(paste0(path_folder,"/",grep("^Sil_group_sil-", list.files(path=path_folder),value=T)[1]))[,1]
  } else {
    SIL_spl_map = ES_spl_map*0
  }
  
  if(length(grep("^Enh_group_enh-",list.files(path=path_folder),value=T))!=0){
    ENH_spl_map = read.csv(paste0(path_folder,"/",grep("^Enh_group_enh-", list.files(path=path_folder),value=T)[1]))[,1]
  } else {
    ENH_spl_map = ES_spl_map*0
  }
  
  tmp_col                 = e_spl_map/ES_spl_map
  tmp_col[is.na(tmp_col)] = 0
  
  return(suppressWarnings(cbind(pr         = name,
                                hw         = params[1],
                                ew         = params[2],
                                rbp        = name,
                                pos        = 1:length(SIL_spl_map),
                                height_enh = ENH_spl_map, 
                                height_sil = SIL_spl_map, 
                                height     = ES_spl_map,
                                col        = tmp_col)))
}))


all_splicingMaps            = as.data.frame(all_splicingMaps)
all_splicingMaps$pos        = as.numeric(all_splicingMaps$pos)
all_splicingMaps$height_enh = as.numeric(all_splicingMaps$height_enh)
all_splicingMaps$height_sil = as.numeric(all_splicingMaps$height_sil)

message(noquote("[*] RNAmotifs splicing maps created and ready for final heatmap."))

# DEseq2 ===============================
message(noquote("\n[*] Importing DESeq2 differential genes..."))
if (file.exists(deseq_file) & deseq_file %like% '.tsv') {
  deseq     = read.delim(deseq_file)
  sub_deseq = subset(deseq, gene_name %in% rbps)
} else if (file.exists(deseq_file) & deseq_file %like% '.rds') {
  deseq     = as.data.frame(readRDS(deseq_file))
  sub_deseq = subset(deseq, gene_name %in% rbps)
} else {
  # baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
  sub_deseq = data.frame(log2FoldChange = rep(0, length(rbps)), padj =  rep(1, length(rbps)), gene_name =rbps)
}



# Import scores ===================
srr_sil         = readRDS( paste0(DIR_OUTPUT, list.files(DIR_OUTPUT, pattern = 'SCORE1_sil') ))
profile_sim_sil = readRDS( paste0(DIR_OUTPUT, list.files(DIR_OUTPUT, pattern = 'SCORE2_sil') ))
srr_enh         = readRDS( paste0(DIR_OUTPUT, list.files(DIR_OUTPUT, pattern = 'SCORE1_enh') ))
profile_sim_enh = readRDS( paste0(DIR_OUTPUT, list.files(DIR_OUTPUT, pattern = 'SCORE2_enh') ))

ENH = lapply(names(profile_sim_enh), function(params) {
  srr_mat = srr_enh[[params]]
  prof_mat = profile_sim_enh[[params]]
  res = prof_mat*srr_mat
  return(res)
})
names(ENH) = names(profile_sim_enh)
ENH = list(ENH)
names(ENH) = name


SIL = lapply(names(profile_sim_sil), function(params) {
  srr_mat = srr_sil[[params]]
  prof_mat = profile_sim_sil[[params]]
  res = prof_mat*srr_mat
  return(res)
})
names(SIL) = names(profile_sim_sil)
SIL = list(SIL)
names(SIL) = name



info = do.call(rbind, lapply(c(5, 15, 25, 35), function(hw){
  do.call(rbind, lapply(c(30, 50, 100, 200, 300), function(ew){  
    params       = paste0('hw_',hw,'_ew_',ew)
 
    summary_file = paste0(DIR_OUTPUT, '/summary_results_', name,'_',params,'.rds')
    if(!file.exists(summary_file)) {return()} 
    else {
      res = readRDS(summary_file)
    }
    res$key = params
    return(res)
  }))
}))


for ( typeAS in c('enh','sil')) {
  if (typeAS =='sil') {MAT = SIL[[name]]} else {MAT = ENH[[name]]}
      
  sub_res_df = subset(info, key %in% unique(PARAMS$params))
  if (typeAS =='sil') {
    sub_res_df = sub_res_df[which(sapply(1:nrow(sub_res_df), function(x){"Silenced" %in% sub_res_df[x,] | "Both" %in% sub_res_df[x,]})),]
  } else {
    sub_res_df = sub_res_df[which(sapply(1:nrow(sub_res_df), function(x){"Enhanced" %in% sub_res_df[x,] | "Both" %in% sub_res_df[x,]})),]
  }
    
  if (dim(sub_res_df)[1]!=0) {
    res = plot_association_heatmap(MAT            = MAT, 
                                   PARAMS         = PARAMS,
                                   PEAK           = PEAK,
                                   sub_res_df     = sub_res_df,
                                   typeAS         = typeAS,
                                   input_name     = name,
                                   ranking_method = "rowmean",
                                   sub_deseq      = sub_deseq
    )
    p           = res[[1]]
    mat         = res[[2]]
    ranking     = res[[3]]
    tet_score   = res[[4]]
    pval_legend = res[[5]]
    
    message(noquote("\n[*] Saving heatmap as rds file..."))
    saveRDS(mat, paste0(DIR_OUTPUT, '/final_mat_',typeAS,'.rds'))
    saveRDS(tet_score, paste0(DIR_OUTPUT, '/final_mat_tet_score_',typeAS,'.rds'))
    saveRDS(ranking, paste0(DIR_OUTPUT, '/final_mat_prot_score_',typeAS,'.rds'))
    
    if ( all(dim(mat) ==c(1,1))) {
      ht_dim = c(10,10)
    } else {
      ht_dim = getdim(draw(p))
    }  
    # Save plot
    message(noquote(paste0("[*] Printing heatmap for ", ifelse(typeAS == 'sil', "silenced", "enhanced"), " exons...")))
    figure_final_path = paste0(DIR_OUTPUT, cell_line_target,'_', name, '_', typeAS,'.pdf')
    pdf(figure_final_path, width = ht_dim[1]+5, height = ht_dim[2])
    draw(p, annotation_legend_list = list(pval_legend))
    add_splicing_maps_to_final_ht(sub_splicingMaps = all_splicingMaps, mat = mat, params = PARAMS)
    dev.off()
  } else {next}
}

