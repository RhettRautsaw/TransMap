#!/usr/bin/env Rscript
# Copyright 2022 Rhett M. Rautsaw
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#=====================#
##### DESCRIPTION #####
#=====================#
#
# TransMap.R is designed to perform cross-species transcriptome expression homology graphing.
# 
#=============#
##### USE #####
#=============#
#
# TransMap.R -t transcriptomes -s species_tree.newick -e expression -m metadata.txt -c 8
#
#======================#
##### REQUIREMENTS #####
#======================#

## Packages
packages = c("argparse", "tidyverse", "treeio","ape", "scales", "ggpubr", "patchwork", "umap")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  print("Installing required packages")
  install.packages(packages[!installed_packages], repos="https://cloud.r-project.org")
}


#=========================#
##### SETUP ARGUMENTS #####
#=========================#

suppressPackageStartupMessages(library("argparse"))
cpu<-parallel::detectCores()

# create parser object
parser <- ArgumentParser(description="TransMap.R is designed to perform cross-species transcriptome expression homology graphing.")

# specify our desired options 
# by default ArgumentParser will add an help option
parser$add_argument("-t", default="transcriptomes", 
                    help="Folder containing fastas for each species transcriptome named as Genus_species.fasta [default: \"%(default)s\"]")
parser$add_argument("-s", default="species_tree.newick", 
                    help="Species tree for creation of homology graph [default: \"%(default)s\"]")
parser$add_argument("-q", default=50,  
                    help="Query hsp coverage for BLAST searches [default: \"%(default)s\"]")
parser$add_argument("-e", default="expression", 
                    help="Folder containing tab-delimited TPM values for individuals of each species named as Genus_species.txt [default: \"%(default)s\"]")
parser$add_argument("-m", default="metadata.txt",
                    help="Tab-delimited file of metadata for each sample with the first column matching the sample names in the expression files. [default: \"%(default)s\"]")
parser$add_argument("-c", type="integer", default=cpu, metavar="number",
                    help="Number of cpus to speed up some steps. [default: \"%(default)s\"]")


#========================#
##### READ ARGUMENTS #####
#========================#

args <- parser$parse_args()
trans <- normalizePath(args$t)
stree <- normalizePath(args$s)
exprs <- normalizePath(args$e)
metad <- normalizePath(args$m)


cat("Starting TransMap\n")
cat(paste0("\t Transcriptomes:\t", trans,"\n"))
cat(paste0("\t Species Tree:\t", stree,"\n"))
cat(paste0("\t Expression Data:\t", exprs,"\n"))
cat(paste0("\t Metadata:\t", metad,"\n"))
cat(paste0("\t CPUs:\t", args$c,"\n"))

#=======================#
##### LOAD PACKAGES #####
#=======================#

cat(paste0("\n", Sys.time(), " ::: Loading packages and preparing directory :::\n"))
suppressWarnings(invisible(suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))))

#=========================#
##### SETUP FUNCTIONS #####
#=========================#

`%!in%` = Negate(`%in%`)

fill_mat = function(cols, m){
  res = `dimnames<-`(m[match(cols,rownames(m)), match(cols,colnames(m))], list(cols, cols))
  ifelse(is.na(res), 0, res)
}

samap_hom<-function(map, exprs, species="Crotalus_scutulatus_scutulatus"){
  # Create empty dataframe for transformed values
  results = data.frame()
  # Read Expression Data
  df = as.data.frame(t(read.table(paste0(exprs,"/",species[1],".txt"),check.names = F, header=1, row.names=1)))
  rownames(df)<-gsub("\\.","-",rownames(df))
  # Log Transform
  df1 = log2(df+1)
  # Pivot Longer
  df0 = df1 %>% rownames_to_column("sample") %>% gather("gene","r",2:ncol(.))
  # Bind to Results
  results = rbind.data.frame(results,df0)
  # Get Map results for missing genes
  map2 = map[colnames(map) %in% colnames(df1),colnames(map) %!in% colnames(df1)]
  # Calculate linear operators to impute each missing gene's expression based on similarlity of genes
  su = colSums(map2)
  lo = map2 %*% diag(1/su)
  colnames(lo) = colnames(map2)
  lo[is.nan(lo)] = 0
  # Multiply known expression by linear operators and sum to get the expression for missing genes
  df2 = as.matrix(df1[,rownames(lo)])%*%lo
  # Make dataframe long rather than wide
  df2 = as.data.frame(df2) %>% rownames_to_column("sample") %>% gather("gene","r",2:ncol(.))
  # Bind the data to the full dataframe
  results = rbind.data.frame(results,df2)
  # Return Results
  return(results)
}

#========================#
##### BEGIN TRANSMAP #####
#========================#

if(!file.exists("blastmap.txt")){
  #### PREPARING FOR BLAST SEARCHES ####
  cat(paste0("\n", Sys.time(), " ::: Reading Species Tree :::\n"))
  
  tree<-read.tree(stree)
  cat(paste0("\n", Sys.time(), " ::: Creating Seq Identity Cutoffs from Phylogenetic Distance :::\n"))
  x=cophenetic(tree)
  x=as.data.frame(x) %>% rownames_to_column("sp1") %>% gather("sp2", "dist", 2:ncol(.))
  x$dist <- round(rescale(x$dist, to = c(97, 72)))
  
  cat(paste0("\n", Sys.time(), " ::: Creating Species Comparisons :::\n"))
  species<-gsub(".fasta","",list.files(trans,".fasta$"))
  write.table(species,"species.list",sep="\t", quote=F,row.names = F, col.names=F)
  comps = expand.grid(species,species) %>% select(sp1=Var1, sp2=Var2) %>% mutate(comp=paste0(sp1,"-",sp2)) %>% mutate(same=sp1==sp2) %>% filter(same=="FALSE") %>% select(-same)
  comps<-left_join(comps,x)
  write.table(comps,"comparisons.txt", sep="\t", quote=F,row.names = F, col.names=F)
  
  
  #### PAIRWISE BLAST SEARCHES ####
  system("mkdir maps")
  
  cat(paste0("\n", Sys.time(), " ::: Creating BLAST Databases :::\n"))
  system(paste0("parallel -a species.list -j ", args$c, " --bar 'cd ",trans, "; makeblastdb -in {}.fasta -dbtype nucl > /dev/null'"))
  
  cat(paste0("\n", Sys.time(), " ::: Running BLAST Searches :::\n"))
  system(paste0("parallel -a comparisons.txt -j ", args$c, " --bar --colsep '\t' 'blastn -query ", trans, "/{1}.fasta -db ", trans, "/{2}.fasta -outfmt 6 -out maps/{3}.txt -num_threads 1 -max_hsps 1 -perc_identity {4} -qcov_hsp_perc ", args$q, " -max_target_seqs 10000 -evalue 1e-10'"))
  
  cat(paste0("\n", Sys.time(), " ::: Concatenating BLAST Results :::\n"))
  system("ls maps | parallel -j 1 --bar 'cat maps/{} >> blastmap.txt'")
}else{
  species=readLines("species.list")
}

#### HOMOLOGY GRAPH PREPARATION ####
cat(paste0("\n", Sys.time(), " ::: Creating Homology Graph :::\n"))
map = read_tsv("blastmap.txt", col_names=FALSE, col_types = cols())
map = as.matrix(map %>% select(X1,X2,X12) %>% spread(X2,X12) %>% replace(is.na(.), 0) %>% column_to_rownames("X1"))
# Add Missing Rows/Columns as 0
cols = sort(union(colnames(map), rownames(map)))
map = fill_mat(cols, map)
## Artificially inflate bitscores so that similar things are more heavily weighted later
map2 = map^10


#### HOMOLOGY GRAPH IMPUTATION ####
if(!file.exists("imputed_results.txt")){
  imputed_results = data.frame()
  pb = txtProgressBar(min = 0, max = length(species), initial = 0) 
  cat(paste0("\n", Sys.time(), " ::: Imputing Expression for Each Species based on Homology Graph :::\n"))
  for(i in 1:length(species)){
    setTxtProgressBar(pb,i)
    res = samap_hom(map2, exprs, species[i])
    imputed_results = rbind.data.frame(imputed_results,res)
  }
  close(pb)
  rm(res, pb)
  
  cat(paste0("\n", Sys.time(), " ::: Writing Results :::\n"))
  # Pivot Wider
  imputed_results = imputed_results %>% spread(gene, r) %>% replace(is.na(.), 0) #%>% #mutate(sample=gsub("\\.","-",sample)) 
    #column_to_rownames("sample")
  
  write.table(imputed_results, "imputed_results.txt", sep="\t", quote=F, row.names = F)
  
}else{
  imputed_results = read.table("imputed_results.txt", sep='\t', header = T, check.names = F)
}

cat(paste0("\n", Sys.time(), " ::: Reading Metadata :::\n"))
metadata<-read.table(metad, sep="\t", header=1, check.names = F)

cat(paste0("\n", Sys.time(), " ::: Performing PCA :::\n"))
pca<-prcomp(imputed_results[,2:ncol(imputed_results)])
pca_df<-cbind.data.frame(sample=imputed_results[,1], pca$x)
results<-merge(metadata, pca_df[,1:6], by="sample")

cat(paste0("\n", Sys.time(), " ::: Performing UMAP :::\n"))
umap_df = umap(pca_df[,2:ncol(pca_df)], random_state=2022, n_neighbors=15)
umap_df2=umap_df$layout
colnames(umap_df2)<-c("UMAP1", "UMAP2")
results<-cbind.data.frame(results, umap_df2)

cat(paste0("\n", Sys.time(), " ::: Performing VAE :::\n"))
system("conda activate popvae; ./TransVAE.py --infile imputed_results.txt --seed 42 --out TransVAE_out; conda deactivate")
vae_df = read.table("TransVAE_out_latent_coords.txt", sep='\t', header=T) %>% select(sample=sampleID, VAE1=mean1, VAE2=mean2)
results<-merge(results, vae_df, by="sample")

cat(paste0("\n", Sys.time(), " ::: Plotting/Writing Results :::\n"))
write.table(results, "imputed_latent.txt", sep="\t", quote=F, row.names = F)
A = ggscatter(results, "PC1", "PC2", color="species", label="sample") + rremove("legend")
B = ggscatter(results, "UMAP1", "UMAP2", color="species", label="sample") + rremove("legend")
C = ggscatter(results, "VAE1", "VAE2", color="species", label="sample") + rremove("legend")

cat(paste0("\n", Sys.time(), " ::: Plotting Results :::\n"))
png("TransMap_plot.png", width=3240, height=1080)
A+B+C
dev.off()


#imputed_results2 = right_join(metadata,imputed_results %>% rownames_to_column("sample"))

# write.table(imputed_results2, "imputed_results.txt", sep="\t", quote=F, row.names = F)

# scripts/tranvae.py --infile data/expression.txt --seed 42 --out out/expression_out

# cat(paste0("\n", Sys.time(), " ::: PCA :::\n"))
# pca<-prcomp(imputed_results[,1:ncol(imputed_results2)])
# pca_df<-cbind.data.frame(imputed_results2[,1:3], pca$x)
# 
# A<-ggscatter(pca_df, "PC1", "PC2", color="species", label="sample") + rremove("legend")
# 
# cat(paste0("\n", Sys.time(), " ::: UMAP :::\n"))
# umap_df = umap(pca_df[,4:ncol(pca_df)], random_state=2022, n_neighbors=15)
# umap_df2=umap_df$layout
# colnames(umap_df2)<-c("UMAP1", "UMAP2")
# umap_df2<-cbind.data.frame(imputed_results2[,1:3], umap_df2)
# B<-ggscatter(umap_df2, "UMAP1", "UMAP2", color="species", label="sample") + rremove("legend")
