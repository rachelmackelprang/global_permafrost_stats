#Permancova analysis

library(tidyverse)
library(readxl)
library(Amelia)
library(factoextra)
library(ggfortify)
library(phyloseq)
library(vegan)


set.seed(42)


###################################
######import and transform metadata
###################################

meta<-read_excel("metadata.xlsx")

#change anything of type 'chr' to a factor. Make the first column into row names
meta<-meta %>% 
  as.data.frame() %>% 
  mutate_if(sapply(meta, is.character), as.factor) %>% 
  column_to_rownames('SampleID')


#transform metadata. All continuous variables except pH were skewed right so they are log (y+1) transformed. 
#pH was skewed left. Box-Cox estimation showed a cubic transformation was best for pH. 

meta$Latitude<-log10(meta$Latitude+1)
meta$Elevation_m<-log10(meta$Elevation_m+1)
meta$KYBP<-log10(meta$KYBP+1)
meta$Depth_m<-log10(meta$Depth_m+1)
meta$TotalC_percent<-log10(meta$TotalC_percent+1)
meta$TotalN_percent<-log10(meta$TotalN_percent+1)
meta$CN<-log10(meta$CN+1)
meta$OrganicC_percent<-log10(meta$OrganicC_percent+1)
meta$`water content gwat per gsoil`<-log10(meta$`water content gwat per gsoil`+1)
meta$electrical_cond<-log10(meta$electrical_cond+1)
meta$pH<-meta$pH^(3)

################################
#import and manipulate KEGG data
################################

#import KEGG data (counts)
kegg_counts<-read_excel("KEGG_allsamples_counts_newnames.xlsx")
kegg_counts<-kegg_counts %>% as.data.frame() %>% column_to_rownames("KEGG")

#create phyloseq object to perform filtering and change to relative abundance
OTU<-otu_table(kegg_counts, taxa_are_rows = TRUE)
physeq<-phyloseq(OTU, meta)

#prune low abundance genes (keep only those observed more than 10 times in at least 10% of samples)
physeq.f<-filter_taxa(physeq, function(x) sum(x>10) > (0.1 * length(x)), prune=TRUE)

#get relative abundance
physeq.f.ra<-transform_sample_counts(physeq.f, function(x) x*100/sum(x))

#export trimmed relative abundance data
kegg_trimmed_ra<-as.data.frame(otu_table(physeq.f.ra))
kegg_trimmed_ra<-t(kegg_trimmed_ra)

###############################
#imputations###################
###############################

#impute missing values of continuous variables (except longitude)
imputed<-amelia(meta, m=1, idvars=c("Site", "Region1", "Type", "Origin", "Region2", "Period1", "Longitude", "Period2"))
imputed_df<-as.data.frame(imputed$imputations)

#normalize by subtracting the mean and dividing by the standard deviation. This also selects only one
#of the imputations output by amelia. Note we are excluding longitude

imputed_norm_cont<-scale(imputed_df[c(8,10:19)])

##############################
#Correlation matrix and PCA###
##############################

#create a correlation matrix
imputed_norm_cont_corr<-cor(imputed_norm_cont)

#PCA of continuous metadata variables from normalized data
meta_norm_pca<-princomp(imputed_norm_cont, scores = TRUE)

#get z-scores 
meta_norm_zscores<-meta_norm_pca$scores
meta_norm_zscores<-as.data.frame(meta_norm_zscores)

#PCA of continuous metadata variables based on correlation matrix
#meta_corr_pca<-princomp(imputed_norm_cont_corr, scores=TRUE)

#get z-scores from correlation matrix PCA
#meta_corr_zscores<-meta_corr_pca$scores




#############################
#data wrangling. 
#############################

#add site category and region1 to PC dataframe
meta_norm_zscores<-cbind(meta_norm_zscores, meta$Site)
meta_norm_zscores<-cbind(meta_norm_zscores, meta$Region1)

#change column name of Site from "meta$Site" to Site
colnames(meta_norm_zscores)[colnames(meta_norm_zscores) == "meta$Site"] <- "Site"
colnames(meta_norm_zscores)[colnames(meta_norm_zscores) == "meta$Region1"] <- "Region1"

#Find only samples that overlap between KEGG data and metadata
kegg_meta_common_ids<-intersect(rownames(meta_norm_zscores), rownames(kegg_trimmed_ra))
kegg_trimmed_ra_common<-kegg_trimmed_ra[kegg_meta_common_ids,]
meta_norm_zscores_common<-meta_norm_zscores[kegg_meta_common_ids,]



##############################
#PERMANCOVA###################
##############################


#square root transform kegg relative abundance data
kegg_sqrt<-sqrt(kegg_trimmed_ra_common[,1:5504])

#get bray-curtis dissimilarities
kegg_bray<-vegdist(kegg_sqrt, "bray")



#permanova

#Comp.1 and Comp.2 are the z-scores from the PCA
kegg_PC1_PC2_adonis<-adonis(formula=kegg_bray~Site*Comp.1+Site*Comp.2, data=meta_norm_zscores_common, permutations=1000)

kegg_PC1_PC2_adonis

