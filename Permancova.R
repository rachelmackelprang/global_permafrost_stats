#Permancova analysis

library(tidyverse)
library(readxl)
library(Amelia)
library(factoextra)
library(ggfortify)
library(phyloseq)
library(vegan)


###################################
######import and transform metadata
###################################

meta<-read_excel("~/Dropbox/GitHub/global_permafrost_stats//metadata.xlsx")

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
meta$pH<-meta$pH^(1/3)

################################
#import and manipulate KEGG data
################################

#import KEGG data (counts)
kegg_counts<-read_excel("~/Dropbox/GitHub/global_permafrost_stats//KEGG_allsamples_counts_newnames.xlsx")
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
imputed<-amelia(meta, m=5, idvars=c("Site", "Region1", "Type", "Origin", "Region2", "Period1", "Longitude", "Period2"))
imputed_matrix<-as.data.frame(imputed$imputations)

#normalize by subtracting the mean and dividing by the standard deviation. This also selects only one
#of the imputations output by amelia 
imputed_norm_matrix<-lapply(imputed_matrix[c(8,10:19)], function(x) c(scale(x)))

#select only continuous variables that we want to analyze (we are excluding longitude)
imputed_norm_cont<-imputed_matrix[c(8,10:19)]

##############################
#Correlation matrix and PCA##
##############################

#create a correlation matrix
imputed_norm_cont_corr<-cor(imputed_norm_cont)

#PCA of continuous metadata variables
meta_pca<-prcomp(imputed_norm_cont_corr)

#Extract Z-scores from PC1 and PC2 and use in permancova?



##############################
#PERMANCOVA###################
##############################


#square root transform kegg relative abundance data
kegg_sqrt<-sqrt(kegg_trimmed_ra[,1:5504])

#get bray-curtis dissimilarities
kegg_bray<-vegdist(kegg_sqrt, "bray")


