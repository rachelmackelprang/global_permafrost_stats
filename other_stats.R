library(tidyverse)
library(readxl)
library(phyloseq)



set.seed(42)

###################################
######import metadata
###################################

meta<-read_excel("metadata.xlsx")

#change anything of type 'chr' to a factor. Make the first column into row names
meta<-meta %>% 
  as.data.frame() %>% 
  mutate_if(sapply(meta, is.character), as.factor) %>% 
  column_to_rownames('SampleID')


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
kegg_physeq.f<-filter_taxa(physeq, function(x) sum(x>10) > (0.1 * length(x)), prune=TRUE)

#get relative abundance
kegg_physeq.f.ra<-transform_sample_counts(physeq.f, function(x) x*100/sum(x))



###############################
#get diversity stats###########
###############################

kegg_shannon<-estimate_richness(kegg_physeq.f, measures = "Shannon")


#is there a relationship between diversity and pH? Looks for linear and polynomial. 
shannon_fit1<-lm(kegg_shannon$Shannon ~ meta$pH)
shannon_fit2<-lm(kegg_shannon$Shannon ~ meta$pH + I(meta$pH^2))

anova(shannon_fit1)
anova(shannon_fit2)

plot(meta$pH,kegg_shannon$Shannon)
abline(shannon_fit1, lwd=3, col="red")
abline(shannon_fit2, lwd=3, col="blue")
