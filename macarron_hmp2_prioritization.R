library(tidyverse)

# Source MACARRoN functions
source("src/prepInput.R")
source("src/makeDisMat.R")
source("src/findMacMod.R")
source("src/calAVA.R")
source("src/prioritize.R")

# This workflow applies MACARRoN to the HMP2 metabolomics data to prioritize
# metabolic features in IBD subtypes - CD-dysbiosis and UC

# Read HMP2 metabolic feature abundances and name features
abun <- read.csv("input/hmp2_intensities.csv")
rownames(abun) <- paste0("F",1:nrow(abun))
# Read HMP2 metabolic feature annotations and rearrange dataframe
anno <- read.csv("input/hmp2_annotations.csv")
rownames(anno) <- paste0("F",1:nrow(anno))
anno <- anno[,c(4,5,3,2,1,6:ncol(anno))] %>%
# Override all standards as primary features and then keep only primary features
anno[anno$Metabolite != "", "prim_feature"] <- "primary"
anno <- anno[anno$prim_feature == "primary",]
abun <- abun[rownames(anno),]
# Read HMP2 metadata and rearrange dataframe
meta <- read.csv("input/hmp2_metadata.csv")
meta <- meta[,c(1,6,2,4,3,5)] %>%
  column_to_rownames("sample")
# Read chemical taxonomy file
chem_tax <- read.csv("input/hmdb_taxonomy.csv")

#--------------------------
# Identification of modules
#--------------------------

# Create summarized experiment object
hmp2_mbx <- prepInput(abun, anno, meta)
# Phenotype-stratified correlations and distance matrix using best observed 
# bicor between a pair of features
dist_mat <- makeDisMat(hmp2_mbx)
# Identify modules
hmp2_modules <- findMacMod(hmp2_mbx, dist_mat, chem_tax)
module_assignments <- hmp2_modules[[1]]
# Calculate AVA of metabolic features
hmp2_ava <- calAVA(hmp2_mbx, hmp2_modules)


# Read effect sizes and q-values 
# q-values are calculated using Maaslin2 (https://bioconductor.org/packages/Maaslin2/) and
# effect-sizes are difference in means of log2 abundances in the two phenotypes
es_CD <- read.csv("input/primary_prevalent_features/effect_sizes_CDdys_vs_CDnondys.csv", row.names = 1)
es_UC <- read.csv("input/primary_prevalent_features/effect_sizes_UC_vs_nonIBD.csv", row.names = 1)
qval_CD <- read.csv("input/primary_prevalent_features/qvalues_CDdys.csv", row.names = 1)
qval_UC <- read.csv("input/primary_prevalent_features/qvalues_UC.csv", row.names = 1)
#--------------------------
# Prioritization in CD
#--------------------------
prioritized_in_CD <- prioritize(hmp2_mbx, 
                                hmp2_modules,
                                hmp2_ava,
                                qval_CD,
                                es_CD)
all_CD_priority_scores <- prioritized_in_CD[[1]] # Sheet 1 in Dataset EV7

#--------------------------
# Prioritization in UC
#--------------------------
prioritized_in_UC <- prioritize(hmp2_mbx, 
                                hmp2_modules,
                                hmp2_ava,
                                qval_UC,
                                es_UC)

all_UC_priority_scores <- prioritized_in_UC[[1]] # Sheet 2 in Dataset EV7

