#Kolmogrov-Smirnov test

CD_stats <- read.csv("../../output/all_CD_priority_scores.csv")
anno_mods_df <- CD_stats[,c("Module","HMDB.ID","m.z","Covaries_with_standard")] %>%
  filter(Covaries_with_standard==1)

# Calculate ∆mzs in real modules
real_dmzs <- NULL
for(m in annotated_modules){
  df <- anno_mods_df[anno_mods_df$Module == m,]
  dmzs <- as.vector(sapply(unique(df[df$HMDB.ID == "","m.z"]), function(x) x - unique(df[df$HMDB.ID != "","m.z"])))
  module_wise_dmzs <- as.data.frame(cbind(m, dmzs))
  real_dmzs <- rbind(real_dmzs,module_wise_dmzs)
}
colnames(real_dmzs) <- c("module","dmz")
real_dmzs$type <- "real"


# Shuffle
all_ks <- NULL
for(n in seq(1:1000)){
  random_dmz <- NULL
  anno_mods_df$shuffled <- sample(anno_mods_df$Module, size=nrow(anno_mods_df), replace=FALSE)
  has_std <- unique(anno_mods_df[anno_mods_df$HMDB.ID != "","shuffled"])
  for(i in has_std){
    df <- anno_mods_df[anno_mods_df$shuffled == i,]
    dmzs <- as.vector(sapply(unique(df[df$HMDB.ID == "","m.z"]), function(x) x - unique(df[df$HMDB.ID != "","m.z"])))
    random_dmz <- append(random_dmz, dmzs)
  }
  ks_result <- ks.test(real_dmzs$dmz, random_dmz) # two-sided KS test
  output <- cbind(as.numeric(ks_result$statistic), as.numeric(ks_result$p.value))
  all_ks <- rbind(all_ks, output)
  rm(random_dmz)
}
all_ks <- as.data.frame(all_ks) %>% 
  rename(stat = 1, pval = 2)

