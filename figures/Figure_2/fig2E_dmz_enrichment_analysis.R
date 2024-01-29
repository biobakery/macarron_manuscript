CD_stats <- read.csv("../../output/all_CD_priority_scores.csv")
anno_mods_df <- CD_stats[,c("Module","HMDB.ID","m.z","Covaries_with_standard")] %>%
  filter(Covaries_with_standard==1)

#Enrichment of ∆mzs
real_dmzs <- NULL
for(m in annotated_modules){
  df <- anno_mods_df[anno_mods_df$Module == m,]
  dmzs <- as.vector(sapply(unique(df[df$HMDB.ID == "","m.z"]), function(x) x - unique(df[df$HMDB.ID != "","m.z"])))
  real_dmzs <- append(real_dmzs,dmzs)
}

# Find frequencies of positive and negative ∆mzs
real_dmzs_freq <- as.data.frame(table(round(real_dmzs))) %>%
  rename(dmz = 1, Freq = 2) %>%
  filter(Freq >= 3) %>%
  mutate(dmz = as.numeric(as.character(dmz))) %>%
  filter(abs(dmz) > 1) %>%
  mutate(type = ifelse(dmz < 0, "neg", "pos"))

pos_dmzs <- real_dmzs_freq[real_dmzs_freq$type == "pos",c(1,2)]
neg_dmzs <- real_dmzs_freq[real_dmzs_freq$type == "neg",c(1,2)]
rownames(pos_dmzs) <- NULL
rownames(neg_dmzs) <- NULL

# Permutation test
empty_df1 <- as.data.frame(t(pos_dmzs$Freq)) # for positive random ∆mzs
empty_df2 <- as.data.frame(t(neg_dmzs$Freq)) # for negative random ∆mzs
rownames(pos_dmzs) <- pos_dmzs$dmz
rownames(neg_dmzs) <- neg_dmzs$dmz
for (n in 1:10000){
  random_dmz <- NULL
  anno_mods_df$shuffled <- sample(anno_mods_df$Module, size=nrow(anno_mods_df), replace=FALSE)
  has_std <- unique(anno_mods_df[anno_mods_df$HMDB.ID != "","shuffled"])
  for(i in has_std){
    df <- anno_mods_df[anno_mods_df$shuffled == i,]
    dmzs <- as.vector(sapply(unique(df[df$HMDB.ID == "","m.z"]), function(x) x - unique(df[df$HMDB.ID != "","m.z"])))
    random_dmz <- append(random_dmz, dmzs)
  }
  random_dmz_counts <- as.data.frame(table(round(random_dmz)))
  random_dmz_counts$Var1 <- as.numeric(as.character(random_dmz_counts$Var1))
  rownames(random_dmz_counts) <- random_dmz_counts$Var1
  new_counts <- random_dmz_counts[rownames(pos_dmzs),2]
  empty_df1 <- rbind(empty_df1, new_counts)
  new_counts <- random_dmz_counts[rownames(neg_dmzs),2]
  empty_df2 <- rbind(empty_df2, new_counts)
}
combined_df <- cbind(empty_df1, empty_df2)
combined_df[is.na(combined_df)] <- 0

# Calculate empirical p-value
emp_pval <- function(f){
  n = 10001
  r = length(combined_df[2:nrow(combined_df),f][which(combined_df[2:nrow(combined_df),f] >= combined_df[1,f])]) + 1
  pval <- r/n
}
real_dmzs_freq <- as.data.frame(rbind(pos_dmzs, neg_dmzs))
real_dmzs_freq$emp_results <- sapply(seq(1:nrow(real_dmzs_freq)), function(f) emp_pval(f))
real_dmzs_freq$emp_results <- as.numeric(as.character(real_dmzs_freq$emp_results))
real_dmzs_freq$adj_p <- p.adjust(real_dmzs_freq$emp_results, method = "BH")
real_dmzs_freq <- real_dmzs_freq[order(real_dmzs_freq$adj_p),]
write.csv(real_dmzs_freq, file="data/fig2E_permutation_dmz_results.csv", row.names = FALSE)