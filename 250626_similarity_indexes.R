################################################################################
### SIMILARITY INDEXES ###
################################################################################


### R script for valuation of different aa seq if they're biased towards 
## some lineage or not

## corresponds to the lab notebook 250626_similarity_indexes_boehms_paper_table

getwd()
## load the data
finalResult <- read_csv("250428_result_filtered_table_3a.csv")
countTable <- read_csv("250627_count_table2_seq_samples.csv")




## add metadata to last few rows
for (i in 2:ncol(countTable)){
  countTable %>% mutate_at(vars(colnames(countTable)[i]),as.numeric)
}

## creating a table with metadata info that can be added to the end of the 
## count table
samples_names <- unique(finalResult$Sample.bio.name)

df_unique <- finalResult %>%
  distinct(Sample.bio.name, .keep_all = TRUE)
df <- data.frame(ID = c(1:22),df_unique[,c(29,26,25,27)])
dfT <- t(df)
dfT <- data.frame(dfT)
dfT <- rownames_to_column(dfT,"aa_seq")
colnames(dfT)[2:ncol(dfT)] <- samples_names
### merging the metadata and the count table ###
#countTable <- count_Table1
#rownames(countTable) <- count_Table1$aa_seq
countTableFinal <- rbind(countTable,dfT)
colnames(dfT)
colnames(countTable)

excludeColList <- c()
for (col in 1:ncol(dfT)){
  if (dfT[[3,col]] == "T"){
    excludeColList[length(excludeColList) +1] = col
  }
}
countTableSpleen <- countTable[,-excludeColList]


# the contingency table together with columns for results
similarityIndexesTable <- data.frame(
  aa_seq = countTable$aa_seq,
  CD4_present = NA,
  CD8_present = NA,
  log_odds_ratio = NA,
  fisher_p = NA,
  jaccard_distance = NA
  #stringsAsFactors = FALSE
)



# separate the samples according to the cell type
group <- as.character(dfT[4,-excludeColList ])  
cd4_samples <- which(group == "CD4")
cd8_samples <- which(group == "CD8")




for (i in 1:nrow(countTableSpleen)) { 
  oneSeq <- countTableSpleen[i,] # the current sequence
  a <- sum(oneSeq[cd4_samples])  
  c <- sum(oneSeq[cd8_samples]) 
  b <- length(cd4_samples) - a
  d <- length(cd8_samples) - c

  print(i)
  # contingency table for Fisher and log odds
  contingTable <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(Group = c("CD4", "CD8"), Status = c("Present", "Absent")))
  print(contingTable)
  test <- fisher.test(contingTable)
  # Log-odds ratio
  # Add 0.5 correction to avoid log(0) if needed
  log_odds <- log((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)))
  
  # jaccard
  jaccard <- sum(cd4_samples & cd8_samples)/sum(cd4_samples | cd8_samples)
  jaccard_distance <- 1-jaccard
  print(jaccard_distance)
  
    
  # Store results
  similarityIndexesTable[i, "CD4_present"] <- a
  similarityIndexesTable[i, "CD8_present"] <- c
  similarityIndexesTable[i, "log_odds_ratio"] <- log_odds
  similarityIndexesTable[i, "fisher_p"] <- test$p.value
  similarityIndexesTable[i, "jaccard_distance"] <- jaccard_distance
}

print(similarityIndexesTable)
write_csv(similarityIndexesTable, "250701_similarity_indexes.csv")
write_csv(similarityIndexesTable, "250701_similarity_indexes_just_spleen.csv")


similarityIndexesTable <- similarityIndexesTable %>% mutate(diffexpr = ifelse(similarityIndexesTable$log_odds_ratio > 0.6 & similarityIndexesTable$fisher_p < 0.05, 
                                                    "CD4 biased",
                                                    ifelse(similarityIndexesTable$log_odds_ratio < -0.6 & similarityIndexesTable$fisher_p < 0.05,"CD8 biased","NO")))

library(ggrepel)
ggplot(data=similarityIndexesTable, aes(x=similarityIndexesTable$log_odds_ratio, 
                                        y=-log10(similarityIndexesTable$fisher_p), 
                                        col = diffexpr,
                                        label = similarityIndexesTable$aa_seq)) + 
  geom_point() + 
  theme_classic() +
  scale_color_manual(values=c("red","blue", "black")) +
  geom_text_repel()

ggsave("250701_volcano_plot_aaseq_bias1.png")

summary(fisher.test(contingTable))


## create contingency table ##

contingTable <- data.frame(present = c(-1,-1),absent = c(-1,-1))
rownames(contingTable) <- c("CD4","CD8")


library(dplyr)

#write_csv(countTableFinal, "250626_count_table_with_metadata.csv")
#write.xlsx(countTableFinal, "250626_count_table_with_metadata.xlsx")


countTable %>% mutate(across(all_of(names(countTable)[2:ncol(countTable)])),as.numeric)
