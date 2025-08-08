## Monte Carlo simulations for null hypothesis
countTable.3a <- read_csv(paste(getwd(),"results_table/MM.xp.3a_results","250702_commit_accuracy_non-unique_seq_spleen.csv",sep="/"))
countTable.3b <- read_csv(paste(getwd(),"results_table/MM.xp.3b_results","250714_3b_commit_accuracy_non-unique_seq_spleen.csv",sep="/"))




# for alpha 780 sequences
simulationsTable <- data.frame(ID = "",Mean_Accurancy = "",Weighted_Mean = "", 
                               Strictly_CD4 = "",Strictly_CD8 = "", Undetermined = "")
simulationsTable <- simulationsTable[-1,]
samples.3a <- countTable.3a[1:779,2:16]
group <- countTable.3a[783,2:16]
cd4_samples <- which(group == "CD4")
cd8_samples <- which(group == "CD8")


## function to work with the simulation table so we get the final mean value
getMeanAccurance <- function(dataTable, aa_seq_col){
  dataTable <- data.frame(dataTable)
  dataTable <- dataTable  %>%
    mutate(Occurrence = rowSums(dataTable)) %>%
    mutate(Occurr_CD4 = 0) %>%
    mutate(Occurr_CD8 = 0) %>%
    mutate(aa_seq = aa_seq_col)
  dataTable$Occurr_CD4 <- rowSums(dataTable[, cd4_samples])
  dataTable$Occurr_CD8 <- rowSums(dataTable[, cd8_samples])
  dataTable <- dataTable %>%
    mutate(Dominant_pheno = ifelse(Occurr_CD4 > Occurr_CD8, "CD4", ifelse(Occurr_CD4 < Occurr_CD8, "CD8","none"))) %>%
    mutate(Commit_accuracy = ifelse(Dominant_pheno == "CD4", Occurr_CD4/Occurrence*100,
                                    ifelse(Dominant_pheno == "CD8", Occurr_CD8/Occurrence*100,50))) %>%
    mutate(Dominant_pheno = factor(Dominant_pheno, levels = c("CD4", "no", "CD8")))
  ncol1 <- ncol(dataTable)
  # we want to have it in order that the aa_seq is first, then cd4 and cd8 samples and then the rest
  colOrder <- c(ncol1-2,cd4_samples,cd8_samples,ncol1-5,ncol1-4,ncol1-3,ncol1-1,ncol1)
  dataTable.result <- dataTable[,colOrder]
  return(dataTable.result)
}


generateSimulationsTable <- function(samples,dataTable) {
  for (i in 1:1000){
    simulation <- permatswap(samples,fixedmar = "both",mtype = "prab", times = 1)$perm[[1]]
    simulation.new <- getMeanAccurance(simulation,dataTable$aa_seq[1:(nrow(dataTable)-5)])
    meanCurrent <- mean(simulation.new$Commit_accuracy)
    simulation.new$Commit_accuracy <- round(simulation.new$Commit_accuracy, digits = 2)
  
    weighted_avg <- weighted.mean(simulation.new$Commit_accuracy, simulation.new$Occurrence, na.rm = TRUE)
    strictly_CD4 = length(which(simulation.new$Occurr_CD8 == 0))
    strictly_CD8 = length(which(simulation.new$Occurr_CD4 == 0))
    both = nrow(simulation.new) - strictly_CD4 - strictly_CD8 - 5
    simulationsTable.newRow <- data.frame(ID = i,Mean_Accurancy = meanCurrent,
                                          Weighted_Mean = weighted_avg,
                                          Strictly_CD4 = strictly_CD4,
                                          Strictly_CD8 = strictly_CD8,
                                          Undetermined = both)
    print(i)
    simulationsTable <- rbind(simulationsTable,simulationsTable.newRow)
    print("*****")
  }
  return(simulationsTable)
}
#data.frame(3,4)
simulationTable.3a <- generateSimulationsTable(samples.3a,countTable.3a)
write_csv(simulationsTable.3a, "250716_3a_simulations_1k_table_result.csv")

#### beta ####
simulationsTable <- data.frame(ID = "",Mean_Accurancy = "",Weighted_Mean = "", 
                               Strictly_CD4 = "",Strictly_CD8 = "", Undetermined = "")
simulationsTable <- simulationsTable[-1,]
samples.3b <- countTable.3b[1:268,2:32]
group <- countTable.3b[272,2:32]
cd4_samples <- which(group == "CD4")
cd8_samples <- which(group == "CD8")

simulationTable.3b <- generateSimulationsTable(samples.3b,countTable.3b)
write_csv(simulationTable.3b, "250716_3b_simulations_1k_table_result.csv")



ggplot(simulationTable.3b, aes(x = Mean_Accurancy)) +
 geom_histogram(bins = 50, fill = "skyblue", color = "black") +
   labs(
     title = "Mean values of commitment accurancy simulations",
     x = "Mean of commit accurancies",
     y = "Abs counts"
   ) +
  # our mean value of commitment accurancy
  geom_vline(xintercept = 85.90598) +
  geom_text(x=90.2, y=30, label="85.90598") +
   theme_classic()
ggsave("250716_3b_simulations_of_means.png")


ggplot(simulationTable.3b, aes(x = Weighted_Mean)) +
  geom_histogram(bins = 50, fill = "#097725", color = "black") +
  labs(
    title = "Weighted average values of commitment accurancy simulations",
    x = "Weighted mean of commit accurancies",
    y = "Abs counts"
  ) +
  geom_vline(xintercept = 82.47451) +
  geom_text(x=86.7, y=30, label="82.47451") +
  theme_classic()
ggsave("250716_3b_simulations_of_weighted_means.png")

#t.test(simulationsTable$Mean_Accurancy)



mean(simulationsTable$Mean_Accurancy >= 90.46551)
(sum(simulationsTable$Mean_Accurancy >= 90.46551) + 1) / (length(simulationsTable$Mean_Accurancy >= 90.46551) + 1)


(sum(simulationsTable$Weighted_Mean >= 86.81582) + 1) / (length(simulationsTable$Weighted_Mean >= 86.81582) + 1)






#runif samples from a uniform distribution
# means <- c()
# for (i in 1:1000){
#   xs <- runif(runs,min=0,max=1)
#   means[i] <- mean(xs)
# }
# ys <- 0:1
# df <- data.frame(mean = means)
# 
# # Plot histogram of the means
# ggplot(df, aes(x = mean)) +
#   geom_histogram(bins = 50, fill = "skyblue", color = "black") +
#   labs(
#     title = "Distribution of Means from Uniform(0,1)",
#     x = "Sample Mean",
#     y = "Frequency"
#   ) +
#   theme_classic()
