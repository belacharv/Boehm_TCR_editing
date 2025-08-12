### SIMULATIONS BOEHM'S PAPER ###
### The goal of this script is to generate simulations where the presence of the sequence is distributed in samples randomly
# the occurence of each sequence is either preserved from the real data or 
 
## Monte Carlo simulations for null hypothesis
countTable.3a <- read_csv(paste(getwd(),"results_table/MM.xp.3a_results","250702_commit_accuracy_non-unique_seq_spleen.csv",sep="/"))
countTable.3b <- read_csv(paste(getwd(),"results_table/MM.xp.3b_results","250714_3b_commit_accuracy_non-unique_seq_spleen.csv",sep="/"))

#countTable.3a <- read_csv("250808_3a_commit_accuracy_same_no_of_cd4_cd8.csv")


# for alpha 780 sequences
simulationsTable <- data.frame(ID = "",Mean_Accurancy = "",Weighted_Mean = "", 
                               Strictly_CD4 = "",Strictly_CD8 = "", Undetermined = "")
simulationsTable <- simulationsTable[-1,]

# extracting the table where are just samples, no metadata rows down the table or 
  # statistics in later columns
samples.3a <- countTable.3a[1:779,2:16] # cut off the metadata rows and cols
group <- countTable.3a[783,2:16] #783 is a row with the information about cell type (CD4/CD8)
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


generateSimulationsTable <- function(samples,dataTable,preserve_all,fixed_occurance) {
  for (i in 1:1000){
    if (fixed_occurance){
      total_ones <- nrow(samples) * 2
      simulation <- matrix(0, nrow = nrow(samples), ncol = ncol(samples))
      cols_chosen <- replicate(nrow(samples), sample(ncol(samples), 2))
      rows_rep <- rep(1:nrow(samples), each = 2)
      simulation[cbind(rows_rep, as.vector(cols_chosen))] <- 1
    }
    else{
      if (preserve_all){
        simulation <- permatswap(samples,fixedmar = "both",mtype = "prab", times = 1)$perm[[1]]
      }
      else{
        simulation <- permatfull(samples,fixedmar = "rows",mtype = "prab", times = 1)$perm[[1]] 
      }
    }
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
simulationTable.3a <- generateSimulationsTable(samples.3a,countTable.3a, preserve_all = TRUE,fixed_occurance = FALSE)
#write_csv(simulationsTable.3a, "250716_3a_simulations_1k_table_result.csv")
#write_csv(simulationTable.3a,"250808_3a_simulations_1k_filtered_for_same_ratio_cd4_cd8.csv")

simulationTable.3a.preservRows <- generateSimulationsTable(samples.3a,countTable.3a, preserve_all = FALSE, fixed_occurance = FALSE)
simulationTable.3a.fixedOccur <- generateSimulationsTable(samples.3a,countTable.3a, preserve_all = FALSE, fixed_occurance = TRUE)
write_csv(simulationTable.3a, "250811_3a_simulations_1k_table_preserved_just_rows_result.csv")

#### beta ####
simulationsTable <- data.frame(ID = "",Mean_Accurancy = "",Weighted_Mean = "", 
                               Strictly_CD4 = "",Strictly_CD8 = "", Undetermined = "")
simulationsTable <- simulationsTable[-1,]
samples.3b <- countTable.3b[1:268,2:32]
group <- countTable.3b[272,2:32]
cd4_samples <- which(group == "CD4")
cd8_samples <- which(group == "CD8")




simulationTable.3b <- generateSimulationsTable(samples.3b,countTable.3b, preserve_all = TRUE, fixed_occurance = FALSE)
#write_csv(simulationTable.3b, "250716_3b_simulations_1k_table_result.csv")

simulationTable.3b.preservRows <- generateSimulationsTable(samples.3b,countTable.3b, preserve_all = FALSE, fixed_occurance = FALSE)
simulationTable.3b.fixedOccur <- generateSimulationsTable(samples.3b,countTable.3b, preserve_all = FALSE, fixed_occurance = TRUE)


################################################################################
#### PLOTS ####
################################################################################




################################################################################
#### PRESERVED BOTH ROWS AND COLS ####
ggplot(simulationTable.3b, aes(x = Mean_Accurancy)) +
  geom_histogram(bins = 50, fill = "#FCD74C", color = "black") +
  labs(
    title = "Beta mean of commit accurancy preservAll",
    x = "Mean of commit accurancies",
    y = "Abs counts"
  ) +
  # our mean value of commitment accurancy
  geom_vline(xintercept = 85.90598) +
  geom_text(x=85.7, y=30, label="85.90598") +
  theme_classic()
#ggsave("250716_3b_simulations_of_means.png", width = 12, height = 9, units = "cm")
ggsave("250811_3b_simulations_of_means.png", width = 16, height = 8, units = "cm")


ggplot(simulationTable.3b, aes(x = Weighted_Mean)) +
  geom_histogram(bins = 50, fill = "#81E374", color = "black") +
  labs(
    title = "Beta weighted mean of commit accurancy preservAll",
    x = "Weighted mean of commit accurancies",
    y = "Abs counts"
  ) +
  geom_vline(xintercept = 82.47451) +
  geom_text(x=82.3, y=30, label="82.47451") +
  theme_classic()
#ggsave("250716_3b_simulations_of_weighted_means.png", width = 12, height = 9, units = "cm")
ggsave("250811_3b_simulations_of_weighted_means.png", width = 16, height = 8, units = "cm")


ggplot(simulationTable.3a, aes(x = Mean_Accurancy)) +
  geom_histogram(bins = 50, fill = "#74DCE3", color = "black") +
  labs(
    title = "Alpha mean of commit accurancy preservAll",
    x = "Mean of commit accurancies",
    y = "Abs counts"
  ) +
  # our mean value of commitment accurancy
  geom_vline(xintercept = 90.46551) +
  geom_text(x=90.1, y=30, label="90.46551") +
  theme_classic()
#ggsave("250808_3a_simulations_of_means_filtered.png")
ggsave("250811_3a_simulations_of_means.png", width = 16, height = 8, units = "cm")

ggplot(simulationTable.3a, aes(x = Weighted_Mean)) +
  geom_histogram(bins = 50, fill = "#E38374", color = "black") +
  labs(
    title = "Alpha Weighted mean val of commit accurancy",
    x = "Weighted mean of commit accurancies",
    y = "Abs counts"
  ) +
  geom_vline(xintercept = 86.81582) +
  geom_text(x=86.6, y=50, label="86.81582") +
  theme_classic()
ggsave("250811_3a_simulations_of_weighted_means.png", width = 16, height = 8, units = "cm")


################################################################################
#### PRESERVED JUST ROWS ####
ggplot(simulationTable.3b.preservRows, aes(x = Mean_Accurancy)) +
 geom_histogram(bins = 50, fill = "#C29A02", color = "black") +
   labs(
     title = "Beta mean of commit accurancy preservRows",
     x = "Mean of commit accurancies",
     y = "Abs counts"
   ) +
  # our mean value of commitment accurancy
  geom_vline(xintercept = 85.90598) +
  geom_text(x=84.2, y=30, label="85.90598") +
   theme_classic()
#ggsave("250716_3b_simulations_of_means.png", width = 12, height = 9, units = "cm")
ggsave("250811_3b_simulations_of_means_just_rows.png", width = 16, height = 8, units = "cm")


ggplot(simulationTable.3b.preservRows, aes(x = Weighted_Mean)) +
  geom_histogram(bins = 50, fill = "#097725", color = "black") +
  labs(
    title = "Beta weighted mean of commit accurancy preservRows",
    x = "Weighted mean of commit accurancies",
    y = "Abs counts"
  ) +
  geom_vline(xintercept = 82.47451) +
  geom_text(x=81, y=30, label="82.47451") +
  theme_classic()
#ggsave("250716_3b_simulations_of_weighted_means.png", width = 12, height = 9, units = "cm")
ggsave("250811_3b_simulations_of_weighted_means_just_rows.png", width = 16, height = 8, units = "cm")


ggplot(simulationTable.3a.preservRows, aes(x = Mean_Accurancy)) +
  geom_histogram(bins = 50, fill = "#2D418A", color = "black") +
  labs(
    title = "Alpha mean of commit accurancy preservRows",
    x = "Mean of commit accurancies",
    y = "Abs counts"
  ) +
  # our mean value of commitment accurancy
  geom_vline(xintercept = 90.46551) +
  geom_text(x=89, y=30, label="90.46551") +
  theme_classic()
#ggsave("250808_3a_simulations_of_means_filtered.png")
ggsave("250811_3a_simulations_of_means_only_rows.png", width = 16, height = 8, units = "cm")

ggplot(simulationTable.3a.preservRows, aes(x = Weighted_Mean)) +
  geom_histogram(bins = 50, fill = "#96220B", color = "black") +
  labs(
    title = "Alpha Weighted mean of commit accurancy preservRows",
    x = "Weighted mean of commit accurancies",
    y = "Abs counts"
  ) +
  geom_vline(xintercept = 86.81582) +
  geom_text(x=85.2, y=50, label="86.81582") +
  theme_classic()
#ggsave("250808_3a_simulations_of_weighted_means_filtered.png")
ggsave("250811_3a_simulations_of_weighted_means_only_rows.png", width = 16, height = 8, units = "cm")



################################################################################
### FIXED OCCURANCY FOR 2 SAMPLES ###
ggplot(simulationTable.3a.fixedOccur, aes(x = Mean_Accurancy)) +
  geom_histogram(bins = 50, fill = "#BE79F2", color = "black") +
  labs(
    title = "Alpha mean of commit accurancy fixedOccur",
    x = "Mean of commit accurancies",
    y = "Abs counts"
  ) +
  # our mean value of commitment accurancy
  geom_vline(xintercept = 90.46551, color = "black", lwd = 2) +
  geom_text(x=87, y=40, label="Real data observed mean")  +
  geom_text(x=89, y=30, label="90.46551") +
  theme_classic()
#ggsave("250808_3a_simulations_of_means_filtered.png")
ggsave("250811_3a_simulations_of_means_fixed_occur.png", width = 16, height = 8, units = "cm")

ggplot(simulationTable.3b.fixedOccur, aes(x = Mean_Accurancy)) +
  geom_histogram(bins = 50, fill = "#F2AF49", color = "black") +
  labs(
    title = "Beta mean of commit accurancy fixedOccur",
    x = "Mean of commit accurancies",
    y = "Abs counts"
  ) +
  # our mean value of commitment accurancy
  geom_vline(xintercept = 85.90598, color = "black", lwd = 2) +
  geom_text(x=83, y=38, label="Real data observed mean") +
  geom_text(x=84, y=30, label="85.90598") +
  theme_classic()
#ggsave("250716_3b_simulations_of_means.png", width = 12, height = 9, units = "cm")
ggsave("250811_3b_simulations_of_means_fixed_occur.png", width = 16, height = 8, units = "cm")

?geom_vline
#t.test(simulationsTable$Mean_Accurancy)

mean(simulationTable.3a$Mean_Accurancy)
mean(simulationTable.3a$Weighted_Mean)

mean(simulationTable.3b$Mean_Accurancy)
mean(simulationTable.3b$Weighted_Mean)

mean(simulationTable.3a.preservRows$Mean_Accurancy)
mean(simulationTable.3a.preservRows$Weighted_Mean)

mean(simulationTable.3b.preservRows$Mean_Accurancy)
mean(simulationTable.3b.preservRows$Weighted_Mean)


mean(simulationTable.3a.fixedOccur$Mean_Accurancy)
mean(simulationTable.3b.fixedOccur$Mean_Accurancy)


mean(simulationsTable$Mean_Accurancy >= 90.46551)
(sum(simulationsTable$Mean_Accurancy >= 90.46551) + 1) / (length(simulationsTable$Mean_Accurancy >= 90.46551) + 1)


(sum(simulationsTable$Weighted_Mean >= 86.81582) + 1) / (length(simulationsTable$Weighted_Mean >= 86.81582) + 1)



simulation <- permatswap(samples.3a,fixedmar = "rows",mtype = "prab", times = 1)$perm[[1]]

simulation <- permatfull(samples.3a,fixedmar = "rows",mtype = "prab", times = 1)$perm[[1]]
















# set.seed(123)
# n_rows <- 779
# n_cols <- 15

# Total number of ones needed
total_ones <- n_rows * 2

# Preallocate matrix
mat <- matrix(0, nrow = n_rows, ncol = n_cols)

# Randomly choose positions for each row
cols_chosen <- replicate(n_rows, sample(n_cols, 2))

# Convert to long vector of row indices and col indices
rows_rep <- rep(1:n_rows, each = 2)

# Assign ones directly
mat[cbind(rows_rep, as.vector(cols_chosen))] <- 1











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
