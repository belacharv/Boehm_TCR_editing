################################################################################
### COUNT TABLE #####
################################################################################

getwd()
library("xlsx")

#setwd(paste(getwd(),"results_table",sep=""))
## Load data

# if we wanna run it on alpha
finalResult <- read_csv(paste(getwd(),"results_table/MM.xp.3a_results","250428_result_filtered_table_3a.csv",sep="/"))
getwd()
# in case we want to deal with beta
finalResult <- read_csv(paste(getwd(),"results_table/MM.xp.3b_results","250714_result_filtered_table_3b.csv",sep="/"))


# prepare count table
count_table <- data.frame(finalResult$CDR3.amino.acid.sequence)
samples_names <- unique(finalResult$Sample.bio.name)
samples_names

library(tibble)
## creating a table with metadata info that can be added to the end of the 
  ## count table
df_unique <- finalResult %>%
  distinct(Sample.bio.name, .keep_all = TRUE)
df <- data.frame(ID = c(1:length(samples_names)),df_unique[,c(29,26,25,27,15)])
dfT <- t(df)
dfT <- data.frame(dfT)
dfT <- rownames_to_column(dfT,"aa_seq")
colnames(dfT)[2:ncol(dfT)] <- samples_names


## MM_438_S and
excludeColList <- c()
for (col in 1:ncol(dfT)){
  # we want to keep only spleen samples and are interested in populations CD4 and CD8
  if (dfT[[3,col]] == "T" | is.na(dfT[[4,col]]) | dfT[[4,col]] == "DN" | dfT[[4,col]] == "DP"){
    excludeColList[length(excludeColList) +1] = col
  }
}


#rownames(count_table) <- count_table$finalResult.CDR3.amino.acid.sequence
sorted_aa <- data.frame(aa_seq = sort(count_table$finalResult.CDR3.amino.acid.sequence))
count_Table1 <- unique(sorted_aa)
for (name in samples_names) {
  # create the template for a count table with only zeros for now
  count_Table1 <- mutate(count_Table1, !!name := 0)
}

## testing
#sample_name = samples_names[1]
#sample_table <- filter(finalResult, finalResult$Sample.bio.name == sample_name)



#Seq_count_sample <- rbind(Seq_count_sample,c("he",0))
# function filling the count table with umi
umiCountSample <- function(sample_name, countTable, umiValue) {
  sample_table <- filter(finalResult, finalResult$Sample.bio.name == sample_name)
  #print(paste(sample_name,nrow(sample_table)))
  #print(nrow(sample_table))
  #Seq_count_sample <- rbind(Seq_count_sample, c(sample_name,nrow(sample_table)))
  
  for (aaseq in 1:nrow(sample_table)) {
    nameSeq <- sample_table$CDR3.amino.acid.sequence[aaseq]
    value <- sample_table$Umi.count[aaseq]
    if (umiValue == TRUE) {
      # if we want the exact number of umi of that seq in that sample
      countTable[countTable$aa_seq == nameSeq,sample_name] <- value
      }
    else{
      # if we just wanna know whether the seq is in the sample (1) or not (0)
      countTable[countTable$aa_seq == nameSeq,sample_name] <- 1}
  }
  return(countTable)
}


#count_table1 = count_Table1
count_Table2 = data.frame(count_Table1)
for (sample in samples_names){
  count_Table1 = umiCountSample(sample,count_Table1,TRUE)
  count_Table2 = umiCountSample(sample,count_Table2,FALSE)
  #print(paste("done",sample))
}


#write_csv(count_Table1, "250505_count_table_seq_samples.csv")
#write.xlsx(count_Table1, "250505_count_table_seq_samples.xlsx")

#write_csv(count_Table2, "250627_count_table2_seq_samples.csv")

#write.xlsx(count_Table2, "250627_count_table_seq_samples.xlsx")


samples_names


### merging the metadata and the count table ###
countTable <- count_Table1
rownames(countTable) <- count_Table1$aa_seq
#countTable <- countTable[-1]
countTableFinal <- rbind(countTable,dfT)
#colnames(dfT)
#colnames(countTable)

#write_csv(countTableFinal, "250626_count_table_with_metadata.csv")
#write.xlsx(countTableFinal, "250626_count_table_with_metadata.xlsx")


################################################################################
### FILTERING OUT ALL THE SEQ IN JUST 1 SAMPLE AND UNPRODUCTIVE SEQ ###
################################################################################
# we want to get rid of the thymus samples
countTableSpleen <- count_Table2[,-excludeColList]

Seq_count_sample <- data.frame(df)
Seq_count_sample <- mutate(Seq_count_sample,seq_count = 0)

excluderows <- c()
for(i in 1:nrow(Seq_count_sample)){
  print(i)
  if(Seq_count_sample[[i,3]] == "T" || 
     is.na(Seq_count_sample[[i,4]])  ||
     Seq_count_sample[[i,4]] == "DN" ||
     Seq_count_sample[[i,4]] == "DP"){
    excluderows[length(excluderows)+1] <- i
    
  }
}
Seq_count_sample <- Seq_count_sample[-excluderows,]

for (sample in Seq_count_sample$Sample.bio.name) {
  sample_table <- filter(finalResult, finalResult$Sample.bio.name == sample)
  print(nrow(sample_table))
  Seq_count_sample$seq_count[Seq_count_sample$Sample.bio.name == sample] <- nrow(sample_table)
}
#Seq_count_sample<- Seq_count_sample[-1,]
#write_csv(Seq_count_sample, "250808_3b_seq_number_spleen_samples.csv")

samples_to_exclude <- Seq_count_sample %>%
 filter(cell == "CD4") %>%
  slice_min(order_by = seq_count,n=5)
  #which(low_seq_count_filtr$tissue == "CD4")
countTableSpleen.filtr <- countTableSpleen[,!(colnames(countTableSpleen) %in% samples_to_exclude$Sample.bio.name)]
countTableSpleen <- countTableSpleen.filtr

countTableSpleen <- countTableSpleen %>% 
  mutate(Occurrence = rowSums(countTableSpleen[,2:ncol(countTableSpleen)])) %>%
  filter(Occurrence > 1) %>% # only seqs that are Occurring in multiple samples
  filter(!grepl("\\*", aa_seq)) # get rid of the unproductive seqs

dfT <- dfT %>% mutate(Occurrence = NA)
dfTSpleen <- dfT[,-excludeColList]
dfTSpleen <- dfTSpleen[,!(colnames(dfTSpleen) %in% samples_to_exclude$Sample.bio.name)]
colnames(dfTSpleen)
colnames(countTableSpleen)
countTableSpleen <- rbind(countTableSpleen,dfTSpleen) # add metadata to the table
No_of_samples <- ncol(dfTSpleen)-2

group <- as.character(dfTSpleen[4,])  
cd4_samples <- which(group == "CD4")
cd8_samples <- which(group == "CD8")

#colnames(dfT)
colnames(countTableSpleen)

countTableSpleen <- countTableSpleen %>% 
  mutate(Occurr_CD4 = 0) %>%
  mutate(Occurr_CD8 = 0)


#No_of_samples <- length(samples_names) - length(excludeColList)
sample_cols <- 2:(No_of_samples+1)
sample_data <- cbind(0,countTableSpleen[, sample_cols])
sample_data <- as.data.frame(lapply(sample_data, as.numeric))

countTableSpleen$Occurr_CD4 <- rowSums(sample_data[, cd4_samples] > 0)
countTableSpleen$Occurr_CD8 <- rowSums(sample_data[, cd8_samples] > 0)


cols <- colnames(countTableSpleen)
non_sample_cols <- c("aa_seq", "Occurrence", "Occurr_CD4", "Occurr_CD8")
sample_cols <- setdiff(cols, non_sample_cols)
sample_groups <- group[1:length(sample_cols)]  
ordered_samples <- sample_cols[order(factor(sample_groups, levels = c("CD4", "CD8", NA)))]
final_col_order <- c("aa_seq", ordered_samples, non_sample_cols[non_sample_cols != "aa_seq"])
countTableSpleen <- countTableSpleen[, final_col_order]


  
countTableSpleen <- countTableSpleen %>%
  mutate(Dominant_pheno = ifelse(Occurr_CD4 > Occurr_CD8, "CD4", ifelse(Occurr_CD4 < Occurr_CD8, "CD8","none"))) %>%
  mutate(Commit_accuracy = ifelse(Dominant_pheno == "CD4", Occurr_CD4/Occurrence*100,
                                  ifelse(Dominant_pheno == "CD8", Occurr_CD8/Occurrence*100,50))) %>%
  mutate(Dominant_pheno = factor(Dominant_pheno, levels = c("CD4", "no", "CD8"))) #%>%
  #arrange(Dominant_pheno)
#countTableSpleen$Commit_accuracy <- round(countTableSpleen$Commit_accuracy, digits = 2)
  
mean(countTableSpleen$Commit_accuracy[1:(nrow(countTableSpleen)-6)])

countTableSpleen$Commit_accuracy <- round(countTableSpleen$Commit_accuracy, digits = 2)

weighted_avg <- weighted.mean(countTableSpleen$Commit_accuracy[1:(nrow(countTableSpleen)-6)], countTableSpleen$Occurrence[1:(nrow(countTableSpleen)-6)], na.rm = TRUE)
weighted_avg
seqCommitCounts <- data.frame("strictly_CD4" = length(which(countTableSpleen$Occurr_CD8 == 0)), 
                              "strictly_CD8" = length(which(countTableSpleen$Occurr_CD4 == 0)), 
                              "both" = nrow(countTableSpleen) - length(which(countTableSpleen$Occurr_CD8 == 0)) - length(which(countTableSpleen$Occurr_CD4 == 0))-6)


#write_csv(countTableSpleen, "250714_3b_commit_accuracy_non-unique_seq_spleen.csv")
#write_csv(countTableSpleen,"250808_3a_commit_accuracy_same_no_of_cd4_cd8.csv")

################################################################################
### WILL NOT BE USED ###
################################################################################
# ## Only repeating sequences
# only_aa_seq <- data.frame(aa_seq="",samples = "", Same_mouse = "")
# samples_having_this_aa <- hash() 
# for (i in 2:nrow(count_Table2)){
#   if (count_Table2$Occurrence[i] > 1) {
#     key = count_Table2$aa_seq[i]
#     values = list()
#     names = ""
#     mice = c()
#     #extract_sample_name = hash()
#     for (column in 2:(ncol(count_Table2)-1)){
#       if (count_Table2[i,column] >0){
#         values = append(values,colnames(count_Table2)[column])
#         names = paste(names,colnames(count_Table2)[column],", ")
#         if (substr(colnames(count_Table2)[column],1,8) %in% mice){
#           print("*******************************************")
#         }
#         else {
#           mice = c(mice,substr(colnames(count_Table2)[column],1,8))
#         }
#       }
#     }
#     print(mice)
#     #same_mouse = 
#     if (length(mice) == 1) {same_mouse = TRUE}
#     else {same_mouse = FALSE}
#     samples_having_this_aa[[key]] = values
#     #print(length(values))
#     #new_row = data.frame()
#     only_aa_seq = rbind(only_aa_seq,c(count_Table2$aa_seq[i],names,same_mouse))
#   }
# }
# result = only_aa_seq[-1,]
# ## will not be used later because it's not numerical
# #write_csv(result, "250518_repeating_aa_seq.csv")
# #write.xlsx(result, "250518_repeating_aa_seq.xlsx")
# 
# 


# 
# ### is this section even needed? ###
# x <- count_Table1[1,]
# shared_aa_seq <- data.frame(aa_seq=c(""),Amount_of_samples=c(""),Samples=c(""))
# x1 <- x[x!=0]
# non_zero <- which(x!=0)
# samples_with_seq = c()
# for (column_index in 2:length(non_zero)){
#   sample_with_this_aa <- colnames(count_Table1[column_index])
#   print(sample_with_this_aa)
#   #samples_with_seq.append(sample_with_this_aa)
#   samples_with_seq = append(samples_with_seq,sample_with_this_aa)
# }
# 
# 
# for (row in 1:nrow(count_Table1)){
#   x <- count_Table1[row,]
#   non_zero <- which(x!=0)
#   if (length(non_zero) >2){
#     #print(count_Table1$aa_seq[row])
#     samples_with_seq = c()
#     for (column_index in 2:length(non_zero)){
#       sample_with_this_aa <- colnames(count_Table1[column_index])
#       print(sample_with_this_aa)
#       #samples_with_seq.append(sample_with_this_aa)
#       samples_with_seq = append(samples_with_seq,sample_with_this_aa)
#     }
#     print(samples_with_seq)
#     shared_aa_seq = rbind(shared_aa_seq,c(count_Table1$aa_seq[row],length(non_zero)-1,samples_with_seq))
#   }
# }
# shared_aa_seq <- shared_aa_seq[-1,]
# shared_aa_seq <- mutate(shared_aa_seq, Amount_of_samples = as.numeric(shared_aa_seq$Amount_of_samples))


#rbind(shared_aa_seq,c(count_Table1$aa_seq[row],length(non_zero)-1,FALSE))
#getwd()


#count_table <- count_table %>%
#length(samples_names) 

