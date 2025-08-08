################################################################################
### FILTERING AA SEQ FROM BOEHMS PAPER DATA ###
################################################################################

## This R script works with data from Boehms publication Origin and evolutionary 
## malleability of T cell receptor α diversity

### The goal is to have table with 
#library(data.table)
library(NLP)
library(tidyverse)
library(hash)


## load data
load("MM.xp.3a") # Counts tables with the sequences from sort experiments (TCR alpha edited mice)

load("MM.xp.3b") # Counts tables with the sequences from sort experiments (TCR beta edited mice)


################################################################################
## AA SEQ ##

# restoring the original conditions for filtrating the AA seq

# 
# resultTable <- data.table::copy(MM.xp.3a)
# resultTable <- mutate(resultTable, Used = FALSE)
# resultTabulaAA <- filter(resultTable,CDR3.nucleotide.sequence == 'a')



################################################################################
# Fuction that returns list of doubles with indexes of the rows that we want to 
# merge together


# This function returns the list of columns we want to check if they are the same
# in row1 and row2 during merging
findIndexOfCols<- function(data) {
  # selected columns that should be the same
  neededMatch <- c("genotype","cell","tissue","tetramer","raton","Sample.bio.name")
  indexOfCols <- c(0,0,0,0,0,0)
  for (name in 1:6){
    for (column in 1:ncol(data)) { 
      if (colnames(data)[column] == neededMatch[name]) { 
        indexOfCols[name] <- column
      }
    }
  }
  # returns indexes of the columns that need to be the same
  return(indexOfCols)
}
#indexOfCols = findIndexOfCols()print(findIndexOfCols(resultTabulaAA) )


# function where we get the filtrated table with seq_name and all the reads present
# output is a table of two columns where we have the indexes of the 
  # pairs of rows we need to merge
getNeededRows <- function(randomData,seq_name) {
  sampleBioNamesHash = hash()
  row = 1
  needToBeMerged = data.frame(row1 = c(""),row2 = c(""))
  while (row <= nrow(randomData)) {
    print(row)
    sample_name = randomData$Sample.bio.name[row]
    # if there are two rows with the same Sample.bio.name ( = from same sample)
    if(!is.null(sampleBioNamesHash[[sample_name]])) {
      row2 = sampleBioNamesHash[[sample_name]]
      print(paste(row, row2, sep = " "))
      # add them to the table of what needs to be merged
      needToBeMerged = rbind(needToBeMerged,c(row,row2))
      } 
    else { # if this the first time we found the seq in the sample
      sampleBioNamesHash[[sample_name]] = row }
    row = row +1
    }
  print(sampleBioNamesHash)
  # getting rid of the first empty row
  needToBeMerged = needToBeMerged[2:nrow(needToBeMerged),]
  print(needToBeMerged)
  return(needToBeMerged)
}




mergeRowsInTable <- function(data, needToBeMerged, matchCols) {
  removedRows = c()
  for (i in 1:nrow(needToBeMerged)) {
    row1 <- as.numeric(needToBeMerged$row1[i])
    row2 <- as.numeric(needToBeMerged$row2[i])
    print(i)
    print(paste(row1,row2,sep=" "))
    
    # Adjust indices if previous rows were deleted
    row1_adj <- row1 - sum(removedRows < row1)
    row2_adj <- row2 - sum(removedRows < row2)
    
    # Check again after adjustment
    if (row1_adj <= 0 || row2_adj <= 0 || 
        row1_adj > nrow(data) || row2_adj > nrow(data)) next
    
    # Merge the rows
    data <- MergeSequences(row1_adj, row2_adj, data, matchCols)
    
    # Track the removed row (row2 was removed)
    removedRows <- c(removedRows, row2)
  }
  
  return(data)
}


### Compare two rows and if they're the same, call merging function ############

# matchCols is output of indexOfCols
MergeSequences <- function(row1,row2,data,matchCols){
  # check if all the parameters are the same so we can merge the two rows together
  for (name in 1:5) 
  {
    print(data[[row1,matchCols[name]]])
    print(data[[row2,matchCols[name]]])
    # if some of the values is NA
    if((is.na(data[[row1,matchCols[name]]])) || (is.na(data[[row2,matchCols[name]]]))) 
    { 
      print(paste("missing values in",name, sep=" "))
      return(data)
    }
    if(data[[row1,matchCols[name]]] == data[[row2,matchCols[name]]])
    { 
      print("this row ok")
    }
    else {
      # it crashed on one of the names we want to check
      print(paste("not a match in", name, sep = " "))
      return(data)}
  }
  # if everything was okay and we went through all the parameters, we can
  # say that the rows fit and merge them together
  data <- mergeTwoRows(row1,row2,data)
  print(paste("now we merged",row1,row2, sep = " "))
  print(paste("number of rows",nrow(data), sep = " "))
  return(data)
}

### Function to merge the two rows together ####################################

mergeTwoRows <- function(row1,row2,data) {
  print("hey there")
  ## The parameters we want to change because of merging the rows
  data$Umi.count[row1] = data$Umi.count[row1] + data$Umi.count[row2]
  data$Read.count[row1] = data$Read.count[row1] + data$Read.count[row2]
  data$Umi.proportion[row1] = data$Umi.proportion[row1] + data$Umi.proportion[row2]
  data$Read.proportion[row1] = data$Read.proportion[row1] + data$Read.proportion[row2]
  # delete the second row
  dataNew <- data[-row2,]
  print("how's it goin?")
  return(dataNew)
}




#resultTab <- sortOnlySameAA(randomData, indexOfCols)
#randomTry <- sortOnlySameAA(randomData, indexOfCols)

################################################################################
## GLOBAL FUNCTION ##

#resultTable <- data.table::copy(MM.xp.3a)
#resultTable <- mutate(resultTable, Used = FALSE)


# function going through all the rows in the table 
# calling all the amino acid sequences found
globalMerge <- function(resultTable, indexOfCols) {
  count <- 0
  resultTabulaAA <- filter(resultTable,CDR3.nucleotide.sequence == 'a')
  for (i in 1:nrow(resultTable)) {
    print(paste(i," row", sep = ""))
    if (resultTable$Used[i] == TRUE) { 
      print('ops') 
      next }
    seq_name <- resultTable$CDR3.amino.acid.sequence[i]
    print(seq_name)
    # mark all the reads with same seq as used
    resultTable$Used[resultTable$CDR3.amino.acid.sequence == seq_name] <- TRUE
    # filter out all the rows with this sequence
    new_table <- filter(resultTable,resultTable$CDR3.amino.acid.sequence == seq_name)
    # if there are more than 1 reads or more UMIs in the only sample
    needToBeMerged = getNeededRows(new_table,seq_name)
    print(needToBeMerged[[1,1]])
    if (is.na(needToBeMerged[[1,1]])) { newTableResult = new_table}
    else {newTableResult <- mergeRowsInTable(new_table,needToBeMerged, indexOfCols)}
    print(nrow(newTableResult))
    # adding to the final table
    # filtering out the aa seq with only one read and umi count 1
    if ((nrow(newTableResult) > 1) | (newTableResult$Umi.count[1] > 1)) { 
      count = count + nrow(newTableResult)
      print('hi')
      # add all the rows to our blank table
      for (row in 1:nrow(newTableResult)) {
        resultTabulaAA[nrow(resultTabulaAA)+1,] <- newTableResult[row,]
      }
    }
  }
  print(count)
  return(resultTabulaAA)
}

# mayberesult = globalMerge(resultTable, indexOfCols)
# write_csv(mayberesult, file = "250426_result_table_filtrated_3a.csv")
# write.xlsx(mayberesult,file = "250426_3a_filtered_table_boehm.xlsx")
# getwd()







# function with input as table for one sample filtreted merged
getUmiCountAndProportion <- function(data) {
  print("Before the proportion edit")
  umiCount <- sum(data$Umi.count)
  readCount <- sum(data$Read.count)
  ## check the sum of umi proportion and read proportion
  print(sum(data$Umi.proportion))
  print(sum(data$Read.proportion))
  for (row in 1:nrow(data)){
    umiRow <- data$Umi.count[row]
    readRow <- data$Read.count[row]
    data$Umi.proportion[row] <- umiRow/umiCount
    data$Read.proportion[row] <- readRow/readCount
  }
  print("After the proportion edit")
  print(sum(data$Umi.proportion))
  print(sum(data$Read.proportion))
  return(data)
}


#checkedTable <- filter(mayberesult,mayberesult$Sample.bio.name == "a" )
checkUmi <- function(dataTable) {
  # empty table with the info about columns
  checkedTable <- filter(dataTable,dataTable$Sample.bio.name == "a" )
  samples <- unique(dataTable$Sample.bio.name)
  for(sample in samples) {
    oneSampleTable <- filter(dataTable, dataTable$Sample.bio.name == sample)
    oneSampleResult <- getUmiCountAndProportion(oneSampleTable)
    checkedTable <- rbind(checkedTable, oneSampleResult)
  }
  return(checkedTable)
}



################################################################################
## FILTRATION RESULT ##
main <- function(){
  # load data
  load("MM.xp.3a") 
  load("MM.xp.3b")
  
  # prepare the data
  resultTable <- data.table::copy(MM.xp.3b)
  resultTable <- mutate(resultTable, Used = FALSE)
  resultTabulaAA <- filter(resultTable,CDR3.nucleotide.sequence == 'a')
  resultTable <- resultTable[resultTable$raton != "MM_327", ]
  
  # The columns we want to keep the same when merging two rows
  indexOfCols = findIndexOfCols(resultTable)
  
  # filtering out rows with 1 umi in 1 sample
  # merging the rows different only in nucleo seq, but with same aa seq
  tableAfterFiltration = globalMerge(resultTable, indexOfCols)
  # editing the final umi and read proportion
  checkedTable <- filter(tableAfterFiltration,tableAfterFiltration$Sample.bio.name == "a" )
  finalResult = checkUmi(tableAfterFiltration)
  #write_csv(finalResult,"250428_result_filtered_table_3a.csv")
  #write.xlsx(finalResult, "250428_result_filtered_table_3a.xlsx")
  return(finalResult)
}


#samples_names
#unique(finalResult$raton)
finalResult <- main()

og_names <- unique(MM.xp.3b$Sample.bio.name)
og_names


#write_csv(finalResult,"250714_result_filtered_table_3b.csv")







################################################################################
### COUNT TABLE #####
# 
# count_table <- data.frame(finalResult$CDR3.amino.acid.sequence)
# samples_names <- unique(finalResult$Sample.bio.name)
# samples_names
# #rownames(count_table) <- count_table$finalResult.CDR3.amino.acid.sequence
# sorted_aa <- data.frame(aa_seq = sort(count_table$finalResult.CDR3.amino.acid.sequence))
# count_Table1 <- unique(sorted_aa)
# for (name in samples_names) {
#   count_Table1 <- mutate(count_Table1, !!name := "0")
# }
# 
# 
# umiCountSample <- function(sample_name) {
#   sample_table <- filter(finalResult, finalResult$Sample.bio.name == sample_name)
#   print(nrow(sample_table))
#   #sample_table$CDR3.amino.acid.sequence
#   #aaseq_list <- unique(sample_table$CDR3.amino.acid.sequence)
#   for (aaseq in 1:nrow(sample_table)) {
#     nameSeq <- sample_table$CDR3.amino.acid.sequence[aaseq]
#     #print(nameSeq)
#     value <- sample_table$Umi.count[aaseq]
#     #print(value)
#     count_Table1[count_Table1$aa_seq == nameSeq,sample_name] <- 1#[count_table$finalResult.CDR3.amino.acid.sequence == !!sample_table$CDR3.amino.acid.sequence[aaseq]] = 
#   }
#   return(count_Table1)
# }
# 
# 
# for (sample in samples_names){
#   count_Table1 = umiCountSample(sample)
#   print(paste("done",sample))
# }
# 
# 
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
# 
# getwd()
# #rbind(shared_aa_seq,c(count_Table1$aa_seq[row],length(non_zero)-1,FALSE))
# #getwd()
# write_csv(count_Table1, "250505_count_table_seq_samples.csv")
# write.xlsx(count_Table1, "250505_count_table_seq_samples.xlsx")
# #count_table <- count_table %>%
# #length(samples_names) 
# samples_names

#new_table <- data.frame(MM.xp.3a$CDR3.amino.acid.sequence,MM.xp.3a$Umi.count)


# indexOfCols = findIndexOfCols(resultTable)
# #resultTabbi = sortOnlySameAA(resultTable,indexOfCols)
# 
# maybeResult = globalMerge(resultTable, indexOfCols)
# write.csv(maybeResult, file = "250407_3a_filtered_table_boehm.csv")
# write.xlsx(maybeResult,file = "250408_3a_filtered_table_boehm.xlsx")
# library("data.frame")
# library("xlsx")
# result <- checkUmi(mayberesult)
# write_csv(result,"250428_filtered_table_3a.csv")
# write.xlsx(result, "250428_filtered_table_3a.xlsx")
























#### DEBUG ####

## How many different genotypes there are?
findOutgeno = hash()
for (i in 1:nrow(resultTable)){
  if(is.null(findOutgeno[[resultTable$genotype[i]]])) {
    findOutgeno[[resultTable$genotype[i]]] = 1
    print(resultTable$genotype[i])
  } 
}

## How many samples have NA value in cell?
nullInCell = hash()
for (i in 1:nrow(resultTable)){
  if(is.na(resultTable$cell[i])) {
    if(is.null(nullInCell[[resultTable$raton[i]]])) {
      nullInCell[[resultTable$raton[i]]] = 1
      print(resultTable$raton[i])
    } 
  }
  # else {
  # sampleBioNamesHash[[sample_name]] = row }
}





## TESTING ON ONLY SELECTED AMINO ACIDS ##
indexOfCols <- findIndexOfCols(resultTable)

testMergingTable <- function(seqName, indexOfCols) {
  testData <- filter(resultTable, resultTable$CDR3.amino.acid.sequence == seqName)
  print(nrow(testData))
  needToBeMerged = getNeededRows(testData,seqName)
  testResult <- mergeRowsInTable(testData,needToBeMerged, indexOfCols)
  return(list(original = testData, merged = testResult, rowsSupposedToBeMerged = needToBeMerged))
}

firstSeqName = "CAANQGGRALIHSF"
firstTest = testMergingTable(firstSeqName, indexOfCols)
firstTestData = firstTest[[1]] # original
firstTestResult = firstTest[[2]] # merged
firstRowsNeededToBeMerged = firstTest[[3]] # rows supposed to be merged

secondSeqName = "CAANRRF"
secondTest = testMergingTable(secondSeqName, indexOfCols)
secondTestData = secondTest[[1]]
secondTestResult = secondTest[[2]]
secondRowsNeededToBeMerged = secondTest[[3]]

thirdSeqName = "CAANQGGRALN"
testData <- filter(resultTable, resultTable$CDR3.amino.acid.sequence == thirdSeqName)
thirdTest = testMergingTable(thirdSeqName, indexOfCols)
thirdTestData = thirdTest[[1]]
thirdTestResult = thirdTest[[2]]
thirdRowsNeededToBeMerged = thirdTest[[3]]




### TESTING UMI PROPORTION
unique(resultTable$Sample.bio.name)
print(sum(resultTable$Umi.proportion))
Sample.MM_3688a <- filter(mayberesult,mayberesult$Sample.bio.name == "MM_5826a_P5_T")
SampleTest <- getUmiCountAndProportion(Sample.MM_3688a)



################################################################################
################################################################################
### OLD ###
################################################################################
################################################################################


# 
# ## NUCLEO SEQ ##
# 
# resultTable <- data.table::copy(MM.xp.3a)
# resultTable <- mutate(resultTable, Used = FALSE)
# # getting blank table for adding newly obtained rows
# resultTabula <- filter(resultTable,CDR3.nucleotide.sequence == 'a')
# 
# 
# # filtering out the 3a that have either more then 1 umi in one sample, or more 
# count <- 0
# for (i in 1:nrow(MM.xp.3a)) {
#   print(i)
#   if (resultTable$Used[i] == TRUE) { 
#     print('ops') 
#     next }
#   seq_name <- MM.xp.3a$CDR3.nucleotide.sequence[i]
#   # mark all the reads with same seq as used
#   resultTable$Used[resultTable$CDR3.nucleotide.sequence == seq_name] <- TRUE
#   # filter out all the rows with this sequence
#   new_table <- filter(resultTable,resultTable$CDR3.nucleotide.sequence == seq_name)
#   # if there are more than 1 reads or more UMIs in the only sample
#   if ((nrow(new_table) > 1) | (MM.xp.3a$Umi.count[i] > 1)) { 
#     count = count + nrow(new_table)
#     print('hi')
#     # add all the rows to our blank table
#     for (row in 1:nrow(new_table)) {
#       resultTabula[nrow(resultTabula)+1,] <- new_table[row,]
#     }
#   }
# }
# 
# 
# # debug formulas
# seq_name <- MM.xp.3a$CDR3.amino.acid.sequence[8]
# resultTable$Used[resultTable$CDR3.amino.acid.sequence == seq_name] <- TRUE
# new_table <- filter(MM.xp.3a,MM.xp.3a$CDR3.amino.acid.sequence == seq_name)
# if ((nrow(new_table) > 1) | (MM.xp.3a$Umi.count[i] > 1)) { 
#   count = count + nrow(new_table)
#   print('hi')
#   for (row in 1:nrow(new_table)) {
#     resultTabulaAA[nrow(resultTabulaAA)+1,] <- new_table[row,]
#   }
# }
# 
# 



# for (i in 1:nrow(MM.xp.3a)) {
#   print(i)
#   if (resultTable$Used[i] == TRUE) { 
#     print('ops') 
#     next }
#   seq_name <- MM.xp.3a$CDR3.amino.acid.sequence[i]
#   resultTable$Used[resultTable$CDR3.amino.acid.sequence == seq_name] <- TRUE
#   new_table <- filter(resultTable,resultTable$CDR3.amino.acid.sequence == seq_name)
#   
#   if ((nrow(new_table) > 1) | (MM.xp.3a$Umi.count[i] > 1)) { 
#     print('hi')
#     for (row in 1:nrow(new_table)) {
#       resultTabulaAA[nrow(resultTabulaAA)+1,] <- new_table[row,]
#     }
#   }
# }

#write.csv(resultTabulaAA, file = "250210_result_aa_seq_filter.csv")

# ## probably wont be used
# sortOnlySameAA <- function(resultTable,matchCols) {
#   row1 = 1
#   SAMPLEBIONAME <- matchCols[6]
#   merged = c()
#   while (row1 < nrow(resultTable)){
#     print(row1)
#     if (resultTable$Used[row1] == TRUE) { 
#       print('ops')
#       row1 = row1 + 1
#       next }
#     seq_name <- resultTable$CDR3.amino.acid.sequence[row1]
#     resultTable$Used[resultTable$CDR3.amino.acid.sequence == seq_name] <- TRUE
#     new_table <- filter(resultTable,resultTable$CDR3.amino.acid.sequence == seq_name)
#     sampleBioNamesHash <- hash()
#     if ((nrow(new_table) > 1) | (MM.xp.3a$Umi.count[row1] > 1)) { 
#       print('hi')
#       for (row in 1:nrow(new_table)) {
#         if(!is.null(sampleBioNamesHash[[seq_name]])) {
#           row2 = sampleBioNamesHash[[seq_name]]
#           new_table = MergeSequences(row,row2,data = resultTable,matchCols)
#         }
#         sampleBioNamesHash[[seq_name]] = 1
#       }
#       for (row in 1:nrow(new_table)) {
#         resultTabulaAA[nrow(resultTabulaAA)+1,] <- new_table[row,]
#       }
#     }
#     row1 = row1 + 1
#   }
#   print("******************************************************************")
#   #print(length(merged))
#   print("******************************************************************")
#   return(resultTabulaAA)
# }

### function comparing each two rows with each other  ###

# editNewTable<- function(dataTable,matchCols) {
#   SAMPLEBIONAME <- matchCols[6] 
#   print(SAMPLEBIONAME)
#   view(dataTable)
#   print(nrow(dataTable))
#   row1 <- 1
#   row2 <- 2
#   while (row1<nrow(dataTable)) {
#     print("row1")
#     print(row1)
#     while (row2 <=(nrow(dataTable))) {
#       print("row2")
#       print(row2)
#       ## if it's the same sample
#       if(dataTable[[row1,SAMPLEBIONAME]] == dataTable[[row2,SAMPLEBIONAME]]) {
#         dataTable <- MergeSequences(row1,row2,dataTable,matchCols)
#         print("*************************************************")
#       }
#       row2 = row2 + 1
#     }
#     row1 = row1 + 1
#     row2 = row1 + 2
#   }
#   return(dataTable)
# }
# 




### New version ####

## function to find the index of columns we want to keep the same ###
# findIndexOfCols<- function(data) {
#   neededMatch <- c("genotype","cell","tissue","tetramer","raton","Sample.bio.name")
#   indexOfCols <- c(0,0,0,0,0,0)
#   for (name in 1:6){
#     for (column in 1:ncol(data)) { 
#       if (colnames(data)[column] == neededMatch[name]) { 
#         indexOfCols[name] <- column
#       }
#     }
#   }
#   return(indexOfCols)
# }
# print(findIndexOfCols(resultTabulaAA) )
# 
# 
# 
# ### Compare two rows and if they're the same, call merging function ###
# 
# # matchCols is output of indexOfCols
# MergeSequences <- function(row1,row2,data,matchCols){
#   # check if all the parameters are the same so we can merge the two rows together
#   for (name in 1:5) 
#     {
#     print(data[[row1,matchCols[name]]])
#     print(data[[row2,matchCols[name]]])
#     if((is.na(data[[row1,matchCols[name]]])) || (is.na(data[[row2,matchCols[name]]]))) 
#       { 
#           print("missing values")
#           return(data)
#       }
#         if(data[[row1,matchCols[name]]] == data[[row2,matchCols[name]]])
#         { 
#           print("this row ok")
#         }
#         else {
#           print("not a match")
#           return(data)}
#   }
#   # if everything was okay and we went through all the parameters, we can
#   # say that the rows fit and merge them together
#   data <- mergeTwoRows(row1,row2,data)
#   print("good")
#   return(data)
# }
# 
# 
# ### Function to merge the two rows together ###
# mergeTwoRows <- function(row1,row2,data) {
#   print("hey sexy")
#   data$Umi.count[row1] = data$Umi.count[row1] + data$Umi.count[row2]
#   data$Read.count[row1] = data$Read.count[row1] + data$Read.count[row2]
#   data$Umi.proportion[row1] = data$Umi.proportion[row1] + data$Umi.proportion[row2]
#   data$Read.proportion[row1] = data$Read.proportion[row1] + data$Read.proportion[row2]
#   dataNew <- data[-row2,]
#   print("how's it goin?")
#   return(dataNew)
# }

### Main function calling all the other functions ####
# main <- function(data) {
#   indexOfCols <- findIndexOfCols(data)
#   print("here")
#   resultTab <- editNewTable(data,indexOfCols)
#   view(resultTab)
#   return(resultTab)
# }
# 
# hee <- main(new_table)
# maybeRightResult <- main(resultTable)








### function to sort the table according to AA seq ###
# sortAAseq <- function(table) {
#   sortedData <- table %>% arrange(CDR3.amino.acid.sequence)
#   
#   return(sortedData)
# }
# 
# sortedTable <- sortAAseq(resultTable)


# editTable <- function(sortedTable,matchCols) {
#   SAMPLEBIONAME <- matchCols[6]
#   resultTabb <- copy(sortedTable)
#   row1 <- 1
#   excluded <- c()
#   while (row1 < nrow(sortedTable)) {
#     seq_name <- sortedTable$CDR3.amino.acid.sequence[row1]
#     if ((sortedTable$CDR3.amino.acid.sequence[row1+1] != seq_name) && (sortedTable$Umi.count[row1] == 1)) {
#       row1 <- row1 + 1
#       excluded = c(excluded,row1)
#       
#       #print(seq_name)
#       next
#     }
#     
#      row2 <- row1 + 1
#       while ((sortedTable$CDR3.amino.acid.sequence[row2] == seq_name) && (row2 < nrow(sortedTable))){
#          
#          if(sortedTable[[row1,SAMPLEBIONAME]] == sortedTable[[row2,SAMPLEBIONAME]]) {
#            sortedTable <- MergeSequences(row1,row2,sortedTable,matchCols)
#            print("*************************************************")
#          }
#          row2 <- row2 + 1
#    
#    }
#     
#     row1 <- row1 + 1
#   }
#       
#   print(excluded)
#   print(length(excluded))
#   return(resultTabb)
# }
# row1 <- 1
# seq_name = "CAAAF"
# sortedTable$CDR3.amino.acid.sequence[row1+1] != seq_name
# excluded <- c()
# lenght
# resultTabbi <- editTable(sortedTable, indexOfCols)



# sortOnlySameAA <- function(resultTable,matchCols) {
#   row1 = 1
#   SAMPLEBIONAME <- matchCols[6]
#   while (row1 < nrow(resultTable)){
#     print(row1)
#     if (resultTable$Used[row1] == TRUE) { 
#       print('ops')
#       row1 = row1 + 1
#       next }
#     seq_name <- resultTable$CDR3.amino.acid.sequence[row1]
#     resultTable$Used[resultTable$CDR3.amino.acid.sequence == seq_name] <- TRUE
#     new_table <- filter(resultTable,resultTable$CDR3.amino.acid.sequence == seq_name)
#     sampleBioNamesHash <- hash()
#     if ((nrow(new_table) > 1) | (MM.xp.3a$Umi.count[row1] > 1)) { 
#       print('hi')
#       for (row in 1:nrow(new_table)) {
#         if(!is.null(sampleBioNamesHash[[seq_name]])) {
#           row2 = sampleBioNamesHash[[seq_name]]
#           MergeSequences(row,row2,data = resultTable,matchCols)
#         }
#         sampleBioNamesHash[[seq_name]] = 1
#         resultTabulaAA[nrow(resultTabulaAA)+1,] <- new_table[row,]
#       }
#     }
#     row1 = row1 + 1
#   }
#   return(resultTabulaAA)
# }
# 
# indexOfCols = findIndexOfCols(resultTable)
# resultTabbi = sortOnlySameAA(resultTable,indexOfCols)


# !is.null(sampleBioNamesHash[["seq_name"]])
# 
# ree <-(paste("resultTabula$",as.String(he[1]),sep=""))
# print(ree)
# for (name in 1:5){
#   for (column in 1:ncol(resultTabulaAA)) { 
#     if (colnames(resultTabulaAA)[column] == neededMatch[name]) { print(column)}
#   }
# }
# 
# re <- read_csv("250203_result_nucleo_seq_filter.csv")
# resultTabulaAA <- read_csv("250210_result_aa_seq_filter.csv")
# resultTabulaAA <- resultTabulaAA[-1]
# 
# resultTabula[nrow(resultTabula)+1,] <- resultTable[1,]
# resultTable <- filter(resultTable, MM.xp.3a$Umi.count > 1)


       