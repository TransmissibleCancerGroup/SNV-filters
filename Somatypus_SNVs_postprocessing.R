######## Somatypus Variant Calling Post-Processing #############
###################### Version 4.1 #############################

# mrs72
# Last Update - 23/11/2017 ##
#############################

library(stringr)
library(GenomicRanges)


## 1. Load SNVs from Somatypus
##############################
All.SNVs <- read.table("/Your/Somatypus/Directory/Output.vcf")[,-c(3,9)]

BAF_from_Platypus <- function(VCF_Input, sample){
  of <- str_split_fixed(VCF_Input[,sample],":", 6)[,6]
  these <- str_split_fixed(VCF_Input[,sample],":", 6)[,5]
  VCF_Input[,sample] <- paste(of, these, sep="/")
  return(VCF_Input)
}

# Extract support/total reads for SNPs in each sample
for (i in colnames(All.SNVs)[8:ncol(All.SNVs)]){
  All.SNVs <- BAF_from_Platypus(All.SNVs,i)
}


## 2. Post-processing
#####################

# a. Reference Errors
hom.reference <- function(VCF_Input,min_threshold, name_of_reference_sample){
  
  # a. obtain baf
  int_baf <- as.numeric(str_split_fixed(VCF_Input[,name_of_reference_sample],"/",2)[,1])/
    as.numeric(str_split_fixed(VCF_Input[,name_of_reference_sample],"/",2)[,2])
  
  # b. threshold
  homozyg <- which(int_baf>=min_threshold)
  cat("\n Number of potential reference errors found and removed:", length(homozyg))
  VCF_Input <- VCF_Input[-homozyg,]
  
  # c. output
  return(VCF_Input)
}
All.SNVs <- hom.reference(All.SNVs, 0.9)

# b. Strand-biases
strand.bias <- function(VCF_Input, variant_support_edge, min_reads, min_percentage){
  
  # 1. Build 3 columns: N/TC, NF/N, NR/N from Info-field
  
  cat("\n Isolating strand-Coverages...")
  N   <- str_split_fixed(as.character(VCF_Input[,"INFO"]), ";", 20)[,18]
  NF  <- str_split_fixed(as.character(VCF_Input[,"INFO"]), ";", 20)[,8]
  NR  <- str_split_fixed(as.character(VCF_Input[,"INFO"]), ";", 20)[,9]
  TC  <- str_split_fixed(as.character(VCF_Input[,"INFO"]), ";", 20)[,15]
  N_TC <- paste(str_split_fixed(N, "=", 2)[,2], str_split_fixed(TC, "=", 2)[,2], sep="/")
  NF_N <- paste(str_split_fixed(NF, "=", 2)[,2], str_split_fixed(N, "=", 2)[,2], sep="/")
  NR_N <- paste(str_split_fixed(NR, "=", 2)[,2], str_split_fixed(N, "=", 2)[,2], sep="/")
  names <- colnames(VCF_Input)
  VCF_Input <- cbind(VCF_Input, N_TC, NF_N, NR_N)
  colnames(VCF_Input) <- c(names, "Support Coverage", "Forward Support", "Reverse Support")
  #VCF_Input <- VCF_Input[,-7]
  
  # 2. Build two position vectors: which total support is ≤ or > variant_support_edge
  cat("\n Filtering...")
  above <- which(as.numeric(str_split_fixed(VCF_Input[,"Support Coverage"], "/", 2)[,1])>variant_support_edge)
  below_or_equal <- which(as.numeric(str_split_fixed(VCF_Input[,"Support Coverage"], "/", 2)[,1])<=variant_support_edge)
  
  # 3. Apply thresholds to positions in either vector, but on the whole set!
  
  ## above: which Forward Support or Reverse Support are < 20% ?
  above.fw <- as.numeric(str_split_fixed(as.character(VCF_Input[above,"Forward Support"]),"/",2)[,1])/as.numeric(str_split_fixed(as.character(VCF_Input[above,"Forward Support"]),"/",2)[,2])
  above.rv <- as.numeric(str_split_fixed(as.character(VCF_Input[above,"Reverse Support"]),"/",2)[,1])/as.numeric(str_split_fixed(as.character(VCF_Input[above,"Reverse Support"]),"/",2)[,2])
  above <- above[-c(which(above.fw<min_percentage | above.rv<min_percentage))]
  
  ## below: which Forward Support or Reverse Support are < 2 reads ?
  below_or_equal.fw <- as.numeric(str_split_fixed(as.character(VCF_Input[below_or_equal,"Forward Support"]),"/",2)[,1])
  below_or_equal.rv <- as.numeric(str_split_fixed(as.character(VCF_Input[below_or_equal,"Reverse Support"]),"/",2)[,1])
  below_or_equal <- below_or_equal[-c(which(below_or_equal.fw<min_reads | below_or_equal.rv<min_reads))]
  
  ## variants which pass the tresholds
  pass <- sort(union(above,below_or_equal))
  
  # 4. Count the numbers of removed variants (and percentage of the total calls)
  cat("\n \n Variants removed due to percentage underpassing: ", length(which(above.fw<min_percentage | above.rv<min_percentage)), "of", length(N))
  cat("\n Variants removed due to read underpassing: ", length(which(below_or_equal.fw<min_reads | below_or_equal.rv<min_reads)), "of", length(N))
  cat("\n Total removed: ", round(c(length(which(above.fw<min_percentage | above.rv<min_percentage))+length(which(below_or_equal.fw<min_reads | below_or_equal.rv<min_reads)))/length(N),4)*100, "%")
  
  # 5. Output
  VCF_Input <- VCF_Input[pass,]
  return(VCF_Input)
}
All.SNVs <- strand.bias(All.SNVs, 10, 2, 0.2)

# c. Simple repeats
simple.repeats <- function(VCF_Input, repeat_distance){
  
  devil.repeats <- read.table("/Users/ms37/Desktop/Data/Info-Files/Devil_Repeats_translated.txt",header=T)
  require(GenomicRanges)
  
  # 1. Build GenomicRanges +- repeat_distance BP
  Input.Ranges <- GRanges(seqnames = Rle(VCF_Input[,1]),
                          ranges = IRanges(start = as.integer(VCF_Input[,2]),
                                           end = as.integer(VCF_Input[,2])))
  
  Repeat.Ranges <- GRanges(seqnames = Rle(devil.repeats[,1]),
                           ranges = IRanges(start = as.integer(devil.repeats[,2]-repeat_distance),
                                            end = as.integer(devil.repeats[,3])+repeat_distance))
  
  # 2. Match with positions in VCF_Input
  Overlaps <- findOverlaps(Input.Ranges, Repeat.Ranges)
  Overlaps <- as.matrix(Overlaps)
  colnames(Overlaps) <- c("Set", "Repeat")
  
  # 3. Remove the ones which match
  cat("\n Total removed: ", round(length(unique(Overlaps[,1]))/nrow(VCF_Input),4)*100, "%")
  VCF_Input <- VCF_Input[-unique(Overlaps[,1]),]
  
  # 4. Output
  rm(devil.repeats)
  return(VCF_Input)
}
All.SNVs <- simple.repeats(All.SNVs, 5)

# d. Within 500bp of both contig-ends
contig.ends <- function(VCF_Input, gap_distance){
  
  devil.contig.gaps <- read.table("/Users/ms37/Desktop/Data/Info-Files/Devil_Contig_Gaps.txt", header=T)
  require(GenomicRanges)
  
  # 1. Build GenomicRanges +- repeat_distance BP
  Input.Ranges <- GRanges(seqnames = Rle(VCF_Input[,1]),
                          ranges = IRanges(start = as.integer(VCF_Input[,2]),
                                           end = as.integer(VCF_Input[,2])))
  
  Gap.Ranges <- GRanges(seqnames = Rle(devil.contig.gaps[,1]),
                           ranges = IRanges(start = as.integer(devil.contig.gaps[,2]-gap_distance),
                                            end = as.integer(devil.contig.gaps[,3])+gap_distance))
  
  # 2. Match with positions in VCF_Input
  Overlaps <- findOverlaps(Input.Ranges, Gap.Ranges)
  Overlaps <- as.matrix(Overlaps)
  colnames(Overlaps) <- c("Set", "Gap")
  
  # 3. Remove the ones which match
  cat("\n Total removed: ", round(length(unique(Overlaps[,1]))/nrow(VCF_Input),4)*100, "%")
  VCF_Input <- VCF_Input[-unique(Overlaps[,1]),]
  
  # 4. Output
  rm(devil.contig.gaps)
  return(VCF_Input)
}
All.SNVs <- contig.ends(All.SNVs, 500)

# e. Within 1000bp of scaffold ends
supercontig.ends <- function(VCF_Input){
  
  supercontig.ends <- read.table("/Users/ms37/Desktop/Data/Info-Files/Region-Files/Devil7.1_supercontigs_1000bp_off_borders.txt", header=F)
  require(GenomicRanges)
  
  # 1. Build GenomicRanges +- repeat_distance BP
  Input.Ranges <- GRanges(seqnames = Rle(VCF_Input[,1]),
                          ranges = IRanges(start = as.integer(VCF_Input[,2]),
                                           end = as.integer(VCF_Input[,2])))
  
  Ends.Ranges <- GRanges(seqnames = Rle(str_split_fixed(supercontig.ends[,1],":",2)[,1]),
                        ranges = IRanges(start = 1000,
                                         end = as.numeric(str_split_fixed(str_split_fixed(supercontig.ends[,1],":",2)[,2],"-",2)[,2])))
  
  # 2. Match with positions in VCF_Input
  Overlaps <- findOverlaps(Input.Ranges, Ends.Ranges)
  Overlaps <- as.matrix(Overlaps)
  colnames(Overlaps) <- c("Set", "Ends")
  
  # 3. Remove the ones which match
  cat("\n Total removed: ", c(1-round(length(unique(Overlaps[,1]))/nrow(VCF_Input),4))*100, "%")
  VCF_Input <- VCF_Input[unique(Overlaps[,1]),]
  
  # 4. Output
  return(VCF_Input)
  
}
All.SNVs <- supercontig.ends(All.SNVs)