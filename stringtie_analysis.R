#!/usr/bin/env Rscript

# To analyse the output of gffcompare from Stringtie. I use the tmap files that contain 
# the information of the transcripts : the expression level (FPKM) and the assigned code.
# For each code, I plot the number of gene containing transcripts of this code.

library("optparse")
library("ggplot2")

option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="file directory where are the output files from Stringtie (ie .tmap files)", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, 
              help="output file", metavar="character"),
  make_option("--samples", type="character", default=NULL, 
              help="list of samples", metavar="character"),
  make_option("--code", type="character", default=NULL, 
              help="select a specific code to show", metavar="character"),
  make_option("--FPKM", type="numeric", default=0.5, 
              help="limite FPKM to use (default 0.5)"),
  make_option("--conditions", type="character", 
              default=NULL, help="list of conditions (the condition name should be in the name of the file)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#parse the parameters
dir = opt$dir
samples = opt$samples
samples = unlist(strsplit(samples, ","))

if (!is.null(opt$condition)){
  condition = opt$condition
  condition = unlist(strsplit(condition, ","))
} else {
  condition = "condition1"
}

code = opt$code

if (!is.null(opt$FPKM)){
  FPKM_lim = opt$FPKM
} else {
  FPKM_lim = 0.5
}

output_file = opt$output_file

setwd(dirname(output_file))


# start the analysis

#S4 class "code" containing the name of the code, the description, and the result
setClass("code", slots=list(name="character", description="character", result="data.frame"))

give_description_code <- function(code){
  switch(code,
         "=" = {desc <- "complete, exact match of intron chain"},
         "c" = {desc <-"contained in reference (intron compatible)"},
         "k" = {desc <-"containment of reference (reverse containment)"},
         "m" = {desc <-"retained intron(s), full intron chain overlap/match"},
         "n" = {desc <-"retained intron(s), partial or no intron chain match"},
         "j" = {desc <-"multi-exon with at least one junction match"},
         "e" = {desc <-"single exon transflag partially covering an intron, possible pre-mRNA fragment"},
         "o" = {desc <-"other same strand overlap with reference exons"},
         "s" = {desc <-"intron match on the opposite strand (likely a mapping error)"},
         "x" = {desc <-"exonic overlap on the opposite strand (like o or e but on the opposite strand)"},
         "i" = {desc <-"fully contained within a reference intron"},
         "y" = {desc <-"contains a reference within its intron(s)"},
         "p" = {desc <-"possible polymerase run-on (no actual overlap)"},
         "r" = {desc <-"repeat (at least 50% bases soft-masked)"},
         "u" = {desc <-"none of the above (unknown, intergenic)"}
  )
  return(desc)
}




if (is.null(code)){
  code_list = c("=", "c", "k", "m", "n", "j", "e", "o", "s", "x", "i", "y", "p", "r", "u")
   } else {
  codes = unlist(strsplit(code, ","))
  code_list <- vector("list", length(code))
}

files = list.files(dir, pattern = ".gtf.tmap", full.names = T, recursive = T)
print("list of files : ")
print(files)

res = vector("list", length(code_list))

for (e in 1:length(code_list)){
  # For each code, we create an object with the name of the code, its description and a
  # dataframe with, for each sample, the number of genes 
  # that have expressed transcripts (FPKM >= lim) 
  
  current_code = code_list[e]
  dataframe = data.frame(sample = samples, count = NA)
  
  current_code_object <- new("code", 
                             name = current_code, 
                             description = give_description_code(current_code))
  
  for (s in 1:length(dataframe$sample)){
      current_sample = dataframe$sample[s]
      index = grep(current_sample, files)
      current_sample = read.table(files[index], header = T)
      current_sample =  current_sample[current_sample$class_code == current_code & current_sample$FPKM >= FPKM_lim,]
      count = length(unique(current_sample$qry_gene_id))
      dataframe[dataframe$sample == dataframe$sample[s],]$count = count
  }
  slot(current_code_object,"result") <- dataframe
  
  res[[e]] = current_code_object
} 


print(paste0("saving graph result in : ", output_file))
pdf(output_file)

for (e in 1:length(res)){
  # For each code, we create a new dataframe with the information
  # that we want plotted (condition added).
  
  current_object = res[[e]]
  
  dataframe = current_object@result
  current_code = current_object@name
  description = current_object@description
  
  dataframe$case = NA
  for (c in 1:length(condition)){
    current_condition = condition[c]
    index  = grep(current_condition,dataframe$sample)
    dataframe[index,]$case = current_condition
  }
  
  
  dataframe <- within(dataframe, 
                      sample <- factor(sample, 
                                        levels=samples))
  
  
  p<-ggplot(dataframe, aes(x=sample, y=count, fill=case)) +
    geom_bar(stat="identity")+theme_minimal()  +
    ylab(paste0("Number of genes with the code '", current_code, "'")) + 
    xlab(paste0("FPKM >=", FPKM_lim)) + 
    labs(subtitle = description) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
} 

dev.off()

print("analysis ended.")
