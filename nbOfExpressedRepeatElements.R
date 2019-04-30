#This code is used to analyze the number of expressed repeat elements. 
#It uses the gene annotation to filter the intergenic regions. 
#If the gene count are supplied, it will normalized using the total read counts of 
#the genes, else using the total read counts of the repeats.
#If lnRNA annotation are supplied, it will remove repeat elements overlapping these elements.
#It will print a graph for each class of repeat elements and for each sample using a cutoff
#for the expression

library(data.table) 
library(rtracklayer)
library(optparse)
library(ggplot2)

## Functions
checkRequiredArguments <- function(opts, parser){
  # This function check that the required arguments are
  # not missing
  for (o in 1:length(parser@options)){
    current_arg = parser@options[[o]]
    if (grepl("[REQUIRED]",current_arg@help) & is.null(unlist(opt[current_arg@dest]))){
      stop("Required arguments missing. Look at the help.")
    }
    if (grepl("(bed format)",current_arg@help) & 
        (!grepl(".bed", unlist(opt[current_arg@dest])))){
      stop("Gene annotation should be a bed format.")
    }
  }
}

filter_file <- function(dir, sam, pat){ 
  # This function will return the path of a file, using
  # of the name of the sample and the file extention
  pattern1 = paste0("*",sam,"*")
  pattern2 = paste0(pat,"$")
  pattern = paste0(pattern1,pattern2)
  file = list.files(dir, full.names=T, recursive=T, pattern=glob2rx(pattern))
  if (length(file) > 1){
    stop(paste0("The name sample ", sam ," has multiples matches files."))
  } 
  if (length(file) == 0){
    stop(paste0("No matched name with pattern ", pattern, "."))
  }
  return(file)
}

find_condition <- function(cond, sam){
  # This function match the condition with the sample name
  # example : 
  # We have two conditions : K36M, WT
  # We give the sample name : s1_K36M
  # It will return the matching condition : K36M
  res = sapply(cond, function (y) sapply(sam, function (x) grepl(y, x)))
  current_cond = names(res[res == TRUE])
  current_cond = unlist(strsplit(current_cond, "[.]"))[1]
  return(current_cond)
}


give_TRC <- function(file, annotation, chrm_to_remove, isBED){
  # This function gives Total Read Count 
  # for genes or repeats depending of the file and annotation given.
  print(file)
  count = read.table(file, header = T)
  nb_input = dim(count)[1]
  if (!is.null(chrm_to_remove)){
    if (isBED == TRUE){
      to_remove = unique(annotation[seqnames(annotation) %in% chrm_to_remove,]$name)
    } else {
      to_remove = unique(annotation[seqnames(annotation) %in% chrm_to_remove,]$gene_id)
    }
    count = count[!(count$Geneid %in% to_remove),]
  }
  nb_output = dim(count)[1]
  print(paste0(nb_input - nb_output, " elements removes from ", file))
  print(sum(as.numeric(count[,length(count)])))
  print(dim(count)[1])
  total = sum(count[,length(count)])/1000000
  return(total)
}

giveResult <- function(sam, repeat_anno, gene_anno, lnRNA_anno, chrm_to_remove){
  # This function process different steps
  # 1. We are exlcuding the gene coding regions 
  # (ie we are performing the analysis only on intergenic regions)
  # 2. If lnRNA annotation is furnished, we are exlcluding elements overlapping lnRNA.
  # 3. If we ask to remove a particular chromosome, we are excluding it
  # 4. We normalize the expression by total read count and the width of the element.
  gene_anno = reduce(gene_anno)
  gene_anno = gene_anno + 10000
  
  repeat_expression = read.table(sam@repeat_path, header = T)
  if (!is.null(chrm_to_remove)){
    repeat_expression = repeat_expression[!(repeat_expression$Chr %in% chrm_to_remove),]
  }
  repeat_expression = GenomicRanges::makeGRangesFromDataFrame(repeat_expression, keep.extra.columns = TRUE)
  overlaps = findOverlaps(repeat_expression, gene_anno, ignore.strand=TRUE)
  repeat_expression = repeat_expression[-queryHits(overlaps)]
  
  if(!is.null(lnRNA_anno)){
    lnRNA_anno = reduce(lnRNA_anno)
    overlaps = findOverlaps(repeat_expression, lnRNA_anno, ignore.strand=TRUE)
    repeat_expression = repeat_expression[-queryHits(overlaps)]
  }
  
  repeat_expression = as.data.frame(repeat_expression)
  repeat_expression[,length(repeat_expression)] = (repeat_expression[,length(repeat_expression)])/sam@total_read_count
  repeat_expression[,length(repeat_expression)] = repeat_expression[,length(repeat_expression)]/(repeat_expression$width * 0.001)
  return(repeat_expression)
}

give_tab <- function(sams,list_elements, lim){
  # Giving a list of elements (SINE, LINE, LTR...), this function
  # return the number of expressed elements (FPKM > lim) for 
  # each sample
  df = data.frame(sample = samples, count_element = NA , condition = NA)
  for (s in 1:length(sams)){
    print(sams[[s]]@name)
    res = sams[[s]]@result
    print(dim(res))
    print(dim(res[res$Geneid %in% list_elements,]))
    res = res[res[,length(res)] > lim ,]
    
    res = res[res$Geneid %in% list_elements,]
    nb = dim(res)[1]
    df[df$sample == sams[[s]]@name,]$count_element = nb
    df[df$sample == sams[[s]]@name,]$condition = sams[[s]]@condition
  }
  return(df)
}



giveGraph <- function(sams, pdffile, class_list, is_boxplot){
  #This function save a graph (barplot or boxplot).
  
  pdf(pdffile)
  for (c in 1:length(class_list)){
    type_of_element = class_list[c]
    list_of_elements = repeat_annotation[repeat_annotation$repClass == type_of_element,]$gene_id
    list_lim = c(0.5, 1, 5)
    for (l in 1:length(list_lim)){
      lim = list_lim[l]
      tab = give_tab(sams,list_of_elements, lim)
      print(tab)
      if (!isTRUE(is_boxplot)){
        
        g <-  ggplot(tab, aes(x=tab$sample, y=tab$count_element, color=tab$condition)) +
          geom_point() + 
          xlab("SAMPLES") + ylab(paste0("Number of ", type_of_element, " expressed elements")) + ggtitle(paste0("FPKM >",lim)) +
          theme(text = element_text(size=10),
                axis.text.x = element_text(angle=45, hjust=1),
                plot.margin = margin(2, 2, 2, 2, "cm")) 
        print(g)
        
      } else {
        g <-ggplot(tab, aes(x=tab$condition, y=tab$count_element, fill=tab$condition)) +
          geom_boxplot() + xlab("CONDITION") + 
          ylab(paste0("Number of ", type_of_element, " expressed elements")) + ggtitle(paste0("FPKM >",lim)) +
          theme(text = element_text(size=10),
                axis.text.x = element_text(angle=45, hjust=1),
                plot.margin = margin(2, 2, 2, 2, "cm")) 
        
        print(g)
      }
    }
  }
  dev.off()
  
}

###

option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="[REQUIRED] file directory (where are the count files).", dest = "dir", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", dest = "output_file", default=NULL, 
              help="[REQUIRED] pdf file.", metavar="character"),
  make_option("--repeat_annotation", type="character", default=NULL, 
              help="[REQUIRED] repeat annotation (gtf format).", metavar="character"),
  make_option("--gene_annotation", type="character", default=NULL, 
              help="[REQUIRED] gene annotation (bed format).", metavar="character"),
  make_option("--lnRNA_annotation", type="character", default=NULL, 
              help="lnRNA annotation (gtf format).", metavar="character"),
  make_option("--TEsuffix", type="character", default=NULL, 
              help="[REQUIRED] suffix for the counts of Transposable Elements.", metavar="character"),
  make_option("--Gsuffix", type="character", default=NULL, 
              help="[REQUIRED] suffix for the counts of genes.", metavar="character"),
  make_option("--samples", type="character", default=NULL, 
              help="[REQUIRED] sample names.", metavar="character"),
  make_option("--removeChrm", type="character", default=NULL, 
              help="chromosome to remove.", metavar="character"),
  make_option("--class", type="character", default="all", 
              help="class of elements (LTR, LINE..). Default is all.", metavar="character"),
  make_option("--conditions", type="character", 
              default=NULL, help="conditions to analyse (should be included in the name of the samples)", metavar="character"),
  make_option("--boxplot", action="store_true", default=FALSE, 
              help="Plot the graph as a boxplot.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



checkRequiredArguments(opts, opt_parser)

#parse the arguments

input_path = opt$dir
output_file = opt$output_file
repeat_annotation = opt$repeat_annotation
gene_annotation = opt$gene_annotation
lnRNA_annotation = opt$lnRNA_annotation


TEsuffix = opt$TEsuffix
Gsuffix = opt$Gsuffix


samples = opt$samples
removeChrm = opt$removeChrm
conditions = opt$conditions
class = opt$class
boxplot = opt$boxplot

if (!is.null(class)){
  class = unlist(strsplit(class, ","))
}

if (!is.null(removeChrm)){
  removeChrm = unlist(strsplit(removeChrm, ","))
} 

samples = unlist(strsplit(samples, ","))

if (!is.null(conditions)){
  conditions = unlist(strsplit(conditions, ","))
}


##################

## 1. Load the annotations
is_bed = grepl(".bed", gene_annotation)
repeat_annotation = import(repeat_annotation) 
gene_annotation = import(gene_annotation)

if (!is.null(lnRNA_annotation)){
  lnRNA_annotation = import(lnRNA_annotation)
}


## 2. create class sample with info for each sample

setClass("sample", slots=list(name="character", 
                              condition="character",
                              repeat_path = "character", 
                              gene_path = "character",  
                              total_read_count="numeric", 
                              result="data.frame"))


## 3. list of sample objects and enter the information inside

list_of_samples = vector("list", length(samples))

for (j in 1:length(samples)){
  current_sample = samples[j]
  print(current_sample)
  current_sample_object <- new("sample", 
                             name = current_sample,
                             condition = find_condition(conditions, current_sample),
                             repeat_path = filter_file(input_path, current_sample, TEsuffix))
  
  if (!is.null(Gsuffix)){
    current_sample_object@gene_path = filter_file(input_path, current_sample, Gsuffix)
    current_sample_object@total_read_count = give_TRC(current_sample_object@gene_path, gene_annotation, removeChrm, is_bed)
  } else {
    current_sample_object@total_read_count = give_TRC(current_sample_object@repeat_path, repeat_annotation, removeChrm, FALSE)
  }
  list_of_samples[[j]] <- current_sample_object
}


for (s in 1:length(list_of_samples)){
  current_sample_object = list_of_samples[[s]]
  current_sample_object@result = giveResult(current_sample_object, repeat_annotation, 
                                            gene_annotation, lnRNA_annotation, removeChrm)
  list_of_samples[[s]] = current_sample_object
}


giveGraph(list_of_samples, output_file, class, boxplot)






