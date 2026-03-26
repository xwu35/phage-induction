#!/usr/bin/env Rscript

# this script is used to select contigs that are identfied by at least two methods
# usage
# Rscript select_viral_contigs.R -g genomad_summary.tsv -c cenote-taker2_summary.tsv -v virsorter_final_score.tsv -h check_quality_summary.tsv -m mmseq2_result.txt -l genomad_plasmid_summary.txt -o contigListNoControl.txt -p genomad_provirus_contig.txt

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))

option_list = list(
  make_option(c("-g", "--file1"), type="character", default=NULL, 
              help="geNomad file name", metavar="character"),
  make_option(c("-c", "--file2"), type="character", default=NULL, 
              help="cenote-taker3 file name", metavar="character"),
  make_option(c("-v", "--file3"), type="character", default=NULL, 
              help="virsorter2 file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output virus name [default= %default]", metavar="character"),
  make_option(c("-p", "--out2"), type="character", default="provirus_out.txt", 
              help="output provirus name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE);
opt = parse_args(opt_parser);

# check geNomad input
if (is.null(opt$file1)){
  print_help(opt_parser)
  stop("path to geNomad virus summary must be supplied (input file).n", call.=FALSE)
}

# check Cenote-taker2 input
if (is.null(opt$file2)){
  print_help(opt_parser)
  stop("path to Cenote-taker2 summary files must be supplied (input file).n", call.=FALSE)
}

# check Virsorter2 input
if (is.null(opt$file3)){
  print_help(opt_parser)
  stop("path to Virsorter2-pass1 final score must be supplied (input file).n", call.=FALSE)
}

# function to read non-empty files
read_if_not_empty <- function(file, header=T) {
  if (file.info(file)$size > 0) {
    df <- read.csv(file, header=header, sep = "\t")
    return(df)
  } else {
    message(paste("Skippig empty file:", file))
    return(NULL)
  }
}

#------------ geNomad  
# read in genomad summary file, all contigs were retained
genomadQuality <- read_if_not_empty(opt$file1) 

if (is.null(genomadQuality)) {
  genomad <- NULL
  genomad_provirus_list <- NULL
  
  # create an empty geNomad provirus contigs file (same as the touch command in bash)
  file.create(opt$out2)
} else {
  # separate seq_name column in case of the presence of provirus
  # fill = "right" to avoid warning when there is no provirus (coz the format is contig_262186|provirus_3384_9653, so fill with missing values on the right)
  genomadQualitySep <- genomadQuality %>% 
    separate(seq_name, into=c("Contig", "temporary"), sep="\\|", remove = FALSE, fill = "right")
  
  # get the contig names
  genomad <- genomadQualitySep %>% 
    select(Contig)
  
  # provirus identified by genomad
  genomad_provirus <- genomadQualitySep %>% filter(topology == "Provirus") 
  genomad_provirus_list <- genomad_provirus$Contig
  
  # write out geNomad provirus contigs (the original one)
  write.table(genomad_provirus$seq_name, opt$out2,
              sep="\t", quote=F, col.names=F, row.names = F)
}

#------------- Cenote-Taker3 
# read in cenote-taker3 summary file, all contigs were retained
cenoteQuality <- read_if_not_empty(opt$file2)

if (is.null(cenoteQuality)) {
  cenote <- NULL
} else {
  # get the contig names
  cenote <- cenoteQuality %>% 
    select(input_name) %>%
    # the contig names are line C003_1_6MO_k141_222 flag=3 multi=79.0153 len=..
    # extract the strings before the first space
    mutate(Contig = str_extract(input_name, "^[^ ]+")) %>%
    select(Contig)
}

#---------- VirSorter2
# use the quality file to select contigs 
virsorter2Quality <- read_if_not_empty(opt$file3)

if (is.null(virsorter2Quality)) {
  virsorter2 <- NULL
} else {
  # select contigs match the criteria
  virsorter2QualitySep= virsorter2Quality %>%
    separate(seqname, into=c("Contig", "temporary"), sep="\\|\\|") # the double pipe symbol needs to be escaped with double backslashes because '|' is a special character in regular expression
  
  # get the contig names
  virsorter2 <- virsorter2QualitySep %>% select(Contig)
}

#-------------- combine all the identified contigs
AllMethods = rbind(genomad, cenote, virsorter2) 

# if there are identified viral contigs, filter out provirus and plasmid
if (is.null(AllMethods)) {
  # create an empty output
  file.create(opt$out)
} else {
  # remove duplicates
  # remove provirus from the contig list
  ViralContigs <- AllMethods%>%
    distinct() %>% 
    filter(!Contig %in% genomad_provirus_list) 
  
  # write out viral contigs without geNomad provirus
  write.table(ViralContigs, opt$out,
              sep="\t", quote=F, col.names=F, row.names = F)
}
