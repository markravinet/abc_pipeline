### ABC PROCESSOR v1.4 ###
# processing script - formats and separates data
# 25/04/2016

rm(list = ls())

suppressMessages(library(psych))
suppressMessages(library(abc))
suppressMessages(library(getopt))

# specify command line options
spec <- matrix(c(
  'in_sim', 'i', 1, 'character', 'specify sum stats filepath',
  'prefix', 'p', 1, 'character', 'specify prefix for outputs',
  'compress', 'c', 0, 'logical', 'compress output',
  'help', 'h', 0, 'logical', 'display helpful help'
), ncol = 5, byrow = T)

# set command line options
opt = getopt(spec)

# show help if asked for
if (!is.null(opt$help)) {
  cat(paste(getopt(spec, usage=T),"\n"));
  q();
}

# sort out compress
compress <- TRUE
if (is.null(opt$compress)) {
  compress <- FALSE
} 

# set variables for call
in_sim <- opt$in_sim
prefix <- opt$prefix

### INTERNAL FUNCTIONS ###
abc_prep <- function(myData){
  # function to prep abc data for analysis
  param <- myData[, grepl("^p_", names(myData))]
  names(param) <- sub("^p_", "", names(param))
  sust <- myData[, grepl("^s_", names(myData))]
  names(sust) <- sub("^s_", "", names(sust))
  out <- list(param = param, sust = sust)
  return(out)
}

### DUMMY VARIABLES ###
# in_sim <- "./test_dir4/I_growth_test.gz"
# prefix <- "./test_dir/I_growth_test"

## Simulated statistics ##
if(substr(in_sim, nchar(in_sim)-1, nchar(in_sim)) == "gz"){
  nlines <- as.numeric(system(paste0("zcat <", in_sim, " | wc -l"), intern = T))
} else{
  nlines <- as.numeric(system(paste0("cat ", in_sim, " | wc -l"), intern = T))
}
# read in data
myData <- read.table(in_sim, header = T, nrows = nlines)
print("Read in summary statistics")
# remove all Thomson statistics
myData <- myData[, grep("thomson",names(myData), invert = TRUE)]
# remove p_npop
myData <- myData[, grep("p_npop",names(myData), invert = TRUE)]
# run abc_prep
myData <- abc_prep(myData)

# sort data 
# summary statistics, remove all NANs
ss <- myData$sust
ss[is.nan(as.matrix(ss))] <- 0

# params, remove all columns with sum of 0
param <- myData$param
head(param)
str(param)
param <- param[, apply(param, 2, function(z) sum(z) != 0)]
rm(myData) # to save memory

## subset summary stats
sumstat <- colnames(ss)
sumstat <- sub("_mean", "", grep("_mean", sumstat, value = T))
# use indices to subset
sumstat_12ss <- sumstat[c(1:5, 10:12, 17:20)]
sumstat_20ss <- sumstat[c(1:12, 16:20, 25, 27, 29)]
sumstat_29ss <- sumstat
# get indices
sumstat_12ss_idx <- grep(paste(sumstat_12ss, collapse = "|"), names(ss))
sumstat_20ss_idx <- grep(paste(sumstat_20ss, collapse = "|"), names(ss))

if(compress == TRUE){
  # write out both params and sumstats
  print("Writing params")
  write.table(param, gzfile(paste0(prefix, ".params")), quote = FALSE)
  # write out each of the three sumstat files
  print("Writing 12ss")
  write.table(ss[, sumstat_12ss_idx], gzfile(paste0(prefix, "_12ss.ss.gz")), quote = FALSE)
  print("Writing 20ss")
  write.table(ss[, sumstat_20ss_idx], gzfile(paste0(prefix, "_20ss.ss.gz")), quote = FALSE)
  print("Writing 29ss")
  write.table(ss, gzfile(paste0(prefix, "_29ss.ss.gz")), quote = FALSE)
} else{
  # write out both params and sumstats
  print("Writing params")
  write.table(param, paste0(prefix, ".params"), quote = FALSE)
  # write out each of the three sumstat files
  print("Writing 12ss")
  write.table(ss[, sumstat_12ss_idx], paste0(prefix, "_12ss.ss"), quote = FALSE)
  print("Writing 20ss")
  write.table(ss[, sumstat_20ss_idx], paste0(prefix, "_20ss.ss"), quote = FALSE)
  print("Writing 29ss")
  write.table(ss, paste0(prefix, "_29ss.ss"), quote = FALSE)
}


