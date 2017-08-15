### ABC posterior probability estimation ###

rm(list = ls())

suppressMessages(library(abc))
suppressMessages(library(psych))
suppressMessages(library(getopt))

# specify command line options
spec <- matrix(c(
  'dir', 'd', 1, 'character', 'specify data dir',
  'in_obs', 'i', 1, 'character', 'specify observed stats',
  'sum_sub', 's', 0, 'integer', 'optional subset for sumstats - should be 12, 20 or 29 - default 29',
  'nsample', 'n', 0, 'integer', 'number of samples to take for PCA',
  'pca_out', 'p', 1, 'character', 'specify prefix for pca',
  'pp_out', 'o', 1, 'character', 'specify prefix for posterior',
  'help', 'h', 0, 'logical', 'display helpful help'
), ncol = 5, byrow = T)

# set command line options
opt = getopt(spec)

# show help if asked for
if (!is.null(opt$help)) {
  cat(paste(getopt(spec, usage=T),"\n"));
  q();
}

# set calls
dir <- opt$dir
in_obs <- opt$in_obs
#pca_sample <- opt$pca_n
pca_out <- opt$pca_out
posterior_out <- opt$pp_out

# subsample of summary statistics 
# check to ensure the same
if(!nss == "ss12" | !nss == "ss20" | !nss == "ss29" | is.null(opt$sum_sub)){
  print("Summary subsets must be 12, 20 or 29")
  q()
}

# set samples for pca
nsample <- opt$nsample
if (!is.null(opt$nsample)) {
  nsample <- 50000
}

# dummy args 
# dir <- "./test_dir2"
# in_obs <- "./obs/2pop_autosomes_030315_obs.txt"
# pca_out <-  "./test_dir/model_test_pca"
# posterior_out <-  "./test_dir/model_test_pp"
# nss <- "ss12"
# nsample <- 50000

# read in all files
files <- switch(nss,
                ss12 = system(paste0("ls ", dir, "/*12ss.ss*"), intern = TRUE),
                ss20 = system(paste0("ls ", dir, "/*20ss.ss.*"), intern = TRUE),
                ss29 = system(paste0("ls ", dir, "/*20ss.ss.*"), intern = TRUE))
# get nrows
nlines <- as.numeric(system(paste0("cat ", files[1], " | wc -l"), intern = T))
data <- lapply(files, read.table, sep = " ", header = T, colClasses = "numeric",
               nrows = nlines, comment.char = "")
# finalise sumstats
sprintf("Read in %s", basename(files))
model_names <- unlist(lapply(strsplit(basename(files), "_"), function(z) paste(z[1], z[2], z[3], sep = "_")))
sprintf("Models: %s", model_names)
names(data) <- model_names
data <- lapply(data, function(z){
  z[is.nan(as.matrix(z))] <- 0
  z
}) 
sumstat <- do.call(rbind, data)
rownames(sumstat) <- NULL
print("Sumstats formatted")

# set index
indx <- sapply(1:length(data), function(z){
  rep(names(data)[z], nrow(data[[z]]))
})
index <- as.vector(indx)
print("Observed summary statistics read in")

# import and process observed data
obs_ss <- read.delim(in_obs, sep = "")
# dfirst remove thomson stats
obs_ss <- obs_ss[, grep("thomson", colnames(obs_ss), invert = T)]

# choose the obs subset
obs_ss <- switch(nss, 
                 ss12 = obs_ss[, c(1:5, 10:12, 17:20)], 
                 ss20 = obs_ss[, c(1:12, 16:20, 25, 27, 29)],
                 ss29 = obs_ss)
names(obs_ss) <- sub("s_", "", names(obs_ss))
# generate mean, variance, skew, kurtosis
obs_ss <- as.vector(apply(obs_ss, 2, function(z){
  c(mean(z, na.rm = T),
    var(z, na.rm = T),
    skew(z, na.rm = T),
    kurtosi(z, na.rm = T))
}))

# perform PCA to check accepted sum stats approach obs
# sample 50 000 datasets from each of the models
print("Performing PCA")
sprintf("Sampling: %s", nsample)
sample_sumstat <- lapply(data, function(z){
  z[sample(1:nrow(z), nsample, replace = FALSE), ]
})
sample_sumstat <- do.call(rbind, sample_sumstat)
rownames(sample_sumstat) <- NULL
pca_data <- rbind(obs_ss, sample_sumstat)
pca_data[is.na(pca_data)] <- 0
sprintf("Number of datasets: %s", nrow(pca_data))
pca_model <- prcomp(pca_data, scale = TRUE)
pca_idx <- c("obs", sort(rep(model_names, nsample)))

pca_data <- data.frame(model = pca_idx, pca_model$x,
                       type = c("obs", rep("sim", length(pca_idx)-1)))
write.table(pca_data, paste0(pca_out, "_data.txt"), quote = FALSE, sep = " ",
            row.names = FALSE)
# set up for screeplot
pca_imp <- t(summary(pca_model)$importance)
pca_imp <- data.frame(PC = rep(1:nrow(pca_imp)), imp = pca_imp[, 2])
rownames(pca_imp) <- NULL
write.table(pca_imp, paste0(pca_out, "_importance.txt"),
            quote = FALSE, sep = " ",
            row.names = FALSE)

# posterior estimation
print("Begin posterior estimation")
# perform posterior probability estimation
tolerance_range <- c(0.001, 0.005, 0.01, 0.03)
tol_pp <- lapply(tolerance_range, function(z){
  test <- postpr(target = obs_ss, index, sumstat, tol = z, method = "neuralnet")
  summ_test <- summary(test)
  # extract from output
  rej <- summ_test[[1]]$Prob
  mnl <- summ_test[[2]]$Prob
  list(rej, mnl)
})
# sort out for output
saveRDS(tol_pp, file =  paste0(posterior_out, "_Robj", ".rds"))
rejection_pp <- do.call(rbind, lapply(tol_pp, function(z) z[[1]]))
multinom_pp <- do.call(rbind, lapply(tol_pp, function(z) z[[2]]))
colnames(rejection_pp) <- paste0(colnames(rejection_pp), "_rej")
colnames(multinom_pp) <- paste0(colnames(multinom_pp), "_nn")
# finalise data.frame
pp <- data.frame(tol = tolerance_range, rejection_pp, multinom_pp)

# write out
write.table(pp, paste0(posterior_out, "_pp.txt"),
            quote = FALSE, sep = " ",
            row.names = FALSE)

print("Done")
