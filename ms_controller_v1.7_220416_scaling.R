#!/usr/bin/env
#### msABC controllet v1 ####
# Mark Ravinet, 27/11/14
# nb: for constant pop size

rm(list = ls())

suppressMessages(library(psych))
suppressMessages(library(getopt))

# Popns: 1 - Japan Sea; 2 - Pacific Ocean;
# T1 is PO and AT split
# T2 is PO+AT and JS split

# available models
div_models <- c("I", "CM", "IAM", "IRM", "IARM")

# available growth models
growth_models <- c("constant", "growth", "bottleneck")

# available migration models
mig_models <- c("uniform", "hetero")

# specify command line options
spec <- matrix(c(
  'outfile', 'o', 1, 'character', 'specify output file path',
  'priorfile', 'p', 1, 'character', 'specify prior file path',
  'nind', 'i', 1, 'integer', 'number of individuals',
  'nloc', 'l', 1, 'integer', 'number of loci',
  'nsim', 'n', 1, 'integer', 'number of simulations to run',
  'div_model', 'd', 1, 'character', 'choose divergence model',
  'growth_model', 'g', 1, 'character', 'choose growth model',
  'mig_model', 'm', 1, 'character', 'choose migration model',
  'mem_lim', 'f', 0, 'logical', 'run in memory limit mode', 
  'verbose', 'v', 0, 'logical', 'run in verbose mode (print command)',
  'nproc', 'x', 0, 'integer', 'specify number of cpus to use',
  'help', 'h', 0, 'logical', 'display helpful help',
  'model_choice', 'c', 0, 'logical', 'display model options'
   ), ncol = 5, byrow = T)

# set command line options
opt = getopt(spec)

# show help if asked for
if (!is.null(opt$help)) {
  cat(paste(getopt(spec, usage=T),"\n"));
  q();
}

# turn on verbose mode
verbose <- FALSE
if (!is.null(opt$verbose)) {
  verbose <- TRUE
}

# turn on mem_lim mode
mem_lim <- FALSE
if (!is.null(opt$mem_lim)) {
  mem_lim <- TRUE
}

# show model options if asked for
if (!is.null(opt$model_choice)) {
  cat("Divergence models:", models, "\n",
      "Growth models:", growth_models);
  q();
}

# set variables for call
outfile <- opt$outfile
priorfile <- opt$priorfile
nind <- opt$nind
nloc <- opt$nloc
nrep <- opt$nsim
div_model <- opt$div_model
growth_model <- opt$growth_model
mig_model <- opt$mig_model
# non-command line
index <- rep(paste(div_model, growth_model, mig_model, sep = "_"), nrep)

## For parallel, calculate the number of processors available
nproc <- as.numeric(system("nproc", intern = T))
# turn on mem_lim mode
if (!is.null(opt$nproc)) {
  nproc <- opt$nproc
}

# write something to the screen
sprintf("Running: %s sims", nrep)
sprintf("For %s individuals", nind)
sprintf("For %s loci", nloc)
sprintf("Divergence model is: %s", div_model)
sprintf("Growth model is: %s", growth_model)
sprintf("Migration model is: %s", mig_model)
sprintf("Number of processors is: %s", nproc)

cat("\n")

# DUMMY OPTIONS
# nind <- 40
# nloc <- 1800
# div_model <- "IARM"
# growth_model <- "growth"
# nrep <- 10
# index <- rep(paste(div_model, growth_model, sep = "_"), nrep)
# mig_model <- "uniform"
# priorfile <- "test_priorfile"
# outfile <- "test_outfile"
# nproc <- 20

######### PRIORS ######### 

### BASIC PRIORS ###
# set N0 to arbitrary unit of 50 000
N0 <- 50000
# sample theta from a uniform distribution
theta <- runif(nrep, 0.05, 20)

### Time since divergence ###
# sample time to T1 from a lognormal distribution
T1_real <- rlnorm(nrep, log(1500000), 0.8)
# sample T2 from a lognormal distribution 
T1 <- (T1_real/(4*N0))

### MIGRATION PRIORS ###
if(div_model != "I"){
  # migration priors
  rate <- 20 # migration varies from 0.001 to approx 0.5

  # nb: for reference matrix form is :
  #matrix(c("m11", "m12", "m13", "m21", "m22", "m23",
  #         "m31", "m32", "m33"), nrow = 3, byrow = TRUE)
  # JS and AT migration always set to zero
 
  if(mig_model == "uniform"){
    m12 <- (4*N0)*rexp(nrep, rate)
    m21 <- (4*N0)*rexp(nrep, rate)
  }
  
  if(mig_model == "hetero"){
    a <- runif(nrep, 0, 5)
    b <- runif(nrep, 0, 200)
    c1 <- runif(nrep, 0, 15)
    c2 <- runif(nrep, 0, 15)
    # create a matrix of mig_hyper priors
    mig_hyper <- cbind(a, b, c1, c2)
    m12 <- m21 <- rep(NA, nrep)
  }
}

### MIGRATION TIME PRIORS ###
if(div_model == "IAM" | div_model == "IARM"){
  # uniform prior
  anc_migT1 <- T1*runif(nrep, 0.51, 1)
}

if(div_model == "IRM" | div_model == "IARM"){
  # uniform prior
  rec_migT1 <- T1*runif(nrep, 0.01, 0.5)
}

### POPN SIZE AND GROWTH PRIORS ###
# subpopulation sizes (fractions of N0)
N1 <- runif(nrep, 0.01, 3)
N2 <- runif(nrep, 0.01, 3)

# timing of population growth
# note that in the two population case, scT1 is for JA, scT2 is for PA
# nb - this is the time of the onset of growth
if(growth_model != 'constant'){
  scT1 <- T1*runif(nrep, 0.01, 1)
  scT2 <- T1*runif(nrep, 0.01, 1)
  
  # fraction of population size
  x1 <- runif(nrep, 0, 0.5)
  x2 <- runif(nrep, 0, 0.5)
  
  # solve the equation to get the right params
  bpop1 <- (N0*N1)*x1
  alpha1 <- -(1/scT1)*log(bpop1/N0)
  bpop2 <- (N0*N2)*x2
  alpha2 <- -(1/scT1)*log(bpop2/N0)
  
}

### WRITE OUT PRIORS ###
# Collate and write out priors
# divergence priors
priors <- switch(div_model,
                 I = data.frame(index, theta, T1),
                 CM = data.frame(index, theta, T1, m12, m21),
                 IAM = data.frame(index, theta, T1, m12, m21, anc_migT1),
                 IRM = data.frame(index, theta, T1, m12, m21, rec_migT1),
                 IARM = data.frame(index, theta, T1, m12, m21, anc_migT1, rec_migT1))

# add migration hyper priors (if mig_model = "hetero")
if(mig_model == "hetero"){
  priors <- cbind(priors, mig_hyper)
}

# growth priors
# evaluate if present, if not then don't write out
if(!is.null(growth_model)){
  growth_priors <- switch(growth_model,
                          constant = data.frame(N1, N2),
                          growth = data.frame(N1, N2, x1, x2, alpha1, alpha2, scT1, scT2),
                          bottleneck = data.frame(N1, N2, x1, alpha1, scT1))
  priors <- cbind(priors, growth_priors)
}
write.table(priors, file = priorfile, quote = FALSE, row.names = FALSE)

### PREPARE MS COMMANDS ###

# set msline based on model (most complicated part!)
msline_base <- switch(paste0(div_model, growth_model),
                      Iconstant = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20 0.0",
                                         " -n 1 ", priors$N1, 
                                         " -n 2 ", priors$N2, 
                                         " -ej ", priors$T1, " 1 2 --verbose | sed 1d'"),
                      Igrowth = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20 0.0",
                                       " -g 1 ", priors$alpha1, 
                                       " -g 2 ", priors$alpha2, 
                                       " -eg ", priors$scT1, " 1 0.0",
                                       " -eg ", priors$scT2, " 2 0.0",
                                       " -ej ", priors$T1, " 1 2 --verbose | sed 1d'"),
                      Ibottleneck = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20 0.0",
                                           " -n 2 ", priors$N2,
                                           " -g 1 ", priors$alpha1, 
                                           " -eg ", priors$scT1, " 1 0.0",
                                           " -ej ", priors$T1, " 1 2 --verbose | sed 1d'"),
                      CMconstant = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20", 
                                          " -m 1 2 ", priors$m12, " -m 2 1 ", priors$m21,
                                          " -n 1 ", priors$N1, 
                                          " -n 2 ", priors$N2, 
                                          " -ej ", priors$T1, " 2 1 --verbose | sed 1d'"),
                      CMgrowth = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20",
                                        " -m 1 2 ", priors$m12, " -m 2 1 ", priors$m21,
                                        " -g 1 ", priors$alpha1, 
                                        " -g 2 ", priors$alpha2, 
                                        " -eg ", priors$scT1, " 1 0.0",
                                        " -eg ", priors$scT2, " 2 0.0",
                                        " -ej ", priors$T1, " 1 2 --verbose | sed 1d'"),
                      CMbottleneck = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20",
                                            " -m 1 2 ", priors$m12, " -m 2 1 ", priors$m21,
                                            " -n 2 ", priors$N2, 
                                            " -g 1 ", priors$alpha1, 
                                            " -eg ", priors$scT1, " 1 0.0",
                                            " -ej ", priors$T1, " 1 2 --verbose | sed 1d'"),
                      IAMconstant = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20 0.0",
                                           " -n 1 ", priors$N1, 
                                           " -n 2 ", priors$N2, 
                                           " -em ", priors$anc_migT1, " 2 1 ", priors$m21,
                                           " -em ", priors$anc_migT1, " 1 2 ", priors$m12,              
                                           " -ej ", priors$T1, " 1 2 --verbose | sed 1d'"),
                      IAMgrowth = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20 0.0",
                                         " -g 1 ", priors$alpha1, 
                                         " -g 2 ", priors$alpha2, 
                                         " -eg ", priors$scT1, " 1 0.0",
                                         " -eg ", priors$scT2, " 2 0.0",
                                         " -em ", priors$anc_migT1, " 2 1 ", priors$m21,
                                         " -em ", priors$anc_migT1, " 1 2 ", priors$m12,              
                                         " -ej ", priors$T1, " 2 1 --verbose | sed 1d'"),
                      IAMbottleneck = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20 0.0",
                                             " -n 2 ", priors$N2, 
                                             " -g 1 ", priors$alpha1, 
                                             " -eg ", priors$scT1, " 1 0.0",
                                             " -em ", priors$anc_migT1, " 2 1 ", priors$m21,
                                             " -em ", priors$anc_migT1, " 1 2 ", priors$m12,              
                                             " -ej ", priors$T1, " 2 1 --verbose | sed 1d'"),
                      IRMconstant = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20",
                                           " -m 1 2 ", priors$m12, " -m 2 1 ", priors$m21,
                                           " -em ", priors$rec_migT1, " 2 1 ", 0,
                                           " -em ", priors$rec_migT1, " 1 2 ", 0,              
                                           " -ej ", priors$T1, " 2 1 --verbose | sed 1d'"),
                      IRMgrowth = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20",
                                         " -m 1 2 ", priors$m12, " -m 2 1 ", priors$m21,
                                         " -g 1 ", priors$alpha1, 
                                         " -g 2 ", priors$alpha2, 
                                         " -eg ", priors$scT2, " 1 0.0",
                                         " -eg ", priors$scT1, " 2 0.0",
                                         " -em ", priors$rec_migT1, " 2 1 ", 0,
                                         " -em ", priors$rec_migT1, " 1 2 ", 0,              
                                         " -ej ", priors$T1, " 2 1 --verbose | sed 1d'"),
                      IRMbottleneck = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20",
                                             " -m 1 2 ", priors$m12, " -m 2 1 ", priors$m21,
                                             " -n 2 ", priors$N2, 
                                             " -g 1 ", priors$alpha1, 
                                             " -eg ", priors$scT1, " 1 0.0 ",  
                                             " -em ", priors$rec_migT1, " 2 1 ", 0,
                                             " -em ", priors$rec_migT1, " 1 2 ", 0,              
                                             " -ej ", priors$T1, " 2 1 --verbose | sed 1d'"),
                      IARMconstant = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20",
                                            " -m 1 2 ", priors$m12, " -m 2 1 ", priors$m21,
                                            " -n 1 ", priors$N1, 
                                            " -n 2 ", priors$N2, 
                                            " -em ", priors$rec_migT1, " 2 1 ", 0,
                                            " -em ", priors$rec_migT1, " 1 2 ", 0,
                                            " -em ", priors$anc_migT1, " 2 1 ", priors$m21,
                                            " -em ", priors$anc_migT1, " 1 2 ", priors$m12, 
                                            " -ej ", priors$T1, " 2 1 --verbose | sed 1d'"),
                      IARMgrowth = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20",
                                          " -m 1 2 ", priors$m12, " -m 2 1 ", priors$m21,
                                          " -g 1 ", priors$alpha1, 
                                          " -g 2 ", priors$alpha2, 
                                          " -eg ", priors$scT1, " 1 0.0",
                                          " -eg ", priors$scT2, " 2 0.0",
                                          " -em ", priors$rec_migT1, " 2 1 ", 0,
                                          " -em ", priors$rec_migT1, " 1 2 ", 0,
                                          " -em ", priors$anc_migT1, " 2 1 ", priors$m21,
                                          " -em ", priors$anc_migT1, " 1 2 ", priors$m12, 
                                          " -ej ", priors$T1, " 2 1 --verbose | sed 1d'"),
                      IARMbottleneck = paste0("'msABC ", nind, " {} -t ", priors$theta, " -I 2 20 20",
                                              " -m 1 2 ", priors$m12, " -m 2 1 ", priors$m21,
                                              " -n 2 ", priors$N2, 
                                              " -g 1 ", priors$alpha1, 
                                              " -eg ", priors$scT1, " 1 0.0",
                                              " -em ", priors$rec_migT1, " 2 1 ", 0,
                                              " -em ", priors$rec_migT1, " 1 2 ", 0,
                                              " -em ", priors$anc_migT1, " 2 1 ", priors$m21,
                                              " -em ", priors$anc_migT1, " 1 2 ", priors$m12, 
                                              " -ej ", priors$T1, " 2 1 --verbose | sed 1d'")
)


# if heterozygous migration model, need to set migration rate per locus to be set
# from piped input file
if(mig_model == "hetero"){
  msline_base <- gsub("NA", "tbs", msline_base)
  msline_base <- sub("--verbose \\| sed 1d", "--verbose < loc_m_rate | sed 1d", msline_base)
}

# prep parallel commands
quart <- nloc*(1/nproc)
nloc_split <- rep(ceiling(quart), nproc)
nloc_split[length(nloc_split)] <- nloc_split[length(nloc_split)] - abs(nloc - sum(nloc_split))
if(sum(nloc_split) != nloc){
  print("Error in number of loci")
  q()
}
parallel_prefix <- "parallel"
if(mem_lim){
  parallel_prefix <- "parallel --memfree 500M"
}
parallel_suffix <- sprintf("::: %s", paste(nloc_split, collapse = " "))


## INITIALISE OUTPUT ###
# generate statnames (now internal)
if(mig_model == "hetero"){
  # create a file for locus specific rates in tbs
  m12 <- rbeta(nloc, a[1], b[1])*c1[1]
  m21 <- rbeta(nloc, a[1], b[1])*c2[1]
  y <- cbind(m12, m21)
  print(head(y))
  # write out heterozygous loci
  write.table(y, "./loc_m_rate", quote = F,
              row.names = F, col.names = F)
}
stat_base <- sub("\\{}", "1", msline_base[1])
stat_base <- sub("sed 1d", "sed -n 1p", stat_base)
stat_base <- gsub("\'", "", stat_base)
stat_name <- system(stat_base, intern = T)
stat_name <- unlist(strsplit(stat_name, "\t"))
# get param/stat indexes
param_ind <- grep("^p_", stat_name)
sust_ind <- grep("^s_", stat_name)
param_names <- stat_name[param_ind]
# add genome migration dist if hetero
if(mig_model == "hetero"){
  param_names <- c(stat_name[param_ind], c("p_a", "p_b", "p_c1", "p_c2"))
}
# get moments
mean <- paste0(stat_name[sust_ind], "_mean")
var <- paste0(stat_name[sust_ind], "_var")
skew <- paste0(stat_name[sust_ind], "_skew")
kurtosis <- paste0(stat_name[sust_ind], "_krt")
moment_names <- as.vector(rbind(mean, var, skew, kurtosis))
# recombine
stat_name <- c(param_names, moment_names)

### RUN MSABC ###
## open a file connection
myOutput <- file(outfile, open = "w")
# initialise output with statnames
write.table(t(stat_name), file = myOutput, append = FALSE, sep = " ", col.names = FALSE, quote = F,
            row.names = FALSE)

print(mig_model)
## run MSABC loop ##
sapply(1:nrep, function(z){
  if(mig_model == "hetero"){
    # create a file for locus specific rates in tbs
    m12 <- rbeta(nloc, a[z], b[z])*c1[z]
    m21 <- rbeta(nloc, a[z], b[z])*c2[z]
    y <- cbind(m12, m21)
    # write out heterozygous loci
    write.table(y, "./loc_m_rate", quote = F,
                row.names = F, col.names = F)
  }
  # construct msABC command
  msline <- paste(parallel_prefix, msline_base[z], parallel_suffix)
  # print line if verbose mode is on
  if(verbose){
    print(msline)
  }
  print(z)
  # run the msABC command - now reading data straight into R
  x <- suppressWarnings(system(msline, intern = T))
  x <- do.call(rbind, strsplit(x, "\t"))
  sim_data <- x[-1, ]
  storage.mode(sim_data) <- "numeric"
  # split into sum stats and param
  sust <- sim_data[, sust_ind]
  param <- sim_data[1 , param_ind]
  # add heterogenous migration params if hetero
  if(mig_model == "hetero"){
    param <- c(param, priors$a[z], priors$b[z], priors$c1[z], priors$c2[z])
  }
  #
  moments <- apply(sust, 2, function(z){
    c(mean(z, na.rm = T),
      var(z, na.rm = T),
      skew(z, na.rm = T),
      kurtosi(z, na.rm = T))
  })
  # recombine
  myOut <- c(param, c(moments))
  # write to connection
  write.table(t(myOut), file = myOutput, append = TRUE, sep = " ", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

})
# close connection
close(myOutput)