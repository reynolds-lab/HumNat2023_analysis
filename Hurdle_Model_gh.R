### Hurdle Model ###
#Austin W Reynolds
#April.7.2021

### User-defined variables ###
outputpath <- "~/Dropbox/Cederberg_project/Datasets/"
datapath <- "~/Dropbox/Cederberg_project/Datasets/df.migration.csv"
modelpath <- "~/Dropbox/Cederberg_project/STAN_models/basemodel_loglik_6MAR2020.stan"

### Load required modules ###
library(rstan) 
library(loo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

### Define functions ###
data.assigner <- function(data,design.matrix){
  # makes data list input for STAN model
  datalist <- list(
    # For all cases
    N = dim(data)[1], # Sample size of the input dataset
    N_subjects = length(unique(data$subjectID)), # Number of probands in the input dataset
    N_birthplaces = length(unique(data$birthplace)), # Number of birthplaces in the input dataset
    migrated = data$migrated, # Binary indicating whether individual migrated (1) or not (0)
    X = design.matrix, # design matrix
    K = dim(design.matrix)[2], # Number of co-variates in the model
    mu_mig = rep(0,dim(design.matrix)[2]), # Prior mean vector for coefficients of migration model
    mu_dis = c(log(100), rep(0,dim(design.matrix)[2]-1)), # Prior mean vector for coefficients of distance model
    subID = data$subjectID.numeric, # Vector of indices for proband (family) random effects
    bp = data$birthplace.numeric, # Vector of indices for birthplaces random effects
    
    # For complete cases (case == 1, 2)  
    N_complete = sum(data$case == 1 | data$case == 2), # Count of complete cases in the input dataset
    complete_cases = which(data$case == 1 | data$case == 2), # Row indices for complete cases in the input dataset
    subID_complete_cases = data$subjectID.numeric[which(data$case == 1 | data$case == 2)], #SubjectID indices for complete cases in the input dataset
    bp_complete_cases = data$birthplace.numeric[which(data$case == 1 | data$case == 2)], #Birthplace indices for complete cases in the input dataset
    nonzero = ifelse(data$distance[which(data$case == 1 | data$case == 2)] > 0, data$distance[which(data$case == 1 | data$case == 2)], 1), # Indicator for the hurdle model for complete cases in the input dataset
    
    # For left-censored cases (case = 3) 
    N_lftcensored = sum(data$case == 3), # Count of left-censored cases in the input dataset
    lftcensored_cases = which(data$case == 3), # Row indices for left-censored cases in the input dataset
    subID_lftcensored_cases = data$subjectID.numeric[which(data$case == 3)], #SubjectID indices for left-censored cases in the input dataset
    bp_lftcensored_cases = data$birthplace.numeric[which(data$case == 3)],#Birthplace indices for left-censored cases in the input dataset
    maxdist_lftcensored = data$maxdist[which(data$case == 3)], # Upper bound on an uncertain migration distance in left-censored cases in the input dataset
    
    # For interval-censored cases (case = 4) 
    N_intcensored = sum(data$case == 4), # Count of interval-censored cases in the input dataset
    intcensored_cases = which(data$case == 4), # Row indices for interval-censored cases in the input dataset
    subID_intcensored_cases = data$subjectID.numeric[which(data$case == 4)], #SubjectID indices for interval-censored cases in the input dataset
    bp_intcensored_cases = data$birthplace.numeric[which(data$case == 4)],#Birthplace indices for interval-censored cases in the input dataset
    maxdist_intcensored = data$maxdist[which(data$case == 4)], # Upper bound on an uncertain migration distance in interval-censored cases in the input dataset
    mindist_intcensored = data$mindist[which(data$case == 4)] # Lower bound on an uncertain migration distance in interval-censored cases in the input dataset
  ) 
  return(datalist)
}

model.runner <- function(datalist, modelpath, filename, iter = 2000, warmup = 1000, chains = 4){
  modelout <- stan(file = modelpath, 
                   data = datalist,
                   iter = iter, warmup = warmup, chains = chains, control = list(adapt_delta = 0.99))
  
  saveRDS(modelout,paste(outputpath,filename,sep = ""))
  return(modelout)
}


### Read in dataset ###
df.migration <- read.csv(datapath)


### Run STAN models ###
# Base Model
x.frame <- model.frame(ones ~  sex * generation * Region + subject_birthyear_centered_centennial, data=df.migration)
base.matrix <- model.matrix(ones ~ sex * generation * Region + subject_birthyear_centered_centennial, data=x.frame)
base.datalist <- data.assigner(data = df.migration, design.matrix = base.matrix)
base.model <- model.runner(base.datalist,modelpath,filename = "harmonization.age.test.rds", iter = 500, warmup = 250, chains = 2)
# Model summary 
summary(base.model,
        pars=c(
          "beta_dis","beta_mig",
          "sd_lnorm", "sd_u_dis", "sd_v_dis", "sd_u_mig", "sd_v_mig"))

# Base Model + Size Category
x.frame <- model.frame(ones ~  sex * generation * Region + subject_birthyear_centered_centennial + size_categorical, data=df.migration)
size.matrix <- model.matrix(ones ~ sex * generation * Region + subject_birthyear_centered_centennial + size_categorical, data=x.frame)
size.datalist <- data.assigner(data = df.migration, design.matrix = size.matrix)
size.model <- model.runner(size.datalist,modelpath,filename = "harmonization.age.test.size.rds", iter = 500, warmup = 250, chains = 2)
# Model summary 
summary(size.model,
        pars=c(
          "beta_dis","beta_mig",
          "sd_lnorm", "sd_u_dis", "sd_v_dis", "sd_u_mig", "sd_v_mig"))

# Base Model + Ancestry
x.frame <- model.frame(ones ~  sex * generation * Region + subject_birthyear_centered_centennial + logratio_KHS_EUR + logratio_BTU_EUR, data=df.migration)
ancestry.matrix <- model.matrix(ones ~ sex * generation * Region + subject_birthyear_centered_centennial + logratio_KHS_EUR + logratio_BTU_EUR, data=x.frame)
ancestry.datalist <- data.assigner(data = df.migration, design.matrix = ancestry.matrix)
ancestry.model <- model.runner(ancestry.datalist, modelpath, filename = "harmonization.age.test.ancestry.rds", iter = 500, warmup = 250, chains = 2)
# Model summary 
summary(ancestry.model,
        pars=c(
          "beta_dis","beta_mig",
          "sd_lnorm", "sd_u_dis", "sd_v_dis", "sd_u_mig", "sd_v_mig"))


### Perform approximate leave-one-out cross-validation using Pareto smoothed importance sampling
# for model comparison ###
base.loo <- loo(x = base.model, cores = 4)
size.loo <- loo(x = size.model, cores = 4)
ancestry.loo <- loo(x = ancestry.model, cores = 4)

# Model comparison using PSIS-LOO cross-validations as information criteria
loo_compare(base.loo, size.loo, ancestry.loo)





