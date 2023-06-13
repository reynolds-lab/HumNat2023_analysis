### Multilogit Model ###
#Austin W Reynolds
#April.7.2021

### User-defined variables ###
outputpath <- "~/Dropbox/Cederberg_project/Datasets/"
datapath <- "~/Dropbox/Cederberg_project/Datasets/df.locality.csv"
modelpath <- "~/Dropbox/Cederberg_project/STAN_models/expanded_multilogit_apr14.stan"

### Load required modules ###
library(rstan) 
library(loo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#functions
data.assigner <- function(data,design.matrix){
  # makes data list input for STAN model
  datalist <- list(
    K = 4, # number of categories of the response
    N = dim(data)[1], # sample size for the model
    y = data$locality.numeric, # integer coding for locality pattern used in the STAN model
    D = ncol(design.matrix), # Number of co-variates in the model
    x = design.matrix, # design matrix
    mu_raw = rep(0,3*dim(design.matrix)[2]), # prior mean vector for coefficients of locality model
    N_vmp = max(data$m_bp.numeric), # Number of unique mother's birthplaces in the dataset
    N_vfp = max(data$f_bp.numeric), # Number of unique father's birthplaces in the dataset
    N_id = max(data$subjectID.numeric), # Number of probands in the dataset
    id = data$subjectID.numeric, # integer vector of proband IDs for indexing proband random effects
    vmp = data$m_bp.numeric, # integer vector of mother's birthplace IDs for indexing birthplace effects
    vfp = data$f_bp.numeric # integer vector of father's birthplace IDs for indexing birthplace effects
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
df.locality <- read.csv(datapath)


### Run STAN models ###
# Base model
x.frame <- model.frame(ones ~ generation * Region + subject_birthyear_centered_centennial, data=df.locality)
base.matrix <- model.matrix(ones ~ generation * Region + subject_birthyear_centered_centennial, data=x.frame)
base.datalist <- data.assigner(data = df.locality, design.matrix = base.matrix)
base.model <- model.runner(base.datalist,modelpath,filename = "harmonization.age.test.multilogit.rds", iter = 500, warmup = 250, chains = 2)
# Model summary
summary(base.model, pars=c("beta_raw", "sigma_id", "sigma_vmp", "sigma_vfp"))


# Base Model + Size Category
x.frame <- model.frame(ones ~ generation * Region + subject_birthyear_centered_centennial + size_categorical, data=df.locality)
size.matrix <- model.matrix(ones ~ generation * Region + subject_birthyear_centered_centennial + size_categorical, data=x.frame)
size.datalist <- data.assigner(data = df.locality, design.matrix = size.matrix)
size.model <- model.runner(size.datalist,modelpath,filename = "harmonization.age.test.multilogit.size.rds", iter = 500, warmup = 250, chains = 2)
# Model summary
summary(size.model, pars=c("beta_raw", "sigma_id", "sigma_vmp", "sigma_vfp"))

# Base Model + Ancestry
x.frame <- model.frame(ones ~ generation * Region + subject_birthyear_centered_centennial + logratio_KHS_EUR + logratio_BTU_EUR, data=df.locality)
ancestry.matrix <- model.matrix(ones ~ generation * Region + subject_birthyear_centered_centennial + logratio_KHS_EUR + logratio_BTU_EUR, data=x.frame)
ancestry.datalist <- data.assigner(data = df.locality, design.matrix = ancestry.matrix)
ancestry.model <- model.runner(ancestry.datalist,modelpath,filename = "harmonization.age.test.multilogit.ancestry.rds", iter = 500, warmup = 250, chains = 2)
# Model summary
summary(ancestry.model, pars=c("beta_raw", "sigma_id", "sigma_vmp", "sigma_vfp"))


### Perform approximate leave-one-out cross-validation using Pareto smoothed importance sampling
# for model comparison ###
base.loo <- loo(x = base.model, cores = 4)
size.loo <- loo(x = size.model, cores = 4)
ancestry.loo <- loo(x = ancestry.model, cores = 4)

# Model comparison using PSIS-LOO cross-validations as information criteria
loo_compare(base.loo, size.loo, ancestry.loo)
