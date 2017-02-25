# install.packages("RMark")
# install.packages("stringr")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("grid")
# install.packages("gridExtra")
# install.packages("reshape2")
# install.packages("RColorBrewer")
# install.packages("extrafont")
# install.packages("stats")
# install.packages("lme4")
# install.packages("parallel")
#library(snowfall)
library(RMark) 
library(parallel) 
library(stringr)
#library(ggplot2)
#library(dplyr)
#library(gridExtra)
#library(grid)
#library(reshape2)
#library(RColorBrewer)
#library(extrafont)
#library(stats)
#library(lme4)
###############################################################################
chick <- 
  read.table("data/chick_survival_data.txt",
             header = TRUE, colClasses = c("factor", "character","factor",
                                           "numeric","factor","factor", "numeric"))

fledgling_adult <- 
  read.table("data/fledgling_adult_survival_data.txt",
             header = TRUE, colClasses = c("factor", "character","factor","factor"))
###############################################################################
MarkPath <- "/usr/local/bin/mark"
MarkViewer <- "nano"
###############################################################################
# Hatching sex ratio (proportion of hatchlings that are male)
HSR <- 0.4856322 
# Sex-specific annual per captia fecundity rates (hatchlings)
# Females:
RF <- 2.06#1.967778
# Males:
RM <- 0#1.847123
###############################################################################
plover_matrix <- 
  function(demographic_rates, mating_function = TRUE)
  {
    if(mating_function)
    {
      # Define plover life-stages of the Ceuta snowy plover matrix model
      stages <- c("F_1st_yr",  "F_Adt",  "M_1st_yr",  "M_Adt")
      # Build the 4x4 matrix
      result <- 
        matrix(c(0, NA, 0, NA, 
                 (demographic_rates$F_Chk_survl*demographic_rates$F_Fdg_survl),
                 demographic_rates$F_Adt_survl, 
                 0, 0,
                 0, NA, 0, NA,
                 0, 0, 
                 (demographic_rates$M_Chk_survl*demographic_rates$M_Fdg_survl),
                 demographic_rates$M_Adt_survl),
               nrow = 4, byrow = TRUE,
               dimnames = list(stages, stages))
      result
    }
    else
    {
      # Define plover life-stages of the Ceuta snowy plover matrix model
      stages <- c("F_1st_yr",  "F_Adt",  "M_1st_yr",  "M_Adt")
      # Build the 4x4 matrix
      result <- 
        matrix(c(0, demographic_rates$RF * (1 - HSR), 0, demographic_rates$RM * (1 - HSR), 
                 (demographic_rates$F_Chk_survl*demographic_rates$F_Fdg_survl),
                 demographic_rates$F_Adt_survl, 
                 0, 0,
                 0, demographic_rates$RF * HSR, 0, demographic_rates$RM * HSR,
                 0, 0, 
                 (demographic_rates$M_Chk_survl*demographic_rates$M_Fdg_survl),
                 demographic_rates$M_Adt_survl),
               nrow = 4, byrow = TRUE,
               dimnames = list(stages, stages))
      result
    }
  }
###############################################################################
freq_dep_SSD_ASR <- 
  function (A, n = rep(10, nrow(A)), h = 1, k = 3, 
            iterations = 30, HSR = 0.5, mating_function = TRUE) 
  {
    # Number of stages in matrix
    x <- length(n) 
    if(mating_function)
    {
      # Number of time steps to simulate
      t <- iterations 
      # an empty t by x matrix
      stage <- matrix(numeric(x * t), nrow = x) 
      # for loop that goes through each of t time steps
      for (i in 1:t) { 
        # stage distribution at time t
        stage[,i] <- n 
        # number of male adults at time t
        M2 <- stage[4, i] 
        # number of female adults at time t
        F2 <- stage[2, i] 
        # Female freq-dep fecundity of Female chicks
        A[1,x/2]        <- (k*M2)/(M2+(F2/h))*HSR 
        # Female freq-dep fecundity of Male chicks
        A[(x/4)*3,x/2]  <- (k*M2)/(M2+(F2/h))*HSR
        # Male freq-dep fecundity of Female chicks
        A[1,x]          <- (k*F2)/(M2+(F2/h))*HSR 
        # Male freq-dep fecundity of Male chicks
        A[(x/4)*3,x]    <- (k*F2)/(M2+(F2/h))*HSR 
        # define the new n (i.e., new stage distribution at time t)
        n <- A %*% n 
        # define rownames of stage matrix
        rownames(stage) <- rownames(A) 
        # define colnames of stage matrix
        colnames(stage) <- 0:(t - 1) 
        # define stable stage as the last stage
        w <- stage[, t] 
      }
      # calculate the proportional stable stage distribution
      stable.stage <- w/sum(w)
      # calc ASR as the proportion of the adult stable stage class that is male
      ASR <- stable.stage[x]/(stable.stage[x/2] + stable.stage[x])
      
      # make a list of results
      pop.proj <- list(ASR = ASR,
                       stable.stage = stable.stage, 
                       stage.vectors = stage,
                       SSD_M2 = stable.stage[4],
                       SSD_F2 = stable.stage[2])
    }
    else
    {
      ev <- eigen(A)
      lmax <- which.max(Re(ev$values))
      W <- ev$vectors
      w <- abs(Re(W[, lmax]))
      names(w) <- colnames(A)
      # calculate the proportional stable stage distribution
      stable.stage <- w/sum(w)
      # calc ASR as the proportion of the adult stable stage class that is male
      ASR <- stable.stage[x]/(stable.stage[x/2] + stable.stage[x])
      
      # make a list of results
      pop.proj <- list(ASR = ASR,
                       stable.stage = stable.stage, 
                       SSD_M2 = stable.stage[4],
                       SSD_F2 = stable.stage[2])
    }
    # print the list as output to the function
    pop.proj 
  }
###############################################################################
bootstrap_data <- function(fledgling_adult, chick) {
  chick_boot <- chick[sample(1:nrow(chick), 
                             size = nrow(chick), 
                             replace = TRUE), ]
  present <- fledgling_adult$ring %in% chick_boot$ring
  fledgling_adult_boot1 <- fledgling_adult[present, ]
  spare_fledgling_adult <- fledgling_adult[!present, ]
  fledgling_adult_boot2 <- 
    spare_fledgling_adult[sample(1:nrow(spare_fledgling_adult), 
                                 size = nrow(fledgling_adult) - 
                                   nrow(fledgling_adult_boot1), 
                                 replace = TRUE), ]
  fledgling_adult_boot <- rbind(fledgling_adult_boot1, fledgling_adult_boot2)
  out <- list(chick_boot = chick_boot, fledgling_adult_boot = fledgling_adult_boot)
}
###############################################################################
bootstrap_survival_ASR <- function(plover_boot_list, num_boot) {
  
  # First steps require some wrangling
  
  # specify the bootstrapped data samples
  chick <- plover_boot_list[["chick_boot"]]
  fledgling_adult <- plover_boot_list[["fledgling_adult_boot"]]
  
  # remove ring column
  fledgling_adult <- fledgling_adult[,-1]
  chick <- chick[,-1]
  
  # Remove capture histories that have no resights (i.e., all zeros in the ch)
  nrow(chick[!which(str_detect(chick[,"ch"],"1") == TRUE),])
  nrow(fledgling_adult[!which(str_detect(fledgling_adult[,"ch"],"1") == TRUE),])
  
  # Create processed RMARK data format as Cormack-Jolly_Seber with 2 groups 
  # (sex and age initally ringed), starting at year 2006, two age groups
  # (first-years and adults) in which the first-year stage only lasts for 
  # one year.
  fledgling_adult.proc <- RMark::process.data(fledgling_adult, model = "CJS",
                                              groups = c("sex", "age"),
                                              begin.time = 2006, age.var = 2, 
                                              initial.age = c(1, 0))
  
  # Create processed RMARK data format as Cormack-Jolly_Seber with 3 groups 
  # (sex, year, and brood ID).
  chick.proc <-  RMark::process.data(chick, model = "CJS",
                                     groups = c("sex", "year", "brood_ID"))
  
  # Create the design matrix from the processed mark-recapture datasets
  fledgling_adult.ddl <- RMark::make.design.data(fledgling_adult.proc)
  chick.ddl <- RMark::make.design.data(chick.proc)
  
  # adds first-year / adult age field to design data in column "Age"
  fledgling_adult.ddl <- RMark::add.design.data(data = fledgling_adult.proc,
                                                ddl = fledgling_adult.ddl, 
                                                parameter = "Phi", 
                                                type = "age",
                                                bins = c(0, 1, 7), right = FALSE,
                                                name = "age", replace = TRUE)
  
  # create a dummy field in the design matrix called marked.as.adult 
  # which is "0" for the group initally ringed as chicks and "1" for the group
  # marked as adults.
  fledgling_adult.ddl$Phi$marked.as.adult = 0
  fledgling_adult.ddl$Phi$marked.as.adult[fledgling_adult.ddl$Phi$initial.age.class == "A"] = 1 
  fledgling_adult.ddl$p$marked.as.adult = 0
  fledgling_adult.ddl$p$marked.as.adult[fledgling_adult.ddl$p$initial.age.class == "A"] = 1
  
  # # check parameter matrices to see if groups were binned correctly
  # PIMS(mark(fledgling_adult.proc,fledgling_adult.ddl,model.parameters=list(Phi=list(formula=~age+sex)),output=F),"Phi")
  # 
  # Create quadratic time variable so that it can be
  # tested for temporal variation in fledgling and adult survival...
  time <- c(0:(fledgling_adult.proc$nocc[1] - 1))
  quadratic <- time^2
  quad_time <- data.frame(time, quadratic)
  quad_time$time <- c(2006:2012)
  fledgling_adult.ddl$p <- 
    RMark::merge_design.covariates(fledgling_adult.ddl$Phi,
                                   quad_time, bygroup = FALSE, bytime = TRUE)
  fledgling_adult.ddl$Phi <- 
    RMark::merge_design.covariates(fledgling_adult.ddl$Phi,
                                   quad_time, bygroup = FALSE, bytime = TRUE)
  
  # ...and chicks
  time <- c(0:(chick.proc$nocc[1] - 1))
  quadratic <- time^2
  quad_time <- data.frame(time, quadratic)
  chick.ddl$p <- 
    RMark::merge_design.covariates(chick.ddl$Phi,
                                   quad_time, bygroup = FALSE, bytime = TRUE)
  chick.ddl$Phi <- 
    RMark::merge_design.covariates(chick.ddl$Phi,
                                   quad_time, bygroup = FALSE, bytime = TRUE)
  
  # Second step is to create the candidate survival models for 
  # fledglings and adults as a function
  fledgling_adult_survival = function() 
  {
    # sex- and stage-specific survival:
    Phi.agexsex = list(formula = ~ age * sex) 
    
    # Models exploring variation in encounter probability
    # constant:
    p.dot = list(formula =  ~ 1)
    # sex-dependent:
    p.sex = list(formula =  ~ sex)
    # age-dependent:
    p.age = list(formula =  ~ age)
    # factorial variation across year:
    p.Time = list(formula =  ~ Time)
    # interaction between sex and factorial year:
    p.sexxTime = list(formula =  ~ sex * Time)
    # interaction between age and factorial year:
    p.agexTime = list(formula =  ~ age * Time)
    # interaction between age and sex:
    p.agexsex = list(formula =  ~ age * sex)
    # additive effects of sex and factorial year:
    p.sex_Time = list(formula =  ~ sex + Time)
    # additive effects of age and factorial year:
    p.age_Time = list(formula =  ~ age + Time)
    # additive effects of age and sex:
    p.age_sex = list(formula =  ~ age + sex)
    #additive effects of sex, age, factorial year:
    p.Time_age_sex = list(formula =  ~ Time + age + sex)
    #additive effect of year and interaction between age and sex:
    p.Time_age_x_sex = list(formula =  ~ Time + age * sex)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml <-  RMark::create.model.list("CJS")
    
    # specify the data, design matrix, the number of threads to use 
    # (if using a server) and run the models in Program MARK
    model.list <-  RMark::mark.wrapper(cml, data = fledgling_adult.proc,
                                       threads = 8, delete = TRUE)
    
    # output the model list and sotre the results
    return(model.list)
  }
  
  # Run the models on the bootstrapped data
  fledgling_adult_survival_run <- 
    fledgling_adult_survival()
  
  # Extract the AIC model table from the model output
  AIC_table_fledgling_adult <- 
    fledgling_adult_survival_run$model.table
  
  # Find the model number for the first ranked model of the AIC table
  model_fledgling_adult_num <- 
    as.numeric(rownames(fledgling_adult_survival_run$model.table[1,]))
  
  # extract and format survival rates from fledgling and adult model output
  fledgling_adult_reals <- 
    fledgling_adult_survival_run[[model_fledgling_adult_num]]$results$real
  
  # format the output to tidy up the sex- and age-specific effects
  Groups <- data.frame(str_split_fixed(rownames(fledgling_adult_reals), " ", n = 5))
  fledgling_adult_reals <- cbind(Groups, fledgling_adult_reals)
  fledgling_adult_reals <- 
    fledgling_adult_reals[which(fledgling_adult_reals$X1 == "Phi"),]
  fledgling_adult_reals$age <- 
    unlist(str_extract_all(fledgling_adult_reals$X2,"[AJ]"))
  fledgling_adult_reals$age <- 
    as.factor(ifelse(fledgling_adult_reals$age == "A","Adult","Fledgling"))
  fledgling_adult_reals$sex <- 
    unlist(str_extract_all(fledgling_adult_reals$X2,"[FM]"))
  fledgling_adult_reals$sex <- 
    as.factor(ifelse(fledgling_adult_reals$sex == "F","Female","Male"))
  fledgling_adult_reals$sex_age <- 
    paste(fledgling_adult_reals$sex,fledgling_adult_reals$age,sep = "_")
  fledgling_adult_survival_real <- 
    fledgling_adult_reals[,c("sex_age", "estimate")]
  row.names(fledgling_adult_survival_real) <- NULL
  
  # Do the same for chicks by creating the candidate survival models for 
  # chicks as a function
  chick_survival = function() 
  {
    # sex- and quadratic age-specific survival:
    Phi.quadratic.x.sex_C = list(formula = ~ sex * quadratic)
    
    # Models exploring variation in encounter probability
    # constant:
    p.dot = list(formula = ~ 1)
    # quadratic across age
    p.quadratic = list(formula = ~ quadratic)
    # annual variation
    p.year = list(formula = ~ year)
    # sex-specific
    p.sex = list(formula = ~ sex)
    # interaction between year and quadratic age
    p.year.x.quadratic = list(formula = ~ year * quadratic)
    # interaction between year and quadratic age
    p.sex.x.quadratic = list(formula = ~ sex * quadratic)
    # additive effects of sex and linear age
    p.sex.quadratic = list(formula = ~ sex + quadratic)
    # additive effects of year and quadratic age
    p.year.quadratic = list(formula = ~ year + quadratic)
    # additive effects of year, sex, and quadratic age
    p.year.quadratic.Sex = list(formula = ~ year + quadratic + sex)
    # additive effect of year and interaction between sex and quadratic age
    p.year.quadratic.x.Sex = list(formula = ~ year + quadratic * sex)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml <-  RMark::create.model.list("CJS")
    
    # specify the data, design matrix, the number of threads to use 
    # (if using a server) and run the models in Program MARK
    model.list <-  RMark::mark.wrapper(cml, data = chick.proc, ddl = chick.ddl,
                                       threads = 8, delete = TRUE)
    
    # output the model list and sotre the results
    return(model.list)
  }
  
  # Run the models on the bootstrapped data
  chick_survival_run <- chick_survival()
  
  # Extract the AIC model table from the model output
  AIC_table_chick <- chick_survival_run$model.table
  
  # Find the model number for the first ranked model of the AIC table
  model_chick_num <- as.numeric(rownames(chick_survival_run$model.table[1,]))
  
  # extract real parameter estimates from top models
  chick_reals <- chick_survival_run[[model_chick_num]]$results$real
  
  # format the output to tidy up the sex- and age-specific effects
  Groups <- data.frame(str_split_fixed(rownames(chick_reals), " ", n = 5))
  chick_reals <- cbind(Groups, chick_reals)
  chick_reals <- chick_reals[which(chick_reals$X1 == "Phi"),]
  chick_reals$sex <- unlist(str_extract_all(chick_reals$X2,"[FM]"))
  chick_reals$sex <- as.factor(ifelse(chick_reals$sex == "F","Female","Male"))
  
  # transform the daily chick survival (DCS) to apparent hatching success.
  # this can either be done by DCS^25 if there is no age effect:
  if(nrow(chick_reals) == 2)
  {
    plover_Survival_to_Fledge_F <- 
      chick_reals[which(chick_reals$sex == "Female"),
                  c("estimate")]^25
    plover_Survival_to_Fledge_M <- 
      chick_reals[which(chick_reals$sex == "Male"),
                  c("estimate")]^25
  }
  # or by calculating the product of all DCS estimates:
  if(nrow(chick_reals) != 2){
    plover_Survival_to_Fledge_F <- 
      prod(chick_reals[which(chick_reals$sex == "Female"),
                       c("estimate")][c(1:26)])
    plover_Survival_to_Fledge_M <- 
      prod(chick_reals[which(chick_reals$sex == "Male"),
                       c("estimate")][c(1:26)])
  }
  
  # tidy up the output and put it in a dataframe.
  estimate <- c(plover_Survival_to_Fledge_F, plover_Survival_to_Fledge_M)
  sex <- c("Female", "Male")
  age <- c("Chick", "Chick")
  sex_age <- paste(sex, age, sep = "_")
  chick_survival_real <- data.frame(sex_age, estimate)
  
  # Bind the fledgling and adult dataframe with the chicks
  survival_rates <- rbind(fledgling_adult_survival_real, chick_survival_real)
  
  # Create a list of demographic rates from the survival analyses above
  demographic_rates <- list(F_Chk_survl = survival_rates[5,2],
                            F_Fdg_survl = survival_rates[3,2],
                            F_Adt_survl = survival_rates[1,2],
                            M_Chk_survl = survival_rates[6,2],
                            M_Fdg_survl = survival_rates[4,2],
                            M_Adt_survl = survival_rates[2,2],
                            # # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                            # h = 1,
                            # k = 3,
                            # Define hatching sex ratio
                            HSR = HSR,
                            # Define the fecundity of males (RM) and females (RF)
                            RF = RF,
                            RM = RM)
  
  # Build matrix based on rates specified in the list above
  demographic_matrix <- plover_matrix(demographic_rates, mating_function = FALSE)
  
  # Determine the ASR at the stable stage distribution
  ASR_SSD <- freq_dep_SSD_ASR(A = demographic_matrix, mating_function = FALSE)
  
  # Extract ASR
  ASR_estimate <- ASR_SSD$ASR
  
  # make a list of all the results from this iteration
  bootstrap_results_list <- 
    list(AIC_table_chick, 
         AIC_table_fledgling_adult, 
         survival_rates,
         ASR_estimate)
}
###############################################################################
run_bootstrap_survival_ASR <- function(num_boot, fledgling_adult, chick)
{
  bootstrap_data_list <- bootstrap_data(fledgling_adult, chick)
  result <- bootstrap_survival_ASR(bootstrap_data_list, num_boot)
}
###############################################################################
niter <- 1
# ncores <- 1
# e1 <- environment()
# cl <- parallel::makeCluster(ncores)
# parallel::clusterExport(cl, c("run_bootstrap_survival_ASR", "bootstrap_data", "bootstrap_survival_ASR", "plover_matrix", "freq_dep_SSD_ASR",
#                               "str_detect", "process.data", "make.design.data",
#                               "add.design.data",
#                               "merge_design.covariates",
#                               "create.model.list",
#                               "mark.wrapper",
#                               "str_split_fixed",
#                               "str_extract_all",
#                               "HSR",
#                               "RF",
#                               "RM"))
# 
# survival_ASR_bootstrap_result <- parallel::parSapply(cl, 1:niter, run_bootstrap_survival_ASR, fledgling_adult, chick)
# parallel::stopCluster(cl)

system.time(survival_ASR_bootstrap_result <- sapply(1:niter, run_bootstrap_survival_ASR, fledgling_adult, chick))

###############################################################################
# extract AIC_table_chick
AIC_table_chick_boot <- do.call(rbind, lapply(seq(from = 1, to = niter * 4, by = 4), function(x) survival_ASR_bootstrap_result[[x]]))
num_mods <- nrow(AIC_table_chick_boot)/niter
AIC_table_chick_boot$iter <- rep(1:niter, each = num_mods)
write.table(AIC_table_chick_boot, file = "output/bootstrap/AIC_table_chick_boot_out.txt")

# extract AIC_table_fledgling_adult
AIC_table_fledgling_adult_boot <- do.call(rbind, lapply(seq(from = 2, to = niter * 4, by = 4), function(x) survival_ASR_bootstrap_result[[x]]))
num_mods <- nrow(AIC_table_fledgling_adult_boot)/niter
AIC_table_fledgling_adult_boot$iter <- rep(1:niter, each = num_mods)
write.table(AIC_table_fledgling_adult_boot, file = "output/bootstrap/AIC_table_fledgling_adult_boot_out.txt")

# extract Survival_rates
survival_rates_boot <- do.call(rbind, lapply(seq(from = 3, to = niter * 4, by = 4), function(x) survival_ASR_bootstrap_result[[x]]))
survival_rates_boot$iter <- rep(1:niter, each = 6)
write.table(survival_rates_boot, file = "output/bootstrap/survival_rates_boot_out.txt")

# extract ASR
ASR_boot <- sapply(seq(from = 4, to = niter * 4, by = 4), function(x) survival_ASR_bootstrap_result[[x]])
ASR_boot <- data.frame(ASR_boot = unname(ASR_boot), iter = 1:niter)
write.table(ASR_boot, file = "output/bootstrap/ASR_boot_out.txt")