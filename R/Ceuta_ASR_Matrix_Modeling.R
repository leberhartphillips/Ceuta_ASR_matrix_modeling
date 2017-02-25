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
library(RMark) 
library(stringr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(reshape2)
library(RColorBrewer)
library(extrafont)
library(stats)
library(lme4)
###############################################################################
chick <- 
  read.table("/external/spazedawgs/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/data/chick_survival_data.txt",
             header = TRUE, colClasses = c("factor", "character","factor",
                                           "numeric","factor","factor", "numeric"))

fledgling_adult <- 
  read.table("/external/spazedawgs/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/data/fledgling_adult_survival_data.txt",
             header = TRUE, colClasses = c("factor", "character","factor","factor"))
###############################################################################
MarkPath <- "/usr/local/bin/mark"
MarkViewer <- "nano"
###############################################################################
# subset the captures so that only chicks captured on the day of hatch are included
# (the "ch" refers to the capture history of an individual on each day of its life as
# a chick.  Thus, if the first charactoer of the ch string is a 1, it was captured
# on the day of hatch and is included in the hatch sex ratio dataset)
caught_at_hatch <- chick[which(substring(chick$ch, 1, 1) == "1"),]

# sum the number of chicks that are included for each hatch ID
brood_ID_count <- 
  caught_at_hatch %>% count(brood_ID)

# join this data to the subset capture data
caught_at_hatch <- left_join(caught_at_hatch, brood_ID_count, by = "brood_ID")

# subset these data so that clutch size equals the number of chicks sampled from each nest
HSR_df <- caught_at_hatch[which(caught_at_hatch$clutch_size == caught_at_hatch$n),]

# make new columns "Male" and "Female" that have a 1 or a 0 to describe the sex of the chick
HSR_df$male <- ifelse(HSR_df$sex == "M", 1, 0)
HSR_df$female <- ifelse(HSR_df$sex == "F", 1, 0)

# define hatch ID as a factor
HSR_df$brood_ID <- as.factor(HSR_df$brood_ID)

# run mixed effects linear regression
# Brood ID is used as a random effect to control for the non-independence of siblings
HSR_model <- lme4::glmer(cbind(male, female) ~ (1| brood_ID),
                         data = HSR_df, family = binomial)

# check out the model results. P = 0.588, therefore hatching sex ratio doesn't
# deviate from parity
summary(HSR_model)

# calculate what the average hatching sex ratio is
# summarize the data so that each row is a nest instead of an individual
HSR_df_summary <- 
  HSR_df %>% 
  group_by(brood_ID) %>%
  summarise(no_males = sum(male),
            hatch_date_season = min(day_of_season),
            clutch_size = mean(n),
            year = first(year))

# calculate the proportion of the brood that was male
HSR_df_summary$prop_male <- HSR_df_summary$no_males/HSR_df_summary$clutch_size

# calculate the average hatching sex ratio across all nests
HSR <- mean(HSR_df_summary$prop_male)

# calculate the 95% confidence interval of the hatching sex ratio
HSR_95CI <- c(mean(HSR_df_summary$prop_male)-((sd(HSR_df_summary$prop_male)/sqrt(length(HSR_df_summary$prop_male)))*1.96),
              mean(HSR_df_summary$prop_male)+((sd(HSR_df_summary$prop_male)/sqrt(length(HSR_df_summary$prop_male)))*1.96))
###############################################################################
plover_matrix <- 
  function(demographic_rates)
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
###############################################################################
freq_dep_SSD_ASR <- 
  function (A, n = rep(10, nrow(A)), h = 1, k = 3, iterations = 30, HSR = 0.5) 
  {
    # Number of stages in matrix
    x <- length(n) 
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
    }
    
    # define rownames of stage matrix
    rownames(stage) <- rownames(A) 
    # define colnames of stage matrix
    colnames(stage) <- 0:(t - 1) 
    # define stable stage as the last stage
    w <- stage[, t] 
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
                                 size = nrow(fledgling_adult) - nrow(fledgling_adult_boot1), 
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
  # chick <- chick[!which(str_detect(chick[,"ch"],"1") == TRUE),]
  
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
    merge_design.covariates(fledgling_adult.ddl$Phi,
                            quad_time, bygroup = FALSE, bytime = TRUE)
  fledgling_adult.ddl$Phi <- 
    merge_design.covariates(fledgling_adult.ddl$Phi,
                            quad_time, bygroup = FALSE, bytime = TRUE)
  
  # ...and chicks
  time <- c(0:(chick.proc$nocc[1] - 1))
  quadratic <- time^2
  quad_time <- data.frame(time, quadratic)
  chick.ddl$p <- 
    merge_design.covariates(chick.ddl$Phi,
                            quad_time, bygroup = FALSE, bytime = TRUE)
  chick.ddl$Phi <- 
    merge_design.covariates(chick.ddl$Phi,
                            quad_time, bygroup = FALSE, bytime = TRUE)
  
  # Second step is to create the candidate survival models for 
  # fledglings and adults as a function
  fledgling_adult_survival = function() 
  {
    setwd("/external/spazedawgs/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/output/temp/") # set wd so that results go to the correct folder
    # sex- and stage-specific survival:
    Phi.agexsex = list(formula = ~ age * sex) 
    
    # Models exploring variation in encounter probability
    # sex-dependent:
    p.sex = list(formula =  ~ sex) 
    # age-dependent:
    p.age = list(formula =  ~ age) 
    # constant:
    p.dot = list(formula =  ~ 1) 
    # linear variation across time:
    p.time = list(formula =  ~ time) 
    # factorial variation across time:
    p.Time = list(formula =  ~ Time) 
    # quadratic variation across time:
    p.quadratic = list(formula =  ~ quadratic) 
    # interaction between sex and linear time:
    p.sexxtime = list(formula =  ~ sex * time) 
    # interaction between age and time:
    p.agextime = list(formula =  ~ age * time)
    # interaction between sex and factorial time:
    p.sexxTime = list(formula =  ~ sex * Time)
    # interaction between age and factorial time:
    p.agexTime = list(formula =  ~ age * Time)
    # interaction between age and sex:
    p.agexsex = list(formula =  ~ age * sex)
    # interaction between quadratic time and sex:
    p.quadraticxsex = list(formula =  ~ quadratic * sex)
    # interaction between quadratic time and age:
    p.quadraticxage = list(formula =  ~ quadratic * age)
    # # interaction between quadratic time and sex:
    # p.quadraticxagexsex = list(formula =  ~ quadratic * age * sex)
    # p.Timexagexsex = list(formula =  ~ Time * age * sex)
    # p.timexagexsex = list(formula =  ~ time * age * sex)
    # additive effects of sex and linear time:
    p.sex_time = list(formula =  ~ sex + time)
    # additive effects of age and linear time:
    p.age_time = list(formula =  ~ age + time)
    # additive effects of sex and factorial time:
    p.sex_Time = list(formula =  ~ sex + Time)
    # additive effects of age and factorial time:
    p.age_Time = list(formula =  ~ age + Time)
    # additive effects of age and sex:
    p.age_sex = list(formula =  ~ age + sex)
    # additive effects of sex and quadratic time:
    p.quadratic_sex = list(formula =  ~ quadratic + sex)
    # additive effects of age and quadratic time:
    p.quadratic_age = list(formula =  ~ quadratic + age)
    # additive effects of sex, age, and quadratic time:
    p.quadratic_age_sex = list(formula =  ~ quadratic + age + sex)
    # additive effects of sex, age, factorial time:
    p.Time_age_sex = list(formula =  ~ Time + age + sex)
    # additive effects of sex, age, linear time:
    p.time_age_sex = list(formula =  ~ time + age + sex)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml = RMark::create.model.list("CJS")
    
    # specify the data, design matrix, the number of threads to use 
    # (if using a server) and run the models in Program MARK
    model.list = RMark::mark.wrapper(cml, data = fledgling_adult.proc, 
                                     ddl = fledgling_adult.ddl,
                                     threads = 20, brief = TRUE, delete = TRUE)
    
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
    setwd("/external/spazedawgs/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/output/temp/") # set wd so that results go to the correct folder    
    # # sex- and linear age-specific survival:
    # Phi.Time.x.sex = list(formula = ~ Sex * Time)
    # sex- and quadratic age-specific survival:
    Phi.quadratic.x.sex = list(formula = ~ sex * quadratic)
    # # sex-specific survival:
    # Phi.Sex = list(formula = ~ Sex)
    
    # Models exploring variation in encounter probability
    # constant:
    p.dot = list(formula = ~ 1)
    # # linear across time
    # p.Time = list(formula = ~ Time)
    # quadratic across time
    p.quadratic = list(formula = ~ quadratic)
    # annual variation
    p.year = list(formula = ~ year)
    # sex-specific
    p.sex = list(formula = ~ sex)
    # # interaction between year and linear age
    # p.year.x.Time = list(formula = ~ Year * Time)
    # interaction between year and quadratic age
    p.year.x.quadratic = list(formula = ~ year * quadratic)
    # # interaction between sex and linear age
    # p.year.x.Time = list(formula = ~ Sex * Time)
    # interaction between year and quadratic age
    p.sex.x.quadratic = list(formula = ~ sex * quadratic)
    # # additive effects of sex and linear age
    # p.sex.Time = list(formula = ~ Sex + Time)
    # additive effects of sex and linear age
    p.sex.quadratic = list(formula = ~ sex + quadratic)
    # # additive effects of year and linear age    
    # p.year.Time = list(formula = ~ Year + Time)
    # additive effects of year and quadratic age    
    p.year.quadratic = list(formula = ~ year + quadratic)
    # # additive effects of year, sex, and linear age
    # p.year.Time.Sex = list(formula = ~ Year + Time + Sex)
    # additive effects of year, sex, and quadratic age
    p.year.quadratic.Sex = list(formula = ~ year + quadratic + sex)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml = RMark::create.model.list("CJS")
    
    # specify the data, design matrix, the number of threads to use 
    # (if using a server) and run the models in Program MARK
    model.list = RMark::mark.wrapper(cml, data = chick.proc, ddl = chick.ddl,
                                     threads  =  20, brief = TRUE, delete = TRUE)
    
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
                            # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                            h = 1,
                            k = 3,
                            # Define primary sex ratio (assumed to be 0.5)
                            HSR = HSR)
  
  # Build matrix based on rates specified in the list above
  demographic_matrix <- plover_matrix(demographic_rates)
  
  # Determine the ASR at the stable stage distribution
  ASR_SSD <- freq_dep_SSD_ASR(A = demographic_matrix, h = 1, k = 3, HSR = HSR)
  
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
niter <- 100

survival_ASR_bootstrap_result <- sapply(1:niter, run_bootstrap_survival_ASR, fledgling_adult, chick)
###############################################################################
# extract AIC_table_chick
AIC_table_chick_boot <- do.call(rbind, lapply(seq(from = 1, to = niter * 4, by = 4), function(x) survival_ASR_bootstrap_result[[x]]))
num_mods <- nrow(AIC_table_chick_boot)/niter
AIC_table_chick_boot$iter <- rep(1:niter, each = num_mods)

# extract AIC_table_fledgling_adult
AIC_table_fledgling_adult_boot <- do.call(rbind, lapply(seq(from = 2, to = niter * 4, by = 4), function(x) survival_ASR_bootstrap_result[[x]]))
num_mods <- nrow(AIC_table_fledgling_adult_boot)/niter
AIC_table_fledgling_adult_boot$iter <- rep(1:niter, each = num_mods)

# extract Survival_rates
survival_rates_boot <- do.call(rbind, lapply(seq(from = 3, to = niter * 4, by = 4), function(x) survival_ASR_bootstrap_result[[x]]))
survival_rates_boot$iter <- rep(1:niter, each = 6)

# extract ASR
ASR_boot <- sapply(seq(from = 4, to = niter * 4, by = 4), function(x) survival_ASR_bootstrap_result[[x]])
ASR_boot <- data.frame(ASR_boot = unname(ASR_boot), iter = 1:niter)
###############################################################################
sex_diff_survival <- function(survival_rates_boot) {
  # make an empty datarame to store the results
  sex_diff_surv_output <- data.frame(Adult = numeric(niter),
                                     Fledgling = numeric(niter),
                                     Chick = numeric(niter))
  # for loop to go through each iteration and calculate the differece between female and male
  # survival rates for each stage.
  for(i in 1:niter){
    Adult <- 
      survival_rates_boot[which(survival_rates_boot$iter == i), 2][2] -
      survival_rates_boot[which(survival_rates_boot$iter == i), 2][1]
    Fledgling <- 
      survival_rates_boot[which(survival_rates_boot$iter == i), 2][4] -
      survival_rates_boot[which(survival_rates_boot$iter == i), 2][3]
    Chick <- 
      survival_rates_boot[which(survival_rates_boot$iter == i), 2][6] -
      survival_rates_boot[which(survival_rates_boot$iter == i), 2][5]
    
    sex_diff_surv_output[i, 1] <- Adult
    sex_diff_surv_output[i, 2] <- Fledgling
    sex_diff_surv_output[i, 3] <- Chick
  }
  # restructure the output and lable columns
  sex_diff_surv_output <- reshape2::melt(data = sex_diff_surv_output)
  colnames(sex_diff_surv_output) <- c("stage", "difference")
  # return the output
  sex_diff_surv_output
}

# run the function on the bootstrap list from above
sex_diff_survival_output <- sex_diff_survival(survival_rates_boot)

# Calculate some summary statistics
sex_diff_survival_summary <- 
  sex_diff_survival_output %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(avg = mean(Difference),
                   median = median(Difference),
                   var = var(Difference))

# specify custom color palette to distingush first-year stages (i.e. chicks and
# fledglings) from adults
cbPalette <- c("#D9D9D9", "#A6A6A6", "#A6A6A6")

# reorder the levels of the stage factors
sex_diff_survival_output$stage <- 
  factor(sex_diff_survival_output$Stage, levels = c("Adult", "Fledgling", "Chick"))

# plot the sex-biases in survival across the three stages
Sex_bias_violin <- 
ggplot2::ggplot(aes(y = difference, x = stage, fill = stage), data = sex_diff_survival_output) + 
  coord_flip() +
  theme_bw() +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  annotate("text", x = 2, y = -0.25,
           label = c("\u2640"), size = 14,
           vjust = c(0.5), hjust = c(0.5)) +
  annotate("text", x = 2, y = 0.25,
           label = c("\u2642"), size = 14,
           vjust = c(0.5), hjust = c(0.5)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.x = element_text(size=12, margin = margin(10, 0, 0, 0)),
        axis.text.x  = element_text(size=10, margin = margin(5, 0, 0, 0)), 
        axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
        axis.text.y  = element_text(size=10, angle = 90, hjust = 0.5, margin = margin(0, 5, 0, 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        panel.border = element_rect(linetype = "solid", colour = "grey"),
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  scale_fill_manual(values = cbPalette) +
  scale_y_continuous(limits=c(-0.25, 0.25)) +
  xlab("Life-stage") + 
  ylab("Sex-bias in apparent survival (\u03D5)")
###############################################################################
# specify the confidence interval bounds
CI <- 0.95

# Calculate the confidence interval, mean, and median of the ASR bootstraps
ASR_boot_95CI <- stats::quantile(ASR_boot$ASR_boot, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
ASR_boot_mean <- mean(ASR_boot$ASR_boot)
ASR_boot_median <- median(ASR_boot$ASR_boot)

# consolidate the results
ASR_boot_summary <- as.data.frame(cbind(ASR_boot_95CI[1], ASR_boot_95CI[2], ASR_boot_mean, ASR_boot_median))
rownames(ASR_boot_summary) <- NULL
colnames(ASR_boot_summary) <- c("lcl", "ucl", "mean", "median")

# plot the ASR histogram
ASR_histogram <- 
ggplot2::ggplot() +
          annotate("text", x = c(-Inf,Inf), y = c(70, 70),
                   label = c("\u2640", "\u2642"), size = 14,
                   vjust = c(0.5,0.5), hjust = c(-0.5,1.5)) +
          geom_histogram(binwidth = 0.02, data = ASR_boot, aes(x = ASR_boot)) +
          geom_errorbarh(data = ASR_boot_summary, aes(y = 125, x = lcl, xmin = lcl, xmax = ucl), color = "black", size = 0.8, linetype = "solid") +
          theme_bw() +
          theme(legend.position="none",
                legend.position = c(0, 1), 
                legend.justification = c(0, 1),
                legend.text=element_text(size=11),
                legend.title=element_blank(),
                legend.key.height=unit(0.8,"line"),
                legend.key.width=unit(0.8,"line"),
                legend.background = element_rect(fill=NA),
                axis.title.x = element_text(size=12, margin = margin(10, 0, 0, 0)),
                axis.text.x  = element_text(size=10, margin = margin(5, 0, 0, 0)), 
                axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
                axis.text.y  = element_text(size=10, angle = 90, hjust = 0.5, margin = margin(0, 5, 0, 0), color = "white"),
                axis.ticks.y = element_line(size = 0.5, colour = "white"),
                axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
                axis.ticks.length = unit(0.2, "cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(linetype = "solid", colour = "grey"),
                plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
                strip.background = element_blank(), 
                strip.text = element_blank(),
                panel.margin = unit(0.75, "lines")) +
          ylab("Frequency") +
          xlab("Adult sex ratio (proportion \u2642)") +
          scale_x_continuous(limits = c(0.0, 1)) +
          scale_y_continuous(limits = c(0, 130))