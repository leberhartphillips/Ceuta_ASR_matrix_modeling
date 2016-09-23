# Script for running analyses AFTER the bootstrap survival analysis

# Author: Luke J. Eberhart-Phillips
# Date: September 14, 2016

# The script takes the results from the bootstrapped survival analysis and 
# runs/creates the following output:
# 1) Summary stats/plot of sex-difference in survival
# 2) Summary stats/plot of ASR for the population
# 3) Life-table response experiment

# POINTS TO DISCUSS:
# 1) Use of the non-linear mating function. I've added an argument to all 
#    functions that allow for the option of using either the frequency-dependent
#    mating function or not (i.e., with a constant fecundity)
# 2) When using the constant fecundity, I'm not sure if I should use just
#    female fecundity, or both male and female fecundity.  Based on our
#    previous discussions, I think you recommended that females are the only
#    sex to produce offspring.  I've set up the script so that male fecundity
#    can easily be set to 0 or to the estimate calculated from our raw field
#    data.

############################ LOAD LIBRARIES ##################################
library(ggplot2)
library(dplyr)
library(reshape2)
library(Rmisc)
library(stringr)
library(RColorBrewer)
library(grid)

############################ DATA IMPORT ######################################
survival_rates_boot <- 
  read.table("data/bootstrap_surv_results.txt",
             header = TRUE, colClasses = c("factor", "numeric","factor"))

################### CONSTANTS TO SET PRIOR TO ANALYSIS ########################
# NOTE: all the following constant parameters were calculated from raw field data.
# These calculations will be included in the script elsewhere.

# number of iterations used in bootstrap
niter <- 1000 

# Hatching sex ratio (proportion of hatchlings that are male)
HSR <- 0.4856322 

# Sex-specific annual per captia fecundity rates (hatchlings)
# Females:
RF <- 1.967045

# Males:
# NOTE: I think that this should be set as "0", but I may be wrong. 
# According to our data, males produced on average "1.839908" hatchlings
# annually per capita. Should this go in for male fecundity?...or should
# it be left as "0"?
#RM <- 0
RM <- 1.839908
############################ FUNCTIONS ########################################
# The following section outlines all the functions used in the analysis.
# NOTE: the "mating_function" TRUE/FALSE argument specifies if the harmonic
# mating function is used (i.e., "TRUE", which is default), or not (i.e.,
# "FALSE")
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

freq_dep_SSD_ASR <- 
  function (A, n = rep(10, nrow(A)), h = 1, k = 3, iterations = 30, HSR = 0.5, mating_function = TRUE) 
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

lower_level_sens_analysis <- 
  function(freq_dep_ASR, VR_list, mating_function = TRUE)
  {
    if(mating_function)
    {
      # make a list of all parameters
      vr <- list(
        F_Chk_survl = VR_list$F_Chk_survl,
        F_Fdg_survl = VR_list$F_Fdg_survl,
        F_Adt_survl = VR_list$F_Adt_survl,
        M_Chk_survl = VR_list$M_Chk_survl,
        M_Fdg_survl = VR_list$M_Fdg_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        h = VR_list$h,
        k = VR_list$k,
        HSR = VR_list$HSR,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
      # make a matrix of the elements
      el <- expression(0, ((k * M2) / (M2 + (F2 / h))) * (1 - HSR), 
                       0, ((k * F2) / (M2 + (F2 / h))) * (1 - HSR),
                       (F_Chk_survl * F_Fdg_survl), F_Adt_survl, 0, 0,
                       0, ((k * M2) / (M2 + (F2 / h))) * HSR, 0, 
                       ((k * F2) / (M2 + (F2 / h))) * HSR,
                       0, 0, (M_Chk_survl * M_Fdg_survl), M_Adt_survl)
      # calculate the effect of proportional changes in vital rates
      n <- length(vr)
      vr_nums <- seq(0, 2, 0.01) # proportional changes in increments of 0.01 between 0 and 2
      # First calculate sensitivites
      vrsen <- matrix(numeric(n * length(vr_nums)), 
                      ncol = n, dimnames = list(vr_nums, names(vr)))
      for (h in 1:n)
      {
        vr2 <- vr
        for (i in 1:length(vr_nums))
        {
          vr2[[h]] <- vr_nums[i]
          A <- 
            matrix(sapply(el, eval, vr2, NULL), nrow = sqrt(length(el)), byrow=TRUE)
          vrsen[i, h] <- 
            eigen(A)$vectors[4, 1] / (eigen(A)$vectors[2, 1] + eigen(A)$vectors[4, 1])
        }
      }
      # Next calculate rescaled elasticities
      vrelas <- matrix(numeric(n * length(vr_nums)), 
                       ncol = n, dimnames = list(vr_nums, names(vr)))
      for (h in 1:n)
      {
        for (i in 1:length(vr_nums))
        {
          vr2 <- vr
          vr2[[h]] <- vr_nums[i] * vr2[[h]]
          A <- matrix(sapply(el, eval, vr2 , NULL), nrow = sqrt(length(el)), byrow = TRUE)
          vrelas[i, h] <- 
            (eigen(A)$vectors[4, 1] / 
               (eigen(A)$vectors[2, 1] + eigen(A)$vectors[4, 1])) / 
            unname(freq_dep_ASR$ASR)
        }
      }
      # tidy up and label results
      colnames(vrsen) <- c("Female chick survival", 
                           "Female fledgling survival", 
                           "Female adult survival", 
                           "Male chick survival", 
                           "Male fledgling survival", 
                           "Male adult survival",
                           "Mating system index (h)", 
                           "Clutch size",
                           "Hatching sex ratio",
                           "Breeding males",
                           "Breeding females")
      colnames(vrelas) <- c("Female chick survival", 
                            "Female fledgling survival", 
                            "Female adult survival", 
                            "Male chick survival", 
                            "Male fledgling survival", 
                            "Male adult survival",
                            "Mating system (h)", 
                            "Clutch size",
                            "Hatching sex ratio",
                            "Breeding males",
                            "Breeding females")
    }
    else
    {
      # make a list of all parameters
      vr <- list(
        F_Chk_survl = VR_list$F_Chk_survl,
        F_Fdg_survl = VR_list$F_Fdg_survl,
        F_Adt_survl = VR_list$F_Adt_survl,
        M_Chk_survl = VR_list$M_Chk_survl,
        M_Fdg_survl = VR_list$M_Fdg_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        RF = VR_list$RF,
        RM = VR_list$RM,
        HSR = VR_list$HSR)
      # make a matrix of the elements
      el <- expression(0, RF * (1 - HSR), 
                       0, RM * (1 - HSR),
                       (F_Chk_survl * F_Fdg_survl), F_Adt_survl, 0, 0,
                       0, RF * HSR,
                       0, RM * HSR,
                       0, 0, (M_Chk_survl * M_Fdg_survl), M_Adt_survl)
      # calculate the effect of proportional changes in vital rates
      n <- length(vr)
      vr_nums <- seq(0, 2, 0.01) # proportional changes in increments of 0.01 between 0 and 2
      # First calculate sensitivites
      vrsen <- matrix(numeric(n * length(vr_nums)), 
                      ncol = n, dimnames = list(vr_nums, names(vr)))
      for (h in 1:n)
      {
        vr2 <- vr
        for (i in 1:length(vr_nums))
        {
          vr2[[h]] <- vr_nums[i]
          A <- 
            matrix(sapply(el, eval, vr2, NULL), nrow = sqrt(length(el)), byrow=TRUE)
          vrsen[i, h] <- 
            eigen(A)$vectors[4, 1] / (eigen(A)$vectors[2, 1] + eigen(A)$vectors[4, 1])
        }
      }
      # Next calculate rescaled elasticities
      vrelas <- matrix(numeric(n * length(vr_nums)), 
                       ncol = n, dimnames = list(vr_nums, names(vr)))
      for (h in 1:n)
      {
        for (i in 1:length(vr_nums))
        {
          vr2 <- vr
          vr2[[h]] <- vr_nums[i] * vr2[[h]]
          A <- matrix(sapply(el, eval, vr2 , NULL), nrow = sqrt(length(el)), byrow = TRUE)
          vrelas[i, h] <- 
            (eigen(A)$vectors[4, 1] / 
               (eigen(A)$vectors[2, 1] + eigen(A)$vectors[4, 1])) / 
            unname(freq_dep_ASR$ASR)
        }
      }
      # tidy up and label results
      colnames(vrsen) <- c("Female chick survival", 
                           "Female fledgling survival", 
                           "Female adult survival", 
                           "Male chick survival", 
                           "Male fledgling survival", 
                           "Male adult survival",
                           "Female Fecundity",
                           "Male Fecundity",
                           "Hatching sex ratio")
      colnames(vrelas) <- c("Female chick survival", 
                            "Female fledgling survival", 
                            "Female adult survival", 
                            "Male chick survival", 
                            "Male fledgling survival", 
                            "Male adult survival",
                            "Female Fecundity",
                            "Male Fecundity",
                            "Hatching sex ratio")
    }
    Sensitivities <- melt(vrsen)
    colnames(Sensitivities) <- c("Perturbation", "Vitalrate", "Sensitivity")
    Elasticities <- melt(vrelas)
    colnames(Elasticities) <- c("Perturbation", "Vitalrate", "Elasticity")
    results <- list(Sensitivities = Sensitivities,
                    Elasticities = Elasticities,
                    Element_expression = el)
    results
  }

ASR_analysis <- 
  function (A, zero = TRUE) 
  {
    ev <- eigen(A) # makes list of the eigen values and eigen vectors of A
    lmax <- which.max(Re(ev$values)) # index of dominant eigen value
    W <- ev$vectors # Eigen vectors
    w <- abs(Re(W[, lmax])) # dominant eigen vector
    stable.stage = w / sum(w) # stable stage distribution
    ASR <- stable.stage[4] / (stable.stage[2] + stable.stage[4]) # SSD ASR
    V <- try(Conj(solve(W)), silent = TRUE) # check if possible to proceed
    if (class(V) == "try-error") {
      ASR.analysis <- list(ASR = ASR, stable.stage = stable.stage, 
                           sensitivities = A * NA, elasticities = A * NA)
    }
    else {
      v <- abs(Re(V[lmax, ])) # solve matrix
      s <- v %o% w # outer product of v and w
      if (zero) {
        s[A == 0] <- 0
      }
      e <- s * A/ASR # calculate elasticities
      x <- dimnames(A) # get vital rate names
      dimnames(s) <- x # assign vital rate names to s
      names(w) <- x[[1]]
      names(v) <- x[[1]]
      ASR.analysis <- list(ASR = ASR, stable.stage = stable.stage, 
                           sensitivities = s, elasticities = e)
    }
    ASR.analysis
  }

# LTRE analysis of ASR
vitalsens_ASR <- 
  function (elements, VR_list, freq_dep_ASR, mating_function = TRUE) 
  {
    if(mating_function)
    {
      # list of parameters in the treatment matrix
      # this contains the observed paramters based on previous analyses
      treatment_matrix <- list(
        F_Chk_survl = VR_list$F_Chk_survl,
        F_Fdg_survl = VR_list$F_Fdg_survl,
        F_Adt_survl = VR_list$F_Adt_survl,
        M_Chk_survl = VR_list$M_Chk_survl,
        M_Fdg_survl = VR_list$M_Fdg_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        h = VR_list$h,
        k = VR_list$k,
        HSR = VR_list$HSR,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
      # list of parameters in the control matrix
      # this contains parameters with no sex-differences
      control_matrix <- list(
        F_Chk_survl = VR_list$M_Chk_survl,
        F_Fdg_survl = VR_list$M_Fdg_survl,
        F_Adt_survl = VR_list$M_Adt_survl,
        M_Chk_survl = VR_list$M_Chk_survl,
        M_Fdg_survl = VR_list$M_Fdg_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        h = VR_list$h,
        k = VR_list$k,
        HSR = 0.5,
        M2 = unname(freq_dep_ASR$SSD_M2),
        F2 = unname(freq_dep_ASR$SSD_F2))
      # check if everything is correctly structured before proceeding
      if (is.vector(treatment_matrix)) {
        treatment_matrix <- as.list(treatment_matrix)
      }
      if (!is.list(treatment_matrix)) {
        stop("Vital rates should be a vector or list")
      }
      if (class(elements) != "expression") {
        stop("Matrix elements should be an expression")
      }
      if (is.vector(control_matrix)) {
        control_matrix <- as.list(control_matrix)
      }
      if (!is.list(control_matrix)) {
        stop("Vital rates should be a vector or list")
      }
      # find the number of stage and sex specific parameters
      n <- sqrt(length(elements))
      if (n%%1 != 0) {
        stop(paste("Length of element expression is", length(elements), 
                   "- Expecting power of 2 like 4, 9, 16 to form a square matrix"))
      }
      # add the mating function parameters to the matrices
      vrs <- try(sapply(elements, eval, treatment_matrix, NULL), silent = TRUE)
      vrs_LTRE <- try(sapply(elements, eval, control_matrix, NULL), silent = TRUE)
      if (class(vrs) == "try-error") {
        vrs <- sub("Error in eval\\(expr, envir, enclos\\) :",
                   "", vrs[1])
        stop(paste("Cannot evaluate element expression using given vital rates:",
                   vrs))
      }
      if (class(vrs_LTRE) == "try-error") {
        vrs_LTRE <- sub("Error in eval\\(expr, envir, enclos\\) :",
                        "", vrs_LTRE[1])
        stop(paste("Cannot evaluate element expression using given vital rates:",
                   vrs_LTRE))
      }
      # make an empty dataframe where all the perturbation stats will go
      res <- data.frame(estimate = unlist(treatment_matrix), sensitivity = 0, 
                        elasticity = 0, LTRE = 0)
      # build the treatment matrix
      A <- matrix(vrs, nrow = n, byrow = TRUE)
      # build the control matrix
      A_LTRE <- matrix(vrs_LTRE, nrow = n, byrow = TRUE)
      # transform the matrix to M-prime (see formula in manuscript)
      Ac <- (A + A_LTRE) / 2
      # run sensitivity analyses on both matrices
      SAc <- ASR_analysis(Ac)
      ASR <- ASR_analysis(A)
      # calculate derivatives of lower-level matrix elements
      deriv.funcs <- sapply(elements, deriv, namevec = names(treatment_matrix), 
                            function.arg = TRUE)
      devs <- lapply(deriv.funcs, function(x) do.call(x, treatment_matrix))
      # run for loop to go through each parameter and estimate elasticity,
      # sensitivity, and LTRE
      for (i in 1:length(treatment_matrix)) {
        derivs <- matrix(as.numeric(lapply(devs, function(x) attr(x, "gradient")[i])), 
                         nrow = n, byrow = TRUE)
        res[i, 2] <- sum(derivs * ASR$sensitivities)
        res[i, 3] <- treatment_matrix[[i]] / ASR$ASR * sum(derivs * ASR$sensitivities)
        # only do LTRE calculations on survival parameters RELATIVE to one sex
        # i.e., don't calculate LTRE on mating system components
        res[i, 4] <- ifelse(i > 3 & i < 6, NA,
                            ifelse(i < 4, (treatment_matrix[[i + 3]] - treatment_matrix[[i]]) * 
                                     sum(derivs * SAc$sensitivities),
                                   ifelse(i == 9, (control_matrix[[i]] - treatment_matrix[[i]]) * 
                                            sum(derivs * SAc$sensitivities),
                                          NA)))
      }
      # consolidate results
      y <- res
      y$Vital_rate <- as.factor(rownames(y))
      colnames(y) <- c("Estimate", "Sensitivity", "Elasticity", "LTRE", "Vital_rate")
      y_melt <- suppressMessages(melt(y[,c(2:5)]))
      y_melt$parameter <- 
        as.factor(ifelse(str_detect(y_melt$Vital_rate,"Chk"), "Chick survival",
                         ifelse(str_detect(y_melt$Vital_rate,"Fdg"), "Fledgling survival",
                                ifelse(str_detect(y_melt$Vital_rate,"Adt"), "Adult survival",
                                       ifelse(str_detect(y_melt$Vital_rate,"HSR"), "Hatching sex ratio", 
                                              ifelse(str_detect(y_melt$Vital_rate,"h"), "Mating System",
                                                     ifelse(str_detect(y_melt$Vital_rate,"k"), "Clutch size", 
                                                            "No. breeding adults")))))))
      y_melt$parameter <- factor(y_melt$parameter, levels = c("Adult survival",
                                                              "Fledgling survival",
                                                              "Chick survival",
                                                              " ",
                                                              "No. breeding adults",
                                                              "Mating System",
                                                              "Clutch size",
                                                              "Hatching sex ratio"))
    }
    else
    {
      # list of parameters in the treatment matrix
      # this contains the observed paramters based on previous analyses
      treatment_matrix <- list(
        F_Chk_survl = VR_list$F_Chk_survl,
        F_Fdg_survl = VR_list$F_Fdg_survl,
        F_Adt_survl = VR_list$F_Adt_survl,
        M_Chk_survl = VR_list$M_Chk_survl,
        M_Fdg_survl = VR_list$M_Fdg_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        RF = VR_list$RF,
        RM = VR_list$RM,
        HSR = VR_list$HSR)
      # list of parameters in the control matrix
      # this contains parameters with no sex-differences
      control_matrix <- list(
        F_Chk_survl = VR_list$M_Chk_survl,
        F_Fdg_survl = VR_list$M_Fdg_survl,
        F_Adt_survl = VR_list$M_Adt_survl,
        M_Chk_survl = VR_list$M_Chk_survl,
        M_Fdg_survl = VR_list$M_Fdg_survl,
        M_Adt_survl = VR_list$M_Adt_survl,
        RF = VR_list$RF,
        RM = VR_list$RM,
        HSR = VR_list$HSR)
      # check if everything is correctly structured before proceeding
      if (is.vector(treatment_matrix)) {
        treatment_matrix <- as.list(treatment_matrix)
      }
      if (!is.list(treatment_matrix)) {
        stop("Vital rates should be a vector or list")
      }
      if (class(elements) != "expression") {
        stop("Matrix elements should be an expression")
      }
      if (is.vector(control_matrix)) {
        control_matrix <- as.list(control_matrix)
      }
      if (!is.list(control_matrix)) {
        stop("Vital rates should be a vector or list")
      }
      # find the number of stage and sex specific parameters
      n <- sqrt(length(elements))
      if (n%%1 != 0) {
        stop(paste("Length of element expression is", length(elements), 
                   "- Expecting power of 2 like 4, 9, 16 to form a square matrix"))
      }
      # add the mating function parameters to the matrices
      vrs <- try(sapply(elements, eval, treatment_matrix, NULL), silent = TRUE)
      vrs_LTRE <- try(sapply(elements, eval, control_matrix, NULL), silent = TRUE)
      if (class(vrs) == "try-error") {
        vrs <- sub("Error in eval\\(expr, envir, enclos\\) :",
                   "", vrs[1])
        stop(paste("Cannot evaluate element expression using given vital rates:",
                   vrs))
      }
      if (class(vrs_LTRE) == "try-error") {
        vrs_LTRE <- sub("Error in eval\\(expr, envir, enclos\\) :",
                        "", vrs_LTRE[1])
        stop(paste("Cannot evaluate element expression using given vital rates:",
                   vrs_LTRE))
      }
      # make an empty dataframe where all the perturbation stats will go
      res <- data.frame(estimate = unlist(treatment_matrix), sensitivity = 0, 
                        elasticity = 0, LTRE = 0)
      # build the treatment matrix
      A <- matrix(vrs, nrow = n, byrow = TRUE)
      # build the control matrix
      A_LTRE <- matrix(vrs_LTRE, nrow = n, byrow = TRUE)
      # transform the matrix to M-prime (see formula in manuscript)
      Ac <- (A + A_LTRE) / 2
      # run sensitivity analyses on both matrices
      SAc <- ASR_analysis(Ac)
      ASR <- ASR_analysis(A)
      # calculate derivatives of lower-level matrix elements
      deriv.funcs <- sapply(elements, deriv, namevec = names(treatment_matrix), 
                            function.arg = TRUE)
      devs <- lapply(deriv.funcs, function(x) do.call(x, treatment_matrix))
      # run for loop to go through each parameter and estimate elasticity,
      # sensitivity, and LTRE
      for (i in 1:length(treatment_matrix)) {
        derivs <- matrix(as.numeric(lapply(devs, function(x) attr(x, "gradient")[i])), 
                         nrow = n, byrow = TRUE)
        res[i, 2] <- sum(derivs * ASR$sensitivities)
        res[i, 3] <- treatment_matrix[[i]] / ASR$ASR * sum(derivs * ASR$sensitivities)
        # only do LTRE calculations on survival parameters RELATIVE to one sex
        # i.e., don't calculate LTRE on mating system components
        res[i, 4] <- ifelse(i > 3 & i < 6, NA,
                            ifelse(i < 4, (treatment_matrix[[i + 3]] - treatment_matrix[[i]]) * 
                                     sum(derivs * SAc$sensitivities),
                                   ifelse(i == 9, (control_matrix[[i]] - treatment_matrix[[i]]) * 
                                            sum(derivs * SAc$sensitivities),
                                          NA)))
      }
      # consolidate results
      y <- res
      y$Vital_rate <- as.factor(rownames(y))
      colnames(y) <- c("Estimate", "Sensitivity", "Elasticity", "LTRE", "Vital_rate")
      y_melt <- suppressMessages(melt(y[,c(2:5)]))
      y_melt$parameter <- 
        as.factor(ifelse(str_detect(y_melt$Vital_rate,"Chk"), "Chick survival",
                         ifelse(str_detect(y_melt$Vital_rate,"Fdg"), "Fledgling survival",
                                ifelse(str_detect(y_melt$Vital_rate,"Adt"), "Adult survival",
                                       ifelse(str_detect(y_melt$Vital_rate,"HSR"), "Hatching sex ratio",
                                              ifelse(str_detect(y_melt$Vital_rate,"RF"), "Female Fecundity", "Male Fecundity"))))))
      y_melt$parameter <- factor(y_melt$parameter, levels = c("Adult survival",
                                                              "Fledgling survival",
                                                              "Chick survival",
                                                              "Female Fecundity",
                                                              "Male Fecundity",
                                                              "Hatching sex ratio"))
    }
    y_melt$Sex <- as.factor(ifelse(str_detect(y_melt$Vital_rate,"F_") & 
                                     y_melt$variable != "LTRE", "Female", 
                                   ifelse(str_detect(y_melt$Vital_rate,"M_") & 
                                            y_melt$variable != "LTRE","Male", "Other")))
    y_melt$Sex <- 
      factor(y_melt$Sex,
             levels = c("Female","Male", "Other"))
    y_melt$value_trans <- ifelse(y_melt$Sex == "Female", 
                                 abs(y_melt$value)*-1, y_melt$value)
    y_melt <- y_melt[,-c(1)]
    results <- list(Sensitivity = subset(y_melt, (variable == "Sensitivity")),
                    Elasticity = subset(y_melt, (variable == "Elasticity")),
                    LTRE = subset(y_melt, (variable == "LTRE" & !is.na(value))))
    results$LTRE$parameter <- 
      factor(results$LTRE$parameter,
             levels = c("Adult survival",
                        "Fledgling survival",
                        "Chick survival",
                        "Hatching sex ratio"))
    row.names(results$Sensitivity) <- NULL
    row.names(results$Elasticity) <- NULL
    row.names(results$LTRE) <- NULL
    results$Sensitivity$value_trans <- 
      as.numeric(results$Sensitivity$value_trans)
    results$Elasticity$value_trans <- 
      as.numeric(results$Elasticity$value_trans)
    results$LTRE$value_trans <- 
      as.numeric(results$LTRE$value_trans)
    results
  }

ASR_extract <- 
  function(survival_rates) 
  {
    # make an empty datarame to store the results
    ASR_output <- data.frame(R = numeric(niter),
                             h2= numeric(niter),
                             h1= numeric(niter),
                             h0.5= numeric(niter),
                             h0.3= numeric(niter))
    # for loop to go through each iteration and calculate the differece between female and male
    # survival rates for each stage.
    for(i in 1:niter){
      # Create a list of demographic rates from the survival analyses above
      demographic_rates <- list(F_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][5],
                                F_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][3],
                                F_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][1],
                                M_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][6],
                                M_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][4],
                                M_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][2],
                                # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                                h = 1,
                                k = 3,
                                # Define primary sex ratio (assumed to be 0.5)
                                HSR = HSR,
                                # Define female and male fecundity
                                RF = RF,
                                RM = RM)
      
      # Build matrix based on rates specified in the list above
      demographic_matrix_R <- plover_matrix(demographic_rates, mating_function = FALSE)
      demographic_matrix_fun <- plover_matrix(demographic_rates, mating_function = TRUE)
      
      # Determine the ASR at the stable stage distribution
      ASR_SSD_R <- freq_dep_SSD_ASR(A = demographic_matrix_R, HSR = HSR, mating_function = FALSE)
      ASR_SSD_h2 <- freq_dep_SSD_ASR(A = demographic_matrix_fun, h = 2, HSR = HSR, mating_function = TRUE)
      ASR_SSD_h1 <- freq_dep_SSD_ASR(A = demographic_matrix_fun, h = 1, HSR = HSR, mating_function = TRUE)
      ASR_SSD_h0.5 <- freq_dep_SSD_ASR(A = demographic_matrix_fun, h = 0.5, HSR = HSR, mating_function = TRUE)
      ASR_SSD_h0.3 <- freq_dep_SSD_ASR(A = demographic_matrix_fun, h = 0.3, HSR = HSR, mating_function = TRUE)
      
      # Extract ASR
      ASR_output[i, 1] <- ASR_SSD_R$ASR
      ASR_output[i, 2] <- ASR_SSD_h2$ASR
      ASR_output[i, 3] <- ASR_SSD_h1$ASR
      ASR_output[i, 4] <- ASR_SSD_h0.5$ASR
      ASR_output[i, 5] <- ASR_SSD_h0.3$ASR
      
    }
    # restructure the output and lable columns
    ASR_output <- suppressMessages(reshape2::melt(data = ASR_output))
    colnames(ASR_output) <- c("ASR_model", "estimate")
    # return the output
    ASR_output
  }

################################ SURVIVAL RATES ###############################
# summarise the sex differences in survival across the three stages for the 
# bootstrap results
sex_diff_survival <- function(survival_rates_boot) {
  # make an empty datarame to store the results
  sex_diff_surv_output <- data.frame(Adult = numeric(niter),
                                     Fledgling = numeric(niter),
                                     Chick = numeric(niter))
  # for loop to go through each iteration and calculate the differece between 
  # female and male survival rates for each stage.
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
  sex_diff_surv_output <- suppressMessages(reshape2::melt(data = sex_diff_surv_output))
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
  dplyr::summarise(avg = mean(difference),
                   median = median(difference),
                   var = var(difference))

# specify custom color palette to distingush first-year stages (i.e. chicks and
# fledglings) from adults
cbPalette <- c("#A6A6A6", "#D9D9D9", "#D9D9D9")
cbPalette <- c("#D9D9D9", "#A6A6A6", "#A6A6A6")

# reorder the levels of the stage factors
sex_diff_survival_output$stage <- 
  factor(sex_diff_survival_output$stage, levels = c("Adult", "Fledgling", "Chick"))

# plot the sex-biases in survival across the three stages
Background <- 
  ggplot(aes(y = difference, x = stage, fill = stage), data = sex_diff_survival_output) + 
  coord_flip() +
  theme_bw() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(2)]) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(1)]) +
  annotate("text", x = c(2,2), y = c(-Inf, Inf),
           label = c("\u2640", "\u2642"), size = 7,
           family="Arial", vjust = c(0.5,0.5), hjust = c(-0.3,1.3)) +
  theme(text = element_text(family="Arial", color = "white"), # set the font as Candara
        legend.position = "none",
        axis.title.x = element_text(size=12, margin = margin(10, 0, 0, 0)),
        axis.text.x  = element_text(size=10, margin = margin(5, 0, 0, 0)), 
        axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
        axis.text.y  = element_text(size=10, angle = 90, hjust = 0.5, margin = margin(0, 5, 0, 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.5, colour = "white"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_blank(),
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  scale_x_continuous(limits=c(0,4),breaks=c(0,1,2), labels=c("Chick", "Fledgling", "Adult")) +
  scale_y_continuous(limits=c(-0.25,0.25)) +
  xlab("Life-stage") + 
  ylab("Sex-bias in apparent survival (\u03D5)")
Background

Bootstrap_sex_diff_VR_plot <- 
  ggplot(aes(y = difference, x = stage, fill = stage), data = sex_diff_survival_output) + 
  coord_flip() +
  theme_bw() +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(width = 0.2) +
  theme(text = element_text(family="Arial"),
        legend.position = "none",
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
  scale_y_continuous(limits=c(-0.25,0.25)) +
  xlab("Life-stage") + 
  ylab("Sex-bias in apparent survival (\u03D5)")
Bootstrap_sex_diff_VR_plot

jpeg(filename = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/figs/Surv_sex_differences_final_arial_grey_switch.jpeg",
     quality = 100,
     width = 4,
     height = 4,
     units = "in",
     res = 300)

grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) )
print( Background + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( Bootstrap_sex_diff_VR_plot + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
dev.off()

####################################### ASR ###################################
# extract the ASR estimate for each iteration of the bootstrapped survival
# analysis
ASR_boot <- ASR_extract(survival_rates_boot)

# Calculate the confidence interval, mean, and median of the ASR bootstraps
ASR_boot_summary <- 
  Rmisc::summarySE(filter(ASR_boot, ASR_model == "R"), 
                   measurevar = "estimate", 
                   groupvars = c("ASR_model"))

CI <- 0.95
Ceuta_ASR_95CI_quan <- stats::quantile(filter(ASR_boot, ASR_model == "R")[2], c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
ASR_boot_summary <- as.data.frame(cbind(ASR_boot_summary, Ceuta_ASR_95CI_quan[1], Ceuta_ASR_95CI_quan[2]))
rownames(ASR_boot_summary) <- NULL
colnames(ASR_boot_summary) <- c("ASR_model", "N", "estimate", "sd", "se", "ci", "lcl", "ucl")

# plot the histogram distributions of the two models (i.e., "...R" vs. "...fun")
# Note: "ASR_R" is the matrix model without the non-linear mating function, and
# "ASR_fun" is the matrix model with the non-linear mating function.

ASR_bootstrap_histogram <- 
  ggplot() +
  annotate("rect", xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(2)]) +
  annotate("rect", xmin=0.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(1)]) +
  annotate("text", x = c(-Inf,Inf), y = c(95, 95),
           label = c("\u2640", "\u2642"), size = 7,
           family="Arial", vjust = c(1.5,1.5), hjust = c(-0.5,1.5)) +
  geom_histogram(binwidth = 0.02, data = filter(ASR_boot, ASR_model == "R"), aes(x = estimate)) +
  geom_errorbarh(data = ASR_boot_summary, aes(y = 155, x = lcl, xmin = lcl, xmax = ucl), color = "black", size = 0.8, linetype = "solid") +
  theme_bw() +
  theme(text=element_text(family="Arial"),
        legend.position="none",
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
  scale_y_continuous(limits = c(0, 160))
ASR_bootstrap_histogram

ggsave(ASR_bootstrap_histogram,
       filename = "ASR_distributions_final_arial.jpg",
       path = "figs/",
       width = 4,
       height = 3, units = "in",
       dpi = 300,
       scale = 1)

############################## LTRE ###########################################
# Summarise the bootstrap stage- and sex-specific survival rates for the
# deterministic matrix
survival_rates_boot_summary <-
  Rmisc::summarySE(survival_rates_boot,
                   measurevar = "estimate",
                   groupvars = c("sex_age"),
                   conf.interval = 0.95)

# Define Ceuta vital rates estimated from mark-recapture analysis:
deterministic_list <- list(F_Chk_survl = survival_rates_boot_summary[2,3],
                           F_Fdg_survl = survival_rates_boot_summary[3,3],
                           F_Adt_survl = survival_rates_boot_summary[1,3],
                           M_Chk_survl = survival_rates_boot_summary[5,3],
                           M_Fdg_survl = survival_rates_boot_summary[6,3],
                           M_Adt_survl = survival_rates_boot_summary[4,3],
                           # Define h (harem size, h = 1 is monogamy) and k (clutch size)
                           h = 1,
                           k = 3,
                           # Define primary sex ratio (assumed to be 0.5)
                           HSR = HSR,
                           RF = RF,
                           RM = RM)

# Ceuta matrix model. "..._R" is the matrix model with a simple annual fecundity
# value for females, whereas the "..._fun" is the matrix model that includes the
# non-linear harmonic mean mating function (Note: it will be added to the matrix
# later in a mathematical expression. Here, the cell is left as "NA". 
# FYI: The "..._R" and "...fun" notation is retained throughout the script.
deterministic_matrix_R <- plover_matrix(deterministic_list, mating_function = FALSE)
deterministic_matrix_fun <- plover_matrix(deterministic_list, mating_function = TRUE)

# Determine the ASR at the stable stage distribution
deterministic_ASR_R <- 
  freq_dep_SSD_ASR(A = deterministic_matrix_R, HSR = HSR, mating_function = FALSE)
deterministic_ASR_R$ASR

deterministic_ASR_fun <- 
  freq_dep_SSD_ASR(A = deterministic_matrix_fun, h = 1, HSR = HSR, mating_function = TRUE, iterations = 1000)
deterministic_ASR_fun$ASR

# Lower-level vital rate sensitivity analysis
deterministic_LLSA_R <- 
  lower_level_sens_analysis(freq_dep_ASR = deterministic_ASR_R, 
                            VR_list = deterministic_list,
                            mating_function = FALSE)

deterministic_LLSA_fun <- 
  lower_level_sens_analysis(freq_dep_ASR = deterministic_ASR_fun, 
                            VR_list = deterministic_list,
                            mating_function = TRUE)

# Calculate vital rate sensitivities, elasticities, and LTRE
deterministic_LTRE_R <- 
  vitalsens_ASR(elements = deterministic_LLSA_R$Element_expression,
                VR_list = deterministic_list, freq_dep_ASR = deterministic_ASR_R,
                mating_function = FALSE)

deterministic_LTRE_fun <- 
  vitalsens_ASR(elements = deterministic_LLSA_fun$Element_expression,
                VR_list = deterministic_list, freq_dep_ASR = deterministic_ASR_fun,
                mating_function = TRUE)

# Custom color palette for the plotting of Fledgling and Adult stats
cbPalette <- c("#A6A6A6", "#D9D9D9", "#D9D9D9", "#D9D9D9")
#cbPalette <- c("#D9D9D9", "#A6A6A6", "#A6A6A6", "#A6A6A6")


# plot the comparative LTRE results
Background_LTRE <- 
  ggplot2::ggplot(data = deterministic_LTRE_R$LTRE,
         aes(x = parameter, y = value, fill = parameter)) + 
  coord_flip() +
  theme_bw() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(2)]) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(1)]) +
  annotate("text", x = c(2,2), y = c(-Inf, Inf),
           label = c("\u2640", "\u2642"), size = 7,
           family="Arial", vjust = c(0.5,0.5), hjust = c(-0.3,1.3)) +
  theme(text = element_text(family="Arial", color = "white"), # set the font as Candara
        legend.position = "none",
        axis.title.x = element_text(size=12, margin = margin(10, 0, 0, 0)),
        axis.text.x  = element_text(size=10, margin = margin(5, 0, 0, 0)), 
        axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
        axis.text.y  = element_text(size=10, angle = 90, hjust = 0.5, margin = margin(0, 1, 0, 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.5, colour = "white"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_blank(),
        plot.margin = unit(c(1,0.5,0.5,1.12), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  ylab("Contribution to adult sex ratio") +
  xlab("Sex-bias in parameter") +
  scale_fill_manual(values = cbPalette) +
  scale_y_continuous(limits = c(-0.04, 0.04)) +
  scale_x_discrete(labels = c("Adult survival" = expression(Adult["\u03D5"]),
                              "Fledgling survival" = expression(Fledgling ["\u03D5"]),
                              "Chick survival" = expression(Chick ["\u03D5"]),
                              "Hatching sex ratio" = "Hatching SR"))
Background_LTRE

LTRE_bar_plot_fun <- 
  ggplot2::ggplot() +
  theme_bw() +
  coord_flip() +
  geom_bar(data = deterministic_LTRE_R$LTRE,
           aes(x = parameter, y = value, fill = parameter), color = "black", stat = "identity", alpha = 0.8) +
  theme(text = element_text(family="Arial"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.x = element_text(size=12, margin = margin(10, 0, 0, 0)),
        axis.text.x  = element_text(size=10, margin = margin(5, 0, 0, 0)), 
        axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
        axis.text.y  = element_text(size=10, angle = 90, hjust = 0.5, margin = margin(0, 1, 0, 0)),
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
  ylab("Contribution to adult sex ratio") +
  xlab("Sex-bias in parameter") +
  scale_fill_manual(values = cbPalette) +
  scale_y_continuous(limits = c(-0.04, 0.04)) +
  scale_x_discrete(labels = c("Adult survival" = expression(Adult["\u03D5"]),
                              "Fledgling survival" = expression(Fledgling ["\u03D5"]),
                              "Chick survival" = expression(Chick ["\u03D5"]),
                              "Hatching sex ratio" = "Hatching SR"))
LTRE_bar_plot_fun

jpeg(filename = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/figs/LTRE_final_arial.jpeg",
     quality = 100,
     width = 4,
     height = 5,
     units = "in",
     res = 300)

grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) )
print( Background_LTRE + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( LTRE_bar_plot_fun + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
dev.off()

################################## BREEDING SYSTEM ############################
breeding_data <- 
  read.table("data/breeding_data.txt",
             header=T, stringsAsFactors = FALSE)
######## MATING BEHAVIOUR #####################################################
# remove any cases in which one mate was not identified (i.e., "NA")
mating_df <- breeding_data[which(!is.na(breeding_data$female) & !is.na(breeding_data$male)),]

# Bind the two mates together to make a unique pair
mating_df$pair <- as.factor(paste(mating_df$female, mating_df$male, sep = "-"))

# Determine how many mating attempts each individual had each year
females <- dcast(mating_df, female  ~ year)
males <- dcast(mating_df, male  ~ year)

# determine how many different mates each individual had over their lifetime in the popualtion
number_males_p_female <- aggregate(male ~ female, mating_df, function(x) length(unique(x)))
number_females_p_male <- aggregate(female ~ male, mating_df, function(x) length(unique(x)))

# Join these two dataframes together and define as numeric
females <- inner_join(females, number_males_p_female)
females[,c(2:8)] <- 
  lapply(females[,c(2:8)], as.numeric)
males <- inner_join(males, number_females_p_male)
males[,c(2:8)] <- 
  lapply(males[,c(2:8)], as.numeric)

# Calculate the total number of mating attempts over each individual's lifetime
females$attempts <- rowSums(females[, c(2:8)])
males$attempts <- rowSums(males[, c(2:8)])

# Calculate the number of years breeding
females$years <- rowSums(females[, c(2:8)] > 0)
males$years <- rowSums(males[, c(2:8)] > 0)

# Filter out all individuals that only had one mating attempt
females_no_1 <- filter(females, male  != 1 | years != 1 | attempts != 1)
males_no_1 <- filter(males, female  != 1 | years != 1 | attempts != 1)

# tidy up dataframes then bind them together
females_no_1$sex <- "Female"
females_no_1$sex <- as.factor(females_no_1$sex)
colnames(females_no_1)[c(1,9)] <- c("focal", "mate")
males_no_1$sex <- "Male"
males_no_1$sex <- as.factor(males_no_1$sex)
colnames(males_no_1)[c(1,9)] <- c("focal", "mate")
mating <- rbind(females_no_1, males_no_1)

# Determine if an individual was either:
# a) monogamous between years (i.e. only 1 mate in lifetime, with the number of attempts equaling the number years mating)
# b) monogamous within years (i.e. only 1 mate in lifetime, with the number of attempts greater the number years mating)
# c) polygamous between years (i.e. more than one mate in lifetime, with the number of attempts equaling the number years mating)
# d) polygamous within years (i.e. more than one mate in lifetime, with the number of attempts greater the number years mating)
mating$status <- ifelse(mating$mate == 1 & mating$years == mating$attempts, "Monogamous between years",
                        ifelse(mating$mate == 1 & mating$years < mating$attempts, "Monogamous within years",
                               ifelse(mating$mate > 1 & mating$years == mating$attempts, "Polygamous between years",
                                      ifelse(mating$mate > 1 & mating$years < mating$attempts, "Polygamous within years", "XXX"))))

# Calculate the number of mates per year
mating$no_mates_per_year <- mating$mate/mating$years

# Run chi-squared test of the sex-differences in polygamy rates
chisq.test(table(mating$sex, mating$status)[,c(3,4)])

# Run chi-squared test of the sex-differences in monogamy rates
chisq.test(table(mating$sex, mating$status)[,c(1,2)])

# Run chi-squared test of the sex-differences in mating behaviour rates
chisq.test(table(mating$sex, mating$status))

# Set the factor levels for plotting
mating$status <- factor(mating$status, 
                        levels = c("Monogamous within years",
                                   "Monogamous between years",
                                   "Polygamous between years",
                                   "Polygamous within years"))

# Determine the number of males and females used in the analysis
sample_sizes_sex <- aggregate(focal ~ sex, data = mating, FUN = function(x){NROW(x)})

# Define the color palatte to use in the plot
custom_pal <- c("#7b3294", "#9E6BB1", "#91bfdb", "#4575b4")

# plot the sex-differences in mating behaviour
matefidelity_plot_by_sex <- 
  ggplot() +
  geom_bar(position = "fill", alpha = 0.75, data = mating, aes(x = sex, fill = status)) +
  geom_text(aes(y = c(1.05, 1.05), x = c(1.11, 2.11), label = focal, family = "Arial"), data = sample_sizes_sex, size = 3) +
  annotate("text", x = c(0.92, 1.92), y = c(1.05, 1.05), label = "n = ", size = 3) +
  theme_bw() +
  theme(text = element_text(family = "Arial"),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.height=unit(0.8,"line"),
        legend.key.width=unit(0.8,"line"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 10), 
        axis.title.y = element_text(size = 12, margin = margin(0, 15, 0, 0)),
        axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=12),
        strip.background = element_blank(),
        strip.text = element_text(vjust = -10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(linetype = "solid", colour = "grey"),
        plot.margin = unit(c(0.2,0.2,-0.2,0.2), "cm")) +
  ylab("Proportion of individuals") +
  scale_fill_manual(values = custom_pal) +
  #facet_grid(. ~ sex) +
  scale_y_continuous(limits = c(0, 1.05)) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))
matefidelity_plot_by_sex

ggplot2::ggsave(matefidelity_plot_by_sex, 
                filename = "matefidelity_plot.jpg", 
                path = "figs/",
                width = 3,
                height = 5, units = "in",
                dpi = 300)

################
# extract the female column, add a sex column.  extract the male colum, add a 
# sex column.  Stack these two dataframes.
Sex <- rep("Female", nrow(breeding_data))
Ring <- breeding_data$female
females <- data.frame(Ring, Sex)
Sex <- rep("Male", nrow(breeding_data))
Ring <- breeding_data$male
males <- data.frame(Ring, Sex)
Individuals <- rbind(males, females)

# Replicate each row by 2 then cbind the stacked dataframe from the previous
# step
reproduction_df <- cbind(breeding_data[rep(row.names(breeding_data), 2), c("no_chicks", "clutch_size", "ID", "year")],
                      Individuals)

# Change the order of the sex levels, so that females are first (for the plot)
reproduction_df$Sex <- factor(reproduction_df$Sex, levels = c("Female", "Male"))

# subset the data to remove entries that have a NA in the Ring column
reproduction_df <- reproduction_df[!is.na(reproduction_df$Ring),]
dcast(reproduction_df, Ring ~ year)

# subset the data to remove entries that have a NA in the Ring column
reproduction_df <- reproduction_df[!is.na(reproduction_df$no_chicks),]

# group data according to Year, Sex, then Ring
reproduction_df <- group_by(reproduction_df, year, Sex, Ring)

# sum the total chicks produced per bird each year
reproduction_df_sum <- ungroup(dplyr::summarise(reproduction_df, total_chicks_p_year = sum(as.numeric(no_chicks))))

# calculate avg total chicks produced per bird in each year
fecundity_annual_summary <- Rmisc::summarySE(reproduction_df_sum, measurevar = "total_chicks_p_year", groupvars = c("Sex", "year"))

# group data according to Sex then Ring
reproduction_df_sum <- group_by(reproduction_df_sum, Sex, Ring)

# calculate avg total chicks produced per bird each year
reproduction_df_sum_avg <- ungroup(dplyr::summarise(reproduction_df_sum, avg_chicks_p_year = mean(as.numeric(total_chicks_p_year))))

# summarize the avg annual no_chicks by sex
fecundity_summary <- Rmisc::summarySE(reproduction_df_sum_avg, measurevar = "avg_chicks_p_year", groupvars = c("Sex"))

# Determine how many individuals were included in the analysis
sample_sizes_sex <- aggregate(Ring ~ Sex, data = reproduction_df_sum_avg, FUN = function(x){NROW(x)})

# specify the color pallete to use for the plot
cbPalette <- brewer.pal(8, "Dark2")[c(2,1)]

# plot distributions of fecundity by sex
Sex_specific_fucund_plot <- 
  ggplot() +
  geom_boxplot(aes(y = avg_chicks_p_year, x = Sex, fill = Sex), data = reproduction_df_sum_avg, size = .3, alpha = 0.6) +
  geom_text(aes(y = c(6.5, 6.5), x = c(1.12, 2.12), label = Ring, family = "Arial"), data = sample_sizes_sex, size = 3) +
  annotate("text", x = c(0.92, 1.92), y = c(6.5, 6.5), label = "n = ", size = 3) +
  theme_bw() +
  theme(text=element_text(size=16, family="Arial"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 10), 
        axis.title.y = element_text(size = 12, margin = margin(0, 15, 0, 0)),
        axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "grey")) +
  scale_fill_manual(values = cbPalette) +
  ylab("Per capita annual hatchlings parented") +
  scale_y_continuous(limits = c(0, 6.5))
Sex_specific_fucund_plot

ggsave(Sex_specific_fucund_plot,
       filename = "sex-specific_reproductive_success_tall.jpg",
       path = "figs/",
       width = 2.8,
       height = 5, units = "in",
       dpi = 300,
       scale = 1)

# Run F-test to assess sex-specific variation in per capita fecundity
reproduction_df_sum_avg <- as.data.frame(reproduction_df_sum_avg)
var.test(reproduction_df_sum_avg[which(reproduction_df_sum_avg$Sex == "Female"),c("avg_chicks_p_year")],
         reproduction_df_sum_avg[which(reproduction_df_sum_avg$Sex == "Male"),c("avg_chicks_p_year")],
         alternative = "greater")