library(ggplot2)
library(grid)

# load the following dataframe from the github repo
survival_rates_boot <-
  read.table("output/bootstrap/survival_rates_boot_out.txt", header = TRUE)

# number of iterations used in bootstrap
niter <- 1000
# Hatching sex ratio (proportion of hatchlings that are male)
HSR <- 0.4856322
# mating function parameters to specify
h <- 1.215808#0.8224983
k <- 3

###############################################################################
####### functions needed in the script ########################################

# this function builds the matrix based on a list of specified vital rates
plover_matrix <-
  function(vital_rates){
    # Define plover life-stages of the Ceuta snowy plover matrix model
    stages <- c("F_1st_yr",  "F_Adt",  "M_1st_yr",  "M_Adt")
    # Build the 4x4 matrix (note: NAs are assigned to cells where the mating
    # function will go later)
    result <-
      matrix(c(0, NA, 0, NA,
               (vital_rates$F_Chk_survl*vital_rates$F_Fdg_survl),
               vital_rates$F_Adt_survl,
               0, 0,
               0, NA, 0, NA,
               0, 0,
               (vital_rates$M_Chk_survl*vital_rates$M_Fdg_survl),
               vital_rates$M_Adt_survl),
             nrow = length(stages), byrow = TRUE,
             dimnames = list(stages, stages))
    result
  }

# this function calculates the ASR of the matrix from the matrix created from
# the previous function (to assess if the SSD is chaotic, write "plot = TRUE")
matrix_ASR <-
  function(A, n = rep(10, nrow(A)), h = 1, k = 3, iterations = 1000, HSR = 0.5, plot = FALSE){
    # Number of stages in matrix
    x <- length(n)
    # Number of time steps to simulate
    t <- iterations
    # an empty t by x matrix to store the stage distributions
    stage <- matrix(numeric(x * t), nrow = x)
    # an empty t vector to store the population sizes
    pop <- numeric(t)
    # for loop that goes through each of t time steps
    for (i in 1:t) {
      # stage distribution at time t
      stage[,i] <- n
      # population size at time t
      pop[i] <- sum(n)
      # number of male adults at time t
      M2 <- stage[4, i]
      # number of female adults at time t
      F2 <- stage[2, i]
      # Female freq-dep fecundity of Female chicks
      A[1,x/2]        <- (k*M2)/(M2+(F2*h))*HSR
      # Female freq-dep fecundity of Male chicks
      A[(x/4)*3,x/2]  <- (k*M2)/(M2+(F2*h))*HSR
      # Male freq-dep fecundity of Female chicks
      A[1,x]          <- (k*F2)/(M2+(F2*h))*HSR
      # Male freq-dep fecundity of Male chicks
      A[(x/4)*3,x]    <- (k*F2)/(M2+(F2*h))*HSR
      # define the new n (i.e., new stage distribution at time t)
      n <- A %*% n
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(t - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, t]
    }
    # calc ASR as the proportion of the adult stable stage class that is male
    ASR <- stable.stage[x]/(stable.stage[x/2] + stable.stage[x])

    if(plot)
    {
      # plot distrubution to assure that it is not chaotic
      matplot(rownames(t(stage)), t(stage), type='l', lwd=2, las=1)
    }
    # make a list of results
    pop.proj <- list(ASR = ASR,
                     lambda = pop[t]/pop[t - 1],
                     stable.stage = stable.stage,
                     stage.vectors = stage,
                     SSD_M2 = stable.stage[4],
                     SSD_F2 = stable.stage[2])
    # print the list as output to the function
    pop.proj
  }

# this function extracts the ASR and lambda estimates from each bootstrap iteration using
# a matrix with the tow-sex mating function
lambda_and_ASR_extract <-
  function(survival_rates, niter = 1000){
    # make an empty datarame to store the results
    output <- data.frame(treat_ASR = numeric(niter),
                         treat_lambda = numeric(niter),
                         contr_ASR = numeric(niter),
                         contr_lambda = numeric(niter))
    # for loop to go through each iteration and calculate the differece between female and male
    # survival rates for each stage.
    for(i in 1:niter){
      # Create a list of demographic rates from the survival analyses above
      VR_treat_boot <- list(F_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][5],
                            F_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][3],
                            F_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][1],
                            M_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][6],
                            M_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][4],
                            M_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][2],
                            # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                            h = h,
                            k = k,
                            # Define primary sex ratio (assumed to be 0.5)
                            HSR = HSR)

      VR_contr_boot <- list(F_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][6],
                            F_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][4],
                            F_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][2],
                            M_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][6],
                            M_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][4],
                            M_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][2],
                            # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                            h = h,
                            k = k,
                            # Define primary sex ratio (assumed to be 0.5)
                            HSR = 0.5)

      # Build matrix based on rates specified in the list above
      treat_matrix_boot <- plover_matrix(VR_treat_boot)
      contr_matrix_boot <- plover_matrix(VR_contr_boot)

      # Determine the ASR at the stable stage distribution
      treat_lambda_and_ASR_analysis <-
        matrix_ASR(A = treat_matrix_boot, h = VR_treat_boot$h, HSR = VR_treat_boot$HSR, iterations = 75)
      contr_lambda_and_ASR_analysis <-
        matrix_ASR(A = contr_matrix_boot, h = VR_contr_boot$h, HSR = VR_contr_boot$HSR, iterations = 75)
      # Extract ASR
      output[i, 1] <- treat_lambda_and_ASR_analysis$ASR
      output[i, 2] <- treat_lambda_and_ASR_analysis$lambda
      output[i, 3] <- contr_lambda_and_ASR_analysis$ASR
      output[i, 4] <- contr_lambda_and_ASR_analysis$lambda
    }
    # restructure the output and lable columns
    output <- suppressMessages(reshape2::melt(data = output))
    colnames(output) <- c("parameter", "estimate")
    # return the output
    output
  }

# this function runs a perturbation analysis on the two-sex non-linear matrix to
# calculate the sensitivites and elasticities of each parameter included in the
# matrix (including the components of the mating function)
sensitivity_analysis <-
  function(vital_rates, matrix_str, h = 1, k = 3, HSR, niter = 1000, ASR, lambda){
    # make a list of all parameters
    vr <-
      list(F_Chk_survl = vital_rates$F_Chk_survl,
           F_Fdg_survl = vital_rates$F_Fdg_survl,
           F_Adt_survl = vital_rates$F_Adt_survl,
           M_Chk_survl = vital_rates$M_Chk_survl,
           M_Fdg_survl = vital_rates$M_Fdg_survl,
           M_Adt_survl = vital_rates$M_Adt_survl)

    # number of stages in the matrix
    no_stages <- sqrt(length(matrix_str))

    # Define plover life-stages of the Ceuta snowy plover matrix model
    stages <- c("F_1st_yr",  "F_Adt",  "M_1st_yr",  "M_Adt")

    # an empty t by x matrix
    stage <- matrix(numeric(no_stages * niter), nrow = no_stages)

    # an empty t vector to store the population sizes
    pop <- numeric(niter)

    # dataframe to store the perturbation results
    ASR_pert_results <-
      data.frame(parameter = c("F_Chk_survl", "F_Fdg_survl", "F_Adt_survl",
                               "M_Chk_survl", "M_Fdg_survl", "M_Adt_survl",
                               "h", "k", "HSR"),
                 sensitivities = numeric(9),
                 elasticities = numeric(9))

    lambda_pert_results <-
      data.frame(parameter = c("F_Chk_survl", "F_Fdg_survl", "F_Adt_survl",
                               "M_Chk_survl", "M_Fdg_survl", "M_Adt_survl",
                               "h", "k", "HSR"),
                 sensitivities = numeric(9),
                 elasticities = numeric(9))

    # specifiy how many survival rates there are
    n <- length(vr)

    # create vectors of perturbations to test on parameters of the matrix model
    vr_nums <- seq(0, 1, 0.01) # proportional changes in survival and HSR (i.e., between 0 an 1)
    h_nums <- seq(0, 2, 0.02) # proportional changes in h index (i.e., between 0 and 2)
    k_nums <- seq(2, 4, 0.02) # proportional changes in k (i.e, between 2 and 4)

    # create empty dataframes to store the perturbation results
    vr_pert_ASR <- matrix(numeric(n * length(vr_nums)),
                      ncol = n, dimnames = list(vr_nums, names(vr)))
    h_pert_ASR <- matrix(numeric(length(h_nums)),
                     ncol = 1, dimnames = list(h_nums, "h"))
    k_pert_ASR <- matrix(numeric(length(k_nums)),
                     ncol = 1, dimnames = list(k_nums, "k"))
    HSR_pert_ASR <- matrix(numeric(length(vr_nums)),
                       ncol = 1, dimnames = list(vr_nums, "HSR"))
    vr_pert_lambda <- matrix(numeric(n * length(vr_nums)),
                          ncol = n, dimnames = list(vr_nums, names(vr)))
    h_pert_lambda <- matrix(numeric(length(h_nums)),
                         ncol = 1, dimnames = list(h_nums, "h"))
    k_pert_lambda <- matrix(numeric(length(k_nums)),
                         ncol = 1, dimnames = list(k_nums, "k"))
    HSR_pert_lambda <- matrix(numeric(length(vr_nums)),
                           ncol = 1, dimnames = list(vr_nums, "HSR"))

    # perturbation of vital rates survival rates
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      vr2 <- vr # reset the vital rates to the original
      for (i in 1:length(vr_nums)) # pick a perturbation level
      {
        vr2[[g]] <- vr_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr2, NULL), nrow = sqrt(length(matrix_str)), byrow=TRUE, dimnames = list(stages, stages)) # build the matrix with the new value
        m <- rep(10, no_stages) # reset the starting stage distribution for simulation (all with 10 individuals)
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          # number of male adults at time t
          M2 <- stage[4, j]
          # number of female adults at time t
          F2 <- stage[2, j]
          # Female freq-dep fecundity of Female chicks
          A[1,no_stages/2]        <- ((k*M2)/(M2+(F2/h)))*HSR
          # Female freq-dep fecundity of Male chicks
          A[(no_stages/4)*3,no_stages/2]  <- ((k*M2)/(M2+(F2/h)))*HSR
          # Male freq-dep fecundity of Female chicks
          A[1,no_stages]          <- ((k*F2)/(M2+(F2/h)))*HSR
          # Male freq-dep fecundity of Male chicks
          A[(no_stages/4)*3,no_stages]    <- ((k*F2)/(M2+(F2/h)))*HSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        vr_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
        # calc lambda as the pop change in the counts of the last two iterations
        vr_pert_lambda[i, g] <- pop[niter]/pop[niter - 1]
      }
      spl_ASR <- smooth.spline(vr_pert_ASR[,g] ~ rownames(vr_pert_ASR))
      ASR_pert_results[g, 2] <- predict(spl_ASR, x=vr[[g]], deriv=1)$y
      ASR_pert_results[g, 3] <- vr[[g]]/ASR * ASR_pert_results[g, 2]

      spl_lambda <- smooth.spline(vr_pert_lambda[,g] ~ rownames(vr_pert_lambda))
      lambda_pert_results[g, 2] <- predict(spl_lambda, x=vr[[g]], deriv=1)$y
      lambda_pert_results[g, 3] <- vr[[g]]/lambda * lambda_pert_results[g, 2]
    }
    # perturbation of the h index parameter
    for (i in 1:length(h_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), byrow=TRUE, dimnames = list(stages, stages)) # build the matrix with the new value
      m <- rep(10, no_stages) # reset the starting stage distribution for simulation (all with 10 individuals)
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        # number of male adults at time t
        M2 <- stage[4, j]
        # number of female adults at time t
        F2 <- stage[2, j]
        # Female freq-dep fecundity of Female chicks
        A[1,no_stages/2]        <- ((k*M2)/(M2+(F2/h_nums[i])))*HSR
        # Female freq-dep fecundity of Male chicks
        A[(no_stages/4)*3,no_stages/2]  <- ((k*M2)/(M2+(F2/h_nums[i])))*HSR
        # Male freq-dep fecundity of Female chicks
        A[1,no_stages]          <- ((k*F2)/(M2+(F2/h_nums[i])))*HSR
        # Male freq-dep fecundity of Male chicks
        A[(no_stages/4)*3,no_stages]    <- ((k*F2)/(M2+(F2/h_nums[i])))*HSR
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      h_pert_ASR[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
      # calc lambda as the pop change in the counts of the last two iterations
      h_pert_lambda[i, ] <- pop[niter]/pop[niter - 1]
    }
    spl_ASR <- smooth.spline(h_pert_ASR[, 1] ~ rownames(h_pert_ASR))
    ASR_pert_results[n+1, 2] <- predict(spl_ASR, x=h, deriv=1)$y
    ASR_pert_results[n+1, 3] <- h/ASR * ASR_pert_results[n+1, 2]

    spl_lambda <- smooth.spline(h_pert_lambda[,1] ~ rownames(h_pert_lambda))
    lambda_pert_results[n+1, 2] <- predict(spl_lambda, x=h, deriv=1)$y
    lambda_pert_results[n+1, 3] <- h/lambda * lambda_pert_results[n+1, 2]

    # perturbation of k parameter
    for (i in 1:length(k_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), byrow=TRUE, dimnames = list(stages, stages)) # build the matrix with the new value
      m <- rep(10, no_stages) # reset the starting stage distribution for simulation (all with 10 individuals)
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        # number of male adults at time t
        M2 <- stage[4, j]
        # number of female adults at time t
        F2 <- stage[2, j]
        # Female freq-dep fecundity of Female chicks
        A[1,no_stages/2]        <- ((k_nums[i]*M2)/(M2+(F2/h)))*HSR
        # Female freq-dep fecundity of Male chicks
        A[(no_stages/4)*3,no_stages/2]  <- ((k_nums[i]*M2)/(M2+(F2/h)))*HSR
        # Male freq-dep fecundity of Female chicks
        A[1,no_stages]          <- ((k_nums[i]*F2)/(M2+(F2/h)))*HSR
        # Male freq-dep fecundity of Male chicks
        A[(no_stages/4)*3,no_stages]    <- ((k_nums[i]*F2)/(M2+(F2/h)))*HSR
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      k_pert_ASR[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
      # calc lambda as the pop change in the counts of the last two iterations
      k_pert_lambda[i, ] <- pop[niter]/pop[niter - 1]
    }
    spl_ASR <- smooth.spline(k_pert_ASR[,1] ~ rownames(k_pert_ASR))
    ASR_pert_results[n+2, 2] <- predict(spl_ASR, x=k, deriv=1)$y
    ASR_pert_results[n+2, 3] <- k/ASR * ASR_pert_results[n+2, 2]

    spl_lambda <- smooth.spline(k_pert_lambda[,1] ~ rownames(k_pert_lambda))
    lambda_pert_results[n+2, 2] <- predict(spl_lambda, x=k, deriv=1)$y
    lambda_pert_results[n+2, 3] <- k/lambda * lambda_pert_results[n+2, 2]
    # perturbation of HSR
    for (i in 1:length(vr_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), byrow=TRUE, dimnames = list(stages, stages)) # build the matrix with the new value
      m <- rep(10, no_stages) # reset the starting stage distribution for simulation (all with 10 individuals)
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        # number of male adults at time t
        M2 <- stage[4, j]
        # number of female adults at time t
        F2 <- stage[2, j]
        # Female freq-dep fecundity of Female chicks
        A[1,no_stages/2]        <- ((k*M2)/(M2+(F2/h)))*vr_nums[i]
        # Female freq-dep fecundity of Male chicks
        A[(no_stages/4)*3,no_stages/2]  <- ((k*M2)/(M2+(F2/h)))*vr_nums[i]
        # Male freq-dep fecundity of Female chicks
        A[1,no_stages]          <- ((k*F2)/(M2+(F2/h)))*vr_nums[i]
        # Male freq-dep fecundity of Male chicks
        A[(no_stages/4)*3,no_stages]    <- ((k*F2)/(M2+(F2/h)))*vr_nums[i]
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      HSR_pert_ASR[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
      # calc lambda as the pop change in the counts of the last two iterations
      HSR_pert_lambda[i, ] <- pop[niter]/pop[niter - 1]
    }
    spl_ASR <- smooth.spline(HSR_pert_ASR[,1] ~ rownames(HSR_pert_ASR))
    ASR_pert_results[n+3, 2] <- predict(spl_ASR, x=HSR, deriv=1)$y
    ASR_pert_results[n+3, 3] <- HSR/ASR * ASR_pert_results[n+3, 2]

    spl_lambda <- smooth.spline(HSR_pert_lambda[,1] ~ rownames(HSR_pert_lambda))
    lambda_pert_results[n+3, 2] <- predict(spl_lambda, x=HSR, deriv=1)$y
    lambda_pert_results[n+3, 3] <- HSR/lambda * lambda_pert_results[n+3, 2]

    result <- list(ASR_pert_results = ASR_pert_results,
                   lambda_pert_results = lambda_pert_results)
  }

# this function estimates the contribution that each vital rate has on ASR bias,
# given the sensitivities calculated in the previous function
# (see formula 8 on page 133 of Veran and Beissinger (2009))
LTRE_analysis <-
  function(Mprime_sens, matrix_str, vital_rates){
    # make an empty dataframe to stroe LTRE results
    LTRE_ASR <-
      data.frame(parameter = c("Chick survival", "Fledgling survival",
                               "Adult survival", "Hatching sex ratio",
                               "Mating system"),
                 contribution = numeric(5))

    LTRE_lambda <-
      data.frame(parameter = c("Chick survival", "Fledgling survival",
                               "Adult survival", "Hatching sex ratio",
                               "Mating system"),
                 contribution = numeric(5))

    # run a for loop to extract the parameter contributions
    for(i in 1:nrow(LTRE_ASR))
    {
      LTRE_ASR[i, 2] <-
        ifelse(i < 4, (vital_rates[[i + 3]] - vital_rates[[i]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 3],
               ifelse(i == 4, ((1-vital_rates[[9]]) - vital_rates[[9]]) * Mprime_sens$ASR_pert_results$sensitivities[9],
                      (1 - vital_rates[[7]]) * Mprime_sens$ASR_pert_results$sensitivities[7]))
    }

    for(i in 1:nrow(LTRE_lambda))
    {
      LTRE_lambda[i, 2] <-
        ifelse(i < 4, (vital_rates[[i + 3]] - vital_rates[[i]]) * Mprime_sens$lambda_pert_results$sensitivities[i + 3],
               ifelse(i == 4, (vital_rates[[9]] - (1-vital_rates[[9]])) * Mprime_sens$lambda_pert_results$sensitivities[9],
                      (vital_rates[[7]] - 1) * Mprime_sens$lambda_pert_results$sensitivities[7]))
    }

    LTRE_ASR$parameter <- factor(LTRE_ASR$parameter, levels = c("Adult survival",
                                                                "Fledgling survival",
                                                                "Chick survival",
                                                                "Hatching sex ratio",
                                                                "Mating system"))

    LTRE_lambda$parameter <- factor(LTRE_lambda$parameter, levels = c("Adult survival",
                                                                      "Fledgling survival",
                                                                      "Chick survival",
                                                                      "Hatching sex ratio",
                                                                      "Mating system"))

    LTRE_results <- list(LTRE_ASR = LTRE_ASR,
                         LTRE_lambda = LTRE_lambda)
  }

bootstrap_matrix_extract <-
  function(survival_rates, niter = 1000){

    # make empty lists for the matricies
    treat_mats <- vector("list", niter)
    contr_mats <- vector("list", niter)

    for(i in 1:niter){
      # Create a list of demographic rates from the survival analyses above
      VR_treat_boot <- list(F_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][5],
                            F_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][3],
                            F_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][1],
                            M_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][6],
                            M_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][4],
                            M_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][2],
                            # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                            h = h,
                            k = k,
                            # Define primary sex ratio (assumed to be 0.5)
                            HSR = HSR)

      VR_contr_boot <- list(F_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][6],
                            F_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][4],
                            F_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][2],
                            M_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][6],
                            M_Fdg_survl = survival_rates[which(survival_rates$iter == i), 2][4],
                            M_Adt_survl = survival_rates[which(survival_rates$iter == i), 2][2],
                            # Define h (harem size, h < 1 is polyandry) and k (clutch size)
                            h = 1,
                            k = 3,
                            # Define primary sex ratio (assumed to be 0.5)
                            HSR = 0.5)

      # Build matrix based on rates specified in the list above
      treat_mats[[i]] <- plover_matrix(VR_treat_boot)
      contr_mats[[i]]  <- plover_matrix(VR_contr_boot)
    }
      # make a list of matrix lists
      output_mats <- list(treat_mats = treat_mats,
                          contr_mats = contr_mats)
      output_mats
  }

###############################################################################
####### ASR estimate with two-sex matrix that includes mating function ########

# extract the ASR estimate for each iteration of the bootstrapped survival
# analysis
lambda_and_ASR_boot <- lambda_and_ASR_extract(survival_rates_boot)

ASR_boot <- data.frame(filter(lambda_and_ASR_boot, parameter == "treat_ASR")[2], c(1:1000))
colnames(ASR_boot) <- c("ASR_boot", "iter")
write.table(ASR_boot,
            file = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/output/bootstrap/ASR_boot_out.txt",
            row.names = FALSE, col.names = TRUE)


CI <- 0.95

# Calculate the mean, and median of the ASR bootstraps
lambda_ASR_boot_summary <-
  lambda_and_ASR_boot%>%
  dplyr::group_by(parameter)%>%
  dplyr::summarise(mean = mean(estimate),
                   lcl = quantile(estimate, (1 - CI)/2),
                   ucl = quantile(estimate, 1 - (1 - CI)/2))

lambda_diff_boot_summary <-
  lambda_diff%>%
  dplyr::summarise(mean = mean(lambda_diff),
                   lcl = quantile(lambda_diff, (1 - CI)/2),
                   ucl = quantile(lambda_diff, 1 - (1 - CI)/2))

# Calculate the lambda difference between the treatment and control matrix
lambda_difference <- 
  function(niter = 1000, lambda_and_ASR_boot)
    {
      lambda_diff <- matrix(numeric(niter * 2), nrow = niter)
      for(i in 1:niter)
        {
        lambda_diff[i,1] <- i
        lambda_diff[i,2] <- lambda_and_ASR_boot$estimate[niter+i]-lambda_and_ASR_boot$estimate[3*niter+i]
        }
      lambda_diff <- as.data.frame(lambda_diff)
      colnames(lambda_diff) <- c("iter", "lambda_diff")
      lambda_diff
    }

lambda_diff <- lambda_difference(lambda_and_ASR_boot = lambda_and_ASR_boot)

# plot the histogram distribution of ASRs calculated from each bootstrap iteration
ASR_bootstrap_histogram <- 
  ggplot() +
  annotate("rect", xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(2)]) +
  annotate("rect", xmin=0.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(1)]) +
  annotate("text", x = c(-Inf,Inf), y = c(80, 80),
           label = c("\u2640", "\u2642"), size = 7,
           family="Arial", vjust = c(1.5,1.5), hjust = c(-0.5,1.5)) +
  geom_histogram(binwidth = 0.02, data = lambda_and_ASR_boot[which(lambda_and_ASR_boot$parameter == "treat_ASR"),], 
                 aes(x = estimate)) +
  geom_errorbarh(data = lambda_ASR_boot_summary[1,], 
                 aes(y = 125, x = lcl, xmin = lcl, xmax = ucl), color = "black", size = 0.8, linetype = "solid") +
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
        axis.text.y  = element_text(size=10, angle = 90, hjust = 0.5, margin = margin(0, 5, 0, 0), color = "black"),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
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
  scale_y_continuous(limits = c(0, 125))
ASR_bootstrap_histogram

ggsave(ASR_bootstrap_histogram,
       filename = "ASR_distribution_mating_function.jpg",
       path = "figs/final/",
       width = 4,
       height = 3, units = "in",
       dpi = 300,
       scale = 1)

# plot the histogram distribution of lambdas calculated from each bootstrap iteration
lambda_bootstrap_histogram <- 
  ggplot() +
  # annotate("rect", xmin=-Inf, xmax=1, ymin=-Inf, ymax=Inf, alpha=0.4,
  #          fill=brewer.pal(8, "Set1")[c(1)]) +
  # annotate("rect", xmin=1, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.4,
  #          fill=brewer.pal(8, "Set1")[c(3)]) +
  # annotate("text", x = c(-Inf,Inf), y = c(105, 105),
  #          label = c("-\u03BB", "+\u03BB"), size = 7,
  #          family="Arial", vjust = c(1.5,1.5), hjust = c(-0.5,1.5)) +
  geom_histogram(binwidth = 0.01, data = lambda_and_ASR_boot[which(lambda_and_ASR_boot$parameter == "contr_lambda"),],
                 aes(x = estimate), fill=brewer.pal(8, "Set1")[c(1)], alpha = 0.7) +
  geom_histogram(binwidth = 0.01, data = lambda_and_ASR_boot[which(lambda_and_ASR_boot$parameter == "treat_lambda"),],
                 aes(x = estimate), alpha = 0.7) +
  geom_errorbarh(data = lambda_ASR_boot_summary[2,],
                 aes(y = 180, x = lcl, xmin = lcl, xmax = ucl), color = "black", size = 0.8, linetype = "solid") +
  geom_errorbarh(data = lambda_ASR_boot_summary[4,],
                 aes(y = 177, x = lcl, xmin = lcl, xmax = ucl), color = "red", size = 0.8, linetype = "solid") +
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
  xlab("Population growth rate (\u03BB)") +
  scale_x_continuous(limits = c(0.75, 1.25)) +
  scale_y_continuous(limits = c(0, 180))
lambda_bootstrap_histogram

# ggsave(lambda_bootstrap_histogram,
#        filename = "lambda_distribution_mating_function.jpg",
#        path = "figs/final/",
#        width = 4,
#        height = 3, units = "in",
#        dpi = 300,
#        scale = 1)

lambda_diff_bootstrap_histogram <- 
  ggplot() +
  geom_histogram(binwidth = 0.01, data = lambda_diff,
                 aes(x = lambda_diff), alpha = 0.7) +
  geom_errorbarh(data = lambda_diff_boot_summary[1,],
                 aes(y = 177, x = lcl, xmin = lcl, xmax = ucl), color = "black", size = 0.8, linetype = "solid") +
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
  xlab("Difference in population growth rate (\u03BB)") +
  scale_y_continuous(limits = c(0, 180))
lambda_diff_bootstrap_histogram

########### Life Table Response Experiment ####################################
# Summarise the bootstrap stage- and sex-specific survival rates for the
# deterministic matrix
survival_rates_boot_summary <-
  Rmisc::summarySE(survival_rates_boot,
                   measurevar = "estimate",
                   groupvars = c("sex_age"),
                   conf.interval = 0.95)

# Define deterministic vital rates estimated from mark-recapture analysis
VR_treat <- list(F_Chk_survl = survival_rates_boot_summary[2,3],
                 F_Fdg_survl = survival_rates_boot_summary[3,3],
                 F_Adt_survl = survival_rates_boot_summary[1,3],
                 M_Chk_survl = survival_rates_boot_summary[5,3],
                 M_Fdg_survl = survival_rates_boot_summary[6,3],
                 M_Adt_survl = survival_rates_boot_summary[4,3],
                 # Define h (harem size, h = 1 is monogamy) and k (clutch size)
                 h = h,
                 k = k,
                 # Define primary sex ratio
                 HSR = HSR)

# Define vital rates of the M prime matrix (i.e., average between a "control 
# matrix" and the "treatment matrix"). The control matrix is a matrix in which
# the female vital rates are set to the male vital rates, and the treatment
# matrix is the matrix containing the sex-specific values estimated from the field
# (see formula 8 on page 133 of Veran and Beissinger (2009))
VR_mprime <- list(F_Chk_survl = (survival_rates_boot_summary[2,3] +
                                   survival_rates_boot_summary[5,3])/2,
                  F_Fdg_survl = (survival_rates_boot_summary[3,3] +
                                   survival_rates_boot_summary[6,3])/2,
                  F_Adt_survl = (survival_rates_boot_summary[1,3] +
                                   survival_rates_boot_summary[4,3])/2,
                  M_Chk_survl = (survival_rates_boot_summary[5,3] +
                                   survival_rates_boot_summary[5,3])/2,
                  M_Fdg_survl = (survival_rates_boot_summary[6,3] +
                                   survival_rates_boot_summary[6,3])/2,
                  M_Adt_survl = (survival_rates_boot_summary[4,3] +
                                   survival_rates_boot_summary[4,3])/2,
                  # Define h (harem size, h = 1 is monogamy) and k (clutch size)
                  h = (h+1)/2,
                  k = k,
                  # Define primary sex ratio
                  HSR = (HSR+0.5)/2)

# define the structure of the two-sex matrix (Note: NA's specify where the mating 
# function will go)
matrix_str <- expression(0, NA, 0, NA,
                         (F_Chk_survl * F_Fdg_survl), F_Adt_survl, 0, 0,
                         0, NA, 0, NA,
                         0, 0, (M_Chk_survl * M_Fdg_survl), M_Adt_survl)

# build the treatment matrix
treatment_matrix <- 
  plover_matrix(VR_treat)

# build the M-Prime matrix 
# (see formula 8 on page 133 of Veran and Beissinger (2009))
M_prime_matrix <- 
  plover_matrix(VR_mprime)

# calculate the ASR of each matrix
treatment_ASR_analysis <- 
  matrix_ASR(A = treatment_matrix, h = h, HSR = VR_treat$HSR, iterations = 1000)

M_prime_ASR_analysis <- 
  matrix_ASR(A = M_prime_matrix, h = 1, HSR = VR_mprime$HSR, iterations = 1000)

# specify the ASR of the treatment and M-prime matricies
ASR_treat <- 
  treatment_ASR_analysis$ASR

ASR_mprime <- 
  M_prime_ASR_analysis$ASR

# specify the lambda of the treatment and the M-prime matricies
lambda_treat <- 
  treatment_ASR_analysis$lambda

lambda_mprime <- 
  M_prime_ASR_analysis$lambda

# conduct a sensitivity analysis on the treatment matrix
treat_sensitivity_analysis <- 
  sensitivity_analysis(vital_rates = VR_treat, 
                       matrix_str = matrix_str, 
                       h = VR_treat$h, 
                       k = VR_treat$k, 
                       HSR = VR_treat$HSR, 
                       niter = 1000, 
                       ASR = ASR_treat,
                       lambda = lambda_treat)

# conduct a sensitivity analysis on the M-Prime matrix
Mprime_sensitivity_analysis <- 
  sensitivity_analysis(vital_rates = VR_mprime, 
                       matrix_str = matrix_str, 
                       h = VR_mprime$h, 
                       k = VR_mprime$k, 
                       HSR = VR_mprime$HSR, 
                       niter = 1000, 
                       ASR = ASR_mprime,
                       lambda = lambda_mprime)

# conduct the LTRE comparing the two matrices
LTRE_plover <- 
  LTRE_analysis(Mprime_sens = Mprime_sensitivity_analysis, 
                matrix_str = matrix_str, 
                vital_rates = VR_treat)

sum(abs(LTRE_plover$LTRE_ASR$contribution))/c(ASR_treat-ASR_mprime)
sum(abs(LTRE_plover$LTRE_lambda$contribution))/(lambda_treat-lambda_mprime)


# specify the color palette to use in LTRE
cbPalette <- c("#A6A6A6", "#D9D9D9", "#D9D9D9", "#D9D9D9", "#A6A6A6")

# plot the LTRE results
# first draw the background (i.e., the colors)
Background_LTRE_ASR <-
  ggplot2::ggplot(data = LTRE_plover$LTRE_ASR,
                  aes(x = parameter, y = contribution, fill = parameter)) +
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
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  scale_x_discrete(labels = c("Adult survival" = expression(Adult["\u03D5"]),
                              "Fledgling survival" = expression(Fledgling ["\u03D5"]),
                              "Chick survival" = expression(Chick ["\u03D5"]),
                              "Hatching sex ratio" = "Hatching SR",
                              "Mating system" = "Mat. sys."))
# second draw the LTRE output
LTRE_ASR <-
  ggplot2::ggplot() +
  theme_bw() +
  coord_flip() +
  geom_bar(data = LTRE_plover$LTRE_ASR,
           aes(x = parameter, y = contribution, fill = parameter), color = "black", stat = "identity") +
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
  scale_y_continuous(limits = c(-0.10, 0.10)) +
  scale_x_discrete(labels = c("Adult survival" = expression(Adult["\u03D5"]),
                              "Fledgling survival" = expression(Fledgling ["\u03D5"]),
                              "Chick survival" = expression(Chick ["\u03D5"]),
                              "Hatching sex ratio" = "Hatching SR",
                              "Mating system" = "Mat. sys."))
# Save the plot
jpeg(filename = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/figs/final/LTRE_mating_function_ASR.jpeg",
     quality = 100,
     width = 4,
     height = 5,
     units = "in",
     res = 300)

# draw the background and the LTRE on top of eachother for the final plot
grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) )
print( Background_LTRE_ASR + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( LTRE_ASR + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
dev.off()

Background_LTRE_lambda <-
  ggplot2::ggplot(data = LTRE_plover$LTRE_lambda,
                  aes(x = parameter, y = contribution, fill = parameter)) +
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
  scale_y_continuous(limits = c(-0.025, 0.025)) +
  scale_x_discrete(labels = c("Adult survival" = expression(Adult["\u03D5"]),
                              "Fledgling survival" = expression(Fledgling ["\u03D5"]),
                              "Chick survival" = expression(Chick ["\u03D5"]),
                              "Hatching sex ratio" = "Hatching SR",
                              "Mating system" = "Mat. sys."))
# second draw the LTRE output
LTRE_lambda <-
  ggplot2::ggplot() +
  theme_bw() +
  coord_flip() +
  geom_bar(data = LTRE_plover$LTRE_lambda,
           aes(x = parameter, y = contribution, fill = parameter), color = "black", stat = "identity") +
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
  ylab("Contribution to population growth rate (\u03BB)") +
  xlab("Sex-bias in parameter") +
  scale_fill_manual(values = cbPalette) +
  scale_y_continuous(limits = c(-0.025, 0.025)) +
  scale_x_discrete(labels = c("Adult survival" = expression(Adult["\u03D5"]),
                              "Fledgling survival" = expression(Fledgling ["\u03D5"]),
                              "Chick survival" = expression(Chick ["\u03D5"]),
                              "Hatching sex ratio" = "Hatching SR",
                              "Mating system" = "Mat. sys."))
#Save the plot
jpeg(filename = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/figs/final/LTRE_mating_function_lambda.jpeg",
     quality = 100,
     width = 4,
     height = 5,
     units = "in",
     res = 300)

# draw the background and the LTRE on top of eachother for the final plot
grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) )
print( Background_LTRE_lambda + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( LTRE_lambda + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
dev.off()

# check to see if the sum of LTRE contributions equals the difference
# between the ASR of the treatment matrix and the M-prime matrix
sum(abs(LTRE_plover$LTRE_ASR$contribution)) # ---> 0.1131799

abs(ASR_treat-ASR_mprime) # ---> 0.07136055

# check the amount of the difference explained by M-prime
ASR_mprime/ASR_treat # ---> 0.89%

sum(abs(LTRE_plover$LTRE_lambda$contribution)) # ---> 0.06294061

abs(lambda_treat-lambda_mprime) # ---> 0.04437934

# check the amount of the difference explained by M-prime
lambda_mprime/lambda_treat # ---> 105.3838%

# calculate relative contributions
# fledglings vs. chicks
LTRE_plover[2,2]/LTRE_plover[1,2]
# fledglings vs. adults
LTRE_plover[2,2]/LTRE_plover[3,2]


############# STEP BY STEP THROUGH THE SENSITIVITY ANALYSIS ###################
# Illustrate how the sensitivity analysis of the non-linear matrix works
# Note: This example only calculates the sensitivity of Female Adult Survival 
# (i.e., "F_Adt_survl") and so some features of the sensitivity_analysis() 
# function above are removed for brevity and clarity in the example

# first specify all the parameters (these are calculated or assigned earlier
# in this script)
vital_rates = VR_treat
matrix_str = matrix_str
h = VR_treat$h
k = VR_treat$k
HSR = VR_treat$HSR
niter = 1000
ASR = ASR_treat
lambda = lambda_treat

# Here's the flow through the function:
# make a list of all survival parameters
vr <- 
  list(F_Chk_survl = vital_rates$F_Chk_survl,
       F_Fdg_survl = vital_rates$F_Fdg_survl,
       F_Adt_survl = vital_rates$F_Adt_survl,
       M_Chk_survl = vital_rates$M_Chk_survl,
       M_Fdg_survl = vital_rates$M_Fdg_survl,
       M_Adt_survl = vital_rates$M_Adt_survl)

# number of stages in the matrix
no_stages <- sqrt(length(matrix_str))

# Define plover life-stages of the Ceuta snowy plover matrix model
stages <- c("F_1st_yr",  "F_Adt",  "M_1st_yr",  "M_Adt")

# an empty t by x matrix
stage <- matrix(numeric(no_stages * niter), nrow = no_stages)

# an empty t vector to store the population sizes
pop <- numeric(niter)

# dataframe to store the perturbation results
ASR_pert_results <- 
  data.frame(parameter = c("F_Fdg_survl"),
             sensitivities = numeric(1),
             elasticities = numeric(1))

lambda_pert_results <-
  data.frame(parameter = c("F_Fdg_survl", "h", "HSR"),
             sensitivities = numeric(3),
             elasticities = numeric(3))

# specifiy how many survival rates there are
n <- length(vr)-5

# create a vector of perturbations to test on parameters of the matrix model, in
# this case this proportional changes in survival (i.e., between 0 an 1)
vr_nums <- seq(0, 1, 0.01) 
h_nums <- seq(0, 2, 0.02) # proportional changes in h index (i.e., between 0 and 2)

# create an empty dataframe to store the perturbation results
vr_pert_ASR <- matrix(numeric(n * length(vr_nums)), 
                  ncol = n, dimnames = list(vr_nums, "F_Fdg_survl"))
vr_pert_lambda <- matrix(numeric(n * length(vr_nums)), 
                      ncol = n, dimnames = list(vr_nums, "F_Fdg_survl"))
h_pert_lambda <- matrix(numeric(length(h_nums)),
                        ncol = 1, dimnames = list(h_nums, "h"))
HSR_pert_lambda <- matrix(numeric(length(vr_nums)),
                          ncol = 1, dimnames = list(vr_nums, "HSR"))

# for loop to simulate the ASR response to each perturbation in the vital rate
for (i in 1:length(vr_nums)) # pick a perturbation level
  {
    # specify the vital rate with the new perturbation level
    vr[["F_Fdg_survl"]] <- vr_nums[i] 
    # build the matrix with the new value
    A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), 
                byrow=TRUE, dimnames = list(stages, stages)) 
    # reset the starting stage distribution for simulation (all with 10 individuals)
    m <- rep(10, no_stages) 
    # project the matrix through 1000 iterations
    for (j in 1:niter) { 
      # stage distribution at time t
      stage[,j] <- m 
      # population size at time t
      pop[j] <- sum(m)
      # number of male adults at time t
      M2 <- stage[4, j] 
      # number of female adults at time t
      F2 <- stage[2, j] 
      # Female freq-dep fecundity of Female chicks
      A[1,no_stages/2] <- ((k*M2)/(M2+(F2/h)))*HSR 
      # Female freq-dep fecundity of Male chicks
      A[(no_stages/4)*3,no_stages/2]  <- ((k*M2)/(M2+(F2/h)))*HSR
      # Male freq-dep fecundity of Female chicks
      A[1,no_stages] <- ((k*F2)/(M2+(F2/h)))*HSR 
      # Male freq-dep fecundity of Male chicks
      A[(no_stages/4)*3,no_stages] <- ((k*F2)/(M2+(F2/h)))*HSR 
      # define the new n (i.e., new stage distribution at time t)
      m <- A %*% m 
    }
    # define rownames of stage matrix
    rownames(stage) <- rownames(A) 
    # define colnames of stage matrix
    colnames(stage) <- 0:(niter - 1)
    # calculate the proportional stable stage distribution
    stage <- apply(stage, 2, function(x) x/sum(x))
    # define stable stage as the last stage
    stable.stage <- stage[, niter]
    # calc ASR as the proportion of the adult stable stage class that is male
    vr_pert_ASR[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
    # calc lambda as the pop change in the counts of the last two iterations
    vr_pert_lambda[i, ] <- pop[niter]/pop[niter - 1]
}
# perturbation of the h index parameter
for (i in 1:length(h_nums)) # pick a perturbation level
  {
  A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), byrow=TRUE, dimnames = list(stages, stages)) # build the matrix with the new value
  m <- rep(10, no_stages) # reset the starting stage distribution for simulation (all with 10 individuals)
  for (j in 1:niter) { # project the matrix through t iteration
    # stage distribution at time t
    stage[,j] <- m
    # population size at time t
    pop[j] <- sum(m)
    # number of male adults at time t
    M2 <- stage[4, j]
    # number of female adults at time t
    F2 <- stage[2, j]
    # Female freq-dep fecundity of Female chicks
    A[1,no_stages/2]        <- ((k*M2)/(M2+(F2/h_nums[i])))*HSR
    # Female freq-dep fecundity of Male chicks
    A[(no_stages/4)*3,no_stages/2]  <- ((k*M2)/(M2+(F2/h_nums[i])))*HSR
    # Male freq-dep fecundity of Female chicks
    A[1,no_stages]          <- ((k*F2)/(M2+(F2/h_nums[i])))*HSR
    # Male freq-dep fecundity of Male chicks
    A[(no_stages/4)*3,no_stages]    <- ((k*F2)/(M2+(F2/h_nums[i])))*HSR
    # define the new n (i.e., new stage distribution at time t)
    m <- A %*% m
  }
  # define rownames of stage matrix
  rownames(stage) <- rownames(A)
  # define colnames of stage matrix
  colnames(stage) <- 0:(niter - 1)
  # calculate the proportional stable stage distribution
  stage <- apply(stage, 2, function(x) x/sum(x))
  # define stable stage as the last stage
  stable.stage <- stage[, niter]
  # calc lambda as the pop change in the counts of the last two iterations
  h_pert_lambda[i, ] <- pop[niter]/pop[niter - 1]
}
# perturbation of HSR
for (i in 1:length(vr_nums)) # pick a perturbation level
  {
  A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), byrow=TRUE, dimnames = list(stages, stages)) # build the matrix with the new value
  m <- rep(10, no_stages) # reset the starting stage distribution for simulation (all with 10 individuals)
  for (j in 1:niter) { # project the matrix through t iteration
    # stage distribution at time t
    stage[,j] <- m
    # population size at time t
    pop[j] <- sum(m)
    # number of male adults at time t
    M2 <- stage[4, j]
    # number of female adults at time t
    F2 <- stage[2, j]
    # Female freq-dep fecundity of Female chicks
    A[1,no_stages/2]        <- ((k*M2)/(M2+(F2/h)))*vr_nums[i]
    # Female freq-dep fecundity of Male chicks
    A[(no_stages/4)*3,no_stages/2]  <- ((k*M2)/(M2+(F2/h)))*vr_nums[i]
    # Male freq-dep fecundity of Female chicks
    A[1,no_stages]          <- ((k*F2)/(M2+(F2/h)))*vr_nums[i]
    # Male freq-dep fecundity of Male chicks
    A[(no_stages/4)*3,no_stages]    <- ((k*F2)/(M2+(F2/h)))*vr_nums[i]
    # define the new n (i.e., new stage distribution at time t)
    m <- A %*% m
  }
  # define rownames of stage matrix
  rownames(stage) <- rownames(A)
  # define colnames of stage matrix
  colnames(stage) <- 0:(niter - 1)
  # calculate the proportional stable stage distribution
  stage <- apply(stage, 2, function(x) x/sum(x))
  # define stable stage as the last stage
  stable.stage <- stage[, niter]
  # calc lambda as the pop change in the counts of the last two iterations
  HSR_pert_lambda[i, ] <- pop[niter]/pop[niter - 1]
}

# find the spline for the results
VR_Perturbations_FS <- as.numeric(rownames(vr_pert_ASR))
ASR_Results_FS <- vr_pert_ASR[,1]
spl_Fdg_survl <- smooth.spline(ASR_Results_FS ~ VR_Perturbations_FS)

VR_Perturbations_FS_lambda <- as.numeric(rownames(vr_pert_lambda))
lambda_Results_FS <- vr_pert_lambda[,1]
spl_Fdg_survl_lambda <- smooth.spline(lambda_Results_FS ~ VR_Perturbations_FS_lambda)

VR_Perturbations_h <- as.numeric(rownames(h_pert_lambda))
lambda_Results_h <- h_pert_lambda[,1]
spl_h <- smooth.spline(lambda_Results_h ~ VR_Perturbations_h)

VR_Perturbations_HSR <- as.numeric(rownames(HSR_pert_lambda))
lambda_Results_HSR <- HSR_pert_lambda[,1]
spl_HSR <- smooth.spline(lambda_Results_HSR ~ VR_Perturbations_HSR)

# Assign the emperical value of the vital rate (i.e., where to estimate the tangent line)
newx_FS <- vital_rates$F_Fdg_survl
newx_FS_lambda <- vital_rates$F_Fdg_survl
newx_h <- vital_rates$h
newx_HSR <- vital_rates$HSR

# The first order derivative
pred0_FS <- predict(spl_Fdg_survl, x=newx_FS, deriv=0)
pred0_FS_lambda <- predict(spl_Fdg_survl_lambda, x=newx_FS_lambda, deriv=0)
pred0_h <- predict(spl_h, x=newx_h, deriv=0)
pred0_HSR <- predict(spl_HSR, x=newx_HSR, deriv=0)

# The second order derivative
pred1_FS <- predict(spl_Fdg_survl, x=newx_FS, deriv=1)
pred1_FS_lambda <- predict(spl_Fdg_survl_lambda, x=newx_FS_lambda, deriv=1)
pred1_h <- predict(spl_h, x=newx_h, deriv=1)
pred1_HSR <- predict(spl_HSR, x=newx_HSR, deriv=1)

# Find the y-intercept of the tangent
yint_FS <- pred0_FS$y - (pred1_FS$y*newx_FS)
yint_FS_lambda <- pred0_FS_lambda$y - (pred1_FS_lambda$y*newx_FS_lambda)
yint_h <- pred0_h$y - (pred1_h$y*newx_h)
yint_HSR <- pred0_HSR$y - (pred1_HSR$y*newx_HSR)

# plot the perturbation results
plot(x = VR_Perturbations_FS, y = ASR_Results_FS, 
     main = "non-linear matrix sensitivity of female fledgling survival",
     xlab = "survival rate value",
     ylab = "adult sex ratio response")
lines(spl_Fdg_survl, col=2)
points(pred0_FS, col=2, pch=19)
lines(ASR_Results_FS, yint_FS + pred1_FS$y*ASR_Results_FS, col=3)

plot(x = VR_Perturbations_FS_lambda, y = lambda_Results_FS, 
     main = "non-linear matrix sensitivity of female fledgling survival",
     xlab = "survival rate value",
     ylab = "lambda response")
lines(spl_Fdg_survl_lambda, col=2)
points(pred0_FS_lambda, col=2, pch=19)
lines(lambda_Results_FS, yint_FS_lambda + pred1_FS_lambda$y*lambda_Results_FS, col=3)

plot(x = VR_Perturbations_h, y = lambda_Results_h, 
     main = "non-linear matrix sensitivity of h index",
     xlab = "h value",
     ylab = "lambda response")
lines(spl_h, col=2)
points(pred0_h, col=2, pch=19)
lines(lambda_Results_h, yint_h + pred1_h$y*lambda_Results_h, col=3)

plot(x = VR_Perturbations_HSR, y = lambda_Results_HSR, 
     main = "non-linear matrix sensitivity of HSR",
     xlab = "HSR value",
     ylab = "lambda response")
lines(spl_HSR, col=2)
points(pred0_HSR, col=2, pch=19)
lines(lambda_Results_HSR, yint_HSR + pred1_HSR$y*lambda_Results_HSR, col=3)

# assign the slope of the tangent line to the dataframe (this is the sensitivity)
ASR_pert_results[, 2] <- predict(spl_Fdg_survl, x=vr[["F_Fdg_survl"]], deriv=1)$y
lambda_pert_results[2, 2] <- predict(spl_h, x=h, deriv=1)$y
lambda_pert_results[3, 2] <- predict(spl_HSR, x=HSR, deriv=1)$y

# calculate the elasticity of the vital rate given the sensitivity
ASR_pert_results[, 3] <- vr[["F_Fdg_survl"]]/ASR * ASR_pert_results[, 2]
lambda_pert_results[, 3] <- h/lambda * lambda_pert_results[, 2]
lambda_pert_results[, 3] <- HSR/lambda * lambda_pert_results[, 2]

############## PVA simulation #################################################
boot_matricies <- bootstrap_matrix_extract(survival_rates = survival_rates_boot, niter = 1000)

stoch_projection_plover <- 
  function (matrices, n0, tmax = 50, nreps = 5000, Quasi_Ex = 0, no_stages = 4, treatment = TRUE, HSR, h = 1, k = 3) 
  {
    x <-  no_stages
    F_Adults <- matrix(numeric(nreps * tmax), nrow = nreps)
    F_Chicks <- matrix(numeric(nreps * tmax), nrow = nreps)
    M_Adults <- matrix(numeric(nreps * tmax), nrow = nreps)
    M_Chicks <- matrix(numeric(nreps * tmax), nrow = nreps)
    Adults <- matrix(numeric(nreps * tmax), nrow = nreps)
    F_Adults[,1] <- n0[2]
    F_Chicks[,1] <- n0[1]
    M_Adults[,1] <- n0[4]
    M_Chicks[,1] <- n0[3]
    Adults[,1] <- n0[2] + n0[4]
    for (i in 1:nreps) {
      if(treatment){
        A <- sample(matrices$treat_mats, (tmax-1), replace = TRUE)
      }
      else{
        A <- sample(matrices$contr_mats, (tmax-1), replace = TRUE)
      }
      for (j in 1:(tmax-1)) {
        n <- c(F_Chicks[i,j], F_Adults[i,j], M_Chicks[i,j], M_Adults[i,j])
        # number of male adults at time t
        M2 <- n[4] 
        # number of female adults at time t
        F2 <- n[2] 
        # Female freq-dep fecundity of Female chicks
        A[[j]][1,x/2]        <- (k*M2)/(M2+(F2/h))*HSR 
        # Female freq-dep fecundity of Male chicks
        A[[j]][(x/4)*3,x/2]  <- (k*M2)/(M2+(F2/h))*HSR
        # Male freq-dep fecundity of Female chicks
        A[[j]][1,x]          <- (k*F2)/(M2+(F2/h))*HSR 
        # Male freq-dep fecundity of Male chicks
        A[[j]][(x/4)*3,x]    <- (k*F2)/(M2+(F2/h))*HSR 
        # define the new n (i.e., new stage distribution at time t)
        B <- A[[j]] %*% n
        F_Chicks[i,(j+1)] <- B[1]
        F_Adults[i,(j+1)] <- B[2]
        M_Chicks[i,(j+1)] <- B[3]
        M_Adults[i,(j+1)] <- B[4]
        Adults[i,(j+1)] <- B[2] + B[4]
      }
    }
    Adults_CV <- apply(Adults, 1, sd)/apply(Adults, 1, mean)
    Adults_Final <- Adults[,tmax]
    Adults_melted <- melt(t(Adults))
    Adults_ExProb <- sum(Adults_Final < Quasi_Ex)/nreps
    colnames(Adults_melted) <- c("Year", "Iteration", "Adults")
    Sim_Avg <- Adults_melted %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(mean(Adults))
    Adults <- list(Adults_CV = Adults_CV, Adults_ExProb = Adults_ExProb, 
                   Adults_Final = Adults_Final, Adults_melted = Adults_melted,
                   Sim_Avg = Sim_Avg)
    Adults
  }

testn <- c(100, 50, 100, 50)
treat_stoch_proj <- stoch_projection_plover(matrices = boot_matricies, n0 = testn, 
                                            tmax = 25, nreps = 1000, Quasi_Ex = 2, HSR = HSR)

contr_stoch_proj <- stoch_projection_plover(matrices = boot_matricies, n0 = testn, 
                                            tmax = 25, nreps = 1000, Quasi_Ex = 2, HSR = HSR,
                                            treatment = FALSE)

# Line plots of simulation trajectories with trend lines illustrating
# the predicted population dynamics post-2012 for 10, 25, and 50 years into
# the future
Stoch_Proj_plot <- {
  ggplot(NULL) +  
    theme_bw() +
    geom_line(data = treat_stoch_proj$Adults_melted, 
              aes(x = (Year-1), y = Adults, group = Iteration), 
              size = 0.25,  alpha = 0.02) +
    geom_line(data = treat_stoch_proj$Sim_Avg, 
              aes(x = (Year-1), y = `mean(Adults)`), 
              size = 1,  alpha = 1) +
    geom_line(data = contr_stoch_proj$Adults_melted, 
              aes(x = (Year-1), y = Adults, group = Iteration), 
              size = 0.25,  alpha = 0.02, color = "red") +
    geom_line(data = contr_stoch_proj$Sim_Avg, 
              aes(x = (Year-1), y = `mean(Adults)`), 
              size = 1,  alpha = 1, color = "red") +
    scale_colour_manual("", values = c("black", "black")) +
    theme(text=element_text(family="Arial"),
          axis.title.x = element_text(size=8, vjust=-0.1),
          axis.text.x  = element_text(size=7), 
          axis.title.y = element_text(size=8, vjust=1.2),
          axis.text.y  = element_text(size=7), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_line(size = 0.2, colour = "grey40"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks.y = element_line(size = 0.2, colour = "grey40")) +
    ylab("Adult population size") +
    xlab("Annual iterations") +
    # annotate("text", x = 10, y = 60, vjust = 1, hjust = 1, 
    #          label = "P(10 year extinction) = 0.166", size = 3.5, family="Arial") +
    scale_y_continuous(limits = c(0, 100)) +
    scale_x_continuous(limits = c(0, 25), breaks = 0:25)
}
Stoch_Proj_plot