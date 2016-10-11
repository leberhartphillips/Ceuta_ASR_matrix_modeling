library(ggplot2)

# load the following dataframe from the github repo
survival_rates_boot <- 
  read.table("output/bootstrap/survival_rates_boot_out.txt", header = TRUE)

# number of iterations used in bootstrap
niter <- 1000 
# Hatching sex ratio (proportion of hatchlings that are male)
HSR <- 0.4856322
h <- 1
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
                     stable.stage = stable.stage, 
                     stage.vectors = stage,
                     SSD_M2 = stable.stage[4],
                     SSD_F2 = stable.stage[2])
    # print the list as output to the function
    pop.proj 
  }

# this function extracts the ASR estimate from each bootstrap iteration using
# a matrix with the tow-sex mating function
ASR_extract <- 
  function(survival_rates, niter = 1000){
    # make an empty datarame to store the results
    ASR_output <- data.frame(ASR_output = numeric(niter))
    # for loop to go through each iteration and calculate the differece between female and male
    # survival rates for each stage.
    for(i in 1:niter){
      # Create a list of demographic rates from the survival analyses above
      vital_rates_boot <- list(F_Chk_survl = survival_rates[which(survival_rates$iter == i), 2][5],
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
      
      # Build matrix based on rates specified in the list above
      matrix_boot <- plover_matrix(vital_rates_boot)

      # Determine the ASR at the stable stage distribution
      treatment_ASR_analysis <- 
        matrix_ASR(A = matrix_boot, h = vital_rates_boot$h, HSR = vital_rates_boot$HSR, iterations = 75)
      
      # Extract ASR
      ASR_output[i, 1] <- treatment_ASR_analysis$ASR
    }
    # restructure the output and lable columns
    ASR_output <- suppressMessages(reshape2::melt(data = ASR_output))[-1]
    colnames(ASR_output) <- c("estimate")
    # return the output
    ASR_output
  }

# this function runs a perturbation analysis on the two-sex non-linear matrix to 
# calculate the sensitivites and elasticities of each parameter included in the 
# matrix (including the components of the mating function)
sensitivity_analysis <- 
  function(vital_rates, matrix_str, h = 1, k = 3, HSR, niter = 1000, ASR){
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
    
    # dataframe to store the perturbation results
    perturbation_results <- 
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
    vr_pert <- matrix(numeric(n * length(vr_nums)), 
                      ncol = n, dimnames = list(vr_nums, names(vr)))
    h_pert <- matrix(numeric(length(h_nums)), 
                     ncol = 1, dimnames = list(h_nums, "h"))
    k_pert <- matrix(numeric(length(k_nums)), 
                     ncol = 1, dimnames = list(k_nums, "k"))
    HSR_pert <- matrix(numeric(length(vr_nums)), 
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
        # calculate the proportional stable stage distribution
        #stable.stage <- w/sum(w)
        # calc ASR as the proportion of the adult stable stage class that is male
        vr_pert[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
      }
      spl <- smooth.spline(vr_pert[,g] ~ rownames(vr_pert))
      perturbation_results[g, 2] <- predict(spl, x=vr[[g]], deriv=1)$y
      perturbation_results[g, 3] <- vr[[g]]/ASR * perturbation_results[g, 2]
    }
    # perturbation of the h index parameter
    for (i in 1:length(h_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), byrow=TRUE, dimnames = list(stages, stages)) # build the matrix with the new value
      m <- rep(10, no_stages) # reset the starting stage distribution for simulation (all with 10 individuals)
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m 
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
      h_pert[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
    }
    spl <- smooth.spline(h_pert[, 1] ~ rownames(h_pert))
    perturbation_results[n+1, 2] <- predict(spl, x=h, deriv=1)$y
    perturbation_results[n+1, 3] <- h/ASR * perturbation_results[n+1, 2]
    
    # perturbation of k parameter
    for (i in 1:length(k_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), byrow=TRUE, dimnames = list(stages, stages)) # build the matrix with the new value
      m <- rep(10, no_stages) # reset the starting stage distribution for simulation (all with 10 individuals)
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m 
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
      k_pert[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
    }
    spl <- smooth.spline(k_pert[,1] ~ rownames(k_pert))
    perturbation_results[n+2, 2] <- predict(spl, x=k, deriv=1)$y
    perturbation_results[n+2, 3] <- k/ASR * perturbation_results[n+2, 2]
    # perturbation of HSR
    for (i in 1:length(vr_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), byrow=TRUE, dimnames = list(stages, stages)) # build the matrix with the new value
      m <- rep(10, no_stages) # reset the starting stage distribution for simulation (all with 10 individuals)
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m 
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
      HSR_pert[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
    }
    spl <- smooth.spline(HSR_pert[,1] ~ rownames(HSR_pert))
    perturbation_results[n+3, 2] <- predict(spl, x=HSR, deriv=1)$y
    perturbation_results[n+3, 3] <- HSR/ASR * perturbation_results[n+3, 2]
    perturbation_results
  }

# this function estimates the contribution that each vital rate has on ASR bias,
# given the sensitivities calculated in the previous function
# (see formula 8 on page 133 of Veran and Beissinger (2009))
LTRE_analysis <- 
  function(treatment_sens, Mprime_sens, matrix_str, vital_rates){
    # make an empty dataframe to stroe LTRE results
    LTRE_results <- 
      data.frame(parameter = c("Chick survival", "Fledgling survival", 
                               "Adult survival", "Hatching sex ratio"),
                 contribution = numeric(4))
    
    # run a for loop to extract the parameter contributions
    for(i in 1:nrow(LTRE_results))
    {
      LTRE_results[i, 2] <- 
        ifelse(i < 4, (vital_rates[[i + 3]] - vital_rates[[i]]) * Mprime_sens$sensitivities[i + 3],
               (vital_rates[[i]] - 0.5) * Mprime_sens$sensitivities[9])
    }
    LTRE_results$parameter <- factor(LTRE_results$parameter, levels = c("Adult survival",
                                                                        "Fledgling survival",
                                                                        "Chick survival",
                                                                        "Hatching sex ratio"))
    LTRE_results
  }

###############################################################################
####### ASR estimate with two-sex matrix that includes mating function ########

# extract the ASR estimate for each iteration of the bootstrapped survival
# analysis
ASR_boot <- ASR_extract(survival_rates_boot)

# Calculate the mean, and median of the ASR bootstraps
ASR_boot_summary <-
  Rmisc::summarySE(ASR_boot,
                   measurevar = "estimate")

# specify the confidence interval to assess
CI <- 0.95
Ceuta_ASR_95CI_quan <- stats::quantile(ASR_boot, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
ASR_boot_summary <- as.data.frame(cbind(ASR_boot_summary, Ceuta_ASR_95CI_quan[1], Ceuta_ASR_95CI_quan[2]))
rownames(ASR_boot_summary) <- NULL
colnames(ASR_boot_summary) <- c(".id","N", "estimate", "sd", "se", "ci", "lcl", "ucl")

# plot the histogram distribution of ASRs calculated from each bootstrap iteration
ASR_bootstrap_histogram <- 
  ggplot() +
  annotate("rect", xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(2)]) +
  annotate("rect", xmin=0.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(1)]) +
  annotate("text", x = c(-Inf,Inf), y = c(70, 70),
           label = c("\u2640", "\u2642"), size = 7,
           family="Arial", vjust = c(1.5,1.5), hjust = c(-0.5,1.5)) +
  geom_histogram(binwidth = 0.02, data = ASR_boot, aes(x = estimate)) +
  geom_errorbarh(data = ASR_boot_summary, aes(y = 125, x = lcl, xmin = lcl, xmax = ucl), color = "black", size = 0.8, linetype = "solid") +
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
  scale_y_continuous(limits = c(0, 125))
ASR_bootstrap_histogram

# ggsave(ASR_bootstrap_histogram,
#        filename = "ASR_distribution_final.jpg",
#        path = "figs/final/",
#        width = 4,
#        height = 3, units = "in",
#        dpi = 300,
#        scale = 1)

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
                  h = h,
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
  matrix_ASR(A = treatment_matrix, h = 1, HSR = VR_treat$HSR, iterations = 1000)
M_prime_ASR_analysis <- 
  matrix_ASR(A = M_prime_matrix, h = 1, HSR = VR_mprime$HSR, iterations = 1000)
ASR_treat <- 
  treatment_ASR_analysis$ASR
ASR_mprime <- 
  M_prime_ASR_analysis$ASR

# conduct a sensitivity analysis on the treatment matrix
treat_sensitivity_analysis <- 
  sensitivity_analysis(vital_rates = VR_treat, 
                       matrix_str = matrix_str, 
                       h = VR_treat$h, 
                       k = VR_treat$k, 
                       HSR = VR_treat$HSR, 
                       niter = 1000, 
                       ASR = ASR_treat)

# conduct a sensitivity analysis on the M-Prime matrix
Mprime_sensitivity_analysis <- 
  sensitivity_analysis(vital_rates = VR_mprime, 
                       matrix_str = matrix_str, 
                       h = VR_treat$h, 
                       k = VR_treat$k, 
                       HSR = VR_mprime$HSR, 
                       niter = 1000, 
                       ASR = ASR_mprime)

# conduct the LTRE comparing the two matrices
LTRE_plover <- 
  LTRE_analysis(treatment_sens = treat_sensitivity_analysis, 
                Mprime_sens = Mprime_sensitivity_analysis, 
                matrix_str = matrix_str, 
                vital_rates = VR_treat)

# specify the color palette to use in LTRE
cbPalette <- c("#737373", "#BDBDBD", "#BDBDBD", "#BDBDBD")

# plot the LTRE results
LTRE_plot <- 
ggplot2::ggplot() +
  theme_bw() +
  coord_flip() +
  geom_bar(data = LTRE_plover,
           aes(x = parameter, y = contribution, fill = parameter), 
           color = "black", stat = "identity", alpha = 0.8) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.x = element_text(size=12, margin = margin(10, 0, 0, 0)),
        axis.text.x  = element_text(size=10, margin = margin(5, 0, 0, 0)), 
        axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
        axis.text.y  = element_text(size=10, angle = 90, hjust = 0.5, 
                                    margin = margin(0, 1, 0, 0)),
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
  scale_y_continuous(limits = c(-0.10, 0.1))
LTRE_plot

# check to see if the sum of LTRE contributions equals the difference
# between the ASR of the treatment matrix and the M-prime matrix
sum(LTRE_plover$contribution) # ---> 0.1131799

ASR_treat-ASR_mprime # ---> 0.07136055

# check the amount of the difference explained by M-prime
ASR_mprime/ASR_treat # ---> 0.89%

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

# dataframe to store the perturbation results
perturbation_results <- 
  data.frame(parameter = c("F_Adt_survl"),
             sensitivities = numeric(1),
             elasticities = numeric(1))

# specifiy how many survival rates there are
n <- length(vr)-5

# create a vector of perturbations to test on parameters of the matrix model, in
# this case this proportional changes in survival (i.e., between 0 an 1)
vr_nums <- seq(0, 1, 0.01) 

# create an empty dataframe to store the perturbation results
vr_pert <- matrix(numeric(n * length(vr_nums)), 
                  ncol = n, dimnames = list(vr_nums, "F_Adt_survl"))

# for loop to simulate the ASR response to each perturbation in the vital rate
for (i in 1:length(vr_nums)) # pick a perturbation level
  {
    # specify the vital rate with the new perturbation level
    vr[["F_Adt_survl"]] <- vr_nums[i] 
    
    # build the matrix with the new value
    A <- matrix(sapply(matrix_str, eval, vr, NULL), nrow = sqrt(length(matrix_str)), 
                byrow=TRUE, dimnames = list(stages, stages)) 
    # reset the starting stage distribution for simulation (all with 10 individuals)
    m <- rep(10, no_stages) 
    
    # project the matrix through 1000 iterations
    for (j in 1:niter) { 
      
      # stage distribution at time t
      stage[,j] <- m 
      
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
    vr_pert[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
  }

# plot the perturbation results
VR_Perturbations <- rownames(vr_pert)
ASR_Results <- vr_pert
plot(x = VR_Perturbations, y = ASR_Results, 
     main = "non-linear matrix sensitivity of female adult survival",
     xlab = "vital rate value",
     ylab = "adult sex ratio response")

# find the spline for the results
spl <- smooth.spline(ASR_Results ~ VR_Perturbations)
lines(spl, col=2)

# Assign the emperical value of the vital rate (i.e., where to estimate the tangent line)
newx <- vital_rates$F_Adt_survl

# The first order derivative
pred0 <- predict(spl, x=newx, deriv=0)

# The second order derivative
pred1 <- predict(spl, x=newx, deriv=1)

# Plot the point where to estimate the tangent
points(pred0, col=2, pch=19)

# Find the y-intercept of the tangent
yint <- pred0$y - (pred1$y*newx)

# plot the tangent
lines(ASR_Results, yint + pred1$y*ASR_Results, col=3)

# assign the slope of the tangent line to the dataframe (this is the sensitivity)
perturbation_results[, 2] <- predict(spl, x=vr[["F_Adt_survl"]], deriv=1)$y

# calculate the elasticity of the vital rate given the sensitivity
perturbation_results[, 3] <- vr[["F_Adt_survl"]]/ASR * perturbation_results[, 2]