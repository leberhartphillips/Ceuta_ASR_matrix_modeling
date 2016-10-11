library(RMark) 
library(stringr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(reshape2)
library(RColorBrewer)
library(Rmisc)
library(stats)
library(lme4)
library(magrittr)
setwd("~/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling")
chick <- 
  read.table("data/chick_mark-recapture_data.txt",
             header = TRUE, colClasses = c("factor", "character","factor",
                                           "numeric","factor","factor", "numeric"))
fledgling_adult <- 
  read.table("data/fledgling_adult_mark-recapture_data.txt",
             header = TRUE, colClasses = c("factor","character","factor","factor"))
breeding_data <- 
  read.table("data/breeding_data.txt", 
             header = TRUE)

# number of iterations used in bootstrap
niter <- 1000 
# Hatching sex ratio (proportion of hatchlings that are male)
HSR <- 0.4856322 
# Sex-specific annual per captia fecundity rates (hatchlings)
# Females:
RF <- 2.061886
Sex <- rep("Female", nrow(breeding_data))
Ring <- breeding_data$female
females <- data.frame(Ring, Sex)
Sex <- rep("Male", nrow(breeding_data))
Ring <- breeding_data$male
males <- data.frame(Ring, Sex)
Individuals <- rbind(males, females)
reproduction_df <- cbind(breeding_data[rep(row.names(breeding_data), 2), 
                                       c("no_chicks", "clutch_size", "brood_ID", "year")],
                         Individuals)
reproduction_df$Sex <- factor(reproduction_df$Sex, levels = c("Female", "Male"))
reproduction_df <- reproduction_df[!is.na(reproduction_df$Ring),]
reproduction_df <- reproduction_df[!is.na(reproduction_df$no_chicks),]
reproduction_df <- dplyr::group_by(reproduction_df, year, Sex, Ring)
reproduction_df_sum <- 
  dplyr::ungroup(dplyr::summarise(reproduction_df, 
                                  total_chicks_p_year = sum(as.numeric(no_chicks))))
fecundity_annual_summary <- 
  Rmisc::summarySE(reproduction_df_sum, measurevar = "total_chicks_p_year", 
                   groupvars = c("Sex", "year"))
reproduction_df_sum <- dplyr::group_by(reproduction_df_sum, Sex, Ring)
reproduction_df_sum_avg <- 
  dplyr::ungroup(dplyr::summarise(reproduction_df_sum, 
                                  avg_chicks_p_year = mean(as.numeric(total_chicks_p_year))))
fecundity_sex_summary <- 
  Rmisc::summarySE(fecundity_annual_summary, 
                   measurevar = "total_chicks_p_year", groupvars = c("Sex"))
sample_sizes_sex <- 
  aggregate(Ring ~ Sex, data = reproduction_df_sum_avg, FUN = function(x){NROW(x)})
cbPalette <- RColorBrewer::brewer.pal(8, "Dark2")[c(2,1)]
Sex_specific_fucund_plot <- 
ggplot2::ggplot() +
  geom_boxplot(aes(y = avg_chicks_p_year, x = Sex, fill = Sex), 
               data = reproduction_df_sum_avg, size = .3, alpha = 0.6) +
  geom_text(data = sample_sizes_sex, size = 3, 
            aes(y = c(6.5, 6.5), x = c(1.12, 2.12), label = Ring)) +
  annotate("text", x = c(0.92, 1.92), y = c(6.5, 6.5), label = "n = ", 
           size = 3) +
  theme_bw() +
  theme(text=element_text(size=16),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 10), 
        axis.title.y = element_text(size = 12, 
                                    margin = margin(0, 15, 0, 0)),
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
# ggsave(Sex_specific_fucund_plot,
#        filename = "Sex-specific_reproductive_success_final.jpg",
#        path = "figs/final/",
#        width = 2.8,
#        height = 5, units = "in",
#        dpi = 300,
#        scale = 1)
###############################################################################

plover_matrix <- 
  function(demographic_rates)
  {
    # Define plover life-stages of the Ceuta snowy plover matrix model
    stages <- c("F_1st_yr",  "F_Adt",  "M_1st_yr",  "M_Adt")
    # Build the 4x4 matrix
    result <- 
      matrix(c(
        # top row of matrix
        0, (demographic_rates$RF * (1 - HSR))/2, 0, (demographic_rates$RF * (1 - HSR))/2, 
        # second row of matrix
        (demographic_rates$F_Chk_survl*demographic_rates$F_Fdg_survl),
        demographic_rates$F_Adt_survl, 
        0, 0,
        # third row of matrix
        0, (demographic_rates$RF * HSR)/2, 0, (demographic_rates$RF * HSR)/2,
        # fourth row of matrix
        0, 0, 
        (demographic_rates$M_Chk_survl*demographic_rates$M_Fdg_survl),
        demographic_rates$M_Adt_survl),
        nrow = 4, byrow = TRUE,
        dimnames = list(stages, stages))
    result
  }

matrix_ASR <- 
  function (A) 
  {
    # Number of stages in the matrix
    x <- ncol(A) 
    {
      # determine the dominant eigen value and vector of the matrix
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

ASR_analysis <- 
  function (A, zero = TRUE) 
  {
    
    # makes list of the eigen values and eigen vectors of A
    ev <- eigen(A) 
    
    # index of dominant eigen value
    lmax <- which.max(Re(ev$values)) 
    
    # Eigen vectors
    W <- ev$vectors 
    
    # dominant eigen vector
    w <- abs(Re(W[, lmax])) 
    
    # stable stage distribution
    stable.stage = w / sum(w) 
    
    # ASR
    ASR <- stable.stage[4] / (stable.stage[2] + stable.stage[4]) 
    
    # check if possible to proceed
    V <- try(Conj(solve(W)), silent = TRUE)
    if (class(V) == "try-error") {
      ASR.analysis <- list(ASR = ASR, stable.stage = stable.stage, 
                           sensitivities = A * NA, elasticities = A * NA)
    }
    else {
      
      # solve matrix
      v <- abs(Re(V[lmax, ])) 
      
      # outer product of v and w
      s <- v %o% w 
      if (zero) {
        s[A == 0] <- 0
      }
      
      # calculate elasticities
      e <- s * A/ASR
      
      # get vital rate names
      x <- dimnames(A) 
      
      # assign vital rate names to s
      dimnames(s) <- x 
      names(w) <- x[[1]]
      names(v) <- x[[1]]
      
      # output a list containing the ASR, SSD, sensitivites, and elasticities
      ASR.analysis <- list(ASR = ASR, stable.stage = stable.stage, 
                           sensitivities = s, elasticities = e)
    }
    ASR.analysis
  }

ASR_perturbation <- 
  function (elements, VR_list, freq_dep_ASR) 
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
      HSR = 0.5)
    
    # list of parameters in the M-prime matrix
    # this contains the average difference between each parameter
    M_prime_matrix <- list(
      F_Chk_survl = (treatment_matrix$F_Chk_survl + control_matrix$F_Chk_survl)/2,
      F_Fdg_survl = (treatment_matrix$F_Fdg_survl + control_matrix$F_Fdg_survl)/2,
      F_Adt_survl = (treatment_matrix$F_Adt_survl + control_matrix$F_Adt_survl)/2,
      M_Chk_survl = (treatment_matrix$M_Chk_survl + control_matrix$M_Chk_survl)/2,
      M_Fdg_survl = (treatment_matrix$M_Fdg_survl + control_matrix$M_Fdg_survl)/2,
      M_Adt_survl = (treatment_matrix$M_Adt_survl + control_matrix$M_Adt_survl)/2,
      RF = (treatment_matrix$RF + control_matrix$RF)/2,
      HSR = (treatment_matrix$HSR + control_matrix$HSR)/2)
    
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
    
    # add lower-level functions to the matrix
    vrs_treatment <- try(sapply(elements, eval, treatment_matrix, NULL), silent = TRUE)
    vrs_control <- try(sapply(elements, eval, control_matrix, NULL), silent = TRUE)
    vrs_M_prime <- try(sapply(elements, eval, M_prime_matrix, NULL), silent = TRUE)
    
    
    # check if its okay to proceed
    if (class(vrs_treatment) == "try-error") {
      vrs_treatment <- sub("Error in eval\\(expr, envir, enclos\\) :",
                           "", vrs_treatment[1])
      stop(paste("Cannot evaluate element expression using given vital rates:",
                 vrs_treatment))
    }
    if (class(vrs_control) == "try-error") {
      vrs_control <- sub("Error in eval\\(expr, envir, enclos\\) :",
                         "", vrs_control[1])
      stop(paste("Cannot evaluate element expression using given vital rates:",
                 vrs_control))
    }
    if (class(vrs_M_prime) == "try-error") {
      vrs_M_prime <- sub("Error in eval\\(expr, envir, enclos\\) :",
                         "", vrs_M_prime[1])
      stop(paste("Cannot evaluate element expression using given vital rates:",
                 vrs_M_prime))
    }
    
    # make an empty dataframe where all the perturbation stats will go
    res <- data.frame(estimate = unlist(treatment_matrix), sensitivity = 0, 
                      elasticity = 0, LTRE = 0)
    
    # build the treatment matrix
    A_treatment <- matrix(vrs_treatment, nrow = n, byrow = TRUE)
    
    # build the control matrix
    A_control <- matrix(vrs_control, nrow = n, byrow = TRUE)
    
    # build the control matrix
    M_prime <- matrix(vrs_M_prime, nrow = n, byrow = TRUE)
    
    # run sensitivity analyses on both matrices
    ASR_M_prime <- ASR_analysis(M_prime)
    ASR_treatment <- ASR_analysis(A_treatment)
    
    # calculate derivatives of lower-level matrix elements
    deriv.funcs <- sapply(elements, deriv, namevec = names(treatment_matrix), 
                          function.arg = TRUE)
    devs_treatment <- lapply(deriv.funcs, function(x) do.call(x, treatment_matrix))
    devs_m_prime <- lapply(deriv.funcs, function(x) do.call(x, M_prime_matrix))
    
    # run for loop to go through each parameter and estimate elasticity,
    # sensitivity, and LTRE
    for (i in 1:length(treatment_matrix)) {
      # first extract the 
      derivs_treatment <- 
        matrix(as.numeric(lapply(devs_treatment, function(x) attr(x, "gradient")[i])), 
               nrow = n, byrow = TRUE)
      derivs_m_prime <- 
        matrix(as.numeric(lapply(devs_m_prime, function(x) attr(x, "gradient")[i])), 
               nrow = n, byrow = TRUE)
      
      res[i, 2] <- 
        sum(derivs_treatment * ASR_treatment$sensitivities)
      res[i, 3] <- 
        treatment_matrix[[i]] / ASR_treatment$ASR * 
        sum(derivs_treatment * ASR_treatment$sensitivities)
      
      # only do LTRE calculations on survival and HSR parameters RELATIVE to one sex
      # i.e., don't calculate LTRE on fecundity
      res[i, 4] <- ifelse(i > 3 & i < 6, NA,
                          # calculate the sex differences of each survival rate, then multiply it
                          # by the sensitivity of that parameter in the M prime matrix
                          ifelse(i < 4, (treatment_matrix[[i + 3]] - treatment_matrix[[i]]) * 
                                   sum(derivs_m_prime * ASR_M_prime$sensitivities),
                                 # calculate the difference in HSR from 0.5 and multiply it by its
                                 # sensitivity
                                 ifelse(i == 8, (control_matrix[[i]] - treatment_matrix[[i]]) * 
                                          sum(derivs_m_prime * ASR_M_prime$sensitivities),
                                        NA)))
    }
    
    # consolidate results
    y <- res
    y$Vital_rate <- as.factor(rownames(y))
    colnames(y) <- c("Estimate", "Sensitivity", "Elasticity", "LTRE", "Vital_rate")
    y_melt <- suppressMessages(reshape2::melt(y[,c(2:5)]))
    y_melt$parameter <- 
      as.factor(ifelse(stringr::str_detect(y_melt$Vital_rate,"Chk"), "Chick survival",
                       ifelse(stringr::str_detect(y_melt$Vital_rate,"Fdg"), "Fledgling survival",
                              ifelse(stringr::str_detect(y_melt$Vital_rate,"Adt"), "Adult survival",
                                     "Hatching sex ratio"))))
    y_melt$parameter <- factor(y_melt$parameter, levels = c("Adult survival",
                                                            "Fledgling survival",
                                                            "Chick survival",
                                                            "Hatching sex ratio"))
    y_melt$Sex <- as.factor(ifelse(stringr::str_detect(y_melt$Vital_rate,"F_") & 
                                     y_melt$variable != "LTRE", "Female", 
                                   ifelse(stringr::str_detect(y_melt$Vital_rate,"M_") & 
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

chick_AIC_tables <- 
  read.table("output/bootstrap/AIC_table_chick_boot_out.txt", header = TRUE)

fldg_ad_AIC_tables <- 
  read.table("output/bootstrap/AIC_table_fledgling_adult_boot_out.txt", header = TRUE)

survival_rates_boot <- 
  read.table("output/bootstrap/survival_rates_boot_out.txt", header = TRUE)

# ASR_boot <- 
#   read.table("output/bootstrap/ASR_boot_out.txt", header = TRUE)

ASR_extract <- 
  function(survival_rates) 
  {
    # make an empty datarame to store the results
    ASR_output <- data.frame(RF = numeric(niter),
                             RFM = numeric(niter))
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
                                RF = RF)
      
      # Build matrix based on rates specified in the list above
      demographic_matrix_RF <- plover_matrix(demographic_rates)
      demographic_matrix_RFM <- plover_matrix(demographic_rates)
      
      # Determine the ASR at the stable stage distribution
      ASR_SSD_RF <- matrix_ASR(A = demographic_matrix_RF)
      ASR_SSD_RFM <- matrix_ASR(A = demographic_matrix_RFM)
      
      # Extract ASR
      ASR_output[i, 1] <- ASR_SSD_RF$ASR
      ASR_output[i, 2] <- ASR_SSD_RFM$ASR
    }
    # restructure the output and lable columns
    ASR_output <- suppressMessages(reshape2::melt(data = ASR_output))
    colnames(ASR_output) <- c("ASR_model", "ASR_boot")
    # return the output
    ASR_output
  }

###############################################################################
### Survival sex differences plot #############################################
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
  sex_diff_surv_output <- reshape2::melt(data = sex_diff_surv_output)
  colnames(sex_diff_surv_output) <- c("stage", "difference")
  # return the output
  sex_diff_surv_output
}
sex_diff_survival_output <- sex_diff_survival(survival_rates_boot)
sex_diff_survival_summary <- 
  sex_diff_survival_output %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(avg = mean(difference),
                   median = median(difference),
                   var = var(difference))
cbPalette <- c("#A6A6A6", "#D9D9D9", "#D9D9D9")
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

# jpeg(filename = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/figs/final/Survival_sex_differences_final.jpeg",
#      quality = 100,
#      width = 4,
#      height = 4,
#      units = "in",
#      res = 300)

grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) )
print( Background + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( Bootstrap_sex_diff_VR_plot + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
# dev.off()

###############################################################################
################ ASR Bootstrap distribution plot ##############################
ASR_boot <- filter(ASR_extract(survival_rates_boot), ASR_model == "RFM")


CI <- 0.95
ASR_boot_95CI <- stats::quantile(ASR_boot$ASR_boot, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
ASR_boot_mean <- mean(ASR_boot$ASR_boot)
ASR_boot_median <- median(ASR_boot$ASR_boot)
ASR_boot_summary <- as.data.frame(cbind(ASR_boot_95CI[1], ASR_boot_95CI[2], 
                                        ASR_boot_mean, ASR_boot_median))
rownames(ASR_boot_summary) <- NULL
colnames(ASR_boot_summary) <- c("lcl", "ucl", "mean", "median")
ASR_bootstrap_histogram <- 
  ggplot() +
  annotate("rect", xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(2)]) +
  annotate("rect", xmin=0.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.6,
           fill=brewer.pal(8, "Dark2")[c(1)]) +
  annotate("text", x = c(-Inf,Inf), y = c(55, 55),
           label = c("\u2640", "\u2642"), size = 7,
           family="Arial", vjust = c(1.5,1.5), hjust = c(-0.5,1.5)) +
  #geom_histogram(binwidth = 0.02, data = filter(ASR_boot, ASR_model == "R"), aes(x = estimate)) +
  geom_histogram(binwidth = 0.02, data = ASR_boot, aes(x = ASR_boot)) +
  geom_errorbarh(data = ASR_boot_summary, aes(y = 85, x = lcl, xmin = lcl, xmax = ucl), color = "black", size = 0.8, linetype = "solid") +
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
  scale_y_continuous(limits = c(0, 110))
ASR_bootstrap_histogram

# ggsave(ASR_bootstrap_histogram,
#        filename = "ASR_distribution_final.jpg",
#        path = "figs/final/",
#        width = 4,
#        height = 3, units = "in",
#        dpi = 300,
#        scale = 1)

###############################################################################
########### AIC table plots ###################################################
# define the model number
chick_AIC_tables$model_number <- as.numeric(chick_AIC_tables$model)
fldg_ad_AIC_tables$model_number <- as.numeric(fldg_ad_AIC_tables$model)
# summarize the average AIC stats for each candidate model across all 1000 iterations
chick_AIC_tables_summary <- 
  chick_AIC_tables %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(avg_Delta = mean(DeltaAICc),
                   IQR_Delta = IQR(DeltaAICc),
                   avg_Weight = mean(weight),
                   IQR_Weight = IQR(weight))
fldg_ad_AIC_tables_summary <- 
  fldg_ad_AIC_tables %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(avg_Delta = mean(DeltaAICc),
                   IQR_Delta = IQR(DeltaAICc),
                   avg_Weight = mean(weight),
                   IQR_Weight = IQR(weight))
# rank the output by delta AIC and determine model number
chick_AIC_tables_summary <- dplyr::arrange(chick_AIC_tables_summary, avg_Delta)
chick_AIC_tables_summary$model_number <- as.numeric(chick_AIC_tables_summary$model)
fldg_ad_AIC_tables_summary <- dplyr::arrange(fldg_ad_AIC_tables_summary, avg_Delta)
fldg_ad_AIC_tables_summary$model_number <- as.numeric(fldg_ad_AIC_tables_summary$model)
# merge the two datasets for plotting
chick_AIC_tables <- 
  dplyr::left_join(chick_AIC_tables_summary, chick_AIC_tables, by = "model_number")
fldg_ad_AIC_tables <- 
  dplyr::left_join(fldg_ad_AIC_tables_summary, fldg_ad_AIC_tables, by = "model_number")
# extract the model structure explaining resighting probability
chick_AIC_tables$p <- 
  factor(chick_AIC_tables$p, 
         levels = str_sub(as.character(chick_AIC_tables_summary$model), 
                          start = 24, end = str_length(chick_AIC_tables_summary$model)-1))
fldg_ad_AIC_tables$p <- 
  factor(fldg_ad_AIC_tables$p,
         levels = str_sub(as.character(fldg_ad_AIC_tables_summary$model), 
                          start = 18, end = str_length(fldg_ad_AIC_tables_summary$model)-1))

Bootstrap_Delta_AIC_plot_C <- 
  ggplot(aes(y = DeltaAICc, x = p), data = chick_AIC_tables) + 
  theme_bw() +
  #geom_violin(fill = "grey40") +
  geom_boxplot(width = 0.3, fill = "grey70", outlier.size = 0.5) +
  theme(text = element_text(family="Arial"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=12, margin = margin(0, 18, 0, 0)),
        axis.text.y  = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  scale_y_continuous(limits=c(0,50)) +
  xlab("Model") + 
  ylab("\u0394 AIC") +
  ggtitle("Chicks resighting model selection")
Bootstrap_Delta_AIC_plot_C

Bootstrap_AIC_weight_plot_C <- 
  ggplot(aes(y = weight, x = p), data = chick_AIC_tables) + 
  theme_bw() +
  geom_boxplot(width = 0.3, fill = "grey70", outlier.size = 0.5) +
  theme(text = element_text(family="Arial"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size=10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size=12, margin = margin(0, 15, 0, 0)),
        axis.text.y  = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  scale_y_continuous(limits=c(0,1)) +
  xlab("Model") + 
  ylab("AIC weight")
Bootstrap_AIC_weight_plot_C

# jpeg(filename = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/figs/AIC_multi_plot_chicks.jpeg",
#      quality = 100,
#      width = 6,
#      height = 9,
#      units = "in",
#      res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(Bootstrap_Delta_AIC_plot_C, vp = vplayout(1:2, 1))  # key is to define vplayout
print(Bootstrap_AIC_weight_plot_C, vp = vplayout(3:5, 1))
# dev.off()

Bootstrap_Delta_AIC_plot_F_A <- 
  ggplot(aes(y = DeltaAICc, x = p), data = fldg_ad_AIC_tables) + 
  theme_bw() +
  #geom_violin(fill = "grey40") +
  geom_boxplot(width = 0.3, fill = "grey70", outlier.size = 0.5) +
  theme(text = element_text(family="Arial"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=12, margin = margin(0, 18, 0, 0)),
        axis.text.y  = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        plot.margin = unit(c(0.5,0.5,0,0.8), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  scale_y_continuous(limits=c(0,50)) +
  xlab("Model") + 
  ylab("\u0394 AIC")  +
  ggtitle("Fledgling and adult resighting model selection")
Bootstrap_Delta_AIC_plot_F_A

Bootstrap_AIC_weight_plot_F_A <- 
  ggplot(aes(y = weight, x = p), data = fldg_ad_AIC_tables) + 
  theme_bw() +
  geom_boxplot(width = 0.3, fill = "grey70", outlier.size = 0.5) +
  theme(text = element_text(family="Arial"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size=10, angle = 45, hjust = 1), 
        axis.title.y = element_text(size=12, margin = margin(0, 18, 0, 0)),
        axis.text.y  = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        panel.margin = unit(0.75, "lines"),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(limits=c(0,1)) +
  xlab("Model") + 
  ylab("AIC weight")
Bootstrap_AIC_weight_plot_F_A

# jpeg(filename = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/figs/AIC_multi_plot_fledgling_adult.jpeg",
#      quality = 100,
#      width = 6,
#      height = 9,
#      units = "in",
#      res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(Bootstrap_Delta_AIC_plot_F_A, vp = vplayout(1:2, 1))  # key is to define vplayout
print(Bootstrap_AIC_weight_plot_F_A, vp = vplayout(3:5, 1))
# dev.off()

###############################################################################
################ LTRE plot ####################################################
survival_rates_boot$iter <- as.factor(survival_rates_boot$iter)
survival_rates_boot_summary <- 
  survival_rates_boot %>%
  dplyr::group_by(sex_age) %>%
  dplyr::summarise(Avg = mean(estimate),
                   lcl = stats::quantile(estimate, (1 - CI)/2, na.rm = TRUE),
                   ucl = stats::quantile(estimate, 1 - (1 - CI)/2, na.rm = TRUE))
survival_rates_boot_summary <- as.data.frame(survival_rates_boot_summary)
deterministic_list <- list(F_Chk_survl = survival_rates_boot_summary[2,2],
                           F_Fdg_survl = survival_rates_boot_summary[3,2],
                           F_Adt_survl = survival_rates_boot_summary[1,2],
                           M_Chk_survl = survival_rates_boot_summary[5,2],
                           M_Fdg_survl = survival_rates_boot_summary[6,2],
                           M_Adt_survl = survival_rates_boot_summary[4,2],
                           HSR = HSR,
                           RF = RF)
# deterministic_list_male <- list(F_Chk_survl = survival_rates_boot_summary[5,2],
#                            F_Fdg_survl = survival_rates_boot_summary[6,2],
#                            F_Adt_survl = survival_rates_boot_summary[4,2],
#                            M_Chk_survl = survival_rates_boot_summary[5,2],
#                            M_Fdg_survl = survival_rates_boot_summary[6,2],
#                            M_Adt_survl = survival_rates_boot_summary[4,2],
#                            HSR = 0.5,
#                            RF = RF)
# deterministic_list_female <- list(F_Chk_survl = survival_rates_boot_summary[2,2],
#                                 F_Fdg_survl = survival_rates_boot_summary[3,2],
#                                 F_Adt_survl = survival_rates_boot_summary[1,2],
#                                 M_Chk_survl = survival_rates_boot_summary[2,2],
#                                 M_Fdg_survl = survival_rates_boot_summary[3,2],
#                                 M_Adt_survl = survival_rates_boot_summary[1,2],
#                                 HSR = 0.5,
#                                 RF = RF)
matrix_structure <- expression(
  # top row of matrix
  0, (RF * (1 - HSR))/2, 0, (RF * (1 - HSR))/2,
  # second row of matrix
  (F_Chk_survl * F_Fdg_survl), F_Adt_survl, 0, 0,
  # third row of matrix
  0, (RF * HSR)/2, 0, (RF * HSR)/2,
  # fourth row of matrix
  0, 0, (M_Chk_survl * M_Fdg_survl), M_Adt_survl
)

deterministic_matrix <- plover_matrix(deterministic_list)

A <- matrix(c(0,0.5302839,0,0.5302839,0.07052423,0.6877241,0,0,0,0.5006591,0,0.5006591,0,0,0.1190305,0.6911536), byrow=TRUE, ncol=4)
B <- A
B[1,4] <- 0
B[3,4] <- 0
B
deterministic_matrix_male <- plover_matrix(deterministic_list_male)
deterministic_matrix_female <- plover_matrix(deterministic_list_female)

eigen(deterministic_matrix)$values[1]
eigen(deterministic_matrix_male)$values[1]
eigen(deterministic_matrix_female)$values[1]

deterministic_ASR <- 
  matrix_ASR(A = deterministic_matrix)
deterministic_ASR$ASR
deterministic_LTRE <- 
  ASR_perturbation(elements = matrix_structure,
                   VR_list = deterministic_list, 
                   freq_dep_ASR = deterministic_ASR)
cbPalette <- c("#A6A6A6", "#D9D9D9", "#D9D9D9", "#D9D9D9")

# plot the comparative LTRE results
Background_LTRE <-
  ggplot2::ggplot(data = deterministic_LTRE$LTRE,
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
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  scale_x_discrete(labels = c("Adult survival" = expression(Adult["\u03D5"]),
                              "Fledgling survival" = expression(Fledgling ["\u03D5"]),
                              "Chick survival" = expression(Chick ["\u03D5"]),
                              "Hatching sex ratio" = "Hatching SR"))
Background_LTRE

LTRE_bar_plot <-
  ggplot2::ggplot() +
  theme_bw() +
  coord_flip() +
  geom_bar(data = deterministic_LTRE$LTRE,
           aes(x = parameter, y = value, fill = parameter), color = "black", stat = "identity") +
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
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  scale_x_discrete(labels = c("Adult survival" = expression(Adult["\u03D5"]),
                              "Fledgling survival" = expression(Fledgling ["\u03D5"]),
                              "Chick survival" = expression(Chick ["\u03D5"]),
                              "Hatching sex ratio" = "Hatching SR"))
LTRE_bar_plot

# jpeg(filename = "/Users/Luke/Dropbox/Luke/R_projects/Ceuta_ASR_Matrix_Modeling/figs/final/LTRE_final.jpeg",
#      quality = 100,
#      width = 4,
#      height = 5,
#      units = "in",
#      res = 300)

grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) )
print( Background_LTRE + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
print( LTRE_bar_plot + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
# dev.off()

###############################################################################
############## Mating system plot #############################################
mating_df <- 
  breeding_data[which(!is.na(breeding_data$female) & !is.na(breeding_data$male)),]
mating_df$pair <- as.factor(paste(mating_df$female, mating_df$male, sep = "-"))
females <- reshape2::dcast(mating_df, female  ~ year)
males <- reshape2::dcast(mating_df, male  ~ year)
number_males_p_female <- 
  stats::aggregate(male ~ female, mating_df, function(x) length(unique(x)))
number_females_p_male <- 
  stats::aggregate(female ~ male, mating_df, function(x) length(unique(x)))
females <- dplyr::inner_join(females, number_males_p_female)
females[,c(2:8)] <- 
  lapply(females[,c(2:8)], as.numeric)
males <- dplyr::inner_join(males, number_females_p_male)
males[,c(2:8)] <- 
  lapply(males[,c(2:8)], as.numeric)
females$attempts <- rowSums(females[, c(2:8)])
males$attempts <- rowSums(males[, c(2:8)])
females$years <- rowSums(females[, c(2:8)] > 0)
males$years <- rowSums(males[, c(2:8)] > 0)
females_no_1 <- dplyr::filter(females, male  != 1 | years != 1 | attempts != 1)
males_no_1 <- dplyr::filter(males, female  != 1 | years != 1 | attempts != 1)
females_no_1$sex <- "Female"
females_no_1$sex <- as.factor(females_no_1$sex)
colnames(females_no_1)[c(1,9)] <- c("focal", "mate")
males_no_1$sex <- "Male"
males_no_1$sex <- as.factor(males_no_1$sex)
colnames(males_no_1)[c(1,9)] <- c("focal", "mate")
mating <- rbind(females_no_1, males_no_1)
mating$status <- ifelse(mating$mate == 1 & mating$years == mating$attempts, 
                        "Monogamous between years",
                        ifelse(mating$mate == 1 & mating$years < mating$attempts, 
                               "Monogamous within years",
                               ifelse(mating$mate > 1 & mating$years == mating$attempts, 
                                      "Polygamous between years",
                                      ifelse(mating$mate > 1 & mating$years < mating$attempts, 
                                             "Polygamous within years", "XXX"))))
mating$no_mates_per_year <- mating$mate/mating$years
chisq.test(table(mating$sex, mating$status)[,c(3,4)])
chisq.test(table(mating$sex, mating$status)[,c(1,2)])
chisq.test(table(mating$sex, mating$status))
mating$status <- factor(mating$status, 
                        levels = c("Monogamous within years",
                                   "Monogamous between years",
                                   "Polygamous between years",
                                   "Polygamous within years"))
sample_sizes_sex <- 
  stats::aggregate(focal ~ sex, data = mating, FUN = function(x){NROW(x)})
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
                filename = "Mating_system_final.jpg",
                path = "figs/final/",
                width = 3,
                height = 5, units = "in",
                dpi = 300)