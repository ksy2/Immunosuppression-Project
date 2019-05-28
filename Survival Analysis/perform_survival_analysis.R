## Survival Analysis Code 
## Kaelan Yu

# remove all variables
rm(list = ls())


# set current working directory to where the liver_dat.sas7bdat file is located
setwd("C:/Users/Kaelan/Desktop/Research Project/Survival Analysis")

# suppress command output by redirecting it to discard.txt
sink(file = "discard.txt")

# to perform survival analysis
library(survival)
library(survminer)

# to perform data manipulation using SQL like operations
library(dplyr)

# to be able to output tables to LaTeX
library(xtable)

# do not treat strings as factors
options(stringsAsFactors = FALSE)

# We will work directly with the cleaned versions of the liver and 
# immuno data sets.

# import cleaned liver data set
liver_dat <- read.csv("liver_dat_cleaned.csv", header = T)

# import cleaned immuno data set
immuno_dat <- read.csv("immuno_dat_cleaned.csv", header = T)

# use the select and inner_join functions from the dplyr library

# select relevant attributes from the immuno_dat data set
immuno_dat <- select(immuno_dat, TRR_ID_CODE, TTG_STEROIDS, SIMULECT_STEROIDS, STEROIDS_ONLY)

# perform an inner join on TRR_ID_CODE (encrypted transplant identifier) 
survival_dat <- inner_join(liver_dat, immuno_dat, by = "TRR_ID_CODE")

# select relevant variables for survival analysis:
# (1) Patient Survival: time = PTIME, event = PSTATUS
# (2) Graft Survival:   time = GTIME, event = GSTATUS
survival_dat <- select(survival_dat, AGE, TRR_ID_CODE, PSTATUS, PTIME,GSTATUS, GTIME, 
                    TTG_STEROIDS, SIMULECT_STEROIDS, STEROIDS_ONLY)

# we want to convert the time variable from days to years
survival_dat[, "PYEARS"] <- survival_dat$PTIME / 365
survival_dat[, "GYEARS"] <- survival_dat$GTIME / 365

# make a new "TYPE" variable for the induction type:
survival_dat[, "TYPE"] <- NA
survival_dat$TYPE[survival_dat$TTG_STEROIDS == 1] <- "TTG"
survival_dat$TYPE[survival_dat$SIMULECT_STEROIDS == 1] <- "SIMULECT"
survival_dat$TYPE[survival_dat$STEROIDS_ONLY == 1] <- "ONLY"
survival_dat$TYPE[survival_dat$TTG_STEROIDS == 0 & survival_dat$SIMULECT_STEROIDS == 0 
               & survival_dat$STEROIDS_ONLY == 0] <- "NO INDUCTION"


# create a variable for the AGE attribute from the data set for convenience
age <- survival_dat$AGE

# split survival_dat data set into 3 data sets for each age group
dat_18_49 <- survival_dat[(age >= 18 & age <= 49), ]
dat_50_64 <- survival_dat[(age >= 50 & age <= 64), ]
dat_above_65 <- survival_dat[(age >= 65), ]

# store the data sets in a list data structure for easy access
dat_vector <- list(dat_18_49, dat_50_64, dat_above_65)

# store the number of age groups as a variable
num_age_groups <- 3


# Takes in a data set, name of age group, and limits for y-axis with default
# of c(0, 1) and plots the patient survival curves for the 3 induction groups.
plot_survival_curves <- function(dataset, age, yrange = c(0, 1)) {
  surv <- Surv(time = dataset$PYEARS, event = dataset$PSTATUS)
  fit_TTG <- survfit(surv ~ dataset$TTG_STEROIDS == 1, type = "kaplan-meier",
                                    subset = dataset$TYPE == "TTG")

  fit_sim <- survfit(surv ~ dataset$SIMULECT_STEROIDS == 1, type = "kaplan-meier",
                                    subset = dataset$TYPE == "SIMULECT")

  fit_only <- survfit(surv ~ dataset$STEROIDS_ONLY == 1, type = "kaplan-meier",
                                     subset = dataset$TYPE == "ONLY")

  fit <- list(Thymoglobulin = fit_TTG, Simulect = fit_sim, Steroids_Only = fit_only)

  ggsurvplot(fit, data = dataset, combine = TRUE, risk.table = "percentage", title = paste("Age", age, "Survival Curves"),
             xlab = "Number of Years", ylab = "Patient Survival Probability", ylim = yrange,
             censor = F, size = 1.5, break.x.by = 1, legend.labs = c("TTG",
                                                                     "Simulect",
                                                                     "Steroids Only"))
}

# Takes in a data set, name of age group, and limits for y-axis with default
# of c(0, 1) and plots the graft survival curves for the 3 induction groups.
plot_graft_curves <- function(dataset, age, yrange = c(0, 1)) {
  surv <- Surv(time = dataset$GYEARS, event = dataset$GSTATUS)
  fit_TTG <- survfit(surv ~ dataset$TTG_STEROIDS == 1, type = "kaplan-meier", 
                              subset = dataset$TYPE == "TTG")
  
  fit_sim <- survfit(surv ~ dataset$SIMULECT_STEROIDS == 1, type = "kaplan-meier", 
                              subset = dataset$TYPE == "SIMULECT")
  
  fit_only <- survfit(surv ~ dataset$STEROIDS_ONLY == 1, type = "kaplan-meier", 
                               subset = dataset$TYPE == "ONLY")
  
  fit <- list(Thymoglobulin = fit_TTG, Simulect = fit_sim, Steroids_Only = fit_only)
  
  ggsurvplot(fit, data = dataset, combine = TRUE, risk.table = "percentage", title = paste("Age", age, "Survival Curves"),
             xlab = "Number of Years", ylab = "Graft Survival Probability", ylim = yrange, 
             censor = F, size = 1.5, break.x.by = 1, legend.labs = c("Steroids + Thymoglobulin", 
                                                                     "Steroids + Simulect", 
                                                                     "Steroids Only"))
}

# Plot Patient Survival Curves

png(filename = "survival_18_49.png", width = 800, height = 480)
plot_survival_curves(dat_18_49, "18 - 49", yrange = c(0.3, 1))
graphics.off()

png(filename = "survival_50_64.png", width = 800, height = 480)
plot_survival_curves(dat_50_64, "50 - 64", yrange = c(0.3, 1))
graphics.off()

png(filename = "survival_above_65.png", width = 800, height = 480)
plot_survival_curves(dat_above_65, "65+", yrange = c(0.3, 1))
graphics.off()


## Plot Graft Survival Curves

png(filename = "graft_18_49.png", width = 800, height = 480)
plot_graft_curves(dat_18_49, "18 - 49", yrange = c(0.3, 1))
graphics.off()

png(filename = "graft_50_64.png", width = 800, height = 480)
plot_graft_curves(dat_50_64, "50 - 64", yrange = c(0.3, 1))
graphics.off()

png(filename = "graft_above_65.png", width = 800, height = 480)
plot_graft_curves(dat_above_65, "65+", yrange = c(0.3, 1))
graphics.off()



# store the name of the age groups in a vector
age_groups <- c("18 - 49", "50 - 64", "65+")


# p-Values from the log-rank test for survival analysis

# create a vector to store the p_values
survival_pvalue_TTG_vec <- rep(0, num_age_groups) # TTG vs steroids
survival_pvalue_sim_vec <- rep(0, num_age_groups) # simulect vs steroids
graft_pvalue_TTG_vec <- rep(0, num_age_groups)    # TTG vs steroids
graft_pvalue_sim_vec <- rep(0, num_age_groups)    # simulect vs steroids

# get the p-values from survival analysis
for (i in 1:num_age_groups) {
  dat <- dat_vector[[i]]
  
  survival_surv <- Surv(time = dat$PYEARS, event = dat$PSTATUS)
  survival_fit_TTG <- survfit(survival_surv ~ dat$TYPE, data = dat, type = "kaplan-meier", subset = dat$TYPE %in% c("TTG", "ONLY"))
  survival_fit_sim <- survfit(survival_surv ~ dat$TYPE, data = dat, type = "kaplan-meier", subset = dat$TYPE %in% c("SIMULECT", "ONLY"))
  survival_pvalue_TTG <- surv_pvalue(survival_fit_TTG)$pval.txt
  survival_pvalue_TTG_vec[i] <- survival_pvalue_TTG
  survival_pvalue_sim <- surv_pvalue(survival_fit_sim)$pval.txt
  survival_pvalue_sim_vec[i] <- survival_pvalue_sim
  print(survival_pvalue_TTG)
  print(survival_pvalue_sim)
  
  graft_surv <- Surv(time = dat$GYEARS, event = dat$GSTATUS)
  graft_fit_TTG <- survfit(graft_surv ~ dat$TYPE, data = dat, type = "kaplan-meier", subset = dat$TYPE %in% c("TTG", "ONLY"))
  graft_fit_sim <- survfit(graft_surv ~ dat$TYPE, data = dat, type = "kaplan-meier", subset = dat$TYPE %in% c("SIMULECT", "ONLY"))
  graft_pvalue_TTG <- surv_pvalue(graft_fit_TTG)$pval.txt
  graft_pvalue_TTG_vec[i] <- graft_pvalue_TTG
  graft_pvalue_sim <- surv_pvalue(graft_fit_sim)$pval.txt
  graft_pvalue_sim_vec[i] <- graft_pvalue_sim
  print(graft_pvalue_TTG)
  print(graft_pvalue_sim)
}

## Input p-value results in a data frame
pvalue_table <- data.frame()

for (i in 1:num_age_groups) {
  pvalue_table["Patient Survival", i] <- ""
  pvalue_table["Patient Survival (TTG v. Steroids Only)", i] <- survival_pvalue_TTG_vec[i]
  pvalue_table["Patient Survival (Simulect v. Steroids Only)", i] <- survival_pvalue_sim_vec[i]
  pvalue_table["Graft Survival", i] <- ""
  pvalue_table["Graft Survival (TTG v. Steroids Only)", i] <- graft_pvalue_TTG_vec[i]
  pvalue_table["Graft Survival (Simulect v. Steroids Only)", i] <- graft_pvalue_sim_vec[i]
}

View(pvalue_table)

names(pvalue_table) <- age_groups

# store LaTeX output in pvalue.txt
sink(file = "pvalue.txt")

print("The LaTeX output of the p-value table is provided below.")
print(xtable(pvalue_table, auto = TRUE))
