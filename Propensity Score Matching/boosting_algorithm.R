## Propensity Score Matching - Boosting
## Kaelan Yu


## 0. Preliminaries

# remove all variables
rm(list = ls())

# set current working directory
setwd("C:/Users/Kaelan/Desktop/Research Project/Propensity Score Matching")

# do not interpret strings as factors
options(stringsAsFactors = FALSE)

# number of age groups
num_age_groups <- 3

## 1. Libraries

# to perform data manipulation using SQL like operations
library(dplyr)

# to be able to output tables to LaTeX
library(xtable)

# to be able to create binary dummy variables from categorical variables
library(fastDummies)

# to gain access to the random forest method used in the 
# propensity score matching reference paper
library(twang)




## 2. Import Data sets

# import liver and immuno data sets (liver_dat_cleaned and dat_small)
liver_dat <- read.csv(file = "liver_dat_large.csv", header = TRUE)
immuno_dat <- read.csv(file = "immuno_dat_small.csv", header = TRUE)

# get the TRR_ID_CODE and induction types from the immuno data set
immuno_dat <- dplyr::select(immuno_dat, TRR_ID_CODE, TYPE)

# perform an inner join on the two data sets on TRR_ID_CODE
dat <- inner_join(liver_dat, immuno_dat, by = "TRR_ID_CODE")




## 3. Prepare Data sets for Propensity Score Matching

# store AGE variable as separate variable for convenience
age <- dat$AGE

# between 18 and 49 logical
between_18_49 <- (age >= 18) & (age <= 49)

# between 50 and 64 logical
between_50_64 <- (age >= 50) & (age <= 64)

# 65 and above logical
above_65 <- (age >= 65)

# split data set into 3 smaller data sets based on age group
dat_18_49 <- dat[between_18_49, ]
dat_50_64 <- dat[between_50_64, ]
dat_above_65 <- dat[above_65, ]

# covariates for propensity score matching
variables <- c("TYPE", "AGE", "GENDER", "ALBUMIN_TX", "ETHCAT", "BMI_CALC", "DIAB", 
               "DIAL_TX", "MELD_PELD_LAB_SCORE", "ASCITES_TX", "CREAT_TX", "ENCEPH_TX", "DIAG", 
               "TBILI_TX", "INR_TX", "BACT_PERIT_TCR", "FUNC_STAT_TRR", "DON_TY", 
               "AGE_DON", "GENDER_DON", "ETHCAT_DON", "CARDARREST_NEURO", "HGT_CM_DON_CALC", 
               "HEP_C_ANTI_DON", "LITYP", "PSTATUS")

dat_vector <- list(dat_18_49, dat_50_64, dat_above_65)


# This function takes in a list of data sets (vec) and cleans the list of 
# data sets as follows:
#    1. assigns patients with a specific UNOS database code to a certain
#       group for categorical variables
#    2. creates binary dummy variables out of categorical variables using the 
#       fastDummies package and eliminates original categorical variable
get_cleaned_datasets <- function(vec) {
  cleaned_vector <- list()
  n <- length(vec)
  # for each data set
  for (i in 1:n) {
    # omit rows with NA values
    cleaned_vector[[i]] <- na.omit(vec[[i]][, variables])
    dataset <- cleaned_vector[[i]]
    
    ## Gender of the recipient - male or female
    dataset[, "GENDER"][dataset[, "GENDER"] == "M"] <- "Male"
    dataset[, "GENDER"][dataset[, "GENDER"] == "F"] <- "Female"
    
    ## Gender of the donor - male or female
    dataset[, "GENDER_DON"][dataset[, "GENDER_DON"] == "M"] <- "Male"
    dataset[, "GENDER_DON"][dataset[, "GENDER_DON"] == "F"] <- "Female"
    
    ## Diabetes - Yes/no
    dataset[, "DIAB"][dataset[, "DIAB"] == "N"] <- "No"
    dataset[, "DIAB"][dataset[, "DIAB"] == "Y"] <- "Yes"
    
    ## Dialysis at Time of Transplant - yes/no
    dataset[, "DIAL_TX"][dataset[, "DIAL_TX"] == "N"] <- "No"
    dataset[, "DIAL_TX"][dataset[, "DIAL_TX"] == "Y"] <- "Yes"
    
    ## Ethnicity - White, Black, Hispanic, Asian, Other
    dataset[, "ETHCAT"][dataset[, "ETHCAT"] == 1] <- "White"
    dataset[, "ETHCAT"][dataset[, "ETHCAT"] == 2] <- "Black"
    dataset[, "ETHCAT"][dataset[, "ETHCAT"] == 4] <- "Hispanic"
    dataset[, "ETHCAT"][dataset[, "ETHCAT"] == 5] <- "Asian"
    dataset[, "ETHCAT"][dataset[, "ETHCAT"] %in% c(6, 7, 9, 998)] <- "Other"
    
    ## Functional Status - No Assist, Some Assist, Total Assist, Other
    dataset[, "FUNC_STAT_TRR"][dataset[, "FUNC_STAT_TRR"] %in% 
                                 c(1, 2070, 2080, 2090, 2100, 4070, 4080, 4090, 4100)] <- "No Assist"
    dataset[, "FUNC_STAT_TRR"][dataset[, "FUNC_STAT_TRR"] %in% 
                                 c(2, 2040, 2050, 2060, 4040, 4050, 4060)] <- "Some Assist"
    dataset[, "FUNC_STAT_TRR"][dataset[, "FUNC_STAT_TRR"] %in% 
                                 c(3, 2010, 2020, 2030, 4010, 4020, 4030)] <- "Total Assist"
    dataset[, "FUNC_STAT_TRR"][dataset[, "FUNC_STAT_TRR"] %in% 
                                 c(996, 998)] <- "Other"
    
    # Diabetes - yes/no
    dataset[, "DIAB"] <- ifelse(dataset[, "DIAB"] %in% c(2, 3, 4, 5), "Yes", "No")
    
    # Ascites at time of transplant - Absent, Slight, MOderate, Other
    dataset[, "ASCITES_TX"][dataset[, "ASCITES_TX"] == 1] <- "Absent"
    dataset[, "ASCITES_TX"][dataset[, "ASCITES_TX"] == 2] <- "Slight"
    dataset[, "ASCITES_TX"][dataset[, "ASCITES_TX"] == 3] <- "Moderate"
    dataset[, "ASCITES_TX"][dataset[, "ASCITES_TX"] == 4] <- "Other"
    
    # Encephelopathy at time of transplant - None, Stage 1-2, Stage 3-4
    dataset[, "ENCEPH_TX"][dataset[, "ENCEPH_TX"] == 1] <- "None"
    dataset[, "ENCEPH_TX"][dataset[, "ENCEPH_TX"] == 2] <- "Stage 1-2"
    dataset[, "ENCEPH_TX"][dataset[, "ENCEPH_TX"] == 3] <- "Stage 3-4"
    
    # Donor Type - Deceased or Living
    dataset[, "DON_TY"] <- ifelse(dataset[, "DON_TY"] == "C", "Deceased", 
                                  "Living")
    
    # Donor Ethnicity - White, Black, Hispanic, Asian, Other
    dataset[, "ETHCAT_DON"][dataset[, "ETHCAT_DON"] == 1] <- "White"
    dataset[, "ETHCAT_DON"][dataset[, "ETHCAT_DON"] == 2] <- "Black"
    dataset[, "ETHCAT_DON"][dataset[, "ETHCAT_DON"] == 4] <- "Hispanic"
    dataset[, "ETHCAT_DON"][dataset[, "ETHCAT_DON"] == 5] <- "Asian"
    dataset[, "ETHCAT_DON"][dataset[, "ETHCAT_DON"] %in% c(6, 7, 9, 998)] <- "Other"
    
    # deceased donor cardiac arrest - yes/no
    dataset[, "CARDARREST_NEURO"] <- ifelse(dataset[, "CARDARREST_NEURO"] == "Y", 
                                            "Yes", "No")
    
    
    # Diagnosis - HCV, HBC, alcoholic disease, PBC, PSC, autoimmune, 
    #             cryptogenic/fatty liver, other
    dataset[, "DIAG"][dataset[, "DIAG"] %in% c(4104, 4204, 4216, 4593)] <- "HCV"
    dataset[, "DIAG"][dataset[, "DIAG"] %in% c(4102, 4106, 4107, 4202, 4592)] <- "HBV"
    dataset[, "DIAG"][dataset[, "DIAG"] %in% c(4215, 4217)] <- "Alcoholic"
    dataset[, "DIAG"][dataset[, "DIAG"] == 4220] <- "PBC"
    
    dataset[, "DIAG"][dataset[, "DIAG"] %in% c(4240, 4241, 4242, 4245)] <- "PSC"
    dataset[, "DIAG"][dataset[, "DIAG"] == 4212] <- "Autoimmune"
    
    dataset[, "DIAG"][dataset[, "DIAG"] %in% c(4208, 4213, 4214)] <- "Cryptogenic/Fatty Liver"
    dataset[, "DIAG"][!(dataset[, "DIAG"] %in% c("HCV", "HBV", "Alcoholic", "PBC", "PSC", 
                                                 "Autoimmune", "Cryptogenic/Fatty Liver"))] <- "Other"
    
    
    # Spontaneous Bacterial Peritonitis - yes/no
    dataset[, "BACT_PERIT_TCR"] <- ifelse(dataset[, "BACT_PERIT_TCR"] == "Y", "Yes", "No")
    
    # result of deceased donor antibody to hepatitis c - yes/no
    dataset[, "HEP_C_ANTI_DON"] <- ifelse(dataset[, "HEP_C_ANTI_DON"] == "P", "Yes", "No")
    
    # split liver - yes/no
    dataset[, "LITYP"] <- ifelse(dataset[, "LITYP"] %in%
                                   c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17), "Yes", "No")
    
    dataset <- dataset[, order(names(dataset))]
    
    cleaned_vector[[i]] <- dataset
    
    # apply fastDummies package to create dummy variables
    cleaned_vector[[i]] <- dummy_cols(cleaned_vector[[i]])
    
    # predictor variable - patient status (dead or alive)
    cleaned_vector[[i]][, "PATIENTSTATUS"] <- ifelse(cleaned_vector[[i]][, "PSTATUS"] == 1, "Dead", "Alive")
    
    
    # remove original categorical variables
    cleaned_vector[[i]] <- subset(cleaned_vector[[i]], select = -c(GENDER, ETHCAT, DIAB, DIAL_TX, ASCITES_TX, 
                                                                   ENCEPH_TX, DIAG, BACT_PERIT_TCR, 
                                                                   FUNC_STAT_TRR, DON_TY, GENDER_DON, 
                                                                   ETHCAT_DON, CARDARREST_NEURO, 
                                                                   HEP_C_ANTI_DON, LITYP))
    cleaned_vector[[i]] <- dplyr::select(cleaned_vector[[i]], -c(`TYPE_NO INDUCTION`, TYPE_TTG, TYPE_ONLY, TYPE_SIMULECT))
    
  }
  
  return(cleaned_vector)
}

induction_cleaned_vector <- get_cleaned_datasets(dat_vector)

cleaned_dat_18_49 <- induction_cleaned_vector[[1]]
cleaned_dat_50_64 <- induction_cleaned_vector[[2]]
cleaned_dat_above_65 <- induction_cleaned_vector[[3]]


# remove TYPE and PATIENTSTATUS to create predictor variables
predictor_variables <- names(cleaned_dat_18_49)[-c(11, 60)]

mnps.induction <- 
mnps(as.factor(TYPE) ~ 
       AGE + AGE_DON + ALBUMIN_TX + BMI_CALC + CREAT_TX + HGT_CM_DON_CALC + 
       INR_TX + MELD_PELD_LAB_SCORE + PSTATUS + TBILI_TX + 
       ASCITES_TX_Slight + ASCITES_TX_Moderate + ASCITES_TX_Absent + ASCITES_TX_Other + 
       BACT_PERIT_TCR_No + BACT_PERIT_TCR_Yes + 
       CARDARREST_NEURO_No + CARDARREST_NEURO_Yes + 
       DIAB_No + DIAB_Yes + 
       DIAL_TX_No + DIAL_TX_Yes + 
       DON_TY_Deceased + DON_TY_Living + 
       ETHCAT_White + ETHCAT_Asian + ETHCAT_Black + ETHCAT_Hispanic + ETHCAT_Other + 
       ETHCAT_DON_White + ETHCAT_DON_Asian + ETHCAT_DON_Black + ETHCAT_DON_Hispanic + ETHCAT_DON_Other + 
       GENDER_Female + GENDER_Male + 
       GENDER_DON_Female + GENDER_DON_Male + 
       HEP_C_ANTI_DON_No + HEP_C_ANTI_DON_Yes + 
       LITYP_No + LITYP_Yes,
     
data = cleaned_dat_18_49, estimand = "ATE", 
stop.method = c("es.mean", "ks.mean"), n.trees = 500)
# n.trees = number of iterations

