## Propensity Score Matching
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

# suppress command output
sink(file = "discard.txt")


## 1. Libraries

# to perform data manipulation using SQL like operations
library(dplyr)

# to be able to output tables to LaTeX
library(xtable)

# to be able to create binary dummy variables from categorical variables
library(fastDummies)

# to gain access to the random forest method used in the 
# propensity score matching reference paper
library(tree)
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

# further split each smaller data set into 3 even smaller data sets
# based on induction drug group (TTG, simulect, steroids only)
dat_TTG_18_49 <- dat_18_49[dat_18_49$TYPE == "TTG", ]
dat_sim_18_49 <- dat_18_49[dat_18_49$TYPE == "SIMULECT", ]
dat_only_18_49 <- dat_18_49[dat_18_49$TYPE == "ONLY", ]

dat_TTG_50_64 <- dat_50_64[dat_50_64$TYPE == "TTG", ]
dat_sim_50_64 <- dat_50_64[dat_50_64$TYPE == "SIMULECT", ]
dat_only_50_64 <- dat_50_64[dat_50_64$TYPE == "ONLY", ]

dat_TTG_above_65 <- dat_above_65[dat_above_65$TYPE == "TTG", ]
dat_sim_above_65 <- dat_above_65[dat_above_65$TYPE == "SIMULECT", ]
dat_only_above_65 <- dat_above_65[dat_above_65$TYPE == "ONLY", ]


# covariates for propensity score matching
variables <- c("AGE", "GENDER", "ALBUMIN_TX", "ETHCAT", "BMI_CALC", "DIAB", 
               "DIAL_TX", "MELD_PELD_LAB_SCORE", "ASCITES_TX", "CREAT_TX", "ENCEPH_TX", "DIAG", 
               "TBILI_TX", "INR_TX", "BACT_PERIT_TCR", "FUNC_STAT_TRR", "DON_TY", 
               "AGE_DON", "GENDER_DON", "ETHCAT_DON", "CARDARREST_NEURO", "HGT_CM_DON_CALC", 
               "HEP_C_ANTI_DON", "LITYP", "PSTATUS")


# store all 9 data sets in a list for convenience
dat_vector <- list(dat_TTG_18_49, dat_sim_18_49, dat_only_18_49, 
                   dat_TTG_50_64, dat_sim_50_64, dat_only_50_64, 
                   dat_TTG_above_65, dat_sim_above_65, dat_only_above_65)


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
  }
  
  return(cleaned_vector)
}



# store the 9 cleaned data sets in a list data structure in 
# preparation for logistic regression:

# 1 - 3: age 18 - 49 (TTG, sim, only)
# 4 - 6: age 50 - 64 (TTG, sim, only)
# 7 - 9: age 65+     (TTG, sim, only)
induction_cleaned_vector <- get_cleaned_datasets(dat_vector)


# remove statistical anomalies from the 6th data set (age 50 - 64, steroids only)
induction_cleaned_vector[[6]] <- induction_cleaned_vector[[6]][-c(2091, 2161, 2236), ]


# information for the following piece of code:

# 1. age vector is a list of 3 lists, i.e. it's 2 dimensional
# 2. format: [[age group]][[induction group (TTG, sim, only)]]
# -  the first list contains the data sets for age 18 - 49 (TTG, sim, only)
# -  the second list contains the data sets for age 50 - 64 (TTG, sim, only)
# -  the third list contains the data sets for age 65+ (TTG, sim, only)

age_vector <- list()

for (i in 1:3) {
  age_vector[[i]] <- list(induction_cleaned_vector[[3 * i - 2]], 
                          induction_cleaned_vector[[3 * i - 1]], 
                          induction_cleaned_vector[[3 * i]])
}


## This is just a test for all the dimensions of the data set - columns should match

# for (i in 1:3) {
#   for (j in 1:3) {
#    print(dim(age_vector[[i]][[j]]))
#  }
# }



## 4. Run Logistic Regression

# propensity score will be a 3D list containing:

# - 3 age groups [[i]]: [[1]] = age 18 - 49, [[2]] = age 50 - 64, [[3]] = age 65+ which each contain
# - 3 total score vectors [[j]]: [[1]] = TTG results, [[2]] = sim results, [[3]] = only results, which each contain
# - 3 score vectors [[k]]: [[1]] = applying model to TTG, [[2]] = applying model to sim, [[3]] = applying model to only


propensity_score <- list()

# i = age group (18 - 49, 50 - 64, 65+)

# j = induction type (TTG, simulect, only)

# k = propensity score from model (TTG on TTG, simulect, only)

for (i in 1:3) {
  
  # make a new total_score_vector for each age group/induction group (which will contain 3 score_vectors)
  total_score_vector <- list()
  
  # for each induction type (1 = TTG, 2 = simulect, 3 = only)
  for (j in 1:3) {
    
    # run logistic regression (add option, family = "binomial") to predict the patient
    # status using the data set from age_vector [[age_group]] [[induction type]]
    glm.fits <- glm(PSTATUS ~ .- PSTATUS - PATIENTSTATUS, data = age_vector[[i]][[j]], 
                    family = "binomial")
    
    # obtain the propensity score from applying the fitted model
    # (ex: i = 1 applies the fitted model to the data sets: 
    # k = 1:3 (TTG, sim, only))
    
    # score_vector gets those 3 propensity score lists
    score_vector <- list()
    
    
    for (k in 1:3) {
      glm.probs <- as.data.frame(predict(glm.fits, age_vector[[i]][[k]], type = "response"))
      score_vector[[k]] <- glm.probs
    }
    
    # fill in the jth entries of the total score vectors
    total_score_vector[[j]] <- score_vector
  }
  
  propensity_score[[i]] <- total_score_vector
}

# get propensity score results for each age group: these are 2D lists containing [[j]] (induction type) and
# [[k]] (the results on each induction group)
age_18_49_propensity_scores <- propensity_score[[1]]
age_50_64_propensity_scores <- propensity_score[[2]]
age_above_65_propensity_scores <- propensity_score[[3]]

# get induction group propensity score results for each age group/induction group: these are 1D lists
# containing [[k]] (a list containing the application of the model to each induction group)
age_18_49_TTG_propensity_scores <- age_18_49_propensity_scores[[1]]
age_18_49_sim_propensity_scores <- age_18_49_propensity_scores[[2]]
age_18_49_only_propensity_scores <- age_18_49_propensity_scores[[3]]

age_50_64_TTG_propensity_scores <- age_50_64_propensity_scores[[1]]
age_50_64_sim_propensity_scores <- age_50_64_propensity_scores[[2]]
age_50_64_only_propensity_scores <- age_50_64_propensity_scores[[3]]

age_above_65_TTG_propensity_scores <- age_above_65_propensity_scores[[1]]
age_above_65_sim_propensity_scores <- age_above_65_propensity_scores[[2]]
age_above_65_only_propensity_scores <- age_above_65_propensity_scores[[3]]

# create a vector of 3 colors for stacked histograms
colors <- c(rgb(1, 0, 0, 0.7), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.3))

# age groups
age_groups <- c("18 - 49", "50 - 64", "65+")

# induction groups
induction_groups <- c("TTG", "Simulect", "Steroids Only")


## 5. Data Visualization


## 5-1. Plot Propensity Scores as Stacked Histograms

# this function plots a stacked histogram of the propensity scores by taking 
# in the following parameters:

# 1. propensity_score_vector = propensity score for each age group/induction group
# 2. age = name of age group
# 3. induction = name of induction group
# 4. ylimit = y-axis limits
# 5. xlimit = x-axis limits
# 6. nbreaks = number of brwaks in the histogram
plot_histograms <- function(propensity_score_vector, age, induction, ylimit, xlimit, nbreaks) {
  hist(as.numeric(unlist(propensity_score_vector[[1]])), col = colors[1], breaks = nbreaks, ylim = ylimit, 
       main = paste("Age ", age, "Propensity Scores from \n Logistic Regression using", 
                    induction, "\n(", nbreaks, "Breaks )"), xlab = "Propensity Score", xlim = xlimit)
  hist(as.numeric(unlist(propensity_score_vector[[2]])), col = colors[2], breaks = nbreaks, add = TRUE)
  hist(as.numeric(unlist(propensity_score_vector[[3]])), col = colors[3], breaks = nbreaks, add = TRUE)
  legend(x = "right", legend = c("TTG", "Simulect", "Steroids Only"), fill = colors)
}

# put all plots in a single pdf file - stacked_histograms.pdf
pdf(file = "stacked_histograms.pdf")

# age 18 - 49
plot_histograms(age_18_49_TTG_propensity_scores, "18 - 49", "TTG", c(0, 800), c(0, 1), 50)
plot_histograms(age_18_49_sim_propensity_scores, "18 - 49", "Simulect", c(0, 800), c(0, 1), 50)
plot_histograms(age_18_49_only_propensity_scores, "18 - 49", "Steroids Only", c(0, 800), c(0, 1), 50)

# age 50 - 64 group
plot_histograms(age_50_64_TTG_propensity_scores, "50 - 64", "TTG", c(0, 2500), c(0, 1), 50)
plot_histograms(age_50_64_sim_propensity_scores, "50 - 64", "Simulect", c(0, 2500), c(0, 1), 50)
plot_histograms(age_50_64_only_propensity_scores, "50 - 64", "Steroids Only", c(0, 2500), c(0, 1), 50)

# age 65+ group
plot_histograms(age_above_65_TTG_propensity_scores, "65+ ", "TTG", c(0, 450), c(0, 1), 50)
plot_histograms(age_above_65_sim_propensity_scores, "65+ ", "Simulect", c(0, 450), c(0, 1), 50)
plot_histograms(age_above_65_only_propensity_scores, "65+ ", "Steroids Only", c(0, 450), c(0, 1), 50)
graphics.off()


## 5-2. Kernel Density Plots

# bandwidth parameter
h = 0.01

pdf(file = "kernel_density_plots.pdf")

# age 18 - 49
plot(density(as.numeric(unlist(age_18_49_TTG_propensity_scores)), bw = h), col = "red", 
            xlab = "Propensity Score", main = paste("Kernel Density Plot of Propensity Scores from Logistic Regression
            Age 18 - 49, bandwidth =", h), ylim = c(0, 5), xlim = c(0, 1), lwd = 3)
lines(density(as.numeric(unlist(age_18_49_sim_propensity_scores)), bw = h), col = "blue", lwd = 3)
lines(density(as.numeric(unlist(age_18_49_only_propensity_scores)), bw = h), col = "green", lwd = 3)
legend(x = "right", legend = c("TTG", "Simulect", "Only"), fill = c("red", "blue", "green"))

# age 50 - 64
plot(density(as.numeric(unlist(age_50_64_TTG_propensity_scores), bw = h)), col = "red", 
     xlab = "Propensity Score", main = paste("Kernel Density Plot of Propensity Scores from Logistic Regression
            Age 50 - 64, bandwidth =", h), ylim = c(0, 7), xlim = c(0, 1), lwd = 3)
lines(density(as.numeric(unlist(age_50_64_sim_propensity_scores)), bw = h), col = "blue", lwd = 3)
lines(density(as.numeric(unlist(age_50_64_only_propensity_scores)), bw = h), col = "green", lwd = 3)
legend(x = "right", legend = c("TTG", "Simulect", "Only"), fill = c("red", "blue", "green"))

# age 65+
plot(density(as.numeric(unlist(age_above_65_TTG_propensity_scores)), bw = h), col = "red", 
     xlab = "Propensity Score", main = paste("Kernel Density Plot of Propensity Scores from Logistic Regression
                                             Age 65+, bandwidth =", h), ylim = c(0, 5), xlim = c(0, 1), lwd = 3)
lines(density(as.numeric(unlist(age_above_65_sim_propensity_scores)), bw = h), col = "blue", lwd = 3)
lines(density(as.numeric(unlist(age_above_65_only_propensity_scores)), bw = h), col = "green", lwd = 3)
legend(x = "right", legend = c("TTG", "Simulect", "Only"), fill = c("red", "blue", "green"))
graphics.off()



## 6. LaTeX Output

# Print out list of covariates for propensity score matching
sink("covariates.txt")
print("The LaTeX output of the list of covariates is provided below.")
x <- xtable(as.data.frame(names(induction_cleaned_vector[[1]])), auto = TRUE)
print(x, tabular.environment = "longtable", floating = FALSE)


# remove all variables
rm(list = ls())