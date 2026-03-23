###############################################
#OC Regression & Survival Analysis
###############################################

#Setup and data import

getwd()  #Check current working directory

#Read dataset (choose CSV interactively)
BS4 <- read.csv("OC_Dataset.csv")

#Drop unwanted ID columns and remove incomplete cases
BS4 <- BS4[, !names(BS4) %in% c("X.1", "X.2")]
BS4 <- na.omit(BS4)

#Load required packages
library(lattice)
library(tidyr)
library(epiDisplay)

attach(BS4)     #For direct column access
names(BS4)      #Inspect variable names


###############################################
# Q1: Linear regression with age as outcome
###############################################

par(mfrow = c(1, 1))

#Scatter plot of age vs meal calories
plot(age, meal.cal,
     xlab = "Age (years)",
     ylab = "Meal calories",
     main = "Age vs Meal Calories")

#Simple linear regression: age ~ meal.cal
age_cal <- lm(age ~ meal.cal)
summary(age_cal)

#Overlay regression line
plot(age ~ meal.cal,
     xlab = "Meal calories",
     ylab = "Age (years)")
abline(coef = coef(age_cal), col = "red")

#Diagnostic plots for linear model assumptions
par(mfrow = c(2, 2))
plot(age_cal)

#R² is low (~5%), so age is poorly explained by meal.cal alone.

#Univariate linear models with age as outcome vs each predictor
model_test1 <- lm(age ~ Survival_months)
summary(model_test1)

model_test2 <- lm(age ~ status)
summary(model_test2)

model_test3 <- lm(age ~ Treatment_type)
summary(model_test3)

model_test4 <- lm(age ~ as.factor(ph.ecog))
summary(model_test4)

model_test5 <- lm(age ~ meal.cal)
summary(model_test5)

model_test6 <- lm(age ~ wt.loss)
summary(model_test6)

#Multivariable models with age as outcome
age_multi  <- lm(age ~ status + as.factor(ph.ecog) + meal.cal)
summary(age_multi)

age_multi2 <- lm(age ~ as.factor(ph.ecog) + meal.cal)
summary(age_multi2)

age_multi3 <- lm(age ~ Treatment_type + as.factor(ph.ecog) + meal.cal)
summary(age_multi3)

#Same model as GLM with Gaussian family (equivalent to lm)
age_multi5 <- glm(age ~ as.factor(ph.ecog) + meal.cal, family = gaussian)
summary(age_multi5)

#Formatted regression table
regress.display(age_multi5)

#Diagnostics for the GLM
par(mfrow = c(2, 2))
plot(age_multi5)


###############################################
# Q2: Logistic regression for treatment type
# Outcome: Treatment_type (Chemotherapy vs Radiation)
###############################################

#Recode Treatment_type to numeric:
#Radiation = 1, Chemotherapy = 0
Treatment_type[which(Treatment_type == "Radiation")]    <- 1
Treatment_type[which(Treatment_type == "Chemotherapy")] <- 0
Treatment_type <- as.numeric(Treatment_type)

Treatment_type  #Check recode

#Univariate logistic regressions
model_test7  <- glm(Treatment_type ~ age, family = binomial)
summary(model_test7)

model_test8  <- glm(Treatment_type ~ Survival_months, family = binomial)
summary(model_test8)

model_test9  <- glm(Treatment_type ~ status, family = binomial)
summary(model_test9)

model_test10 <- glm(Treatment_type ~ as.factor(ph.ecog), family = binomial)
summary(model_test10)

model_test11 <- glm(Treatment_type ~ meal.cal, family = binomial)
summary(model_test11)

model_test12 <- glm(Treatment_type ~ wt.loss, family = binomial)
summary(model_test12)

#Multivariable logistic model for treatment allocation
treat_multi4 <- glm(Treatment_type ~ status + meal.cal + wt.loss,
                    family = binomial)
summary(treat_multi4)

#Odds ratios and CI
logistic.display(treat_multi4)

#Example: manually computed odds for a specific covariate profile
#(status = 1, meal.cal = 875, wt.loss = 4), using rounded coefficients
OR <- exp(-1.88 + 0.99 + 0.001 * 875 + 0.04 * 4)
OR

#ROC curve and AUC for the logistic model
curve <- lroc(treat_multi4)
curve$auc

#Likelihood ratio test vs null model
pchisq(treat_multi4$null.deviance - treat_multi4$deviance,
       treat_multi4$df.null - treat_multi4$df.residual,
       lower.tail = FALSE)

#Check overall fit via residual deviance per df
deviance(treat_multi4) / df.residual(treat_multi4)

#Influence diagnostics using Cook's distance
cooksD <- cooks.distance(treat_multi4)
plot(cooksD, ylab = "Cook's distance",
     main = "Influence diagnostics: logistic model")
n <- nrow(BS4)
abline(h = 4 / n, lty = 2, col = "red")  #Common rule-of-thumb threshold


###############################################
# Q3: Survival analysis by treatment type
# Question: Does survival differ by treatment type,
# and which covariates impact survival?
###############################################

library(survival)

#Convert Treatment_type back to factor labels for survival analysis
Treatment_type[which(Treatment_type == "1")] <- "Radiation"
Treatment_type[which(Treatment_type == "0")] <- "Chemotherapy"
Treatment_type <- as.character(Treatment_type)

#Kaplan–Meier curves by treatment type
model_fit_surv <- survfit(Surv(Survival_months, status) ~ Treatment_type)
model_fit_surv

plot(model_fit_surv,
     col  = 1:2,
     lwd  = 2,
     xlab = "Months",
     ylab = "Survival probability",
     main = "Kaplan–Meier Survival by Treatment Type")

legend("topright",
       c("Chemotherapy", "Radiation"),
       col = 1:2,
       lwd = 2,
       bty = "n")

#Log-rank test for equality of survivor functions
survdiff(Surv(Survival_months, status) ~ Treatment_type)

#Univariate Cox models for each covariate
model_cox  <- coxph(Surv(Survival_months, status) ~ Treatment_type)
summary(model_cox)

model_cox2 <- coxph(Surv(Survival_months, status) ~ age)
summary(model_cox2)

model_cox3 <- coxph(Surv(Survival_months, status) ~ as.factor(ph.ecog))
summary(model_cox3)

model_cox4 <- coxph(Surv(Survival_months, status) ~ meal.cal)
summary(model_cox4)

model_cox5 <- coxph(Surv(Survival_months, status) ~ wt.loss)
summary(model_cox5)

#Multivariable Cox models with different covariate sets
cox_multi  <- coxph(Surv(Survival_months, status) ~ Treatment_type + age + as.factor(ph.ecog))
summary(cox_multi)

cox_multi2 <- coxph(Surv(Survival_months, status) ~ Treatment_type + as.factor(ph.ecog))
summary(cox_multi2)

#Nicely formatted Cox model output
cox.display(cox_multi2, simplified = TRUE)

#Proportional hazards (PH) assumption diagnostics
library(survminer)

model_cox.zph <- cox.zph(cox_multi2)
print(model_cox.zph)

#Graphical check of PH assumption
ggcoxzph(model_cox.zph)
