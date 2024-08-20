library("survminer")
library(cutpointr)
library("ggplot2")
require("survival")
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(edgeR)
library(pROC)
library(verification)

cellFrac = read.csv("Tumor_Cell_Prop_4.csv")
cellFrac = cellFrac[which(cellFrac$TMA %in% c("YTMA_471")),]

train_data = cellFrac[,c(38:51)] # for tumor cells 
nCT = dim(train_data)[2]

#y_df = Surv(as.numeric(cellFrac$PFS_2Years_months, cellFrac$PFS_2Years_Index)) ###
y_df = Surv(as.numeric(cellFrac$PFS_5Years_months, cellFrac$PFS_5Years_Index)) ###

cell_counts = rep(0,nCT)
alpha = 1

for (cSeed in 1:100){
  
  print(cSeed)
  set.seed(cSeed)
  
  lambdas_to_try <- 10^seq(-3, 2, length.out = 100)
  ##
  lasso_cv <- cv.glmnet(as.matrix(train_data), y_df, alpha = alpha, lambda = lambdas_to_try,
                        standardize = FALSE, nfolds = 10, family = "cox", lower.limits = 0) ###
  #plot(lasso_cv)
  lambda_cv <- lasso_cv$lambda.min
  model_cv <- glmnet(as.matrix(train_data), y_df, alpha = alpha, lambda = lambda_cv, standardize = FALSE, family = "cox", lower.limits = 0) 
  coeffs <- predict(model_cv, as.matrix(train_data), type = "coefficients")
  print(coeffs)
  
  nonzero = which(abs(coeffs)>0)
  cell_counts[nonzero] = cell_counts[nonzero] + 1
  
}

names(cell_counts) = colnames(train_data)
cell_counts

## refit model
nonzero = which(cell_counts>=50) 

cSeed = 2001
print(cSeed)
set.seed(cSeed)

if (length(nonzero)>0){
  coxph_fit <- coxph(y_df ~ as.matrix(train_data[,nonzero])) # fit regular cox reg model, without lasso penalty
  coeffs = coef(coxph_fit)
  names(coeffs) = colnames(train_data[,nonzero])
}else{
  print("failed")
  return()
}
print(coeffs)
coeffs = coeffs[which(coeffs >= 2)]
print(coeffs)
write.csv(coeffs,'tumor_protein_coeffs_resistance_PFS.csv')

## hzr on training

coeffs = read.csv('tumor_protein_coeffs_resistance_PFS.csv')
coeffs1 = coeffs$x
names(coeffs1) = coeffs$X
coeffs = coeffs1
nonzero = which(colnames(train_data) %in% names(coeffs))

scores = matrix(0,dim(train_data)[1],1)
for(i in 1:dim(train_data)[1]){
  scores[i] = sum(as.matrix(train_data[i,nonzero])*t(abs(coeffs)))
}  

binary_score = (scores>=quantile(scores, probs = 0.66))

#y_df2 = Surv(as.numeric(cellFrac$PFS_5Years_months, cellFrac$PFS_5Years_Index)) # PFS-5years test
y_df2 = Surv(as.numeric(cellFrac$PFS_2Years_months, cellFrac$PFS_2Years_Index)) # PFS-2years test
train_data_2 = train_data
train_data_2$binary_score = binary_score
fit <- survfit(y_df2 ~ binary_score, data = train_data_2)

cox <- coxph(y_df2 ~ binary_score, data = train_data_2)
cox
hzr = exp(cox$coefficients)
hzr

# # Plot informative survival curves
A = ggsurvplot(fit, data = train_data_2,
               #title = "Tumor (Training Cohort:YTMA-471)",
               pval = FALSE, pval.method = FALSE,     # Add p-value &  method name
               surv.median.line = "hv",               # Add median survival lines
               legend.title = "Proliferating tumor+Granulocytes+Vessels",               # Change legend titles
               legend.labs = c("Low", "High"),        # Change legend labels
               palette = "lancet",                    # Use JCO journal color palette
               risk.table = TRUE,                    # Add No at risk table
               cumevents = FALSE,                     # Add cumulative No of events table
               tables.height = 0.25,                  # Specify tables height
               tables.theme = theme_cleantable(),     # Clean theme for tables
               tables.y.text = TRUE,
              xlab = "PFS-2Year (months)"   # Hide tables y axis text
)
A

cox <- coxph(y_df2 ~ binary_score, data = train_data_2)
summary(cox)
ggforest(cox)
A$plot <- A$plot+
  ggplot2::annotate("text",
                    x = 15, y = 0.75, # x and y coordinates of the text
                    label = "HR = 3.8(1.5-9.1) \n p = 0.004** (Log-Rank-2-sided) \n cutpoint = tertile", size = 4)
A

## hzr on test

cellFrac_test = read.csv("Tumor_Cell_Prop_4.csv")
cellFrac_test = cellFrac_test[which(cellFrac_test$TMA %in% c("4301")),]
# cellFrac_test = cellFrac_test[which(cellFrac_test$Histology.1 %in% c("Squamous")),]
# cellFrac_test = cellFrac_test[which(cellFrac_test$ImmunoTx %in% c("Nivo")),]
# cellFrac_test = cellFrac_test[which(cellFrac_test$ImmunoTx %in% c("Pembro")),]

y_df_test = Surv(as.numeric(cellFrac_test$PFS_2Years_months, cellFrac_test$PFS_2Years_Index)) 


test_data = cellFrac_test[,c(38:51)] 

scores = matrix(0,dim(test_data)[1],1)
for(i in 1:dim(test_data)[1]){
  scores[i] = sum(as.matrix(test_data[i,nonzero])*t(abs(coeffs)))
}  

binary_score = (scores>=quantile(scores, probs = 0.66))
test_data_2 = test_data
test_data_2$scores = scores # to test continuous scores
test_data_2$binary_score = binary_score
fit <- survfit(y_df_test ~ binary_score, data = test_data_2)

cox <- coxph(y_df_test ~ binary_score, data = test_data_2)
cox
hzr = exp(cox$coefficients)
hzr

# # Plot informative survival curves
B = ggsurvplot(fit, data = test_data_2,
               #title = "Tumor (Validation Cohort:UQ-4301)",
               pval = FALSE, pval.method = FALSE,              # Add p-value &  method name
               surv.median.line = "hv",                        # Add median survival lines
               legend.title = "Proliferating tumor+Granulocytes+Vessels",               # Change legend titles
               legend.labs = c("Low", "High"),                 # Change legend labels
               palette = "lancet",                             # Use JCO journal color palette
               risk.table = TRUE,                             # Add No at risk table
               cumevents = FALSE,                              # Add cumulative No of events table
               tables.height = 0.25,                           # Specify tables height
               tables.theme = theme_cleantable(),              # Clean theme for tables
               tables.y.text = TRUE,
              xlab = "PFS-2Years (months)"                     # Hide tables y axis text
)
B

cox <- coxph(y_df_test ~ binary_score, data = test_data_2)
cox

# Get two-sided p-value
p_value_two_sided <- summary(cox)$coefficients[5]

# Convert to one-sided p-value (assuming we are testing group1 > group0)
p_value_one_sided <- ifelse(summary(cox)$coefficients[2] > 1, p_value_two_sided / 2, 1 - p_value_two_sided / 2)
p_value_one_sided

ggforest(cox)
B$plot <- B$plot+
  ggplot2::annotate("text",
                    x = 9, y = 0.75, # x and y coordinates of the text
                    label = "HR = 1.8(0.87-3.7) \n p = 0.05* (Log-Rank-1-sided) \n cutpoint = tertile", size = 4)
B

