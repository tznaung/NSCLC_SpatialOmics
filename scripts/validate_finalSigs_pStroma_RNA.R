library("survminer")
library("ggplot2")
require("survival")
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(edgeR)
library(pROC)
library(ggplot2)
library(verification)

pStroma = read.csv("UQ_Stroma_large_PFS_OS_only_1st_2nd_line.csv", row.names = 2)
#pStroma = read.csv("UQ_stroma_PFS_OS.csv", row.names = 2)

data_pStroma = pStroma[, 17:1843]
#data_pStroma = pStroma[, 3:1829]

###############################################################

patIds = rownames(data_pStroma)
nPat = length(patIds)

maxN_coeff = 10 # HR = 0.305, p = 0.255 

nSig = 1

hzrs_all = matrix(0,maxN_coeff,nSig)
pvals_all = matrix(0,maxN_coeff,nSig)

nLargestCoeffs = maxN_coeff

print(nLargestCoeffs)

sigScores = matrix(0,nPat,nSig)
sigScores = matrix(0,nPat,nSig)

for (cSig in 1:nSig){
  
  if (cSig==1){
    coeffs = read.csv("final_coeffs_PFS_5Years_3celltypes_Stroma_glmnet_prop0p05_CW.csv")
    c_data = data_pStroma
    resps = pStroma$PFS_2Years_index
  }
  if (cSig==2){
  }
  if (cSig==3){
  }
  if (cSig==4){
  }
  if (cSig==5){
  }
  
  c_nPat = dim(c_data)[1]
  c_resps = rep(0,c_nPat)
  for (i in 1:c_nPat){
    idx = which(rownames(c_data)[i] == patIds)
    c_resps[i] = resps[idx]
  }
  
  genes = coeffs$X
  nGene = length(genes)
  cfs = coeffs$V1
  print(cSig)
  print("before:")
  print(genes)

  ### sort and step down
  cfs_sort = sort(abs(cfs),decreasing = TRUE)
  if ((min(nLargestCoeffs,nGene)+1) <= nGene){
    cfs_sort = cfs_sort[(min(nLargestCoeffs,nGene)+1):nGene]
    cf_idx = which(abs(cfs) %in% cfs_sort)
    cfs[cf_idx] = 0
  }
  
  names(cfs) = genes
  print(cfs)

  for (cPat in 1:nPat){
    vals = rep(0,nGene)
    pat_idx = which(rownames(c_data)==patIds[cPat])
    if(length(pat_idx)>0){
      for (i in 1:nGene){
        idx = which(colnames(c_data)==genes[i])
        if (length(idx)>0){
          if (is.finite(c_data[pat_idx,idx])){
            vals[i] = c_data[pat_idx,idx]
          }
        }else{
          print("not found")
          print(genes[i])
        }
      }
      sigScores[cPat,cSig] = sum((as.matrix(vals))*(as.matrix(cfs)))
    }else{
      sigScores[cPat,cSig] = NA
    }
  } 
}

#y_df_test = Surv(pStroma$PFS_2Years_months, pStroma$PFS_2Years_index) ###
y_df_test = Surv(pStroma$OS_5Years_months, pStroma$OS_5Years_index) ###

pStroma$binary_score = (sigScores[,1]>=quantile(sigScores[,1],probs=0.66))
#pStroma$binary_score = (sigScores[,1]>=quantile(sigScores[,1],probs=0.67))
cox <- coxph(y_df_test ~ binary_score, data = pStroma)

hzr = exp(cox$coefficients)
pval = summary(cox)$waldtest[3]

# Plot

fit <- survfit(y_df_test ~ binary_score, data = pStroma)
fit$n.risk

col = c("deeppink","slateblue1")
# # Plot informative survival curves
A = ggsurvplot(fit, data = pStroma,
               #title = "Stroma (Validation Cohort:UQ-4301)",
               pval = FALSE, pval.method = FALSE,    # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "sig scores (response)",               # Change legend titles
               legend.labs = c("Low", "High"),  # Change legend labels
               #palette = "Dark2",                    # Use JCO journal color palette
               palette = "jco",                    # Use JCO journal color palette
               risk.table = TRUE,                  # Add No at risk table
               cumevents = FALSE,                  # Add cumulative No of events table
               tables.height = 0.25,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE,              # Hide tables y axis text
               xlab = "OS (months)"
)
A

print(hzr)
# 
print(pval)

cox <- coxph(y_df_test ~ binary_score, data = pStroma)
cox

# find the hazard ratio with this function.
cox <- coxph(y_df_test ~ binary_score, data = pStroma)
cox

# Get two-sided p-value
p_value_two_sided <- summary(cox)$coefficients[5]

# Convert to one-sided p-value (assuming we are testing group1 > group0)
p_value_one_sided <- ifelse(summary(cox)$coefficients[2] < 1, p_value_two_sided / 2, 1 - p_value_two_sided / 2)
p_value_one_sided

ggforest(cox)
A$plot <- A$plot+
  ggplot2::annotate("text",
                    x = 15, y = 0.75, # x and y coordinates of the text
                    label = "HR = 1.1 (0.43-3) \n p = 0.6 (Log-Rank-1-sided) \n cutpoint = tertile", size = 4)
A

