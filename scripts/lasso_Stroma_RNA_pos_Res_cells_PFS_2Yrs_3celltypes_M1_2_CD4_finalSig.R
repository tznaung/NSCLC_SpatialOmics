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


cell_genes1 = read.csv("Macrophages.M1_cell_genes.csv")
cell_genes2 = read.csv("Macrophages.M2_cell_genes.csv")
cell_genes3 = read.csv("CD4_cell_genes.csv")
cell_genes4 = read.csv("Macrophage_cell_genes_ChinWee.csv")
cell_genes5 = read.csv("CD4_cell_genes_ChinWee.csv")

cell_genes1 = cell_genes1$x
cell_genes2 = cell_genes2$x
cell_genes3 = cell_genes3$x
cell_genes4 = cell_genes4$x
cell_genes5 = cell_genes5$x

cell_genes = union(cell_genes1,cell_genes2)
cell_genes = union(cell_genes,cell_genes3)
cell_genes = union(cell_genes,cell_genes4)
cell_genes = union(cell_genes,cell_genes5)

UQ_Stroma = read.csv("UQ_Stroma_large_PFS_OS_only_1st_2nd_line.csv", row.names = 2)
cell_genes = intersect(cell_genes, colnames(UQ_Stroma))

# get Sig DE genes for training set 

Stroma_df = read.csv("Stroma_av_3B.csv", header = T)
Stroma_df$Spot_ID <- Stroma_df$X
Res_Stroma = read.csv("Yale_firstLine_pre_patients.csv", header = T)
Stroma_2 = merge(Stroma_df, Res_Stroma, by = "Spot_ID", all = FALSE)
Stroma_2$group <- Stroma_2$PFS_5Years_Index # define the group you want to compare
data = Stroma_2[,3:6095]
prop_test = 0.20

idx0 = which(Stroma_2$group==1) 
idx1 = which(Stroma_2$group==0)

hzrs = read.csv("hzrs_Stroma_PFS_5Years_3celltypes_response_CW.csv")

#good_models = which((hzrs$x)<=0.7 & (hzrs$x)>0)
good_models = which((hzrs$x)<=0.7)
#good_models = which((hzrs$x)>=1.5)
nGoodModels = length(good_models)
goodModelGeneCounts = rep(0,dim(data)[2])
names(goodModelGeneCounts) = colnames(data)

for (cSplit in good_models){

  print(cSplit)
  set.seed(cSplit)
  idx0_shuffled = sample(idx0, length(idx0), replace = FALSE)
  idx1_shuffled = sample(idx1, length(idx1), replace = FALSE)
  
  test_idx = c(idx0_shuffled[1:(length(idx0)*prop_test)],idx1_shuffled[1:(length(idx1)*prop_test)])
  train_idx = setdiff(1:dim(data)[1],test_idx)
  train_data = Stroma_2[train_idx,]
  train_data2 = train_data[,3:6095]
  test_data = Stroma_2[test_idx,]
  test_data2 = test_data[,3:6095]
  
  ## ADD
  train_data2 = train_data2[,which(colnames(train_data2) %in% cell_genes)]
  test_data2 = test_data2[,which(colnames(test_data2) %in% cell_genes)]
  
  # make the pseudo counts of the data
  
  d = DGEList(t(train_data2), group= train_data$PFS_5Years_Index)
  d = calcNormFactors(d)    
  d = estimateCommonDisp(d)
  d = estimateTagwiseDisp(d)
  
  #Find significant genes
  response <- exactTest(d, pair = c("1", "0"))
  #FDR <- p.adjust(response$table$PValue, method="BH")
  FDR <- response$table$PValue
  sum(FDR < 0.05)
  topTags(response)
  plotMD(response)
  abline(h=c(-1,1), col="blue")
  response$table$FDR = FDR
  #write.csv(response$table, "Stroma_train_set_Sig_Genes_PFS.csv")
  
  dex_df = response$table
  dex_sel_genes = dex_df[which(dex_df$FDR < 0.05),] # take only genes that have FDR < 0.1
  dex_df$gene = row.names(dex_df)
  do_lasso = 1
  
  if (dim(dex_sel_genes)[1]<=1){
    next
  }
  
  train_data3 = train_data[,which(colnames(train_data) %in% rownames(dex_sel_genes))]
  # train_data3$response = train_data$response
  # y_df = Surv(train_data$PFS_Days, train_data$PFS_Index) ###
  y_df = Surv(train_data$PFS_5Years_months, train_data$PFS_5Years_Index) ###
  
  if (do_lasso){
    nGene = dim(train_data3)[2] - 1
    dim2 = dim(train_data3)[2]
    var_counts = matrix(0,1,nGene)
    mean_coeffs = matrix(0,nGene,1)
    nSeeds = 100
    
    for(cSeed in 1:nSeeds){
      
      print(cSeed)
      set.seed(cSeed)
      
      # lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
      lambdas_to_try <- 10^seq(-3, 2, length.out = 100)
      ##
      lasso_cv <- cv.glmnet(as.matrix(train_data3[,1:nGene]), y_df, alpha = 0.1, lambda = lambdas_to_try,
                            standardize = FALSE, nfolds = 10, family = "cox", upper.limits = 0) ###
      #plot(lasso_cv)
      lambda_cv <- lasso_cv$lambda.min
      model_cv <- glmnet(as.matrix(train_data3[,1:nGene]), y_df, alpha = 0.1, lambda = lambda_cv, standardize = FALSE, family = "cox", upper.limits = 0) ###
      coeffs <- predict(model_cv, as.matrix(train_data3[,1:nGene]), type = "coefficients")
      
      var_counts[which(coeffs!=0)] = var_counts[which(coeffs!=0)] + 1
      mean_coeffs = mean_coeffs + (coeffs/nSeeds)
      
    }
    
    mean_coeffs_weighted = mean_coeffs * (nSeeds/t(var_counts));
    
    var_props = var_counts / nSeeds
    colnames(var_props) = colnames(train_data3)[1:nGene]
    rownames(mean_coeffs) = colnames(train_data3)[1:nGene]
    rownames(mean_coeffs_weighted) = colnames(train_data3)[1:nGene]
    t(var_props)
    #mean_coeffs
    mean_coeffs_weighted #[2:nGene+1,]
    
    ## refit model
    nonzero = which(var_props>0)
    
    cSeed = 2001
    print(cSeed)
    set.seed(cSeed)
    
    if (length(nonzero)>1){
      #coxph_fit <- coxph(y_df ~ as.matrix(train_data3[,nonzero])) # fit regular cox reg model, without lasso penalty
      #coeffs = coef(coxph_fit)
      #names(coeffs) = colnames(train_data3[,nonzero])
      
      model_cv <- glmnet(as.matrix(train_data3[,nonzero]), y_df, alpha = 0.1, lambda = 0, standardize = FALSE, family = "cox", upper.limits = 0) ###
      coeffs <- predict(model_cv, as.matrix(train_data3[,nonzero]), type = "coefficients")  
      coeffs = coeffs[,1]
      
    }else{
      print("failed")
      next     
    }
  }else{
    nGene = dim(train_data3)[2] - 1
    dim2 = dim(train_data3)[2]  
    nonzero = (1:nGene)
  }
  
  print(coeffs)
  
  # Testing set 
  test_data = Stroma_2[test_idx,]
  
  test_data3 = test_data[,which(colnames(test_data) %in% rownames(dex_sel_genes))]
  # test_data3$response = test_data$response  
  # y_df_test = Surv(test_data$PFS_Days, test_data$PFS_Index) ###
  y_df_test = Surv(test_data$PFS_5Years_months, test_data$PFS_5Years_Index) ###
  test_data3$lasso_scores = matrix(0,dim(test_data3)[1],1)
  for(i in 1:dim(test_data3)[1]){
    test_data3$lasso_scores[i] = sum(as.matrix(test_data3[i,nonzero])*t(coeffs))
  }  
  

  # test_data3$binary_score = (test_data3$lasso_scores>=median(test_data3$lasso_scores))
  test_data3$binary_score = (test_data3$lasso_scores>=quantile(test_data3$lasso_scores,probs=0.66))
  # test_data3$binary_resp = (test_data3$response == "yes")
  cox <- coxph(y_df_test ~ binary_score, data = test_data3)
  print(c(exp(cox$coefficients),hzrs$x[cSplit]))
  print(names(coeffs))
  
  gene_idxs = which(names(goodModelGeneCounts) %in% names(coeffs))
  goodModelGeneCounts[gene_idxs] = goodModelGeneCounts[gene_idxs] + 1
  
  computed_hr = exp(cox$coefficients[1])  # 
  expected_hr = hzrs$x[cSplit]
  
  print(paste("Computed HR: ", computed_hr))
  print(paste("Expected HR from hzrs: ", expected_hr))
  
  # Check if they are the same (within a small numerical tolerance)
  if (abs(computed_hr - expected_hr) < .Machine$double.eps^0.5) {
    print("Values are the same within numerical tolerance.")
  } else {
    print("Values are different.")
  }

}

write.csv(goodModelGeneCounts,"goodModelCounts_Stroma_PFS_5Years_3celltypes_CW.csv")

goodModelGeneCounts = read.csv("goodModelCounts_Stroma_PFS_5Years_3celltypes_CW.csv",row.names = 1)
goodModelGeneCounts_ = as.matrix(goodModelGeneCounts)
names(goodModelGeneCounts_) = row.names(goodModelGeneCounts)
goodModelGeneCounts = goodModelGeneCounts_

goodModelGeneProps = goodModelGeneCounts / nGoodModels
names(goodModelGeneProps) = row.names(goodModelGeneCounts)
goodModelGeneProps = goodModelGeneProps[goodModelGeneProps>0]
sort(goodModelGeneProps)

final_genes = names(goodModelGeneProps)[which(goodModelGeneProps>0.05)]
#final_genes = names(goodModelGeneProps)[which(goodModelGeneProps>0.01)]

train_data = Stroma_2
train_data2 = train_data[,3:6095]
train_data3 = train_data[,which(colnames(train_data) %in% final_genes)]
train_data3$response = train_data$response
y_df = Surv(train_data$PFS_5Years_months, train_data$PFS_5Years_Index) ###
# y_df = Surv(train_data$OS_Days_2Yrs, train_data$OS_Index_2Yrs) ###

cSeed = 4001
print(cSeed)
set.seed(cSeed)

model_cv <- glmnet(as.matrix(train_data3[,1:length(final_genes)]), y_df, alpha = 0, lambda = 0, standardize = FALSE, family = "cox", upper.limits = 0) ###
coeffs <- predict(model_cv, as.matrix(train_data3[,1:length(final_genes)]), type = "coefficients")  
final_coeffs = coeffs[,1]
names(final_coeffs) = rownames(coeffs)
coeffs

#coxph_fit <- coxph(y_df ~ as.matrix(train_data3[,1:length(final_genes)])) # fit regular cox reg model, without lasso penalty
#final_coeffs = coef(coxph_fit)
#names(final_coeffs) = colnames(train_data3[,1:length(final_genes)])

#print(final_coeffs)
#final_coeffs = final_coeffs[which(final_coeffs <= 0)]
#print(final_coeffs)
#final_coeffs = final_coeffs[which(final_coeffs <= -0.0002)]
#print(final_coeffs)
write.csv(as.matrix(final_coeffs),"final_coeffs_PFS_5Years_3celltypes_Stroma_glmnet_prop0p05_CW.csv")

## hzr on training--------------------------------------------------------------

Stroma_df = read.csv("Stroma_av_3B.csv", header = T)
Stroma_df$Spot_ID <- Stroma_df$X
Res_Stroma = read.csv("Yale_firstLine_pre_patients.csv", header = T)
Stroma_2 = merge(Stroma_df, Res_Stroma, by = "Spot_ID", all = FALSE)
Stroma_2$group <- Stroma_2$PFS_5Years_Index # define the group you want to compare
data = Stroma_2[,3:6095]
prop_test = 0.20
train_data = Stroma_2

coeffs = read.csv('final_coeffs_PFS_5Years_3celltypes_Stroma_glmnet_prop0p05_CW.csv')
coeffs1 = coeffs$V1
names(coeffs1) = coeffs$X
coeffs = coeffs1
nonzero = which(colnames(train_data) %in% names(coeffs))

scores = matrix(0,dim(train_data)[1],1)
for(i in 1:dim(train_data)[1]){
  scores[i] = sum(as.matrix(train_data[i,nonzero])*t(abs(coeffs)))
}  


# binary_score = (scores>=median(scores))
binary_score = (scores>=quantile(scores, probs = 0.66))

#y_df2 = Surv(train_data$PFS_2Years_months, train_data$PFS_2Years_Index)
#y_df2 = Surv(train_data$PFS_5Years_months, train_data$PFS_5Years_Index)
#y_df2 = Surv(train_data$OS_2Years_months, train_data$OS_2Years_Index)
y_df2 = Surv(train_data$OS_5Years_months, train_data$OS_5Years_Index)

train_data_2 = train_data
train_data_2$binary_score = binary_score
fit <- survfit(y_df2 ~ binary_score, data = train_data_2)

cox <- coxph(y_df2 ~ binary_score, data = train_data_2)
cox
hzr = exp(cox$coefficients)
hzr

# # Plot informative survival curves
A = ggsurvplot(fit, data = train_data_2,
               # title = "Tumor (Training Cohort)",
               # title = "Stroma (Training Cohort:YTMA-471)",
               pval = FALSE, pval.method = FALSE,    # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "sig scores (response)",               # Change legend titles
               #legend.title = "M-1 + M-2 + CD4",
               legend.labs = c("Low", "High"),  # Change legend labels
               palette = "jco",                    # Use JCO journal color palette
               risk.table = TRUE,                  # Add No at risk table
               cumevents = FALSE,                  # Add cumulative No of events table
               tables.height = 0.25,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = TRUE,
               xlab = "OS-5Year (months)"  # Hide tables y axis text
)
A

# find the hazard ratio with this function.
cox <- coxph(y_df2 ~ binary_score, data = train_data_2)
summary(cox)
ggforest(cox)
A$plot <- A$plot+
  ggplot2::annotate("text",
                    x = 40, y = 0.5, # x and y coordinates of the text
                    label = "HR = 0.34 (0.096-1.2) \n p =0.06 (Log-Rank-2-sided) \n cutpoint = tertile", size = 4)
A



## heatmap

# pdf(file = paste("output_split_Stroma_final.pdf", sep=""))
# 
# # Training set 
# train_data4 = train_data3[,1:length(final_genes)]
# train_data4 <- log2(train_data4)
# rownames(train_data4) = train_data$X
# #  sapply(rownames(train_data4),function(x) 
# #  strsplit(as.character(x),split = "\\\\")[[1]][1])
# train_data4 <- train_data4[order(train_data$response),]
# train_data5 = data.frame("Response" = sort(train_data$response))
# rownames(train_data5) = rownames(train_data4)
# #rownames(train_data4) = train_data5$Response # name matching
# pheatmap(train_data4,annotation_row = train_data5, main = "Signature genes Stroma Comp",
#          cluster_rows = T, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
#          # cutree_cols = 2,
#          cutree_rows = 2, 
#          fontsize = 8)
# pheatmap(train_data4,annotation_row = train_data5, main = "Signature genes for CK comp",
#          cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
#          cutree_cols = 2,
#          #cutree_rows = 3, 
#          fontsize = 8)
# 
# 
# #dev.off()



