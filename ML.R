library(readxl)
library(MuMIn)
library(tidyverse)
library(randomForest)
library(varSelRF)
library(pROC)
library(ggplot2)
library(caret)
library(rfPermute)
library(neuralnet)
library(lightgbm)
library(xgboost)
library(shapviz)
library(PRROC)



simulated_data=get_simdata(n_samples=10000,n_snps=100000,c(miu=0.3,A=0.8,B=-0.9,ia=0.9,ib=-0.8, sd=0.1))

new_dat <- as.data.frame(simulated_data)
train <- sample(nrow(new_dat), nrow(new_dat)*0.8)

traindata <- new_dat[train,]
testdata <- new_dat[-train,]
# traindata <- xgb.DMatrix(data = as.matrix(traindata[,-1]),label=traindata$Subtype)
# testdata <- xgb.DMatrix(data = as.matrix(testdata[,-1]),label=testdata$Subtype)
traindata <- list(data = as.matrix(traindata[,-ncol(new_dat)]),label=traindata$Y)
testdata <- list(data = as.matrix(testdata[,-ncol(new_dat)]),label=testdata$Y)

res_xgboost = xgboost(
  data = traindata$data,
  label = traindata$label,
  max_depth = 5, 
  eta = 0.3, 
  
  nrounds = 100,
  verbose = 0,
  objective = "binary:logistic",
  eval_metric="auc")
xgb_pred <- predict(res_xgboost,testdata$data)

# pr_xgboost<- pr.curve(scores.class0 = xgb_pred[testdata$label == 1],
#                   scores.class1 = xgb_pred[testdata$label == 0],
#                   curve = TRUE)
pr_xgboost<-pr.curve(scores.class0 = xgb_pred,
                  weights.class0 = testdata$label,
                  curve = TRUE)

roc_xgboost1<-roc(as.ordered(testdata$label) ,as.ordered(xgb_pred),smoth=TRUE)
plot(roc_xgboost1)


new_dat <- as.data.frame(simulated_data)
train <- sample(nrow(new_dat), nrow(new_dat)*0.8)

traindata <- new_dat[train,]
testdata <- new_dat[-train,]
# traindata <- xgb.DMatrix(data = as.matrix(traindata[,-1]),label=traindata$Subtype)
# testdata <- xgb.DMatrix(data = as.matrix(testdata[,-1]),label=testdata$Subtype)
traindata <- list(data = as.matrix(traindata[,1:2]),label=traindata$Y)
testdata <- list(data = as.matrix(testdata[,1:2]),label=testdata$Y)

res_xgboost = xgboost(
  data = traindata$data,
  label = traindata$label,
  max_depth = 5, 
  eta = 0.3, 
  
  nrounds = 100,
  verbose = 0,
  objective = "binary:logistic",
  eval_metric="auc")
xgb_pred <- predict(res_xgboost,testdata$data)

# pr_xgboost<- pr.curve(scores.class0 = xgb_pred[testdata$label == 1],
#                   scores.class1 = xgb_pred[testdata$label == 0],
#                   curve = TRUE)
pr_xgboost<-pr.curve(scores.class0 = xgb_pred,
                     weights.class0 = testdata$label,
                     curve = TRUE)

roc_xgboost2<-roc(as.ordered(testdata$label) ,as.ordered(xgb_pred),smoth=TRUE)
roc_xgboost2


new_dat <- as.data.frame(simulated_data)
train <- sample(nrow(new_dat), nrow(new_dat)*0.8)

traindata <- new_dat[train,]
testdata <- new_dat[-train,]
# traindata <- xgb.DMatrix(data = as.matrix(traindata[,-1]),label=traindata$Subtype)
# testdata <- xgb.DMatrix(data = as.matrix(testdata[,-1]),label=testdata$Subtype)
traindata <- list(data = as.matrix(traindata[,1:4]),label=traindata$Y)
testdata <- list(data = as.matrix(testdata[,1:4]),label=testdata$Y)

res_xgboost = xgboost(
  data = traindata$data,
  label = traindata$label,
  max_depth = 5, 
  eta = 0.3, 
  
  nrounds = 100,
  verbose = 0,
  objective = "binary:logistic",
  eval_metric="auc")
xgb_pred <- predict(res_xgboost,testdata$data)

# pr_xgboost<- pr.curve(scores.class0 = xgb_pred[testdata$label == 1],
#                   scores.class1 = xgb_pred[testdata$label == 0],
#                   curve = TRUE)
pr_xgboost<-pr.curve(scores.class0 = xgb_pred,
                     weights.class0 = testdata$label,
                     curve = TRUE)

roc_xgboost3<-roc(as.ordered(testdata$label) ,as.ordered(xgb_pred),smoth=TRUE)
roc_xgboost3




new_dat <- as.data.frame(simulated_data)

new_dat <- new_dat  %>%
  mutate(
    SNP_1 = as.numeric(SNP_1),
    SNP_2 = as.numeric(SNP_2),
    SNP_3 = as.numeric(SNP_3),
    SNP_4 = as.numeric(SNP_4),
    # 显式构建交互项（方案4专用）
    "Interaction1_3" = case_when(
      SNP_1 == 0 & SNP_3== 0 ~ 1,  # (00)
      SNP_1 == 1 & SNP_3 == 1 ~ 1,  # (11)
      TRUE ~ 0                         # (01, 10)
    ),
    "Interaction2_4"= case_when(
      SNP_2 == 0 & SNP_4== 0 ~ 1,  # (00)
      SNP_2 == 1 & SNP_4 == 1 ~ 1,  # (11)
      TRUE ~ 0                         # (01, 10)
    ))
train <- sample(nrow(new_dat), nrow(new_dat)*0.8)

traindata <- new_dat[train,]
testdata <- new_dat[-train,]
# traindata <- xgb.DMatrix(data = as.matrix(traindata[,-1]),label=traindata$Subtype)
# testdata <- xgb.DMatrix(data = as.matrix(testdata[,-1]),label=testdata$Subtype)
traindata <- list(data = as.matrix(traindata[,c(1:4,100002:100003)]),label=traindata$Y)
testdata <- list(data = as.matrix(testdata[,c(1:4,100002:100003)]),label=testdata$Y)

res_xgboost = xgboost(
  data = traindata$data,
  label = traindata$label,
  max_depth = 5, 
  eta = 0.3, 
  
  nrounds = 100,
  verbose = 0,
  objective = "binary:logistic",
  eval_metric="auc")
xgb_pred <- predict(res_xgboost,testdata$data)

# pr_xgboost<- pr.curve(scores.class0 = xgb_pred[testdata$label == 1],
#                   scores.class1 = xgb_pred[testdata$label == 0],
#                   curve = TRUE)
pr_xgboost<-pr.curve(scores.class0 = xgb_pred,
                     weights.class0 = testdata$label,
                     curve = TRUE)

roc_xgboost4<-roc(as.ordered(testdata$label) ,as.ordered(xgb_pred),smoth=TRUE)

plot(roc_xgboost4)

c(roc_xgboost1$auc,roc_xgboost2$auc,roc_xgboost3$auc,roc_xgboost4$auc)


roc_data <- data.frame(
  fpr = c(1 - roc_xgboost1$specificities,
          1 - roc_xgboost2$specificities,
          1 - roc_xgboost3$specificities,
          1 - roc_xgboost4$specificities),
  tpr = c(roc_xgboost1$sensitivities,
          roc_xgboost2$sensitivities,
          roc_xgboost3$sensitivities,
          roc_xgboost4$sensitivities),
  model = rep(c("Model1", "Model2","Model3","Model4"), 
              c(length(roc_xgboost1$sensitivities),
                length(roc_xgboost2$sensitivities),
                length(roc_xgboost3$sensitivities),
                length(roc_xgboost4$sensitivities)))
)

# 设置因子水平（重要！）
roc_data$model <- factor(roc_data$model, 
                         levels = c("Model1","Model2","Model3","Model4"))

# 生成AUC标签（确保顺序与因子水平一致）
auc_labels <- c(
  paste0("Model1 (AUC = ", round(auc_xgboost1, 4), ")"),
  paste0("Model2 (AUC = ", round(auc_xgboost2, 4), ")"),
  paste0("Model3 (AUC = ", round(auc_xgboost3, 4), ")"),
  paste0("Model4 (AUC = ", round(auc_xgboost4, 4), ")")
)

# 绘图代码修正
tt <- ggplot(roc_data, aes(x = fpr, y = tpr, color = model)) +
  geom_line(size = 1.5, alpha = 0.8) +  # 调整alpha增加可见性
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("darkgreen", "royalblue3","brown2","orange"),
                     labels = auc_labels) +
  labs(x = NULL,#"False Positive Rate", 
       y = NULL,#"True Positive Rate",
       title = NULL)+#"ROC Curve Comparison") +
  theme_bw() +
  theme(
    text = element_text(size = 14, family = "serif"),
    legend.position = c(0.75, 0.2),  # 调整图例位置
    legend.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", size = 1)
  ) +
  coord_equal() +

  # scale_x_continuous(expand = c(0,0),breaks  =seq(0,1,0.25),limits = c(0,1))+
  # scale_y_continuous(expand = c(0,0),breaks  =seq(0,1,0.25),limits = c(0,1))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill ="transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(#face="bold",
          color="black", size=16),
        axis.text.y =element_text(#face="bold",
          color="black", size=16),plot.margin=unit(c(1,1,1,1),"cm"))+
  theme(axis.ticks.length.y = unit(0.4,"cm"))+    theme(axis.ticks.length.x = unit(0.4,"cm"))+
  theme(text=element_text(size=16,  family="serif"),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid")) + theme(legend.title = element_blank(), legend.position = c(0.65, 0.25)) # 确保坐标轴比例一致
  pdf("roc_PLOT.pdf",width = 6,height = 6)
  print(tt)
  
  dev.off()
