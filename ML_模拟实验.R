library(xgboost)
library(shapviz)
library(pROC)
library(caret)
library(tidyverse)
library(doParallel)

# 启用并行计算
registerDoParallel(cores = 4)

# 定义交叉验证函数

run_kfold <- function(data, features, target, params, nfold=5) {
  folds <- createFolds(data[[target]], k = nfold)
  
  results <- foreach(i = 1:nfold, 
                     .combine = rbind, 
                     .packages = c("xgboost", "pROC")) %dopar% {  # 传递包到子进程
                       # 划分训练测试集
                       train_data <- data[-folds[[i]], ]
                       test_data <- data[folds[[i]], ]
                       
                       # 转换数据格式（此时xgb.DMatrix可用）
                       dtrain <- xgb.DMatrix(
                         data = as.matrix(train_data[, features]),
                         label = train_data[[target]]
                       )
                       dtest <- xgb.DMatrix(
                         data = as.matrix(test_data[, features]),
                         label = test_data[[target]]
                       )
                       
                       # 训练模型
                       model <- xgb.train(
                         params = params,
                         data = dtrain,
                         nrounds = 100,
                         # early_stopping_rounds = 10,
                         # watchlist = list(test = dtest),
                         verbose = 0
                       )
                       
                       # 预测并计算AUC
                       pred <- predict(model, dtest)
                       auc_val <- auc(response = test_data[[target]], predictor = pred)
                       
                       list(fold = i, auc = auc_val,model=model)
                     }
  
  return(results)
}
# 定义各方案特征集
define_features <- function(data) {
  # 方案1: 所有SNP
  scheme1 <- grep("^SNP_", colnames(data), value = TRUE)
  
  # 方案2: 显著SNP (假设已知SNP1和SNP2)
  scheme2 <- c("SNP_1", "SNP_2")
  
  # 方案3: 显著SNP + 互作SNP (假设SNP3和SNP4)
  scheme3 <- c(scheme2, "SNP_3", "SNP_4")
  
  # 方案4: 方案3 + 显式交互项
  scheme4 <- c(
    scheme3,
    "Interaction1_3" = "Interaction1_3",
    "Interaction2_4" = "Interaction2_4"
  )
  
  list(
    scheme1 = scheme1,
    scheme2 = scheme2,
    scheme3 = scheme3,
    scheme4 = scheme4
  )
}

# 生成带交互项的数据
data_prep <- function(data) {
  data %>%
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
}
data_prep <- function(data) {
  data %>%
    mutate(
      SNP_1 = as.numeric(SNP_1),
      SNP_2 = as.numeric(SNP_2),
      SNP_3 = as.numeric(SNP_3),
      SNP_4 = as.numeric(SNP_4),
      # 显式构建交互项（方案4专用）
      "Interaction1_3" = SNP_1 * SNP_3,
      "Interaction2_4"= SNP_2 * SNP_4
    )
}


simulated_data=get_simdata(n_samples=1000,n_snps=100000,c(miu=0.3,A=0.8,B=-0.9,ia=0.9,ib=-0.8, sd=0.1))

prepped_data <- data_prep(simulated_data)


feature_sets <- define_features(prepped_data)

# 模型参数
params <- list(
  objective = "binary:logistic",
  max_depth = 5,
  eta = 0.3,
  eval_metric="auc"
  # subsample =0.8,
  # colsample_bytree = 0.5,
  # reg_alpha = 1,
  # reg_lambda = 1
)

# 并行运行各方案
results <- list()
for (scheme in names(feature_sets)) {
  cat("model:", scheme, "\n")
  res <- run_kfold(
    data = prepped_data,
    features = feature_sets[[scheme]],
    target = "Y",
    params = params
  )
  results[[scheme]] <- res
}

# 计算各方案性能指标
performance <- sapply(results,function(res){
 nres <-  as.numeric(res[,2])
  return(c(mean_auc= mean(nres),sd_auc= sd(nres)))
})
performance
# 可视化比较
library(ggplot2)
library(tidyr)
library(dplyr)
df_auc <- as.data.frame(performance)
colnames(df_auc) <- c("scheme1", "scheme2", "scheme3", "scheme4")
df_auc$metric <- rownames(performance)

df_long <- df_auc %>%
  pivot_longer(cols = starts_with("scheme"), 
               names_to = "scheme", 
               values_to = "value") %>%
  pivot_wider(names_from = metric, values_from = value)

df_long$design <- factor(df_long$scheme, levels = c("scheme1", "scheme2", "scheme3", "scheme4"),
                         labels = c("Whole Genome", 
                                    "GWAS Significant Loci", 
                                    "GWAS + Epistatic Interactions", 
                                    "GWAS + Epistasis + Constructed Interactions"))



pdf("auc_plot.pdf",width = 10,height = 6)
ggplot((df_long), aes(x = scheme, y = mean_auc, fill = design)) +
  scale_fill_manual(values = c("darkgreen", "royalblue3","brown2","orange"))+
  geom_bar(stat = "identity", width = 0.8) +  # 绘制柱状图
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), width = 0.2) +  # 添加误差条
  labs(x = NULL, y = NULL)  +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
   # 不显示图例
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )+  coord_cartesian(ylim = c(0.5, 0.8))+  geom_text(aes(label = round(mean_auc, 3)), vjust = -0.5, nudge_y = 0.03,size = 5) + 
  scale_x_discrete(breaks=NULL,labels=NULL)+
  scale_y_continuous(breaks = seq(0.5,0.8,0.1),position = "right")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()
  )+theme(axis.text.x = element_text( color="black", size=18),axis.text.y =element_text( color="black", size=18),
          panel.border = element_rect(color = "black", size = 1.5, fill = NA),
          axis.ticks.length.y = unit(0.4,"cm"),axis.ticks.length.x = unit(0.4,"cm") ,axis.title.x = element_text(size = 18),  # 修改横坐标标题的字体大小
          axis.title.y = element_text(size = 18))
dev.off()



# 确保正确生成特征矩阵
X_pred <- as.matrix(prepped_data[, feature_sets$scheme3])
colnames(X_pred) <- feature_sets$scheme3  # 保留特征名

# 重新训练模型（确保特征名传递）
dtrain <- xgb.DMatrix(X_pred, label = prepped_data$Y)
xgb_model <- xgb.train(params, dtrain, nrounds = 100)

# 正确调用shapviz（使用X_pred参数）
shap_inter <- shapviz(
  xgb_model, 
  X_pred = X_pred,  # 必须传递矩阵，而非数据框
  interactions = TRUE
)

# 提取交互值矩阵
inter_values <- get_shap_interactions(shap_inter)

# 验证交互值维度
print(dim(inter_values))  # 应输出 [n_samples, n_features, n_features]
# 检查矩阵特征名和维度
if ("SNP_1" %in% colnames(X_pred) && "SNP_3" %in% colnames(X_pred)) {
  avg_inter <- mean(inter_values[, "SNP_1", "SNP_3"])
  cat("SNP1-SNP3平均交互强度:", avg_inter, "\n")
} else {
  cat("未找到预设互作特征，当前特征名为:", colnames(X_pred))
}
