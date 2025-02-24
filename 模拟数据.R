
rm(list = ls())
library(tidyr)
library(dplyr)
get_simdata <- function(n_samples=1000,n_snps=100,par=c(miu=1,A=1,B=-1.2,ia=0.8,ib=-0.5,sd=1)){
#   par=c(miu=0.5,A=0.6,B=-0.8,ia=0.3,ib=-0.2)
# n_samples <- 1000
# n_snps <- 10  # 10 个 SNP

# 生成 SNP 数据（0/1 编码）
SNP_data <- as.data.frame(matrix(sample(0:1, n_samples * n_snps, replace = TRUE), 
                                 nrow = n_samples, ncol = n_snps))
colnames(SNP_data) <- paste0("SNP_", 1:n_snps)

# 定义 GWAS 显著 SNP 和 互作 SNP
SNP_A <- SNP_data$SNP_1  # 显著 GWAS SNP
SNP_B <- SNP_data$SNP_2  # 显著 GWAS SNP
SNP_X <- SNP_data$SNP_3  # 互作 SNP
SNP_Y <- SNP_data$SNP_4  # 互作 SNP
SNP_data <- SNP_data%>%
  mutate(
    Interaction_A_X = case_when(
      SNP_A == 0 & SNP_X== 0 ~ 1,  # (00)
      SNP_A == 1 & SNP_X == 1 ~ 1,  # (11)
      TRUE ~ 0                         # (01, 10)
    )
)
SNP_data <- SNP_data%>%
  mutate(
    Interaction_B_Y = case_when(
      SNP_B == 0 & SNP_Y== 0 ~ 1,  # (00)
      SNP_B == 1 & SNP_Y == 1 ~ 1,  # (11)
      TRUE ~ 0                         # (01, 10)
    )
  )
# SNP_data <- SNP_data%>%
#   mutate(
#     Interaction_A_X = 
#       SNP_A *SNP_X
#   )
# SNP_data <- SNP_data%>%
#   mutate(
#     Interaction_B_Y = 
#       SNP_B *SNP_Y
#   )


random_SNPs <- SNP_data[, 5:n_snps]  # 其他 6 个随机 SNP

# 生成 表型（包含 GWAS SNP + 互作效应 + 噪声）
interaction_term1 <- SNP_data$Interaction_A_X  # SNP_A 和 SNP_X 的互作
interaction_term2 <- SNP_data$Interaction_B_Y  # SNP_B 和 SNP_Y 的互作
miu=par[1]
# 计算线性得分
linear_score <- miu+par[2] * SNP_A - par[3] * SNP_B + par[4] * interaction_term1 - par[5] * interaction_term2

# 添加随机噪声
noise <- rnorm(n_samples, mean = 0, sd =par[6])
linear_score <- (linear_score + noise)
linear_score <- scale(linear_score)  # 归一化
#linear_score <- 2*(linear_score - min(linear_score)) / (max(linear_score) - min(linear_score))-1

prob <- 1 / (1 + exp(-linear_score))

Y <- rbinom(n_samples, 1, prob)  # 以 prob 作为概率生成 0/1 标签

# 合并数据
simulated_data <- cbind(SNP_data[,-c((ncol(SNP_data)-1):ncol(SNP_data))], Y)



return(simulated_data)}

simulated_data=get_simdata()

# 输出数据概览
head(simulated_data)


table(simulated_data$Y)
