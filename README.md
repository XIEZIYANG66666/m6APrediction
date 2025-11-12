---
title: "README.Rmd"
output: html_document
date: "2025-11-11"
---

# m6APrediction

`m6APrediction` is an R package based on machine learning for predicting m6A methylation sites in RNA sequences. This package includes functions such as model training, site prediction, and performance evaluation (e.g., drawing ROC and PRC curves), and supports common algorithms like random forest, helping researchers efficiently identify potential m6A modification sites.

# Installation

You can install this package from GitHub. Use the following code to install:

```{r}
# Install the devtools package (if you haven't installed it yet) install.packages("devtools")

# Install the m6APrediction package using devtools
devtools::install_github("[XIEZIYANG66666]/m6APrediction")
```

If you use the remotes package, you can use the following command:

```{r}
remotes::install_github("[XIEZIYANG66666]/m6APrediction")
```

#Usage example 1. Data loading and preprocessing Load the dataset and perform feature engineering (take m6A_unbalanced.csv as an example) :

```{r}
# Loading package
data_df <- read.csv("m6A_unbalanced.csv")

# 对5-mer序列进行独热编码（生成nt_pos1至nt_pos5）
data_df <- encode_5mer(data_df, seq_col = "5mer_sequence")  # 假设encode_5mer是包中的编码函数

# 将目标变量转为因子（负类为参考水平）
data_df$m6A_status <- factor(data_df$m6A_status, levels = c("Negative", "Positive"))

# 查看类别比例
table(data_df$m6A_status)
```

2.  Model training and prediction

```{r}
# 设置随机种子，划分训练集（80%）和测试集（20%）
set.seed(123)
train_idx <- sample(nrow(data_df), 0.8 * nrow(data_df))
train_df <- data_df[train_idx, ]
test_df <- data_df[-train_idx, ]

# 训练随机森林模型
rf_fit <- train_rf_model(
  data = train_df,
  target = "m6A_status",
  features = c("gc_content", "RNA_type", "RNA_region", "exon_length", 
               "distance_to_junction", "evolutionary_conservation", 
               "nt_pos1", "nt_pos2", "nt_pos4", "nt_pos5")
)

# 查看模型信息
rf_fit
```

3.  Model Prediction and Evaluation

```{r}
# 预测测试集概率
prob_test <- predict(rf_fit, test_df, type = "prob")

# 以0.5为阈值的混淆矩阵
cm <- confusion_matrix(test_df$m6A_status, prob_test[, "Positive"], cutoff = 0.5)
cm
# 计算评估指标（阈值0.5）
metrics_05 <- calculate_metrics(cm)  # 假设calculate_metrics是包中的函数
metrics_05
# 阈值0.1时的指标
metrics_01 <- calculate_metrics(cm, cutoff = 0.1)
metrics_01

```

This package offers robust model performance. The following is a visual display based on ROC and PRC curves, illustrating the model's strong predictive capabilities. #Visualization of model performance

```{r}
# 生成不同阈值下的性能指标
n <- 100
thresholds <- seq(0.05, 0.95, length.out = n)
plot_df <- calculate_threshold_metrics(test_df$m6A_status, prob_test[, "Positive"], thresholds)

# 绘制ROC曲线（红色）
ggplot(plot_df, aes(x = FPR, y = recall)) +
  geom_point(color = "red") +
  geom_line(color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(title = "ROC曲线", x = "假阳性率（FPR）", y = "召回率（Recall）") +
  theme_classic()

# 绘制PRC曲线（蓝色）
ggplot(plot_df, aes(x = recall, y = precision)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  labs(title = "PRC曲线", x = "召回率（Recall）", y = "精确率（Precision）") +
  theme_classic()
```

ROC and PRC curve description: The left side of the curve corresponds to a larger decision threshold (such as close to 1), representing a stricter positive judgment criterion (only high-probability samples are classified as m6A loci). The right side corresponds to a smaller threshold, and the judgment criteria are more lenient. Complete ROC/PRC curve and AUC value (based on precrec)

```{r}
# 使用precrec包生成完整曲线
library(precrec)
curves <- evalmod(scores = prob_test[, "Positive"], labels = test_df$m6A_status)

# 绘制曲线
autoplot(curves) + theme_classic()

# 提取AUC值
auc(curves)
```

#ROC curve

The ROC curve (Receiver Operating Characteristic curve) is used to evaluate the classification performance of the model. By comparing the True Positive Rate with the False Positive Rate, the performance of the model at different thresholds can be evaluated.

![ROC Curve](roc_curve.png)

#PRC curve

The PRC curve (Precision-Recall curve) pays more attention to evaluating the performance of the model in cases where the class distribution is imbalanced. It illustrates the relationship between the model's precision and recall rates.

![PRC Curve](prc_curve.png)
