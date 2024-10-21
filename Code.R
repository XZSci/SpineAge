#LOAD PACKAGE----
pacman::p_load(
  nhanesR, survey, ggplot2, forcats, dplyr, extrafont,
  recipes, gridExtra, rms, patchwork,
  grid, forestploter, tidyverse, forcats,caret,rsample,ranger
)
create_roc_df <- function(roc_object, model_name, time_point) {
  data.frame(
    FP = as.vector(t(roc_object$FP[, paste0("t=", time_point)])),
    TP = as.vector(t(roc_object$TP[, paste0("t=", time_point)])),
    Model = model_name,
    AUC = roc_object$AUC[paste0("t=", time_point)]
  )
}
getwd()
setwd("D:/desktop/database/nhanes/nhanesR/vertebrae4/202405/aging cell")
# extract the data----
demo <- nhs_tsv("demo.", years = 2013:2014)

T4 <- nhs_tsv("DXXT4.", years = 2013:2014)
T5 <- nhs_tsv("DXXT5.", years = 2013:2014)
T6 <- nhs_tsv("DXXT6.", years = 2013:2014)
T7 <- nhs_tsv("DXXT7.", years = 2013:2014)
T8 <- nhs_tsv("DXXT8.", years = 2013:2014)
T9 <- nhs_tsv("DXXT9.", years = 2013:2014)
T10 <- nhs_tsv("DXXT10.", years = 2013:2014)
T11 <- nhs_tsv("DXXT11.", years = 2013:2014)
T12 <- nhs_tsv("DXXT12.", years = 2013:2014)
L1 <- nhs_tsv("DXXL1.", years = 2013:2014)
L2 <- nhs_tsv("DXXL2.", years = 2013:2014)
L3 <- nhs_tsv("DXXL3.", years = 2013:2014)
L4 <- nhs_tsv("DXXL4.", years = 2013:2014)
cancer <- nhs_tsv("mcq_.", years = 2013:2014)
Xverte <- nhs_read(
  T4, "DXDLSPST:status", "DXXT4CC", "DXXPOSTH:T4P", "DXXANTEH:T4A", "DXXMIDH:T4M",
  T5, "DXXT5CC", "DXXPOSTH:T5P", "DXXANTEH:T5A", "DXXMIDH:T5M",
  T6, "DXXT6CC", "DXXPOSTH:T6P", "DXXANTEH:T6A", "DXXMIDH:T6M",
  T7, "DXXT7CC", "DXXPOSTH:T7P", "DXXANTEH:T7A", "DXXMIDH:T7M",
  T8, "DXXT8CC", "DXXPOSTH:T8P", "DXXANTEH:T8A", "DXXMIDH:T8M",
  T9, "DXXT9CC", "DXXPOSTH:T9P", "DXXANTEH:T9A", "DXXMIDH:T9M",
  T10, "DXXT10CC", "DXXPOSTH:T10P", "DXXANTEH:T10A", "DXXMIDH:T10M",
  T11, "DXXT11CC", "DXXPOSTH:T11P", "DXXANTEH:T11A", "DXXMIDH:T11M",
  T12, "DXXT12CC", "DXXPOSTH:T12P", "DXXANTEH:T12A", "DXXMIDH:T12M",
  L1, "DXXL1CC", "DXXPOSTH:L1P", "DXXANTEH:L1A", "DXXMIDH:L1M",
  L2, "DXXL2CC", "DXXPOSTH:L2P", "DXXANTEH:L2A", "DXXMIDH:L2M",
  L3, "DXXL3CC", "DXXPOSTH:L3P", "DXXANTEH:L3A", "DXXMIDH:L3M",
  L4, "DXXL4CC", "DXXPOSTH:L4P", "DXXANTEH:L4A", "DXXMIDH:L4M"
)
exm <- nhs_tsv("bmx.|bpx.", years = 2013:2014)
xtotal <- nhs_read(
  demo, "sdmvpsu", "sdmvstra", "RIAGENDR:sex", "RIDAGEYR:age",
  "RIDRETH1:eth1", "RIDRETH3:eth2",
  "RIDEXPRG:pregnacy", "WTINT2YR",
  "WTMEC2YR",
  exm, "bmxht:height", "bmxwaist", "bmxwt:weight",
  "bmxbmi:bmi"
)


Xcancer <- nhs_read(cancer, "mcq220:cancer")
smk <- diag_smoke(years = 2013:2014)
dm <- diag_DM(years = 2013:2014)
HBP <- diag_Hypertension(years = 2013:2014)
ASCVD <- diag_ASCVD(years = 2013:2014)
HyLipid <- diag_Hyperlipidemia(years = 2013:2014)
Mortality <- db_mort(years = 2013:2014)

df <- Full_Join(xtotal, Xverte, Xcancer, by = "seqn") |>
  diag_smoke(years = 2013:2014) |>
  diag_DM(years = 2013:2014) |>
  diag_Hypertension(years = 2013:2014) |>
  diag_ASCVD(years = 2013:2014) |>
  diag_Hyperlipidemia(years = 2013:2014) |>
  dex_eGFR() |>
  db_mort(years = 2013:2014)
# save original data----
write.table(df, "result.csv", row.names = FALSE, col.names = TRUE, sep = ",")

# datawashing----
df1 <- read.csv("result.csv")
## Recode(df$sex)
df1$sex <- Recode(df1$sex,
  "Male::1",
  "Female::0",
  to.numeric = FALSE
)
# Recode(df1$DM)
df1$DMnew <- Recode(df1$DM,
  "no::0",
  "DM::1",
  "IFG::0",
  "IGT::0",
  to.numeric = FALSE
)
# Recode(df1$eth1)
df1$eth1 <- Recode(df1$eth1,
  "Non-Hispanic Black::",
  "Other Race - Including Multi-Racial::Other race",
  "Non-Hispanic White::",
  "Mexican American::Other race",
  "Other Hispanic::Other race",
  to.numeric = FALSE
)
# Recode(df1$smoke)
df1$smoke <- Recode(df1$smoke,
  "never::0",
  "former::1",
  "now::1",
  "NA::",
  to.numeric = FALSE
)
# Recode(df1$caner)
df1$cancer <- Recode(df1$cancer,
  "No::0",
  "Yes::1",
  "NA::",
  to.numeric = FALSE
)
# Recode(df1$Hyperlipidemia)
df1$Hyperlipidemia <- Recode(df1$Hyperlipidemia,
  "yes::1",
  "no::0",
  to.numeric = FALSE
)
# Recode(df1$ASCVD)
df1$ASCVD <- Recode(df1$ASCVD,
  "no::0",
  "yes::1",
  "NA::",
  to.numeric = FALSE
)
# Recode(df1$Hypertension)
df1$Hypertension <- Recode(df1$Hypertension,
  "no::0",
  "yes::1",
  "NA::",
  to.numeric = FALSE
)
# Recode(df3$mortstat)
df1$mortstat <- Recode(df1$mortstat,
  "Assumed deceased::1",
  "Assumed alive::0",
  "NA::",
  to.numeric = FALSE
)



write.table(df1, "result1.csv", row.names = FALSE, col.names = TRUE, sep = ",")
## inclusion and exclusion----
df1 <- read.csv("result1.csv")
df2 <- drop_row(df1, df1$age < 40)
df2 <- drop_row(df2, is.na(df2$status))
df2 <- drop_row(df2, df2$status =='IVA Lateral Spine not scanned, pregnancy')
df2 <- drop_row(df2, df2$status =='IVA Lateral Spine not scanned, weight > 450 lbs')
df2 <- drop_row(df2, df2$status =='IVA Lateral Spine not scanned, other reason')
df2 <- drop_row(df2, df2$status =='IVA Lateral Spine scan completed, but one or more vertebrae are invalid')
df2 <- drop_row(df2, is.na(df2$bmi))
df2 <- drop_row(df2, is.na(df2$bmxwaist))
df2 <- drop_row(df2, is.na(df2$CKD_EPI_Scr_2009))
df2 <- drop_row(df2, is.na(df2$smoke))
df2 <- drop_row(df2, is.na(df2$Hypertension))
df2 <- drop_row(df2, is.na(df2$Hyperlipidemia))
df2 <- drop_row(df2, is.na(df2$DM))
df2 <- drop_row(df2, is.na(df2$ASCVD))
df2 <- drop_row(df2, df2$cancer == 1)
df2 <- drop_row(df2, df2$eligstat == "Ineligible")
df2$SumWt2y <- (df2$WTMEC2YR)
write.table(df2, "result2.csv", row.names = FALSE, col.names = TRUE, sep = ",")
df <- read.csv("result2.csv")

df$T4MPR <- df$T4M / df$T4P
df$T4APR <- df$T4A / df$T4P
df$T4PPR <- df$T4P / df$T5P


df$T5MPR <- df$T5M / df$T5P
df$T5APR <- df$T5A / df$T5P
df$T5PPR <- 2 * df$T5P / (df$T4P + df$T6P)

df$T6MPR <- df$T6M / df$T6P
df$T6APR <- df$T6A / df$T6P
df$T6PPR <- 2 * df$T6P / (df$T5P + df$T7P)

df$T7MPR <- df$T7M / df$T7P
df$T7APR <- df$T7A / df$T7P
df$T7PPR <- 2 * df$T7P / (df$T6P + df$T8P)


df$T8MPR <- df$T8M / df$T8P
df$T8APR <- df$T8A / df$T8P
df$T8PPR <- 2 * df$T8P / (df$T7P + df$T9P)


df$T9MPR <- df$T9M / df$T9P
df$T9APR <- df$T9A / df$T9P
df$T9PPR <- 2 * df$T9P / (df$T8P + df$T10P)


df$T10MPR <- df$T10M / df$T10P
df$T10APR <- df$T10A / df$T10P
df$T10PPR <- 2 * df$T10P / (df$T9P + df$T11P)

df$T11MPR <- df$T11M / df$T11P
df$T11APR <- df$T11A / df$T11P
df$T11PPR <- 2 * df$T11P / (df$T10P + df$T12P)


df$T12MPR <- df$T12M / df$T12P
df$T12APR <- df$T12A / df$T12P
df$T12PPR <- 2 * df$T12P / (df$T11P + df$L1P)

df$L1MPR <- df$L1M / df$L1P
df$L1APR <- df$L1A / df$L1P
df$L1PPR <- 2 * df$L1P / (df$T12P + df$L2P)

df$L2MPR <- df$L2M / df$L2P
df$L2APR <- df$L2A / df$L2P
df$L2PPR <- 2 * df$L2P / (df$L1P + df$L3P)

df$L3MPR <- df$L3M / df$L3P
df$L3APR <- df$L3A / df$L3P
df$L3PPR <- 2 * df$L3P / (df$L2P + df$L4P)

df$L4MPR <- df$L4M / df$L4P
df$L4APR <- df$L4A / df$L4P
df$L4PPR <- df$L4P / df$L3P


df <- select_col(df, c(
  "seqn", "sdmvpsu", "sdmvstra", "sex", "age",
  "eth1", "bmxwaist", "CKD_EPI_Scr_2009",
  "bmi", "smoke", "DMnew", "Hypertension", "ASCVD", "Hyperlipidemia",
  "T4MPR", "T4APR", "T4PPR", "T5MPR", "T5APR", "T5PPR", "T6MPR", "T6APR",
  "T6PPR", "T7MPR", "T7APR", "T7PPR", "T8MPR", "T8APR", "T8PPR", "T9MPR",
  "T9APR", "T9PPR", "T10MPR", "T10APR", "T10PPR", "T11MPR", "T11APR",
  "T11PPR", "T12MPR", "T12APR", "T12PPR", "L1MPR", "L1APR", "L1PPR", "L2MPR",
  "L2APR", "L2PPR", "L3MPR", "L3APR", "L3PPR", "L4MPR", "L4APR",
  "L4PPR", "seqn", "SumWt2y", "morstat", "eligstat", "mortstat", "ucod_leading", "diabetes", "hyperten", "permth_int", "permth_exm"
))
## save final data----
write.table(df, "result3.csv", row.names = FALSE, col.names = TRUE, sep = ",")


# analysis----
dff <- read.csv("result3.csv")

dff$sex <- as.factor(dff$sex)
dff$eth1 <- as.factor(dff$eth1)
dff$smoke <- as.factor(dff$smoke)
dff$DMnew <- as.factor(dff$DMnew)
dff$Hypertension <- as.factor(dff$Hypertension)
dff$ASCVD <- as.factor(dff$ASCVD)
dff$Hyperlipidemia <- as.factor(dff$Hyperlipidemia)
# table 1----
nhs <- svy_design(dff, weights = "SumWt2y")
svy_tableone(
  design = nhs, gv = c("sex", "DMnew", "eth1", "smoke", "Hypertension", "ASCVD", "Hyperlipidemia"),
  cv.nn = c("age", "bmi", "bmxwaist", "CKD_EPI_Scr_2009"),
  total = T, round = 2, by = "mortstat", g_perSQse = T
)

t_test_result <- svyttest(age ~ mortstat, design =nhs)
svyranktest(CKD_EPI_Scr_2009 ~ mortstat, design =nhs, test = "wilcoxon")$p.value
svyranktest(CKD_EPI_Scr_2009 ~ mortstat, design =nhs, test = "wilcoxon")$statistic

svychisq(~ASCVD+mortstat, design =nhs)



# LASSO----
## LASSO1P----
df1 <- dff[, c(15:53)]
y <- as.matrix(dff[, "age"])
x <- as.matrix(df1)
x <- scale(x)
summary(x)

set.seed(1)
lassoCV <- glmnet::cv.glmnet(x, y, family = "gaussian", maxit = 100000, nfolds = 10) # 拟合模型
lassoCV$lambda.min
lassoCV$lambda.1se
coef1 <- coef(lassoCV, s = "lambda.min")
names(which(coef(lassoCV, s = "lambda.min")[-1, ] != 0))
coef2 <- coef(lassoCV, s = "lambda.1se")
names(which(coef(lassoCV, s = "lambda.1se")[-1, ] != 0))
coef <- cbind2(coef1, coef2) # 合并显示
df_coef <- as.data.frame(as.matrix(coef))



plot.cvlasso <- data.frame(
  lambda = lassoCV$lambda,
  mse = lassoCV$cvm,
  sd = lassoCV$cvsd
)
(lasso1p <- ggplot(plot.cvlasso, aes(log(lambda), mse)) +
  geom_errorbar(aes(ymin = mse - sd, ymax = mse + sd), width = 0.1, color = "#bababa") + # 添加误差棒
  geom_point(color = "#404040", size = 1) + # 点图
  geom_vline(xintercept = log(lassoCV$lambda.min), linetype = "dashed", color = "#0571b0", linewidth = 1) +
  geom_vline(xintercept = log(lassoCV$lambda.1se), linetype = "dashed", color = "#de2d26", linewidth = 1) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "sans", size = 12, face = "bold"),
    axis.title = element_text(family = "sans", size = 12, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold")
  ) +
  xlim(-7, 1) +
  labs(
    title = "Lasso Regression with 10-fold Cross-validation",
    x = "Log-transformed Lambda", y = "Mean Squared Error"
  )


)

## LASSO2P----
fit.lasso <- glmnet::glmnet(x, y, family = "gaussian", alpha = 1, nlambda = 1000)
plot.lasso <- data.frame(as.matrix(fit.lasso$lambda), as.matrix(t(fit.lasso$beta)))
plot.lasso |>
  tidyr::pivot_longer(
    cols = 2:ncol(plot.lasso),
    names_to = "variable",
    values_to = "coef"
  ) -> plot.lasso

colors <- viridis::viridis(n = 39, option = "C")
(lasso2p <- ggplot(plot.lasso, aes(x = log(as.matrix.fit.lasso.lambda.), y = coef, color = variable)) +
  geom_line(linewidth = 0.6, alpha = 1) +
  scale_color_manual(values = colors) +
  labs(
    title = "LASSO Regression Path of Each Parameter",
    x = "Log-transformed Lambda",
    y = "Coefficients of Parameters"
  ) +
  xlim(-7, 1) +
  theme_classic() +
  theme(
    legend.position = "",
    plot.title = element_text(hjust = 0.5, family = "sans", size = 12, face = "bold"),
    axis.title = element_text(family = "sans", size = 12, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold")
  ) +
  geom_vline(xintercept = log(lassoCV$lambda.min), linetype = "dashed", color = "#0571b0", linewidth = 1) +
  geom_vline(xintercept = log(lassoCV$lambda.1se), linetype = "dashed", color = "#ca0020", linewidth = 1)
)


## LASSOTOTAL----
(lassototal <- lasso1p + lasso2p + plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(
    axis.title = element_text(family = "sans", size = 12, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold"),
    legend.title = element_text(family = "sans", size = 12, face = "bold"),
    legend.text = element_text(family = "sans", size = 12, face = "bold")
  )
)

ggsave("Figure 3.png", plot = lassototal, dpi = 600, width = 180, height = 110, units = "mm")
ggsave("Figure 3.pdf", plot = lassototal, dpi = 600, width = 180, height = 110, units = "mm")
ggsave("Figure 3.tiff", plot = lassototal, dpi = 600, width = 180, height = 110, units = "mm")


# Random Forest----
dfrf <- dff[, c("age", "T5APR", "T6APR", "T7APR", "T8APR", "T9MPR", "T9APR", "T10APR", "T11MPR", "L2APR", "L4MPR")]
baked_rf <- recipe(age ~ ., dfrf) |>
  prep(, training = dfrf) |>
  bake(, new_data = dfrf)
## RF model----
# 创建 10 折交叉验证，重复 5 次
set.seed(1)
hyper_grid <- expand.grid(
  ntrees=c(250,500,750,1000,1250,1500,1750,2000),
  mtry = c(1,2,3,4,5,6,7,8,9,10),
  min.node.size = c(1,2,3,4,5),
  replace = c(TRUE),
  sample.fraction = c(1),
  rmse = NA
)
cv_splits <- vfold_cv(baked_rf, v = 10, repeats = 5)

for (j in seq_len(nrow(hyper_grid))) {
  cv_rmse <- c()  
  
  for (i in seq_along(cv_splits$splits)) {
    train_data <- analysis(cv_splits$splits[[i]])
    test_data <- assessment(cv_splits$splits[[i]])
    fit <- ranger::ranger(
      formula = age ~ .,
      data = train_data,
      num.trees = hyper_grid$ntrees[j],
      mtry = hyper_grid$mtry[j],
      min.node.size = hyper_grid$min.node.size[j],
      replace = hyper_grid$replace[j],
      sample.fraction = hyper_grid$sample.fraction[j],
      seed = 1,
      respect.unordered.factors = 'order'
    )
    preds <- predict(fit, test_data)$predictions
    rmse <- sqrt(mean((test_data$age - preds)^2))
    cv_rmse <- c(cv_rmse, rmse)
  }
  hyper_grid$rmse[j] <- mean(cv_rmse)
}
print(hyper_grid)

rffit <- ranger::ranger(
  formula = age ~ .,
  data = baked_rf,
  num.trees = 500,
  mtry = 3,
  min.node.size = 5,
  replace = T,
  sample.fraction = 1,
  verbose = F,
  seed = 1,
  respect.unordered.factors = "order",
  importance = "impurity",
  oob.error = T
)

(MSE <- sqrt(rffit$prediction.error))
dff$pred_age <- predict(rf_final, newdata=dfrf) 
dff$pred_age <- predict(rffit, data = baked_rf)$predictions
dff$deltaage <- dff$pred_age - dff$age
dff$age2 <- dff$age
summary(dff$pred_age)
(MAE <- (sum(abs(dff$pred_age - dff$age))) / nrow(dff))
lm(age ~ pred_age, data = dff) |> summary()

# Linear association between SpineAge and Chronological Age----
(rfp1 <- ggplot(dff, aes(x = age, y = pred_age, color = sex)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("#de2d26", "#0571b0"), labels = c("Female", "Male")) + # 绘制散点图
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Chronological Age", y = "SpineAge", color = "Sex") +
  xlim(35, 85) +
  ylim(35, 85) +
  theme_classic() +
  theme(
    axis.title = element_text(family = "sans", size = 12, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold"),
    legend.title = element_text(family = "sans", size = 10, face = "bold"),
    legend.text = element_text(family = "sans", size = 10, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "right",
    legend.position.inside = c(0.1, 0.8)
  ) +
  annotate("text",
    x = 60, y = 80, size = 4, family = "sans", fontface = "bold",
    label = expression(paste(R^2, " = 0.956, MAE = 3.91"))
  )
)


## RF importance----
rfimp <- as.data.frame(rffit$variable.importance)
rfimp <- cbind(Parameters = rownames(rfimp), rfimp) |> dplyr::rename("InclineNodePurity" = "rffit$variable.importance")
# Recode(rfimp$Parameters)
rfimp$Parameters <- Recode(rfimp$Parameters,
  "T5APR::T5-APR",
  "T6APR::T6-APR",
  "T7APR::T7-APR",
  "T8APR::T8-APR",
  "T9MPR::T9-MPR",
  "T9APR::T9-APR",
  "T10APR::T10-APR",
  "T11MPR::T11-MPR",
  "L2APR::L2-APR",
  "L4MPR::L4-MPR",
  to.numeric = FALSE
)

(rfp2 <- rfimp |> mutate(Parameters = fct_reorder(Parameters, InclineNodePurity)) %>%
  ggplot(aes(x = Parameters, y = InclineNodePurity)) +
  geom_segment(aes(x = Parameters, xend = Parameters, y = 25000, yend = InclineNodePurity), linewidth = 1, color = "grey", linetype = "longdash") +
  geom_point(size = 3, color = "black", fill = alpha("#de2d26", 1), alpha = 1, shape = 21, stroke = 1) +
  theme_classic() +
  labs(x = "Parameters", y = "Importance (Increase of Node Purity)") +
  scale_y_continuous(limits = c(25000, NA)) +
  coord_flip() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "serif", size = 12, face = "bold"),
    axis.title = element_text(family = "serif", size = 12, face = "bold"),
    axis.text = element_text(family = "serif", size = 10, face = "bold")
  )
)
ggsave("RFimportance.png", plot = p3, dpi = 600, width = 10, height = 6.18, units = "in")
ggsave("RFimportance.tiff", plot = p3, dpi = 600, width = 10, height = 6.18, units = "in")
ggsave("RFimportance.pdf", plot = p3, dpi = 600, width = 10, height = 6.18, units = "in")
## SHAP----
explain_rf <- DALEX::explain(
  model = rffit,
  data = dplyr::select(as.data.frame(baked_rf), -"age"),
  y = baked_rf$age,
  label = "randomforest", type = "regression"
)
shap_rf <- DALEX::predict_parts(
  explainer = explain_rf,
  new_observation = baked_rf,
  type = "shap", B = 1000
)

shap1 <- shap_rf |> as.data.frame()
shap1$variable <- Recode(shap1$variable,
  "L2APR = 0.8902::L2-APR",
  "L4MPR = 0.945::L4-MPR",
  "T10APR = 0.9275::T10-APR",
  "T11MPR = 0.8722::T11-MPR",
  "T5APR = 0.9297::T5-APR",
  "T6APR = 0.8119::T6-APR",
  "T7APR = 0.8713::T7-APR",
  "T8APR = 0.872::T8-APR",
  "T9APR = 0.8659::T9-APR",
  "T9MPR = 0.8608::T9-MPR",
  to.numeric = FALSE
)

(rfp3 <- shap1 |>
  mutate(mean_con = mean(contribution), .by = variable) %>%
  mutate(variable = fct_reorder(variable, abs(mean_con))) %>%
  ggplot() +
  geom_bar(
    data = \(x) dplyr::distinct(x, variable, mean_con),
    aes(mean_con, variable, fill = mean_con > 0), alpha = 0.5,
    stat = "identity"
  ) +
  geom_boxplot(aes(contribution, variable, fill = mean_con > 0), width = 0.4) +
  scale_fill_manual(values = c("#0571b0", "#de2d26")) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, family = "serif", size = 12, face = "bold"),
    axis.title = element_text(family = "serif", size = 12, face = "bold"),
    axis.text = element_text(family = "serif", size = 10, face = "bold")
  ) +
  labs(
    x = "SHAP Value (Impact on SpineAge Calculation)", y = "Parameters"
  ))



ggsave("SHAP.png", plot = pshap, dpi = 600, width = 10, height = 6.18, units = "in")
ggsave("SHAP.tiff", plot = pshap, dpi = 600, width = 10, height = 6.18, units = "in")
ggsave("SHAP.pdf", plot = pshap, dpi = 600, width = 10, height = 6.18, units = "in")

## breakdown fig----
breakdown_rf <- DALEX::predict_parts(
  explainer = explain_rf,
  new_observation = baked_rf[5, ],
  type = "break_down",
  B = 1000 # 选择多少个排列组合
)
breakdown_rf$variable <- Recode(breakdown_rf$variable,
  "intercept::",
  "T8APR = 0.97741935483871::T8-APR",
  "T6APR = 0.862838915470494::T6-APR",
  "L4MPR = 0.911492734478203::L4-MPR",
  "T11MPR = 0.914285714285714::T11-MPR",
  "L2APR = 0.908972691807542::L2-APR",
  "T5APR = 0.934819897084048::T5-APR",
  "T9MPR = 0.928343949044586::T9-MPR",
  "T7APR = 0.888188976377953::T7-APR",
  "T9APR = 0.979299363057325::T9-APR",
  "T10APR = 0.919354838709677::T10-APR",
  "prediction::",
  to.numeric = FALSE
)
breakdown_rf$label <- Recode(breakdown_rf$label,
  "randomforest::SpineAge Calculation",
  to.numeric = FALSE
)


(rfp4 <- plot(breakdown_rf) +
  scale_fill_manual(values = c("#0571b0", "#de2d26", "#fdae61")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(family = "serif", size = 12, face = "bold"),
    axis.text = element_text(family = "serif", size = 10, face = "bold"),
    strip.text.x = element_text(family = "serif", size = 12, face = "bold"),
    
  ) +
  labs(
    title = "",
    x = "SpineAge", y = "Parameters"
  ) 
)
ggsave("Figure 2c.png", plot = rfp3, dpi = 600, width = 10, height = 6.18, units = "in")
ggsave("Figure 2c.tiff", plot = rfp3, dpi = 600, width = 10, height = 6.18, units = "in")
ggsave("Figure 2c.pdf", plot = rfp3, dpi = 600, width = 10, height = 6.18, units = "in")
## RFtotal----
(rftotal <- rfp2 + rfp3 + rfp4 + rfp1 + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(
    axis.title = element_text(family = "sans", size = 12, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold"),
    legend.title = element_text(family = "sans", size = 12, face = "bold"),
    legend.text = element_text(family = "sans", size = 10, face = "bold")
  )
)

ggsave("Figure 4.png", plot = rftotal, dpi = 600, width = 330, height = 150, units = "mm")
ggsave("Figure 4.tiff", plot = rftotal, dpi = 600, width = 330, height = 150, units = "mm")
ggsave("Figure 4.pdf", plot = rftotal, dpi = 600, width = 330, height = 240, units = "mm")

# timeroc----
# timeROC:
roc.predage <- timeROC::timeROC(
  T = dff$permth_int,
  delta = dff$mortstat,
  marker = dff$pred_age,
  cause = 1,
  times = c(12, 24, 36, 48, 60, 72),
  ROC = TRUE,
  iid = TRUE,
  weighting = "marginal"
)
roc.predage$AUC

roc.deltaage <- timeROC::timeROC(
    T = dff$permth_int,
    delta = dff$mortstat,
    marker = -dff$deltaage,
    cause = 1,
    times = c(12, 24, 36, 48, 60, 72),
    ROC = TRUE,
    iid = TRUE,
    weighting = "marginal"
)
roc.deltaage$AUC
roc.age <- timeROC::timeROC(
  T = dff$permth_int,
  delta = dff$mortstat,
  marker = dff$age,
  cause = 1,
  times = c(12, 24, 36, 48, 60, 72),
  ROC = TRUE,
  iid = TRUE,
  weighting = "marginal"
)
roc.age$AUC
(auc <- timeROC::compare(roc.age, roc.predage, adjusted = F))
##ROCSPINEAGE----
rocpredage <-data.frame(
    TP12 = roc.predage$TP[,1],
    TP24= roc.predage$TP[,2],
    TP36= roc.predage$TP[,3],
    TP48= roc.predage$TP[,4],
    TP60= roc.predage$TP[,5],
    FP12= roc.predage$FP[,1],
    FP24= roc.predage$FP[,2],
    FP36= roc.predage$FP[,3],
    FP48= roc.predage$FP[,4],
    FP60= roc.predage$FP[,5])
rocpredage$yd12 <- rocpredage$TP12-rocpredage$FP12
rocpredage$yd24 <- rocpredage$TP24-rocpredage$FP24
rocpredage$yd36 <- rocpredage$TP36-rocpredage$FP36
rocpredage$yd48 <- rocpredage$TP48-rocpredage$FP48
rocpredage$yd60 <- rocpredage$TP60-rocpredage$FP60


sort(unique(dff$pred_age),decreasing = T)[which.max(rocpredage$TP12-rocpredage$FP12)] 
sort(unique(dff$pred_age),decreasing = T)[which.max(rocpredage$TP24-rocpredage$FP24)]
sort(unique(dff$pred_age),decreasing = T)[which.max(rocpredage$TP36-rocpredage$FP36)]
sort(unique(dff$pred_age),decreasing = T)[which.max(rocpredage$TP48-rocpredage$FP48)]
sort(unique(dff$pred_age),decreasing = T)[which.max(rocpredage$TP60-rocpredage$FP60)]

sens12 <- rocpredage$TP12[which.max(rocpredage$TP12-rocpredage$FP12)]
speci12 <- (1-rocpredage$FP12[which.max(rocpredage$TP12-rocpredage$FP12)]) 
inci12 <-roc.predage$CumulativeIncidence['t=12']
#inci12 <- length(which(dff$permth_int<12))/length(dff$permth_int)
PPV12 <- (sens12*inci12)/(sens12*inci12+(1-inci12)*(1-speci12)) 
NPV12 <- (speci12*(1-inci12))/(speci12*(1-inci12)+(1-sens12)*inci12) 
print(c(sens12,speci12,inci12,PPV12,NPV12))

sens24 <- rocpredage$TP24[which.max(rocpredage$TP24-rocpredage$FP24)] 
speci24 <- (1-rocpredage$FP24[which.max(rocpredage$TP24-rocpredage$FP24)]) 
inci24 <-roc.predage$CumulativeIncidence['t=24']
#inci24 <- length(which(dff$permth_int<24))/length(dff$permth_int)
PPV24 <- (sens24*inci24)/(sens24*inci24+(1-inci24)*(1-speci24)) 
NPV24 <- (speci24*(1-inci24))/(speci24*(1-inci24)+(1-sens24)*inci24) 
print(c(sens24,speci24,inci24,PPV24,NPV24))


sens36 <- rocpredage$TP36[which.max(rocpredage$TP36-rocpredage$FP36)] 
speci36 <- (1-rocpredage$FP36[which.max(rocpredage$TP36-rocpredage$FP36)]) 
inci36 <-roc.predage$CumulativeIncidence['t=36']
#inci36 <- length(which(dff$permth_int<36))/length(dff$permth_int)
PPV36 <- (sens36*inci36)/(sens36*inci36+(1-inci36)*(1-speci36)) 
NPV36 <- (speci36*(1-inci36))/(speci36*(1-inci36)+(1-sens36)*inci36) 
print(c(sens36,speci36,inci36,PPV36,NPV36))

sens48 <- rocpredage$TP48[which.max(rocpredage$TP48-rocpredage$FP48)] 
speci48 <- (1-rocpredage$FP48[which.max(rocpredage$TP48-rocpredage$FP48)]) 
inci48 <-roc.predage$CumulativeIncidence['t=48']
#inci48 <- length(which(dff$permth_int<48))/length(dff$permth_int)
PPV48 <- (sens48*inci48)/(sens48*inci48+(1-inci48)*(1-speci48)) 
NPV48 <- (speci48*(1-inci48))/(speci48*(1-inci48)+(1-sens48)*inci48) 
print(c(sens48,speci48,inci48,PPV48,NPV48))

sens60 <- rocpredage$TP60[which.max(rocpredage$TP60-rocpredage$FP60)] 
speci60 <- (1-rocpredage$FP60[which.max(rocpredage$TP60-rocpredage$FP60)]) 
inci60 <-roc.predage$CumulativeIncidence['t=60']
#inci60 <- length(which(dff$permth_int<60))/length(dff$permth_int)
PPV60 <- (sens60*inci60)/(sens60*inci60+(1-inci60)*(1-speci60)) 
NPV60 <- (speci60*(1-inci60))/(speci60*(1-inci60)+(1-sens60)*inci60) 
print(c(sens60,speci60,inci60,PPV60,NPV60))


##ROCAGE----
rocage <-data.frame(
    TP12 = roc.age$TP[,1],
    TP24= roc.age$TP[,2],
    TP36= roc.age$TP[,3],
    TP48= roc.age$TP[,4],
    TP60= roc.age$TP[,5],
    FP12= roc.age$FP[,1],
    FP24= roc.age$FP[,2],
    FP36= roc.age$FP[,3],
    FP48= roc.age$FP[,4],
    FP60= roc.age$FP[,5])
rocage$yd12 <- rocage$TP12-rocage$FP12
rocage$yd24 <- rocage$TP24-rocage$FP24
rocage$yd36 <- rocage$TP36-rocage$FP36
rocage$yd48 <- rocage$TP48-rocage$FP48
rocage$yd60 <- rocage$TP60-rocage$FP60
sort(unique(dff$age),decreasing = T)[which.max(rocage$TP12-rocage$FP12)] 
sort(unique(dff$age),decreasing = T)[which.max(rocage$TP24-rocage$FP24)]
sort(unique(dff$age),decreasing = T)[which.max(rocage$TP36-rocage$FP36)]
sort(unique(dff$age),decreasing = T)[which.max(rocage$TP48-rocage$FP48)]
sort(unique(dff$age),decreasing = T)[which.max(rocage$TP60-rocage$FP60)]


sens12 <- rocage$TP12[which.max(rocage$TP12-rocage$FP12)] 
speci12 <- (1-rocage$FP12[which.max(rocage$TP12-rocage$FP12)]) 
inci12 <-roc.predage$CumulativeIncidence['t=12']
#inci12 <- length(which(dff$permth_int<12))/length(dff$permth_int)
PPV12 <- (sens12*inci12)/(sens12*inci12+(1-inci12)*(1-speci12)) 
NPV12 <- (speci12*(1-inci12))/(speci12*(1-inci12)+(1-sens12)*inci12) 
print(c(sens12,speci12,inci12,PPV12,NPV12))

sens24 <- rocage$TP24[which.max(rocage$TP24-rocage$FP24)] 
speci24 <- (1-rocage$FP24[which.max(rocage$TP24-rocage$FP24)]) 
inci24 <-roc.predage$CumulativeIncidence['t=24']
#inci24 <- length(which(dff$permth_int<24))/length(dff$permth_int)
PPV24 <- (sens24*inci24)/(sens24*inci24+(1-inci24)*(1-speci24)) 
NPV24 <- (speci24*(1-inci24))/(speci24*(1-inci24)+(1-sens24)*inci24) 
print(c(sens24,speci24,inci24,PPV24,NPV24))


sens36 <- rocage$TP36[which.max(rocage$TP36-rocage$FP36)] 
speci36 <- (1-rocage$FP36[which.max(rocage$TP36-rocage$FP36)]) 
inci36 <-roc.predage$CumulativeIncidence['t=36']
#inci36 <- length(which(dff$permth_int<36))/length(dff$permth_int)
PPV36 <- (sens36*inci36)/(sens36*inci36+(1-inci36)*(1-speci36)) 
NPV36 <- (speci36*(1-inci36))/(speci36*(1-inci36)+(1-sens36)*inci36) 
print(c(sens36,speci36,inci36,PPV36,NPV36))

sens48 <- rocage$TP48[which.max(rocage$TP48-rocage$FP48)] 
speci48 <- (1-rocage$FP48[which.max(rocage$TP48-rocage$FP48)]) 
inci48 <-roc.predage$CumulativeIncidence['t=48']
#inci48 <- length(which(dff$permth_int<48))/length(dff$permth_int)
PPV48 <- (sens48*inci48)/(sens48*inci48+(1-inci48)*(1-speci48)) 
NPV48 <- (speci48*(1-inci48))/(speci48*(1-inci48)+(1-sens48)*inci48) 
print(c(sens48,speci48,inci48,PPV48,NPV48))

sens60 <- rocage$TP60[which.max(rocage$TP60-rocage$FP60)] 
speci60 <- (1-rocage$FP60[which.max(rocage$TP60-rocage$FP60)]) 
inci60 <-roc.predage$CumulativeIncidence['t=60']
#inci60 <- length(which(dff$permth_int<60))/length(dff$permth_int)
PPV60 <- (sens60*inci60)/(sens60*inci60+(1-inci60)*(1-speci60)) 
NPV60 <- (speci60*(1-inci60))/(speci60*(1-inci60)+(1-sens60)*inci60) 
print(c(sens60,speci60,inci60,PPV60,NPV60))







dff[which(dff$permth_int<12),]
age80 <- dff[which(dff$age>56.60487),]

rocage12 <- create_roc_df(roc.age, "Chronological Age", 12)
rocpredage12 <- create_roc_df(roc.predage, "SpineAge", 12)
roc12 <- rbind(rocage12, rocpredage12)
roc12$Model <- factor(roc12$Model, levels = c("SpineAge", "Chronological Age"))


(p12 <- ggplot(data = roc12, aes(FP, TP, color = Model)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, color = "black", size = 1, linetype = 2) +
  scale_color_manual(
    values = c("SpineAge" = "#de2d26", "Chronological Age" = "#0571b0")
  ) +
  annotate("text",
    x = 0.8, y = 0.25, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 12-month: ",
      sprintf("%.3f", roc12$AUC[roc12$Model == "SpineAge"])
    ),
    color = "#de2d26"
  ) +
  annotate("text",
    x = 0.8, y = 0.15, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 12-month:",
      sprintf("%.3f", roc12$AUC[roc12$Model == "Chronological Age"])
    ),
    color = "#0571b0"
  ) +
  annotate("text",
    x = 0.8, y = 0.05, size = 2.5, family = "sans",
    label = paste0(
      "P-value:",
      sprintf("%.3f", auc$p_values_AUC["t=12"])
    ),
    color = "black"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_classic() +
  theme(
    axis.title = element_text(family = "sans", size = 10, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85)
  ))
rocage24 <- create_roc_df(roc.age, "Chronological Age", 24)
rocpredage24 <- create_roc_df(roc.predage, "SpineAge", 24)
roc24 <- rbind(rocage24, rocpredage24)
roc24$Model <- factor(roc24$Model, levels = c("SpineAge", "Chronological Age"))
(p24 <- ggplot(data = roc24, aes(FP, TP, color = Model)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, color = "black", size = 1, linetype = 2) +
  scale_color_manual(
    values = c("SpineAge" = "#de2d26", "Chronological Age" = "#0571b0")
  ) +
  annotate("text",
    x = 0.8, y = 0.25, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 24-month: ",
      sprintf("%.3f", roc24$AUC[roc24$Model == "SpineAge"])
    ),
    color = "#de2d26"
  ) +
  annotate("text",
    x = 0.8, y = 0.15, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 24-month:",
      sprintf("%.3f", roc24$AUC[roc24$Model == "Chronological Age"])
    ),
    color = "#0571b0"
  ) +
  annotate("text",
    x = 0.8, y = 0.05, size = 2.5, family = "sans",
    label = paste0(
      "P-value:",
      sprintf("%.3f", auc$p_values_AUC["t=24"])
    ),
    color = "black"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_classic() +
  theme(
    axis.title = element_text(family = "sans", size = 10, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85)
  ))
rocage36 <- create_roc_df(roc.age, "Chronological Age", 36)
rocpredage36 <- create_roc_df(roc.predage, "SpineAge", 36)
roc36 <- rbind(rocage36, rocpredage36)
roc36$Model <- factor(roc36$Model, levels = c("SpineAge", "Chronological Age"))


(p36 <- ggplot(data = roc36, aes(FP, TP, color = Model)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, color = "black", size = 1, linetype = 2) +
  scale_color_manual(
    values = c("SpineAge" = "#de2d26", "Chronological Age" = "#0571b0")
  ) +
  annotate("text",
    x = 0.8, y = 0.25, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 36-month: ",
      sprintf("%.3f", roc36$AUC[roc36$Model == "SpineAge"])
    ),
    color = "#de2d26"
  ) +
  annotate("text",
    x = 0.8, y = 0.15, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 36-month:",
      sprintf("%.3f", roc36$AUC[roc36$Model == "Chronological Age"])
    ),
    color = "#0571b0"
  ) +
  annotate("text",
    x = 0.8, y = 0.05, size = 2.5, family = "sans",
    label = paste0(
      "P-value:",
      sprintf("%.3f", auc$p_values_AUC["t=36"])
    ),
    color = "black"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_classic() +
  theme(
    axis.title = element_text(family = "sans", size = 10, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85)
  ))
rocage48 <- create_roc_df(roc.age, "Chronological Age", 48)
rocpredage48 <- create_roc_df(roc.predage, "SpineAge", 48)
roc48 <- rbind(rocage48, rocpredage48)
roc48$Model <- factor(roc48$Model, levels = c("SpineAge", "Chronological Age"))


(p48 <- ggplot(data = roc48, aes(FP, TP, color = Model)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, color = "black", size = 1, linetype = 2) +
  scale_color_manual(
    values = c("SpineAge" = "#de2d26", "Chronological Age" = "#0571b0")
  ) +
  annotate("text",
    x = 0.8, y = 0.25, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 48-month: ",
      sprintf("%.3f", roc48$AUC[roc48$Model == "SpineAge"])
    ),
    color = "#de2d26"
  ) +
  annotate("text",
    x = 0.8, y = 0.15, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 48-month:",
      sprintf("%.3f", roc48$AUC[roc48$Model == "Chronological Age"])
    ),
    color = "#0571b0"
  ) +
  annotate("text",
    x = 0.8, y = 0.05, size = 2.5, family = "sans",
    label = paste0(
      "P-value:",
      sprintf("%.3f", auc$p_values_AUC["t=48"])
    ),
    color = "black"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_classic() +
  theme(
    axis.title = element_text(family = "sans", size = 10, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85)
  ))
rocage60 <- create_roc_df(roc.age, "Chronological Age", 60)
rocpredage60 <- create_roc_df(roc.predage, "SpineAge", 60)
roc60 <- rbind(rocage60, rocpredage60)
roc60$Model <- factor(roc60$Model, levels = c("SpineAge", "Chronological Age"))


(p60 <- ggplot(data = roc60, aes(FP, TP, color = Model)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, color = "black", size = 1, linetype = 2) +
  scale_color_manual(
    values = c("SpineAge" = "#de2d26", "Chronological Age" = "#0571b0")
  ) +
  annotate("text",
    x = 0.8, y = 0.25, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 60-month: ",
      sprintf("%.3f", roc60$AUC[roc60$Model == "SpineAge"])
    ),
    color = "#de2d26"
  ) +
  annotate("text",
    x = 0.8, y = 0.15, size = 2.5, family = "sans",
    label = paste0(
      "AUC at 60-month:",
      sprintf("%.3f", roc60$AUC[roc60$Model == "Chronological Age"])
    ),
    color = "#0571b0"
  ) +
  annotate("text",
    x = 0.8, y = 0.05, size = 2.5, family = "sans",
    label = paste0(
      "P-value:",
      sprintf("%.3f", auc$p_values_AUC["t=60"])
    ),
    color = "black"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_classic() +
  theme(
    axis.title = element_text(family = "sans", size = 10, face = "bold"),
    axis.text = element_text(family = "sans", size = 10, face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85)
  ))

(roctotal <- p12 + p24 + p36 + p48 + p60 + guide_area() + plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  labs(x = "False positive rate", y = "True positive rate") &
  theme(
    axis.title = element_text(family = "sans", size = 10, face = "bold"),
    axis.text = element_text(family = "sans", size = 8, face = "bold"),
    legend.title = element_text(family = "sans", size = 10, face = "bold"),
    legend.text = element_text(family = "sans", size = 8, face = "bold")
  )
)


ggsave("Figure 5.png", plot = roctotal, dpi = 600, width = 180, height = 180, units = "mm")
ggsave("Figure 5.tiff", plot = roctotal, dpi = 600, width = 180, height = 180, units = "mm")
ggsave("Figure 5.pdf", plot = roctotal, dpi = 600, width = 180, height = 180, units = "mm")

ggsave("pic24.emf", plot = p24,device = devEMF::emf, width = 10, height = 6, units = "cm")
ggsave("pic60.emf", plot = p60,device = devEMF::emf, width = 10, height = 6, units = "cm")
# Cox regresion----
dff$deltaageQ <- quant(dff$deltaage, n = 4, Q = TRUE, round = 3)
dff$deltaageQ.median <- quant.median(dff$deltaage, n = 4, round = 3)

nhs <- svy_design(dff, weights = "SumWt2y")

svycoxph(Surv(permth_int, mortstat) ~ sex + age + eth1 + CKD_EPI_Scr_2009 + bmxwaist + bmi + DMnew + smoke + Hypertension +
  ASCVD + Hyperlipidemia + deltaage, nhs) %>% reg_table(round = 3)
svycoxph(Surv(permth_int, mortstat) ~ sex + age + eth1 + CKD_EPI_Scr_2009 + bmxwaist + bmi + DMnew + smoke + Hypertension +
             ASCVD + Hyperlipidemia + deltaage, nhs) |> cox.zph()

svycoxph(Surv(permth_int, mortstat) ~ sex + age + eth1 + CKD_EPI_Scr_2009 + bmxwaist + bmi + DMnew + smoke + Hypertension +
  ASCVD + Hyperlipidemia + deltaageQ, nhs) %>% reg_table(round = 3)
svycoxph(Surv(permth_int, mortstat) ~ sex + age + eth1 + CKD_EPI_Scr_2009 + bmxwaist + bmi + DMnew + smoke + Hypertension +
             ASCVD + Hyperlipidemia + deltaageQ, nhs) |> cox.zph()
svycoxph(Surv(permth_int, mortstat) ~ sex + age + eth1 + CKD_EPI_Scr_2009 + bmxwaist + bmi + DMnew + smoke + Hypertension +
  ASCVD + Hyperlipidemia + deltaageQ.median, nhs) %>% reg_table(round = 3)
svycoxph(Surv(permth_int, mortstat) ~ sex + age + eth1 + CKD_EPI_Scr_2009 + bmxwaist + bmi + DMnew + smoke + Hypertension +
             ASCVD + Hyperlipidemia + deltaageQ.median, nhs) |> cox.zph()
# Forestplot----


fore <- read.csv("forestplot.csv", fileEncoding = "GBK")
fore$Variable <- ifelse(fore$coef == "",
  fore$Variable,
  paste0("    ", fore$Variable)
)
fore$` ` <- paste(rep("   ", 10), collapse = " ")
fore$`Hazards Ratio (95% CI)` <- ifelse(is.na(fore$Hazards.Ratio), "",
  sprintf(
    "%.3f (%.3f to %.3f)",
    fore$Hazards.Ratio, fore$lower, fore$upper
  )
)
fore <- fore |> dplyr::rename("Coefficient" = "coef")
fore <- fore |> dplyr::rename("   SE" = "se.coef.")
fore <- fore |> dplyr::rename("Robust SE" = "robust.se")
fore <- fore |> dplyr::rename("    Z" = "z")
fore <- fore |> dplyr::rename("P-value" = "P.value")
fore <- fore |> dplyr::rename("P-value for Trend" = "P.for.trend")
my_theme <- forest_theme(
  core = list(bg_params = list(fill = c("white")), padding = unit(c(4, 3), "mm")),
  base_size = 10,
  base_family = "sans",
  refline_col = "black", # 参考线颜色
  refline_lty = 2, # 参考线类型
  panel_bg = "white", # 面板背景颜色
  ci_col = "#de2d26", # 置信区间颜色
  ci_fill = "#0571b0", # 置信区间填充颜色
  ci_alpha = 1,
  ci_Theight = 0.3,
  ci_lwd = 2,
  ci_lty = 1
)
(pfore <- forest(fore[, c(1:5, 12, 11, 9, 10)],
  est = fore$Hazards.Ratio,
  lower = fore$lower,
  upper = fore$upper,
  ci_column = 7, # CI图所在的列
  ref_line = 1,
  arrow_lab = c("Lower Risk", "Higher Risk"),
  xlim = c(0.5, 64),
  ticks_at = c(0.5, 1, 2, 4, 8, 16, 32), x_trans = "log2",
  theme = my_theme
))
(pfore1 <- edit_plot(pfore,
  row = c(1, 3),
  gp = gpar(fontface = "bold")
) |>
  edit_plot(,
    row = c(2, 5, 6, 7), col = 8,
    gp = gpar(fontface = "bold")
  ) |>
  edit_plot(,
    col = 1:9,
    row = c(2, 4, 6),
    which = "background",
    gp = gpar(fill = "#e0f3f8")
  ) |>
  edit_plot(,
    col = 2:9,
    which = "text",
    hjust = unit(0.5, "npc"),
    x = unit(0.5, "npc")
  ) |>
  add_border(,
    part = "header",
    row = 1,
    where = "top",
    gp = gpar(lwd = 2)
  ) |>
  add_border(,
    part = "header",
    row = 8,
    where = "bottom",
    gp = gpar(lwd = 2)
  ) |>
  add_border(,
    part = "header",
    row = 2,
    where = "top",
    gp = gpar(lwd = 1)
  ) |>
  edit_plot(,
    row = c(1, 3),
    part = "body",
    hjust = unit(0.5, "npc"),
    x = unit(0.5, "npc")
  )
)
ggsave("Figure 6.png", plot = pfore1, dpi = 600, width = 250, height = 100, units = "mm")
ggsave("Figure 6.tiff", plot = pfore1, dpi = 600, width = 250, height = 100, units = "mm")
ggsave("Figure 6.pdf", plot = pfore1, dpi = 600, width = 250, height = 100, units = "mm")
