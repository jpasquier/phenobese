library(ggplot2)
library(parallel)
library(pROC)
library(writexl)

options(mc.cores = detectCores() - 1)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Phenobese)")

# Data
load("data/bl.rda")

# Output directory
outdir <- paste0("results/analyses_Q2_", format(Sys.Date(), "%Y%m%d"))

# Selection of obese patients
ob <- bl[grepl("obese", bl$Group), ]

# Define the variable "Poor metabolic profile"
ncrit <- 2
ob$PoorMetabolicProfile <-
  ob$`Glucose tolerance OGTT` %in% c("DT2", "IGT") &
  ob$`HOMA-IR` > 4
crit <- paste("Poor Metabolic Profile IF Glucose tolerance OGTT in {DT2, IGT}",
              "& HOMA-IR > 4")
if (ncrit == 2) {
  outdir <- file.path(outdir, "Two_criteria")
} else {
  outdir <- file.path(outdir, "Three_criteria")
  ob$PoorMetabolicProfile <- ob$PoorMetabolicProfile & ob$`hs-CRP` > 5
  crit <- paste(crit, "& hs-CRP > 5")
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
cat(crit, file = file.path(outdir, "criteria.txt"))
rm(crit, ncrit)

# Predictors
ob$Testosterone <- ob$`Testosterone (nmol/l)`
ob$TestosteroneLibre <- ob$`Testosterone libre (pmol/l)`
ob$TE <- ob$`T/E`
m <- glm(PoorMetabolicProfile ~ Testosterone + TestosteroneLibre + TE,
         family = binomial, data = ob)
ob$LinearCombination <- predict(m, newdata = ob)

# ROC analyses
X <- c(`Testosterone (nmol/l)` = "Testosterone",
       `Testosterone libre (pmol/l)` = "TestosteroneLibre",
       `T/E` = "TE",
       `Linear combination` = "LinearCombination")
R <- lapply(setNames(names(X), X), function(u) {
  y <- "PoorMetabolicProfile"
  x <- X[u]
  ob <- na.omit(ob[c(x, y)])
  r <- roc(as.formula(paste(y, "~", x)), data = ob)
  roc_curve <- ggroc(r) +
    geom_abline(intercept = 1, slope = 1, color = "grey60") +
    annotate("text", x = .25, y = .25,
             label = paste("Area under the curve:", signif(r$auc, 4))) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()) +
    labs(title = u, subtitle = "ROC curve",
         caption = paste("Number of observations:", nrow(ob)))
  roc_table <- data.frame(
    Threshold = r$thresholds,
    Sensitivity = r$sensitivities,
    Specificity = r$specificities
  )
  list(curve = roc_curve, table = roc_table)
})

# Export tables and curves
tmp <- as.data.frame(with(ob, table(PoorMetabolicProfile, useNA = "ifany")),
                     stringsAsFactors = FALSE)
tmp[is.na(tmp)] <- "NA"
write_xlsx(tmp, file.path(outdir, "Q2_PoorMetabolicProfile.xlsx"))
sink(file.path(outdir, "Q2_linear_combination.txt"))
print(summary(m))
sink()
write_xlsx(lapply(R, function(r) r$table),
           file.path(outdir, "Q2_roc_tables.xlsx"))
pdf(file.path(outdir, "Q2_roc_curves.pdf"))
for(r in R) print(r$curve)
dev.off()
rm(tmp, r)

# SessionInfo
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
