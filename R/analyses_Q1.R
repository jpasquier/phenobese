library(ggplot2)
library(parallel)
library(writexl)

kable <- knitr::kable

options(mc.cores = detectCores() - 1)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Hypogonadisme)")

# Data
load("data/bl.rda")

# Output directory
outdir <- "results/analyses_Q1_20210927"
if (!dir.exists(outdir)) dir.create(outdir)

# Recode `Group` as factor
bl$Group <- factor(bl$Group, c("Lean control", "Non-HH obese", "HH obese"))

# Q1 - Compare all parameters at baseline - Numeric variables
V <- names(bl)[sapply(bl, class) == "numeric"]
rows <- mclapply(setNames(V, V), function(v) {
  x <- list(
    All = bl[[v]],
    Control = bl[bl$Group %in% "Lean control", v],
    Obese = bl[grepl("obese", bl$Group), v],
    NonHH = bl[bl$Group %in% "Non-HH obese", v],
    HH = bl[bl$Group %in% "HH obese", v]
  )
  x <- lapply(x, function(z) z[!is.na(z)])
  n <- length
  q25 <- function(z) unname(quantile(z, 0.25))
  q75 <- function(z) unname(quantile(z, 0.75))
  fct <- c("n", "mean", "sd", "min", "q25", "median", "q75", "max")
  r <- do.call(c, lapply(names(x), function(u) {
    r <- sapply(fct, function(f) {
      if (length(x[[u]]) > 0) get(f)(x[[u]]) else NA
    })
    names(r) <- paste(u, names(r), sep = ".")
    return(r)
  }))
  lg <- do.call(rbind, lapply(names(x), function(u) {
    if (length(x[[u]]) > 0) data.frame(grp = u, y = x[[u]]) else NULL
  }))
  lg$grp <- droplevels(factor(
    lg$grp, c("All", "Control", "Obese", "NonHH", "HH")))
  fig <- ggplot(lg, aes(x = grp, y = y)) +
    geom_boxplot() +
    labs(title = v, x = "", y = v)
  r <- cbind(
    data.frame(Variable = v),
    t(r),
    Control_Obese_t.test_p.value =
      tryCatch(t.test(x$Control, x$Obese)$p.value, error = function(err) NA),
    Control_Obese_wilcox_p.value =
      tryCatch(wilcox.test(x$Control, x$Obese, exact = FALSE)$p.value,
               error = function(err) NA),
    HH_NonHH_t.test_p.value =
      tryCatch(t.test(x$HH, x$NonHH)$p.value, error = function(err) NA),
    HH_NonHH_wilcox_p.value =
      tryCatch(wilcox.test(x$HH, x$NonHH, exact = FALSE)$p.value,
               error = function(err) NA)
  )
  attr(r, "fig") <- fig
  return(r)
})
Q1_num <- do.call(rbind, rows)
Q1_boxplots <- lapply(rows, function(r) attr(r, "fig"))
rm(V, rows)

# Q1 - Compare all parameters at baseline - Categorical variables
V <- names(bl)[sapply(bl, class) == "character"]
V <- V[grepl("^(ADAM|Glucose|NAFLD)", V)]
Q1_cat <- do.call(rbind, mclapply(setNames(V, V), function(v) {
  lvls <- unique(na.omit(bl[[v]]))
  x <- list(
    All = bl[[v]],
    Control = bl[bl$Group %in% "Lean control", v],
    Obese = bl[grepl("obese", bl$Group), v],
    NonHH = bl[bl$Group %in% "Non-HH obese", v],
    HH = bl[bl$Group %in% "HH obese", v]
  )
  x <- lapply(x, function(z) factor(z[!is.na(z)], lvls))
  Merge <- function(x1, x2) merge(x1, x2, by = "value", all = TRUE)
  r <- Reduce(Merge, lapply(names(x), function(u) {
    r <- table(x[[u]])
    r <- cbind(n = r, prop = prop.table(r))
    colnames(r) <- paste(u, colnames(r), sep = ".")
    cbind(data.frame(value = rownames(r)), r)
  }))
  r <- r[match(lvls, r$value), ]
  tbl1 <- rbind(data.frame(g = 1, v = x$Control),
                data.frame(g = 2, v = x$Obese))
  tbl1 <- table(tbl1$g, tbl1$v)
  tbl2 <- rbind(data.frame(g = 1, v = x$HH),
                data.frame(g = 2, v = x$NonHH))
  tbl2 <- table(tbl2$g, tbl2$v)
  na <- rep(NA, nrow(r) - 1)
  cbind(
    variable = c(v, na),
    r,
    Control_Obese_fisher.test_p.value = c(fisher.test(tbl1)$p.value, na),
    HH_NonHH_fisher.test_p.value = c(fisher.test(tbl2)$p.value, na)
  )
}))
rm(V)

# Q1 - Area under the curve - Glucose/Insuline
Y <- c("Glucose", "Insuline")
auc <- lapply(setNames(Y, Y), function(y) {
  V <- grep(paste0("^T[0-9]+ ", y, "$"), names(bl), value = TRUE)
  sapply(1:nrow(bl), function(i) {
    x <- (0:4)*30
    y <- bl[i, V]
    if (is.na(y[1]) | is.na(y[length(y)])) {
      auc <- NA
    } else {
      x <- x[!is.na(y)]
      y <- y[!is.na(y)]
      n <- length(y)
      auc <- sum((x[-1] - x[-n]) * (y[-1] + y[-n]) / 2)
    }
    return(auc)
  })
})
bl$AUC_Glucose <- auc$Glucose
bl$AUC_Insuline <- auc$Insuline
Y1 <- paste0("AUC_", Y)
Q1_auc_tbls <- lapply(setNames(Y1, Y), function(y) {
  fml <- as.formula(paste(y, "~ Group"))
  Merge <- function(x1, x2) merge(x1 ,x2, by = "Group", all = TRUE, sort = FALSE)
  n <- function(z, na.rm = FALSE) sum(!is.na(z))
  q25 <- function(z, na.rm = FALSE) unname(quantile(z, 0.25, na.rm = na.rm))
  q75 <- function(z, na.rm = FALSE) unname(quantile(z, 0.75, na.rm = na.rm))
  fct <- c("n", "mean", "sd", "min", "q25", "median", "q75", "max")
  Reduce(Merge, lapply(fct, function(f) {
    setNames(aggregate(fml, bl, get(f), na.rm = TRUE), c("Group", f))
  }))
})
rm(auc, Y, Y1)

# Q1 - Area under the curve - Glucose/Insuline - Boxplots
Y <- c("Glucose", "Insuline")
Q1_auc_boxplots <- lapply(setNames(Y, Y), function(y) {
  y1 <- paste0("AUC_", y)
  d <- na.omit(bl[c("Group", y1)])
  ggplot(d, aes_string(x = "Group", y = y1)) +
  geom_boxplot() +
  labs(title = y, x = "")
})
rm(Y)

# Q1 - Export tables
write_xlsx(list(numeric = Q1_num, categorical = Q1_cat),
           file.path(outdir, "Q1_tables.xlsx"))

# Q1 - Export boxplots
f1 <- file.path(outdir, "Q1_boxplots.pdf")
pdf(f1)
for (fig in Q1_boxplots) print(fig)
dev.off()
f2 <- "/tmp/tmp_bookmarks.txt"
for(i in 1:length(Q1_boxplots)) {
  write("BookmarkBegin", file = f2, append = (i!=1))
  s <- paste("BookmarkTitle:", names(Q1_boxplots)[i])
  write(s, file = f2, append = TRUE)
  write("BookmarkLevel: 1", file = f2, append = TRUE)
  write(paste("BookmarkPageNumber:", i), file = f2, append = TRUE)
}
system(paste("pdftk", f1, "update_info", f2, "output /tmp/tmp_fig.pdf"))
system(paste("mv /tmp/tmp_fig.pdf", f1, "&& rm", f2))
rm(fig, f1, f2, i, s)

# Q1 - Area under the curve - Glucose/Insuline - Export tables and analyses
for (y in c("Glucose", "Insuline")) {
  z <- file(file.path(outdir, paste0("Q1_AUC_", y, ".txt")), open = "wt")
  sink(z)
  sink(z, type = "message")
  cat(paste0(y, "\n========\n"))
  print(kable(Q1_auc_tbls[[y]]))
  cat("\n\nLinear regression\n-----------------\n")
  fml <- as.formula(paste0("AUC_", y, " ~ Group"))
  fit <- do.call("lm", list(formula = fml, data = quote(bl)))
  print(summary(fit))
  cat("\nANOVA Table\n-----------\n")
  print(anova(fit))
  sink(type = "message")
  sink()
}
rm(y, z, fml, fit)
write_xlsx(Q1_auc_tbls, file.path(outdir, "Q1_AUC.xlsx"))

# Q1 - Area under the curve - Glucose/Insuline - Export boxplots
pdf(file.path(outdir, "Q1_AUC_boxplots.pdf"))
for (fig in Q1_auc_boxplots) print(fig)
dev.off()
rm(fig)
