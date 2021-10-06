library(ggplot2)
library(gridExtra)
library(parallel)
library(writexl)

options(mc.cores = detectCores() - 1)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Phenobese)")

# Output directory
outdir <- paste0("results/analyses_Q5_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# Data
for (z in list.files("data")) load(file.path("data", z))
rm(z)

# Independent variables
load("data/metabolic_variables.dta")
load("data/exploratory_variables.dta")
X <- c(metabolic_variables, exploratory_variables)

# lapply(list(bl = bl, j1 = j1, j28 = j28, m3 = m3, m6 = m6, m12 = m12), function(d) {
#   X[!(X %in% names(d))]
# })

# metabolic_variables <- gsub(" +", " ", metabolic_variables)
# save(metabolic_variables, file = "data/metabolic_variables.dta", compress = "xz")
# exploratory_variables <- gsub(" +", " ", exploratory_variables)
# save(exploratory_variables, file = "data/exploratory_variables.dta", compress = "xz")


dta <- mclapply(c("bl", "j1", "j28", "m3", "m6", "m12"), function(z) {
  d <- get(z)
  x <- names(d)[sapply(d, class) == "numeric"]
  x <- if (z == "bl") c("Subject ID", "Group", x) else c("Subject ID", x)
  d <- d[x]
  b <- !(names(d) %in% c("Subject ID", "Group"))
  names(d)[b] <- paste0(names(d)[b], " (", z, ")")
  return(d)
})
dta <- Reduce(function(x, y) merge(x, y, by = "Subject ID", all = TRUE), dta)

W <- names(dta)[!(names(dta) %in% c("Subject ID", "Group"))]
W <- W[!grepl("\\(bl\\)$", W)]
for (x in W) {
  x0 <- sub("\\((j(1|28)|m(3|6|12))\\)$", "(bl)", x)
  if (x0 %in% names(dta)) {
    y <- paste("Δ", x)
    dta[[y]] <- dta[[x]] - dta[[x0]]
  }
}
rm(W, x, x0, y)

# Univariate regressions - Table
uni_reg_tbl <- do.call(rbind, mclapply(X, function(x) {
  Y <- paste("Δ Testosterone", c("(nmol/l)", "libre (pmol/l)"), "(m12)")
  Z <- c("j1", "j28", "m3", "m6", "m12")
  do.call(rbind, lapply(Y, function(y) {
    do.call(rbind, lapply(Z, function(z) {
      x1 <- paste0("Δ ", x, " (", z, ")")
      if (!(x1 %in% names(dta))) return(NULL)
      i <- !is.na(dta[[x1]]) & !is.na(dta[[y]])
      n <- sum(i)
      if (sum(i) >= 2) {
        u <- dta[i, x1]
        v <- dta[i, y]
        fit <- lm(v ~ u)
        beta <- coef(fit)
        if (sum(i) > 2) {
          ci <- confint(fit)
          pv <- anova(fit)$`Pr(>F)`[1]
        } else {
          ci <- matrix(NA, nrow = 2, ncol = 2)
          pv <- NA
        }
      } else {
        beta <- NA
        ci <- matrix(NA, nrow = 2, ncol = 2)
        pv <- NA
      }
      r <- c(beta, as.vector(t(ci)))[c(1, 3:4, 2, 5:6)]
      names(r) <- c(paste0("Intercept", c("", ".lwr", ".upr")),
                    paste0("Slope", c("", ".lwr", ".upr")))
      r <- c(r, p.value = pv)
      s <- data.frame(
        outcome = y,
        predictor = x1,
        nobs = n
      )
      cbind(s, t(r))
    }))
  }))
}))
write_xlsx(uni_reg_tbl, file.path(outdir, "Q5_univariable_regressions.xlsx"))

# Univariate regressions - Figures
uni_reg_figs <- mclapply(setNames(X, X), function(x) {
  Y <- paste("Δ Testosterone", c("(nmol/l)", "libre (pmol/l)"), "(m12)")
  Z <- c("j1", "j28", "m3", "m6", "m12")
  figs <- unlist(recursive = FALSE, lapply(Y, function(y) {
    lapply(Z, function(z) {
      x1 <- paste0("Δ ", x, " (", z, ")")
      i <- !is.na(dta[[x1]]) & !is.na(dta[[y]])
      n <- sum(i)
      if (n >= 2) {
        d <- data.frame(u = dta[i, x1], v = dta[i, y])
        fig <- ggplot(d, aes(y = v, x = u)) +
          geom_point() +
          geom_smooth(method = lm, formula = y ~ x, se = FALSE) +
          labs(x = x1, y = y) +
          theme(axis.title = element_text(size = rel(.8)))
      } else {
        fig <- NULL
      }
    })
  }))
  figs[!sapply(figs, is.null)]
})
f1 <- file.path(outdir, "Q5_univariable_regressions.pdf")
cairo_pdf(f1, width = 14, height = 14, onefile = TRUE)
for (figs in uni_reg_figs) {
  do.call(grid.arrange, append(figs, list(ncol = round(sqrt(length(figs))))))
}
dev.off()
f2 <- "/tmp/tmp_bookmarks.txt"
for(i in 1:length(uni_reg_figs)) {
  write("BookmarkBegin", file = f2, append = (i!=1))
  s <- paste("BookmarkTitle:", names(uni_reg_figs)[i])
  write(s, file = f2, append = TRUE)
  write("BookmarkLevel: 1", file = f2, append = TRUE)
  write(paste("BookmarkPageNumber:", i), file = f2, append = TRUE)
}
system(paste("pdftk", f1, "update_info", f2, "output /tmp/tmp_fig.pdf"))
system(paste("mv /tmp/tmp_fig.pdf", f1, "&& rm", f2))
rm(figs, f1, f2, i, s)

# Adjusted regressions - Table
adj_reg_tbl <- do.call(rbind, mclapply(X, function(x) {
  Y <- paste("Testosterone", c("(nmol/l)", "libre (pmol/l)"), "(m12)")
  Z <- c("j1", "j28", "m3", "m6", "m12")
  do.call(rbind, lapply(Y, function(y) {
    do.call(rbind, lapply(Z, function(z) {
      x1 <- paste0("Δ ", x, " (", z, ")")
      x2 <- sub("\\(m12\\)$", "(bl)", y)
      if (!(x1 %in% names(dta)) | !(x2 %in% names(dta))) return(NULL)
      i <- !is.na(dta[[x1]]) & !is.na(dta[[x2]]) & !is.na(dta[[y]])
      n <- sum(i)
      if (sum(i) >= 3) {
        u1 <- dta[i, x1]
        u2 <- dta[i, x2]
        v <- dta[i, y]
        fit <- lm(v ~ u1 + u2)
        beta <- coef(fit)
        if (sum(i) > 3) {
          ci <- confint(fit)
          pv <- coef(summary(fit))[, "Pr(>|t|)"]
        } else {
          ci <- matrix(NA, nrow = 3, ncol = 2)
          pv <- rep(NA, 3)
        }
      } else {
        beta <- NA
        ci <- matrix(NA, nrow = 3, ncol = 2)
        pv <- rep(NA, 3)
      }
      r <- c(beta, as.vector(t(ci)), pv)[c(1, 4:5, 2, 6:7, 11, 3, 8:9, 12)]
      names(r) <- c(paste0("Intercept", c("", ".lwr", ".upr")),
                    paste0("Slope", c("", ".lwr", ".upr", ".pval")),
                    paste0("Slope2", c("", ".lwr", ".upr", ".pval")))
      s <- data.frame(
        outcome = y,
        predictor = x1,
        nobs = n
      )
      cbind(s, t(r))
    }))
  }))
}))
write_xlsx(adj_reg_tbl, file.path(outdir, "Q5_adjusted_regressions.xlsx"))
