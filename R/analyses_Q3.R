library(car)
library(ggplot2)
library(gridExtra)
library(parallel)
library(writexl)

options(mc.cores = detectCores() - 1)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Phenobese)")

# Data
load("data/bl.rda")

# Output directory
outdir <- paste0("results/analyses_Q3_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# Rename variables
names(bl) <- trimws(names(bl))
names(bl) <- gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "",
                  names(bl), perl = TRUE)
names(bl) <- gsub(" |/|-|:|\\+", "_", names(bl))
names(bl) <- gsub("_+", "_", names(bl))
names(bl) <- sub("^1_", "One_", names(bl))
names(bl) <- sub("^2_", "Two_", names(bl))
names(bl) <- sub("^3Oxo", "ThreeOxo", names(bl))
names(bl) <- sub("^7Keto", "SevenKeto", names(bl))
names(bl) <- sub("^12a_", "TwelveA_", names(bl))
grep("^[0-9]", names(bl), value = TRUE)

# Predictors
X <- c(
  "Age", "BMI", "FGF21_AM", "Leptin_AM", "ASAT", "ALAT", "NAFLD_score",
  "Glucose", "Insuline", "HOMA_IR", "Hb1Ac", "Cholesterol", "HDL",
  "Triglycerides", "LDL_calculated", "hs_CRP", "Android_Gynoid_ratio",
  "Fat_percent", "Intra_visceral_fat_mass", "ALMI", "Fat_mass_index",
  "Fat_Mass_Total", "Fat_mass_Android", "Lean_mass_index", "Lean_mass_total",
  "b_hydroxybutyrate", "Acetoacetate", "FGF21_PM", "NEFA_AM", "NEFA_PM",
  "Adiponectine", "FGF19_AM", "FGF19_PM", "Leptin_PM", "ThreeOxo_CA_Fast",
  "SevenKeto_DCA_Fast", "SevenKeto_LCA_Fast", "CA_Fast", "CDCA_Fast",
  "DCA_Fast", "GCA_Fast", "GCDCA_Fast", "GLCA_Fast", "GUDCA_Fast", "LCA_Fast",
  "TCA_Fast", "TCDCA_Fast", "TDCA_Fast", "TLCA_Fast", "TUDCA_Fast",
  "UDCA_Fast", "Total_Bile_Acids_Fast", "Total_unconjugated_Fast",
  "Glycine_conjugated_Fast", "Taurine_conjgated_Fast", "Total_conjugated_Fast",
  "Primary_unconjugated_Fast", "Primary_ALL_Fast",
  "Secondary_unconjugated_Fast", "Secondary_ALL_Fast", "TwelveA_OH_Fast",
  "Non_12a_OH_Fast", "TwelveA_Non_12a_Fast", "ThreeOxo_CA_Prand",
  "SevenKeto_DCA_Prand", "SevenKeto_LCA_Prand", "CA_Prand", "CDCA_Prand",
  "DCA_Prand", "GCA_Prand", "GCDCA_Prand", "GLCA_Prand", "GUDCA_Prand",
  "LCA_Prand", "TCA_Prand", "TCDCA_Prand", "TDCA_Prand", "TLCA_Prand",
  "TUDCA_Prand", "UDCA_Prand", "Total_Bile_Acids_Prand",
  "Total_unconjugated_Prand", "Glycine_conjugated_Prand",
  "Taurine_conjgated_Prand", "Total_conjugated_Prand",
  "Primary_unconjugated_Prand", "Primary_ALL_Prand",
  "Secondary_unconjugated_Prand", "Secondary_ALL_Prand", "TwelveA_OH_Prand",
  "Non_12a_OH_Prand", "TwelveA_Non_12a_Prand", "Cer1P_C16_0", "Cer_C14_0",
  "Cer_C16_0", "Cer_C18_1", "Cer_C18_0","Cer_C20_0", "Cer_C22_0", "Cer_C24_1",
  "Cer_C10_0", "Cer_C24_0", "Ceramides_PURE", "Ceramides_C16_18_20",
  "One_DeoxyCer_C24_1", "One_Deoxy_Spa", "DhCer_C16_0", "DhCer_C18_0",
  "DhCer_C20", "DhCer_C22_1_DeoxyCer_C22", "DhCer_C24_1", "DhCer_C24_0",
  "DihydroCeramides_TOTAL", "Dihydroceramides_C16_18_20", "Cer_DhCer_TOTAL",
  "HexCer_C16_0", "HexCer_C18_1", "HexCer_C18_0", "HexCer_C24_1",
  "HexocylCeramides_TOTAL", "LacCer_C12_0", "LacCer_C16_0", "LacCer_C24_1",
  "LacCer_C24_0", "LactocylCeramides_TOTAL", "Sph", "Spa",
  "N_Nervonoyl_1_Deoxy_Spa", "Sphingoid_base_total", "Sph_1P", "SM_C12_0",
  "SM_C16_0", "SM_C18_1", "SM_C18_0", "SM_C24_1", "SM_C24_0",
  "Sphingomyelins_TOTAL", "Sphingomyelins_C16_18_20", "creatinine",
  "phenylalanine", "leucine", "kynurenine", "tryptophan",
  "isoleucine_alloisoleucine", "methionine", "taurine", "valine", "proline",
  "pipecolate", "tyrosine", "Alpha_aminobutyrate", "beta_alanine", "creatine",
  "alanine", "hydroxyproline", "guanidinoacetate", "threonine",
  "Two_aminoadipate", "glycine","glutamate", "serine", "glutamine",
  "asparagine", "homocitrulline", "citrulline", "aspartate", "arginine",
  "histidine", "lysine", "ornithine"
)
all(X %in% names(bl))

# Outcomes
Y <- c("Testosterone", "Testosterone_libre", "T_E")

# Univariable regressions - Table
uni_reg_tbl <- do.call(rbind, mclapply(X, function(x) {
  do.call(rbind, lapply(Y, function(y) {
    do.call(rbind, lapply(1:2, function(k) {
      d <- if (k == 1) bl else subset(bl, grepl("obese", Group))
      fit <- lm(as.formula(paste(y, "~", x)), d)
      r <- c(coef(fit), as.vector(t(confint(fit))))[c(1, 3:4, 2, 5:6)]
      names(r) <- c(paste0("Intercept", c("", ".lwr", ".upr")),
                    paste0("Slope", c("", ".lwr", ".upr")))
      r <- c(r, p.value = anova(fit)$`Pr(>F)`[1])
      s <- data.frame(
        outcome = y,
        predictor = x,
        subset = c("All", "Obese")[k],
        nobs = nrow(fit$model)
      )
      cbind(s, t(r))
    }))
  }))
}))

# Univariate regressions - Figures
uni_reg_figs <- mclapply(setNames(X, X), function(x) {
  unlist(recursive = FALSE, lapply(setNames(Y, Y), function(y) {
    lapply(c(All = 1, Obese = 2), function(k) {
      d <- if (k == 1) bl else subset(bl, grepl("obese", Group))
      d <- na.omit(d[c(x, y)])
      fit <- lm(as.formula(paste(y, "~", x)), d)
      b <- signif(coef(fit)[[2]], 3)
      ci <- signif(confint(fit)[2, ], 3)
      p <- coef(summary(fit))[2, "Pr(>|t|)" ]
      p <- if (p >= 0.001) paste0("p=", round(p, 3)) else "p<0.001"
      r2 <- paste0("R2=", round(summary(fit)$r.squared, 3))
      cap <- paste0("b = ", b, " (", ci[1], ",", ci[2], "), ", p, ", ", r2)
      fig <- ggplot(d, aes_string(y = y, x = x)) +
        geom_point() +
        geom_smooth(method = lm, formula = y ~ x, se = FALSE) +
        labs(subtitle = c("All", "Obese")[k], caption = cap) +
        theme(axis.title=element_text(size = rel(.75)))
    })
  }))
})

# Multivariate regressions
multi_reg <- mclapply(setNames(Y, Y), function(y) {
  lapply(c(All = 1, Obese = 2), function(k) {
    dta <- if (k == 1) bl else subset(bl, grepl("obese", Group))
    lapply(1:4, function(l) {
      fml <- as.formula(paste(y, "~ Age + HOMA_IR + hs_CRP +",
                              "Intra_visceral_fat_mass"))
      if (l %in% c(2, 4)) fml <- update(fml, . ~ . + BMI)
      if (l %in% 3:4) fml <- update(fml, . ~ . + FGF21_AM)
      fit <- do.call("lm", list(formula = fml, data = quote(dta)))
      tbl <- cbind(data.frame(variable = names(coef(fit)), beta = coef(fit)),
                   confint(fit), `p-value` = coef(summary(fit))[, 4])
      tbl0 <- do.call(rbind, lapply(names(fit$model)[-1], function(x) {
        fit0 <- lm(paste(y, "~", x), fit$model)
        cbind(data.frame(variable = names(coef(fit0)), beta = coef(fit0)),
              confint(fit0), `p-value` = coef(summary(fit0))[, 4])[-1, ]
      }))
      names(tbl0)[-1] <- paste(names(tbl0)[-1], "(univar)")
      vif <- vif(fit)
      vif <- data.frame(variable = names(vif), vif = vif)
      tbl$dummy_row_number <- 1:nrow(tbl)
      tbl <- merge(tbl0, tbl, by = "variable", all = TRUE)
      tbl <- merge(tbl, vif, by = "variable", all = TRUE)
      tbl <- tbl[order(tbl$dummy_row_number), ]
      tbl$dummy_row_number <- NULL
      names(tbl)[names(tbl) == "variable"] <- "independent_variable"
      tbl <- cbind(dependent_variable = c(y, rep(NA, nrow(tbl) - 1)),
                   group = c(c("All", "Obese")[k], rep(NA, nrow(tbl) - 1)),
                   nobs = c(nrow(fit$model), rep(NA, nrow(tbl) - 1)),
                   tbl)
      list(fit = fit, tbl = tbl, n = nrow(fit$model))
    })
  })
})

# Export table of univariate regressions
write_xlsx(uni_reg_tbl, file.path(outdir, "Q3_univariable_regressions.xlsx"))

# Export figures of univariate regressions
f1 <- file.path(outdir, "Q3_univariable_regressions.pdf")
pdf(f1)
for (figs in uni_reg_figs) do.call(grid.arrange, append(figs, list(ncol = 2)))
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

# Export tables and diagnostic plots of multivariable regressions
get_fml <- function(fit) {
  paste(as.character(fit$call$formula)[c(2, 1, 3)], collapse = " ")
}
for (y in names(multi_reg)) {
  for (s in names(multi_reg[[y]])) {
    f <- paste0("Q3_multivariable_regression_", y, "_", s, ".xlsx")
    f <- file.path(outdir, f)
    tbls <- lapply(multi_reg[[y]][[s]], function(z) z$tbl)
    names(tbls) <- paste("Model", 1:length(tbls))
    write_xlsx(tbls, f)
    f <- sub("xlsx$", "pdf", f)
    fits <- lapply(multi_reg[[y]][[s]], function(z) z$fit)
    pdf(f)
    for(fit in fits) {
      par(mfrow = c(2, 2))
      for (i in 1:4) plot(fit, i)
      par(mfrow = c(1, 1))
      mtext(get_fml(fit), outer = TRUE, line = -1.8, cex = 1)
    }
    dev.off()
  }
}
rm(get_fml, y, s, f, tbls, fits, fit)

# SessionInfo
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
