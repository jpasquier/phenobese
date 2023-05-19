library(broom)
library(broom.mixed)
library(dplyr)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(here)
library(lme4)
library(lmerTest)
library(parallel)
library(sjPlot)
library(tidyr)
library(writexl)

options(mc.cores = detectCores() - 1)

# Declare script location
i_am("R/longitudinal_analysis.R")

# Outcomes
y_vars <- c("BMI (kg/m2)", "HOMA-IR", "hs-CRP", "Fat percent",
            "Intra-visceral fat mass", "FGF21 AM", "NEFA AM",
            "isoleucine/alloisoleucine", "Leptin AM (ng/ml)",
            "b-hydroxybutyrate", "Adiponectine", "Ceramides - PURE",
            "leucine", "valine", "Testosterone (nmol/l)", "SHBG",
            "Testosterone libre (pmol/l)")
y_vars <- setNames(
  y_vars, tolower(gsub(" - |-| |/", "_", sub(" \\(.+\\)", "", y_vars))))

# Data
dfs <- c("bl", "j1", "j28", "m3", "m6", "m12")
for (df in dfs) {
  load(here("data", paste0(df, ".rda")))
}
dta <- do.call(bind_rows, mclapply(dfs, function(df) {
  get(df) %>%
    select(id = `Subject ID`, any_of(unname(y_vars))) %>%
    mutate(time = df)
})) %>%
  mutate(time = factor(time, dfs)) %>%
  group_by(id) %>%
  filter(n() > 1) %>%
  ungroup()
rm(df, dfs)

# Output directory
outdir <- here("results", paste0("longitudinal_analysis_",
                                 format(Sys.Date(), "%Y%m%d")))
if (!dir.exists(outdir)) dir.create(outdir)

# Longitudinal analyses using mixed linear models
R <- mclapply(y_vars, function(y) {
  d <- dta %>%
    select(id, time, y = !!sym(y)) %>%
    drop_na(everything())
  fit <- lmer(y ~ time + (1 | id), d)
  tbl <- tidy(fit, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(
      term = sub("^time", "", term),
      term = sub("sd__", "SD ", term),
    )
  figs <- list()
  dgp <- plot_model(fit, type = "diag")
  dgp[[2]] <- dgp[[2]]$id
  figs$diag <- ggarrange(plotlist = dgp, nrow = 2, ncol = 2) %>%
    annotate_figure(top = text_grob("Diagnostic plots",
                    face = "bold", size = 16)) %>%
    suppressMessages()
  figs$pred <- emmeans(fit, "time") %>%
    as_tibble() %>%
    ggplot(aes(x = time, y = emmean, group = 1)) +
    geom_point() +
    geom_line(linetype = "dashed") +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) +
    labs(x = NULL, y = y) +
    theme_bw()
  list(fit = fit, tbl = tbl, figs = figs)
})

# Any singular fit?
if (any(sapply(R, function(x) isSingular(x$fit)))) {
  warning("Singular fit(s)")
}

# Export results
mclapply(R, function(r) r$tbl) %>%
  write_xlsx(file.path(outdir, "regression_coefficients.xlsx"))
mclapply(R, function(r) {
  emmeans(r$fit, "time") %>%
    as_tibble() %>%
    rename(mean = emmean) %>%
    select(!df)
}) %>%
  write_xlsx(file.path(outdir, "means.xlsx"))
for (y in names(R)) {
  svg(file.path(outdir, paste0(y, "_diag_plots.svg")))
  print(R[[y]]$figs$diag)
  dev.off()
  svg(file.path(outdir, paste0(y, "_means.svg")))
  print(R[[y]]$figs$pred)
  dev.off()
}

# SessionInfo
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
