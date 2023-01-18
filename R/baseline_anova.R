library(broom)
library(dplyr)
library(here)
library(parallel)
library(tidyr)
library(writexl)

options(mc.cores = detectCores() - 1)

# Declare script location
i_am("R/baseline_anova.R")

# Data
load(here("data/bl.rda"))

# Output directory
outdir <- here("results", paste0("baseline_anova_",
                                 format(Sys.Date(), "%Y%m%d")))
if (!dir.exists(outdir)) dir.create(outdir)

# Anova - Numeric variables
V <- names(bl)[sapply(bl, class) == "numeric"]
V <- V[V != "Calculated age"]
tab_num <- do.call(bind_rows, mclapply(V, function(v) {
  dta <- bl %>%
    select(group = Group, y = !!sym(v)) %>%
    mutate(group = factor(group, c("Lean control", "Non-HH obese", "HH obese"),
                          c("Control", "NonHH", "HH"))) %>%
    drop_na(everything()) 
  d <- dta %>%
    group_by(group) %>%
    summarise(n = n(), mean = mean(y), sd = sd(y)) %>%
    pivot_wider(names_from = group, values_from = c(n, mean, sd),
                names_vary = "slowest") %>%
    mutate(
      anova_pval = anova(lm(y ~ group, dta))$`Pr(>F)`[1],
      kruskal_pval = kruskal.test(y ~ group, data = dta)$p.value
    )
  cbind(variable = v, d)
}))
rm(V)

# Fisher test - Categorical variables
V <- names(bl)[sapply(bl, class) == "factor"]
V <- V[V != "Group"]
tab_cat <- do.call(bind_rows, mclapply(V, function(v) {
  dta <- bl %>%
    select(group = Group, y = !!sym(v)) %>%
    mutate(group = factor(group, c("Lean control", "Non-HH obese", "HH obese"),
                          c("Control", "NonHH", "HH"))) %>%
    drop_na(everything())
  d <- dta %>%
    group_by(group, y) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(prop = n / sum(n)) %>%
    pivot_wider(names_from = group, values_from = c(n, prop),
                names_vary = "slowest") %>%
    rename(value = y)
  d[is.na(d)] <- 0
  p <- fisher.test(dta$y, dta$group)$p.value
  na <- rep(NA, nrow(d) - 1)
  cbind(variable = c(v, na), d, fisher_test_pval = c(p, na))
}))
rm(V)

# Export results
write_xlsx(list(numeric = tab_num, categorial = tab_cat),
           file.path(outdir, "baseline_anova.xlsx"))

# SessionInfo
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
