library(broom)
library(dplyr)
library(here)
library(parallel)
library(writexl)

options(mc.cores = detectCores() - 1)

# Declare script location
i_am("R/adjusted_baseline_differences.R")

# Data
load(here("data/bl.rda"))

# Output directory
outdir <- here("results", paste0("adjusted_baseline_differences_",
                                 format(Sys.Date(), "%Y%m%d")))
if (!dir.exists(outdir)) dir.create(outdir)

# Convert yes/no variables to logical
for (j in grep("\\(1=yes, 0=no\\)", names(bl))) {
  if (!all(bl[[j]] %in% c(0:1, NA))) stop("not a logical var")
  if (class(bl[[j]]) == "factor") bl[[j]] <- as.numeric(as.character(bl[[j]]))
  bl[[j]] <- as.logical((bl[[j]]))
}
for (j in which(sapply(bl, class) == "factor")) {
  x <- bl[[j]]
  x <- as.character(x)
  if (all(x %in% c("Negative (-)", "Positive (+)", "no", "yes", NA))) {
    x <- recode(x, `Negative (-)` = "0", `no` = "0",  `Positive (+)` = "1",
                `yes` = "1")
    bl[[j]] <- as.logical(as.numeric(x))
  }
}
rm(j, x)

# Adjusted differences and adjusted odds ratios
V <- names(bl)[sapply(bl, class) %in% c("numeric", "logical")]
V <- V[!(V %in% c("Calculated age", "Age (yrs)", "BMI (kg/m2)"))]
tab <- do.call(bind_rows, mclapply(V, function(v) {
  type <- class(bl[[v]])
  dta <- bl %>%
    filter(Group %in% c("Non-HH obese", "HH obese")) %>%
    mutate(Group = droplevels(Group)) %>%
    select(group = Group, age = `Age (yrs)`, bmi = `BMI (kg/m2)`, y = !!sym(v))
  fml <- y ~ group
  do.call(bind_rows, lapply(1:3, function(k) {
    if (k == 2) {
      fml <- update(fml, . ~ . + age)
    } else if (k == 3) {
      fml <- update(fml, . ~ . + bmi)
    }
    fit <- if (type == "numeric") {
      lm(fml, dta)
    } else {
      glm(fml, dta, family = "binomial")
    }
    tidy(fit, conf.int = TRUE, exponentiate = (type == "logical")) %>%
      mutate(
        variable = v,
        type = type,
        adjustment = c("none", "age", "bmi")[k],
        n = nrow(dta)
      ) %>%
      filter(grepl("^group", term)) %>%
      select(variable, type, adjustment, n, estimate, std.error, conf.low,
             conf.high, p.value)
  }))
}))
rm(V)

# Export results
write_xlsx(tab, file.path(outdir, "adjusted_baseline_differences.xlsx"))

# SessionInfo
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
