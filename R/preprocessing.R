library(readxl)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Phenobese)")

# Import and clean data
p <- "data-raw/BaseDesDonnées_CopieMrPasquier_Sept2.xlsx"
S <- excel_sheets(p)
for (s in S) {
  d <- as.data.frame(read_xlsx(path = p, sheet = s))
  d <- d[apply(!is.na(d), 1, any), ]
  names(d) <- trimws(gsub(" +", " ", gsub("\r\n", " ", names(d))))
  for(v in names(d)[sapply(d, class) == "character"]) {
    x <- d[[v]]
    na <- c("-", "ND")
    b <- is.na(x) | grepl("^[0-9]+(\\.[0-9]+)?$", x) | x %in% na
    if (all(b)) {
      x[x %in% na] <- NA
      d[[v]] <- as.numeric(x)
    }
    if (v == "NAFLD score Interpretation") {
      d[d[[v]] %in% "No Fibrosis", v] <- "No fibrosis"
    }
  }
  if (s == "Baseline") {
    u <- "bl"
  } else {
    u <- tolower(sub("^PostRYGB-", "", s))
  }
  assign(u, d)
}
rm(b, d, na, p, s, S, u, v, x)

# Harmonization of variable names (the names of the baseline are used as a
# reference)
for (z in c("j1", "j28", "m3", "m6", "m12")) {
  d <- get(z)
  names(d)[names(d) == "T/E ratio"] <- "T/E"
  names(d)[names(d) == "LH"] <- "LH (U/l)"
  names(d)[names(d) == "FSH"] <- "FSH (U/l)"
  names(d)[names(d) == "Estradiol"] <- "Estradiol (nmol/l)"
  names(d)[names(d) == "Insulin"] <- "Insuline"
  names(d)[names(d) == "hsCRP"] <- "hs-CRP"
  names(d)[names(d) == "FGF19 AM corr"] <- "FGF19 AM"
  names(d)[names(d) == "FGF19 PM corr"] <- "FGF19 PM"
  names(d)[names(d) == "Hb1AC"] <- "Hb1Ac"
  assign(z, d)
}
rm(z, d)

# Recode `Group` as factor
bl$Group <- factor(bl$Group, c("Lean control", "Non-HH obese", "HH obese"))

# Categorical variables
V <- names(bl)[sapply(bl, class) == "character"]
V <- V[grepl("^(ADAM|Glucose|NAFLD)", V)]
V <- c(V, grep("\\(1=yes, 0=no\\)$", names(bl), value = TRUE))
Z <- c("bl", "j1", "j28", "m3", "m6", "m12")
for (v in V) {
  lvls <- do.call(c, lapply(Z, function(z) unique(get(z)[[v]])))
  lvls <- sort(unique(lvls))
  if (v == "NAFLD score Interpretation") {
    lvls <- c("No fibrosis", "Fibrosis", "Indeterminate")
  }
  for (z in Z) {
    d <- get(z)
    if (v %in% names(d)) d[[v]] <- factor(d[[v]], lvls)
    assign(z, d)
  }
}
rm(V, Z, d, lvls, v, z)

# Save data
if (!dir.exists("data")) dir.create("data")
for (d in c("bl", "j1", "j28", "m3", "m6", "m12")) {
  save(list = d, file = paste0("data/", d, ".rda"), compress = "xz")
}
rm(d)
