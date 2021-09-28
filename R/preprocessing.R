library(readxl)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Hypogonadisme)")

# Import and clean data
p <- "data-raw/BaseDesDonnées_CopieMrPasquier_Sept2.xlsx"
S <- excel_sheets(p)
for (s in S) {
  d <- as.data.frame(read_xlsx(path = p, sheet = s))
  d <- d[apply(!is.na(d), 1, any), ]
  names(d) <- gsub("\r\n", " ", names(d))
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

# Check character variables
table(sapply(bl, class))
lapply(bl[sapply(bl, class) == "character"], table, useNA = "ifany")
lapply(m12[sapply(m12, class) == "character"], table, useNA = "ifany")

# Save data
if (!dir.exists("data")) dir.create("data")
for (d in c("bl", "j1", "j28", "m3", "m6", "m12")) {
  save(list = d, file = paste0("data/", d, ".rda"), compress = "xz")
}
rm(d)
