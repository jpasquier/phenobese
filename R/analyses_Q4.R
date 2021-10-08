library(parallel)
library(writexl)

options(mc.cores = detectCores() - 1)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Phenobese)")

# Data
load("data/bl.rda")
load("data/m12.rda")

# Output directory
outdir <- paste0("results/analyses_Q4_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# Check that groups are well defined
tmp <- merge(bl[c("Subject ID", "Group")], m12[c("Subject ID", "Group")],
             by = "Subject ID")
with(tmp, all(Group.x == Group.y))
rm(tmp)

# Define the set of variables to compare
V_num <- intersect(names(bl)[sapply(bl, class) == "numeric"],
                   names(m12)[sapply(m12, class) == "numeric"])
V_cat <- intersect(names(bl)[sapply(bl, class) == "factor"],
                   names(m12)[sapply(m12, class) == "factor"])

# Join the two datasets
dta <- merge(bl[c("Subject ID", "Group", V_num, V_cat)],
             m12[c("Subject ID", V_num, V_cat)],
             by = "Subject ID", suffixes = c(" bl", " m12"))
dta$Group <- droplevels(dta$Group)

# Table (num)
tbl_num <- do.call(rbind, mclapply(V_num, function(v) {
  d <- na.omit(dta[c("Group", paste(v, c("bl", "m12")))])
  d[[paste(v, "diff")]] <- d[[paste(v, "m12")]] - d[[paste(v, "bl")]]
  d[[paste(v, "rdiff")]] <- d[[paste(v, "m12")]] / d[[paste(v, "bl")]] - 1
  Z <- c("n", "bl", "m12", "diff", "rdiff")
  Z <- setNames(Z, Z)
  G <- c("Obese", "Non-HH obese", "HH obese")
  G <- setNames(G, gsub("-| ", "", G))
  r <- unlist(recursive = TRUE, lapply(Z, function(z) {
    lapply(G, function(g) {
      if (g != "Obese") d <- d[d$Group == g, ]
      if (z == "n") {
        r <- nrow(d)
      } else {
        x <- d[[paste(v, z)]]
        r <- list(mean = mean(x), sd = sd(x))
        if (z %in% c("diff", "rdiff")) {
          if (length(unique(x)) >= 2) {
            tt.pv <- t.test(x)$p.value
            wt.pv <- wilcox.test(x, exact = FALSE)$p.value
          } else {
            tt.pv = NA
            wt.pv = NA
          }
          r <- append(r, list(paired.t.test.pv = tt.pv,
                              paired.wilcox.test.pv = wt.pv))
        }
      }
      return(r)
    })
  }))
  for (u in c("diff", "rdiff")) {
    if (all(table(d$Group) >= 2)) {
      tt.pv <- t.test(d[[paste(v, u)]] ~ d$Group)$p.value
      wt.pv <- wilcox.test(d[[paste(v, u)]] ~ d$Group, exact = FALSE)$p.value
    } else {
      tt.pv = NA
      wt.pv = NA
    }
    s <- c(
      mean = r[[paste0(u, ".HHobese.mean")]] -
        r[[paste0(u, ".NonHHobese.mean")]],
      se = sqrt(r[[paste0(u, ".HHobese.sd")]]^2 / r[["n.HHobese"]] +
                  r[[paste0(u, ".NonHHobese.sd")]]^2 / r[["n.NonHHobese"]]),
      t.test.pv = tt.pv,
      wilcox.test.pv = wt.pv
    )
    names(s) <- paste("diff", u, names(s), sep = ".")
    r <- c(r, s)
  }
  cbind(data.frame(variable = v), t(r))
}))

# Table (cat)
tbl_cat <- do.call(rbind, mclapply(V_cat, function(v) {
  d <- na.omit(dta[c("Group", paste(v, c("bl", "m12")))])
  if (nrow(d) == 0) return(NULL)
  G <- c("Obese", "Non-HH obese", "HH obese")
  G <- setNames(G, gsub("-| ", "", G))
  Merge <- function(x, y) merge(x, y, by = "value", all = TRUE, sort = FALSE)
  r <- Reduce(Merge, mclapply(G, function(g) {
    if (g != "Obese") d <- d[d$Group == g, ]
    r <- Reduce(Merge, lapply(c("bl", "m12"), function(z) {
      u <- table(d[[paste(v, z)]])
      u <- data.frame(value = names(u), n = as.vector(u), prop = as.vector(prop.table(u)))
      names(u)[2:3] <- paste(names(u)[2:3], z, sep = ".")
      u
    }))
    pv <- mcnemar.test(d[[paste(v, "bl")]], d[[paste(v, "m12")]])$p.value
    r <- cbind(r, mcnemar.test.pv = c(pv, rep(NA, nrow(r) - 1)))
    names(r)[-1] <- paste(gsub("-| ", "", g), names(r)[-1], sep = ".")
    return(r)
  }))
  r <- cbind(variable = c(v, rep(NA, nrow(r) - 1)), r)
}))

# Export
write_xlsx(list(numeric = tbl_num, categorical = tbl_cat), file.path(outdir, "Q4_tables.xlsx"))
