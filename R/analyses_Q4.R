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
V <- intersect(names(bl)[sapply(bl, class) == "numeric"],
               names(m12)[sapply(m12, class) == "numeric"])

# Join the to dataset
dta <- merge(bl[c("Subject ID", "Group", V)], m12[c("Subject ID", V)],
             by = "Subject ID", suffixes = c(" bl", " m12"))

# Table
tbl <- do.call(rbind, mclapply(V, function(v) {
  d <- na.omit(dta[c("Group", paste(v, c("bl", "m12")))])
  d[[paste(v, "diff")]] <- d[[paste(v, "m12")]] - d[[paste(v, "bl")]]
  Z <- c("n", "bl", "m12", "diff")
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
        if (z == "diff") {
          if (length(x) >= 2) {
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
  if (all(table(d$Group) >= 2)) {
    tt.pv <- t.test(d[[paste(v, "diff")]] ~ d$Group)$p.value
    wt.pv <- wilcox.test(d[[paste(v, "diff")]] ~ d$Group,
                         exact = FALSE)$p.value
  } else {
    tt.pv = NA
    wt.pv = NA
  }
  sqrt(r[["diff.HHobese.sd"]]^2 / r[["n.HHobese"]] +
               r[["diff.NonHHobese.sd"]]^2 / r[["n.NonHHobese"]])
  r <- c(r, c(
    diff.diff.mean = r[["diff.HHobese.mean"]] - r[["diff.NonHHobese.mean"]],
    diff.diff.se = sqrt(r[["diff.HHobese.sd"]]^2 / r[["n.HHobese"]] +
                          r[["diff.NonHHobese.sd"]]^2 / r[["n.NonHHobese"]]),
    diff.diff.t.test.pv = tt.pv,
    diff.diff.wilcox.test.pv = wt.pv
  ))
  cbind(data.frame(variable = v), t(r))
}))
write_xlsx(tbl, file.path(outdir, "Q4_table.xlsx"))
