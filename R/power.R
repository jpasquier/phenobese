n_welch_t_test <- function(m1, s1, m2, s2, r = 1, alpha = .05, pwr = .8) {
  root_func <- function(n1) {
    n2 <- r * n1
    d <- m2 - m1
    s <- sqrt(s1^2 / n1 + s2^2 / n2)
    ncp = abs(d) / s
    nu <- (s1^2 / n1 + s2^2 / n2)^2 /
      (s1^4 / (n1^2 * (n1 - 1)) + s2^4 / (n2^2 * (n2 - 1)))
    qu <- qt(1 - alpha / 2, nu)
    1 - pt(qu, nu, ncp = ncp) + pt(-qu, nu, ncp = ncp) - pwr
  }
  n1 <- uniroot(root_func, c(2, 10^7))$root
  n2 <- r * n1
  c(n1 = n1, n2 = n2)
}
n_welch_t_test(m1 = 3.49, s1 = 2.15, m2 = 5.72, s2 = 2.29, r = 1.5)
#       n1       n2
# 13.73042 20.59563
