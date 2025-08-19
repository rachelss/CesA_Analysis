# LRT comparing M0 (null model) and M1a (nearly neutral)

lnL_M0 <- -59147.35
lnL_M1a <- -58502.65

LR <- -2* (lnL_M0 - lnL_M1a)

degrees_of_freedom <- 1

p_value <- pchisq(LR, df = degrees_of_freedom, lower.tail = FALSE)

cat("Likelihood Ratio:", LR, "\n")
cat("P-value:", p_value, "\n")

# LRT comparing M1a (nearly neutral) and M2a (positive selection)

lnL_M1a <- -58502.65
lnL_M2a <- -58502.65

LR <- -2* (lnL_M1a - lnL_M2a)

degrees_of_freedom <- 2

p_value <- pchisq(LR, df = degrees_of_freedom, lower.tail = FALSE)

cat("Likelihood Ratio:", LR, "\n")
cat("P-value:", p_value, "\n")

curve(dchisq(x, df = degrees_of_freedom), from = 0, to = 10, col = "blue", lwd = 2)
abline(v = LR, col = "red", lwd = 2)  # Add vertical line for LR value

# comparing M7 (neutral/purifying) vs M8 (positive)

lnL_M7 <- -57030.17
lnL_M8 <- -57030.18

LR <- -2* (lnL_M7 - lnL_M8)

degrees_of_freedom <- 2

p_value <- pchisq(LR, df = degrees_of_freedom, lower.tail = FALSE)

cat("Likelihood Ratio:", LR, "\n")
cat("P-value:", p_value, "\n")

curve(dchisq(x, df = degrees_of_freedom), from = 0, to = 10, col = "blue", lwd = 2)
abline(v = LR, col = "red", lwd = 2)  # Add vertical line for LR value
