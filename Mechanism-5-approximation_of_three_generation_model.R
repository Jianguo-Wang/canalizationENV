############################################################
## Bias-corrected simulation:
## GDL process vs 3-generation approximation
############################################################

#set.seed(123)
############################################################
## Simulation: test OLS and NLS estimation using only 
## three-generation trios from full GDL model
############################################################

set.seed(123)
library(ggplot2)

## ---------------------------
## True parameters
## ---------------------------
phi_true <- (sqrt(5) - 1) / 2
alpha <- phi_true
beta  <- phi_true

true_c1 <- -alpha * beta
true_c2 <- -alpha * beta^2

cat("Theoretical coefficients:\n")
cat("V(t-1):", round(true_c1, 3), "\n")
cat("V(t-2):", round(true_c2, 3), "\n\n")

## ---------------------------
## Simulation parameters
## ---------------------------
n_embryos <- 105
lineages_per_embryo <- 30
max_gen <- 10          # finite lineage length
burn_in <- 20          # to stabilize GDL kernel
noise_sd <- 0.2

## ---------------------------
## Function to simulate one lineage (full GDL)
## ---------------------------
simulate_lineage <- function(phi, T, burn_in, noise_sd) {
  total_T <- T + burn_in
  V <- numeric(total_T)
  V[1] <- rnorm(1)
  
  for (t in 2:total_T) {
    feedback <- sum(phi^(1:(t - 1)) * V[(t - 1):1])
    V[t] <- -phi * feedback + rnorm(1, 0, noise_sd)
  }
  
  # Return only biologically observable window
  V[(burn_in + 1):total_T]
}

## ---------------------------
## Generate pooled data (all embryos)
## ---------------------------
V_t <- numeric()
V_t1 <- numeric()
V_t2 <- numeric()

for (emb in 1:n_embryos) {
  for (lin in 1:lineages_per_embryo) {
    V <- simulate_lineage(phi_true, max_gen, burn_in, noise_sd)
    
    # Only extract 3-generation trios for regression
    for (t in 3:max_gen) {
      V_t  <- c(V_t,  V[t])
      V_t1 <- c(V_t1, V[t - 1])
      V_t2 <- c(V_t2, V[t - 2])
    }
  }
}

df <- data.frame(V_t, V_t1, V_t2)
cat("Total three-generation trios:", nrow(df), "\n\n")

## ---------------------------
## 1) OLS estimation (linear model)
## ---------------------------
fit_ols <- lm(V_t ~ V_t1 + V_t2 - 1, data = df)
summary(fit_ols)

est_ols <- coef(fit_ols)

## ---------------------------
## 2) NLS estimation (α = β = phi)
## ---------------------------
# Define residual function using only V_t, V_t1, V_t2
nls_fit <- nls(
  V_t ~ -phi^2 * V_t1 - phi^3 * V_t2,
  start = list(phi = 0.5),
  algorithm = "port",
  lower = 0,
  upper = 1,
  data = df
)

phi_est <- coef(nls_fit)["phi"]

## ---------------------------
## Compare results
## ---------------------------
comparison <- data.frame(
  Method = c("OLS V(t-1)", "OLS V(t-2)", "NLS phi"),
  Estimated = c(round(est_ols[1],3), round(est_ols[2],3), round(phi_est,3)),
  Theoretical = c(round(true_c1,3), round(true_c2,3), round(phi_true,3))
)
print(comparison)

## ---------------------------
## Optional plot: OLS vs theory
## ---------------------------
plot_df <- data.frame(
  Term = factor(c("V(t-1)", "V(t-2)")),
  Estimated = est_ols,
  Theoretical = c(true_c1, true_c2)
)

ggplot(plot_df, aes(x = Theoretical, y = Estimated, label = Term)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(vjust = -1) +
  theme_classic() +
  labs(
    title = "Three-generation trios: OLS recovers coefficients, NLS recovers phi",
    subtitle = paste("NLS phi =", round(phi_est,3), "| True phi =", round(phi_true,3)),
    x = "Theoretical",
    y = "Estimated"
  )


#####################################################################################

############################################################
## Three-generation simulation WITHOUT burn-in
## Mimicking real 10-generation embryo lineages
############################################################

set.seed(123)
library(ggplot2)

## ---------------------------
## True parameters
## ---------------------------
phi_true <- (sqrt(5) - 1) / 2
alpha <- phi_true
beta  <- phi_true

true_c1 <- -alpha * beta
true_c2 <- -alpha * beta^2

cat("Theoretical coefficients:\n")
cat("V(t-1):", round(true_c1, 3), "\n")
cat("V(t-2):", round(true_c2, 3), "\n\n")

## ---------------------------
## Simulation parameters
## ---------------------------
n_embryos <- 105
lineages_per_embryo <- 30
max_gen <- 10       # real embryo lineage length
noise_sd <- 0.5     # realistic noise

## ---------------------------
## Function to simulate one lineage (full GDL)
## ---------------------------
simulate_lineage <- function(phi, T, noise_sd) {
  V <- numeric(T)
  V[1] <- rnorm(1)
  
  for (t in 2:T) {
    feedback <- sum(phi^(1:(t - 1)) * V[(t - 1):1])
    V[t] <- -phi * feedback + rnorm(1, 0, noise_sd)
  }
  
  V
}

## ---------------------------
## Generate pooled data (all embryos)
## ---------------------------
V_t <- numeric()
V_t1 <- numeric()
V_t2 <- numeric()

for (emb in 1:n_embryos) {
  for (lin in 1:lineages_per_embryo) {
    V <- simulate_lineage(phi_true, max_gen, noise_sd)
    
    # Only extract three-generation trios for regression
    for (t in 6:9) {
      V_t  <- c(V_t,  V[t])
      V_t1 <- c(V_t1, V[t - 1])
      V_t2 <- c(V_t2, V[t - 2])
    }
  }
}

df <- data.frame(V_t, V_t1, V_t2)
cat("Total three-generation trios:", nrow(df), "\n\n")

## ---------------------------
## 1) OLS estimation (linear model)
## ---------------------------
fit_ols <- lm(V_t ~ V_t1 + V_t2 - 1, data = df)
summary(fit_ols)

est_ols <- coef(fit_ols)

## ---------------------------
## 2) NLS estimation (α = β = phi)
## ---------------------------
nls_fit <- nls(
  V_t ~ -phi^2 * V_t1 - phi^3 * V_t2,
  start = list(phi = 0.618),
  algorithm = "port",
  lower = 0,
  upper = 1,
  data = df
)

phi_est <- coef(nls_fit)["phi"]

## ---------------------------
## Compare results
## ---------------------------
comparison <- data.frame(
  Method = c("OLS V(t-1)", "OLS V(t-2)", "NLS phi"),
  Estimated = c(round(est_ols[1],3), round(est_ols[2],3), round(phi_est,3)),
  Theoretical = c(round(true_c1,3), round(true_c2,3), round(phi_true,3))
)
print(comparison)

## ---------------------------
## Optional plot: OLS vs theory
## ---------------------------
plot_df <- data.frame(
  Term = factor(c("V(t-1)", "V(t-2)")),
  Estimated = est_ols,
  Theoretical = c(true_c1, true_c2)
)

ggplot(plot_df, aes(x = Theoretical, y = Estimated, label = Term)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(vjust = -1) +
  theme_classic() +
  labs(
    title = "Three-generation trios: OLS and NLS without burn-in",
    subtitle = paste("NLS phi =", round(phi_est,3), "| True phi =", round(phi_true,3)),
    x = "Theoretical",
    y = "Estimated"
  )



