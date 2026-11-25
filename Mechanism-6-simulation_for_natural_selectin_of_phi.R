



############################################################
## Natural selection of φ via variance minimization
## Finite lineage, population-level selection
## Consistent with conservation-law simulations
############################################################

### only considering the initial velocity deviation of a focal cell 

library(ggplot2)

set.seed(123)

## ---------------------------
## Biological parameters
## ---------------------------
lineage_length <- 10          # C. elegans scale
n_samples      <- 400         # initial noise realizations
population     <- 300         # embryos per generation
n_generations  <- 120

## noise
X0_sd <- 1

## φ mutation (bounded, conservative)
mutation_sd <- 0.01
phi_min     <- 0.2
phi_max     <- 0.85

## selection strength
selection_beta <- 10

## golden ratio
phi_golden <- (sqrt(5) - 1) / 2

## ---------------------------
## Deterministic residual function
## ---------------------------
## This is exactly your proven formula
compute_residual <- function(phi, X0, n) {
  if (n == 0) return(X0)
  X0 * (1 - sum(phi^(2:(n + 1))))
}

## ---------------------------
## Variance-based embryo fitness
## ---------------------------
## Key: variance across different initial noises
variance_fitness <- function(phi) {
  X0 <- rnorm(n_samples, 0, X0_sd)

  residuals <- vapply(
    X0,
    compute_residual,
    numeric(1),
    phi = phi,
    n   = lineage_length
  )

  var(residuals)
}

## ---------------------------
## Initialize population
## ---------------------------
phi_population <- rnorm(population, mean = 0.55, sd = 0.05)
phi_population <- pmin(pmax(phi_population, phi_min), phi_max)

history <- data.frame(
  generation = integer(),
  mean_phi   = numeric(),
  sd_phi     = numeric()
)

## ---------------------------
## Evolution loop
## ---------------------------
for (gen in 1:n_generations) {

  # developmental instability
  variances <- sapply(phi_population, variance_fitness)

  # fitness: lower variance = higher fitness
  fitness <- exp(-selection_beta * variances)
  fitness <- fitness / sum(fitness)

  # selection
  parents <- sample(
    phi_population,
    size = population,
    replace = TRUE,
    prob = fitness
  )

  # mutation
  phi_population <- parents + rnorm(population, 0, mutation_sd)
  phi_population <- pmin(pmax(phi_population, phi_min), phi_max)

  # record
  history <- rbind(
    history,
    data.frame(
      generation = gen,
      mean_phi   = mean(phi_population),
      sd_phi     = sd(phi_population)
    )
  )

  if (gen %% 10 == 0) {
    cat("Gen", gen,
        "Mean φ =", round(mean(phi_population), 4), "\n")
  }
}

## ---------------------------
## Visualization
## ---------------------------
ggplot(history, aes(generation, mean_phi)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = mean_phi - sd_phi,
        ymax = mean_phi + sd_phi),
    fill = "grey70",
    alpha = 0.3
  ) +
  geom_hline(
    yintercept = phi_golden,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  labs(
    x = "Selection generation",
    y = expression(phi),
    title = "Golden ratio φ is uniquely selected by variance minimization",
    subtitle = "Finite lineage instability penalizes both under- and over-correction"
  ) +
  theme_classic()


### considering the noise injectin into the velocity dynamics 
############################################################
## Natural selection of φ under repeated noise injection
## Finite lineage (~10 generations), variance-based fitness
############################################################

library(ggplot2)

set.seed(123)

## ---------------------------
## Biological parameters
## ---------------------------
lineage_length <- 10            # realistic C. elegans scale
population     <- 300
n_generations  <- 200

## noise injection
noise_sd       <- 0.1           # amplitude of injected noise
inject_every   <- 1             # inject noise every N generations
n_replicates   <- 400           # embryos evaluated per φ

## φ mutation
mutation_sd <- 0.01
phi_min     <- 0.2
phi_max     <- 0.85

## selection strength
selection_beta <- 10

## golden ratio
phi_golden <- (sqrt(5) - 1) / 2

## ---------------------------
## Feedback kernel (your proven form)
## ---------------------------
feedback_kernel <- function(phi, noise, dt) {
  if (dt == 0) return(noise)
  noise * (1 - sum(phi^(2:(dt + 1))))
}

## ---------------------------
## Development with noise injection
## ---------------------------
## Returns final residual at generation T
simulate_embryo <- function(phi) {

  noise_gens <- seq(0, lineage_length, by = inject_every)
  noise_vals <- rnorm(length(noise_gens), 0, noise_sd)

  T <- lineage_length
  residual <- 0

  for (k in seq_along(noise_gens)) {
    tk <- noise_gens[k]
    if (tk <= T) {
      residual <- residual +
        feedback_kernel(phi, noise_vals[k], T - tk)
    }
  }

  residual
}

## ---------------------------
## Variance-based fitness
## ---------------------------
variance_fitness <- function(phi) {
  finals <- replicate(n_replicates, simulate_embryo(phi))
  var(finals)
}

## ---------------------------
## Initialize population
## ---------------------------
phi_population <- rnorm(population, mean = 0.55, sd = 0.05)
phi_population <- pmin(pmax(phi_population, phi_min), phi_max)

history <- data.frame(
  generation = integer(),
  mean_phi   = numeric(),
  sd_phi     = numeric()
)

## ---------------------------
## Evolution loop
## ---------------------------
for (gen in 1:n_generations) {

  variances <- sapply(phi_population, variance_fitness)

  fitness <- exp(-selection_beta * variances)
  fitness <- fitness / sum(fitness)

  parents <- sample(
    phi_population,
    size = population,
    replace = TRUE,
    prob = fitness
  )

  phi_population <- parents + rnorm(population, 0, mutation_sd)
  phi_population <- pmin(pmax(phi_population, phi_min), phi_max)

  history <- rbind(
    history,
    data.frame(
      generation = gen,
      mean_phi   = mean(phi_population),
      sd_phi     = sd(phi_population)
    )
  )

  if (gen %% 10 == 0) {
    cat("Gen", gen,
        "Mean φ =", round(mean(phi_population), 4), "\n")
  }
}

## ---------------------------
## Visualization
## ---------------------------
ggplot(history, aes(generation, mean_phi)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = mean_phi - sd_phi,
        ymax = mean_phi + sd_phi),
    fill = "grey70",
    alpha = 0.3
  ) +
  geom_hline(
    yintercept = phi_golden,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  labs(
    x = "Selection generation",
    y = expression(phi),
    title = "Golden ratio φ is selected under repeated noise injection",
    subtitle = paste(
      "Finite lineage (", lineage_length,
      " generations), noise every ", inject_every, " generation",
      sep = ""
    )
  ) +
  theme_classic()



############################################################
## Unified simulation:
## Natural selection of φ via variance minimization
## Finite lineage (~10 generations)
## Supports:
##   (1) initial-noise-only
##   (2) repeated noise injection
############################################################

library(ggplot2)
set.seed(1)

## =========================================================
## Global biological parameters
## =========================================================
lineage_length <- 10            # realistic C. elegans lineage
population     <- 300           # embryos per selection generation
n_generations  <- 200           # evolutionary generations
n_replicates   <- 400           # embryos evaluated per φ

## noise
X0_sd          <- 1.0           # initial noise SD
noise_sd       <- 0.1           # injected noise SD
inject_every   <- 1             # generations between injections

## mutation
mutation_sd <- 0.01
phi_min     <- 0.2
phi_max     <- 0.85

## selection strength
selection_beta <- 10

## golden ratio
phi_golden <- (sqrt(5) - 1) / 2

## =========================================================
## Feedback kernel (proven conservation-law form)
## =========================================================
feedback_kernel <- function(phi, noise, dt) {
  if (dt == 0) return(noise)
  noise * (1 - sum(phi^(2:(dt + 1))))
}

## =========================================================
## Developmental simulation
## mode = "initial"  → noise only at generation 0
## mode = "inject"   → repeated noise injection
## =========================================================
simulate_embryo <- function(phi, mode = c("initial", "inject")) {

  mode <- match.arg(mode)
  T <- lineage_length
  residual <- 0

  if (mode == "initial") {
    X0 <- rnorm(1, 0, X0_sd)
    residual <- feedback_kernel(phi, X0, T)

  } else {
    noise_gens <- seq(1, T, by = inject_every)
    noise_vals <- rnorm(length(noise_gens), 0, noise_sd)

    for (k in seq_along(noise_gens)) {
      tk <- noise_gens[k]
      residual <- residual +
        feedback_kernel(phi, noise_vals[k], T - tk)
    }
  }

  residual
}

## =========================================================
## Variance-based fitness
## =========================================================
variance_fitness <- function(phi, mode) {
  finals <- replicate(n_replicates, simulate_embryo(phi, mode))
  var(finals)
}

## =========================================================
## Evolutionary simulation
## =========================================================
run_selection <- function(mode) {

  phi_population <- rnorm(population, mean = 0.3, sd = 0.05)
  phi_population <- pmin(pmax(phi_population, phi_min), phi_max)

  history <- data.frame(
    generation = integer(),
    mean_phi   = numeric(),
    sd_phi     = numeric()
  )

  for (gen in 1:n_generations) {

    variances <- sapply(phi_population, variance_fitness, mode = mode)

    fitness <- exp(-selection_beta * variances)
    fitness <- fitness / sum(fitness)

    parents <- sample(
      phi_population,
      size = population,
      replace = TRUE,
      prob = fitness
    )

    phi_population <- parents + rnorm(population, 0, mutation_sd)
    phi_population <- pmin(pmax(phi_population, phi_min), phi_max)

    history <- rbind(
      history,
      data.frame(
        generation = gen,
        mean_phi   = mean(phi_population),
        sd_phi     = sd(phi_population)
      )
    )

    if (gen %% 20 == 0) {
      cat("Mode:", mode,
          "| Gen", gen,
          "| Mean φ =", round(mean(phi_population), 4), "\n")
    }
  }

  history
}

## =========================================================
## Run both scenarios for generation = 10
## =========================================================
history_initial <- run_selection("initial")
history_inject  <- run_selection("inject")

history_initial$mode <- "Initial noise only"
history_inject$mode  <- "Repeated noise injection"

history_all <- rbind(history_initial, history_inject)

## =========================================================
## Visualization
## =========================================================
ggplot(history_all, aes(generation, mean_phi, color = mode, fill = mode)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(
    aes(ymin = mean_phi - sd_phi,
        ymax = mean_phi + sd_phi),
    alpha = 0.25,
    color = NA
  ) +
  geom_hline(
    yintercept = phi_golden,
    linetype = "dashed",
    color = "black",
    linewidth = 1
  ) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  labs(
    x = "Selection generation",
    y = expression(phi),
    title = "Golden ratio φ emerges under finite-lineage selection",
    subtitle = "Variance minimization penalizes both under- and over-correction"
  ) +
  theme_classic()



## =========================================================
## Run both scenarios for generation = 13
## =========================================================
lineage_length <- 13            # realistic C. elegans lineage

history_initial_13 <- run_selection("initial")
history_inject_13  <- run_selection("inject")

history_initial_13$mode <- "Initial noise only"
history_inject_13$mode  <- "Repeated noise injection"

history_all_13 <- rbind(history_initial_13, history_inject_13)

## =========================================================
## Visualization
## =========================================================
ggplot(history_all_13, aes(generation, mean_phi, color = mode, fill = mode)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(
    aes(ymin = mean_phi - sd_phi,
        ymax = mean_phi + sd_phi),
    alpha = 0.25,
    color = NA
  ) +
  geom_hline(
    yintercept = phi_golden,
    linetype = "dashed",
    color = "black",
    linewidth = 1
  ) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  labs(
    x = "Selection generation",
    y = expression(phi),
    title = "Golden ratio φ emerges under finite-lineage selection",
    subtitle = "Variance minimization penalizes both under- and over-correction"
  ) +
  theme_classic()



##################################################################



## =========================================================
## Run selection for lineage length = 10 (inject only)
## =========================================================
lineage_length <- 10
history_inject_10 <- run_selection("inject")
history_inject_10$lineage <- "Lineage length = 10"

## =========================================================
## Run selection for lineage length = 13 (inject only)
## =========================================================
lineage_length <- 13
history_inject_13 <- run_selection("inject")
history_inject_13$lineage <- "Lineage length = 13"

## =========================================================
## Combine histories
## =========================================================
history_compare <- rbind(history_inject_10, history_inject_13)

## =========================================================
## Visualization: lineage-length comparison
## =========================================================
ggplot(history_compare, aes(generation, mean_phi, color = lineage, fill = lineage)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(
    aes(ymin = mean_phi - sd_phi, ymax = mean_phi + sd_phi), alpha = 0.25, color = NA) +
  geom_hline(yintercept = phi_golden, linetype = "dashed", color = "black", linewidth = 1) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,0.618,1))+
  labs(
    x = "Selection generation",
    y = expression(phi),
    title = "Golden ratio φ under repeated noise injection",
    subtitle = "Comparison between finite lineage lengths (noise injected every generation)"
  )

ggplot(history_compare, aes(generation, mean_phi, color = lineage, fill = lineage)) +
  geom_line(linewidth = 0.3) +
  geom_ribbon(
    aes(ymin = mean_phi - sd_phi, ymax = mean_phi + sd_phi), alpha = 0.15, color = NA) +
  geom_hline(yintercept = phi_golden, linetype = "dashed", color = "black", linewidth = 0.3) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,0.618,1))+
  labs(
    x = "Selection generation",
    y = "Average phi in population",
    title = "Evolutionary simulations of golden ratio φ"
  )+
  theme(legend.position=c(0.6,0.2))

ggsave("Evolution_optimal_phi_finite_lineage.pdf",width=unit(2.2,"in"),height=unit(1.7,"in"))


### smaller or larger lineage length 


## =========================================================
## Run selection for lineage length = 5 (inject only)
## =========================================================
lineage_length <- 5
history_inject_5 <- run_selection("inject")
history_inject_5$lineage <- "Lineage length = 5"

## =========================================================
## Run selection for lineage length = 20 (inject only)
## =========================================================
lineage_length <- 20
history_inject_20 <- run_selection("inject")
history_inject_20$lineage <- "Lineage length = 20"


## =========================================================
## Combine histories
## =========================================================
history_compare_abnormal <- rbind(history_inject_5, history_inject_20)

## =========================================================
## Visualization: lineage-length comparison
## =========================================================

ggplot(history_compare_abnormal, aes(generation, mean_phi, color = lineage, fill = lineage)) +
  geom_line(linewidth = 0.3) +
  geom_ribbon(
    aes(ymin = mean_phi - sd_phi, ymax = mean_phi + sd_phi), alpha = 0.15, color = NA) +
  geom_hline(yintercept = phi_golden, linetype = "dashed", color = "black", linewidth = 0.3) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,0.618,1))+
  labs(
    x = "Selection generation",
    y = "Average phi in population",
    title = "Evolution of golden ratio phi"
  )+
  theme(legend.position=c(0.6,0.2))

ggsave("Evolution_optimal_phi_finite_lineage_abnormal.pdf",width=unit(2,"in"),height=unit(2,"in"))




#####################################################################################
## ===============================
## Exact finite-lineage variance V_n(phi) — Eq. (7)
## ===============================
Vn <- function(phi, n, sigma2 = 1) {

  if (phi <= 0 || phi >= 1) return(Inf)

  A <- 1 - phi^2 / (1 - phi)

  term1 <- (n + 1) * A^2
  term2 <- (2 * A * phi^2 / (1 - phi)^2) * (1 - phi^(n + 1))
  term3 <- (phi^4 / ((1 - phi)^2 * (1 - phi^2))) *
           (1 - phi^(2 * (n + 1)))

  sigma2 * (term1 + term2 + term3)
}

## ===============================
## Optimization for optimal phi
## ===============================
opt_phi <- function(n, phi_min = 0.2, phi_max = 0.85) {

  opt <- optimize(
    f = function(phi) Vn(phi, n),
    interval = c(phi_min, phi_max)
  )

  data.frame(
    n = n,
    phi_opt = opt$minimum,
    Vmin = opt$objective
  )
}

## ===============================
## Scan lineage lengths
## ===============================
n_values <- 1:25
results <- do.call(rbind, lapply(n_values, opt_phi))

## Golden ratio
phi_golden <- (sqrt(5) - 1) / 2

## Percent of optimal phi relative to golden ratio
results$phi_percent_of_golden <- 100 * results$phi_opt / phi_golden

## ===============================
## Visualization
## ===============================
library(ggplot2)

ggplot(results, aes(x = n, y = phi_opt)) +
  ## C. elegans lineage shadow (n = 10–13)
  annotate(
    "rect",
    xmin = 10, xmax = 13,
    ymin = 0, ymax = 1,
    alpha = 0.15,
    fill = "gray"
  ) +
  ## optimal phi curve
  geom_line(linewidth = 0.3) +
  geom_point(size = 0.5) +
  ## golden ratio reference
  geom_hline(
    yintercept = phi_golden,
    linetype = "dashed",
    color = "red",
    linewidth = 0.3
  ) +
  ## axis settings
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, phi_golden, 1),
    labels = c("0", "0.5", expression(phi), "1")
  ) +
  labs(
    x = "Lineage length",
    y = "Optimum phi",
    title = "Analytic proving and numerical results"
  )

ggsave("finite_lineage_optimal_phi_converges_to_the_golden_ratio.pdf",width=unit(2.2,"in"),height=unit(1.7,"in"))






