
# working model 
library(ggplot2)

set.seed(1)

## ---------------------------
## Parameters
## ---------------------------
Tmax <- 50
phi  <- (sqrt(5) - 1) / 2
noise_sd <- 0.1

eps <- rnorm(Tmax, 0, noise_sd)

## ---------------------------
## (1) Geometric distributed lag
## ---------------------------
V_gdl <- numeric(Tmax)

for (t in 2:Tmax) {
  lags    <- V_gdl[1:(t - 1)]
  weights <- phi^(seq_len(t - 1))
  V_gdl[t] <- -phi * sum(weights * rev(lags)) + eps[t]
}

## ---------------------------
## (2) Iterative inherited form
## ---------------------------
V_new <- numeric(Tmax)
V_inh <- numeric(Tmax)

for (t in 2:Tmax) {
  V_new[t] <- -phi * V_inh[t] + eps[t]
  V_inh[t + 1] <- phi * (V_new[t] + V_inh[t])
}

## ---------------------------
## Prepare data for ggplot
## ---------------------------
df <- data.frame(
  t     = rep(1:Tmax, 2),
  V     = c(V_gdl, V_new),
  model = factor(
    rep(c("Geometric distributed lag",
          "Iterative form"),
        each = Tmax),
    levels = c("Geometric distributed lag",
               "Iterative form")
  )
)

## ---------------------------
## Visualization (two panels)
## ---------------------------
ggplot(df, aes(t, V,col=model)) +
  geom_line(linewidth = 0.5,show.legend=F) +
  facet_wrap(~ model, ncol = 1, scales = "free_y") +
  theme(
    strip.background = element_blank()
  )+  
  labs(
    x = "Generation (t)",
    y = expression(V[t]),
    title = "Micro-macro transformation",
  )

ggsave("micro_macro_transformation.pdf",width=unit(3,"in"),height=unit(2.2,"in"))


