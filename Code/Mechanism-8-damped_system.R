

library(ggplot2)

theme_set(theme_classic())
theme_update(panel.border=element_blank(),
             plot.title=element_text(size=rel(8/11),hjust=0.5),
             plot.title.position="plot",
             axis.text=element_text(size=rel(8/11)),
             axis.text.y=element_text(angle=90,hjust=0.5),
             axis.title=element_text(size=rel(9/11)),
             axis.title.x=element_text(hjust=0.5),
             axis.title.y=element_text(hjust=0.5),
             plot.margin=unit(c(0.05,0.15,0.15,0.15),"in"),
             legend.title=element_blank(),
             legend.background = element_rect(colour = NA,fill = NA),
             legend.key.size = unit(8, "pt"),
             legend.key.height = unit(8, "pt"),
             legend.text=element_text(size=rel(7/11))
)


###############################################################

library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(grid)

set.seed(20230703)

#----------------------------------------------------------
# Simulation parameters
#----------------------------------------------------------

n         <- 10
n_sim     <- 100
phi_grid  <- seq(0.05,0.95,by=0.02)
sigma_vec <- c(0.5,1,2)

results <- expand.grid(
  phi_true = phi_grid,
  sigma_eps = sigma_vec,
  phi_est = NA,
  r_sq = NA
)

for(s_idx in seq_along(sigma_vec)){

  sigma_eps <- sigma_vec[s_idx]

  for(i in seq_along(phi_grid)){

    phi_true <- phi_grid[i]

    pool_df <- data.frame(
      V3=numeric(0),
      V2=numeric(0),
      V1=numeric(0)
    )

    for(sim in 1:n_sim){

      V <- numeric(n)

      eps <- rnorm(n,sd=sigma_eps)

      V[1] <- eps[1]
      V[2] <- -phi_true^2*V[1] + eps[2]

      for(t in 3:n){
        V[t] <- -phi_true^2*V[t-1] -
          phi_true^3*V[t-2] +
          eps[t]
      }

      pool_df <- rbind(
        pool_df,
        data.frame(
          V3=V[3:n],
          V2=V[2:(n-1)],
          V1=V[1:(n-2)]
        )
      )
    }

    model <- tryCatch(

      nls(
        V3 ~ -phi^2*V2 - phi^3*V1,
        data=pool_df,
        start=list(phi=phi_true),
        control=nls.control(
          maxiter=200,
          warnOnly=TRUE)
      ),

      error=function(e) NULL
    )

    idx <- which(results$phi_true==phi_true &
                   results$sigma_eps==sigma_eps)

    if(is.null(model)){

      results$phi_est[idx] <- NA
      results$r_sq[idx] <- NA

    }else{

      phi_est <- coef(model)[["phi"]]

      fitted_vals <- predict(model)

      r_sq <- cor(pool_df$V3,fitted_vals)^2

      results$phi_est[idx] <- phi_est
      results$r_sq[idx] <- r_sq
    }
  }
}

#----------------------------------------------------------
# Factors
#----------------------------------------------------------

results$sigma_eps <- factor(
  results$sigma_eps,
  levels=sigma_vec,
  labels=c(
    "sigma[epsilon]==0.5",
    "sigma[epsilon]==1",
    "sigma[epsilon]==2"
  )
)

#----------------------------------------------------------
# Colors
#----------------------------------------------------------

cols <- c(
  "#1b9e77",
  "#d95f02",
  "#7570b3"
)

#----------------------------------------------------------
# Common theme
#----------------------------------------------------------

base_theme <- theme_classic(base_size=14) +
  theme(
    plot.title=element_text(
      hjust=.5,
      size=14
    ),
    axis.title=element_text(size=14),
    axis.text=element_text(size=12)
  )

#----------------------------------------------------------
# Panel A
#----------------------------------------------------------

p1 <- ggplot(
  results,
  aes(phi_true,
      phi_est,
      colour=sigma_eps)
)+

  geom_abline(
    slope=1,
    intercept=0,
    linetype="dashed",
    colour="grey50"
  )+

  geom_line(linewidth=1)+

  geom_point(size=2)+

  scale_color_manual(
    values=cols,
    labels=label_parse()
  )+

  coord_fixed(
    xlim=c(0,1),
    ylim=c(0,1)
  )+

  labs(
    title=" ",
    x=expression(True~phi),
    y=expression(Estimated~hat(phi)),
    colour="Noise level"
  )+

  base_theme+

  theme(

    legend.position=c(.72,.18),

    legend.background=
      element_rect(
        fill=alpha("white",0.8),
        colour="grey70"
      ),

    legend.title=
      element_text(face="bold"),

    legend.key.height=
      unit(.45,"cm"),

    legend.key.width=
      unit(.6,"cm")
  )

#----------------------------------------------------------
# Panel B
#----------------------------------------------------------

p2 <- ggplot(
  results,
  aes(phi_true,
      r_sq,
      colour=sigma_eps,
      linetype=sigma_eps)
)+

  geom_hline(
    yintercept=seq(0,1,.25),
    colour="grey80",
    linetype="dotted",
    linewidth=.5
  )+

  geom_line(linewidth=1)+

  geom_point(size=2)+

  geom_vline(
    xintercept=.618,
    colour="darkred",
    linetype="dashed",
    linewidth=.7
  )+

  annotate(
    "text",
    x=.618,
    y=.95,
    label="phi == 0.618",
    parse=TRUE,
    colour="darkred",
    size=4,
    hjust=-0.1
  )+

  scale_color_manual(
    values=cols,
    labels=label_parse()
  )+

  scale_linetype_discrete(
    labels=label_parse()
  )+

  coord_fixed(
    ylim=c(0,1)
  )+

  labs(
    title=expression(
      V[t]==-phi^2*V[t-1]-phi^3*V[t-2]+epsilon[t]
    ),
    x=expression(True~phi),
    y=expression(Squared~Pearson~correlation~(R^2))
  )+

  base_theme+

  theme(
    legend.position="none"
  )

#----------------------------------------------------------
# Combine
#----------------------------------------------------------

final_plot <- plot_grid(
  p1,
  p2,
  ncol=1,
  align="hv",
  axis="lr",
  rel_heights=c(1,1),
  labels=c("a","b"),
  label_size=16,
  hjust=-0.3
)

#----------------------------------------------------------
# Save
#----------------------------------------------------------

ggsave(
  "E:/canalizationENV/figs/damped_system_three_generation.pdf",
  final_plot,
  width=5,
  height=9.5,
  device="pdf"
)

print(final_plot)


###############################################################

library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
set.seed(20230703)

#----------------------------------------------------------
# Simulation parameters
#----------------------------------------------------------
n         <- 10
n_sim     <- 100
phi_grid  <- seq(0.05, 0.95, by = 0.02)
sigma_vec <- c(0.5, 1, 2)

results <- expand.grid(phi_true = phi_grid, sigma_eps = sigma_vec,
                       phi_est = NA, R2 = NA, R = NA)

for (s_idx in seq_along(sigma_vec)) {
  sigma_eps <- sigma_vec[s_idx]
  
  for (i in seq_along(phi_grid)) {
    phi_true <- phi_grid[i]
    
    pool_df <- data.frame(V3 = numeric(0), V2 = numeric(0), V1 = numeric(0))
    
    for (sim in 1:n_sim) {
      V <- numeric(n)
      eps <- rnorm(n, sd = sigma_eps)
      
      V[1] <- eps[1]
      if (n >= 2) {
        V[2] <- -phi_true^2 * V[1] + eps[2]
      }
      if (n >= 3) {
        for (t in 3:n) {
          coeffs <- phi_true^(2:t)
          V[t] <- -sum(coeffs * rev(V[1:(t-1)])) + eps[t]
        }
      }
      
      triples <- data.frame(
        V3 = V[3:n],
        V2 = V[2:(n-1)],
        V1 = V[1:(n-2)]
      )
      pool_df <- rbind(pool_df, triples)
    }
    
    model <- tryCatch(
      nls(V3 ~ -phi^2 * V2 - phi^3 * V1,
          data = pool_df,
          start = list(phi = phi_true),
          control = nls.control(maxiter = 200, warnOnly = TRUE)),
      error = function(e) NULL
    )
    
    idx <- which(results$phi_true == phi_true & results$sigma_eps == sigma_eps)
    
    if (is.null(model)) {
      results$phi_est[idx] <- NA
      results$R2[idx]       <- NA
      results$R[idx]        <- NA
    } else {
      phi_est <- coef(model)[["phi"]]
      fitted_vals <- predict(model)
      r_val <- cor(pool_df$V3, fitted_vals)
      results$phi_est[idx] <- phi_est
      results$R2[idx]       <- r_val^2
      results$R[idx]        <- r_val
    }
  }
}

results$sigma_eps <- factor(results$sigma_eps,
                            levels = sigma_vec,
                            labels = c("sigma[epsilon]==0.5",
                                       "sigma[epsilon]==1",
                                       "sigma[epsilon]==2"))

cols <- c("#1b9e77", "#d95f02", "#7570b3")
base_theme <- theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))

# ---- Panel A: Parameter recovery ----
pA <- ggplot(results, aes(phi_true, phi_est, colour = sigma_eps)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0.618, colour = "darkred", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 0.618, y = 0.95, label = "phi == 0.618",
           parse = TRUE, colour = "darkred", size = 4, hjust = -0.1) +
  geom_line(linewidth = 0.6) +           # thinner lines
  geom_point(size = 1.0) +               # smaller points
  scale_color_manual(values = cols, labels = label_parse()) +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = " ",
       x = expression(True~phi),
       y = expression(Estimated~hat(phi)),
       colour = "Noise level") +
  theme(legend.position = c(0.25, 0.82),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = "grey70"),
        legend.title = element_text(face = "bold"),
        legend.key.height = unit(0.45, "cm"),
        legend.key.width  = unit(0.6, "cm"))

# ---- Panel B: Squared Pearson correlation (R²) ----
pB <- ggplot(results, aes(phi_true, R2, colour = sigma_eps, linetype = sigma_eps)) +
  geom_hline(yintercept = c(0, 0.25, 0.5, 0.75, 1),
             colour = "grey80", linetype = "dotted", linewidth = 0.5) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.0) +
  geom_vline(xintercept = 0.618, colour = "darkred", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 0.618, y = 0.95, label = "phi == 0.618",
           parse = TRUE, colour = "darkred", size = 4, hjust = -0.1) +
  scale_color_manual(values = cols, labels = label_parse()) +
  scale_linetype_discrete(labels = label_parse()) +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = expression(V[t] == -phi^2 * V[t-1] - phi^3 * V[t-2] + epsilon[t] ~~ "(fitted model)"),
       x = expression(True~phi),
       y = expression(Squared~Pearson~correlation~(R^2))) +
  theme(legend.position = "none")

# ---- Combine ----
final_plot <- plot_grid(pA, pB, 
                        ncol = 1,
                        align = "hv",
                        axis = "lr",
                        rel_heights = c(1, 1),
						rel_widths=c(1,1),
                        labels = c("a", "b"),
                        label_size = 16,
                        hjust = -0.3)

ggsave("E:/canalizationENV/figs/damped_system_all_generation.pdf",
       final_plot, width = 4.5, height = 9, device = "pdf")
print(final_plot)

#################################################

library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(grid)

#----------------------------------------------------------
# THEME
#----------------------------------------------------------
theme_set(theme_classic())

theme_update(
  panel.border = element_blank(),
  plot.title = element_text(size = rel(8/11), hjust = 0.5),
  plot.title.position = "plot",
  axis.text = element_text(size = rel(8/11)),
  axis.text.y = element_text(angle = 0),   # ✔ ensures horizontal y-axis text (global)
  axis.title = element_text(size = rel(9/11)),
  axis.title.x = element_text(hjust = 0.5),
  axis.title.y = element_text(hjust = 0.5),
  plot.margin = unit(c(0.05,0.15,0.15,0.15),"in"),
  legend.background = element_blank(),
  legend.key.size = unit(8, "pt"),
  legend.key.height = unit(8, "pt"),
  legend.text = element_text(size = rel(7/11))
)

set.seed(20230703)



cols <- c("#1b9e77", "#d95f02", "#7570b3")

x_breaks <- c(0, 0.5, 1)
y_breaks <- c(0, 0.25, 0.5, 0.75, 1)

#----------------------------------------------------------
# PANEL A
#----------------------------------------------------------
p1 <- ggplot(results,
             aes(phi_true, phi_est, colour = sigma_eps)) +

  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              colour = "grey50") +

  # ✔ added vertical reference line (same as panel B)
  geom_vline(
    xintercept = 0.618,
    colour = "darkred",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 0.95,
    label = "italic(phi)==0.618",
    parse = TRUE,
    colour = "darkred",
    size = 2.6,
    hjust = 1
  ) +
  geom_line(linewidth = 0.22) +
  geom_point(size = 0.65, stroke = 0.2) +

  scale_color_manual(
    values = cols,
    name = "Noise level",
    labels = c(
      expression(sigma[epsilon] == 0.5),
      expression(sigma[epsilon] == 1),
      expression(sigma[epsilon] == 2)
    )
  ) +

  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +

  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +

  labs(
    x = expression(phi~"(simulation setting)"),
    y = expression("estimated "~phi)
  ) +

  theme(
    legend.position = c(-0.05, 1.15),
    legend.justification = c(0, 1),
    legend.title = element_text(size = rel(6/11)),
    legend.spacing.y = unit(0.12, "cm")
  )

#----------------------------------------------------------
# PANEL B
#----------------------------------------------------------
p2 <- ggplot(results,
             aes(phi_true, R2,
                 colour = sigma_eps,
                 linetype = sigma_eps)) +

  geom_line(linewidth = 0.22) +
  geom_point(size = 0.65, stroke = 0.2) +

  scale_color_manual(values = cols) +
  scale_linetype_discrete() +

  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +

  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +

  geom_vline(
    xintercept = 0.618,
    colour = "darkred",
    linetype = "dashed",
    linewidth = 0.5
  ) +

  annotate(
    "text",
    x = 0.52,
    y = 0.95,
    label = "italic(phi)==0.618",
    parse = TRUE,
    colour = "darkred",
    size = 2.6,
    hjust = 1
  ) +

  labs(
    x = expression(phi~"(simulation setting)"),
    y = expression(R^2),
    title = expression(V[t] == -phi^2*V[t-1] - phi^3*V[t-2] + epsilon[t])
  ) +

  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = rel(5/11), hjust = 0.68)
  )

#----------------------------------------------------------
# COMBINE
#----------------------------------------------------------
final_plot <- plot_grid(
  p1, p2,
  ncol = 2,
  align = "hv",
  axis = "tb",
  labels = c("a", "b"),
  label_size = 12,
  label_fontface = "plain"
)

#----------------------------------------------------------
# SAVE
#----------------------------------------------------------
ggsave(
  "E:/canalizationENV/figs/damped_system_for_publish.pdf",
  final_plot,
  width = 4,
  height = 2.1
)

print(final_plot)

