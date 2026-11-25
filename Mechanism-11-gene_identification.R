# ========================================================================
# 0. Global ggplot2 theme settings
# ========================================================================
library(ggplot2)

theme_set(theme_classic())
theme_update(
  panel.border          = element_blank(),
  plot.title            = element_text(size = rel(8 / 11), hjust = 0.5),
  plot.title.position   = "plot",
  axis.text             = element_text(size = rel(8 / 11)),
  axis.text.y           = element_text(angle = 90, hjust = 0.5),
  axis.title            = element_text(size = rel(9 / 11)),
  axis.title.x          = element_text(hjust = 0.5),
  axis.title.y          = element_text(hjust = 0.5),
  plot.margin           = unit(c(0.05, 0.15, 0.15, 0.15), "in"),
  legend.title          = element_blank(),
  legend.background     = element_rect(colour = NA, fill = NA),
  legend.key.size       = unit(8, "pt"),
  legend.key.height     = unit(8, "pt"),
  legend.text           = element_text(size = rel(7 / 11))
)

# ========================================================================
# 1. Robust phi estimation function (tries multiple start values)
# ========================================================================
estimate_phi <- function(GD, D, M, phi_start = 0.6) {
  df <- data.frame(GD = GD, D = D, M = M)

  # Candidate starting values – you can adjust this list
  start_vals <- c(phi_start, 0.01, 0.1, 0.2,0.3, 0.4,0.5, 0.6,0.7, 0.8,0.9, 1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1)

  for (st in start_vals) {
    fit <- tryCatch({
      nls(
        GD ~ -phi^2 * D - phi^3 * M,
        data = df,
        start = list(phi = st),
        algorithm = "port",
        lower = 0,
        upper = 2
      )
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      return(list(
        phi_hat   = coef(fit),
        se        = sqrt(diag(vcov(fit))),
        fitted    = fitted(fit),
        residuals = resid(fit)
      ))
    }
  }

  warning("No starting value converged for embryo.")
  return(list(phi_hat = NA, se = NA, fitted = NA, residuals = NA))
}

# ========================================================================
# 2. Meta‑analysis function (returns both fixed and random effects)
# ========================================================================
meta_phi <- function(phi, se) {
  if (length(phi) != length(se)) stop("phi and se must have the same length")

  k <- length(phi)

  # Fixed‑effect
  w_fixed    <- 1 / (se^2)
  phi_fixed  <- sum(w_fixed * phi) / sum(w_fixed)
  se_fixed   <- sqrt(1 / sum(w_fixed))
  ci_fixed   <- c(phi_fixed - 1.96 * se_fixed, phi_fixed + 1.96 * se_fixed)

  # Heterogeneity
  Q  <- sum(w_fixed * (phi - phi_fixed)^2)
  df <- k - 1
  C  <- sum(w_fixed) - sum(w_fixed^2) / sum(w_fixed)
  tau2 <- max(0, (Q - df) / C)

  # Random‑effects
  w_random   <- 1 / (se^2 + tau2)
  phi_random <- sum(w_random * phi) / sum(w_random)
  se_random  <- sqrt(1 / sum(w_random))
  ci_random  <- c(phi_random - 1.96 * se_random, phi_random + 1.96 * se_random)

  I2 <- max(0, (Q - df) / Q)

  return(list(
    fixed = list(
      estimate = phi_fixed,
      se       = se_fixed,
      CI95     = ci_fixed
    ),
    random = list(
      estimate = phi_random,
      se       = se_random,
      CI95     = ci_random,
      tau2     = tau2
    ),
    heterogeneity = list(
      Q  = Q,
      df = df,
      I2 = I2
    )
  ))
}

# ========================================================================
# 3. Prepare RNAi and Env combined data frames (velocity derivatives)
# ========================================================================

## RNAi perturbation (embryos 106 – 1928)
# x
RNAi_combined_df_x_all <- na.omit(do.call(rbind, lapply(106:1928, function(emb) {
  lineage_velocity_emb <- (matrix(
    velocity[emb, match(paste(c(lineage_info), "ccl", sep = "_"), colnames(ccl_raw_nonNA_scaling))],
    nrow = dim(lineage_info)[1]
  ) - lineage_velocity_x_avg) / lineage_velocity_x_sd

  M_D_GD_df <- data.frame(
    M        = c(lineage_velocity_emb[1:8, ]),
    D        = c(lineage_velocity_emb[2:9, ]),
    GD       = c(lineage_velocity_emb[3:10, ]),
    embryoID = rep(emb, length(c(lineage_velocity_emb[1:8, ])))
  )
  return(M_D_GD_df)
})))

# y
RNAi_combined_df_y_all <- na.omit(do.call(rbind, lapply(106:1928, function(emb) {
  lineage_velocity_emb <- (matrix(
    velocity[emb, match(paste(c(lineage_info), "ccl", sep = "_"), colnames(ccl_raw_nonNA_scaling)) + 353],
    nrow = dim(lineage_info)[1]
  ) - lineage_velocity_y_avg) / lineage_velocity_y_sd

  M_D_GD_df <- data.frame(
    M        = c(lineage_velocity_emb[1:8, ]),
    D        = c(lineage_velocity_emb[2:9, ]),
    GD       = c(lineage_velocity_emb[3:10, ]),
    embryoID = rep(emb, length(c(lineage_velocity_emb[1:8, ])))
  )
  return(M_D_GD_df)
})))

# z
RNAi_combined_df_z_all <- na.omit(do.call(rbind, lapply(106:1928, function(emb) {
  lineage_velocity_emb <- (matrix(
    velocity[emb, match(paste(c(lineage_info), "ccl", sep = "_"), colnames(ccl_raw_nonNA_scaling)) + 353 * 2],
    nrow = dim(lineage_info)[1]
  ) - lineage_velocity_z_avg) / lineage_velocity_z_sd

  M_D_GD_df <- data.frame(
    M        = c(lineage_velocity_emb[1:8, ]),
    D        = c(lineage_velocity_emb[2:9, ]),
    GD       = c(lineage_velocity_emb[3:10, ]),
    embryoID = rep(emb, length(c(lineage_velocity_emb[1:8, ])))
  )
  return(M_D_GD_df)
})))

# Combine RNAi x, y, z
RNAi_combined_df_xyz_all <- rbind(
  RNAi_combined_df_x_all,
  RNAi_combined_df_y_all,
  RNAi_combined_df_z_all
)

## Env perturbation (embryos 1929 – 2039)
# x
Env_combined_df_x_all <- na.omit(do.call(rbind, lapply(1929:2039, function(emb) {
  lineage_velocity_emb <- (matrix(
    velocity[emb, match(paste(c(lineage_info), "ccl", sep = "_"), colnames(ccl_raw_nonNA_scaling))],
    nrow = dim(lineage_info)[1]
  ) - lineage_velocity_x_avg) / lineage_velocity_x_sd

  M_D_GD_df <- data.frame(
    M        = c(lineage_velocity_emb[1:8, ]),
    D        = c(lineage_velocity_emb[2:9, ]),
    GD       = c(lineage_velocity_emb[3:10, ]),
    embryoID = rep(emb, length(c(lineage_velocity_emb[1:8, ])))
  )
  return(M_D_GD_df)
})))

# y
Env_combined_df_y_all <- na.omit(do.call(rbind, lapply(1929:2039, function(emb) {
  lineage_velocity_emb <- (matrix(
    velocity[emb, match(paste(c(lineage_info), "ccl", sep = "_"), colnames(ccl_raw_nonNA_scaling)) + 353],
    nrow = dim(lineage_info)[1]
  ) - lineage_velocity_y_avg) / lineage_velocity_y_sd

  M_D_GD_df <- data.frame(
    M        = c(lineage_velocity_emb[1:8, ]),
    D        = c(lineage_velocity_emb[2:9, ]),
    GD       = c(lineage_velocity_emb[3:10, ]),
    embryoID = rep(emb, length(c(lineage_velocity_emb[1:8, ])))
  )
  return(M_D_GD_df)
})))

# z
Env_combined_df_z_all <- na.omit(do.call(rbind, lapply(1929:2039, function(emb) {
  lineage_velocity_emb <- (matrix(
    velocity[emb, match(paste(c(lineage_info), "ccl", sep = "_"), colnames(ccl_raw_nonNA_scaling)) + 353 * 2],
    nrow = dim(lineage_info)[1]
  ) - lineage_velocity_z_avg) / lineage_velocity_z_sd

  M_D_GD_df <- data.frame(
    M        = c(lineage_velocity_emb[1:8, ]),
    D        = c(lineage_velocity_emb[2:9, ]),
    GD       = c(lineage_velocity_emb[3:10, ]),
    embryoID = rep(emb, length(c(lineage_velocity_emb[1:8, ])))
  )
  return(M_D_GD_df)
})))

# Combine Env x, y, z
Env_combined_df_xyz_all <- rbind(
  Env_combined_df_x_all,
  Env_combined_df_y_all,
  Env_combined_df_z_all
)

## Combine all embryos: WT (already in combined_df_xyz), RNAi, Env
embryo_combined_df_xyz <- rbind(
  combined_df_xyz,
  RNAi_combined_df_xyz_all,
  Env_combined_df_xyz_all
)

# ========================================================================
# 4. Run phi estimation for all embryos (1 – 2039)
# ========================================================================
population_embryo_xyz_SE <- data.frame(t(sapply(1:2039, function(emb) {
  print(emb)
  sub <- embryo_combined_df_xyz[embryo_combined_df_xyz$embryoID == emb, ]
  result <- tryCatch({
    temp <- estimate_phi(sub$GD, sub$D, sub$M, phi_start = 0.618)
    c(temp$phi_hat, temp$se)
  }, error = function(e) {
    message(paste("Error for embryo", emb, ":", e$message))
    c(NA, NA)
  })
  return(result)
})))
names(population_embryo_xyz_SE) <- c("phi", "SE")
cat("Number of embryos with NA phi:", sum(is.na(population_embryo_xyz_SE$phi)), "\n")

# ========================================================================
# 5. Compute MAE and Z‑scores for each embryo (x, y, z)
#     and identify genes violating self‑similarity
# ========================================================================
coef_M_expected <- -phi_golden^3
coef_D_expected <- -phi_golden^2

M_D_GD_mae_x <- sapply(1:dim(velocity)[1], function(emb) {
  lineage_velocity_emb <- (matrix(
    velocity[emb, match(paste(c(lineage_info), "ccl", sep = "_"), colnames(ccl_raw_nonNA_scaling))],
    nrow = dim(lineage_info)[1]
  ) - lineage_velocity_x_avg) / lineage_velocity_x_sd

  M_D_GD_df <- data.frame(
    M  = c(lineage_velocity_emb[1:8, ]),
    D  = c(lineage_velocity_emb[2:9, ]),
    GD = c(lineage_velocity_emb[3:10, ])
  )

  df <- na.omit(data.frame(
    obs  = M_D_GD_df$GD,
    pred = M_D_GD_df$M * coef_M_expected + M_D_GD_df$D * coef_D_expected
  ))
  mae <- sum(abs(df[, 1] - df[, 2])) / dim(df)[1]
  return(mae)
})

M_D_GD_mae_y <- sapply(1:dim(velocity)[1], function(emb) {
  lineage_velocity_emb <- (matrix(
    velocity[emb, match(paste(c(lineage_info), "ccl", sep = "_"), colnames(ccl_raw_nonNA_scaling)) + 353],
    nrow = dim(lineage_info)[1]
  ) - lineage_velocity_y_avg) / lineage_velocity_y_sd

  M_D_GD_df <- data.frame(
    M  = c(lineage_velocity_emb[1:8, ]),
    D  = c(lineage_velocity_emb[2:9, ]),
    GD = c(lineage_velocity_emb[3:10, ])
  )

  df <- na.omit(data.frame(
    obs  = M_D_GD_df$GD,
    pred = M_D_GD_df$M * coef_M_expected + M_D_GD_df$D * coef_D_expected
  ))
  mae <- sum(abs(df[, 1] - df[, 2])) / dim(df)[1]
  return(mae)
})

M_D_GD_mae_z <- sapply(1:dim(velocity)[1], function(emb) {
  lineage_velocity_emb <- (matrix(
    velocity[emb, match(paste(c(lineage_info), "ccl", sep = "_"), colnames(ccl_raw_nonNA_scaling)) + 353 * 2],
    nrow = dim(lineage_info)[1]
  ) - lineage_velocity_z_avg) / lineage_velocity_z_sd

  M_D_GD_df <- data.frame(
    M  = c(lineage_velocity_emb[1:8, ]),
    D  = c(lineage_velocity_emb[2:9, ]),
    GD = c(lineage_velocity_emb[3:10, ])
  )

  df <- na.omit(data.frame(
    obs  = M_D_GD_df$GD,
    pred = M_D_GD_df$M * coef_M_expected + M_D_GD_df$D * coef_D_expected
  ))
  mae <- sum(abs(df[, 1] - df[, 2])) / dim(df)[1]
  return(mae)
})

# Z‑score relative to first 105 control embryos
M_D_GD_mae_x_Zscore <- sapply(1:length(gene_extend), function(i) {
  (M_D_GD_mae_x[i] - mean(M_D_GD_mae_x[1:105])) / sd(M_D_GD_mae_x[1:105])
})
M_D_GD_mae_y_Zscore <- sapply(1:length(gene_extend), function(i) {
  (M_D_GD_mae_y[i] - mean(M_D_GD_mae_y[1:105])) / sd(M_D_GD_mae_y[1:105])
})
M_D_GD_mae_z_Zscore <- sapply(1:length(gene_extend), function(i) {
  (M_D_GD_mae_z[i] - mean(M_D_GD_mae_z[1:105])) / sd(M_D_GD_mae_z[1:105])
})

# ---- Gene identification (single occurrence) ----
# x
M_D_GD_mae_x_df <- data.frame(
  gene    = gene_extend,
  Zscore  = M_D_GD_mae_x_Zscore
)
M_D_GD_mae_x_df$pvalue <- 1 - pnorm(M_D_GD_mae_x_df$Zscore)
M_D_GD_mae_x_df$padj   <- p.adjust(M_D_GD_mae_x_df$pvalue, method = "BH")
M_D_GD_mae_x_df$sig_count   <- sapply(M_D_GD_mae_x_df$gene, function(x) {
  length(which(M_D_GD_mae_x_df$gene == x & M_D_GD_mae_x_df$padj < 0.05))
})
M_D_GD_mae_x_df$total_count <- sapply(M_D_GD_mae_x_df$gene, function(x) {
  length(which(M_D_GD_mae_x_df$gene == x))
})

# y
M_D_GD_mae_y_df <- data.frame(
  gene    = gene_extend,
  Zscore  = M_D_GD_mae_y_Zscore
)
M_D_GD_mae_y_df$pvalue <- 1 - pnorm(M_D_GD_mae_y_df$Zscore)
M_D_GD_mae_y_df$padj   <- p.adjust(M_D_GD_mae_y_df$pvalue, method = "BH")
M_D_GD_mae_y_df$sig_count   <- sapply(M_D_GD_mae_y_df$gene, function(x) {
  length(which(M_D_GD_mae_y_df$gene == x & M_D_GD_mae_y_df$padj < 0.05))
})
M_D_GD_mae_y_df$total_count <- sapply(M_D_GD_mae_y_df$gene, function(x) {
  length(which(M_D_GD_mae_y_df$gene == x))
})

# z
M_D_GD_mae_z_df <- data.frame(
  gene    = gene_extend,
  Zscore  = M_D_GD_mae_z_Zscore
)
M_D_GD_mae_z_df$pvalue <- 1 - pnorm(M_D_GD_mae_z_df$Zscore)
M_D_GD_mae_z_df$padj   <- p.adjust(M_D_GD_mae_z_df$pvalue, method = "BH")
M_D_GD_mae_z_df$sig_count   <- sapply(M_D_GD_mae_z_df$gene, function(x) {
  length(which(M_D_GD_mae_z_df$gene == x & M_D_GD_mae_z_df$padj < 0.05))
})
M_D_GD_mae_z_df$total_count <- sapply(M_D_GD_mae_z_df$gene, function(x) {
  length(which(M_D_GD_mae_z_df$gene == x))
})

# Extract significant genes (padj < 0.05, sig_count > 1) in each axis
sigGene_x <- unique(M_D_GD_mae_x_df$gene[
  intersect(which(M_D_GD_mae_x_df$sig_count > 1), 106:1928)
])
sigGene_y <- unique(M_D_GD_mae_y_df$gene[
  intersect(which(M_D_GD_mae_y_df$sig_count > 1), 106:1928)
])
sigGene_z <- unique(M_D_GD_mae_z_df$gene[
  intersect(which(M_D_GD_mae_z_df$sig_count > 1), 106:1928)
])

# Genes significant in at least two axes (for GO and phi analysis)
sig_gene_phi <- tolower(names(table(c(sigGene_x, sigGene_y, sigGene_z))[
  which(table(c(sigGene_x, sigGene_y, sigGene_z)) > 1)
]))

# Intersection/differences for reporting
setdiff(sigGene_x, sigGene_y)
setdiff(sigGene_x, sigGene_z)
setdiff(sigGene_y, sigGene_x)
setdiff(sigGene_y, sigGene_z)
setdiff(sigGene_z, sigGene_x)
setdiff(sigGene_z, sigGene_y)


# ========================================================================
# 6. Combine data frames for meta‑analysis of phi
# ========================================================================
data_summary <- data.frame(
  gene          = M_D_GD_mae_x_df$gene,
  Zscore_x      = M_D_GD_mae_x_df$Zscore,
  pvalue_x      = M_D_GD_mae_x_df$pvalue,
  padj_x        = M_D_GD_mae_x_df$padj,
  sig_count_x   = M_D_GD_mae_x_df$sig_count,
  total_count_x = M_D_GD_mae_x_df$total_count,
  Zscore_y      = M_D_GD_mae_y_df$Zscore,
  pvalue_y      = M_D_GD_mae_y_df$pvalue,
  padj_y        = M_D_GD_mae_y_df$padj,
  sig_count_y   = M_D_GD_mae_y_df$sig_count,
  total_count_y = M_D_GD_mae_y_df$total_count,
  Zscore_z      = M_D_GD_mae_z_df$Zscore,
  pvalue_z      = M_D_GD_mae_z_df$pvalue,
  padj_z        = M_D_GD_mae_z_df$padj,
  sig_count_z   = M_D_GD_mae_z_df$sig_count,
  total_count_z = M_D_GD_mae_z_df$total_count,
  phi           = population_embryo_xyz_SE$phi,
  SE            = population_embryo_xyz_SE$SE
)
data_summary$embryoID <- 1:nrow(data_summary)

# ========================================================================
# 7. (Re‑use) Identify significant genes from MAE (already in sig_gene_phi)
#      No need to recompute – we already have sig_gene_phi from Section 5.
# ========================================================================
cat("Number of significant genes (from MAE):", length(sig_gene_phi), "\n")

# ========================================================================
# 8. Extract phi and SE for control and significant embryos
# ========================================================================
control_idx <- 1:105
control_phi <- data_summary$phi[control_idx]
control_se  <- data_summary$SE[control_idx]
valid_ctrl <- !is.na(control_phi) & !is.na(control_se)
control_phi <- control_phi[valid_ctrl]
control_se  <- control_se[valid_ctrl]

sig_idx <- which(tolower(data_summary$gene) %in% sig_gene_phi &
                   !is.na(data_summary$phi) & !is.na(data_summary$SE))
sig_phi   <- data_summary$phi[sig_idx]
sig_se    <- data_summary$SE[sig_idx]
sig_genes <- data_summary$gene[sig_idx]

# Convert gene names to lower case for x-axis labels
sig_genes_lower <- tolower(sig_genes)

plot_df <- data.frame(
  phi   = c(control_phi, sig_phi),
  group = c(rep("Control", length(control_phi)), sig_genes_lower)
)

# Keep genes with at least 2 embryos
keep <- names(table(plot_df$group))[table(plot_df$group) >= 2]
plot_df <- plot_df[plot_df$group %in% keep, ]
plot_df$group <- factor(plot_df$group,
                        levels = c("Control", sort(setdiff(unique(plot_df$group), "Control"))))

# ========================================================================
# 9. Meta‑analysis per group using RANDOM‑EFFECTS (I² > 30%)
# ========================================================================
groups <- levels(plot_df$group)
meta_df <- data.frame(
  group    = groups,
  phi_meta = NA,
  se_meta  = NA,
  I2       = NA,
  tau2     = NA,
  Q        = NA,
  Q_df     = NA,
  pvalue   = NA,
  padj     = NA
)

for (g in groups) {
  if (g == "Control") {
    phi_vec <- control_phi
    se_vec  <- control_se
  } else {
    idx <- which(tolower(data_summary$gene) == g &
                   !is.na(data_summary$phi) & !is.na(data_summary$SE))
    phi_vec <- data_summary$phi[idx]
    se_vec  <- data_summary$SE[idx]
  }
  if (length(phi_vec) > 0) {
    res <- meta_phi(phi_vec, se_vec)
    meta_df$phi_meta[meta_df$group == g] <- res$random$estimate
    meta_df$se_meta[meta_df$group == g]  <- res$random$se
    meta_df$I2[meta_df$group == g]       <- res$heterogeneity$I2
    meta_df$tau2[meta_df$group == g]     <- res$random$tau2
    meta_df$Q[meta_df$group == g]        <- res$heterogeneity$Q
    meta_df$Q_df[meta_df$group == g]     <- res$heterogeneity$df
  }
}

# Ensure factor order (Control first)
meta_df$group <- factor(meta_df$group, levels = levels(plot_df$group))
meta_df <- meta_df[order(meta_df$group), ]

# ========================================================================
# 10. Statistical test (compare random‑effects estimates to control)
# ========================================================================
ctrl_est <- meta_df$phi_meta[meta_df$group == "Control"]
ctrl_se  <- meta_df$se_meta[meta_df$group == "Control"]

gene_groups <- setdiff(meta_df$group, "Control")
pvals <- sapply(gene_groups, function(g) {
  phi_g <- meta_df$phi_meta[meta_df$group == g]
  se_g  <- meta_df$se_meta[meta_df$group == g]
  if (is.na(phi_g) || is.na(se_g)) return(NA)
  z <- (phi_g - ctrl_est) / sqrt(se_g^2 + ctrl_se^2)
  2 * pnorm(-abs(z))
})
padj <- p.adjust(pvals, method = "BH")
names(padj) <- gene_groups

for (g in gene_groups) {
  meta_df$pvalue[meta_df$group == g] <- pvals[g]
  meta_df$padj[meta_df$group == g]   <- padj[g]
}

sig_diff_genes <- names(padj)[padj < 0.05 & !is.na(padj)]
cat("Genes with significant phi difference (padj < 0.05):\n",
    paste(sig_diff_genes, collapse = ", "), "\n")

# ========================================================================
# 11. Barplot with legend (subtitle removed)
# ========================================================================
# Add a significance category column to meta_df
meta_df$significance <- "Not significant"
meta_df$significance[meta_df$group == "Control"] <- "Control"
meta_df$significance[meta_df$group %in% sig_diff_genes] <- "Significant"
meta_df$significance <- factor(meta_df$significance,
                               levels = c("Control", "Not significant", "Significant"))

# Define colors for the three categories
fill_colors <- c(
  "Control"         = "gray75",
  "Not significant" = "#2E86AB",   # soft blue
  "Significant"     = "#D62828"    # strong red
)

ggplot() +
  geom_bar(data = meta_df,
           aes(x = group, y = phi_meta, fill = significance),
           stat = "identity", width = 0.6, alpha = 0.8) +
  geom_point(data = plot_df,
             aes(x = group, y = phi),
             position = position_jitter(width = 0.15),
             size = 1.5, alpha = 0.6, col = "grey40") +
  geom_errorbar(data = meta_df,
                aes(x = group,
                    ymin = phi_meta - 1.96 * se_meta,
                    ymax = phi_meta + 1.96 * se_meta),
                width = 0.2, col = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0.618, lty = 2, col = "grey40") +
  labs(
    title = "Phi estimates: Control vs. significant genes (random-effects meta-analysis)",
    x     = "",
    y     = "Phi"
  ) +
  theme(
    axis.text.x      = element_text(face = "italic", angle = 45, hjust = 1, size = 9),
    legend.position  = "right",   # or "bottom"
    legend.title     = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = fill_colors)

ggsave("phi_estimate_for_individual_significant_genes_meta_analysis.pdf",width=unit(6.4,"in"),height=unit(3,"in"))

# ========================================================================
# 12. Print meta_df for reporting (includes all variables)
# ========================================================================
print(meta_df)

# ========================================================================
# 13. new version (pooled estimate)
# ========================================================================
# ========================================================================
# Alternative version: pooled phi estimation with fallback for unstable fits
# ========================================================================

# ---- 13.1. Extract raw GD, D, M for each group ----
control_embryos <- 1:105
ctrl_data <- embryo_combined_df_xyz[embryo_combined_df_xyz$embryoID %in% control_embryos, ]
ctrl_fit <- estimate_phi(ctrl_data$GD, ctrl_data$D, ctrl_data$M, phi_start = 0.618)
ctrl_phi <- ctrl_fit$phi_hat
ctrl_se  <- ctrl_fit$se

gene_list <- sig_gene_phi
pooled_df <- data.frame(
  group      = c("Control", gene_list),
  phi        = NA,
  se         = NA,
  n_embryos  = NA,
  fit_status = "success"   # "success" or "failed"
)

pooled_df$phi[pooled_df$group == "Control"] <- ctrl_phi
pooled_df$se[pooled_df$group == "Control"]  <- ctrl_se
pooled_df$n_embryos[pooled_df$group == "Control"] <- length(control_embryos)

for (g in gene_list) {
  emb_ids <- data_summary$embryoID[tolower(data_summary$gene) == g]
  if (length(emb_ids) == 0) next
  sub_data <- embryo_combined_df_xyz[embryo_combined_df_xyz$embryoID %in% emb_ids, ]
  fit <- tryCatch({
    estimate_phi(sub_data$GD, sub_data$D, sub_data$M, phi_start = 0.618)
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    # Check for numerical instability
    if (is.na(fit$phi_hat) || is.na(fit$se)) {
      pooled_df$fit_status[pooled_df$group == g] <- "failed"
      next
    }
    if (abs(fit$phi_hat) < 1e-6 && fit$se > 1e3) {
      pooled_df$fit_status[pooled_df$group == g] <- "failed"
      warning(paste("Unstable fit for", g, "- will show as dashed bar"))
      next
    }
    pooled_df$phi[pooled_df$group == g] <- fit$phi_hat
    pooled_df$se[pooled_df$group == g]  <- fit$se
    pooled_df$n_embryos[pooled_df$group == g] <- length(emb_ids)
  } else {
    pooled_df$fit_status[pooled_df$group == g] <- "failed"
  }
}

# ---- Set factor levels: Control first, then genes alphabetically ----
pooled_df$group <- factor(pooled_df$group,
                          levels = c("Control", sort(setdiff(pooled_df$group, "Control"))))
pooled_df <- pooled_df[order(pooled_df$group), ]

# ---- 13.2. Statistical test (only for successful fits) ----
success_df <- pooled_df[pooled_df$fit_status == "success", ]
ctrl_phi <- success_df$phi[success_df$group == "Control"]
ctrl_se  <- success_df$se[success_df$group == "Control"]

gene_groups <- setdiff(success_df$group, "Control")
pvals <- sapply(gene_groups, function(g) {
  phi_g <- success_df$phi[success_df$group == g]
  se_g  <- success_df$se[success_df$group == g]
  if (is.na(phi_g) || is.na(se_g)) return(NA)
  z <- (phi_g - ctrl_phi) / sqrt(se_g^2 + ctrl_se^2)
  2 * pnorm(-abs(z))
})
padj <- p.adjust(pvals, method = "BH")
names(padj) <- gene_groups

pooled_df$pvalue <- NA
pooled_df$padj   <- NA
for (g in gene_groups) {
  pooled_df$pvalue[pooled_df$group == g] <- pvals[g]
  pooled_df$padj[pooled_df$group == g]   <- padj[g]
}

sig_diff_genes <- names(padj)[padj < 0.05 & !is.na(padj)]
cat("Genes with significant phi difference (pooled, padj < 0.05):\n",
    paste(sig_diff_genes, collapse = ", "), "\n")

# ---- 13.3. Define plotting categories ----
pooled_df$category <- "Not significant"
pooled_df$category[pooled_df$group == "Control"] <- "Control"
pooled_df$category[pooled_df$group %in% sig_diff_genes] <- "Significant"
pooled_df$category[pooled_df$fit_status == "failed"] <- "Unstable fit"
pooled_df$category <- factor(pooled_df$category,
                             levels = c("Control", "Not significant", "Significant", "Unstable fit"))

# For unstable fits, set phi to NA so no bar is drawn (we'll add a dashed bar separately)
pooled_df$phi_plot <- pooled_df$phi
pooled_df$phi_plot[pooled_df$fit_status == "failed"] <- NA

# ---- 13.4. Prepare individual phi values for points ----
population_embryo_xyz_SE$embryoID <- 1:nrow(population_embryo_xyz_SE)

embryo_groups <- data_summary[, c("embryoID", "gene")]
embryo_groups$group <- ifelse(embryo_groups$embryoID %in% 1:105, "Control",
                              tolower(embryo_groups$gene))

valid_groups <- pooled_df$group
embryo_groups <- embryo_groups[tolower(embryo_groups$group) %in% valid_groups, ]
embryo_phi <- merge(embryo_groups, population_embryo_xyz_SE, by = "embryoID")
embryo_phi <- embryo_phi[!is.na(embryo_phi$phi), ]
embryo_phi$group <- factor(embryo_phi$group, levels = levels(pooled_df$group))

# ---- 13.5. Barplot with dashed bars for unstable fits ----
library(ggplot2)

# Colors for the four categories (unstable uses transparent fill)
fill_colors <- c(
  "Control"         = "gray75",
  "Not significant" = "#2E86AB",
  "Significant"     = "#D62828",
  "Unstable fit"    = "transparent"   # no fill
)

# Separate unstable groups for dashed border
unstable_df <- pooled_df[pooled_df$fit_status == "failed", ]

ggplot() +
  # 1. Solid bars for all groups (fill mapped to category)
  geom_bar(data = pooled_df,
           aes(x = group, y = phi_plot, fill = category),
           stat = "identity", width = 0.6, alpha = 0.8) +
  # 2. Dashed border bars for unstable fits (drawn on top, no fill)
  geom_bar(data = unstable_df,
           aes(x = group, y = 0),   # height 0 to keep the baseline
           stat = "identity", width = 0.6,
           fill = NA, color = "black", linetype = "dashed", linewidth = 0.5) +
  # 3. Individual points
  geom_point(data = embryo_phi,
             aes(x = group, y = phi),
             position = position_jitter(width = 0.15),
             pch=16, size = 0.5, stroke=0.2,alpha = 0.6, col = "grey40") +
  # 4. Error bars for successful fits
  geom_errorbar(data = pooled_df[pooled_df$fit_status == "success", ],
                aes(x = group,
                    ymin = phi - 1.96 * se,
                    ymax = phi + 1.96 * se),
                width = 0.3, col = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0.618, lty = 2, col = "grey40") +
  labs(
    title = "φ estimates: Control vs. significant genes (pooled phi estimation)",
    x     = "",
    y     = "φ"
  ) +
  theme(
    axis.text.x      = element_text(face = "italic", angle = 45, hjust = 1, size = 9),
    legend.position  = "right",
    legend.title     = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(
    values = fill_colors,
    breaks = c("Control", "Not significant", "Significant", "Unstable fit"),
    labels = c("Control", "Not significant", "Significant", "Unstable fit")
  ) +
  # Override legend key for "Unstable fit" to show dashed border and transparent fill
  guides(fill = guide_legend(
    override.aes = list(
      linetype = c("solid", "solid", "solid", "dashed"),
      color = c("gray75", "#2E86AB", "#D62828", "black"),
      fill = c("gray75", "#2E86AB", "#D62828", "transparent")
    )
  ))

ggsave("phi_estimate_for_individual_significant_genes.pdf",width=unit(6.4,"in"),height=unit(3,"in"))

hist(-log10(pooled_df$padj),xlab="-log10(padj)",1000)

sort(pooled_df$padj[pooled_df$category=="Significant"])

# ---- 13.6. Print pooled results ----
print(pooled_df)



# ========================================================================
# 14. Comparison between lethal and hatched embryos (MAE Z‑scores)
# ========================================================================
dataTemp<-data.frame(
	mae=c(M_D_GD_mae_x,M_D_GD_mae_y,M_D_GD_mae_z),
	type=factor(rep(c(rep("WT",105),rep("RNAi",1928-105),rep("Env",2039-1928)),3),level=c("WT","RNAi","Env")),
	coordinate=c(rep("x",length(M_D_GD_mae_x)),rep("y",length(M_D_GD_mae_y)),rep("z",length(M_D_GD_mae_z))),
	lethal=factor(rep(embryo_lethal,3),
	              levels=c(0,1),
	              labels=c("Hatched","Lethal"))
)

ggplot(dataTemp,aes(x=type,y=mae,col=coordinate))+
geom_boxplot()

ggplot(dataTemp,aes(x=coordinate,y=mae,col=lethal))+
geom_boxplot(size=0.35,width=0.75,lwd=0.2,outlier.size=0.1)+
labs(x="coordinate",y="Deviation level",title="")+
theme(legend.position=c(0.8,0.8))+
scale_y_continuous(limit=c(0,3))

ggsave("compare_mae_lethal_and_normal.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


# Wilcoxon tests (one‑sided: lethal > normal)
wilcox.test(dataTemp$mae[dataTemp$lethal == 1 & dataTemp$coordinate == "x"],
            dataTemp$mae[dataTemp$lethal == 0 & dataTemp$coordinate == "x"],
            alternative = "greater")
wilcox.test(dataTemp$mae[dataTemp$lethal == 1 & dataTemp$coordinate == "y"],
            dataTemp$mae[dataTemp$lethal == 0 & dataTemp$coordinate == "y"],
            alternative = "greater")
wilcox.test(dataTemp$mae[dataTemp$lethal == 1 & dataTemp$coordinate == "z"],
            dataTemp$mae[dataTemp$lethal == 0 & dataTemp$coordinate == "z"],
            alternative = "greater")

dataTemp <- data.frame(
  mae = c(M_D_GD_mae_x_Zscore, M_D_GD_mae_y_Zscore, M_D_GD_mae_z_Zscore),
  type = factor(
    rep(c(rep("WT", 105), rep("RNAi", 1928 - 105), rep("Env", 2039 - 1928)), 3),
    levels = c("WT", "RNAi", "Env")
  ),
  coordinate = c(
    rep("x", length(M_D_GD_mae_x_Zscore)),
    rep("y", length(M_D_GD_mae_y_Zscore)),
    rep("z", length(M_D_GD_mae_z_Zscore))
  ),
  lethal = rep(embryo_lethal, 3)
)

ggplot(dataTemp, aes(x = coordinate, y = mae, col = lethal)) +
  geom_boxplot(width = 0.5, lwd = 0.1, outlier.size = 0.1, outlier.stroke = 0.2) +
  labs(x = "coordinate", y = "Z-score of MAE", title = "") +
  theme(legend.position = c(0.6, 0.8)) +
  scale_y_continuous(breaks = c(0, 20, 40)) +
  theme(
    legend.background      = element_blank(),
    legend.frame           = element_blank(),
    legend.box.background  = element_blank(),
    legend.spacing.x       = unit(0.1, 'cm'),
    legend.spacing.y       = unit(0.1, 'cm'),
    legend.key.width       = unit(0.3, 'cm'),
    legend.key.height      = unit(0.3, 'cm')
  )

#ggsave("compare_mae_lethal_and_normal.pdf", width = unit(1.4, "in"), height = unit(1.5, "in"))

# ========================================================================
# 15. Logistic regression and ROC for lethality prediction
# ========================================================================
library(pROC)
library(caret)
library(tidyverse)

dataTemp <- data.frame(
  x = M_D_GD_mae_x_Zscore,
  y = M_D_GD_mae_y_Zscore,
  z = M_D_GD_mae_z_Zscore,
  lethal = embryo_lethal
)

# Split data (70% training)
set.seed(123)
train_index <- sample(1:dim(dataTemp)[1], round(dim(dataTemp)[1] * 0.7))
train_data <- dataTemp[train_index, ]
test_data  <- dataTemp[-train_index, ]

# Fit logistic regression
logit_model <- glm(as.numeric(lethal) ~ x + y + z,
                   data = train_data,
                   family = binomial(link = "logit"))
summary(logit_model)

# Predict on test set
test_predictions <- predict(logit_model, newdata = test_data, type = "response")
predicted_class <- ifelse(test_predictions > 0.5, 1, 0)
predicted_class <- factor(predicted_class, levels = c(0, 1))

# ROC and AUC
roc_result <- roc(response = as.numeric(test_data$lethal),
                  predictor = test_predictions)
cat("AUC:", auc(roc_result), "\n")

# Basic ROC plot
plot(roc_result,
     main = paste("ROC Curve (AUC =", round(auc(roc_result), 3), ")"),
     col = "blue", lwd = 2,
     print.auc = TRUE,
     print.auc.x = 0.6, print.auc.y = 0.4)
abline(a = 0, b = 1, lty = 2, col = "gray")

# Enhanced ROC with ggplot
roc_df <- data.frame(
  sensitivity = roc_result$sensitivities,
  specificity = roc_result$specificities,
  thresholds  = roc_result$thresholds
)
roc_df$one_minus_specificity <- 1 - roc_df$specificity

ggplot(roc_df, aes(x = one_minus_specificity, y = sensitivity)) +
  geom_line(color = "blue", size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_polygon(aes(x = one_minus_specificity, y = sensitivity),
               fill = "skyblue", alpha = 0.3) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "False Positive Rate", y = "True Positive Rate",
       title = "Lethality prediction") +
  annotate("text", x = 0.7, y = 0.3,
           label = paste("AUC =", round(auc(roc_result), 3)),
           size = 2, color = "darkred", fontface = "bold")

ggsave("lethal_prediction_by_phi_model.pdf", width = unit(1.4, "in"), height = unit(1.5, "in"))

# ========================================================================
# 16. GO enrichment analysis (using sig_gene_phi)
# ========================================================================
library(org.Ce.eg.db)
library(clusterProfiler)

# Map symbols to ENTREZID (first round)
gene_mapped <- mapIds(
  org.Ce.eg.db,
  keys    = sig_gene_phi,
  keytype = "SYMBOL",
  column  = "ENTREZID"
)

# Replace missing with corrected names (manually)
which(is.na(gene_mapped))
sig_gene_phi[which(is.na(gene_mapped))] <- c("hrpk-1", "ints-1", "ints-4")

# Second round mapping
gene_mapped <- mapIds(
  org.Ce.eg.db,
  keys    = sig_gene_phi,
  keytype = "SYMBOL",
  column  = "ENTREZID"
)

# GO enrichment
ego_for_sig_gene_phi <- enrichGO(
  gene          = gene_mapped,
  OrgDb         = org.Ce.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

df_go <- data.frame(ego_for_sig_gene_phi)

# Helper: compute GeneRatio as numeric
calc_ratio <- function(x) {
  v <- as.numeric(unlist(strsplit(x, "/")))
  v[1] / v[2]
}

df_go$Description <- factor(df_go$Description, levels = rev(df_go$Description))
df_go$GeneRatio   <- sapply(df_go$GeneRatio, calc_ratio)

# Dot plot
ggplot(df_go, aes(x = Count, y = Description, col = p.adjust, size = Count)) +
  geom_point() +
  scale_colour_gradient(low = "#F8766D", high = "#C77CFF") +
  labs(x = "Gene count", y = "", title = "") +
  theme(
    axis.text.y          = element_text(angle = 0, hjust = 1, vjust = 0.5, size = rel(10 / 11)),
    legend.title         = element_text(size = rel(8 / 11)),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major.y   = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.y   = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.major.x   = element_blank(),
    panel.grid.minor.x   = element_line(color = "grey90", linewidth = 0.2),
    legend.position      = "bottom"
  ) +
  scale_x_continuous(breaks = c(2, 8)) +
  scale_size(range = c(0.5, 2))

ggsave("E:/canalizationENV/figs/sig_gene_phi_GO_dotplot_all.pdf",
       width = unit(3.6, "in"), height = unit(5, "in"))

# ========================================================================
# 17. Gene–GO matrix plot
# ========================================================================
gene_GO_sparse <- data.frame(
  GO = unlist(lapply(1:dim(df_go)[1], function(i) {
    rep(df_go$Description[i], length(unlist(strsplit(df_go$geneID[i], "/"))))
  })),
  Gene = tolower(unlist(lapply(1:dim(df_go)[1], function(i) {
    unlist(strsplit(df_go$geneID[i], "/"))
  })))
)
gene_GO_sparse$value <- rep(1, dim(gene_GO_sparse)[1])

gene_GO_sparse$Gene <- factor(gene_GO_sparse$Gene, levels = sig_gene_phi)
gene_GO_sparse$GO   <- factor(gene_GO_sparse$GO, levels = rev(df_go$Description))

# Expand to all combinations
gene_GO_full <- merge(
  expand.grid(
    Gene = levels(gene_GO_sparse$Gene),
    GO   = levels(gene_GO_sparse$GO)
  ),
  gene_GO_sparse,
  by = c("Gene", "GO"),
  all.x = TRUE
)


# Plot
ggplot(gene_GO_full, aes(x = Gene, y = GO, color = value)) +
  geom_point(show.legend = FALSE, na.rm = FALSE) +
  scale_color_gradient(low = "#4CAF50", high = "#4CAF50", na.value = NA) +
  scale_x_discrete(position = "bottom") +
  labs(x = "Gene") +
  theme(
    axis.text.x        = element_text(face = "italic", angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y        = element_text(angle = 0, hjust = 1, vjust = 0.5),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.y = element_blank()
  )

ggsave("E:/canalizationENV/figs/sig_gene_phi_GO_relationship.pdf",
       width = unit(6.3, "in"), height = unit(4.52, "in"))



############################################################################
#  END
############################################################################















############################################################################
#  old 
############################################################################



### linear regression 
lm_M_D_GD_coef_x<-t(sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-na.omit(data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	))
	
	M_D_GD_model<-lm(GD~D+M+0, data=M_D_GD_df)
	return(c(coefficients(summary(M_D_GD_model))[1,1],coefficients(summary(M_D_GD_model))[2,1]))

}))

plot(lm_M_D_GD_coef_x)

summary(lm(lm_M_D_GD_coef_x[,2]~lm_M_D_GD_coef_x[,1]+0))

apply(lm_M_D_GD_coef_x,2,mean)
apply(lm_M_D_GD_coef_x,2,median)


apply(lm_M_D_GD_coef_x[1:105,],2,mean)
apply(lm_M_D_GD_coef_x[1:105,],2,mean)[2]/apply(lm_M_D_GD_coef_x[1:105,],2,mean)[1]

hist(lm_M_D_GD_coef_x[1:105,2]/lm_M_D_GD_coef_x[1:105,1],100)

M_D_GD<-na.omit(data.frame(
	M=mo_sis_pair_ccl[,1],
	D=mo_sis_pair_ccl[,2],
	GD=mo_sis_pair_ccl[match(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,1]),2]
))

# 
plot(combined_result_z[[3]],fitted(lm(GD~D+M, data=combined_df_z)))
cor(combined_result_z[[3]],fitted(lm(GD~D+M, data=combined_df_z)))^2
abline(0,1)

plot(combined_result_x[[3]],fitted(lm(GD~D+M, data=combined_df_x)))
cor(combined_result_x[[3]],fitted(lm(GD~D+M, data=combined_df_x)))^2
abline(0,1)

plot(combined_result_y[[3]],fitted(lm(GD~D+M, data=combined_df_y)))
cor(combined_result_y[[3]],fitted(lm(GD~D+M, data=combined_df_y)))^2
abline(0,1)


### gene identification 

## for MAE
golden_phi = (sqrt(5)-1)/2
coef_D_expected=-golden_phi^2
coef_M_expected=-golden_phi^3

M_D_GD_mae_x<-sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	
	df<-na.omit(data.frame(obs=M_D_GD_df$GD,pred=M_D_GD_df$M*coef_M_expected+M_D_GD_df$D*coef_D_expected))
	mae=sum(abs(df[,1]-df[,2]))/dim(df)[1]
	return(mae)	

})


M_D_GD_mae_y<-sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	
	df<-na.omit(data.frame(obs=M_D_GD_df$GD,pred=M_D_GD_df$M*coef_M_expected+M_D_GD_df$D*coef_D_expected))
	mae=sum(abs(df[,1]-df[,2]))/dim(df)[1]
	return(mae)	

})

M_D_GD_mae_z<-sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	
	df<-na.omit(data.frame(obs=M_D_GD_df$GD,pred=M_D_GD_df$M*coef_M_expected+M_D_GD_df$D*coef_D_expected))
	mae=sum(abs(df[,1]-df[,2]))/dim(df)[1]
	return(mae)	

})

cor(data.frame(M_D_GD_mae_x,M_D_GD_mae_y,M_D_GD_mae_z))

### prediction 

dataTemp<-data.frame(
	mae=c(M_D_GD_mae_x,M_D_GD_mae_y,M_D_GD_mae_z),
	type=factor(rep(c(rep("WT",105),rep("RNAi",1928-105),rep("Env",2039-1928)),3),level=c("WT","RNAi","Env")),
	coordinate=c(rep("x",length(M_D_GD_mae_x)),rep("y",length(M_D_GD_mae_y)),rep("z",length(M_D_GD_mae_z))),
	lethal=rep(embryo_lethal,3)
)

ggplot(dataTemp,aes(x=type,y=mae,col=coordinate))+
geom_boxplot()

ggplot(dataTemp,aes(x=coordinate,y=mae,col=lethal))+
geom_boxplot(size=0.35,width=0.75,lwd=0.2,outlier.size=0.1)+
labs(x="coordinate",y="Deviation level from phi-model",title="")+
theme(legend.position=c(0.8,0.8))+
scale_y_continuous(limit=c(0,3))

ggsave("compare_mae_lethal_and_normal.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


### prediction 
dataTemp<-data.frame(
	x=M_D_GD_mae_x,
	y=M_D_GD_mae_y,
	z=M_D_GD_mae_z,
	lethal=embryo_lethal
)


# Load required libraries
library(pROC)      # For ROC/AUC calculations
library(ggplot2)   # For enhanced plotting
library(caret)     # For data splitting (optional but recommended)
library(tidyverse) # For data manipulation

# ==================== STEP 1: Prepare Data ====================
# Use mtcars dataset, convert am to factor for binary classification

# Split data into training and testing sets (80/20 split)
set.seed(123)  # For reproducibility
train_index <- sample(1:dim(dataTemp)[1],400)
train_data <- dataTemp[train_index, ]
test_data <- dataTemp[-train_index, ]

# ==================== STEP 2: Fit Logistic Regression Model ====================
# Model predicting transmission type (am) based on horsepower (hp) and weight (wt)
logit_model <- glm(as.numeric(lethal) ~ x+y+z, 
                   data = train_data, 
                   family = binomial(link = "logit"))

# View model summary
summary(logit_model)

# Interpretation of coefficients:
# - Coefficients are in log-odds
# - Positive coefficient = increases probability of am=1 (manual)
# - Negative coefficient = decreases probability of am=1

# ==================== STEP 3: Make Predictions ====================
# Get predicted probabilities for the test set
test_predictions <- predict(logit_model, 
                            newdata = test_data, 
                            type = "response")

# Convert to factor with 0.5 cutoff (standard threshold)
predicted_class <- ifelse(test_predictions > 0.5, 1, 0)
predicted_class <- factor(predicted_class, levels = c(0, 1))

# ==================== STEP 4: Calculate and Plot ROC/AUC ====================
# Calculate ROC curve
roc_result <- roc(response = as.numeric(test_data$lethal), 
                  predictor = test_predictions)

# Print AUC
cat("AUC:", auc(roc_result), "\n")

# Basic ROC plot
plot(roc_result, 
     main = paste("ROC Curve (AUC =", round(auc(roc_result), 3), ")"),
     col = "blue", 
     lwd = 2,
     print.auc = TRUE,
     print.auc.x = 0.6, 
     print.auc.y = 0.4)

# Add diagonal reference line
abline(a = 0, b = 1, lty = 2, col = "gray")

# ==================== OPTIONAL: Enhanced ROC Plot with ggplot2 ====================
# Create dataframe for ROC curve
roc_df <- data.frame(
  sensitivity = roc_result$sensitivities,
  specificity = roc_result$specificities,
  thresholds = roc_result$thresholds
)

# Calculate 1 - specificity
roc_df$one_minus_specificity <- 1 - roc_df$specificity

# Create enhanced ROC plot
ggplot(roc_df, aes(x = one_minus_specificity, y = sensitivity)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_polygon(aes(x = one_minus_specificity, y = sensitivity), 
               fill = "skyblue", alpha = 0.3) +
  labs(
    title = paste("ROC Curve for Logistic Regression Model"),
    subtitle = paste("AUC =", round(auc(roc_result), 3)),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  annotate("text", x = 0.7, y = 0.3, 
           label = paste("AUC =", round(auc(roc_result), 3)),
           size = 5, color = "darkred", fontface = "bold") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("lethal_prediction_by_phi_model.pdf",width=unit(2,"in"),height=unit(2.1,"in"))

# ==================== STEP 5: Model Evaluation Metrics ====================
# Confusion matrix
conf_matrix <- confusionMatrix(data = predicted_class, 
                               reference = test_data$am, 
                               positive = "1")
print(conf_matrix)

# gene identification 

#x
M_D_GD_mae_Zscore<-sapply(1:length(gene_extend),function(i){
	return((M_D_GD_mae_x[i]-mean(M_D_GD_mae_x[1:105]))/sd(M_D_GD_mae_x[1:105]))
})

M_D_GD_mae_z_df<-data.frame(gene=gene_extend,z=M_D_GD_mae_Zscore)
M_D_GD_mae_z_df$adjp<-p.adjust(1-pnorm(M_D_GD_mae_z_df$z),"BH")
M_D_GD_mae_z_df$count<-sapply(M_D_GD_mae_z_df$gene,function(x){length(which(M_D_GD_mae_z_df$gene==x&M_D_GD_mae_z_df$adjp<0.05))})

M_D_GD_mae_z_df_coX<-M_D_GD_mae_z_df
#y
M_D_GD_mae_Zscore<-sapply(1:length(gene_extend),function(i){
	return((M_D_GD_mae_y[i]-mean(M_D_GD_mae_y[1:105]))/sd(M_D_GD_mae_y[1:105]))
})

M_D_GD_mae_z_df<-data.frame(gene=gene_extend,z=M_D_GD_mae_Zscore)
M_D_GD_mae_z_df$adjp<-p.adjust(1-pnorm(M_D_GD_mae_z_df$z),"BH")
M_D_GD_mae_z_df$count<-sapply(M_D_GD_mae_z_df$gene,function(x){length(which(M_D_GD_mae_z_df$gene==x&M_D_GD_mae_z_df$adjp<0.05))})

M_D_GD_mae_z_df_coY<-M_D_GD_mae_z_df
#z
M_D_GD_mae_Zscore<-sapply(1:length(gene_extend),function(i){
	return((M_D_GD_mae_z[i]-mean(M_D_GD_mae_z[1:105]))/sd(M_D_GD_mae_z[1:105]))
})

M_D_GD_mae_z_df<-data.frame(gene=gene_extend,z=M_D_GD_mae_Zscore)
M_D_GD_mae_z_df$adjp<-p.adjust(1-pnorm(M_D_GD_mae_z_df$z),"BH")
M_D_GD_mae_z_df$count<-sapply(M_D_GD_mae_z_df$gene,function(x){length(which(M_D_GD_mae_z_df$gene==x&M_D_GD_mae_z_df$adjp<0.05))})

M_D_GD_mae_z_df_coZ<-M_D_GD_mae_z_df

sigTemp<-M_D_GD_mae_z_df_coX$gene[M_D_GD_mae_z_df_coX$adjp<0.05]
sigGene_temp<-setdiff(names(which(table(sigTemp)>1)),"ctr")
sigGene_temp_x<-tolower(sigGene_temp)

sigTemp<-M_D_GD_mae_z_df_coY$gene[M_D_GD_mae_z_df_coY$adjp<0.05]
sigGene_temp<-setdiff(names(which(table(sigTemp)>1)),"ctr")
sigGene_temp_y<-tolower(sigGene_temp)

sigTemp<-M_D_GD_mae_z_df_coZ$gene[M_D_GD_mae_z_df_coZ$adjp<0.05]
sigGene_temp<-setdiff(names(which(table(sigTemp)>1)),"ctr")
sigGene_temp_z<-tolower(sigGene_temp)

setdiff(sigGene_temp_x,sigGene_temp_y)
setdiff(sigGene_temp_x,sigGene_temp_z)
setdiff(sigGene_temp_y,sigGene_temp_x)
setdiff(sigGene_temp_y,sigGene_temp_z)
setdiff(sigGene_temp_z,sigGene_temp_x)
setdiff(sigGene_temp_z,sigGene_temp_y)

sig_gene_phi<-names(table(c(sigGene_temp_x,sigGene_temp_y,sigGene_temp_z)))[which(table(c(sigGene_temp_x,sigGene_temp_y,sigGene_temp_z))>1)]



#write.csv(M_D_GD_mae_z_df,"M_D_GD_mae_z_df.csv")

gene_mapped_temp <- mapIds(
  org.Ce.eg.db,
  keys = sig_gene_phi,
  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
  column = "ENTREZID"
)



sig_gene_phi[which(is.na(gene_mapped_temp))]<-c(
"C48B6.2",
"eif-3.E",
"hrpk-1",
"hs25m",
"ints-1",
"sa600",
"ints-4",
"slc-30A5",
"Y106G6E.4",
"edc-4",
"Y47D9A.1",
"cox-4",
"24c",
"B0041.5",
"D2030.2",
"kbrl-1",
"fe1mm",
"fe2mm",    
"fe4mm",
"K12C11.1",
"tcab-1"
)





gene_mapped_temp <- mapIds(
  org.Ce.eg.db,
  keys = sig_gene_phi,
  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
  column = "ENTREZID"
)

ego_temp <- enrichGO(
  gene = gene_mapped_temp,
  OrgDb = org.Ce.eg.db,
  ont = "ALL",  # or BP, CC, MF
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

data.frame(ego_temp)

barplot(ego_temp, showCategory = 10, label_format = 100)

ggsave("GO_phi.pdf",width=unit(6,"in"),height=unit(2.5,"in"))



################### explore displacement noise and CCL noise 

### displacement 

#x
lineage_displacement_x_avg<-matrix(apply(displacement[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_displacement_x_sd<-matrix(apply(displacement[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,sd,na.rm=T),nrow=dim(lineage_info)[1])

combined_displacement_x<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_displacement_x_avg)/lineage_displacement_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_displacement_emb[1:8,]),
		D=c(lineage_displacement_emb[2:9,]),
		GD=c(lineage_displacement_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_result_displacement_x<-estimate_phi(combined_displacement_x$GD, combined_displacement_x$D, combined_displacement_x$M, phi_start = 0.618)
combined_result_displacement_x[[1]]
combined_result_displacement_x[[2]]

#y
lineage_displacement_y_avg<-matrix(apply(displacement[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_displacement_y_sd<-matrix(apply(displacement[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],2,sd,na.rm=T),nrow=dim(lineage_info)[1])

combined_displacement_y<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_displacement_y_avg)/lineage_displacement_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_displacement_emb[1:8,]),
		D=c(lineage_displacement_emb[2:9,]),
		GD=c(lineage_displacement_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_result_displacement_y<-estimate_phi(combined_displacement_y$GD, combined_displacement_y$D, combined_displacement_y$M, phi_start = 0.618)
combined_result_displacement_y[[1]]
combined_result_displacement_y[[2]]

#z
lineage_displacement_z_avg<-matrix(apply(displacement[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_displacement_z_sd<-matrix(apply(displacement[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],2,sd,na.rm=T),nrow=dim(lineage_info)[1])

combined_displacement_z<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_displacement_z_avg)/lineage_displacement_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_displacement_emb[1:8,]),
		D=c(lineage_displacement_emb[2:9,]),
		GD=c(lineage_displacement_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_result_displacement_z<-estimate_phi(combined_displacement_z$GD, combined_displacement_z$D, combined_displacement_z$M, phi_start = 0.618)
combined_result_displacement_z[[1]]
combined_result_displacement_z[[2]]


### ccl
plot(matrix(apply(tim[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1]),
	matrix(apply(tim[106:2039,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1]))

lineage_tim_avg<-matrix(apply(tim[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_tim_sd<-matrix(apply(tim[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,sd,na.rm=T),nrow=dim(lineage_info)[1])

combined_tim<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_tim_emb<-(matrix(tim[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_tim_avg)#/lineage_tim_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_tim_emb[1:8,]),
		D=c(lineage_tim_emb[2:9,]),
		GD=c(lineage_tim_emb[3:10,])
	)
	return(M_D_GD_df)
})))

summary(lm(GD~D+M+0,data=combined_tim))

combined_result_tim<-estimate_phi(combined_tim$GD, combined_tim$D, combined_tim$M, phi_start = 0.618)
combined_result_tim[[1]]
combined_result_tim[[2]]


lineage_tim_emb<-(matrix(tim[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_tim_avg)/lineage_tim_sd
M_D_GD_df<-data.frame(
	M=c(lineage_tim_emb[1:8,]),
	D=c(lineage_tim_emb[2:9,]),
	GD=c(lineage_tim_emb[3:10,])
)
summary(lm(GD~D+M+0,data=M_D_GD_df))

### timing 

lineage_cdt_avg<-matrix(apply(cdt_raw_nonNA_scaling[1:105,match(paste(c(lineage_info),"cdt",sep="_"),colnames(cdt_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_cdt_sd<-matrix(apply(cdt_raw_nonNA_scaling[1:105,match(paste(c(lineage_info),"cdt",sep="_"),colnames(cdt_raw_nonNA_scaling))],2,sd,na.rm=T),nrow=dim(lineage_info)[1])

combined_cdt<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_cdt_emb<-(matrix(cdt_raw_nonNA_scaling[emb,match(paste(c(lineage_info),"cdt",sep="_"),colnames(cdt_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_cdt_avg)/lineage_cdt_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_cdt_emb[1:8,]),
		D=c(lineage_cdt_emb[2:9,]),
		GD=c(lineage_cdt_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_result_cdt<-estimate_phi(combined_cdt$GD, combined_cdt$D, combined_cdt$M, phi_start = 0.618)
combined_result_cdt[[1]]
combined_result_cdt[[2]]


### early cell predict later cell
# start_pos 
lineage_start_pos_x_avg<-matrix(apply(start_pos[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_start_pos_x_sd<-matrix(apply(start_pos[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,sd,na.rm=T),nrow=dim(lineage_info)[1])

combined_start_pos_x<-do.call(rbind,lapply(1:2039,function(emb){

	lineage_start_pos_emb<-(matrix(start_pos[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_start_pos_x_avg)/lineage_start_pos_x_sd
	flat_cell_noise<-c(t(lineage_start_pos_emb))
	return(flat_cell_noise)
}))
generation<-c(t(matrix(rep(1:10,361),nrow=10)))[which(!is.na(apply(combined_start_pos_x,2,mean)))]
combined_start_pos_x_refine<-combined_start_pos_x[,which(!is.na(apply(combined_start_pos_x,2,mean)))]
combined_start_pos_x_refine_g9<-combined_start_pos_x_refine[,generation==9]
combined_start_pos_x_refine_g5to8<-combined_start_pos_x_refine[,generation!=9]

library(glmnet)
predR<-c()
for(i in 1:dim(combined_start_pos_x_refine_g9)[2]){
	cvglm<-cv.glmnet(combined_start_pos_x_refine_g5to8,combined_start_pos_x_refine_g9[,i],nlambda=200,nfolds=10,alpha=1,keep=TRUE)

	predR[i]<-cor(combined_start_pos_x_refine_g9[,i],cvglm$fit.preval[,which(cvglm$lambda==cvglm$lambda.min)])
}

quantile(R_MDD_ccl_within)
hist(apply(ccl_scaled,1,sd),100)
plot(apply(ccl_scaled[,na.omit(match(paste(unique(c(lineage_info[1:7,])),"ccl",sep="_"),colnames(ccl_scaled)))],1,sd),apply(ccl_scaled[,na.omit(match(paste(unique(c(lineage_info[8:10,])),"ccl",sep="_"),colnames(ccl_scaled)))],1,sd))
abline(0,1)
cor(apply(ccl_scaled[,na.omit(match(paste(unique(c(lineage_info[1:7,])),"ccl",sep="_"),colnames(ccl_scaled)))],1,sd),apply(ccl_scaled[,na.omit(match(paste(unique(c(lineage_info[8:10,])),"ccl",sep="_"),colnames(ccl_scaled)))],1,sd))














# x
M_D_GD_mae_pvalue<-sapply(unique(gene_extend),function(x){
	geneInd<-which(gene_extend==x)
	if(length(geneInd)>1){
		return(t.test(M_D_GD_mae_x[geneInd],M_D_GD_mae_x[1:105],alternative="greater")$p.value)
	}else{
		return(NA)
	}
})

M_D_GD_mae_padj<-p.adjust(M_D_GD_mae_pvalue,"BH")
M_D_GD_mae_df_x<-data.frame(gene=unique(gene_extend),padj=M_D_GD_mae_padj)
M_D_GD_mae_df_x[which(M_D_GD_mae_df_x$padj<0.05),]
#y
M_D_GD_mae_pvalue<-sapply(unique(gene_extend),function(x){
	geneInd<-which(gene_extend==x)
	if(length(geneInd)>1){
		return(t.test(M_D_GD_mae_y[geneInd],M_D_GD_mae_y[1:105],alternative="greater")$p.value)
	}else{
		return(NA)
	}
})

M_D_GD_mae_padj<-p.adjust(M_D_GD_mae_pvalue,"BH")
M_D_GD_mae_df_y<-data.frame(gene=unique(gene_extend),padj=M_D_GD_mae_padj)
M_D_GD_mae_df_y[which(M_D_GD_mae_df_y$padj<0.05),]
#y
M_D_GD_mae_pvalue<-sapply(unique(gene_extend),function(x){
	geneInd<-which(gene_extend==x)
	if(length(geneInd)>1){
		return(t.test(M_D_GD_mae_z[geneInd],M_D_GD_mae_z[1:105],alternative="greater")$p.value)
	}else{
		return(NA)
	}
})

M_D_GD_mae_padj<-p.adjust(M_D_GD_mae_pvalue,"BH")
M_D_GD_mae_df_z<-data.frame(gene=unique(gene_extend),padj=M_D_GD_mae_padj)
M_D_GD_mae_df_z[which(M_D_GD_mae_df_z$padj<0.05),]






















cor(combined_result[[3]],combined_df$GD,use="pairwise.complete.obs")^2
#plot(combined_result[[3]],combined_df$GD)
dataTemp<-data.frame(obs=combined_df$GD,pred=combined_result[[3]])
ggplot(dataTemp,aes(x=obs,y=pred))+
geom_point(alpha=0.2,stroke=0.1)+
geom_abline(intercept=0,slope=1,col="red")+
scale_x_continuous(limits=c(-5,5))+
scale_y_continuous(limits=c(-5,5))

1-sum((combined_result[[3]]-combined_df$GD)^2)/sum(combined_df$GD^2)















# x 
M_D_GD_coef_x<-t(sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	M_D_GD_model<-estimate_phi(M_D_GD_df$GD, M_D_GD_df$D, M_D_GD_df$M, phi_start = 0.618)
	return(c(M_D_GD_model[[1]],M_D_GD_model[[2]]))

}))

apply(M_D_GD_coef_x,2,mean)
apply(M_D_GD_coef_x,2,median)

# combine results

combined_df<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_result<-estimate_phi(combined_df$GD, combined_df$D, combined_df$M, phi_start = 0.8)
combined_result[[1]]
combined_result[[2]]
cor(combined_result[[3]],combined_df$GD,use="pairwise.complete.obs")^2
#plot(combined_result[[3]],combined_df$GD)
dataTemp<-data.frame(obs=combined_df$GD,pred=combined_result[[3]])
ggplot(dataTemp,aes(x=obs,y=pred))+
geom_point(alpha=0.2,stroke=0.1)+
geom_abline(intercept=0,slope=1,col="red")+
scale_x_continuous(limits=c(-5,5))+
scale_y_continuous(limits=c(-5,5))

1-sum((combined_result[[3]]-combined_df$GD)^2)/sum(combined_df$GD^2)

# y
M_D_GD_coef_y<-t(sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	M_D_GD_model<-estimate_phi(M_D_GD_df$GD, M_D_GD_df$D, M_D_GD_df$M, phi_start = 0.618)
	return(c(M_D_GD_model[[1]],M_D_GD_model[[2]]))

}))

apply(M_D_GD_coef_y,2,mean)
apply(M_D_GD_coef_y,2,median)

#
combined_df<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_result<-estimate_phi(combined_df$GD, combined_df$D, combined_df$M, phi_start = 0.8)
combined_result[[1]]
combined_result[[2]]
cor(combined_result[[3]],combined_df$GD,use="pairwise.complete.obs")^2
#plot(combined_result[[3]],combined_df$GD)
dataTemp<-data.frame(obs=combined_df$GD,pred=combined_result[[3]])
ggplot(dataTemp,aes(x=obs,y=pred))+
geom_point(alpha=0.2,stroke=0.1)+
geom_abline(intercept=0,slope=1,col="red")+
scale_x_continuous(limits=c(-5,5))+
scale_y_continuous(limits=c(-5,5))

1-sum((combined_result[[3]]-combined_df$GD)^2)/sum(combined_df$GD^2)


# z
M_D_GD_coef_z<-t(sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	if(length(intersect(emb,c(169)))==0){
		lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
		M_D_GD_df<-data.frame(
			M=c(lineage_velocity_emb[1:8,]),
			D=c(lineage_velocity_emb[2:9,]),
			GD=c(lineage_velocity_emb[3:10,])
		)
		M_D_GD_model<-estimate_phi(M_D_GD_df$GD, M_D_GD_df$D, M_D_GD_df$M, phi_start = 0.618)
		return(c(M_D_GD_model[[1]],M_D_GD_model[[2]]))	
	}else{
		return(c(NA,NA))
	}

}))

apply(M_D_GD_coef_z,2,mean,na.rm=T)
apply(M_D_GD_coef_z,2,median,na.rm=T)

#
combined_df<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_result<-estimate_phi(combined_df$GD, combined_df$D, combined_df$M, phi_start = 0.618)
combined_result[[1]]
combined_result[[2]]
cor(combined_result[[3]],combined_df$GD,use="pairwise.complete.obs")^2

cor(combined_df$GD, combined_df$D)^2
cor(combined_df$GD, combined_df$M)^2

summary(lm(GD~D+0,data=combined_df))
summary(lm(GD~M+0,data=combined_df))
summary(lm(GD~D+M+0,data=combined_df))
cor(combined_df$GD, -combined_df$D*0.618^2-combined_df$M*0.618^3)^2


#plot(combined_result[[3]],combined_df$GD)
dataTemp<-data.frame(obs=combined_df$GD,pred=combined_result[[3]])
ggplot(dataTemp,aes(x=obs,y=pred))+
geom_point(alpha=0.2,stroke=0.1)+
geom_abline(intercept=0,slope=1,col="red")+
scale_x_continuous(limits=c(-5,5))+
scale_y_continuous(limits=c(-5,5))

1-sum((combined_result[[3]]-combined_df$GD)^2)/sum(combined_df$GD^2)



dataTemp<-data.frame(
	ratio=c(M_D_GD_coef_x[1:105,1],M_D_GD_coef_y[1:105,1],M_D_GD_coef_z[1:105,1]),
	coordinate=c(rep("x",105),rep("y",105),rep("z",105))
)
ggplot(dataTemp,aes(x=coordinate,y=ratio,col=coordinate))+
geom_boxplot()+
geom_jitter(width=0.15)+
geom_hline(yintercept=0.618,lty=2)

ggplot(dataTemp,aes(x=ratio,col=coordinate))+
geom_density()+
geom_vline(xintercept=0.618,lty=2)+
geom_vline(xintercept=c(mean(M_D_GD_coef_x[1:105,1]),mean(M_D_GD_coef_y[1:105,1]),mean(M_D_GD_coef_z[1:105,1])),lty=2,col="red")

dataTemp<-data.frame(
	ratio=c(mean(M_D_GD_coef_x[c(1:28,30:105),1]),mean(M_D_GD_coef_y[c(1:28,30:105),1]),mean(M_D_GD_coef_z[c(1:28,30:105),1])),
	SE=c(mean(M_D_GD_coef_x[c(1:28,30:105),2]),mean(M_D_GD_coef_y[c(1:28,30:105),2]),mean(M_D_GD_coef_z[c(1:28,30:105),2])),
	coordinate=c("x","y","z")
)
ggplot(dataTemp,aes(x=coordinate,y=ratio,ymax=ratio+SE*1.96,ymin=ratio-SE*1.96))+
geom_point()+
geom_errorbar(width=0.2)+
geom_hline(yintercept=0.618,col="red",lty=2)+
scale_y_continuous(limits=c(0,1))


dataTemp<-data.frame(
	ratio=c(mean(M_D_GD_coef_x[which(M_D_GD_coef_x[,2]<1),1]),mean(M_D_GD_coef_y[which(M_D_GD_coef_y[,2]<1),1]),mean(M_D_GD_coef_z[which(M_D_GD_coef_z[,2]<1),1])),
	SE=c(mean(M_D_GD_coef_x[which(M_D_GD_coef_x[,2]<1),2]),mean(M_D_GD_coef_y[which(M_D_GD_coef_y[,2]<1),2]),mean(M_D_GD_coef_z[which(M_D_GD_coef_z[,2]<1),2])),
	coordinate=c("x","y","z")
)
ggplot(dataTemp,aes(x=coordinate,y=ratio,ymax=ratio+SE*1.96,ymin=ratio-SE*1.96))+
geom_point()+
geom_errorbar(width=0.2)+
geom_hline(yintercept=0.618,col="red",lty=2)+
scale_y_continuous(limits=c(0,1))





# x 
M_D_GD_coef_x<-t(sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	M_D_GD_model<-lm(GD~M+D+0,data=M_D_GD_df)
	return(coefficients(summary(M_D_GD_model))[,1])	

}))
apply(M_D_GD_coef_x,2,mean)
apply(M_D_GD_coef_x,2,median)
mean(M_D_GD_coef_x[,1]/M_D_GD_coef_x[,2])
median(M_D_GD_coef_x[,1]/M_D_GD_coef_x[,2])

# y
M_D_GD_coef_y<-t(sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	M_D_GD_model<-lm(GD~M+D+0,data=M_D_GD_df)
	return(coefficients(summary(M_D_GD_model))[,1])	

}))
apply(M_D_GD_coef_y,2,mean)
apply(M_D_GD_coef_y,2,median)
mean(M_D_GD_coef_y[,1]/M_D_GD_coef_y[,2])
median(M_D_GD_coef_y[,1]/M_D_GD_coef_y[,2])

# z
M_D_GD_coef_z<-t(sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	M_D_GD_model<-lm(GD~M+D+0,data=M_D_GD_df)
	return(coefficients(summary(M_D_GD_model))[,1])	

}))
apply(M_D_GD_coef_z,2,mean)
apply(M_D_GD_coef_z,2,median)
mean(M_D_GD_coef_z[,1]/M_D_GD_coef_z[,2])
median(M_D_GD_coef_z[,1]/M_D_GD_coef_z[,2])

# xyz
M_D_GD_coef_xyz<-t(sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb_x<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df_x<-data.frame(
		M=c(lineage_velocity_emb_x[1:8,]),
		D=c(lineage_velocity_emb_x[2:9,]),
		GD=c(lineage_velocity_emb_x[3:10,])
	)
	lineage_velocity_emb_y<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df_y<-data.frame(
		M=c(lineage_velocity_emb_y[1:8,]),
		D=c(lineage_velocity_emb_y[2:9,]),
		GD=c(lineage_velocity_emb_y[3:10,])
	)
	lineage_velocity_emb_z<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df_z<-data.frame(
		M=c(lineage_velocity_emb_z[1:8,]),
		D=c(lineage_velocity_emb_z[2:9,]),
		GD=c(lineage_velocity_emb_z[3:10,])
	)
	
	M_D_GD_df<-rbind(M_D_GD_df_x,M_D_GD_df_y,M_D_GD_df_z)
	
	M_D_GD_model<-lm(GD~M+D+0,data=M_D_GD_df)
	return(coefficients(summary(M_D_GD_model))[,1])	

}))
apply(M_D_GD_coef_xyz,2,mean)
apply(M_D_GD_coef_xyz,2,median)
mean(M_D_GD_coef_xyz[,1]/M_D_GD_coef_xyz[,2])
median(M_D_GD_coef_xyz[,1]/M_D_GD_coef_xyz[,2])




M_D_GD_coef<-t(sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	M_D_GD_model<-lm(GD~M+D+0,data=M_D_GD_df)
	return(coefficients(summary(M_D_GD_model))[,1])	

}))


par(mfrow=c(2,3))

plot(M_D_GD_coef[,1])
abline(h=0)

plot(M_D_GD_coef[,2])
abline(h=0)

plot(M_D_GD_coef)
abline(0,1,col="red")
abline(0,1.26,col="blue")

summary(lm(D~M+0,data.frame(M_D_GD_coef)))
summary(lm(M~D+0,data.frame(M_D_GD_coef)))

plot(M_D_GD_coef[,2:1])
abline(0,1,col="red")
abline(0,0.7,col="blue")


lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd

N=1000
eps<-rnorm(1000,0,10)
x<-c()
x[1]=1.3
x[2]=0.2
rou=0.618
rou=0.9
rou=0.5
rou=0.5
rou=0.7
rou=0.2

for(i in 3:N){
	x[i]<--rou*(rou^2*x[i-2]+x[i-1]*rou)+eps[i]
}

plot(x,type="b")
abline(h=0,lty=2)
mean(x)
var(x)

#

N=100
eps<-rnorm(1000,0,10)
x<-c()
x[1]=1.3
x[2]=0.2
rou=0.618
c=(1-rou)/rou


for(i in 3:N){
	x[i]<--c*(rou^2*x[i-2]+x[i-1]*rou)+eps[i]
}

plot(x,type="b")
abline(h=0,lty=2)
mean(x)
var(x)



cor(x[1:(N-1)],x[2:N])
cor(x[1:(N-2)],x[3:N])
cor(x[1:(N-3)],x[4:N])
self_corr_x[,1]


x[i]=-0.382*(0.618*x[i-2]+x[i-1])=0

x[i]+φ^2*(x[i-1]+φ*x[i-2])=0

x[i]+φ^2*x[i-1]+φ^3*x[i-2]=0

(1+φ)*x[i]+φ*x[i-1]+φ^2*x[i-2]=0


φ*x[i]=-(x[i]+φ*x[i-1]+φ^2*x[i-2])

(1+1/φ)*x[i]=-(x[i-1]+φ*x[i-2])


X[t]=-0.382*(x[t-1]+0.618*x[t-2])

boxplot(list((M_D_GD_coef[,1]/M_D_GD_coef[,2])[embryo_lethal=="0"],(M_D_GD_coef[,1]/M_D_GD_coef[,2])[embryo_lethal=="1"]))


### 

coef_M_expected<-apply(M_D_GD_coef,2,median)[1]
coef_D_expected<-apply(M_D_GD_coef,2,median)[2]


M_D_GD_R<-sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	
	R<-cor(M_D_GD_df$GD,M_D_GD_df$M*coef_M_expected+M_D_GD_df$D*coef_D_expected,use="pairwise.complete.obs")
	return(R)	

})

#x 
M_D_GD_R_x<-sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	
	R<-cor(M_D_GD_df$GD,M_D_GD_df$M*(-0.618*0.618)+M_D_GD_df$D*(-0.618),use="pairwise.complete.obs")
	return(R)	

})

#y
M_D_GD_R_y<-sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	
	R<-cor(M_D_GD_df$GD,M_D_GD_df$M*(-0.618*0.618)+M_D_GD_df$D*(-0.618),use="pairwise.complete.obs")
	return(R)	

})

#z
M_D_GD_R_z<-sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	
	R<-cor(M_D_GD_df$GD,M_D_GD_df$M*(-0.618*0.618)+M_D_GD_df$D*(-0.618),use="pairwise.complete.obs")
	return(R)	

})

dataTemp<-data.frame(
R=c(M_D_GD_R_x,M_D_GD_R_y,M_D_GD_R_z),
Coordinate=c(rep("x",length(M_D_GD_R_x)),rep("y",length(M_D_GD_R_y)),rep("z",length(M_D_GD_R_z))),
type=rep(c(rep("WT",105),rep("Perturbation",2039-105)),3)
)
ggplot(dataTemp,aes(x=R^2,col=Coordinate))+
geom_density()

ggplot(dataTemp,aes(x=R^2,col=type))+
geom_density()



## for R 

dataTemp<-data.frame(
	R=M_D_GD_R,
	type=c(rep("WT",105),rep("RNAi",1928-105),rep("Env",2039-1928)),
	lethal=embryo_lethal,
	ifsig=rep("insig",length(M_D_GD_R)))
sig_gene_extend_index<-unlist(lapply(sig_gene_extend,function(x){which(gene_extend==x)}))
dataTemp$ifsig[sig_gene_extend_index]<-"sig"

ggplot(dataTemp,aes(x=type,y=R))+
geom_boxplot()


ggplot(dataTemp,aes(x=lethal,y=R))+
geom_boxplot()

t.test(dataTemp$R[dataTemp$lethal=="0"],dataTemp$R[dataTemp$lethal=="1"])

ggplot(dataTemp,aes(x=ifsig,y=R))+
geom_boxplot()


ggplot(dataTemp,aes(x=R,col=lethal))+
geom_density()


ggplot(dataTemp,aes(x=R,col=ifsig))+
geom_density()


## for MAE
M_D_GD_mae<-sapply(1:dim(velocity)[1],function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	
	df<-na.omit(data.frame(obs=M_D_GD_df$GD,pred=M_D_GD_df$M*coef_M_expected+M_D_GD_df$D*coef_D_expected))
	mae=sum(abs(df[,1]-df[,2]))/dim(df)[1]
	return(mae)	

})

dataTemp<-data.frame(
	mae=M_D_GD_mae,
	type=c(rep("WT",105),rep("RNAi",1928-105),rep("Env",2039-1928)),
	lethal=embryo_lethal,
	ifsig=rep("insig",length(M_D_GD_R)))
sig_gene_extend_index<-unlist(lapply(sig_gene_extend,function(x){which(gene_extend==x)}))
dataTemp$ifsig[sig_gene_extend_index]<-"sig"

ggplot(dataTemp,aes(x=type,y=mae))+
geom_boxplot()

# 
M_D_GD_mae_pvalue<-sapply(unique(gene_extend),function(x){
	geneInd<-which(gene_extend==x)
	if(length(geneInd)>1){
		return(t.test(M_D_GD_mae[geneInd],M_D_GD_mae[1:105],alternative="greater")$p.value)
	}else{
		return(NA)
	}
})

M_D_GD_mae_padj<-p.adjust(M_D_GD_mae_pvalue,"BH")
M_D_GD_mae_df<-data.frame(gene=unique(gene_extend),padj=M_D_GD_mae_padj)
M_D_GD_mae_df[which(M_D_GD_mae_df$padj<0.05),]

#
M_D_GD_mae_z<-sapply(1:length(gene_extend),function(i){
	return((M_D_GD_mae[i]-mean(M_D_GD_mae[1:105]))/sd(M_D_GD_mae[1:105]))
})

M_D_GD_mae_z_df<-data.frame(gene=gene_extend,z=M_D_GD_mae_z)
M_D_GD_mae_z_df$adjp<-p.adjust(1-pnorm(M_D_GD_mae_z_df$z),"BH")
M_D_GD_mae_z_df$count<-sapply(M_D_GD_mae_z_df$gene,function(x){length(which(M_D_GD_mae_z_df$gene==x&M_D_GD_mae_z_df$adjp<0.05))})

M_D_GD_mae_z_df[M_D_GD_mae_z_df$z>qnorm(0.98),]
paste(names(table(M_D_GD_mae_z_df$gene[M_D_GD_mae_z_df$z>qnorm(0.98)])),collapse=",")
paste(names(table(M_D_GD_mae_z_df$gene[M_D_GD_mae_z_df$adjp<0.05])),collapse="\t")
sigTemp<-M_D_GD_mae_z_df$gene[M_D_GD_mae_z_df$adjp<0.05]
sigGene_temp<-setdiff(names(which(table(sigTemp)>1)),"ctr")
sigGene_temp<-tolower(sigGene_temp)

#write.csv(M_D_GD_mae_z_df,"M_D_GD_mae_z_df.csv")

gene_mapped_temp <- mapIds(
  org.Ce.eg.db,
  keys = sigGene_temp,
  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
  column = "ENTREZID"
)

sigGene_temp[which(is.na(gene_mapped_temp))]<-c(
"C48B6.2",
"eif-3.E",
"hrpk-1",
"hs25m",
"ints-1",
"sa600",
"ints-4",
"slc-30A5",
"Y106G6E.4",
"edc-4",
"Y47D9A.1"
)

gene_mapped_temp <- mapIds(
  org.Ce.eg.db,
  keys = sigGene_temp,
  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
  column = "ENTREZID"
)

ego_temp <- enrichGO(
  gene = gene_mapped_temp,
  OrgDb = org.Ce.eg.db,
  ont = "ALL",  # or BP, CC, MF
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

data.frame(ego_temp)

barplot(ego_temp, showCategory = 100, label_format = 100)



ggplot(M_D_GD_mae_z_df,aes(x=M_D_GD_mae_z_df$gene=="ctr",y=z))+
geom_boxplot()


ggplot(dataTemp,aes(x=lethal,y=mae))+
geom_boxplot()




t.test(dataTemp$mae[dataTemp$lethal=="0"],dataTemp$mae[dataTemp$lethal=="1"])

ggplot(dataTemp,aes(x=ifsig,y=mae))+
geom_boxplot()


ggplot(dataTemp,aes(x=mae,col=lethal))+
geom_density()


ggplot(dataTemp,aes(x=mae,col=ifsig))+
geom_density()

# Load required libraries
library(pROC)      # For ROC/AUC calculations
library(ggplot2)   # For enhanced plotting
library(caret)     # For data splitting (optional but recommended)
library(tidyverse) # For data manipulation

# ==================== STEP 1: Prepare Data ====================
# Use mtcars dataset, convert am to factor for binary classification

# Split data into training and testing sets (80/20 split)
set.seed(123)  # For reproducibility
train_index <- sample(1:dim(dataTemp)[1],400)
train_data <- dataTemp[train_index, ]
test_data <- dataTemp[-train_index, ]

# ==================== STEP 2: Fit Logistic Regression Model ====================
# Model predicting transmission type (am) based on horsepower (hp) and weight (wt)
logit_model <- glm(as.numeric(lethal) ~ mae, 
                   data = train_data, 
                   family = binomial(link = "logit"))

# View model summary
summary(logit_model)

# Interpretation of coefficients:
# - Coefficients are in log-odds
# - Positive coefficient = increases probability of am=1 (manual)
# - Negative coefficient = decreases probability of am=1

# ==================== STEP 3: Make Predictions ====================
# Get predicted probabilities for the test set
test_predictions <- predict(logit_model, 
                            newdata = test_data, 
                            type = "response")

# Convert to factor with 0.5 cutoff (standard threshold)
predicted_class <- ifelse(test_predictions > 0.5, 1, 0)
predicted_class <- factor(predicted_class, levels = c(0, 1))

# ==================== STEP 4: Calculate and Plot ROC/AUC ====================
# Calculate ROC curve
roc_result <- roc(response = as.numeric(test_data$lethal), 
                  predictor = test_predictions)

# Print AUC
cat("AUC:", auc(roc_result), "\n")

# Basic ROC plot
plot(roc_result, 
     main = paste("ROC Curve (AUC =", round(auc(roc_result), 3), ")"),
     col = "blue", 
     lwd = 2,
     print.auc = TRUE,
     print.auc.x = 0.6, 
     print.auc.y = 0.4)

# Add diagonal reference line
abline(a = 0, b = 1, lty = 2, col = "gray")

# ==================== OPTIONAL: Enhanced ROC Plot with ggplot2 ====================
# Create dataframe for ROC curve
roc_df <- data.frame(
  sensitivity = roc_result$sensitivities,
  specificity = roc_result$specificities,
  thresholds = roc_result$thresholds
)

# Calculate 1 - specificity
roc_df$one_minus_specificity <- 1 - roc_df$specificity

# Create enhanced ROC plot
ggplot(roc_df, aes(x = one_minus_specificity, y = sensitivity)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_polygon(aes(x = one_minus_specificity, y = sensitivity), 
               fill = "skyblue", alpha = 0.3) +
  labs(
    title = paste("ROC Curve for Logistic Regression Model"),
    subtitle = paste("AUC =", round(auc(roc_result), 3)),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  annotate("text", x = 0.7, y = 0.3, 
           label = paste("AUC =", round(auc(roc_result), 3)),
           size = 5, color = "darkred", fontface = "bold") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# ==================== STEP 5: Model Evaluation Metrics ====================
# Confusion matrix
conf_matrix <- confusionMatrix(data = predicted_class, 
                               reference = test_data$am, 
                               positive = "1")
print(conf_matrix)

# Calculate


