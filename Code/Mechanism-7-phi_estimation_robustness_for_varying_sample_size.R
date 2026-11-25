
### nls 
estimate_phi <- function(GD, D, M, phi_start = 0.6) {
  df <- data.frame(GD = GD, D = D, M = M)

  fit <- nls(
    GD ~ -phi^2 * D - phi^3 * M,
    data = df,
    start = list(phi = phi_start),
    algorithm = "port",
    lower = 0,
    upper = 2
  )

  list(
    phi_hat = coef(fit),
    se = sqrt(diag(vcov(fit))),
    fitted = fitted(fit),
    residuals = resid(fit)
  )
}


phi_with_sample<-function(sample_size, n_rep){
	set.seed(1)
	phi_results<-do.call(rbind,
		lapply(sample_size,function(ss){
			print(paste0("sample size: ",ss))
			temp_sample_results<-do.call(rbind,
				lapply(1:n_rep,function(nrep){
					print(paste0("sample size: ",ss," replicate: ", nrep))
					sampleIndex<-sample(1:105,ss)
					## x
					combined_df_x<-na.omit(do.call(rbind,lapply(sampleIndex,function(emb){
						
						lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
						M_D_GD_df<-data.frame(
							M=c(lineage_velocity_emb[1:8,]),
							D=c(lineage_velocity_emb[2:9,]),
							GD=c(lineage_velocity_emb[3:10,]),
							embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
						)
						return(M_D_GD_df)
					})))
					
					#y
					combined_df_y<-na.omit(do.call(rbind,lapply(sampleIndex,function(emb){
						
						lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
						M_D_GD_df<-data.frame(
							M=c(lineage_velocity_emb[1:8,]),
							D=c(lineage_velocity_emb[2:9,]),
							GD=c(lineage_velocity_emb[3:10,]),
							embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
						)
						return(M_D_GD_df)
					})))
					
					#z
					combined_df_z<-na.omit(do.call(rbind,lapply(sampleIndex,function(emb){
						
						lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
						M_D_GD_df<-data.frame(
							M=c(lineage_velocity_emb[1:8,]),
							D=c(lineage_velocity_emb[2:9,]),
							GD=c(lineage_velocity_emb[3:10,]),
							embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
						)
						return(M_D_GD_df)
					})))	
					
					#xyz 
					combined_df_xyz<-rbind(combined_df_x,combined_df_y,combined_df_z)	
					
					## NLS
					
					combined_result_x<-estimate_phi(combined_df_x$GD, combined_df_x$D, combined_df_x$M, phi_start = 0.618)
					combined_result_y<-estimate_phi(combined_df_y$GD, combined_df_y$D, combined_df_y$M, phi_start = 0.618)
					combined_result_z<-estimate_phi(combined_df_z$GD, combined_df_z$D, combined_df_z$M, phi_start = 0.618)
					combined_result_xyz<-estimate_phi(combined_df_xyz$GD, combined_df_xyz$D, combined_df_xyz$M, phi_start = 0.618)
					
					dataTemp<-data.frame(replicates=rep(nrep,4),
						sample_size=rep(ss,4),
						coordinate=c("x","y","z","xyz"),
						phi_estimate=c(combined_result_x[[1]],combined_result_y[[1]],combined_result_z[[1]],combined_result_xyz[[1]]),
						phi_se=c(combined_result_x[[2]],combined_result_y[[2]],combined_result_z[[2]],combined_result_xyz[[2]]))
						
					return(dataTemp)
				})
			)
			
			return(temp_sample_results)
			
		})
	
	)

	return(phi_results)
	
}

robust_phi_sampleSize<-phi_with_sample(c(5,10,20,30,40,50,60,70,80,90,105),10)


robust_phi_sampleSize$coordinate<-factor(robust_phi_sampleSize$coordinate,level=c("x","y","z","xyz"))

###

library(dplyr)
library(ggplot2)


## create one shared jitter position so points and error bars stay aligned
jit <- position_jitter(width = 2, height = 0, seed = 1)

p4 <- ggplot(
  robust_phi_sampleSize,
  aes(
    x = sample_size,
    y = phi_estimate,
    ymin = phi_estimate - phi_se,
    ymax = phi_estimate + phi_se,
    color = factor(replicates)
  )
) +
  geom_errorbar(
    position = jit,
    width = 1,
    alpha = 0.45,
    color = "grey45"
  ) +
  geom_point(
    position = jit,
    size = 1
  ) +
  geom_hline(yintercept = 0.618, lty = 2, lwd = 0.5) +
  facet_wrap(~coordinate, ncol = 2, scales = "free_y") +
  scale_x_continuous(
    breaks = sort(unique(robust_phi_sampleSize$sample_size))
  ) +
  labs(
    x = "Sample size",
    y = expression(hat(phi) %+-% SE),
    color = "Replicate",
    title = expression("Estimated " * phi * " ± standard error for each replicate"),
    subtitle = "Separate panels show x, y, z, and xyz coordinates"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "right"
  )

p4


