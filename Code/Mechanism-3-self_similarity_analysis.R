

# V+ -αβ(1-β^n)/(1-β)
### self-similarity process through continuous negative feedbacks
x<-c()
x[1]=1
x[2]=1-0.618^2
for(i in 3:11){
	x[i]=1-sum(0.618^c(2:i))
}

plot(x,type="b")

library(ggplot2)
dataTemp<-data.frame(generation=c(0:10),Vnoise=x)

ggplot(dataTemp,aes(x=generation,y=Vnoise))+
geom_point(size=0.5)+
geom_line(lwd=0.3)+
geom_hline(yintercept=c(0.236,0.382,0.618,1),lty=2,lwd=0.3,col="grey")+
scale_x_continuous(breaks=c(0:10))+
scale_y_continuous(breaks=c(0,0.236,0.382,0.618,1))+
labs(x="Generation interval",y="Residual of initial velocity noise",title="")+
theme(axis.text.y=element_text(angle=0,vjust=0.5,hjust=1))

ggsave("self-similarity_process.pdf",width=unit(2.2,"in"),height=unit(1.7,"in"))


### estimate ψ and φ through ordinary linear regression 

# combine wild types
combined_df<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	return(M_D_GD_df)
})))

summary(lm(GD~D+M+0,data=combined_df))

combine_D<-coefficients(summary(lm(GD~D+M+0,data=combined_df)))[1,1]
combine_M<-coefficients(summary(lm(GD~D+M+0,data=combined_df)))[2,1]
combine_r<-combine_M/combine_D
combine_c<--combine_D/combine_r

individual_D_mean<-apply(lm_M_D_GD_coef_x[1:105,],2,mean)[1]
individual_M_mean<-apply(lm_M_D_GD_coef_x[1:105,],2,mean)[2]
individual_r<-individual_M_mean/individual_D_mean
individual_c<--individual_D_mean/individual_r


# dataTemp<-data.frame(
	# coeff=c(combine_D,combine_M,combine_r,combine_c,individual_D_mean,individual_M_mean,individual_r,individual_c),
	# name=factor(c("combine_D","combine_M","combine_r","combine_c","individual_D_mean","individual_M_mean","individual_r","individual_c"),
		# level=c("combine_D","combine_M","combine_r","combine_c","individual_D_mean","individual_M_mean","individual_r","individual_c"))
# )

dataTemp<-data.frame(
	coeff=c(combine_D,combine_M,combine_r,combine_c),
	name=factor(c("coeffD","coeffM","gamma","psi"),
		level=c("coeffD","coeffM","gamma","psi"))
)
phi=(sqrt(5)-1)/2

ggplot(dataTemp,aes(x=name,y=coeff))+
geom_bar(stat="identity",position="dodge",width=0.5,col="#619CFF",fill="#619CFF")+
geom_hline(yintercept=c(phi,-phi^2,-phi^3),col="#F8766D",lty=2,lwd=0.5)+
geom_hline(yintercept=0,col="black")+
#theme(axis.text.y=element_text(angle=0,vjust=0.5,hjust=1))+
scale_y_continuous(breaks=c(0,-0.236,-0.382,0.618,1))+
labs(x="",y="Coefficient values (OLS)",title="")

ggsave("coefficient_estimation_by_OLS.pdf",width=unit(2.2,"in"),height=unit(2.1,"in"))

# 

# combine wild types
combined_df_x<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_df_y<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_df_z<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	return(M_D_GD_df)
})))

summary(lm(GD~D+M+0,data=combined_df_x))
summary(lm(GD~D+M+0,data=combined_df_y))
summary(lm(GD~D+M+0,data=combined_df_z))

combined_df_xyz<-rbind(combined_df_x,combined_df_y,combined_df_z)
summary(lm(GD~D+M+0,data=combined_df_xyz))







### coefficient estimation by nls
### nls 
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
      # Success: return the result
      return(list(
        phi_hat = coef(fit),
        se = sqrt(diag(vcov(fit))),
        fitted = fitted(fit),
        residuals = resid(fit)
      ))
    }
  }
  
  # If all starting values fail, return NA (should be rare)
  warning("No starting value converged for this embryo.")
  return(list(phi_hat = NA, se = NA, fitted = NA, residuals = NA))
}


## x
combined_df_x<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

combined_result_x<-estimate_phi(combined_df_x$GD, combined_df_x$D, combined_df_x$M, phi_start = 0.618)
combined_result_x[[1]]
combined_result_x[[2]]




population_control_x<-sapply(1:105,function(emb){
	estimate_phi(combined_df_x[combined_df_x$embryoID==emb,]$GD, combined_df_x[combined_df_x$embryoID==emb,]$D, combined_df_x[combined_df_x$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_control_x,0))
mean(setdiff(population_control_x,0))
median(setdiff(population_control_x,0))
hist(population_control_x)
mean(population_control_x)
median(population_control_x)





## x_unstandardized
combined_df_x_unstandardized<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	return(M_D_GD_df)
})))

combined_result_x_unstandardized<-estimate_phi(combined_df_x_unstandardized$GD, combined_df_x_unstandardized$D, combined_df_x_unstandardized$M, phi_start = 0.618)
combined_result_x_unstandardized[[1]]
combined_result_x_unstandardized[[2]]



#y
combined_df_y<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

combined_result_y<-estimate_phi(combined_df_y$GD, combined_df_y$D, combined_df_y$M, phi_start = 0.618)
combined_result_y[[1]]
combined_result_y[[2]]


population_control_y<-sapply(1:105,function(emb){
	print(emb)
	estimate_phi(combined_df_y[combined_df_y$embryoID==emb,]$GD, combined_df_y[combined_df_y$embryoID==emb,]$D, combined_df_y[combined_df_y$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_control_y,0))
mean(setdiff(population_control_y,0))
median(setdiff(population_control_y,0))
hist(population_control_y)
mean(population_control_y)
median(population_control_y)


#z
combined_df_z<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

combined_result_z<-estimate_phi(combined_df_z$GD, combined_df_z$D, combined_df_z$M, phi_start = 0.618)
combined_result_z[[1]]
combined_result_z[[2]]


population_control_z<-sapply(1:105,function(emb){
	print(emb)
	estimate_phi(combined_df_z[combined_df_z$embryoID==emb,]$GD, combined_df_z[combined_df_z$embryoID==emb,]$D, combined_df_z[combined_df_z$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_control_z,0))
mean(setdiff(population_control_z,0))
median(setdiff(population_control_z,0))
hist(population_control_z)
mean(population_control_z)
median(population_control_z)


#xyz 
combined_df_xyz<-rbind(combined_df_x,combined_df_y,combined_df_z)
combined_result_xyz<-estimate_phi(combined_df_xyz$GD, combined_df_xyz$D, combined_df_xyz$M, phi_start = 0.618)
combined_result_xyz[[1]]
combined_result_xyz[[2]]


population_control_xyz<-sapply(1:105,function(emb){
	print(emb)
	estimate_phi(combined_df_xyz[combined_df_xyz$embryoID==emb,]$GD, combined_df_xyz[combined_df_xyz$embryoID==emb,]$D, combined_df_xyz[combined_df_xyz$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_control_xyz,0))
mean(setdiff(population_control_xyz,0))
median(setdiff(population_control_xyz,0))
hist(population_control_xyz)
mean(population_control_xyz)
median(population_control_xyz)





dataTemp<-data.frame(
coeff=c(combined_result_x[[1]],combined_result_y[[1]],combined_result_z[[1]]),
se=c(combined_result_x[[2]],combined_result_y[[2]],combined_result_z[[2]]),
coordinate=c("x","y","z"))

ggplot(dataTemp,aes(x=coordinate,y=coeff,ymin=coeff-1.96*se,ymax=coeff+1.96*se,col=coordinate))+
geom_point(size=0.1)+
geom_errorbar(width=0.2,lwd=0.1)+
geom_hline(yintercept=0.618,col="grey",lty=2,lwd=0.1)+
scale_y_continuous(limits=c(0,1),breaks=c(0,0.618,1))+
theme(legend.position=c(0.8,0.2))+
labs(x="Coordinate",y="Parameter estimation",title="NLS")+
guides(color = guide_legend(
    override.aes = list(size = 2,alpha=1),        # Increase key size (points in legend)
    keywidth = 0.5,                       # Reduce width of key box
    keyheight = 0.5                       # Reduce height of key box
  )) +
  theme(
	legend.background=element_blank(),
	legend.frame=element_blank(),
	legend.box.background=element_blank(),
    legend.spacing.x = unit(0.1, 'cm'),   # Reduce horizontal spacing between keys
    legend.spacing.y = unit(0.1, 'cm'),   # Reduce vertical spacing between keys
    legend.key.width = unit(0.3, 'cm'),   # Reduce spacing between key and text
    legend.key.height = unit(0.3, 'cm')   # Reduce height of key background
  )


ggsave("coefficient_estimation_by_NLS.pdf",width=unit(2.2,"in"),height=unit(2.1,"in"))


### comparison between OLS and NLS 
dataTemp<-data.frame(
	fitted_by_phi=c(combined_result_x[[3]],combined_result_y[[3]],combined_result_z[[3]]),
	fitted_by_lm=c(fitted(lm(GD~D+M, data=combined_df_x)),fitted(lm(GD~D+M, data=combined_df_y)),fitted(lm(GD~D+M, data=combined_df_z))),
	coordinate=c(rep("x",length(combined_result_x[[3]])),rep("y",length(combined_result_y[[3]])),rep("z",length(combined_result_z[[3]])))
)

ggplot(dataTemp,aes(x=fitted_by_lm,y=fitted_by_phi,col=coordinate))+
geom_point(size=0.15,stroke=0.1,alpha=0.3)+
geom_abline(slope=1,intercept=0,lwd=0.3)+
labs(x="Fitted values (OLS)",y="Fitted values (NLS)",title="OLS Vs NLS")+
theme(legend.position=c(0.8,0.2))+
guides(color = guide_legend(
    override.aes = list(size = 2,alpha=1),        # Increase key size (points in legend)
    keywidth = 0.5,                       # Reduce width of key box
    keyheight = 0.5                       # Reduce height of key box
  )) +
  theme(
	legend.background=element_blank(),
	legend.frame=element_blank(),
	legend.box.background=element_blank(),
    legend.spacing.x = unit(0.1, 'cm'),   # Reduce horizontal spacing between keys
    legend.spacing.y = unit(0.1, 'cm'),   # Reduce vertical spacing between keys
    legend.key.width = unit(0.3, 'cm'),   # Reduce spacing between key and text
    legend.key.height = unit(0.3, 'cm')   # Reduce height of key background
  )

ggsave("fitted_value_comparision_OLS_Vs_NLS.pdf",width=unit(2,"in"),height=unit(2.1,"in"))


### simulation for the conservation law of phi, smaller or larger value 
library(ggplot2)
library(dplyr)

## main figure to prove the phi's conservation law 
set.seed(123)

# parameters
n_samples <- 400
X0_samples <- rnorm(n_samples, mean = 0, sd = 1)

generations <- seq(0, 10, by = 1)

phi_low    <- 0.4
phi_golden <- (sqrt(5) - 1) / 2
phi_high   <- 0.8

phis   <- c(phi_low, phi_golden, phi_high)
labels <- c("Low", "Golden", "High")

# residual function (unchanged, correct)
compute_residual <- function(phi, X0, n) {
  if (n == 0) return(X0)
  X0 * (1 - sum(phi^(2:(n+1))))
}

# simulate
df2 <- expand.grid(
  X0 = X0_samples,
  generation = generations,
  phi = phis
)

df2$residual <- mapply(
  compute_residual,
  df2$phi,
  df2$X0,
  df2$generation
)

df2$phi_label <- factor(
  df2$phi,
  levels = phis,
  labels = labels
)


ggplot(df2, aes(x = generation, y = residual,group=generation)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 0.15,col="grey",stroke=0.1) +
  geom_boxplot(width = 0.5, lwd=0.15,outlier.shape = NA, fill = NA, color = "#F8766D") +
  facet_wrap(~ phi_label, nrow = 1) +
  labs(x = "Generation interval (t)", y = "Velocity residuals")+
  scale_y_continuous(breaks=c(-5,0,5))+
  scale_x_continuous(breaks=c(0,5,10))+
  theme(strip.placement="outside", strip.background=element_blank())

ggsave("simulation_for_golden_ratio_smaller_larger.pdf",width=unit(4,"in"),height=unit(1.7,"in"))


## only considering initial noise 

set.seed(123)

X0    <- 1
n_gen <- 10

df1 <- expand.grid(
  generation = 0:n_gen,
  phi = phis
)

df1$residual <- mapply(
  compute_residual,
  df1$phi,
  X0,
  df1$generation
)

df1$phi_label <- factor(df1$phi, levels = phis, labels = labels)

ggplot(df1, aes(generation, residual, color = phi_label)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.3) +
  geom_point(size = 0.5) +
  scale_x_continuous(breaks=c(0,5,10))+
  labs(
    x = "Generation",
    y = "Residual after cumulative feedback",
    color = "Φ",
    title = "Same initial noise decreases differently under distinct Φ"
  ) 
  
ggsave("simulation_for_golden_ratio_smaller_larger_only_initial noise.pdf",width=unit(2.4,"in"),height=unit(1.7,"in"))


## variance under different φ 

n_samples   <- 500
X0_samples  <- rnorm(n_samples, mean = 0, sd = 1)
generations <- 0:10

df2 <- expand.grid(
  X0 = X0_samples,
  generation = generations,
  phi = phis
)

df2$residual <- mapply(
  compute_residual,
  df2$phi,
  df2$X0,
  df2$generation
)

## compute variance using base R
var_list <- list()

for (p in phis) {
  for (g in generations) {
    vals <- df2$residual[df2$phi == p & df2$generation == g]
    var_list[[length(var_list) + 1]] <- data.frame(
      phi = p,
      generation = g,
      variance = var(vals)
    )
  }
}

var_df <- do.call(rbind, var_list)
var_df$phi_label <- factor(var_df$phi, levels = phis, labels = labels)

ggplot(var_df, aes(generation, variance, color = phi_label)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 0.5) +
  scale_x_continuous(breaks=c(0,5,10))+
  labs(
    x = "Generation",
    y = "Variance of residuals",
    color = "Φ",
    title = "Golden ratio uniquely drives variance to zero"
  )

ggsave("simulation_for_golden_ratio_smaller_larger_variance.pdf",width=unit(2.4,"in"),height=unit(1.7,"in"))



### noise or minor perturbations 
coef_M_expected=-phi_golden^3
coef_D_expected=-phi_golden^2
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

M_D_GD_mae_x_Zscore<-sapply(1:length(gene_extend),function(i){
	return((M_D_GD_mae_x[i]-mean(M_D_GD_mae_x[1:105]))/sd(M_D_GD_mae_x[1:105]))
})

M_D_GD_mae_y_Zscore<-sapply(1:length(gene_extend),function(i){
	return((M_D_GD_mae_y[i]-mean(M_D_GD_mae_y[1:105]))/sd(M_D_GD_mae_y[1:105]))
})

M_D_GD_mae_z_Zscore<-sapply(1:length(gene_extend),function(i){
	return((M_D_GD_mae_z[i]-mean(M_D_GD_mae_z[1:105]))/sd(M_D_GD_mae_z[1:105]))
})

# outlier embryos 
outlier_embryoInd<-which((M_D_GD_mae_x_Zscore+M_D_GD_mae_y_Zscore+M_D_GD_mae_z_Zscore)/3>qnorm(0.95))
outlier_embryoInd

length(setdiff(1:105,outlier_embryoInd))
length(setdiff(106:1928,outlier_embryoInd))
length(setdiff(1929:2039,outlier_embryoInd))

lethal_embryoInd<-which(embryo_lethal=="1")

## RNAi perturbation 
# x
RNAi_combined_df_x<-na.omit(do.call(rbind,lapply(setdiff(106:1928,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

RNAi_combined_result_x<-estimate_phi(RNAi_combined_df_x$GD, RNAi_combined_df_x$D, RNAi_combined_df_x$M, phi_start = 0.618)
RNAi_combined_result_x[[1]]
RNAi_combined_result_x[[2]]


population_RNAi_x<-sapply(setdiff(106:1928,lethal_embryoInd),function(emb){
	print(emb)
	estimate_phi(RNAi_combined_df_x[RNAi_combined_df_x$embryoID==emb,]$GD, RNAi_combined_df_x[RNAi_combined_df_x$embryoID==emb,]$D, RNAi_combined_df_x[RNAi_combined_df_x$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_RNAi_x,0))
mean(setdiff(population_RNAi_x,0))
median(setdiff(population_RNAi_x,0))
hist(population_RNAi_x)
mean(population_RNAi_x)
median(population_RNAi_x)




#y
RNAi_combined_df_y<-na.omit(do.call(rbind,lapply(setdiff(106:1928,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

RNAi_combined_result_y<-estimate_phi(RNAi_combined_df_y$GD, RNAi_combined_df_y$D, RNAi_combined_df_y$M, phi_start = 0.618)
RNAi_combined_result_y[[1]]
RNAi_combined_result_y[[2]]


population_RNAi_y<-sapply(setdiff(106:1928,lethal_embryoInd),function(emb){
	print(emb)
	estimate_phi(RNAi_combined_df_y[RNAi_combined_df_y$embryoID==emb,]$GD, RNAi_combined_df_y[RNAi_combined_df_y$embryoID==emb,]$D, RNAi_combined_df_y[RNAi_combined_df_y$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_RNAi_y,0))
mean(setdiff(population_RNAi_y,0))
median(setdiff(population_RNAi_y,0))
hist(population_RNAi_y)
mean(population_RNAi_y)
median(population_RNAi_y)


#z
RNAi_combined_df_z<-na.omit(do.call(rbind,lapply(setdiff(106:1928,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

RNAi_combined_result_z<-estimate_phi(RNAi_combined_df_z$GD, RNAi_combined_df_z$D, RNAi_combined_df_z$M, phi_start = 0.618)
RNAi_combined_result_z[[1]]
RNAi_combined_result_z[[2]]


population_RNAi_z<-sapply(setdiff(106:1928,lethal_embryoInd),function(emb){
	print(emb)
	estimate_phi(RNAi_combined_df_z[RNAi_combined_df_z$embryoID==emb,]$GD, RNAi_combined_df_z[RNAi_combined_df_z$embryoID==emb,]$D, RNAi_combined_df_z[RNAi_combined_df_z$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_RNAi_z,0))
mean(setdiff(population_RNAi_z,0))
median(setdiff(population_RNAi_z,0))
hist(population_RNAi_z)
mean(population_RNAi_z)
median(population_RNAi_z)



#xyz 
RNAi_combined_df_xyz<-rbind(RNAi_combined_df_x,RNAi_combined_df_y,RNAi_combined_df_z)
RNAi_combined_result_xyz<-estimate_phi(RNAi_combined_df_xyz$GD, RNAi_combined_df_xyz$D, RNAi_combined_df_xyz$M, phi_start = 0.618)
RNAi_combined_result_xyz[[1]]
RNAi_combined_result_xyz[[2]]


population_RNAi_xyz<-sapply(setdiff(106:1928,lethal_embryoInd),function(emb){
	print(emb)
	estimate_phi(RNAi_combined_df_xyz[RNAi_combined_df_xyz$embryoID==emb,]$GD, RNAi_combined_df_xyz[RNAi_combined_df_xyz$embryoID==emb,]$D, RNAi_combined_df_xyz[RNAi_combined_df_xyz$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_RNAi_xyz,0))
mean(setdiff(population_RNAi_xyz,0))
median(setdiff(population_RNAi_xyz,0))
hist(population_RNAi_xyz)
mean(population_RNAi_xyz)
median(population_RNAi_xyz)





## Env perturbation 
# x
Env_combined_df_x<-na.omit(do.call(rbind,lapply(setdiff(1929:2039,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

Env_combined_result_x<-estimate_phi(Env_combined_df_x$GD, Env_combined_df_x$D, Env_combined_df_x$M, phi_start = 0.618)
Env_combined_result_x[[1]]
Env_combined_result_x[[2]]


population_Env_x<-sapply(setdiff(1929:2039,lethal_embryoInd),function(emb){
	print(emb)
	estimate_phi(Env_combined_df_x[Env_combined_df_x$embryoID==emb,]$GD, Env_combined_df_x[Env_combined_df_x$embryoID==emb,]$D, Env_combined_df_x[Env_combined_df_x$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_Env_x,0))
mean(setdiff(population_Env_x,0))
median(setdiff(population_Env_x,0))
hist(population_Env_x)
mean(population_Env_x)
median(population_Env_x)



#y
Env_combined_df_y<-na.omit(do.call(rbind,lapply(setdiff(1929:2039,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

Env_combined_result_y<-estimate_phi(Env_combined_df_y$GD, Env_combined_df_y$D, Env_combined_df_y$M, phi_start = 0.618)
Env_combined_result_y[[1]]
Env_combined_result_y[[2]]


population_Env_y<-sapply(setdiff(1929:2039,lethal_embryoInd),function(emb){
	print(emb)
	estimate_phi(Env_combined_df_y[Env_combined_df_y$embryoID==emb,]$GD, Env_combined_df_y[Env_combined_df_y$embryoID==emb,]$D, Env_combined_df_y[Env_combined_df_y$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_Env_y,0))
mean(setdiff(population_Env_y,0))
median(setdiff(population_Env_y,0))
hist(population_Env_y)
mean(population_Env_y)
median(population_Env_y)



#z
Env_combined_df_z<-na.omit(do.call(rbind,lapply(setdiff(1929:2039,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

Env_combined_result_z<-estimate_phi(Env_combined_df_z$GD, Env_combined_df_z$D, Env_combined_df_z$M, phi_start = 0.618)
Env_combined_result_z[[1]]
Env_combined_result_z[[2]]


population_Env_z<-sapply(setdiff(1929:2039,lethal_embryoInd),function(emb){
	print(emb)
	estimate_phi(Env_combined_df_z[Env_combined_df_z$embryoID==emb,]$GD, Env_combined_df_z[Env_combined_df_z$embryoID==emb,]$D, Env_combined_df_z[Env_combined_df_z$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_Env_z,0))
mean(setdiff(population_Env_z,0))
median(setdiff(population_Env_z,0))
hist(population_Env_z)
mean(population_Env_z)
median(population_Env_z)


#xyz 
Env_combined_df_xyz<-rbind(Env_combined_df_x,Env_combined_df_y,Env_combined_df_z)
Env_combined_result_xyz<-estimate_phi(Env_combined_df_xyz$GD, Env_combined_df_xyz$D, Env_combined_df_xyz$M, phi_start = 0.618)
Env_combined_result_xyz[[1]]
Env_combined_result_xyz[[2]]


population_Env_xyz<-sapply(setdiff(1929:2039,lethal_embryoInd),function(emb){
	print(emb)
	estimate_phi(Env_combined_df_xyz[Env_combined_df_xyz$embryoID==emb,]$GD, Env_combined_df_xyz[Env_combined_df_xyz$embryoID==emb,]$D, Env_combined_df_xyz[Env_combined_df_xyz$embryoID==emb,]$M, phi_start = 0.618)[[1]]
})

hist(setdiff(population_Env_xyz,0))
mean(setdiff(population_Env_xyz,0))
median(setdiff(population_Env_xyz,0))
hist(population_Env_xyz)
mean(population_Env_xyz)
median(population_Env_xyz)


#
dataTemp<-data.frame(
coeff=c(RNAi_combined_result_x[[1]],RNAi_combined_result_y[[1]],RNAi_combined_result_z[[1]],
	Env_combined_result_x[[1]],Env_combined_result_y[[1]],Env_combined_result_z[[1]]),
se=c(RNAi_combined_result_x[[2]],RNAi_combined_result_y[[2]],RNAi_combined_result_z[[2]],
	Env_combined_result_x[[2]],Env_combined_result_y[[2]],Env_combined_result_z[[2]]),
coordinate=rep(c("x","y","z"),2),
perturbation=factor(c(rep("RNAi",3),rep("Env",3)),level=c("RNAi","Env")))

ggplot(dataTemp,aes(x=coordinate,y=coeff,ymin=coeff-1.96*se,ymax=coeff+1.96*se,col=coordinate))+
geom_point(size=0.1)+
geom_errorbar(width=0.2,lwd=0.1)+
geom_hline(yintercept=0.618,col="grey",lty=2,lwd=0.1)+
scale_y_continuous(limits=c(0,1),breaks=c(0,0.618,1))+
theme(legend.position=c(0.8,0.3))+
labs(x="Coordinate",y="Coefficient values (NLS)",title="Robustness to perturbation")+
facet_wrap(~perturbation)+
theme(strip.placement="outside", strip.background=element_blank())


ggsave("coefficient_estimation_by_NLS_for_RNAi_Env.pdf",width=unit(2.4,"in"),height=unit(1.7,"in"))


### noise injection 

library(ggplot2)

# Phi values
phi_low    <- 0.4
phi_golden <- (sqrt(5) - 1) / 2
phi_high   <- 0.8

phis   <- c(phi_low, phi_golden, phi_high)
labels <- c("Low Φ", "Golden Φ", "High Φ")

n_gen <- 40

# generations at which noise is injected
noise_gens <- c(0, 10, 20, 30)

# noise amplitudes (can be fixed or random)
set.seed(3)
noise_values <- rnorm(length(noise_gens), mean = 0, sd = 0.1)

# residual kernel for a single impulse
feedback_kernel <- function(phi, noise, dt) {
  if (dt == 0) return(noise)
  noise * (1 - sum(phi^(2:(dt + 1))))
}

# simulation
df1 <- expand.grid(
  generation = 0:n_gen,
  phi = phis
)

df1$residual <- NA

for (i in seq_len(nrow(df1))) {

  t   <- df1$generation[i]
  phi <- df1$phi[i]

  res <- 0

  for (k in seq_along(noise_gens)) {
    tk <- noise_gens[k]
    if (tk <= t) {
      res <- res + feedback_kernel(
        phi,
        noise_values[k],
        t - tk
      )
    }
  }

  df1$residual[i] <- res
}

df1$phi_label <- factor(df1$phi, levels = phis, labels = labels)

ggplot(df1, aes(generation, residual, color = phi_label)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.3) +
  geom_point(size = 0.5) +
  geom_vline(
    xintercept = noise_gens,
    linetype = "dotted",
    color = "grey50"
  ) +
  labs(
    x = "Generation",
    y = "Velocity residual",
    color = "Φ",
    title = "Noise injection"
  ) 


ggsave("simulation_for_golden_ratio_smaller_larger_noise_injection.pdf",width=unit(3.4,"in"),height=unit(1.7,"in"))


### simulation for generation-skipping 
library(ggplot2)
library(grid)
library(gridExtra)

### -------------------------
### Data generation
### -------------------------
set.seed(123)

n <- 1000
A <- rnorm(n)
B <- rnorm(n)
C <- rnorm(n)

X <- A - B
Y <- B - C
Z <- C - A

df <- data.frame(X = X, Y = Y, Z = Z)

### -------------------------
### Panel 1: correlation triangle
### -------------------------

nodes <- data.frame(
  x = c(0.2, 0.8, 0.5),
  y = c(0.2, 0.2, 0.8),
  label = c("X", "Y", "Z")
)

edges <- data.frame(
  x    = c(0.2, 0.8, 0.5),
  y    = c(0.2, 0.2, 0.8),
  xend = c(0.8, 0.5, 0.2),
  yend = c(0.2, 0.8, 0.2)
)

edge_labels <- data.frame(
  x = c(0.5, 0.67, 0.33),
  y = c(0.15, 0.55, 0.55),
  label = c("corr(X,Y) < 0",
            "corr(Y,Z) < 0",
            "corr(X,Z) < 0")
)

p_triangle <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.8
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y),
    size = 6,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = label),
    size = 4
  ) +
  geom_text(
    data = edge_labels,
    aes(x = x, y = y, label = label),
    size = 2
  ) +
  coord_fixed() +
  theme_void() +
  labs(title = "Pairwise negative correlations")

### -------------------------
### Panels 2–4: scatter plots
### -------------------------

p_xy <- ggplot(df, aes(X, Y)) +
  geom_point(alpha = 0.35, size = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "X vs Y", x = "X", y = "Y")

p_yz <- ggplot(df, aes(Y, Z)) +
  geom_point(alpha = 0.35, size = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Y vs Z", x = "Y", y = "Z")

p_xz <- ggplot(df, aes(X, Z)) +
  geom_point(alpha = 0.35, size = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "X vs Z (generation-skipping)", x = "X", y = "Z")

### -------------------------
### Output to PDF
### -------------------------

pdf("pairwise_negative_correlation_ggplot.pdf", width = 3, height = 3)

grid.arrange(
  p_triangle, p_xy,
  p_yz,       p_xz,
  ncol = 2
)

dev.off()


###

dataTemp<-data.frame(
	embryoPhi=c(population_control_xyz,population_RNAi_xyz,population_Env_xyz),
	condition=c(rep("Control",length(population_control_xyz)),rep("RNAi",length(population_RNAi_xyz)),rep("Env",length(population_Env_xyz)))
)

ggplot(dataTemp,aes(x=embryoPhi,col=condition))+
geom_density(lwd=0.5)+
geom_vline(xintercept=(sqrt(5)-1)/2,lwd=0.5,lty=2)
#geom_vline(xintercept=c(mean(population_control_xyz),mean(population_RNAi_xyz),mean(population_Env_xyz)),col=c("red","blue","green"))


ggsave("embryo_wise_phi_estimation_distribution.pdf",width=unit(2.0,"in"),height=unit(1.7,"in"))

meta_phi <- function(phi, se) {
  
  if(length(phi) != length(se)) {
    stop("phi and se must have the same length")
  }
  
  k <- length(phi)
  
  # Fixed-effect weights
  w_fixed <- 1 / (se^2)
  
  # Fixed-effect estimate
  phi_fixed <- sum(w_fixed * phi) / sum(w_fixed)
  se_fixed  <- sqrt(1 / sum(w_fixed))
  ci_fixed  <- c(phi_fixed - 1.96 * se_fixed,
                 phi_fixed + 1.96 * se_fixed)
  
  # Cochran's Q
  Q <- sum(w_fixed * (phi - phi_fixed)^2)
  
  # DerSimonian-Laird tau^2
  df <- k - 1
  C  <- sum(w_fixed) - sum(w_fixed^2) / sum(w_fixed)
  tau2 <- max(0, (Q - df) / C)
  
  # Random-effects weights
  w_random <- 1 / (se^2 + tau2)
  
  # Random-effects estimate
  phi_random <- sum(w_random * phi) / sum(w_random)
  se_random  <- sqrt(1 / sum(w_random))
  ci_random  <- c(phi_random - 1.96 * se_random,
                  phi_random + 1.96 * se_random)
  
  # I^2
  I2 <- max(0, (Q - df) / Q)
  
  return(list(
    fixed = list(
      estimate = phi_fixed,
      se = se_fixed,
      CI95 = ci_fixed
    ),
    random = list(
      estimate = phi_random,
      se = se_random,
      CI95 = ci_random,
      tau2 = tau2
    ),
    heterogeneity = list(
      Q = Q,
      df = df,
      I2 = I2
    )
  ))
}


# meta for control 
population_control_xyz_SE<-data.frame(t(sapply(setdiff(1:105,lethal_embryoInd),function(emb){
	print(emb)
	temp<-estimate_phi(combined_df_xyz[combined_df_xyz$embryoID==emb,]$GD, combined_df_xyz[combined_df_xyz$embryoID==emb,]$D, combined_df_xyz[combined_df_xyz$embryoID==emb,]$M, phi_start = 0.618)
	return(c(temp[[1]],temp[[2]]))
})))
names(population_control_xyz_SE)<-c("phi","SE")

meta_phi_control<-meta_phi(population_control_xyz_SE[,1], population_control_xyz_SE[,2])


# meta for RNAi 
population_RNAi_xyz_SE<-data.frame(t(sapply(setdiff(106:1928,lethal_embryoInd),function(emb){
	print(emb)
	temp<-estimate_phi(RNAi_combined_df_xyz[RNAi_combined_df_xyz$embryoID==emb,]$GD, RNAi_combined_df_xyz[RNAi_combined_df_xyz$embryoID==emb,]$D, RNAi_combined_df_xyz[RNAi_combined_df_xyz$embryoID==emb,]$M, phi_start = 0.618)
	return(c(temp[[1]],temp[[2]]))
})))
names(population_RNAi_xyz_SE)<-c("phi","SE")

meta_phi_RNAi<-meta_phi(population_RNAi_xyz_SE[,1], population_RNAi_xyz_SE[,2])


# meta for Env 
population_Env_xyz_SE<-data.frame(t(sapply(setdiff(1929:2039,lethal_embryoInd),function(emb){
	print(emb)
	temp<-estimate_phi(Env_combined_df_xyz[Env_combined_df_xyz$embryoID==emb,]$GD, Env_combined_df_xyz[Env_combined_df_xyz$embryoID==emb,]$D, Env_combined_df_xyz[Env_combined_df_xyz$embryoID==emb,]$M, phi_start = 0.618)
	return(c(temp[[1]],temp[[2]]))
})))
names(population_Env_xyz_SE)<-c("phi","SE")

meta_phi_Env<-meta_phi(population_Env_xyz_SE[,1], population_Env_xyz_SE[,2])


dataTemp<-data.frame(
	coeff=c(meta_phi_control$random$estimate,meta_phi_RNAi$random$estimate,meta_phi_Env$random$estimate),
	se=c(meta_phi_control$random$se,meta_phi_RNAi$random$se,meta_phi_Env$random$se),
	condition=c("Control","RNAi","Env")
)


ggplot(dataTemp,aes(x=condition,y=coeff,ymin=coeff-1.96*se,ymax=coeff+1.96*se,col=condition))+
geom_point(size=0.1)+
geom_errorbar(width=0.2,lwd=0.1)+
geom_hline(yintercept=0.618,col="grey",lty=2,lwd=0.1)+
scale_y_continuous(limits=c(0,1),breaks=c(0,0.618,1))+
theme(legend.position=c(0.8,0.2))+
labs(x="Condition",y="Parameter estimation",title="Across-embryo meta-analysis")+
guides(color = guide_legend(
    override.aes = list(size = 2,alpha=1),        # Increase key size (points in legend)
    keywidth = 0.5,                       # Reduce width of key box
    keyheight = 0.5                       # Reduce height of key box
  )) +
  theme(
	legend.background=element_blank(),
	legend.frame=element_blank(),
	legend.box.background=element_blank(),
    legend.spacing.x = unit(0.1, 'cm'),   # Reduce horizontal spacing between keys
    legend.spacing.y = unit(0.1, 'cm'),   # Reduce vertical spacing between keys
    legend.key.width = unit(0.3, 'cm'),   # Reduce spacing between key and text
    legend.key.height = unit(0.3, 'cm')   # Reduce height of key background
  )

ggsave("embryo_wise_phi_estimation_meta.pdf",width=unit(2.0,"in"),height=unit(1.7,"in"))


### without considering sd 
## x 
combined_df_x_withoutSD<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

combined_result_x_withoutSD<-estimate_phi(combined_df_x_withoutSD$GD, combined_df_x_withoutSD$D, combined_df_x_withoutSD$M, phi_start = 0.618)
combined_result_x_withoutSD[[1]]
combined_result_x_withoutSD[[2]]

## y 
combined_df_y_withoutSD<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

combined_result_y_withoutSD<-estimate_phi(combined_df_y_withoutSD$GD, combined_df_y_withoutSD$D, combined_df_y_withoutSD$M, phi_start = 0.618)
combined_result_y_withoutSD[[1]]
combined_result_y_withoutSD[[2]]

## z 
combined_df_z_withoutSD<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,]),
		embryoID=rep(emb,length(c(lineage_velocity_emb[1:8,])))
	)
	return(M_D_GD_df)
})))

combined_result_z_withoutSD<-estimate_phi(combined_df_z_withoutSD$GD, combined_df_z_withoutSD$D, combined_df_z_withoutSD$M, phi_start = 0.618)
combined_result_z_withoutSD[[1]]
combined_result_z_withoutSD[[2]]


focalMD<-na.omit(data.frame(
	MIndex=match(paste(cell_mo_da1_da2[,1],"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling)),
	DIndex=match(paste(cell_mo_da1_da2[,2],"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling)),
	generation=match(cell_mo_da1_da2[,1],lineage_info)%%10
))
sd_t1_2_sd_t<-apply(displacement[1:105,focalMD[,2]],2,sd)/apply(displacement[1:105,focalMD[,1]],2,sd)
hist(sd_t1_2_sd_t)
mean(sd_t1_2_sd_t)
plot(focalMD$generation,sd_t1_2_sd_t)

cor_MD<-sapply(1:dim(focalMD)[1],function(i){cor(start_pos[,focalMD[i,1]],start_pos[,focalMD[i,2]])})
plot(focalMD$generation,cor_MD)


### 4-generation or 5-generation 
## control 
## x4
estimate_phi4 <- function(GD, D, M, GM, phi_start = 0.6) {
  df <- data.frame(GD = GD, D = D, M = M)

  fit <- nls(
    GD ~ -phi^2 * D - phi^3 * M - phi^4 * GM,
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


## x4
combined_df_x4<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		GM=c(lineage_velocity_emb[1:7,]),
		M=c(lineage_velocity_emb[2:8,]),
		D=c(lineage_velocity_emb[3:9,]),
		GD=c(lineage_velocity_emb[4:10,])
	)
	return(M_D_GD_df)
})))

combined_result_x4<-estimate_phi4(combined_df_x4$GD, combined_df_x4$D, combined_df_x4$M, combined_df_x4$GM, phi_start = 0.618)
combined_result_x4[[1]]
combined_result_x4[[2]]

combined_result_x4_x3<-estimate_phi(combined_df_x4$GD, combined_df_x4$D, combined_df_x4$M, phi_start = 0.618)
combined_result_x4_x3[[1]]
combined_result_x4_x3[[2]]

plot(combined_result_x4[[3]],combined_result_x4_x3[[3]])
abline(0,1,"red")
cor(combined_result_x4[[3]],combined_result_x4_x3[[3]])^2


combined_df_y4<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		GM=c(lineage_velocity_emb[1:7,]),
		M=c(lineage_velocity_emb[2:8,]),
		D=c(lineage_velocity_emb[3:9,]),
		GD=c(lineage_velocity_emb[4:10,])
	)
	return(M_D_GD_df)
})))

combined_result_y4<-estimate_phi4(combined_df_y4$GD, combined_df_y4$D, combined_df_y4$M, combined_df_y4$GM, phi_start = 0.618)
combined_result_y4[[1]]
combined_result_y4[[2]]

combined_result_y4_y3<-estimate_phi(combined_df_y4$GD, combined_df_y4$D, combined_df_y4$M, phi_start = 0.618)
combined_result_y4_y3[[1]]
combined_result_y4_y3[[2]]

plot(combined_result_y4[[3]],combined_result_y4_y3[[3]])
abline(0,1,"red")
cor(combined_result_y4[[3]],combined_result_y4_y3[[3]])^2


combined_df_z4<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		GM=c(lineage_velocity_emb[1:7,]),
		M=c(lineage_velocity_emb[2:8,]),
		D=c(lineage_velocity_emb[3:9,]),
		GD=c(lineage_velocity_emb[4:10,])
	)
	return(M_D_GD_df)
})))

combined_result_z4<-estimate_phi4(combined_df_z4$GD, combined_df_z4$D, combined_df_z4$M, combined_df_z4$GM, phi_start = 0.618)
combined_result_z4[[1]]
combined_result_z4[[2]]

combined_result_z4_z3<-estimate_phi(combined_df_z4$GD, combined_df_z4$D, combined_df_z4$M, phi_start = 0.618)
combined_result_z4_z3[[1]]
combined_result_z4_z3[[2]]

plot(combined_result_z4[[3]],combined_result_z4_z3[[3]])
abline(0,1,"red")
cor(combined_result_z4[[3]],combined_result_z4_z3[[3]])^2

combined_df_xyz4<-rbind(combined_df_x4,combined_df_y4,combined_df_z4)
combined_result_xyz4<-estimate_phi4(combined_df_xyz4$GD, combined_df_xyz4$D, combined_df_xyz4$M, combined_df_xyz4$GM, phi_start = 0.618)
combined_result_xyz4[[1]]
combined_result_xyz4[[2]]


## x5

estimate_phi5 <- function(GD, D, M, GM, GGM, phi_start = 0.6) {
  df <- data.frame(GD = GD, D = D, M = M)

  fit <- nls(
    GD ~ -phi^2 * D - phi^3 * M - phi^4 * GM - phi^5 * GM,
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


combined_df_x5<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		GGM=c(lineage_velocity_emb[1:6,]),
		GM=c(lineage_velocity_emb[2:7,]),
		M=c(lineage_velocity_emb[3:8,]),
		D=c(lineage_velocity_emb[4:9,]),
		GD=c(lineage_velocity_emb[5:10,])
	)
	return(M_D_GD_df)
})))

combined_result_x5<-estimate_phi5(combined_df_x5$GD, combined_df_x5$D, combined_df_x5$M, combined_df_x5$GM, combined_df_x5$GGM, phi_start = 0.618)
combined_result_x5[[1]]
combined_result_x5[[2]]

combined_df_y5<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		GGM=c(lineage_velocity_emb[1:6,]),
		GM=c(lineage_velocity_emb[2:7,]),
		M=c(lineage_velocity_emb[3:8,]),
		D=c(lineage_velocity_emb[4:9,]),
		GD=c(lineage_velocity_emb[5:10,])
	)
	return(M_D_GD_df)
})))

combined_result_y5<-estimate_phi5(combined_df_y5$GD, combined_df_y5$D, combined_df_y5$M, combined_df_y5$GM, combined_df_y5$GGM, phi_start = 0.618)
combined_result_y5[[1]]
combined_result_y5[[2]]

combined_df_z5<-na.omit(do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		GGM=c(lineage_velocity_emb[1:6,]),
		GM=c(lineage_velocity_emb[2:7,]),
		M=c(lineage_velocity_emb[3:8,]),
		D=c(lineage_velocity_emb[4:9,]),
		GD=c(lineage_velocity_emb[5:10,])
	)
	return(M_D_GD_df)
})))

combined_result_z5<-estimate_phi5(combined_df_z5$GD, combined_df_z5$D, combined_df_z5$M, combined_df_z5$GM, combined_df_z5$GGM, phi_start = 0.618)
combined_result_z5[[1]]
combined_result_z5[[2]]

combined_df_xyz5<-rbind(combined_df_x5,combined_df_y5,combined_df_z5)
combined_result_xyz5<-estimate_phi5(combined_df_xyz5$GD, combined_df_xyz5$D, combined_df_xyz5$M, combined_df_xyz5$GM, combined_df_xyz5$GGM, phi_start = 0.618)
combined_result_xyz5[[1]]
combined_result_xyz5[[2]]

## RNAi 

## x4
RNAi_combined_df_x4<-na.omit(do.call(rbind,lapply(setdiff(106:1928,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		GM=c(lineage_velocity_emb[1:7,]),
		M=c(lineage_velocity_emb[2:8,]),
		D=c(lineage_velocity_emb[3:9,]),
		GD=c(lineage_velocity_emb[4:10,])
	)
	return(M_D_GD_df)
})))

RNAi_combined_result_x4<-estimate_phi4(RNAi_combined_df_x4$GD, RNAi_combined_df_x4$D, RNAi_combined_df_x4$M, RNAi_combined_df_x4$GM, phi_start = 0.618)
RNAi_combined_result_x4[[1]]
RNAi_combined_result_x4[[2]]


RNAi_combined_df_y4<-na.omit(do.call(rbind,lapply(setdiff(106:1928,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		GM=c(lineage_velocity_emb[1:7,]),
		M=c(lineage_velocity_emb[2:8,]),
		D=c(lineage_velocity_emb[3:9,]),
		GD=c(lineage_velocity_emb[4:10,])
	)
	return(M_D_GD_df)
})))

RNAi_combined_result_y4<-estimate_phi4(RNAi_combined_df_y4$GD, RNAi_combined_df_y4$D, RNAi_combined_df_y4$M, RNAi_combined_df_y4$GM, phi_start = 0.618)
RNAi_combined_result_y4[[1]]
RNAi_combined_result_y4[[2]]


RNAi_combined_df_z4<-na.omit(do.call(rbind,lapply(setdiff(106:1928,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		GM=c(lineage_velocity_emb[1:7,]),
		M=c(lineage_velocity_emb[2:8,]),
		D=c(lineage_velocity_emb[3:9,]),
		GD=c(lineage_velocity_emb[4:10,])
	)
	return(M_D_GD_df)
})))

RNAi_combined_result_z4<-estimate_phi4(RNAi_combined_df_z4$GD, RNAi_combined_df_z4$D, RNAi_combined_df_z4$M, RNAi_combined_df_z4$GM, phi_start = 0.618)
RNAi_combined_result_z4[[1]]
RNAi_combined_result_z4[[2]]


RNAi_combined_df_xyz4<-rbind(RNAi_combined_df_x4,RNAi_combined_df_y4,RNAi_combined_df_z4)
RNAi_combined_result_xyz4<-estimate_phi4(RNAi_combined_df_xyz4$GD, RNAi_combined_df_xyz4$D, RNAi_combined_df_xyz4$M, RNAi_combined_df_xyz4$GM, phi_start = 0.618)
RNAi_combined_result_xyz4[[1]]
RNAi_combined_result_xyz4[[2]]


## x5

RNAi_combined_df_x5<-na.omit(do.call(rbind,lapply(setdiff(106:1928,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		GGM=c(lineage_velocity_emb[1:6,]),
		GM=c(lineage_velocity_emb[2:7,]),
		M=c(lineage_velocity_emb[3:8,]),
		D=c(lineage_velocity_emb[4:9,]),
		GD=c(lineage_velocity_emb[5:10,])
	)
	return(M_D_GD_df)
})))

RNAi_combined_result_x5<-estimate_phi5(RNAi_combined_df_x5$GD, RNAi_combined_df_x5$D, RNAi_combined_df_x5$M, RNAi_combined_df_x5$GM, RNAi_combined_df_x5$GGM, phi_start = 0.618)
RNAi_combined_result_x5[[1]]
RNAi_combined_result_x5[[2]]

RNAi_combined_df_y5<-na.omit(do.call(rbind,lapply(setdiff(106:1928,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		GGM=c(lineage_velocity_emb[1:6,]),
		GM=c(lineage_velocity_emb[2:7,]),
		M=c(lineage_velocity_emb[3:8,]),
		D=c(lineage_velocity_emb[4:9,]),
		GD=c(lineage_velocity_emb[5:10,])
	)
	return(M_D_GD_df)
})))

RNAi_combined_result_y5<-estimate_phi5(RNAi_combined_df_y5$GD, RNAi_combined_df_y5$D, RNAi_combined_df_y5$M, RNAi_combined_df_y5$GM, RNAi_combined_df_y5$GGM, phi_start = 0.618)
RNAi_combined_result_y5[[1]]
RNAi_combined_result_y5[[2]]

RNAi_combined_df_z5<-na.omit(do.call(rbind,lapply(setdiff(106:1928,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		GGM=c(lineage_velocity_emb[1:6,]),
		GM=c(lineage_velocity_emb[2:7,]),
		M=c(lineage_velocity_emb[3:8,]),
		D=c(lineage_velocity_emb[4:9,]),
		GD=c(lineage_velocity_emb[5:10,])
	)
	return(M_D_GD_df)
})))

RNAi_combined_result_z5<-estimate_phi5(RNAi_combined_df_z5$GD, RNAi_combined_df_z5$D, RNAi_combined_df_z5$M, RNAi_combined_df_z5$GM, RNAi_combined_df_z5$GGM, phi_start = 0.618)
RNAi_combined_result_z5[[1]]
RNAi_combined_result_z5[[2]]

RNAi_combined_df_xyz5<-rbind(RNAi_combined_df_x5,RNAi_combined_df_y5,RNAi_combined_df_z5)
RNAi_combined_result_xyz5<-estimate_phi5(RNAi_combined_df_xyz5$GD, RNAi_combined_df_xyz5$D, RNAi_combined_df_xyz5$M, RNAi_combined_df_xyz5$GM, RNAi_combined_df_xyz5$GGM, phi_start = 0.618)
RNAi_combined_result_xyz5[[1]]
RNAi_combined_result_xyz5[[2]]


## Env

## x4
Env_combined_df_x4<-na.omit(do.call(rbind,lapply(setdiff(1929:2039,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		GM=c(lineage_velocity_emb[1:7,]),
		M=c(lineage_velocity_emb[2:8,]),
		D=c(lineage_velocity_emb[3:9,]),
		GD=c(lineage_velocity_emb[4:10,])
	)
	return(M_D_GD_df)
})))

Env_combined_result_x4<-estimate_phi4(Env_combined_df_x4$GD, Env_combined_df_x4$D, Env_combined_df_x4$M, Env_combined_df_x4$GM, phi_start = 0.618)
Env_combined_result_x4[[1]]
Env_combined_result_x4[[2]]


Env_combined_df_y4<-na.omit(do.call(rbind,lapply(setdiff(1929:2039,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		GM=c(lineage_velocity_emb[1:7,]),
		M=c(lineage_velocity_emb[2:8,]),
		D=c(lineage_velocity_emb[3:9,]),
		GD=c(lineage_velocity_emb[4:10,])
	)
	return(M_D_GD_df)
})))

Env_combined_result_y4<-estimate_phi4(Env_combined_df_y4$GD, Env_combined_df_y4$D, Env_combined_df_y4$M, Env_combined_df_y4$GM, phi_start = 0.618)
Env_combined_result_y4[[1]]
Env_combined_result_y4[[2]]


Env_combined_df_z4<-na.omit(do.call(rbind,lapply(setdiff(1929:2039,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		GM=c(lineage_velocity_emb[1:7,]),
		M=c(lineage_velocity_emb[2:8,]),
		D=c(lineage_velocity_emb[3:9,]),
		GD=c(lineage_velocity_emb[4:10,])
	)
	return(M_D_GD_df)
})))

Env_combined_result_z4<-estimate_phi4(Env_combined_df_z4$GD, Env_combined_df_z4$D, Env_combined_df_z4$M, Env_combined_df_z4$GM, phi_start = 0.618)
Env_combined_result_z4[[1]]
Env_combined_result_z4[[2]]


Env_combined_df_xyz4<-rbind(Env_combined_df_x4,Env_combined_df_y4,Env_combined_df_z4)
Env_combined_result_xyz4<-estimate_phi4(Env_combined_df_xyz4$GD, Env_combined_df_xyz4$D, Env_combined_df_xyz4$M, Env_combined_df_xyz4$GM, phi_start = 0.618)
Env_combined_result_xyz4[[1]]
Env_combined_result_xyz4[[2]]


## x5

Env_combined_df_x5<-na.omit(do.call(rbind,lapply(setdiff(1929:2039,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		GGM=c(lineage_velocity_emb[1:6,]),
		GM=c(lineage_velocity_emb[2:7,]),
		M=c(lineage_velocity_emb[3:8,]),
		D=c(lineage_velocity_emb[4:9,]),
		GD=c(lineage_velocity_emb[5:10,])
	)
	return(M_D_GD_df)
})))

Env_combined_result_x5<-estimate_phi5(Env_combined_df_x5$GD, Env_combined_df_x5$D, Env_combined_df_x5$M, Env_combined_df_x5$GM, Env_combined_df_x5$GGM, phi_start = 0.618)
Env_combined_result_x5[[1]]
Env_combined_result_x5[[2]]

Env_combined_df_y5<-na.omit(do.call(rbind,lapply(setdiff(1929:2039,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_D_GD_df<-data.frame(
		GGM=c(lineage_velocity_emb[1:6,]),
		GM=c(lineage_velocity_emb[2:7,]),
		M=c(lineage_velocity_emb[3:8,]),
		D=c(lineage_velocity_emb[4:9,]),
		GD=c(lineage_velocity_emb[5:10,])
	)
	return(M_D_GD_df)
})))

Env_combined_result_y5<-estimate_phi5(Env_combined_df_y5$GD, Env_combined_df_y5$D, Env_combined_df_y5$M, Env_combined_df_y5$GM, Env_combined_df_y5$GGM, phi_start = 0.618)
Env_combined_result_y5[[1]]
Env_combined_result_y5[[2]]

Env_combined_df_z5<-na.omit(do.call(rbind,lapply(setdiff(1929:2039,lethal_embryoInd),function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_D_GD_df<-data.frame(
		GGM=c(lineage_velocity_emb[1:6,]),
		GM=c(lineage_velocity_emb[2:7,]),
		M=c(lineage_velocity_emb[3:8,]),
		D=c(lineage_velocity_emb[4:9,]),
		GD=c(lineage_velocity_emb[5:10,])
	)
	return(M_D_GD_df)
})))

Env_combined_result_z5<-estimate_phi5(RNAi_combined_df_z5$GD, RNAi_combined_df_z5$D, RNAi_combined_df_z5$M, RNAi_combined_df_z5$GM, RNAi_combined_df_z5$GGM, phi_start = 0.618)
Env_combined_result_z5[[1]]
Env_combined_result_z5[[2]]

Env_combined_df_xyz5<-rbind(Env_combined_df_x5,Env_combined_df_y5,Env_combined_df_z5)
Env_combined_result_xyz5<-estimate_phi5(Env_combined_df_xyz5$GD, Env_combined_df_xyz5$D, Env_combined_df_xyz5$M, Env_combined_df_xyz5$GM, Env_combined_df_xyz5$GGM, phi_start = 0.618)
Env_combined_result_xyz5[[1]]
Env_combined_result_xyz5[[2]]


#
dataTemp<-data.frame(
coeff=c(
	combined_result_xyz[[1]],combined_result_xyz4[[1]],combined_result_xyz5[[1]],
	RNAi_combined_result_xyz[[1]],RNAi_combined_result_xyz4[[1]],RNAi_combined_result_xyz5[[1]],
	Env_combined_result_xyz[[1]],Env_combined_result_xyz4[[1]],Env_combined_result_xyz5[[1]]
	),
se=c(
	combined_result_xyz[[2]],combined_result_xyz4[[2]],combined_result_xyz5[[2]],
	RNAi_combined_result_xyz[[2]],RNAi_combined_result_xyz4[[2]],RNAi_combined_result_xyz5[[2]],
	Env_combined_result_xyz[[2]],Env_combined_result_xyz4[[2]],Env_combined_result_xyz5[[2]]
	),
generation=rep(c("3","4","5"),3),
Condition=factor(c(rep("control",3),rep("RNAi",3),rep("Env",3)),level=c("control","RNAi","Env")))

ggplot(dataTemp,aes(x=generation,y=coeff,ymin=coeff-1.96*se,ymax=coeff+1.96*se,col=generation))+
geom_point(size=0.1)+
geom_errorbar(width=0.2,lwd=0.1)+
geom_hline(yintercept=0.618,col="grey",lty=2,lwd=0.1)+
scale_y_continuous(limits=c(0,1),breaks=c(0,0.618,1))+
theme(legend.position=c(0.8,0.3))+
labs(x="Approximated generation",y="Coefficient values (NLS)",title="Robustness to approximated generation")+
facet_wrap(~Condition)+
theme(strip.placement="outside", strip.background=element_blank())


ggsave("coefficient_estimation_by_NLS_for_approximated_generation.pdf",width=unit(4.4,"in"),height=unit(2.1,"in"))











#####################################################################################
# END
#####################################################################################



# for x 
#lineage_velocity_x_avg<-matrix(apply(velocity[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
#lineage_velocity_x_sd<-matrix(apply(velocity[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,sd,na.rm=T),nrow=dim(lineage_info)[1])


# for x 
lineage_displacement_x_avg<-matrix(apply(displacement[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_displacement_x_sd<-matrix(apply(displacement[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,sd,na.rm=T),nrow=dim(lineage_info)[1])



# for x 
lineage_tim_x_avg<-matrix(apply(tim[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_tim_x_sd<-matrix(apply(tim[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,sd,na.rm=T),nrow=dim(lineage_info)[1])


par(mfrow=c(1,2))
plot(lineage_velocity_x_avg,lineage_displacement_x_avg)
plot(lineage_velocity_x_sd,lineage_displacement_x_sd)
cor(as.numeric(lineage_velocity_x_avg),as.numeric(lineage_displacement_x_avg),use="pairwise.complete.obs")
cor(as.numeric(lineage_velocity_x_sd),as.numeric(lineage_displacement_x_sd),use="pairwise.complete.obs")

hist(sapply(1:361,function(i){cor(lineage_velocity_x_sd[8:10,i],lineage_displacement_x_sd[8:10,i],use="pairwise.complete.obs")}),xlab=NA)

## x
combined_df_x_velocity<-do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_velocity_emb[1:8,]),
		D=c(lineage_velocity_emb[2:9,]),
		GD=c(lineage_velocity_emb[3:10,])
	)
	return(M_D_GD_df)
}))


combined_df_x_displacement<-do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_displacement_x_avg)/lineage_displacement_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_displacement_emb[1:8,]),
		D=c(lineage_displacement_emb[2:9,]),
		GD=c(lineage_displacement_emb[3:10,])
	)
	return(M_D_GD_df)
}))




plot(as.numeric(as.matrix(combined_df_x_velocity)), as.numeric(as.matrix(combined_df_x_displacement)))

# alignment between Zscore of velocity and Zscore of displacement 
combined_df_x_velocity<-do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	return(lineage_velocity_emb)
}))

combined_df_x_displacement<-do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_displacement_x_avg)/lineage_displacement_x_sd
	return(lineage_displacement_emb)
}))

#y
combined_df_y_velocity<-do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	return(lineage_velocity_emb)
}))

combined_df_y_displacement<-do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_displacement_y_avg)/lineage_displacement_y_sd
	return(lineage_displacement_emb)
}))

#z
combined_df_z_velocity<-do.call(rbind,lapply(1:105,function(emb){

	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	return(lineage_velocity_emb)
}))

combined_df_z_displacement<-do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_displacement_z_avg)/lineage_displacement_z_sd
	return(lineage_displacement_emb)
}))



combined_df_x_tim<-do.call(rbind,lapply(1:105,function(emb){

	lineage_tim_emb<-(matrix(tim[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_tim_x_avg)/lineage_tim_x_sd
	return(lineage_tim_emb)
}))


### comparison between Z-scores of velocity and displacement 
# Zscore(disp) = (disp-E(disp))/SD(disp)
# Zscore(vel) = (vel-E(vel))/SD(vel)

par(mfrow=c(3,3))
plot(lineage_displacement_x_avg,lineage_velocity_x_avg)
plot(lineage_displacement_y_avg,lineage_velocity_y_avg)
plot(lineage_displacement_z_avg,lineage_velocity_z_avg)

plot(lineage_displacement_x_sd,lineage_velocity_x_sd)
plot(lineage_displacement_y_sd,lineage_velocity_y_sd)
plot(lineage_displacement_z_sd,lineage_velocity_z_sd)

plot(lineage_displacement_x_avg/lineage_displacement_x_sd,lineage_velocity_x_avg/lineage_velocity_x_sd)
plot(lineage_displacement_y_avg/lineage_displacement_y_sd,lineage_velocity_y_avg/lineage_velocity_y_sd)
plot(lineage_displacement_z_avg/lineage_displacement_z_sd,lineage_velocity_z_avg/lineage_velocity_z_sd)



dataTemp<-data.frame(
	Zscore_displacement=c(combined_df_x_displacement,combined_df_y_displacement,combined_df_z_displacement),
	Zscore_velocity=c(combined_df_x_velocity,combined_df_y_velocity,combined_df_z_velocity),
	coordinate=c(rep("x",length(combined_df_x_displacement)),rep("y",length(combined_df_y_displacement)),rep("z",length(combined_df_z_displacement)))
)

ggplot(dataTemp,aes(x=Zscore_velocity,y=Zscore_displacement,col=coordinate))+
geom_point(size=0.15,stroke=0.1,alpha=0.3)+
geom_abline(slope=1,intercept=0,lwd=0.3)+
labs(x="Z-score of velocity",y="Z-score of displacement",title="")+
theme(legend.position=c(0.8,0.2))+
guides(color = guide_legend(
    override.aes = list(size = 2,alpha=1),        # Increase key size (points in legend)
    keywidth = 0.5,                       # Reduce width of key box
    keyheight = 0.5                       # Reduce height of key box
  )) +
  theme(
	legend.background=element_blank(),
	legend.frame=element_blank(),
	legend.box.background=element_blank(),
    legend.spacing.x = unit(0.1, 'cm'),   # Reduce horizontal spacing between keys
    legend.spacing.y = unit(0.1, 'cm'),   # Reduce vertical spacing between keys
    legend.key.width = unit(0.3, 'cm'),   # Reduce spacing between key and text
    legend.key.height = unit(0.3, 'cm')   # Reduce height of key background
  )

ggsave("Zscore_comparision_velocity_Vs_displacement.pdf",width=unit(2,"in"),height=unit(2.1,"in"))




### self-similarity equation for displacement 
# x 
combined_df_x_displacement<-do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_displacement_x_avg)/lineage_displacement_x_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_displacement_emb[1:8,]),
		D=c(lineage_displacement_emb[2:9,]),
		GD=c(lineage_displacement_emb[3:10,])
	)
	return(M_D_GD_df)
}))


combined_result_x_displacement<-estimate_phi(combined_df_x_displacement$GD, combined_df_x_displacement$D, combined_df_x_displacement$M, phi_start = 0.618)
combined_result_x_displacement[[1]]
combined_result_x_displacement[[2]]

# y 
combined_df_y_displacement<-do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_displacement_y_avg)/lineage_displacement_y_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_displacement_emb[1:8,]),
		D=c(lineage_displacement_emb[2:9,]),
		GD=c(lineage_displacement_emb[3:10,])
	)
	return(M_D_GD_df)
}))


combined_result_y_displacement<-estimate_phi(combined_df_y_displacement$GD, combined_df_y_displacement$D, combined_df_y_displacement$M, phi_start = 0.618)
combined_result_y_displacement[[1]]
combined_result_y_displacement[[2]]

# z 
combined_df_z_displacement<-do.call(rbind,lapply(1:105,function(emb){

	lineage_displacement_emb<-(matrix(displacement[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_displacement_z_avg)/lineage_displacement_z_sd
	M_D_GD_df<-data.frame(
		M=c(lineage_displacement_emb[1:8,]),
		D=c(lineage_displacement_emb[2:9,]),
		GD=c(lineage_displacement_emb[3:10,])
	)
	return(M_D_GD_df)
}))


combined_result_z_displacement<-estimate_phi(combined_df_z_displacement$GD, combined_df_z_displacement$D, combined_df_z_displacement$M, phi_start = 0.618)
combined_result_z_displacement[[1]]
combined_result_z_displacement[[2]]


# nls for displacement

dataTemp<-data.frame(
coeff=c(combined_result_x_displacement[[1]],combined_result_y_displacement[[1]],combined_result_z_displacement[[1]]),
se=c(combined_result_x_displacement[[2]],combined_result_y_displacement[[2]],combined_result_z_displacement[[2]]),
coordinate=c("x","y","z"))

ggplot(dataTemp,aes(x=coordinate,y=coeff,ymin=coeff-1.96*se,ymax=coeff+1.96*se,col=coordinate))+
geom_point(size=0.1)+
geom_errorbar(width=0.2,lwd=0.1)+
geom_hline(yintercept=0.618,col="grey",lty=2,lwd=0.1)+
scale_y_continuous(limits=c(0,1),breaks=c(0,0.618,1))+
theme(legend.position=c(0.8,0.2))+
labs(x="Coordinate",y="Parameter estimation",title="NLS")+
guides(color = guide_legend(
    override.aes = list(size = 2,alpha=1),        # Increase key size (points in legend)
    keywidth = 0.5,                       # Reduce width of key box
    keyheight = 0.5                       # Reduce height of key box
  )) +
  theme(
	legend.background=element_blank(),
	legend.frame=element_blank(),
	legend.box.background=element_blank(),
    legend.spacing.x = unit(0.1, 'cm'),   # Reduce horizontal spacing between keys
    legend.spacing.y = unit(0.1, 'cm'),   # Reduce vertical spacing between keys
    legend.key.width = unit(0.3, 'cm'),   # Reduce spacing between key and text
    legend.key.height = unit(0.3, 'cm')   # Reduce height of key background
  )


ggsave("coefficient_estimation_by_NLS_for_displacement.pdf",width=unit(2.2,"in"),height=unit(2.1,"in"))



