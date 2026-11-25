

### load Rscripts and RData 

loadResources<-function(resourceVector,projectFolderPath=focal_projectFolderPath){
	
	library(tools)
	
	projectFolderPath<-gsub("\\\\", "/", projectFolderPath)
	projectFolder<-basename(projectFolderPath)
	
	resourceVector<-sapply(resourceVector,function(rs){
		rs_folders<-unlist(strsplit(rs,"/"))
		return(gsub("//","/",paste0(projectFolderPath,"/",paste0(rs_folders[(which(rs_folders==projectFolder)+1):length(rs_folders)],collapse="/"))))
	})
	
	sapply(1:length(resourceVector),function(i){
		if(file_ext(resourceVector[i])=="R"){
			source(resourceVector[i],chdir=T)
			return(resourceVector[i])
		}
		if(file_ext(resourceVector[i])=="RData"){
			load(resourceVector[i],verbose=T,envir=globalenv())
			return(resourceVector[i])
		}

	})


}



jackknife_cor_parallel <- function(x, y, n_cores = NULL) {

  n <- length(x)
  
  # Full sample correlation
  r_full <- cor(x, y, method = "pearson", use = "pairwise.complete.obs")
  
  # Create cluster
  cl <- makeCluster(n_cores)
  
  # Run jackknife in parallel
  r_jack <- foreach(i=1:n, packages=c("parallel","doParallel","foreach"),.combine=c)%dopar%{
    x_jack <- x[-i]
    y_jack <- y[-i]
    return(cor(x_jack, y_jack, method = "pearson", use = "pairwise.complete.obs"))
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Calculate jackknife statistics
  pseudovalues <- n * r_full - (n - 1) * r_jack
  r_jack_est <- mean(pseudovalues)
  se_jack <- sd(pseudovalues) / sqrt(n)
  bias_jack <- (n - 1) * (mean(r_jack) - r_full)
  ci_lower <- r_jack_est - 1.96 * se_jack
  ci_upper <- r_jack_est + 1.96 * se_jack
  
  return(list(
    full_cor = r_full,
    jackknife_est = r_jack_est,
    jackknife_se = se_jack,
    bias = bias_jack,
    ci_lower=ci_lower, 
	ci_upper=ci_upper)
  )
}


bootstrap_cor_parallel <- function(x, y, n_boot = 1000, n_cores = NULL, seed = 123) {
  n <- length(x)
  
  # Full sample correlation
  r_full <- cor(x, y, method = "pearson", use = "pairwise.complete.obs")
  
  # Set up parallel computing
  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Generate bootstrap indices
  boot_indices <- replicate(n_boot, sample(1:n, n, replace = TRUE))
  
  # Register parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export data to workers
  clusterExport(cl, c("x", "y"), envir = environment())
  
  # Bootstrap in parallel
  r_boot <- foreach(i = 1:n_boot, 
                    .combine = c,
                    .packages = c("foreach")) %dopar% {
    idx <- boot_indices[, i]
    x_boot <- x[idx]
    y_boot <- y[idx]
    cor(x_boot, y_boot, method = "pearson", use = "pairwise.complete.obs")
  }
  
  # Stop cluster
  stopCluster(cl)
  registerDoSEQ()
  
  # Calculate bootstrap statistics
  # Bootstrap estimate (mean of bootstrap replicates)
  r_boot_est <- mean(r_boot)
  
  # Bias correction (bootstrap bias)
  bias_boot <- mean(r_boot) - r_full
  
  # Standard error (standard deviation of bootstrap estimates)
  se_boot <- sd(r_boot)
  
  # Normal confidence interval
  ci_lower <- r_full - 1.96 * se_boot
  ci_upper <- r_full + 1.96 * se_boot
  
  # Alternative: Percentile CI (more commonly used in bootstrap)
  ci_perc <- quantile(r_boot, probs = c(0.025, 0.975))
  
  # Alternative: Bias-corrected confidence interval
  z0 <- qnorm(mean(r_boot < r_full))  # bias correction factor
  ci_bc <- quantile(r_boot, probs = pnorm(2*z0 + c(-1.96, 1.96)))
  
  return(list(
    full_cor = r_full,
    boot_est = r_boot_est,           # Corresponds to jackknife_est
    boot_se = se_boot,               # Corresponds to jackknife_se
    bias = bias_boot,                # Corresponds to bias_jack
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    ci_percentile = ci_perc,
    ci_bias_corrected = ci_bc,
    n_boot = n_boot,
    boot_replicates = r_boot
  ))
}

#######################################################################

#######################################################################

library(ggplot2)
library(reshape)
library(gplots)
library(parallel)
library(doParallel)
library(foreach)

focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"E:/papers/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_R_shuffle.RData",	
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_pvalue.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_R_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_pvalue_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_R_shuffle.RData",	
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_pvalue.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_R_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_pvalue_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_R_shuffle.RData",	
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_pvalue.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_R_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_pvalue_shuffle.RData"	
	),projectFolderPath=focal_projectFolderPath)


#


### mother-sister paradox 

## velocity: cross-embryo
ccl_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling),"_ccl"))

temp1<-unlist(lapply(ccl_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(ccl_names,function(x){which(cell_mo_da1_da2$da1==x)}))

mo_sis_pair_ccl<-cell_mo_da1_da2[intersect(temp1,temp2),]

# plot 
R_MDD_velocity_across_xyz<-t(sapply(1:dim(mo_sis_pair_ccl)[1],function(i){
	M<-mo_sis_pair_ccl[i,1]
	D1<-mo_sis_pair_ccl[i,2]
	D2<-mo_sis_pair_ccl[i,3]
	
	M_index<-match(paste0(M,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D1_index<-match(paste0(D1,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D2_index<-match(paste0(D2,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	
	return(c(cor(velocity[,M_index],velocity[,D1_index],use="pairwise.complete.obs"),cor(velocity[,M_index],velocity[,D2_index],use="pairwise.complete.obs"),cor(velocity[,D1_index],velocity[,D2_index],use="pairwise.complete.obs"),
			 cor(velocity[,M_index+353],velocity[,D1_index+353],use="pairwise.complete.obs"),cor(velocity[,M_index+353],velocity[,D2_index+353],use="pairwise.complete.obs"),cor(velocity[,D1_index+353],velocity[,D2_index+353],use="pairwise.complete.obs"),
			 cor(velocity[,M_index+353*2],velocity[,D1_index+353*2],use="pairwise.complete.obs"),cor(velocity[,M_index+353*2],velocity[,D2_index+353*2],use="pairwise.complete.obs"),cor(velocity[,D1_index+353*2],velocity[,D2_index+353*2],use="pairwise.complete.obs")
			))
}))


dataTemp<-data.frame(
	R=c(R_MDD_velocity_across_xyz),
	coordinate=c(rep("x",dim(R_MDD_velocity_across_xyz)[1]*3),rep("y",dim(R_MDD_velocity_across_xyz)[1]*3),rep("z",dim(R_MDD_velocity_across_xyz)[1]*3)),
	pairGroup=rep(c(rep("MD1",dim(R_MDD_velocity_across_xyz)[1]),rep("MD2",dim(R_MDD_velocity_across_xyz)[1]),rep("D1D2",dim(R_MDD_velocity_across_xyz)[1])),3),
	pairName=c(paste(c(paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,2],sep="_"),paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,3],sep="_"),paste(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3],sep="_")),"x",sep="_"),
			   paste(c(paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,2],sep="_"),paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,3],sep="_"),paste(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3],sep="_")),"y",sep="_"),
			   paste(c(paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,2],sep="_"),paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,3],sep="_"),paste(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3],sep="_")),"z",sep="_")
			   )
)

dataTemp<-dataTemp[which(!duplicated(dataTemp$R)),]

dataTemp$pairGroup2<-ifelse(dataTemp$pairGroup=="D1D2","D1D2","MD")


ggplot(dataTemp,aes(x=pairGroup2,y=R,col=coordinate))+
geom_boxplot(width=0.6,size=0.2,outlier.size=0.1)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(x="Cell pair group",y="Pearson's R",title="Cross-embryo velocity")+
theme(legend.position=c(0.8,0.8))

ggsave("E:/canalizationENV/figs/mother_sister_paradox_cross_embryo_velocity_correlation.pdf",width=unit(1.6,"in"),height=unit(1.4,"in"))


## velocity: within-embryo
velocity_avg<-apply(velocity,2,mean,na.rm=T)
velocity_sd<-apply(velocity,2,sd,na.rm=T)
velocity_scaled<-t(apply(velocity,1,function(x){(x-velocity_avg)/velocity_sd}))

# plot 
R_MDD_velocity_within_xyz<-t(sapply(1:dim(velocity)[1],function(i){
	M_index<-match(paste0(mo_sis_pair_ccl[,1],"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D1_index<-match(paste0(mo_sis_pair_ccl[,2],"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D2_index<-match(paste0(mo_sis_pair_ccl[,3],"_ccl"),colnames(ccl_raw_nonNA_scaling))
	
	temp<-data.frame(
		M=rep(mo_sis_pair_ccl[,1],2),
		D=c(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3]),
		MD=paste(rep(mo_sis_pair_ccl[,1],2),c(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3]),sep="_")
	)
	temp<-temp[which(!duplicated(temp$MD)),]
	M_index_adj<-match(paste0(temp$M,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D_index_adj<-match(paste0(temp$D,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	
	return(c(cor(velocity_scaled[i,M_index_adj],velocity_scaled[i,D_index_adj],use="pairwise.complete.obs"),cor(velocity_scaled[i,D1_index],velocity_scaled[i,D2_index],use="pairwise.complete.obs"),
			 cor(velocity_scaled[i,M_index_adj+353],velocity_scaled[i,D_index_adj+353],use="pairwise.complete.obs"),cor(velocity_scaled[i,D1_index+353],velocity_scaled[i,D2_index+353],use="pairwise.complete.obs"),
			 cor(velocity_scaled[i,M_index_adj+353*2],velocity_scaled[i,D_index_adj+353*2],use="pairwise.complete.obs"),cor(velocity_scaled[i,D1_index+353*2],velocity_scaled[i,D2_index+353*2],use="pairwise.complete.obs")
			))
}))


dataTemp<-data.frame(
	R=c(R_MDD_velocity_within_xyz),
	coordinate=c(rep("x",dim(R_MDD_velocity_within_xyz)[1]*2),rep("y",dim(R_MDD_velocity_within_xyz)[1]*2),rep("z",dim(R_MDD_velocity_within_xyz)[1]*2)),
	pairName=rep(c(rep("MD",dim(R_MDD_velocity_within_xyz)[1]),rep("D1D2",dim(R_MDD_velocity_within_xyz)[1])),3)
)

ggplot(dataTemp,aes(x=pairName,y=R,col=coordinate))+
geom_boxplot(width=0.6,size=0.2,outlier.size=0.1)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(x="Cell pair group",y="Pearson's R",title="Within-embryo velocity")+
theme(legend.position=c(0.8,0.8))

ggsave("E:/canalizationENV/figs/mother_sister_paradox_within_embryo_velocity_correlation.pdf",width=unit(1.6,"in"),height=unit(1.4,"in"))


## negative-feedback k (jackknife)
# for M-D 
feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

k_cells<-do.call(rbind,strsplit(names(xyz_x_model_mother_coef),"_"))[,2]

kMD<-data.frame(
	kMx=rep(xyz_x_model_mother_coef[match(mo_sis_pair[,1],k_cells)],2),
	kDx=c(xyz_x_model_mother_coef[match(mo_sis_pair[,2],k_cells)],xyz_x_model_mother_coef[match(mo_sis_pair[,3],k_cells)]),
	kMy=rep(xyz_y_model_mother_coef[match(mo_sis_pair[,1],k_cells)],2),
	kDy=c(xyz_y_model_mother_coef[match(mo_sis_pair[,2],k_cells)],xyz_y_model_mother_coef[match(mo_sis_pair[,3],k_cells)]),
	kMz=rep(xyz_z_model_mother_coef[match(mo_sis_pair[,1],k_cells)],2),
	kDz=c(xyz_z_model_mother_coef[match(mo_sis_pair[,2],k_cells)],xyz_z_model_mother_coef[match(mo_sis_pair[,3],k_cells)]),
	D=c(rep(mo_sis_pair[,1],2),mo_sis_pair[,2],mo_sis_pair[,3])
)
kMD<-kMD[which(!duplicated(kMD$D)),]


jackknife_cor_x_MD<-jackknife_cor_parallel(kMD$kMx,kMD$kDx)
jackknife_cor_y_MD<-jackknife_cor_parallel(kMD$kMy,kMD$kDy)
jackknife_cor_z_MD<-jackknife_cor_parallel(kMD$kMz,kMD$kDz)

# for D1-D2 

kD1D2<-data.frame(
	kD1x=xyz_x_model_mother_coef[match(mo_sis_pair[,2],k_cells)],
	kD2x=xyz_x_model_mother_coef[match(mo_sis_pair[,3],k_cells)],
	kD1y=xyz_y_model_mother_coef[match(mo_sis_pair[,2],k_cells)],
	kD2y=xyz_y_model_mother_coef[match(mo_sis_pair[,3],k_cells)],
	kD1z=xyz_z_model_mother_coef[match(mo_sis_pair[,2],k_cells)],
	kD2z=xyz_z_model_mother_coef[match(mo_sis_pair[,3],k_cells)]
)


jackknife_cor_x_D1D2<-jackknife_cor_parallel(kD1D2[,1],kD1D2[,2])
jackknife_cor_y_D1D2<-jackknife_cor_parallel(kD1D2[,3],kD1D2[,4])
jackknife_cor_z_D1D2<-jackknife_cor_parallel(kD1D2[,5],kD1D2[,6])


# combine 

dataTemp<-data.frame(
	R=c(jackknife_cor_x_MD$jackknife_est,jackknife_cor_x_D1D2$jackknife_est,
		jackknife_cor_y_MD$jackknife_est,jackknife_cor_y_D1D2$jackknife_est,
		jackknife_cor_z_MD$jackknife_est,jackknife_cor_z_D1D2$jackknife_est),
	SE=c(jackknife_cor_x_MD$jackknife_se,jackknife_cor_x_D1D2$jackknife_se,
		jackknife_cor_y_MD$jackknife_se,jackknife_cor_y_D1D2$jackknife_se,
		jackknife_cor_z_MD$jackknife_se,jackknife_cor_z_D1D2$jackknife_se),		
	coordinate=c(rep("x",2),rep("y",2),rep("z",2)),
	pairName=rep(c("MD","D1D2"),3)
)


ggplot(dataTemp,aes(x=pairName,y=R,col=coordinate,fill=coordinate))+
geom_bar(width=0.6,size=0.2,stat = "identity", position = "dodge")+
labs(x="Cell pair name",y="Pearson's R",title="Negative-feedback strength (k)")+
theme(legend.position=c(0.8,0.8))


ggsave("E:/canalizationENV/figs/mother_sister_paradox_k_correlation_jackknife.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



## negative-feedback k (bootstrap)
# for M-D 

bootstrap_cor_x_MD<-bootstrap_cor_parallel(kMD$kMx,kMD$kDx)
bootstrap_cor_y_MD<-bootstrap_cor_parallel(kMD$kMy,kMD$kDy)
bootstrap_cor_z_MD<-bootstrap_cor_parallel(kMD$kMz,kMD$kDz)

# for D1-D2 

bootstrap_cor_x_D1D2<-bootstrap_cor_parallel(kD1D2[,1],kD1D2[,2])
bootstrap_cor_y_D1D2<-bootstrap_cor_parallel(kD1D2[,3],kD1D2[,4])
bootstrap_cor_z_D1D2<-bootstrap_cor_parallel(kD1D2[,5],kD1D2[,6])


# combine 

dataTemp<-data.frame(
	R=c(bootstrap_cor_x_MD$boot_est,bootstrap_cor_x_D1D2$boot_est,
		bootstrap_cor_y_MD$boot_est,bootstrap_cor_y_D1D2$boot_est,
		bootstrap_cor_z_MD$boot_est,bootstrap_cor_z_D1D2$boot_est),
	SE=c(bootstrap_cor_x_MD$boot_se,bootstrap_cor_x_D1D2$boot_se,
		bootstrap_cor_y_MD$boot_se,bootstrap_cor_y_D1D2$boot_se,
		bootstrap_cor_z_MD$boot_se,bootstrap_cor_z_D1D2$boot_se),		
	coordinate=c(rep("x",2),rep("y",2),rep("z",2)),
	pairName=rep(c("MD","D1D2"),3)
)


ggplot(dataTemp,aes(x=pairName,y=R,col=coordinate,fill=coordinate))+
geom_bar(width=0.6,size=0.2,stat = "identity", position = "dodge", alpha=0.3)+
geom_errorbar(aes(ymin = R - 1.96*SE, ymax = R + 1.96*SE), 
                width = 0.2, 
                position = position_dodge(width = 0.6),
                size = 0.3, color = "black") +
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x="Cell pair group",y="Pearson's R",title="Negative-feedback strength (k)")+
theme(legend.position=c(0.8,0.8))


ggsave("E:/canalizationENV/figs/mother_sister_paradox_k_correlation_bootstrap.pdf",width=unit(1.6,"in"),height=unit(1.4,"in"))

## displacement 

## displacement: cross-embryo
ccl_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling),"_ccl"))

temp1<-unlist(lapply(ccl_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(ccl_names,function(x){which(cell_mo_da1_da2$da1==x)}))

mo_sis_pair_ccl<-cell_mo_da1_da2[intersect(temp1,temp2),]

# plot 
R_MDD_displacement_across_xyz<-t(sapply(1:dim(mo_sis_pair_ccl)[1],function(i){
	M<-mo_sis_pair_ccl[i,1]
	D1<-mo_sis_pair_ccl[i,2]
	D2<-mo_sis_pair_ccl[i,3]
	
	M_index<-match(paste0(M,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D1_index<-match(paste0(D1,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D2_index<-match(paste0(D2,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	
	return(c(cor(displacement[,M_index],displacement[,D1_index],use="pairwise.complete.obs"),cor(displacement[,M_index],displacement[,D2_index],use="pairwise.complete.obs"),cor(displacement[,D1_index],displacement[,D2_index],use="pairwise.complete.obs"),
			 cor(displacement[,M_index+353],displacement[,D1_index+353],use="pairwise.complete.obs"),cor(displacement[,M_index+353],displacement[,D2_index+353],use="pairwise.complete.obs"),cor(displacement[,D1_index+353],displacement[,D2_index+353],use="pairwise.complete.obs"),
			 cor(displacement[,M_index+353*2],displacement[,D1_index+353*2],use="pairwise.complete.obs"),cor(displacement[,M_index+353*2],displacement[,D2_index+353*2],use="pairwise.complete.obs"),cor(displacement[,D1_index+353*2],displacement[,D2_index+353*2],use="pairwise.complete.obs")
			))
}))


dataTemp<-data.frame(
	R=c(R_MDD_displacement_across_xyz),
	coordinate=c(rep("x",dim(R_MDD_displacement_across_xyz)[1]*3),rep("y",dim(R_MDD_displacement_across_xyz)[1]*3),rep("z",dim(R_MDD_displacement_across_xyz)[1]*3)),
	pairGroup=rep(c(rep("MD1",dim(R_MDD_displacement_across_xyz)[1]),rep("MD2",dim(R_MDD_displacement_across_xyz)[1]),rep("D1D2",dim(R_MDD_displacement_across_xyz)[1])),3),
	pairName=c(paste(c(paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,2],sep="_"),paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,3],sep="_"),paste(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3],sep="_")),"x",sep="_"),
			   paste(c(paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,2],sep="_"),paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,3],sep="_"),paste(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3],sep="_")),"y",sep="_"),
			   paste(c(paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,2],sep="_"),paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,3],sep="_"),paste(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3],sep="_")),"z",sep="_")
			   )
)

dataTemp<-dataTemp[which(!duplicated(dataTemp$R)),]

dataTemp$pairGroup2<-ifelse(dataTemp$pairGroup=="D1D2","D1D2","MD")


ggplot(dataTemp,aes(x=pairGroup2,y=R,col=coordinate))+
geom_boxplot(width=0.6,size=0.2,outlier.size=0.1)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(x="Cell pair group",y="Pearson's R",title="Cross-embryo displacement")+
theme(legend.position=c(0.8,0.8))

ggsave("E:/canalizationENV/figs/mother_sister_paradox_cross_embryo_displacement_correlation.pdf",width=unit(1.6,"in"),height=unit(1.4,"in"))


## displacement: within-embryo
displacement_avg<-apply(displacement,2,mean,na.rm=T)
displacement_sd<-apply(displacement,2,sd,na.rm=T)
displacement_scaled<-t(apply(displacement,1,function(x){(x-displacement_avg)/displacement_sd}))

# plot 
R_MDD_displacement_within_xyz<-t(sapply(1:dim(displacement)[1],function(i){
	M_index<-match(paste0(mo_sis_pair_ccl[,1],"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D1_index<-match(paste0(mo_sis_pair_ccl[,2],"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D2_index<-match(paste0(mo_sis_pair_ccl[,3],"_ccl"),colnames(ccl_raw_nonNA_scaling))
	
	temp<-data.frame(
		M=rep(mo_sis_pair_ccl[,1],2),
		D=c(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3]),
		MD=paste(rep(mo_sis_pair_ccl[,1],2),c(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3]),sep="_")
	)
	temp<-temp[which(!duplicated(temp$MD)),]
	M_index_adj<-match(paste0(temp$M,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D_index_adj<-match(paste0(temp$D,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	
	return(c(cor(displacement_scaled[i,M_index_adj],displacement_scaled[i,D_index_adj],use="pairwise.complete.obs"),cor(displacement_scaled[i,D1_index],displacement_scaled[i,D2_index],use="pairwise.complete.obs"),
			 cor(displacement_scaled[i,M_index_adj+353],displacement_scaled[i,D_index_adj+353],use="pairwise.complete.obs"),cor(displacement_scaled[i,D1_index+353],displacement_scaled[i,D2_index+353],use="pairwise.complete.obs"),
			 cor(displacement_scaled[i,M_index_adj+353*2],displacement_scaled[i,D_index_adj+353*2],use="pairwise.complete.obs"),cor(displacement_scaled[i,D1_index+353*2],displacement_scaled[i,D2_index+353*2],use="pairwise.complete.obs")
			))
}))


dataTemp<-data.frame(
	R=c(R_MDD_displacement_within_xyz),
	coordinate=c(rep("x",dim(R_MDD_displacement_within_xyz)[1]*2),rep("y",dim(R_MDD_displacement_within_xyz)[1]*2),rep("z",dim(R_MDD_displacement_within_xyz)[1]*2)),
	pairName=rep(c(rep("MD",dim(R_MDD_displacement_within_xyz)[1]),rep("D1D2",dim(R_MDD_displacement_within_xyz)[1])),3)
)

ggplot(dataTemp,aes(x=pairName,y=R,col=coordinate))+
geom_boxplot(width=0.6,size=0.2,outlier.size=0.1)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(x="Cell pair group",y="Pearson's R",title="Within-embryo displacement")+
theme(legend.position=c(0.8,0.8))

ggsave("E:/canalizationENV/figs/mother_sister_paradox_within_embryo_displacement_correlation.pdf",width=unit(1.6,"in"),height=unit(1.4,"in"))



## position_noise  

## position_noise: cross-embryo
feaature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise),"_x"))

temp1<-unlist(lapply(feaature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feaature_names,function(x){which(cell_mo_da1_da2$da1==x)}))

mo_sis_pair_xyz<-cell_mo_da1_da2[intersect(temp1,temp2),]

# plot 
R_MDD_position_noise_across_xyz<-t(sapply(1:dim(mo_sis_pair_xyz)[1],function(i){
	M<-mo_sis_pair_xyz[i,1]
	D1<-mo_sis_pair_xyz[i,2]
	D2<-mo_sis_pair_xyz[i,3]
	
	M_index<-match(paste0(M,"_x"),colnames(q1_cell350_alignment_xyz_noise))
	D1_index<-match(paste0(D1,"_x"),colnames(q1_cell350_alignment_xyz_noise))
	D2_index<-match(paste0(D2,"_x"),colnames(q1_cell350_alignment_xyz_noise))
	
	return(c(cor(q1_cell350_alignment_xyz_noise[,M_index],q1_cell350_alignment_xyz_noise[,D1_index],use="pairwise.complete.obs"),cor(q1_cell350_alignment_xyz_noise[,M_index],q1_cell350_alignment_xyz_noise[,D2_index],use="pairwise.complete.obs"),cor(q1_cell350_alignment_xyz_noise[,D1_index],q1_cell350_alignment_xyz_noise[,D2_index],use="pairwise.complete.obs"),
			 cor(q1_cell350_alignment_xyz_noise[,M_index+714],q1_cell350_alignment_xyz_noise[,D1_index+714],use="pairwise.complete.obs"),cor(q1_cell350_alignment_xyz_noise[,M_index+714],q1_cell350_alignment_xyz_noise[,D2_index+714],use="pairwise.complete.obs"),cor(q1_cell350_alignment_xyz_noise[,D1_index+714],q1_cell350_alignment_xyz_noise[,D2_index+714],use="pairwise.complete.obs"),
			 cor(q1_cell350_alignment_xyz_noise[,M_index+714*2],q1_cell350_alignment_xyz_noise[,D1_index+714*2],use="pairwise.complete.obs"),cor(q1_cell350_alignment_xyz_noise[,M_index+714*2],q1_cell350_alignment_xyz_noise[,D2_index+714*2],use="pairwise.complete.obs"),cor(q1_cell350_alignment_xyz_noise[,D1_index+714*2],q1_cell350_alignment_xyz_noise[,D2_index+714*2],use="pairwise.complete.obs")
			))
}))


dataTemp<-data.frame(
	R=c(R_MDD_position_noise_across_xyz),
	coordinate=c(rep("x",dim(R_MDD_position_noise_across_xyz)[1]*3),rep("y",dim(R_MDD_position_noise_across_xyz)[1]*3),rep("z",dim(R_MDD_position_noise_across_xyz)[1]*3)),
	pairGroup=rep(c(rep("MD1",dim(R_MDD_position_noise_across_xyz)[1]),rep("MD2",dim(R_MDD_position_noise_across_xyz)[1]),rep("D1D2",dim(R_MDD_position_noise_across_xyz)[1])),3),
	pairName=c(paste(c(paste(mo_sis_pair_xyz[,1],mo_sis_pair_xyz[,2],sep="_"),paste(mo_sis_pair_xyz[,1],mo_sis_pair_xyz[,3],sep="_"),paste(mo_sis_pair_xyz[,2],mo_sis_pair_xyz[,3],sep="_")),"x",sep="_"),
			   paste(c(paste(mo_sis_pair_xyz[,1],mo_sis_pair_xyz[,2],sep="_"),paste(mo_sis_pair_xyz[,1],mo_sis_pair_xyz[,3],sep="_"),paste(mo_sis_pair_xyz[,2],mo_sis_pair_xyz[,3],sep="_")),"y",sep="_"),
			   paste(c(paste(mo_sis_pair_xyz[,1],mo_sis_pair_xyz[,2],sep="_"),paste(mo_sis_pair_xyz[,1],mo_sis_pair_xyz[,3],sep="_"),paste(mo_sis_pair_xyz[,2],mo_sis_pair_xyz[,3],sep="_")),"z",sep="_")
			   )
)

dataTemp<-dataTemp[which(!duplicated(dataTemp$R)),]

dataTemp$pairGroup2<-ifelse(dataTemp$pairGroup=="D1D2","D1D2","MD")


ggplot(dataTemp,aes(x=pairGroup2,y=R,col=coordinate))+
geom_boxplot(width=0.6,size=0.2,outlier.size=0.1)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(x="Cell pair group",y="Pearson's R",title="Cross-embryo position noise")+
theme(legend.position=c(0.8,0.8))

ggsave("E:/canalizationENV/figs/mother_sister_paradox_cross_embryo_position_noise_correlation.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


## displacement: within-embryo
position_noise_avg<-apply(q1_cell350_alignment_xyz_noise,2,mean,na.rm=T)
position_noise_sd<-apply(q1_cell350_alignment_xyz_noise,2,sd,na.rm=T)
position_noise_scaled<-t(apply(q1_cell350_alignment_xyz_noise,1,function(x){(x-displacement_avg)/displacement_sd}))

# plot 
R_MDD_position_noise_within_xyz<-t(sapply(1:dim(displacement)[1],function(i){
	M_index<-match(paste0(mo_sis_pair_ccl[,1],"_x"),colnames(q1_cell350_alignment_xyz_noise))
	D1_index<-match(paste0(mo_sis_pair_ccl[,2],"_x"),colnames(q1_cell350_alignment_xyz_noise))
	D2_index<-match(paste0(mo_sis_pair_ccl[,3],"_x"),colnames(q1_cell350_alignment_xyz_noise))
	
	temp<-data.frame(
		M=rep(mo_sis_pair_ccl[,1],2),
		D=c(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3]),
		MD=paste(rep(mo_sis_pair_ccl[,1],2),c(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3]),sep="_")
	)
	temp<-temp[which(!duplicated(temp$MD)),]
	M_index_adj<-match(paste0(temp$M,"_x"),colnames(q1_cell350_alignment_xyz_noise))
	D_index_adj<-match(paste0(temp$D,"_x"),colnames(q1_cell350_alignment_xyz_noise))
	
	return(c(cor(position_noise_scaled[i,M_index_adj],position_noise_scaled[i,D_index_adj],use="pairwise.complete.obs"),cor(position_noise_scaled[i,D1_index],position_noise_scaled[i,D2_index],use="pairwise.complete.obs"),
			 cor(position_noise_scaled[i,M_index_adj+714],position_noise_scaled[i,D_index_adj+714],use="pairwise.complete.obs"),cor(position_noise_scaled[i,D1_index+714],position_noise_scaled[i,D2_index+714],use="pairwise.complete.obs"),
			 cor(position_noise_scaled[i,M_index_adj+714*2],position_noise_scaled[i,D_index_adj+714*2],use="pairwise.complete.obs"),cor(position_noise_scaled[i,D1_index+714*2],position_noise_scaled[i,D2_index+714*2],use="pairwise.complete.obs")
			))
}))


dataTemp<-data.frame(
	R=c(R_MDD_position_noise_within_xyz),
	coordinate=c(rep("x",dim(R_MDD_position_noise_within_xyz)[1]*2),rep("y",dim(R_MDD_position_noise_within_xyz)[1]*2),rep("z",dim(R_MDD_position_noise_within_xyz)[1]*2)),
	pairName=rep(c(rep("MD",dim(R_MDD_position_noise_within_xyz)[1]),rep("D1D2",dim(R_MDD_position_noise_within_xyz)[1])),3)
)

ggplot(dataTemp,aes(x=pairName,y=R,col=coordinate))+
geom_boxplot(width=0.6,size=0.2,outlier.size=0.1)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(x="Cell pair group",y="Pearson's R",title="Within-embryo position noise")+
theme(legend.position=c(0.8,0.8))

ggsave("E:/canalizationENV/figs/mother_sister_paradox_within_embryo_position_noise_correlation.pdf",width=unit(1.6,"in"),height=unit(1.4,"in"))


## ccl: cross-embryo
ccl_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling),"_ccl"))

temp1<-unlist(lapply(ccl_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(ccl_names,function(x){which(cell_mo_da1_da2$da1==x)}))

mo_sis_pair_ccl<-cell_mo_da1_da2[intersect(temp1,temp2),]

# plot 
R_MDD_ccl_across<-t(sapply(1:dim(mo_sis_pair_ccl)[1],function(i){
	M<-mo_sis_pair_ccl[i,1]
	D1<-mo_sis_pair_ccl[i,2]
	D2<-mo_sis_pair_ccl[i,3]
	
	M_index<-match(paste0(M,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D1_index<-match(paste0(D1,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D2_index<-match(paste0(D2,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	
	return(c(cor(ccl_raw_nonNA_scaling[,M_index],ccl_raw_nonNA_scaling[,D1_index],use="pairwise.complete.obs"),cor(ccl_raw_nonNA_scaling[,M_index],ccl_raw_nonNA_scaling[,D2_index],use="pairwise.complete.obs"),cor(ccl_raw_nonNA_scaling[,D1_index],ccl_raw_nonNA_scaling[,D2_index],use="pairwise.complete.obs")
			))
}))


dataTemp<-data.frame(
	R=c(R_MDD_ccl_across),
	pairGroup=c(rep("MD1",dim(R_MDD_ccl_across)[1]),rep("MD2",dim(R_MDD_ccl_across)[1]),rep("D1D2",dim(R_MDD_ccl_across)[1])),
	pairName=c(paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,2],sep="_"),paste(mo_sis_pair_ccl[,1],mo_sis_pair_ccl[,3],sep="_"),paste(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3],sep="_"))
)

dataTemp<-dataTemp[which(!duplicated(dataTemp$R)),]

dataTemp$pairGroup2<-ifelse(dataTemp$pairGroup=="D1D2","D1D2","MD")


ggplot(dataTemp,aes(x=pairGroup2,y=R))+
geom_boxplot(width=0.6,size=0.2,outlier.size=0.1)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(x="Cell pair group",y="Pearson's R",title="Cross-embryo ccl")+
theme(legend.position=c(0.8,0.8))

ggsave("E:/canalizationENV/figs/mother_sister_paradox_cross_embryo_ccl_correlation.pdf",width=unit(1.6,"in"),height=unit(1.4,"in"))


## ccl: within-embryo
ccl_avg<-apply(ccl_raw_nonNA_scaling,2,mean,na.rm=T)
ccl_sd<-apply(ccl_raw_nonNA_scaling,2,sd,na.rm=T)
ccl_scaled<-t(apply(ccl_raw_nonNA_scaling,1,function(x){(x-ccl_avg)/ccl_sd}))

# plot 
R_MDD_ccl_within<-t(sapply(1:dim(ccl_raw_nonNA_scaling)[1],function(i){
	M_index<-match(paste0(mo_sis_pair_ccl[,1],"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D1_index<-match(paste0(mo_sis_pair_ccl[,2],"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D2_index<-match(paste0(mo_sis_pair_ccl[,3],"_ccl"),colnames(ccl_raw_nonNA_scaling))
	
	temp<-data.frame(
		M=rep(mo_sis_pair_ccl[,1],2),
		D=c(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3]),
		MD=paste(rep(mo_sis_pair_ccl[,1],2),c(mo_sis_pair_ccl[,2],mo_sis_pair_ccl[,3]),sep="_")
	)
	temp<-temp[which(!duplicated(temp$MD)),]
	M_index_adj<-match(paste0(temp$M,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D_index_adj<-match(paste0(temp$D,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	
	return(c(cor(ccl_scaled[i,M_index_adj],ccl_scaled[i,D_index_adj],use="pairwise.complete.obs"),cor(ccl_scaled[i,D1_index],ccl_scaled[i,D2_index],use="pairwise.complete.obs")			))
}))


dataTemp<-data.frame(
	R=c(R_MDD_ccl_within),
	pairName=c(rep("MD",dim(R_MDD_ccl_within)[1]),rep("D1D2",dim(R_MDD_ccl_within)[1]))
)

ggplot(dataTemp,aes(x=pairName,y=R))+
geom_boxplot(width=0.6,size=0.2,outlier.size=0.1)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(x="Cell pair group",y="Pearson's R",title="Within-embryo ccl")+
theme(legend.position=c(0.8,0.8))

ggsave("E:/canalizationENV/figs/mother_sister_paradox_within_embryo_ccl_correlation.pdf",width=unit(1.6,"in"),height=unit(1.4,"in"))


