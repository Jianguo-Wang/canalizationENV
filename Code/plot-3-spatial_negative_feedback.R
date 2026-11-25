



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


#######################################################################

#######################################################################

library(ggplot2)
library(reshape)
library(gplots)

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


feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


###########################################################################################


# delta ~ mo : case

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

xyz_x_model_mother_coef[1]

dataTemp<-data.frame(Mother=X_mo[,1],Delta=delta_X[,1],Daughter=X_da1[,1])
ggplot(dataTemp,aes(x=Mother,y=Delta))+geom_point(size=0.2,stroke=0.3,col="#F8766D")+
geom_vline(xintercept=0,lty=2,lwd=0.25,col="grey")+
geom_hline(yintercept=0,lty=2,lwd=0.25,col="grey")+
geom_abline(slope=-0.68,intercept=0,lwd=0.5,col="#F8766D")+
scale_x_continuous(breaks=c(-20,0,20))+
scale_y_continuous(breaks=c(-20,0,20))+
labs(x=expression(paste(italic(X)[M]," (ABal)")),y=expression(paste(Delta,italic(X)," (ABala-ABal)")),title=expression(paste("Slope ",paste(italic(k)[x]," = -",0.68))))


ggsave("spatial_xyz_x_delta_mother_case.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


par(mfrow=c(1,2))
hist(dataTemp$Mother)
hist(dataTemp$Daughter)

#139~2
dataTemp<-data.frame(Mother=X_mo[,1],Daughter=X_da1[,1])

ggplot(dataTemp,aes(x=Mother))+
geom_histogram(fill="grey",col="grey")+
geom_vline(xintercept=0,lty=2,lwd=0.25,col="#F8766D")+
scale_y_continuous(breaks=c(0,200,400))+
scale_x_continuous(limits=c(-30,30),breaks=c(-30,0,30))+
labs(y="Count",x="Position",title="Mother")

ggsave("spatial_xyz_x_case_mother_distribution.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


ggplot(dataTemp,aes(x=Daughter))+
geom_histogram(fill="grey",col="grey")+
geom_vline(xintercept=0,lty=2,lwd=0.25,col="#F8766D")+
scale_y_continuous(breaks=c(0,200,400))+
scale_x_continuous(limits=c(-30,30),breaks=c(-30,0,30))+
labs(y="Count",x="Position",title="Daughter")

ggsave("spatial_xyz_x_case_daughter_distribution.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


ggplot(dataTemp,aes(x=Mother,y=Daughter))+
geom_point(size=0.2,stroke=0.3,col="#F8766D",alpha=0.1)+
geom_vline(xintercept=0,lty=2,lwd=0.25,col="grey")+
geom_hline(yintercept=0,lty=2,lwd=0.25,col="grey")+
geom_abline(slope=1,intercept=0,lwd=0.5,col="black")+
scale_x_continuous(breaks=c(-20,0,20))+
scale_y_continuous(breaks=c(-20,0,20))+
labs(x="Mother noise",y="Daughter noise",title="")


ggsave("spatial_xyz_x_delta_mother_case.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### distribution of kx, ky and kz 


dataTemp<-data.frame(Coef=c(xyz_x_model_mother_coef,xyz_y_model_mother_coef,xyz_z_model_mother_coef),
					 Coo=factor(c(rep("x",length(xyz_x_model_mother_coef)),
								  rep("y",length(xyz_x_model_mother_coef)),
								  rep("z",length(xyz_x_model_mother_coef))),level=c("x","y","z")))

ggplot(dataTemp,aes(x=Coef,col=Coo))+geom_density(show.legend=TRUE)+
geom_vline(xintercept=c(-1,0),lty=2,col="grey",lwd=0.25)+
scale_y_continuous(breaks=c(0,2))+
scale_x_continuous(limits=c(-1.25,0),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)),y="Density",title="")+
scale_color_manual(values=c("#F8766D","#619CFF","#00BFC4"))+
theme(legend.position=c(0.1,0.8))

ggsave("spatial_xyz_k.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### shuffled cell pairs for x (main)


#R comparison

dataTemp<-data.frame(R=c(xyz_x_model_mo_da1_R,xyz_x_model_mo_da1_R_shuffle),type=c(rep("Observed",706),rep("Shuffled",706)))
ggplot(dataTemp,aes(x=R,col=type,fill=type))+
geom_density(alpha=0.5,lwd=0.1)+
geom_vline(xintercept=0,lty=2,lwd=0.25,col="grey")+
scale_y_continuous(limits=c(0,3.5),breaks=c(0,3))+
labs(x=expression(paste("Pearson's ",italic(R))),y="Density",title="Mother-daughter correlation (x)")+
theme(legend.position=c(0.35,0.9))

ggsave("spatial_xyz_x_shuffle_observe_R.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#k comparison 

dataTemp<-data.frame(kx=c(xyz_x_model_mother_coef,xyz_x_model_mother_coef_shuffle),
ky=c(xyz_y_model_mother_coef,xyz_y_model_mother_coef_shuffle),
kz=c(xyz_z_model_mother_coef,xyz_z_model_mother_coef_shuffle),
type=c(rep("Observed",706),rep("Shuffled",706)))

ggplot(dataTemp,aes(x=kx,col=type,fill=type))+
geom_density(alpha=0.5,lwd=0.1)+
geom_vline(xintercept=c(-1,0),lty=2,lwd=0.25,col=c("#619CFF","#619CFF"))+
scale_y_continuous()+
labs(x=expression(italic(k)[x]),y="Density", title="")+
theme(legend.position=c(0.15,0.9))

ggsave("spatial_xyz_x_shuffle_observe_k.pdf",width=unit(3.2,"in"),height=unit(1.4,"in"))

#> smaller

dataTemp<-data.frame(kx=c(xyz_x_model_mother_coef,xyz_x_model_mother_coef_shuffle),
ky=c(xyz_y_model_mother_coef,xyz_y_model_mother_coef_shuffle),
kz=c(xyz_z_model_mother_coef,xyz_z_model_mother_coef_shuffle),
type=c(rep("Observed",706),rep("Shuffled",706)))

ggplot(dataTemp,aes(x=kx,col=type,fill=type))+
geom_density(alpha=0.5,lwd=0.1)+
geom_vline(xintercept=c(-1,0),lty=2,lwd=0.25,col=c("#619CFF","#619CFF"))+
scale_y_continuous()+
labs(x=expression(italic(k)[x]),y="Density", title="")+
theme(legend.position=c(0.15,0.9))

ggsave("spatial_xyz_x_shuffle_observe_k_smaller.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### comparison between kx, ky and kz 

cor(xyz_x_model_mother_coef,xyz_y_model_mother_coef)
cor(xyz_x_model_mother_coef,xyz_z_model_mother_coef)
cor(xyz_y_model_mother_coef,xyz_z_model_mother_coef)

dataTemp<-data.frame(kx=rep(xyz_x_model_mother_coef,2),kyz=c(xyz_y_model_mother_coef,xyz_z_model_mother_coef),yz=c(rep("y",706),rep("z",706)))

ggplot(dataTemp,aes(x=kx,y=kyz,col=yz))+
geom_point(size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=1,lwd=0.3,lty=2)+
scale_x_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[x]),y=expression(paste(italic(k)[y]," or ",italic(k)[z])),title="")+
theme(legend.position=c(0.1,0.9))

ggsave("spatial_kx_ky_kz.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### comparison between k_Contr, k_RNAi, k_Envir for x (main)

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo


k_wt_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1:105,i]~X_mo[1:105,i])$coefficients[2]})
k_mt_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[106:1928,i]~X_mo[106:1928,i])$coefficients[2]})
k_env_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1929:2039,i]~X_mo[1929:2039,i])$coefficients[2]})


cor.test(k_wt_x,k_mt_x)
cor.test(k_wt_x,k_env_x)

dataTemp<-data.frame(k_wt_x,k_mt_env=c(k_mt_x,k_env_x),mt_env=c(rep("RNAi",706),rep("Env",706)))

ggplot(dataTemp,aes(x=k_wt_x,y=k_mt_env,col=mt_env))+
geom_point(size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=1,lwd=0.3,lty=2)+
scale_x_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[Control]),y=expression(paste(italic(k)[RNAi]," or ",italic(k)[Environment])),title=expression(italic(k)[x]))+
theme(legend.position=c(0.1,0.9))

ggsave("spatial_kContr_kRNAi_kENV_x.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### compare sister 

uniqueInd<-match(names(table(mo_sis_pair[,1]))[table(mo_sis_pair[,1])==2],mo_sis_pair[,1])

coef_mo_da1<-xyz_x_model_mother_coef[match(paste0(mo_sis_pair[,1],"_",mo_sis_pair[,2]),names(xyz_x_model_mother_coef))]
coef_mo_da2<-xyz_x_model_mother_coef[match(paste0(mo_sis_pair[,1],"_",mo_sis_pair[,3]),names(xyz_x_model_mother_coef))]
dataTemp<-data.frame(coef_mo_da1=coef_mo_da1[uniqueInd],coef_mo_da2=coef_mo_da2[uniqueInd],coo=rep("x",length(coef_mo_da2[uniqueInd])))

coef_mo_da1<-xyz_y_model_mother_coef[match(paste0(mo_sis_pair[,1],"_",mo_sis_pair[,2]),names(xyz_y_model_mother_coef))]
coef_mo_da2<-xyz_y_model_mother_coef[match(paste0(mo_sis_pair[,1],"_",mo_sis_pair[,3]),names(xyz_y_model_mother_coef))]
dataTemp<-rbind(dataTemp,data.frame(coef_mo_da1=coef_mo_da1[uniqueInd],coef_mo_da2=coef_mo_da2[uniqueInd],coo=rep("y",length(coef_mo_da2[uniqueInd]))))

coef_mo_da1<-xyz_z_model_mother_coef[match(paste0(mo_sis_pair[,1],"_",mo_sis_pair[,2]),names(xyz_z_model_mother_coef))]
coef_mo_da2<-xyz_z_model_mother_coef[match(paste0(mo_sis_pair[,1],"_",mo_sis_pair[,3]),names(xyz_z_model_mother_coef))]
dataTemp<-rbind(dataTemp,data.frame(coef_mo_da1=coef_mo_da1[uniqueInd],coef_mo_da2=coef_mo_da2[uniqueInd],coo=rep("z",length(coef_mo_da2[uniqueInd]))))

ggplot(dataTemp,aes(x=coef_mo_da1,y=coef_mo_da2,col=coo))+
geom_point(size=0.1,stroke=0.1)+
geom_abline(slope=1,intercept=0,lwd=0.3,lty=2)+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])),breaks=c(-1,0))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])),breaks=c(-1,0))+
scale_color_manual(values=c("#F8766D","#619CFF","#00BFC4"))+
labs(x=expression(italic(k)[Daughter[1]]),y=expression(italic(k)[Daughter[2]]),title="")+
theme(legend.position=c(0.2,0.8))

ggsave("spatial_sister_xyz.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




##########################################################################################################################
								supplementary 
##########################################################################################################################



### shuffled cell pairs for y


#R comparison

dataTemp<-data.frame(R=c(xyz_y_model_mo_da1_R,xyz_y_model_mo_da1_R_shuffle),type=c(rep("Observed",706),rep("Shuffled",706)))
ggplot(dataTemp,aes(x=R,col=type,fill=type))+
geom_density(alpha=0.5,lwd=0.1)+
geom_vline(xintercept=0,lty=2,lwd=0.25,col="grey")+
scale_y_continuous(limits=c(0,3.5),breaks=c(0,3))+
labs(x=expression(paste("Pearson's ",italic(R))),y="Density",title="Mother-daughter correlation (y)")+
theme(legend.position=c(0.35,0.9))

ggsave("spatial_xyz_y_shuffle_observe_R.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#k comparison 

dataTemp<-data.frame(kx=c(xyz_x_model_mother_coef,xyz_x_model_mother_coef_shuffle),
ky=c(xyz_y_model_mother_coef,xyz_y_model_mother_coef_shuffle),
kz=c(xyz_z_model_mother_coef,xyz_z_model_mother_coef_shuffle),
type=c(rep("Observed",706),rep("Shuffled",706)))

ggplot(dataTemp,aes(x=ky,col=type,fill=type))+
geom_density(alpha=0.5,lwd=0.1)+
geom_vline(xintercept=c(-1,0),lty=2,lwd=0.25,col=c("#619CFF","#619CFF"))+
scale_y_continuous()+
labs(x=expression(italic(k)[y]),y="Density", title="")+
theme(legend.position=c(0.15,0.9))

ggsave("spatial_xyz_y_shuffle_observe_k.pdf",width=unit(3.2,"in"),height=unit(1.4,"in"))




### shuffled cell pairs for z


#R comparison

dataTemp<-data.frame(R=c(xyz_z_model_mo_da1_R,xyz_z_model_mo_da1_R_shuffle),type=c(rep("Observed",706),rep("Shuffled",706)))
ggplot(dataTemp,aes(x=R,col=type,fill=type))+
geom_density(alpha=0.5,lwd=0.1)+
geom_vline(xintercept=0,lty=2,lwd=0.25,col="grey")+
scale_y_continuous(limits=c(0,4),breaks=c(0,3))+
labs(x=expression(paste("Pearson's ",italic(R))),y="Density",title="Mother-daughter correlation (z)")+
theme(legend.position=c(0.35,0.9))

ggsave("spatial_xyz_z_shuffle_observe_R.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#k comparison 

dataTemp<-data.frame(kx=c(xyz_x_model_mother_coef,xyz_x_model_mother_coef_shuffle),
ky=c(xyz_y_model_mother_coef,xyz_y_model_mother_coef_shuffle),
kz=c(xyz_z_model_mother_coef,xyz_z_model_mother_coef_shuffle),
type=c(rep("Observed",706),rep("Shuffled",706)))

ggplot(dataTemp,aes(x=kz,col=type,fill=type))+
geom_density(alpha=0.5,lwd=0.1)+
geom_vline(xintercept=c(-1,0),lty=2,lwd=0.25,col=c("#619CFF","#619CFF"))+
scale_y_continuous()+
labs(x=expression(italic(k)[z]),y="Density", title="")+
theme(legend.position=c(0.15,0.9))

ggsave("spatial_xyz_z_shuffle_observe_k.pdf",width=unit(3.2,"in"),height=unit(1.4,"in"))


######################

######################

### distribution of p-value for kx, ky and kz 

dataTemp<-data.frame(k=c(xyz_x_model_mother_coef,xyz_y_model_mother_coef,xyz_z_model_mother_coef),
	pvalue=c(xyz_x_model_mo_da1_pvalue,xyz_y_model_mo_da1_pvalue,xyz_z_model_mo_da1_pvalue),
	coo=c(rep("x",706),rep("y",706),rep("z",706)))

ggplot(dataTemp,aes(x=k,y=-log10(pvalue),col=coo))+geom_point(size=0.2,stroke=0.2)+
geom_hline(yintercept=-log10(0.05/706),lty=2,lwd=0.3)+
labs(x=expression(italic(k)),y=expression(paste("-",log[10],"(",italic(p),")")),title=expression(paste("96% of ",italic(k)," with ",italic(p),"<0.05/706")))+
theme(legend.position=c(0.2,0.8))

length(which(dataTemp$pvalue<0.05/706))/706/3

ggsave("spatial_k_pvalue.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### comparison between k_Contr, k_RNAi, k_Envir for x

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo


k_wt_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1:105,i]~X_mo[1:105,i])$coefficients[2]})
k_mt_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[106:1928,i]~X_mo[106:1928,i])$coefficients[2]})
k_env_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1929:2039,i]~X_mo[1929:2039,i])$coefficients[2]})


cor.test(k_wt_x,k_mt_x)
cor.test(k_wt_x,k_env_x)
cor.test(k_wt_x-k_env_x,k_wt_x-k_mt_x)



dataTemp<-data.frame(k_wt_x,k_mt_env=c(k_mt_x,k_env_x),mt_env=c(rep("RNAi",706),rep("Env",706)))

ggplot(dataTemp,aes(x=k_wt_x,y=k_mt_env,col=mt_env))+
geom_point(size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=1,lwd=0.3,lty=2)+
scale_x_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[Control]),y=expression(paste(italic(k)[RNAi]," or ",italic(k)[Environment])),title=expression(italic(k)[x]))+
theme(legend.position=c(0.1,0.9))

ggsave("spatial_kContr_kRNAi_kENV_x.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### comparison between k_Contr, k_RNAi, k_Envir for y

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)+714]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)+714]

delta_X<-X_da1-X_mo


k_wt_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1:105,i]~X_mo[1:105,i])$coefficients[2]})
k_mt_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[106:1928,i]~X_mo[106:1928,i])$coefficients[2]})
k_env_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1929:2039,i]~X_mo[1929:2039,i])$coefficients[2]})


cor.test(k_wt_x,k_mt_x)
cor.test(k_wt_x,k_env_x)
cor.test(k_wt_x-k_env_x,k_wt_x-k_mt_x)



dataTemp<-data.frame(k_wt_x,k_mt_env=c(k_mt_x,k_env_x),mt_env=c(rep("RNAi",706),rep("Env",706)))

ggplot(dataTemp,aes(x=k_wt_x,y=k_mt_env,col=mt_env))+
geom_point(size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=1,lwd=0.3,lty=2)+
scale_x_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[Control]),y=expression(paste(italic(k)[RNAi]," or ",italic(k)[Environment])),title=expression(italic(k)[y]))+
theme(legend.position=c(0.1,0.9))

ggsave("spatial_kContr_kRNAi_kENV_y.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### comparison between k_Contr, k_RNAi, k_Envir for z

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)+714*2]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)+714*2]

delta_X<-X_da1-X_mo


k_wt_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1:105,i]~X_mo[1:105,i])$coefficients[2]})
k_mt_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[106:1928,i]~X_mo[106:1928,i])$coefficients[2]})
k_env_x<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1929:2039,i]~X_mo[1929:2039,i])$coefficients[2]})


cor.test(k_wt_x,k_mt_x)
cor.test(k_wt_x,k_env_x)
cor.test(k_wt_x-k_env_x,k_wt_x-k_mt_x)


dataTemp<-data.frame(k_wt_x,k_mt_env=c(k_mt_x,k_env_x),mt_env=c(rep("RNAi",706),rep("Env",706)))

ggplot(dataTemp,aes(x=k_wt_x,y=k_mt_env,col=mt_env))+
geom_point(size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=1,lwd=0.3,lty=2)+
scale_x_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=c(-1.53,0.13),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[Control]),y=expression(paste(italic(k)[RNAi]," or ",italic(k)[Environment])),title=expression(italic(k)[z]))+
theme(legend.position=c(0.1,0.9))

ggsave("spatial_kContr_kRNAi_kENV_z.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


















