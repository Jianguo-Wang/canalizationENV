



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



focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_duration_prcomp.RData"
	),projectFolderPath=focal_projectFolderPath)


###

q1_cell350_alignment_xyz_noise_x<-q1_cell350_alignment_xyz_noise[,1:714]
q1_cell350_alignment_xyz_noise_y<-q1_cell350_alignment_xyz_noise[,(1+714):(714*2)]
q1_cell350_alignment_xyz_noise_z<-q1_cell350_alignment_xyz_noise[,(1+714*2):(714*3)]

#

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


####################################################################################

### reestimate k at different levels of scaling factors

###
	
###x

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da2,feature_names)]

#

k_quantile_x<-list()

for(i in 1:10){

	quantileInd<-which(embryo_size[,1]>quantile(embryo_size[,1],0.5-0.1*i/2)&embryo_size[,1]<quantile(embryo_size[,1],0.5+0.1*i/2))


	X_mo_scale_quantile<-scale(X_mo[quantileInd,],scale=FALSE)
	X_da1_scale_quantile<-scale(X_da1[quantileInd,],scale=FALSE)
	delta_X_scale_quantile<-X_da1_scale_quantile-X_mo_scale_quantile


	k_quantile_x[[i]]<-sapply(1:dim(X_mo_scale_quantile)[2],function(i){summary(lm(delta_X_scale_quantile[,i]~X_mo_scale_quantile[,i]))$coefficients[2,1]})

}

names(k_quantile_x)<-paste0("p",c(1:10)*10)


###y

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da2,feature_names)]

#

k_quantile_y<-list()

for(i in 1:10){

	quantileInd<-which(embryo_size[,2]>quantile(embryo_size[,2],0.5-0.1*i/2)&embryo_size[,2]<quantile(embryo_size[,2],0.5+0.1*i/2))


	X_mo_scale_quantile<-scale(X_mo[quantileInd,],scale=FALSE)
	X_da1_scale_quantile<-scale(X_da1[quantileInd,],scale=FALSE)
	delta_X_scale_quantile<-X_da1_scale_quantile-X_mo_scale_quantile


	k_quantile_y[[i]]<-sapply(1:dim(X_mo_scale_quantile)[2],function(i){summary(lm(delta_X_scale_quantile[,i]~X_mo_scale_quantile[,i]))$coefficients[2,1]})

}

names(k_quantile_y)<-paste0("p",c(1:10)*10)


###z

X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da2,feature_names)]

#

k_quantile_z<-list()

for(i in 1:10){

	quantileInd<-which(embryo_size[,3]>quantile(embryo_size[,3],0.5-0.1*i/2)&embryo_size[,3]<quantile(embryo_size[,3],0.5+0.1*i/2))


	X_mo_scale_quantile<-scale(X_mo[quantileInd,],scale=FALSE)
	X_da1_scale_quantile<-scale(X_da1[quantileInd,],scale=FALSE)
	delta_X_scale_quantile<-X_da1_scale_quantile-X_mo_scale_quantile


	k_quantile_z[[i]]<-sapply(1:dim(X_mo_scale_quantile)[2],function(i){summary(lm(delta_X_scale_quantile[,i]~X_mo_scale_quantile[,i]))$coefficients[2,1]})

}

names(k_quantile_z)<-paste0("p",c(1:10)*10)



save(k_quantile_x,file=paste0(focal_projectFolderPath,"/data/confounding/k_quantile_x.RData"))
save(k_quantile_y,file=paste0(focal_projectFolderPath,"/data/confounding/k_quantile_y.RData"))
save(k_quantile_z,file=paste0(focal_projectFolderPath,"/data/confounding/k_quantile_z.RData"))


### CCL 


feature_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


# reestimate kCCL at different levels of scaling factors


X_mo<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da1,feature_names)]
X_da2<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da2,feature_names)]

#

k_quantile_CCL<-list()

for(i in 1:10){

	quantileInd<-which(embryo_duration_prcomp>quantile(embryo_duration_prcomp,0.5-0.1*i/2)&embryo_duration_prcomp<quantile(embryo_duration_prcomp,0.5+0.1*i/2))


	X_mo_scale_quantile<-scale(X_mo[quantileInd,],scale=FALSE)
	X_da1_scale_quantile<-scale(X_da1[quantileInd,],scale=FALSE)
	delta_X_scale_quantile<-X_da1_scale_quantile-X_mo_scale_quantile


	k_quantile_CCL[[i]]<-sapply(1:dim(X_mo_scale_quantile)[2],function(i){summary(lm(delta_X_scale_quantile[,i]~X_mo_scale_quantile[,i]))$coefficients[2,1]})

}

names(k_quantile_CCL)<-paste0("p",c(1:10)*10)

save(k_quantile_CCL,file=paste0(focal_projectFolderPath,"/data/confounding/k_quantile_CCL.RData"))




####################################################################################

# distribution of scaling factor(x)

dataTemp<-data.frame(scaling=embryo_size[,1])
ggplot(dataTemp,aes(x=scaling))+
geom_density(lwd=0.3)+
geom_vline(xintercept=median(embryo_size[,1]),col="red",lty=2,lwd=0.2)+
geom_vline(xintercept=c(quantile(embryo_size[,1],0.5-0.1*3/2),quantile(embryo_size[,1],0.5+0.1*3/2)),col="blue",lty=2,lwd=0.2)+
geom_vline(xintercept=c(quantile(embryo_size[,1],0.5-0.1*7/2),quantile(embryo_size[,1],0.5+0.1*7/2)),col="blue",lty=2,lwd=0.2)+
scale_y_continuous(breaks=c(0,7.5))+
scale_x_continuous(limits=c(0.7,1),breaks=c(0.7,1))+
labs(x="Scaling factor (x)",y="Density",title="")

ggsave("Range_of_scaling_factor_x.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))


# distribution of scaling factor(y)

dataTemp<-data.frame(scaling=embryo_size[,2])
ggplot(dataTemp,aes(x=scaling))+
geom_density(lwd=0.3)+
geom_vline(xintercept=median(embryo_size[,2]),col="red",lty=2,lwd=0.2)+
geom_vline(xintercept=c(quantile(embryo_size[,2],0.5-0.1*3/2),quantile(embryo_size[,2],0.5+0.1*3/2)),col="blue",lty=2,lwd=0.2)+
geom_vline(xintercept=c(quantile(embryo_size[,2],0.5-0.1*7/2),quantile(embryo_size[,2],0.5+0.1*7/2)),col="blue",lty=2,lwd=0.2)+
scale_y_continuous(breaks=c(0,7.5))+
scale_x_continuous(limits=c(0.8,1.2),breaks=c(0.8,1.2))+
labs(x="Scaling factor (y)",y="Density",title="")

ggsave("Range_of_scaling_factor_y.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))


# distribution of scaling factor(z)

dataTemp<-data.frame(scaling=embryo_size[,3])
ggplot(dataTemp,aes(x=scaling))+
geom_density(lwd=0.3)+
geom_vline(xintercept=median(embryo_size[,3]),col="red",lty=2,lwd=0.2)+
geom_vline(xintercept=c(quantile(embryo_size[,3],0.5-0.1*3/2),quantile(embryo_size[,3],0.5+0.1*3/2)),col="blue",lty=2,lwd=0.2)+
geom_vline(xintercept=c(quantile(embryo_size[,3],0.5-0.1*7/2),quantile(embryo_size[,3],0.5+0.1*7/2)),col="blue",lty=2,lwd=0.2)+
scale_y_continuous(breaks=c(0,7.5))+
scale_x_continuous(limits=c(0.7,1.2),breaks=c(0.7,1.2))+
labs(x="Scaling factor (z)",y="Density",title="")

ggsave("Range_of_scaling_factor_z.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))


# kx reestimate

temp<-k_quantile_x
names(temp)<-c(1:10)*10
dataTemp<-stack(temp)

ggplot(dataTemp,aes(x=ind,y=values))+
geom_boxplot(outlier.size=0.2,outlier.stroke=0.2,width=0.5)+
#geom_hline(yintercept=c(-1,0),lty=2,col="red")+
labs(x="Range of scaling factor",y=expression(italic(k)[x]),title="")


ggsave("effect_scaling_factor_on_kx.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#ky reestimate

temp<-k_quantile_y
names(temp)<-c(1:10)*10
dataTemp<-stack(temp)

ggplot(dataTemp,aes(x=ind,y=values))+
geom_boxplot(outlier.size=0.2,outlier.stroke=0.2,width=0.5)+
#geom_hline(yintercept=c(-1,0),lty=2,col="red")+
labs(x="Range of scaling factor",y=expression(italic(k)[y]),title="")


ggsave("effect_scaling_factor_on_ky.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#kz reestimate

temp<-k_quantile_z
names(temp)<-c(1:10)*10
dataTemp<-stack(temp)

ggplot(dataTemp,aes(x=ind,y=values))+
geom_boxplot(outlier.size=0.2,outlier.stroke=0.2,width=0.5)+
#geom_hline(yintercept=c(-1,0),lty=2,col="red")+
labs(x="Range of scaling factor",y=expression(italic(k)[z]),title="")


ggsave("effect_scaleing_factor_on_kz.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



# distribution of scaling factor(CCL)

dataTemp<-data.frame(scaling=embryo_duration_prcomp)
ggplot(dataTemp,aes(x=scaling))+
geom_density(lwd=0.3)+
geom_vline(xintercept=median(embryo_duration_prcomp),col="red",lty=2,lwd=0.2)+
geom_vline(xintercept=c(quantile(embryo_duration_prcomp,0.5-0.1*3/2),quantile(embryo_duration_prcomp,0.5+0.1*3/2)),col="blue",lty=2,lwd=0.2)+
geom_vline(xintercept=c(quantile(embryo_duration_prcomp,0.5-0.1*7/2),quantile(embryo_duration_prcomp,0.5+0.1*7/2)),col="blue",lty=2,lwd=0.2)+
scale_y_continuous(breaks=c(0,7.5))+
scale_x_continuous(limits=c(0.9,1.5),breaks=c(0.9,1.5))+
labs(x="Scaling factor (CCL)",y="Density",title="")

ggsave("Range_of_scaling_factor_CCL.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))


# kx reestimate

temp<-k_quantile_CCL
names(temp)<-c(1:10)*10
dataTemp<-stack(temp)

ggplot(dataTemp,aes(x=ind,y=values))+
geom_boxplot(outlier.size=0.2,outlier.stroke=0.2,width=0.5)+
#geom_hline(yintercept=c(-1,0),lty=2,col="red")+
labs(x="Range of scaling factor",y=expression(italic(k)[CCL]),title="")


ggsave("effect_scaling_factor_on_kCCL.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))
















