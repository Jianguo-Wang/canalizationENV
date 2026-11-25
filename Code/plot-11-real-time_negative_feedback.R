

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
	"E:/papers/canalizationENV/data/realTime/xyz_x_model_mother_coef_M1D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_y_model_mother_coef_M1D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_z_model_mother_coef_M1D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_x_model_mother_coef_M1D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_y_model_mother_coef_M1D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_z_model_mother_coef_M1D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_x_model_mother_coef_M11D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_y_model_mother_coef_M11D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_z_model_mother_coef_M11D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_x_model_mother_coef_M11D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_y_model_mother_coef_M11D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_z_model_mother_coef_M11D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_x_model_mo_da1_R_M1D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_y_model_mo_da1_R_M1D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_z_model_mo_da1_R_M1D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_x_model_mo_da1_R_M1D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_y_model_mo_da1_R_M1D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_z_model_mo_da1_R_M1D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_x_model_mo_da1_R_M11D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_y_model_mo_da1_R_M11D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_z_model_mo_da1_R_M11D11.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_x_model_mo_da1_R_M11D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_y_model_mo_da1_R_M11D1.RData",
	"E:/papers/canalizationENV/data/realTime/xyz_z_model_mo_da1_R_M11D1.RData"
),projectFolderPath=focal_projectFolderPath)



############################################################################3333


dataTemp<-data.frame(
	kx=c(
	xyz_x_model_mother_coef_M11D1,
	xyz_x_model_mother_coef_M11D11,
	xyz_x_model_mother_coef_M1D1,
	xyz_x_model_mother_coef_M1D11),
	ky=c(
	xyz_y_model_mother_coef_M11D1,
	xyz_y_model_mother_coef_M11D11,
	xyz_y_model_mother_coef_M1D1,
	xyz_y_model_mother_coef_M1D11),
	kz=c(
	xyz_z_model_mother_coef_M11D1,
	xyz_z_model_mother_coef_M11D11,
	xyz_z_model_mother_coef_M1D1,
	xyz_z_model_mother_coef_M1D11),
	Rx=c(
	xyz_x_model_mo_da1_R_M11D1,
	xyz_x_model_mo_da1_R_M11D11,
	xyz_x_model_mo_da1_R_M1D1,
	xyz_x_model_mo_da1_R_M1D11),
	Ry=c(
	xyz_y_model_mo_da1_R_M11D1,
	xyz_y_model_mo_da1_R_M11D11,
	xyz_y_model_mo_da1_R_M1D1,
	xyz_y_model_mo_da1_R_M1D11),
	Rz=c(
	xyz_z_model_mo_da1_R_M11D1,
	xyz_z_model_mo_da1_R_M11D11,
	xyz_z_model_mo_da1_R_M1D1,
	xyz_z_model_mo_da1_R_M1D11),	
	type=c(rep("M11D1",706),rep("M11D11",706),rep("M1D1",706),rep("M1D11",706)))

	
### compare k 
# x 

ggplot(dataTemp,aes(x=type,y=kx))+
geom_boxplot(width=0.5,outlier.size=0.2)+
scale_y_continuous(breaks=c(-1,-0.5,0))+
labs(x="",y=expression(italic(k)[x]),title="")

ggsave("realTime_negative_feedback_kx.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


# y

ggplot(dataTemp,aes(x=type,y=ky))+
geom_boxplot(width=0.5,outlier.size=0.2)+
scale_y_continuous(breaks=c(-1,-0.5,0))+
labs(x="",y=expression(italic(k)[y]),title="")

ggsave("realTime_negative_feedback_ky.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


# z

ggplot(dataTemp,aes(x=type,y=kz))+
geom_boxplot(width=0.5,outlier.size=0.2)+
scale_y_continuous(breaks=c(-1,-0.5,0))+
labs(x="",y=expression(italic(k)[z]),title="")

ggsave("realTime_negative_feedback_kz.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

	
### compare R 
# x 

ggplot(dataTemp,aes(x=type,y=Rx))+
geom_boxplot(width=0.5,outlier.size=0.2)+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x="",y=expression(italic(R)[x]),title="")

ggsave("realTime_negative_feedback_Rx.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


# y

ggplot(dataTemp,aes(x=type,y=Ry))+
geom_boxplot(width=0.5,outlier.size=0.2)+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x="",y=expression(italic(R)[y]),title="")

ggsave("realTime_negative_feedback_Ry.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


# z

ggplot(dataTemp,aes(x=type,y=Rz))+
geom_boxplot(width=0.5,outlier.size=0.2)+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x="",y=expression(italic(R)[z]),title="")

ggsave("realTime_negative_feedback_Rz.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### compare k and R for x 

dataTemp<-data.frame(
	k_R=c(xyz_x_model_mother_coef_M11D1,
	xyz_x_model_mother_coef_M11D11,
	xyz_x_model_mother_coef_M1D1,
	xyz_x_model_mother_coef_M1D11,
	xyz_x_model_mo_da1_R_M11D1,
	xyz_x_model_mo_da1_R_M11D11,
	xyz_x_model_mo_da1_R_M1D1,
	xyz_x_model_mo_da1_R_M1D11),
	feature=factor(c(rep("k",706*4),rep("Pearson's R",706*4)),level=c("k","Pearson's R")),
	type=factor(rep(c(rep("M11D1",706),rep("M11D11",706),rep("M1D1",706),rep("M1D11",706)),2),level=c("M1D11","M1D1","M11D11","M11D1")))


ggplot(dataTemp,aes(y=type,x=k_R))+
geom_boxplot(width=0.5,outlier.size=0.2)+
scale_x_continuous(n.breaks=3)+
labs(y="",x="value",title="")+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))+
facet_wrap(vars(feature),nrow=1,scales="free_x")

ggsave("realTime_negative_feedback_kx_Rx.pdf",width=unit(3.2,"in"),height=unit(1.8,"in"))


### compare k and R for y

dataTemp<-data.frame(
	k_R=c(xyz_y_model_mother_coef_M11D1,
	xyz_y_model_mother_coef_M11D11,
	xyz_y_model_mother_coef_M1D1,
	xyz_y_model_mother_coef_M1D11,
	xyz_y_model_mo_da1_R_M11D1,
	xyz_y_model_mo_da1_R_M11D11,
	xyz_y_model_mo_da1_R_M1D1,
	xyz_y_model_mo_da1_R_M1D11),
	feature=factor(c(rep("k",706*4),rep("Pearson's R",706*4)),level=c("k","Pearson's R")),
	type=rep(c(rep("M11D1",706),rep("M11D11",706),rep("M1D1",706),rep("M1D11",706)),2),level=c("M1D11","M1D1","M11D11","M11D1"))


ggplot(dataTemp,aes(y=type,x=k_R))+
geom_boxplot(width=0.5,outlier.size=0.2)+
scale_x_continuous(n.breaks=3)+
labs(y="",x="value",title="")+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))+
facet_wrap(vars(feature),nrow=1,scales="free_x")

ggsave("realTime_negative_feedback_ky_Ry.pdf",width=unit(3.2,"in"),height=unit(1.8,"in"))


### compare k and R for z

dataTemp<-data.frame(
	k_R=c(xyz_z_model_mother_coef_M11D1,
	xyz_z_model_mother_coef_M11D11,
	xyz_z_model_mother_coef_M1D1,
	xyz_z_model_mother_coef_M1D11,
	xyz_z_model_mo_da1_R_M11D1,
	xyz_z_model_mo_da1_R_M11D11,
	xyz_z_model_mo_da1_R_M1D1,
	xyz_z_model_mo_da1_R_M1D11),
	feature=factor(c(rep("k",706*4),rep("Pearson's R",706*4)),level=c("k","Pearson's R")),
	type=rep(c(rep("M11D1",706),rep("M11D11",706),rep("M1D1",706),rep("M1D11",706)),2),level=c("M1D11","M1D1","M11D11","M11D1"))



ggplot(dataTemp,aes(y=type,x=k_R))+
geom_boxplot(width=0.5,outlier.size=0.2)+
scale_x_continuous(n.breaks=3)+
labs(y="",x="value",title="")+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))+
facet_wrap(vars(feature),nrow=1,scales="free_x")

ggsave("realTime_negative_feedback_kz_Rz.pdf",width=unit(3.2,"in"),height=unit(1.8,"in"))






