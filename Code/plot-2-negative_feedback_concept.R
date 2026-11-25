

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
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"E:/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData"
	),projectFolderPath=focal_projectFolderPath)



#######################################################################################################

### noise & distribution in an example embryo 
dataTemp<-data.frame(
	noise=as.numeric(q1_cell350_alignment_xyz_noise[,1]),
	embryoType=factor(c(rep("WT",105),rep("RNAi",1823),rep("Env",111)),level=c("WT","RNAi","Env")),
	embryoNo=c(1:2039)
)

ggplot(dataTemp,aes(x=embryoNo,y=noise,col=embryoType))+
geom_point(stroke=0.3,size=0.35,alpha=0.9)+
theme(legend.position=c(0.8,0.9))+
labs(x="embryos (n=2039)",y="x-coordinate noise of cell ABal",title="")

ggsave("cell_noise_for_x_ABal.pdf",width=unit(3,"in"),height=unit(2.9,"in"))

ggplot(dataTemp,aes(y=noise,col=embryoType,fill=embryoType))+
geom_density(width=0.25,alpha=0.25,show.legend=F)+
scale_x_continuous(breaks=c(0,0.04,0.08))+
labs(x="Density",y="",title="")

ggsave("cell_noise_for_x_ABal_density.pdf",width=unit(1.2,"in"),height=unit(2.9,"in"))



###negative feedback 

xm<-rnorm(200)
xd<-0.3*xm+rnorm(200)*0.3

dataTemp<-data.frame(xm,xd)

#

ggplot(dataTemp,aes(x=xm,y=xd))+
geom_point(col="#F8766D",size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=1,col="black",lwd=0.3,lty=2)+
geom_abline(intercept=0,slope=0.3,col="red",lwd=0.2)+
scale_x_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
labs(x=expression(italic(X)[M]),y=expression(italic(X)[D]),title=expression(paste(italic(a),"=0.3")))

ggsave("ideograph_negative_feedback_XM_XD.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#

ggplot(dataTemp,aes(x=xm,y=xd-xm))+
geom_point(col="#F8766D",size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=-0.7,col="red",lwd=0.2)+
scale_x_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
labs(x=expression(italic(X)[M]),y=expression(paste(italic(X)[D],"-",italic(X)[M])),title=expression(paste(italic(k),"=-0.7")))

ggsave("ideograph_negative_feedback_XM_deltaX.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


###no feedback 

xm<-rnorm(200)
xd<-xm+rnorm(200)*0.3

dataTemp<-data.frame(xm,xd)


#

ggplot(dataTemp,aes(x=xm,y=xd))+
geom_point(col="#F8766D",size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=1,col="red",lwd=0.2)+
scale_x_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
labs(x=expression(italic(X)[M]),y=expression(italic(X)[D]),title=expression(paste(italic(a),"=1")))

ggsave("ideograph_no_feedback_XM_XD.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#

ggplot(dataTemp,aes(x=xm,y=xd-xm))+
geom_point(col="#F8766D",size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=0,col="red",lwd=0.2)+
scale_x_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
labs(x=expression(italic(X)[M]),y=expression(paste(italic(X)[D],"-",italic(X)[M])),title=expression(paste(italic(k),"=0")))

ggsave("ideograph_no_feedback_XM_deltaX.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))

#positive feedback 

xm<-rnorm(200)
xd<-1.3*xm+rnorm(200)*0.3

dataTemp<-data.frame(xm,xd)

#

ggplot(dataTemp,aes(x=xm,y=xd))+
geom_point(col="#F8766D",size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=1,col="black",lwd=0.3,lty=2)+
geom_abline(intercept=0,slope=1.3,col="red",lwd=0.2)+
scale_x_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
labs(x=expression(italic(X)[M]),y=expression(italic(X)[D]),title=expression(paste(italic(a),"=1.3")))

ggsave("ideograph_positive_feedback_XM_XD.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


#

ggplot(dataTemp,aes(x=xm,y=xd-xm))+
geom_point(col="#F8766D",size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=0.3,col="red",lwd=0.2)+
scale_x_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-2.5,2.5),breaks=c(-2,0,2))+
labs(x=expression(italic(X)[M]),y=expression(paste(italic(X)[D],"-",italic(X)[M])),title=expression(paste(italic(k),"=0.3")))

ggsave("ideograph_positive_feedback_XM_deltaX.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


####################################################################################################################

########## noise accumulation 

generation<-c(1:20)

#positive feedback
a_pos<-1.5
var_pos<-(1-a_pos^generation)/(1-a_pos)

#negative feedback
a_neg<-0.5
var_neg<-(1-a_neg^generation)/(1-a_neg)

### no feedback

a_no<-1
var_no<-generation

###
data<-data.frame(generation,var_no,var_pos,var_neg)

ggplot(data,aes(x=generation,y=var_no))+
geom_line(col="red")+
scale_x_continuous(breaks=c(0,10,20))+
scale_y_continuous(breaks=c(0,10,20))+
labs(x="Generation",y="Accumulated variance",title="a=1")

ggsave("ideograph_var_accumulation_no_feedback.pdf",width=unit(1.2,"in"),height=unit(1.6,"in"))


ggplot(data,aes(x=generation,y=var_pos))+
geom_line(col="red")+
scale_x_continuous(breaks=c(0,10,20))+
scale_y_continuous(breaks=c(0,6000))+
labs(x="Generation",y="Accumulated variance",title="a=1.22")

ggsave("ideograph_var_accumulation_positive_feedback.pdf",width=unit(1.2,"in"),height=unit(1.6,"in"))


ggplot(data,aes(x=generation,y=var_neg))+
geom_line(col="red")+
geom_hline(yintercept=1/(1-a_neg),col="grey",lty=2,lwd=0.25)+
scale_x_continuous(breaks=c(0,10,20))+
scale_y_continuous(limits=c(0,1/(1-a_neg)),breaks=c(0,1,2))+
labs(x="Generation",y="Accumulated variance",title="a=sqrt(0.5)")

ggsave("ideograph_var_accumulation_negative_feedback.pdf",width=unit(1.2,"in"),height=unit(1.6,"in"))


