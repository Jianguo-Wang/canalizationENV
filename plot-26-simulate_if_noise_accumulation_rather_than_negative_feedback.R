


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


library(ggplot2)

focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R"),projectFolderPath=focal_projectFolderPath)



####################################################################################333

set.seed(10)
M<-rnorm(1000)
Var=c(0:100)*0.1

# positive feedback

k_initial<-0.5
set.seed(1)
D<-sapply(Var,function(t){M*(k_initial+1)+t*rnorm(1000)})

k_positive<-sapply(1:length(Var),function(i){coefficients(summary(lm(D[,i]~M)))[2,1]-1})
R_positive<-sapply(1:length(Var),function(i){cor(D[,i],M)})

par(mfrow=c(1,2))
plot(Var,k_positive)
plot(Var,R_positive)


# negative feedback

k_initial<--0.5
set.seed(1)
D<-sapply(Var,function(t){M*(k_initial+1)+t*rnorm(1000)})

k_negative<-sapply(1:length(Var),function(i){coefficients(summary(lm(D[,i]~M)))[2,1]-1})
R_negative<-sapply(1:length(Var),function(i){cor(D[,i],M)})

plot(Var,k_negative)
plot(Var,R_negative)


# no feedback

k_initial<-0
set.seed(1)
D<-sapply(Var,function(t){M*(k_initial+1)+t*rnorm(1000)})

k_no<-sapply(1:length(Var),function(i){coefficients(summary(lm(D[,i]~M)))[2,1]-1})
R_no<-sapply(1:length(Var),function(i){cor(D[,i],M)})

plot(Var,k_no)
plot(Var,R_no)


dataTemp<-data.frame(Variance=rep(Var,6),
	value=c(k_positive,k_negative,k_no,R_positive,R_negative,R_no),
	para=c(rep("k (observed)",length(k_positive)*3),rep("R (observed)",length(k_positive)*3)),
	type=c(rep("positive feedback (k=0.5)",length(k_positive)),rep("negative feedback (k=-0.5)",length(k_negative)),rep("no feedback (k=0)",length(k_no))))


ggplot(dataTemp,aes(x=Variance,y=value))+
geom_point(col="#F8766D",size=0.5,stroke=0.5)+
geom_hline(yintercept=c(-0.5,0,0.5,1),lty=2,lwd=0.3)+
facet_grid(vars(para),vars(type))+
labs(x="Daughter's accumulated noise variance",y="observed k or R",title="")


ggsave("simulate_if_accumulated_noise_rather_than_negative_feedback.pdf",width=unit(6.4,"in"),height=unit(4,"in"))
	
	
	
	
