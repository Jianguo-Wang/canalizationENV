
pdf("E:/canalizationENV/figs/Isodisplacement_curve.pdf")

curve(1/x, from = 0, to = 7, col = "blue", lwd = 2, 
      ylab = "Velocity (v)", xlab = "Time (t)", main = "Isodisplacement curve", 
      xlim = c(0, 7), ylim = c(0, 7))
	  
curve(10/x, from = 0, to = 7, col = "blue", lwd = 2, 
      xlim = c(0, 7), ylim = c(0, 7),add=T)

dev.off()


pdf("E:/canalizationENV/figs/Isovelocity curve.pdf")

curve(-1*x+5, from = 0, to = 15, col = "blue", lwd = 2, 
      ylab = "displacement (d)", xlab = "Time (t)", main = "Isovelocity curve", 
      xlim = c(0, 15), ylim = c(-12, 7))
	  
curve(-1*x, from = 0, to = 10, col = "blue", lwd = 2, 
      xlim = c(0, 15), ylim = c(-12, 7),add=T)
	  
curve(-1*x-5, from = 0, to = 5, col = "blue", lwd = 2, 
      xlim = c(0, 15), ylim = c(-12, 7),add=T)

abline(h=-10,lty=2)

dev.off()


pdf("E:/canalizationENV/figs/Isochronous curve.pdf")

curve(-0.5*x-5, from = 0, to = 10, col = "blue", lwd = 2, 
      ylab = "displacement (d)", xlab = "Time (t)", main = "Isochronous curve", 
      xlim = c(0, 12), ylim = c(-12, 7))
	  
curve(-1*x, from = 0, to = 10, col = "blue", lwd = 2, 
      xlim = c(0, 12), ylim = c(-12, 7),add=T)
	  
curve(-1.5*x+5, from = 0, to = 10, col = "blue", lwd = 2, 
      xlim = c(0, 12), ylim = c(-12, 7),add=T)

abline(h=-10,lty=2)
abline(v=10,lty=2)

dev.off()


##########################################################################################

#######################################################################


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

library(foreach)
library(glmnet)
library(parallel)
library(doParallel)
library(ggplot2)

library(dplyr)
library(tidyr)


focal_projectFolderPath<-r"(E:\canalizationENV)"

###

setwd(paste0(focal_projectFolderPath,"/data/velocity"))

###

loadResources(c(
	"E:/canalizationENV/functions/set_theme.R",
	"E:/canalizationENV/functions/rowColRecovery.R",
	"E:/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"E:/canalizationENV/data/basic/q11_cell350_alignment_xyz_noise.RData",
	"E:/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q2_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q3_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q4_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q5_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q6_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q7_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q8_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q9_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q10_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/q11_cell350_alignment_xyz.RData",
	"E:/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData", 
	# CCL is defined as the cell cycle length from the birth to the division of the cell / 
	# the division of its mother and the division of the cell 
	"E:/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"E:/canalizationENV/data/basic/q1_cell350_alignment_xyz_mean.RData",
	"E:/canalizationENV/data/basic/lineage_info_short.RData",
	"E:/canalizationENV/data/prediction/embryo_lethal.RData",
	"E:/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/canalizationENV/data/temporal/ccl_model_mother_coef.RData",
	"E:/canalizationENV/data/basic/cdt_raw_nonNA_scaling.RData"
),projectFolderPath=focal_projectFolderPath)



### cell714 & cell353 
cell714<-do.call(rbind,strsplit(colnames(q1_cell350_alignment_xyz)[1:714],"_x"))[,1]
cell353<-do.call(rbind,strsplit(colnames(ccl_raw_nonNA_scaling),"_ccl"))[,1]

cell714_mo<-cell_mo_da1_da2[match(cell714,cell_mo_da1_da2[,2]),1]
cell353_mo<-cell_mo_da1_da2[match(cell353,cell_mo_da1_da2[,2]),1]

### focal M-D pair 
cell714_ind<-match(paste0(cell714,"_x"),colnames(q1_cell350_alignment_xyz))
cell353_ind<-match(paste0(cell353,"_x"),colnames(q1_cell350_alignment_xyz))

cell714_mo_ind<-match(paste0(cell714_mo,"_x"),colnames(q1_cell350_alignment_xyz))
cell353_mo_ind<-match(paste0(cell353_mo,"_x"),colnames(q1_cell350_alignment_xyz))

### initial and end position (mother birth -> daughter birth), and tim 
start_pos<-q11_cell350_alignment_xyz[,c(cell353_mo_ind,cell353_mo_ind+714,cell353_mo_ind+714*2)]
end_pos<-q11_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]
tim<-cbind(ccl_raw_nonNA_scaling,ccl_raw_nonNA_scaling,ccl_raw_nonNA_scaling)

colnames(start_pos)<-c(paste(paste(cell353_mo,cell353,sep="_"),"x",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"y",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"z",sep="_"))
colnames(end_pos)<-c(paste(paste(cell353_mo,cell353,sep="_"),"x",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"y",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"z",sep="_"))
colnames(tim)<-c(paste(paste(cell353_mo,cell353,sep="_"),"x",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"y",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"z",sep="_"))

### displacement and velocity 
displacement<-end_pos-start_pos
colnames(displacement)<-c(paste(paste(cell353_mo,cell353,sep="_"),"x",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"y",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"z",sep="_"))
velocity<-displacement/tim

### distance and speed, and start_pos refinement 
distance<-abs(q11_cell350_alignment_xyz[,c(cell353_mo_ind,cell353_mo_ind+714,cell353_mo_ind+714*2)]-q1_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q1_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q2_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q2_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q3_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q3_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q4_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q4_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q5_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q5_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q6_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q6_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q7_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q7_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q8_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q8_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q9_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q9_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q10_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])+
		  abs(q10_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]-q11_cell350_alignment_xyz[,c(cell353_ind,cell353_ind+714,cell353_ind+714*2)])

colnames(distance)<-c(paste(paste(cell353_mo,cell353,sep="_"),"x",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"y",sep="_"),paste(paste(cell353_mo,cell353,sep="_"),"z",sep="_"))
speed<-distance/tim

start_pos_refine<-do.call(rbind,lapply(1:dim(q1_cell350_alignment_xyz)[1],function(i){
	ifelse(start_pos[i,]<end_pos[i,],start_pos[i,]-end_pos[i,],-start_pos[i,]+end_pos[i,])
}))

### M-D example: the relationship between displacement, tim and velocity 
exampleIndex<-which.min(sapply(1:353,function(i){
	cor(start_pos[,i],displacement[,i])
}))

#mo<-"ABarpa"
#da1<-"ABarpap"
#coordinate = "x" 

featureName<-colnames(displacement)[exampleIndex]

dataTemp<-data.frame(
	embryoID=c(1:dim(displacement)[1]),
	start_pos=start_pos[,featureName],
	displacement=displacement[,featureName],
	tim=tim[,featureName],
	velocity=velocity[,featureName]
)

cor(dataTemp)

# case: start_pos Vs displacement 
ggplot(dataTemp,aes(x=start_pos,y=displacement))+
geom_point(size=0.3,stroke=0.2,col=colors()[182])+
labs(x="Start position",y="Migration displacement",title=featureName)

ggsave("E:/canalizationENV/figs/physical_principle_case_displacement.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# case: start_pos Vs tim 
ggplot(dataTemp,aes(x=start_pos,y=tim))+
geom_point(size=0.3,stroke=0.2,col=colors()[182])+
labs(x="Start position",y="Migration time",title=featureName)

ggsave("E:/canalizationENV/figs/physical_principle_case_ccl.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# case: start_pos Vs velocity
ggplot(dataTemp,aes(x=start_pos,y=velocity))+
geom_point(size=0.3,stroke=0.2,col=colors()[182])+
labs(x="Start position",y="Migration velocity",title=featureName)

ggsave("E:/canalizationENV/figs/physical_principle_case_velocity.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### start_pos-velocity Vs satrt_pos-tim for all M-D pairs 
dataTemp<-data.frame(
	value=c(
		R_pos_tim=sapply(1:dim(displacement)[2],function(i){cor(start_pos[,i],tim[,i])}),
		R_pos_velocity=sapply(1:dim(displacement)[2],function(i){cor(start_pos[,i],velocity[,i])})
	),
	cause=c(rep("migration time",dim(displacement)[2]),rep("migration velocity",dim(displacement)[2])),	
	coordinate=rep(c(rep("x",353),rep("y",353),rep("z",353)),2)
)

# x coordinate
t.test(dataTemp[(dataTemp$coordinate=="x")&(dataTemp$cause=="migration time"),"value"],dataTemp[(dataTemp$coordinate=="x")&(dataTemp$cause=="migration velocity"),"value"])

ggplot(dataTemp[dataTemp$coordinate=="x",],aes(x=cause,y=value))+
geom_jitter(height=0,width=0.3,size=0.2,stroke=0.1,col="grey")+
geom_boxplot(outlier.size=0.5,col="#F8766D",fill=NA)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(y=expression(paste(italic(R)," with start position (x)")),title="")

ggsave("E:/canalizationENV/figs/physical_principle_velocity_ccl_correlation_x.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# y coordinate
t.test(dataTemp[(dataTemp$coordinate=="y")&(dataTemp$cause=="migration time"),"value"],dataTemp[(dataTemp$coordinate=="y")&(dataTemp$cause=="migration velocity"),"value"])

ggplot(dataTemp[dataTemp$coordinate=="y",],aes(x=cause,y=value))+
geom_jitter(height=0,width=0.3,size=0.2,stroke=0.1,col="grey")+
geom_boxplot(outlier.size=0.5,col="#F8766D",fill=NA)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(y=expression(paste(italic(R)," with start position (y)")),title="")

ggsave("E:/canalizationENV/figs/physical_principle_velocity_ccl_correlation_y.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# z coordinate
t.test(dataTemp[(dataTemp$coordinate=="z")&(dataTemp$cause=="migration time"),"value"],dataTemp[(dataTemp$coordinate=="z")&(dataTemp$cause=="migration velocity"),"value"])

ggplot(dataTemp[dataTemp$coordinate=="z",],aes(x=cause,y=value))+
geom_jitter(height=0,width=0.3,size=0.2,stroke=0.1,col="grey")+
geom_boxplot(outlier.size=0.5,col="#F8766D",fill=NA)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(y=expression(paste(italic(R)," with start position (z)")),title="")

ggsave("E:/canalizationENV/figs/physical_principle_velocity_ccl_correlation_z.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

### start_pos-speed Vs satrt_pos-tim for all M-D pairs 
dataTemp<-data.frame(
	value=c(
		R_pos_tim=sapply(1:dim(displacement)[2],function(i){cor(start_pos_refine[,i],tim[,i])}),
		R_pos_speed=sapply(1:dim(displacement)[2],function(i){cor(start_pos_refine[,i],speed[,i])})
	),
	cause=factor(c(rep("migration time",dim(displacement)[2]),rep("migration speed",dim(displacement)[2])),level=c("migration time","migration speed")),	
	coordinate=rep(c(rep("x",353),rep("y",353),rep("z",353)),2)
)

# x coordinate
t.test(dataTemp[(dataTemp$coordinate=="x")&(dataTemp$cause=="migration time"),"value"],dataTemp[(dataTemp$coordinate=="x")&(dataTemp$cause=="migration speed"),"value"])

ggplot(dataTemp[dataTemp$coordinate=="x",],aes(x=cause,y=value))+
geom_jitter(height=0,width=0.3,size=0.2,stroke=0.1,col="grey")+
geom_boxplot(outlier.size=0.5,col="#F8766D",fill=NA)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(y=expression(paste(italic(R)," with start position (x)")),title="")

ggsave("E:/canalizationENV/figs/physical_principle_speed_ccl_correlation_x.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# y coordinate
t.test(dataTemp[(dataTemp$coordinate=="y")&(dataTemp$cause=="migration time"),"value"],dataTemp[(dataTemp$coordinate=="y")&(dataTemp$cause=="migration speed"),"value"])

ggplot(dataTemp[dataTemp$coordinate=="y",],aes(x=cause,y=value))+
geom_jitter(height=0,width=0.3,size=0.2,stroke=0.1,col="grey")+
geom_boxplot(outlier.size=0.5,col="#F8766D",fill=NA)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(y=expression(paste(italic(R)," with start position (y)")),title="")

ggsave("E:/canalizationENV/figs/physical_principle_speed_ccl_correlation_y.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# z coordinate
t.test(dataTemp[(dataTemp$coordinate=="z")&(dataTemp$cause=="migration time"),"value"],dataTemp[(dataTemp$coordinate=="z")&(dataTemp$cause=="migration speed"),"value"])

ggplot(dataTemp[dataTemp$coordinate=="z",],aes(x=cause,y=value))+
geom_jitter(height=0,width=0.3,size=0.2,stroke=0.1,col="grey")+
geom_boxplot(outlier.size=0.5,col="#F8766D",fill=NA)+
geom_hline(yintercept=0,lty=2,lwd=0.5)+
labs(y=expression(paste(italic(R)," with start position (z)")),title="")

ggsave("E:/canalizationENV/figs/physical_principle_speed_ccl_correlation_z.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

### three example embryos for case M-D pair 
exampleIndex<-which.min(sapply(1:353,function(i){
	cor(start_pos[,i],displacement[,i])
}))

featureName<-colnames(displacement)[exampleIndex]
dataTemp<-data.frame(
	embryoID=rownames(q1_cell350_alignment_xyz_noise),
	start_pos=start_pos[,featureName],
	end_pos=end_pos[,featureName],
	displacement=displacement[,featureName],
	tim=tim[,featureName],
	velocity=velocity[,featureName],
	speed=speed[,featureName],
	start_pos_refine=start_pos_refine[,featureName]
)

dataTemp2<-dataTemp[which(abs(dataTemp$end_pos-mean(dataTemp$end_pos))<1&abs(dataTemp$tim-mean(dataTemp$tim))<1),]

dataTemp2[order(dataTemp2$start_pos),]

par(mfrow=c(1,3))
plot(dataTemp2$start_pos,dataTemp2$velocity)
plot(dataTemp2$start_pos,dataTemp2$tim)
plot(dataTemp2$tim,dataTemp2$velocity)

selectIndex<-sapply(quantile(dataTemp2$start_pos,c(0.1,0.86,0.98)),function(x){which.min(abs(x-dataTemp2$start_pos))})

dataTemp2[selectIndex,]

# theoretical pattern (β)
i=c(1:3)
ggplot() +
  geom_segment(aes(x = 0, y = dataTemp2[selectIndex[i],"start_pos"], xend = dataTemp2[selectIndex[i],"tim"], yend = dataTemp2[selectIndex[i],"end_pos"]), 
	color = "blue", size = 1.5,arrow = arrow(length = unit(0.15, "cm"), type = "closed"),lwd=0.5) +
  labs(title = "Three example embryos", x = "Migration time", y = "Position")
  #xlim(0,30)+ylim(25,55)

ggsave("E:/canalizationENV/figs/physical_principle_case_three_embryos.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# embryo 1
i=1
focalIndex<-match(setdiff(unlist(strsplit(featureName,"_")),c("x","y","z")),cell714)
col1<-rep("grey",714)
col1[focalIndex]<-c("#F8766D","#00BFC4")
shape1<-rep(1,714)
shape1[focalIndex]<-17
size1<-rep(0.2,714)
size1[focalIndex]<-1
stroke1<-rep(0.1,714)
stroke1[focalIndex]<-1

dataTemp3<-data.frame(
	x=q11_cell350_alignment_xyz[match(rownames(dataTemp2[selectIndex[i],]),rownames(q1_cell350_alignment_xyz)),1:714],
	y=q11_cell350_alignment_xyz[match(rownames(dataTemp2[selectIndex[i],]),rownames(q1_cell350_alignment_xyz)),c(1:714)+714]
)

ggplot(dataTemp3,aes(x=x,y=y))+
geom_point(size=size1,stroke=0.1,col=col1,shape=shape1)+
labs(x="Coordinate x",y="Coordinate y",title=rownames(dataTemp2[selectIndex[i],]))+
scale_y_continuous(breaks=c(-40,0,40))+
scale_x_continuous(breaks=c(-100,-50,0,50,100))+
geom_vline(xintercept=dataTemp3$x[focalIndex[2]],lty=2,lwd=0.2)


ggsave("E:/canalizationENV/figs/physical_principle_case_1_embryo.pdf",width=unit(1.7,"in"),height=unit(1.275,"in"))


# embryo 2
i=2
focalIndex<-match(setdiff(unlist(strsplit(featureName,"_")),c("x","y","z")),cell714)
col1<-rep("grey",714)
col1[focalIndex]<-c("#F8766D","#00BFC4")
shape1<-rep(1,714)
shape1[focalIndex]<-17
size1<-rep(0.2,714)
size1[focalIndex]<-1
stroke1<-rep(0.1,714)
stroke1[focalIndex]<-1

dataTemp3<-data.frame(x=q11_cell350_alignment_xyz[match(rownames(dataTemp2[selectIndex[i],]),rownames(q1_cell350_alignment_xyz)),1:714],
y=q11_cell350_alignment_xyz[match(rownames(dataTemp2[selectIndex[i],]),rownames(q1_cell350_alignment_xyz)),c(1:714)+714])

ggplot(dataTemp3,aes(x=x,y=y))+
geom_point(size=size1,stroke=0.1,col=col1,shape=shape1)+
labs(x="Coordinate x",y="Coordinate y",title=rownames(dataTemp2[selectIndex[i],]))+
scale_y_continuous(breaks=c(-40,0,40))+
scale_x_continuous(breaks=c(-100,-50,0,50,100))+
geom_vline(xintercept=dataTemp3$x[focalIndex[2]],lty=2,lwd=0.2)


ggsave("E:/canalizationENV/figs/physical_principle_case_2_embryo.pdf",width=unit(1.7,"in"),height=unit(1.275,"in"))


# embryo 3
i=3
focalIndex<-match(setdiff(unlist(strsplit(featureName,"_")),c("x","y","z")),cell714)
col1<-rep("grey",714)
col1[focalIndex]<-c("#F8766D","#00BFC4")
shape1<-rep(1,714)
shape1[focalIndex]<-17
size1<-rep(0.2,714)
size1[focalIndex]<-1
stroke1<-rep(0.1,714)
stroke1[focalIndex]<-1

dataTemp3<-data.frame(x=q11_cell350_alignment_xyz[match(rownames(dataTemp2[selectIndex[i],]),rownames(q1_cell350_alignment_xyz)),1:714],
y=q11_cell350_alignment_xyz[match(rownames(dataTemp2[selectIndex[i],]),rownames(q1_cell350_alignment_xyz)),c(1:714)+714])

ggplot(dataTemp3,aes(x=x,y=y))+
geom_point(size=size1,stroke=0.1,col=col1,shape=shape1)+
labs(x="Coordinate x",y="Coordinate y",title=rownames(dataTemp2[selectIndex[i],]))+
scale_y_continuous(breaks=c(-40,0,40))+
scale_x_continuous(breaks=c(-100,-50,0,50,100))+
geom_vline(xintercept=dataTemp3$x[focalIndex[2]],lty=2,lwd=0.2)


ggsave("E:/canalizationENV/figs/physical_principle_case_3_embryo.pdf",width=unit(1.7,"in"),height=unit(1.275,"in"))


### cross-embryo correlation analysis for the velocities of each M-D pair

ccl_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling),"_ccl"))

temp1<-unlist(lapply(ccl_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(ccl_names,function(x){which(cell_mo_da1_da2$da1==x)}))

mo_sis_pair_ccl<-cell_mo_da1_da2[intersect(temp1,temp2),]


R_MD_xyz<-t(sapply(1:dim(mo_sis_pair_ccl)[1],function(i){
	M<-mo_sis_pair_ccl[i,1]
	D1<-mo_sis_pair_ccl[i,2]
	M_index<-match(paste0(M,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	D1_index<-match(paste0(D1,"_ccl"),colnames(ccl_raw_nonNA_scaling))
	return(c(cor(velocity[,M_index],velocity[,D1_index]),cor(velocity[,M_index+353],velocity[,D1_index+353]),cor(velocity[,M_index+353*2],velocity[,D1_index+353*2])))
}))

dataTemp<-data.frame(
	R=c(R_MD_xyz),
	coordinate=c(rep("x",dim(R_MD_xyz)[1]),rep("y",dim(R_MD_xyz)[1]),rep("z",dim(R_MD_xyz)[1]))
)

ggplot(dataTemp,aes(x=R,col=coordinate))+
geom_density(alpha=0.25)+
labs(x="M-D velocity correlation",y="Number of M-D pairs",title="Cross embryos")+
scale_x_continuous(breaks=c(-0.8,-0.4,-0.0))+
geom_vline(xintercept=0,lty=2,lwd=0.5)+
theme(legend.position=c(0.14,0.8))

ggsave("E:/canalizationENV/figs/physical_principle_cross_embryo_R.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### within-embryo autocorrelation analysis on the velocities of all M-D pairs for each embryo

# for x 
lineage_velocity_x_avg<-matrix(apply(velocity[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_velocity_x_sd<-matrix(apply(velocity[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],2,sd,na.rm=T),nrow=dim(lineage_info)[1])

self_corr_x<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	return(sapply(1:4,function(i){cor(c(lineage_velocity_emb[1:(10-i),]),c(lineage_velocity_emb[(i+1):10,]),use="pairwise.complete.obs")}))
})

self_corr_x_p<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	return(sapply(1:4,function(i){cor.test(c(lineage_velocity_emb[1:(10-i),]),c(lineage_velocity_emb[(i+1):10,]),use="pairwise.complete.obs")$p.value}))
})

self_corr_x_random<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd

	AR<-c()
	for(i in 1:4){
		ancestral<-c(lineage_velocity_emb[1:(10-i),])
		progenitor<-c(t(apply(lineage_velocity_emb[(i+1):10,],1,sample,replace=FALSE)))
		AR[i]<-cor(ancestral,progenitor,use="pairwise.complete.obs")
	}

	return(AR)
})

self_corr_x_p_random<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd

	AR<-c()
	for(i in 1:4){
		ancestral<-c(lineage_velocity_emb[1:(10-i),])
		progenitor<-c(t(apply(lineage_velocity_emb[(i+1):10,],1,sample,replace=FALSE)))
		AR[i]<-cor.test(ancestral,progenitor,use="pairwise.complete.obs")$p.value
	}

	return(AR)
})

dataTemp<-data.frame(
	AR1=c(self_corr_x[1,]),
	type=c(rep("WT",105),rep("RNAi",1823),rep("Env",111))
)

ggplot(dataTemp,aes(x=AR1,col=type))+geom_density()


# for y 
lineage_velocity_y_avg<-matrix(apply(velocity[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_velocity_y_sd<-matrix(apply(velocity[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],2,sd,na.rm=T),nrow=dim(lineage_info)[1])

self_corr_y<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	return(sapply(1:4,function(i){cor(c(lineage_velocity_emb[1:(10-i),]),c(lineage_velocity_emb[(i+1):10,]),use="pairwise.complete.obs")}))
})

self_corr_y_p<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	return(sapply(1:4,function(i){cor.test(c(lineage_velocity_emb[1:(10-i),]),c(lineage_velocity_emb[(i+1):10,]),use="pairwise.complete.obs")$p.value}))
})

self_corr_y_random<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd

	AR<-c()
	for(i in 1:4){
		ancestral<-c(lineage_velocity_emb[1:(10-i),])
		progenitor<-c(t(apply(lineage_velocity_emb[(i+1):10,],1,sample,replace=FALSE)))
		AR[i]<-cor(ancestral,progenitor,use="pairwise.complete.obs")
	}

	return(AR)
})

self_corr_y_p_random<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd

	AR<-c()
	for(i in 1:4){
		ancestral<-c(lineage_velocity_emb[1:(10-i),])
		progenitor<-c(t(apply(lineage_velocity_emb[(i+1):10,],1,sample,replace=FALSE)))
		AR[i]<-cor.test(ancestral,progenitor,use="pairwise.complete.obs")$p.value
	}

	return(AR)
})


dataTemp<-data.frame(
	AR1=c(self_corr_y[1,]),
	type=c(rep("WT",105),rep("RNAi",1823),rep("Env",111))
)

ggplot(dataTemp,aes(x=AR1,col=type))+geom_density()



# for z 
lineage_velocity_z_avg<-matrix(apply(velocity[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],2,mean,na.rm=T),nrow=dim(lineage_info)[1])
lineage_velocity_z_sd<-matrix(apply(velocity[1:105,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],2,sd,na.rm=T),nrow=dim(lineage_info)[1])

self_corr_z<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	return(sapply(1:4,function(i){cor(c(lineage_velocity_emb[1:(10-i),]),c(lineage_velocity_emb[(i+1):10,]),use="pairwise.complete.obs")}))
})

self_corr_z_p<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	return(sapply(1:4,function(i){cor.test(c(lineage_velocity_emb[1:(10-i),]),c(lineage_velocity_emb[(i+1):10,]),use="pairwise.complete.obs")$p.value}))
})

self_corr_z_random<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd

	AR<-c()
	for(i in 1:4){
		ancestral<-c(lineage_velocity_emb[1:(10-i),])
		progenitor<-c(t(apply(lineage_velocity_emb[(i+1):10,],1,sample,replace=FALSE)))
		AR[i]<-cor(ancestral,progenitor,use="pairwise.complete.obs")
	}

	return(AR)
})

self_corr_z_p_random<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd

	AR<-c()
	for(i in 1:4){
		ancestral<-c(lineage_velocity_emb[1:(10-i),])
		progenitor<-c(t(apply(lineage_velocity_emb[(i+1):10,],1,sample,replace=FALSE)))
		AR[i]<-cor.test(ancestral,progenitor,use="pairwise.complete.obs")$p.value
	}

	return(AR)
})


dataTemp<-data.frame(
	AR1=c(self_corr_z[1,]),
	type=c(rep("WT",105),rep("RNAi",1823),rep("Env",111))
)

ggplot(dataTemp,aes(x=AR1,col=type))+geom_density()



# within-embryo velocity correlation 
dataTemp<-data.frame(
	R=c(self_corr_x[1,],self_corr_y[1,],self_corr_z[1,]),
	coordinate=c(rep("x",dim(self_corr_x)[2]),rep("y",dim(self_corr_x)[2]),rep("z",dim(self_corr_x)[2]))
)

ggplot(dataTemp,aes(x=R,col=coordinate))+
geom_density(alpha=0.25)+
geom_vline(xintercept=0,lty=2,lwd=0.5)+
labs(x="M-D velocity correlation",y="Number of embryos",title="Within embryos")+
scale_x_continuous(breaks=c(-0.6,-0.3,0.0))+
theme(legend.position=c(0.14,0.8))

ggsave("E:/canalizationENV/figs/physical_principle_within_embryo_R.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# net proportion of negative correlation 
dataTemp<-data.frame(
	tim=rep(1:4,6),
	net_neg=c(
		sapply(1:4,function(i){length(which(self_corr_x_p[i,]<0.05&self_corr_x[i,]<0))-length(which(self_corr_x_p[i,]<0.05&self_corr_x[i,]>0))}),
		sapply(1:4,function(i){length(which(self_corr_y_p[i,]<0.05&self_corr_y[i,]<0))-length(which(self_corr_y_p[i,]<0.05&self_corr_y[i,]>0))}),
		sapply(1:4,function(i){length(which(self_corr_z_p[i,]<0.05&self_corr_z[i,]<0))-length(which(self_corr_z_p[i,]<0.05&self_corr_z[i,]>0))}),
		sapply(1:4,function(i){length(which(self_corr_x_p_random[i,]<0.05&self_corr_x_random[i,]<0))-length(which(self_corr_x_p_random[i,]<0.05&self_corr_x_random[i,]>0))}),
		sapply(1:4,function(i){length(which(self_corr_y_p_random[i,]<0.05&self_corr_y_random[i,]<0))-length(which(self_corr_y_p_random[i,]<0.05&self_corr_y_random[i,]>0))}),
		sapply(1:4,function(i){length(which(self_corr_z_p_random[i,]<0.05&self_corr_z_random[i,]<0))-length(which(self_corr_z_p_random[i,]<0.05&self_corr_z_random[i,]>0))})
	),
	coordinate=rep(c(rep("x",4),rep("y",4),rep("z",4)),2),
	type=c(rep("Observed",12),rep("Shuffled",12))
)

total_embryos=dim(self_corr_x)[2]

ggplot(dataTemp, aes(x = tim, y = net_neg/total_embryos, color = coordinate, linetype = type)) +
	geom_line(size = 0.5) +
	geom_point(size = 0.75) +
	labs(x = "Generation interval (t)", y = "Probability of negative Rt", title = "Autocorrelation (Rt)")+
	scale_y_continuous(breaks=c(0,0.5,1.0))+
	theme(legend.position=c(0.85,0.8))  

ggsave("E:/canalizationENV/figs/physical_principle_AR_difference_with_lineage_depth.pdf",width=unit(2.2,"in"),height=unit(1.7,"in"))

### lineage patterns 

## summary of all embryos 

# summary of all embryos for x 

all_types_x<-unique(unlist(lapply(1:dim(self_corr_x)[2],function(emb){
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-
					       lineage_velocity_x_avg)/lineage_velocity_x_sd
	# sub-types 
	emb_lineage_type<-data.frame(t(lineage_velocity_emb[5:9,]))
	colnames(emb_lineage_type)<-paste("g",5:9,sep="")
	rownames(emb_lineage_type)<-colnames(lineage_info)
	emb_lineage_type$type<-apply(ifelse(emb_lineage_type[,1:4]>emb_lineage_type[,2:5],"N","P"),1,function(x){paste(x[!is.na(x)],collapse="")})
	focal_table<-table(emb_lineage_type$type)
	return(names(focal_table))
})))
all_types_x<-setdiff(all_types_x,"")

all_embryo_summary_x<-as.matrix(do.call(cbind,
	lapply(1:dim(self_corr_x)[2],function(emb){
		lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-
						       lineage_velocity_x_avg)/lineage_velocity_x_sd
		# sub-types 
		emb_lineage_type<-data.frame(t(lineage_velocity_emb[5:9,]))
		colnames(emb_lineage_type)<-paste("g",5:9,sep="")
		rownames(emb_lineage_type)<-colnames(lineage_info)
		emb_lineage_type$type<-apply(ifelse(emb_lineage_type[,1:4]>emb_lineage_type[,2:5],"N","P"),1,function(x){paste(x[!is.na(x)],collapse="")})
		focal_table<-table(emb_lineage_type$type)
		focal_count<-data.frame(type=names(focal_table),count=as.numeric(focal_table))
		return(focal_count[match(all_types_x,focal_count$type),2])
	})
))
colnames(all_embryo_summary_x)<-rownames(q1_cell350_alignment_xyz)
rownames(all_embryo_summary_x)<-all_types_x
all_embryo_summary_x[which(is.na(all_embryo_summary_x))]<-0

boxplot(t(all_embryo_summary_x))

dataTemp<-data.frame(
	count=c(all_embryo_summary_x),
	type=factor(rep(all_types_x,dim(all_embryo_summary_x)[2]),level=all_types_x[order(apply(all_embryo_summary_x,1,mean),decreasing=T)])
)

ggplot(dataTemp,aes(x=type,y=count))+
geom_boxplot(width=0.3,size=0.3,outlier.size=0.15)+
labs(x="Lineage patterns",y="Number of lineages",title="Summary of all embryos (x)")+
theme(axis.text.x=element_text(angle=90,hjust=1))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_all_embryos_summary_x.pdf",width=unit(3.2,"in"),height=unit(2.1,"in"))

# summary of all embryos for y 

all_types_y<-unique(unlist(lapply(1:dim(self_corr_y)[2],function(emb){
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-
					       lineage_velocity_y_avg)/lineage_velocity_y_sd
	# sub-types 
	emb_lineage_type<-data.frame(t(lineage_velocity_emb[5:9,]))
	colnames(emb_lineage_type)<-paste("g",5:9,sep="")
	rownames(emb_lineage_type)<-colnames(lineage_info)
	emb_lineage_type$type<-apply(ifelse(emb_lineage_type[,1:4]>emb_lineage_type[,2:5],"N","P"),1,function(x){paste(x[!is.na(x)],collapse="")})
	focal_table<-table(emb_lineage_type$type)
	return(names(focal_table))
})))
all_types_y<-setdiff(all_types_y,"")

all_embryo_summary_y<-as.matrix(do.call(cbind,
	lapply(1:dim(self_corr_y)[2],function(emb){
		lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-
						       lineage_velocity_y_avg)/lineage_velocity_y_sd
		# sub-types 
		emb_lineage_type<-data.frame(t(lineage_velocity_emb[5:9,]))
		colnames(emb_lineage_type)<-paste("g",5:9,sep="")
		rownames(emb_lineage_type)<-colnames(lineage_info)
		emb_lineage_type$type<-apply(ifelse(emb_lineage_type[,1:4]>emb_lineage_type[,2:5],"N","P"),1,function(x){paste(x[!is.na(x)],collapse="")})
		focal_table<-table(emb_lineage_type$type)
		focal_count<-data.frame(type=names(focal_table),count=as.numeric(focal_table))
		return(focal_count[match(all_types_y,focal_count$type),2])
	})
))
colnames(all_embryo_summary_y)<-rownames(q1_cell350_alignment_xyz)
rownames(all_embryo_summary_y)<-all_types_y
all_embryo_summary_y[which(is.na(all_embryo_summary_y))]<-0

boxplot(t(all_embryo_summary_y))

dataTemp<-data.frame(
	count=c(all_embryo_summary_y),
	type=factor(rep(all_types_y,dim(all_embryo_summary_y)[2]),level=all_types_y[order(apply(all_embryo_summary_y,1,mean),decreasing=T)])
)

ggplot(dataTemp,aes(x=type,y=count))+
geom_boxplot(width=0.3,size=0.3,outlier.size=0.15)+
labs(x="Lineage patterns",y="Number of lineages",title="Summary of all embryos (y)")+
theme(axis.text.x=element_text(angle=90,hjust=1))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_all_embryos_summary_y.pdf",width=unit(3.2,"in"),height=unit(2.1,"in"))

# summary of all embryos for z 

all_types_z<-unique(unlist(lapply(1:dim(self_corr_z)[2],function(emb){
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-
					       lineage_velocity_z_avg)/lineage_velocity_z_sd
	# sub-types 
	emb_lineage_type<-data.frame(t(lineage_velocity_emb[5:9,]))
	colnames(emb_lineage_type)<-paste("g",5:9,sep="")
	rownames(emb_lineage_type)<-colnames(lineage_info)
	emb_lineage_type$type<-apply(ifelse(emb_lineage_type[,1:4]>emb_lineage_type[,2:5],"N","P"),1,function(x){paste(x[!is.na(x)],collapse="")})
	focal_table<-table(emb_lineage_type$type)
	return(names(focal_table))
})))
all_types_z<-setdiff(all_types_z,"")

all_embryo_summary_z<-as.matrix(do.call(cbind,
	lapply(1:dim(self_corr_z)[2],function(emb){
		lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-
						       lineage_velocity_z_avg)/lineage_velocity_z_sd
		# sub-types 
		emb_lineage_type<-data.frame(t(lineage_velocity_emb[5:9,]))
		colnames(emb_lineage_type)<-paste("g",5:9,sep="")
		rownames(emb_lineage_type)<-colnames(lineage_info)
		emb_lineage_type$type<-apply(ifelse(emb_lineage_type[,1:4]>emb_lineage_type[,2:5],"N","P"),1,function(x){paste(x[!is.na(x)],collapse="")})
		focal_table<-table(emb_lineage_type$type)
		focal_count<-data.frame(type=names(focal_table),count=as.numeric(focal_table))
		return(focal_count[match(all_types_z,focal_count$type),2])
	})
))
colnames(all_embryo_summary_z)<-rownames(q1_cell350_alignment_xyz)
rownames(all_embryo_summary_z)<-all_types_z
all_embryo_summary_z[which(is.na(all_embryo_summary_z))]<-0

boxplot(t(all_embryo_summary_z))

dataTemp<-data.frame(
	count=c(all_embryo_summary_z),
	type=factor(rep(all_types_z,dim(all_embryo_summary_z)[2]),level=all_types_z[order(apply(all_embryo_summary_z,1,mean),decreasing=T)])
)

ggplot(dataTemp,aes(x=type,y=count))+
geom_boxplot(width=0.3,size=0.3,outlier.size=0.15)+
labs(x="Lineage patterns",y="Number of lineages",title="Summary of all embryos (z)")+
theme(axis.text.x=element_text(angle=90,hjust=1))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_all_embryos_summary_z.pdf",width=unit(3.2,"in"),height=unit(2.1,"in"))

# summary of all embryos (combined)

dataTemp<-data.frame(
	count=c(c(all_embryo_summary_x),c(all_embryo_summary_y),c(all_embryo_summary_z)),
	coordinate=c(rep("x",length(c(all_embryo_summary_x))),rep("y",length(c(all_embryo_summary_y))),rep("z",length(c(all_embryo_summary_x)))),
	type=c(rep(all_types_x,dim(all_embryo_summary_x)[2]),rep(all_types_y,dim(all_embryo_summary_y)[2]),rep(all_types_z,dim(all_embryo_summary_z)[2]))
)
type_table<-data.frame(type=unique(dataTemp$type),avg=sapply(unique(dataTemp$type),function(x){sum(dataTemp[dataTemp[,3]==x,1])}))
type_table<-type_table[order(type_table$avg,decreasing=T),]

dataTemp$type=factor(dataTemp$type,level=type_table$type)

ggplot(dataTemp,aes(x=type,y=count,col=coordinate))+
geom_boxplot(width=0.7,size=0.2,outlier.size=0.15)+
labs(x="Lineage patterns",y="Number of lineages",title="Summary of all embryos")+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_all_embryos_summary_xyz.pdf",width=unit(6.4,"in"),height=unit(2.1,"in"))


## all lineages in an example embryo 

emb=which.min(self_corr_x[1,]) # 1649; 
lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
sapply(1:4,function(i){cor(c(lineage_velocity_emb[1:(10-i),]),c(lineage_velocity_emb[(i+1):10,]),use="pairwise.complete.obs")})

temp_fun<-function(i){
plot(lineage_velocity_emb[5:9,i],type="b")
abline(h=0,lty=2)
}

pdf("E:/canalizationENV/figs/example_embryo_lineage_velocity.pdf")
for(i in 1:dim(lineage_velocity_emb)[2]){
	print(i)
	temp_fun(i)
} 
dev.off()


# sub-types 
emb_lineage_type<-data.frame(
t(lineage_velocity_emb[5:9,])
)
colnames(emb_lineage_type)<-paste("g",5:9,sep="")
rownames(emb_lineage_type)<-colnames(lineage_info)
emb_lineage_type$type<-apply(ifelse(emb_lineage_type[,1:4]>emb_lineage_type[,2:5],"N","P"),1,function(x){
	paste(x[!is.na(x)],collapse="")
})
sort(table(emb_lineage_type$type),decreasing=T)

emb_lineage_type_MD<-data.frame(
	M=c(as.matrix(emb_lineage_type[,1:4])),
	D=c(as.matrix(emb_lineage_type[,2:5])),
	type=rep(emb_lineage_type$type,4),
	Group=rep(c(rep("ABal",64),rep("ABar",64),rep("ABpl",64),rep("ABpr",64),rep("EMS",64),rep("P2",41)),4)
)

ggplot(emb_lineage_type_MD,aes(x=M,y=D,col=type))+
geom_point()

ggplot(emb_lineage_type_MD,aes(x=M,y=D,col=Group))+
geom_point()

cor(emb_lineage_type_MD[emb_lineage_type_MD$Group=="ABal",1:2],use="pairwise.complete.obs")

# total 
sampleIndex<-c(1:dim(lineage_velocity_emb)[2])
dataTemp<-data.frame(
	velocity=c(lineage_velocity_emb[5:9,sampleIndex]),
	generation=rep(5:9,length(sampleIndex)),
	lineage=rep(colnames(lineage_info)[sampleIndex],each=5)
)

ggplot(dataTemp,aes(x=generation,y=velocity,group=lineage,col=lineage))+
geom_line(show.legend=F,lwd=0.1)+
labs(x="Generation",y="Velocity",title="Velocity dynamics in an embryo")

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_example_embryo_all.pdf",width=unit(2.2,"in"),height=unit(2.1,"in"))


# sub-types - 1
typeIndex<-1
typeName=names(sort(table(emb_lineage_type$type),decreasing=T))[typeIndex]
sampleIndex<-which(emb_lineage_type$type==typeName)
typeRatio<-length(sampleIndex)/dim(lineage_velocity_emb)[2]
dataTemp<-data.frame(
	velocity=c(lineage_velocity_emb[5:9,sampleIndex]),
	generation=rep(5:9,length(sampleIndex)),
	lineage=rep(colnames(lineage_info)[sampleIndex],each=5)
)

ggplot(dataTemp,aes(x=generation,y=velocity,group=lineage,col=lineage))+
geom_line(show.legend=F,lwd=0.1)+
labs(x="Generation",y="Standardized velocity",
	title=paste0(typeName," (",round(typeRatio,4)*100," %)"))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_1.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))


# sub-types - 2
typeIndex<-2
typeName=names(sort(table(emb_lineage_type$type),decreasing=T))[typeIndex]
sampleIndex<-which(emb_lineage_type$type==typeName)
typeRatio<-length(sampleIndex)/dim(lineage_velocity_emb)[2]
dataTemp<-data.frame(
	velocity=c(lineage_velocity_emb[5:9,sampleIndex]),
	generation=rep(5:9,length(sampleIndex)),
	lineage=rep(colnames(lineage_info)[sampleIndex],each=5)
)

ggplot(dataTemp,aes(x=generation,y=velocity,group=lineage,col=lineage))+
geom_line(show.legend=F,lwd=0.1)+
labs(x="Generation",y="Standardized velocity",
	title=paste0(typeName," (",round(typeRatio,4)*100," %)"))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_2.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))


# sub-types - 3
typeIndex<-3
typeName=names(sort(table(emb_lineage_type$type),decreasing=T))[typeIndex]
sampleIndex<-which(emb_lineage_type$type==typeName)
typeRatio<-length(sampleIndex)/dim(lineage_velocity_emb)[2]
dataTemp<-data.frame(
	velocity=c(lineage_velocity_emb[5:9,sampleIndex]),
	generation=rep(5:9,length(sampleIndex)),
	lineage=rep(colnames(lineage_info)[sampleIndex],each=5)
)

ggplot(dataTemp,aes(x=generation,y=velocity,group=lineage,col=lineage))+
geom_line(show.legend=F,lwd=0.1)+
labs(x="Generation",y="Standardized velocity",
	title=paste0(typeName," (",round(typeRatio,4)*100," %)"))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_3.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))


# sub-types - 4
typeIndex<-4
typeName=names(sort(table(emb_lineage_type$type),decreasing=T))[typeIndex]
sampleIndex<-which(emb_lineage_type$type==typeName)
typeRatio<-length(sampleIndex)/dim(lineage_velocity_emb)[2]
dataTemp<-data.frame(
	velocity=c(lineage_velocity_emb[5:9,sampleIndex]),
	generation=rep(5:9,length(sampleIndex)),
	lineage=rep(colnames(lineage_info)[sampleIndex],each=5)
)

ggplot(dataTemp,aes(x=generation,y=velocity,group=lineage,col=lineage))+
geom_line(show.legend=F,lwd=0.1)+
labs(x="Generation",y="Standardized velocity",
	title=paste0(typeName," (",round(typeRatio,4)*100," %)"))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_4.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))


# sub-types - 5
typeIndex<-5
typeName=names(sort(table(emb_lineage_type$type),decreasing=T))[typeIndex]
sampleIndex<-which(emb_lineage_type$type==typeName)
typeRatio<-length(sampleIndex)/dim(lineage_velocity_emb)[2]
dataTemp<-data.frame(
	velocity=c(lineage_velocity_emb[5:9,sampleIndex]),
	generation=rep(5:9,length(sampleIndex)),
	lineage=rep(colnames(lineage_info)[sampleIndex],each=5)
)

ggplot(dataTemp,aes(x=generation,y=velocity,group=lineage,col=lineage))+
geom_line(show.legend=F,lwd=0.1)+
labs(x="Generation",y="Standardized velocity",
	title=paste0(typeName," (",round(typeRatio,4)*100," %)"))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_5.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))


# sub-types - 6
typeIndex<-6
typeName=names(sort(table(emb_lineage_type$type),decreasing=T))[typeIndex]
sampleIndex<-which(emb_lineage_type$type==typeName)
typeRatio<-length(sampleIndex)/dim(lineage_velocity_emb)[2]
dataTemp<-data.frame(
	velocity=c(lineage_velocity_emb[5:9,sampleIndex]),
	generation=rep(5:9,length(sampleIndex)),
	lineage=rep(colnames(lineage_info)[sampleIndex],each=5)
)

ggplot(dataTemp,aes(x=generation,y=velocity,group=lineage,col=lineage))+
geom_line(show.legend=F,lwd=0.1)+
labs(x="Generation",y="Standardized velocity",
	title=paste0(typeName," (",round(typeRatio,4)*100," %)"))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_6.pdf",width=unit(1.6,"in"),height=unit(1.2,"in"))

### multivariate regression

# for x 
multivariate_x<-t(sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_collapse<-c(lineage_velocity_emb[1:9,])
	D_collapse<-c(lineage_velocity_emb[2:10,])
	
	lineage_start_pos_emb<-matrix(start_pos[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])
	MD_pos<-c(lineage_start_pos_emb[2:10,])
	
	focal_lm<-lm(D_collapse~M_collapse+MD_pos)
	return(c(coefficients(summary(focal_lm))[,c(1,4)]))
}))
colnames(multivariate_x)<-c("beta_0","beta_D","beta_pos","p_0","p_D","p_pos")
rownames(multivariate_x)<-rownames(ccl_raw_nonNA_scaling)
multivariate_x<-data.frame(multivariate_x)
multivariate_x$AR1<-c(self_corr_x[1,])

# for y
multivariate_y<-t(sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_collapse<-c(lineage_velocity_emb[1:9,])
	D_collapse<-c(lineage_velocity_emb[2:10,])
	
	lineage_start_pos_emb<-matrix(start_pos[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])
	MD_pos<-c(lineage_start_pos_emb[2:10,])
	
	focal_lm<-lm(D_collapse~M_collapse+MD_pos)
	return(c(coefficients(summary(focal_lm))[,c(1,4)]))
}))
colnames(multivariate_y)<-c("beta_0","beta_D","beta_pos","p_0","p_D","p_pos")
rownames(multivariate_y)<-rownames(ccl_raw_nonNA_scaling)
multivariate_y<-data.frame(multivariate_y)
multivariate_y$AR1<-c(self_corr_y[1,])

# for z
multivariate_z<-t(sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_collapse<-c(lineage_velocity_emb[1:9,])
	D_collapse<-c(lineage_velocity_emb[2:10,])
	
	lineage_start_pos_emb<-matrix(start_pos[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])
	MD_pos<-c(lineage_start_pos_emb[2:10,])
	
	focal_lm<-lm(D_collapse~M_collapse+MD_pos)
	return(c(coefficients(summary(focal_lm))[,c(1,4)]))
}))
colnames(multivariate_z)<-c("beta_0","beta_D","beta_pos","p_0","p_D","p_pos")
rownames(multivariate_z)<-rownames(ccl_raw_nonNA_scaling)
multivariate_z<-data.frame(multivariate_z)
multivariate_z$AR1<-c(self_corr_z[1,])



ggplot(multivariate_x,aes(x=beta_D,y=-log10(p_D)))+
geom_point()

ggplot(multivariate_x,aes(x=AR1,y=beta_D))+
geom_point()+
geom_hline(yintercept=0,lty=2)+
geom_vline(xintercept=0,lty=2)


### partial correlation analysis

# for x 
partial_x<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])-lineage_velocity_x_avg)/lineage_velocity_x_sd
	M_collapse<-c(lineage_velocity_emb[1:9,])
	D_collapse<-c(lineage_velocity_emb[2:10,])
	
	lineage_start_pos_emb<-matrix(start_pos[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))],nrow=dim(lineage_info)[1])
	MD_pos<-c(lineage_start_pos_emb[2:10,])
	
	dataTemp<-na.omit(cbind(D_collapse,M_collapse,MD_pos))
	res<-(cor(dataTemp[,1],dataTemp[,2])-cor(dataTemp[,1],dataTemp[,3])*cor(dataTemp[,2],dataTemp[,3]))/
	sqrt((1-cor(dataTemp[,1],dataTemp[,3])^2)*(1-cor(dataTemp[,2],dataTemp[,3])^2))
	
	return(res)
})

partial_y<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])-lineage_velocity_y_avg)/lineage_velocity_y_sd
	M_collapse<-c(lineage_velocity_emb[1:9,])
	D_collapse<-c(lineage_velocity_emb[2:10,])
	
	lineage_start_pos_emb<-matrix(start_pos[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353],nrow=dim(lineage_info)[1])
	MD_pos<-c(lineage_start_pos_emb[2:10,])
	
	dataTemp<-na.omit(cbind(D_collapse,M_collapse,MD_pos))
	res<-(cor(dataTemp[,1],dataTemp[,2])-cor(dataTemp[,1],dataTemp[,3])*cor(dataTemp[,2],dataTemp[,3]))/
	sqrt((1-cor(dataTemp[,1],dataTemp[,3])^2)*(1-cor(dataTemp[,2],dataTemp[,3])^2))
	
	return(res)
})

partial_z<-sapply(1:dim(velocity)[1],function(emb){
	print(emb)
	lineage_velocity_emb<-(matrix(velocity[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])-lineage_velocity_z_avg)/lineage_velocity_z_sd
	M_collapse<-c(lineage_velocity_emb[1:9,])
	D_collapse<-c(lineage_velocity_emb[2:10,])
	
	lineage_start_pos_emb<-matrix(start_pos[emb,match(paste(c(lineage_info),"ccl",sep="_"),colnames(ccl_raw_nonNA_scaling))+353*2],nrow=dim(lineage_info)[1])
	MD_pos<-c(lineage_start_pos_emb[2:10,])
	
	dataTemp<-na.omit(cbind(D_collapse,M_collapse,MD_pos))
	res<-(cor(dataTemp[,1],dataTemp[,2])-cor(dataTemp[,1],dataTemp[,3])*cor(dataTemp[,2],dataTemp[,3]))/
	sqrt((1-cor(dataTemp[,1],dataTemp[,3])^2)*(1-cor(dataTemp[,2],dataTemp[,3])^2))
	
	return(res)
})


corr_partial_xyz<-data.frame(
	pearson_R=c(c(self_corr_x[1,]),c(self_corr_y[1,]),c(self_corr_z[1,])),
	partial_R=c(partial_x,partial_y,partial_z),
	coordinate=c(rep("x",length(partial_x)),rep("y",length(partial_y)),rep("z",length(partial_z)))
)

ggplot(corr_partial_xyz,aes(x=partial_R,col=coordinate))+
geom_density()


ggplot(corr_partial_xyz,aes(x=pearson_R,y=partial_R,col=coordinate))+
geom_point(size=0.2,stroke=0.15)+
#geom_smooth(method="lm")+
geom_vline(xintercept=0,lty=2)+
geom_hline(yintercept=0,lty=2)+
scale_x_continuous(breaks=c(-0.5,0))+
scale_y_continuous(breaks=c(-0.5,0))+
labs(x="Pearson's R",y="Partial R",title="Partial correlation")+
theme(legend.position=c(0.14,0.85))

ggsave("E:/canalizationENV/figs/physical_principle_lineage_patterns_partial_correlation.pdf",width=unit(2,"in"),height=unit(2.1,"in"))


save.image("velocity.RData")



#######################################################################################
# END
#######################################################################################
