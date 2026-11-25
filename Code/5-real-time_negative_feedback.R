




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

library(foreach)
library(glmnet)
library(parallel)
library(doParallel)



focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/realTime"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q11_cell350_alignment_xyz_noise.RData"
	),projectFolderPath=focal_projectFolderPath)

###

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

###

q1_cell350_alignment_xyz_noise_x<-q1_cell350_alignment_xyz_noise[,1:714]
q1_cell350_alignment_xyz_noise_y<-q1_cell350_alignment_xyz_noise[,(1+714):(714*2)]
q1_cell350_alignment_xyz_noise_z<-q1_cell350_alignment_xyz_noise[,(1+714*2):(714*3)]


q11_cell350_alignment_xyz_noise_x<-q11_cell350_alignment_xyz_noise[,1:714]
q11_cell350_alignment_xyz_noise_y<-q11_cell350_alignment_xyz_noise[,(1+714):(714*2)]
q11_cell350_alignment_xyz_noise_z<-q11_cell350_alignment_xyz_noise[,(1+714*2):(714*3)]



###############################################################################################
							M1-D1
###############################################################################################
	
##############
	x
##############

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_x_model_mother_coef_M1D1<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_x_model_mother_coef_M1D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mother_coef_M1D1,file="xyz_x_model_mother_coef_M1D1.RData")


xyz_x_model_mo_da1_R_M1D1<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_x_model_mo_da1_R_M1D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mo_da1_R_M1D1,file="xyz_x_model_mo_da1_R_M1D1.RData")


xyz_x_model_mo_da1_pvalue_M1D1<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_x_model_mo_da1_pvalue_M1D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mo_da1_pvalue_M1D1,file="xyz_x_model_mo_da1_pvalue_M1D1.RData")


##########
    y
##########

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_y_model_mother_coef_M1D1<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_y_model_mother_coef_M1D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mother_coef_M1D1,file="xyz_y_model_mother_coef_M1D1.RData")


xyz_y_model_mo_da1_R_M1D1<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_y_model_mo_da1_R_M1D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mo_da1_R_M1D1,file="xyz_y_model_mo_da1_R_M1D1.RData")


xyz_y_model_mo_da1_pvalue_M1D1<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_y_model_mo_da1_pvalue_M1D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mo_da1_pvalue_M1D1,file="xyz_y_model_mo_da1_pvalue_M1D1.RData")



##########
    z
##########

X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_z_model_mother_coef_M1D1<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_z_model_mother_coef_M1D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mother_coef_M1D1,file="xyz_z_model_mother_coef_M1D1.RData")


xyz_z_model_mo_da1_R_M1D1<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_z_model_mo_da1_R_M1D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mo_da1_R_M1D1,file="xyz_z_model_mo_da1_R_M1D1.RData")


xyz_z_model_mo_da1_pvalue_M1D1<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_z_model_mo_da1_pvalue_M1D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mo_da1_pvalue_M1D1,file="xyz_z_model_mo_da1_pvalue_M1D1.RData")




###############################################################################################
							M1-D11
###############################################################################################
	
##############
	x
##############

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q11_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_x_model_mother_coef_M1D11<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_x_model_mother_coef_M1D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mother_coef_M1D11,file="xyz_x_model_mother_coef_M1D11.RData")


xyz_x_model_mo_da1_R_M1D11<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_x_model_mo_da1_R_M1D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mo_da1_R_M1D11,file="xyz_x_model_mo_da1_R_M1D11.RData")


xyz_x_model_mo_da1_pvalue_M1D11<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_x_model_mo_da1_pvalue_M1D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mo_da1_pvalue_M1D11,file="xyz_x_model_mo_da1_pvalue_M1D11.RData")



##########
    y
##########

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q11_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_y_model_mother_coef_M1D11<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_y_model_mother_coef_M1D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mother_coef_M1D11,file="xyz_y_model_mother_coef_M1D11.RData")


xyz_y_model_mo_da1_R_M1D11<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_y_model_mo_da1_R_M1D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mo_da1_R_M1D11,file="xyz_y_model_mo_da1_R_M1D11.RData")


xyz_y_model_mo_da1_pvalue_M1D11<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_y_model_mo_da1_pvalue_M1D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mo_da1_pvalue_M1D11,file="xyz_y_model_mo_da1_pvalue_M1D11.RData")




##########
    z
##########

X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q11_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_z_model_mother_coef_M1D11<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_z_model_mother_coef_M1D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mother_coef_M1D11,file="xyz_z_model_mother_coef_M1D11.RData")


xyz_z_model_mo_da1_R_M1D11<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_z_model_mo_da1_R_M1D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mo_da1_R_M1D11,file="xyz_z_model_mo_da1_R_M1D11.RData")


xyz_z_model_mo_da1_pvalue_M1D11<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_z_model_mo_da1_pvalue_M1D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mo_da1_pvalue_M1D11,file="xyz_z_model_mo_da1_pvalue_M1D11.RData")




###############################################################################################
							M11-D11
###############################################################################################
	
##############
	x
##############

X_mo<-q11_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q11_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_x_model_mother_coef_M11D11<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_x_model_mother_coef_M11D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mother_coef_M11D11,file="xyz_x_model_mother_coef_M11D11.RData")


xyz_x_model_mo_da1_R_M11D11<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_x_model_mo_da1_R_M11D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mo_da1_R_M11D11,file="xyz_x_model_mo_da1_R_M11D11.RData")


xyz_x_model_mo_da1_pvalue_M11D11<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_x_model_mo_da1_pvalue_M11D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mo_da1_pvalue_M11D11,file="xyz_x_model_mo_da1_pvalue_M11D11.RData")


##########
    y
##########

X_mo<-q11_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q11_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_y_model_mother_coef_M11D11<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_y_model_mother_coef_M11D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mother_coef_M11D11,file="xyz_y_model_mother_coef_M11D11.RData")


xyz_y_model_mo_da1_R_M11D11<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_y_model_mo_da1_R_M11D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mo_da1_R_M11D11,file="xyz_y_model_mo_da1_R_M11D11.RData")


xyz_y_model_mo_da1_pvalue_M11D11<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_y_model_mo_da1_pvalue_M11D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mo_da1_pvalue_M11D11,file="xyz_y_model_mo_da1_pvalue_M11D11.RData")


##########
    z
##########

X_mo<-q11_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q11_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_z_model_mother_coef_M11D11<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_z_model_mother_coef_M11D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mother_coef_M11D11,file="xyz_z_model_mother_coef_M11D11.RData")


xyz_z_model_mo_da1_R_M11D11<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_z_model_mo_da1_R_M11D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mo_da1_R_M11D11,file="xyz_z_model_mo_da1_R_M11D11.RData")


xyz_z_model_mo_da1_pvalue_M11D11<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_z_model_mo_da1_pvalue_M11D11)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mo_da1_pvalue_M11D11,file="xyz_z_model_mo_da1_pvalue_M11D11.RData")





###############################################################################################
							M11-D1
###############################################################################################
	
##############
	x
##############

X_mo<-q11_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_x_model_mother_coef_M11D1<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_x_model_mother_coef_M11D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mother_coef_M11D1,file="xyz_x_model_mother_coef_M11D1.RData")


xyz_x_model_mo_da1_R_M11D1<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_x_model_mo_da1_R_M11D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mo_da1_R_M11D1,file="xyz_x_model_mo_da1_R_M11D1.RData")


xyz_x_model_mo_da1_pvalue_M11D1<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_x_model_mo_da1_pvalue_M11D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mo_da1_pvalue_M11D1,file="xyz_x_model_mo_da1_pvalue_M11D1.RData")

##########
    y
##########

X_mo<-q11_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_y_model_mother_coef_M11D1<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_y_model_mother_coef_M11D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mother_coef_M11D1,file="xyz_y_model_mother_coef_M11D1.RData")


xyz_y_model_mo_da1_R_M11D1<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_y_model_mo_da1_R_M11D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mo_da1_R_M11D1,file="xyz_y_model_mo_da1_R_M11D1.RData")


xyz_y_model_mo_da1_pvalue_M11D1<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_y_model_mo_da1_pvalue_M11D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mo_da1_pvalue_M11D1,file="xyz_y_model_mo_da1_pvalue_M11D1.RData")



##########
    z
##########

X_mo<-q11_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo

#

xyz_z_model_mother_coef_M11D1<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
names(xyz_z_model_mother_coef_M11D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mother_coef_M11D1,file="xyz_z_model_mother_coef_M11D1.RData")


xyz_z_model_mo_da1_R_M11D1<-sapply(1:dim(X_mo)[2],function(i){cor(X_da1[,i],X_mo[,i])})
names(xyz_z_model_mo_da1_R_M11D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mo_da1_R_M11D1,file="xyz_z_model_mo_da1_R_M11D1.RData")


xyz_z_model_mo_da1_pvalue_M11D1<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_da1[,i],X_mo[,i])$p.value})
names(xyz_z_model_mo_da1_pvalue_M11D1)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mo_da1_pvalue_M11D1,file="xyz_z_model_mo_da1_pvalue_M11D1.RData")




