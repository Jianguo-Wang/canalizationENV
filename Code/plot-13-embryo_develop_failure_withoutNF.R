



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
	"F:/wangj/Desktop/temp/canalizationENV/functions/RowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info_short.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_stage_df.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/coo_list_q11_cell350_translation_rotation_scaling.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData"
),projectFolderPath=focal_projectFolderPath)


#

q1_cell350_alignment_xyz_noise_x<-q1_cell350_alignment_xyz_noise[,1:714]
q1_cell350_alignment_xyz_noise_y<-q1_cell350_alignment_xyz_noise[,(1+714):(714*2)]
q1_cell350_alignment_xyz_noise_z<-q1_cell350_alignment_xyz_noise[,(1+714*2):(714*3)]


q1_cell350_alignment_xyz_x<-q1_cell350_alignment_xyz[,1:714]
q1_cell350_alignment_xyz_y<-q1_cell350_alignment_xyz[,(1+714):(714*2)]
q1_cell350_alignment_xyz_z<-q1_cell350_alignment_xyz[,(1+714*2):(714*3)]


cell_name<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))
cell_gen<-sapply(match(cell_name,lineage_info),function(x){rowColRec(x,dim(lineage_info)[1],dim(lineage_info)[2])[1]})





############################################## 
               without negative feedback
##############################################

### x 

x_without_nf<-q1_cell350_alignment_xyz_noise_x

for(i in 5:10){
	focal_da_s<-cell_name[cell_gen==i]
	for(j in 1:length(focal_da_s)){
		
		focal_da<-focal_da_s[j]
		focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da),1]
		
		focal_mo_ind<-match(focal_mo,cell_name)
		focal_da_ind<-match(focal_da,cell_name)
		
		set.seed(1)
		x_without_nf[,focal_da_ind]<-x_without_nf[,focal_mo_ind]+residuals(lm(q1_cell350_alignment_xyz_noise_x[,focal_da_ind]~q1_cell350_alignment_xyz_noise_x[,focal_mo_ind]))
		
	}
	
}

x_without_nf_pos<-x_without_nf

for(i in 1:714){
	x_without_nf_pos[,i]<-x_without_nf[,i]+mean(q1_cell350_alignment_xyz_x[,i])
}


### y 

y_without_nf<-q1_cell350_alignment_xyz_noise_y

for(i in 5:10){
	focal_da_s<-cell_name[cell_gen==i]
	for(j in 1:length(focal_da_s)){
		
		focal_da<-focal_da_s[j]
		focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da),1]
		
		focal_mo_ind<-match(focal_mo,cell_name)
		focal_da_ind<-match(focal_da,cell_name)
		
		set.seed(1)
		y_without_nf[,focal_da_ind]<-y_without_nf[,focal_mo_ind]+residuals(lm(q1_cell350_alignment_xyz_noise_y[,focal_da_ind]~q1_cell350_alignment_xyz_noise_y[,focal_mo_ind]))
		
	}
	
}

y_without_nf_pos<-y_without_nf

for(i in 1:714){
	y_without_nf_pos[,i]<-y_without_nf[,i]+mean(q1_cell350_alignment_xyz_y[,i])
}


### z  

z_without_nf<-q1_cell350_alignment_xyz_noise_z

for(i in 5:10){
	focal_da_s<-cell_name[cell_gen==i]
	for(j in 1:length(focal_da_s)){
		
		focal_da<-focal_da_s[j]
		focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da),1]
		
		focal_mo_ind<-match(focal_mo,cell_name)
		focal_da_ind<-match(focal_da,cell_name)
		
		set.seed(1)
		z_without_nf[,focal_da_ind]<-z_without_nf[,focal_mo_ind]+residuals(lm(q1_cell350_alignment_xyz_noise_z[,focal_da_ind]~q1_cell350_alignment_xyz_noise_z[,focal_mo_ind]))
		
	}
	
}


z_without_nf_pos<-z_without_nf

for(i in 1:714){
	z_without_nf_pos[,i]<-z_without_nf[,i]+mean(q1_cell350_alignment_xyz_z[,i])
}


#################################################################################################3



### without nf 
j=1
cell_stage_df[,1]<-x_without_nf_pos[j,]
cell_stage_df[,2]<-y_without_nf_pos[j,]
cell_stage_df[,3]<-z_without_nf_pos[j,]

dataTemp<-cell_stage_df


i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("withoutNF_develop_failure_contr_emb1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))



### with nf 
j=1

cell_stage_df[,1]<-q1_cell350_alignment_xyz[j,1:714]
cell_stage_df[,2]<-q1_cell350_alignment_xyz[j,c(1:714)+714]
cell_stage_df[,3]<-q1_cell350_alignment_xyz[j,c(1:714)+714*2]

dataTemp<-cell_stage_df


i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("withoutNF_develop_failure_contr_emb1_normal.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


### RNAi

### without nf 
j=106

cell_stage_df[,1]<-x_without_nf_pos[j,]
cell_stage_df[,2]<-y_without_nf_pos[j,]
cell_stage_df[,3]<-z_without_nf_pos[j,]

dataTemp<-cell_stage_df

i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("withoutNF_develop_failure_ABHD-11.2_emb1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))

### with nf 

j=106

cell_stage_df[,1]<-q1_cell350_alignment_xyz[j,1:714]
cell_stage_df[,2]<-q1_cell350_alignment_xyz[j,c(1:714)+714]
cell_stage_df[,3]<-q1_cell350_alignment_xyz[j,c(1:714)+714*2]

dataTemp<-cell_stage_df

i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("withoutNF_develop_failure_ABHD-11.2_emb1_normal.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))



### Env

### without nf 
j=1929

cell_stage_df[,1]<-x_without_nf_pos[j,]
cell_stage_df[,2]<-y_without_nf_pos[j,]
cell_stage_df[,3]<-z_without_nf_pos[j,]

dataTemp<-cell_stage_df

i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("withoutNF_develop_failure_XL_RW10434_16C_20230216_3_s1_emb1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))

### with nf 

j=1929

cell_stage_df[,1]<-q1_cell350_alignment_xyz[j,1:714]
cell_stage_df[,2]<-q1_cell350_alignment_xyz[j,c(1:714)+714]
cell_stage_df[,3]<-q1_cell350_alignment_xyz[j,c(1:714)+714*2]

dataTemp<-cell_stage_df

i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("withoutNF_develop_failure_XL_RW10434_16C_20230216_3_s1_emb1_normal.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))



