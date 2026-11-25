install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Ce.eg.db", "AnnotationDbi","ReactomePA","enrichplot","STRINGdb"))

library(ggplot2)

library(AnnotationDbi)
library(clusterProfiler)
library(org.Ce.eg.db)
library(enrichplot)

library(dplyr)
library(stringr)
library(ggplot2)


##########################################################################################

### colnames(ccl_raw_nonNA_scaling)
### velocity 

rownames(velocity)<-rownames(q1_cell350_alignment_xyz)


### gene name 

gene<-do.call(rbind,strsplit(rownames(q1_cell350_alignment_xyz)[1:1928],"_"))[,1]
env<-c(rep(c(c("16C","18C","20C","22C","24C"),
	         c("BW25113","DA1877"),
			 c("Fe1mM","Fe2mM","Fe4mM","Fe6mM","Fe8mM"),
			 c("GLU100","GLU200","GLU400","GLU50","GLU600"),
			 c("HB101","HS10m","HS15m","HS20m","HS25m","HS30m"),
			 c("HT115","NaCl100","NaCl20","NaCl40")),each=3),
	   rep("NaCl60",2),
	   rep("NaCl80",3),
	   rep("OP50",2),
	   rep(c("PQ1mM","PQ2mM"),each=3),
	   c("PQ4mM","PQ6mM"),
	   rep(c("SA100","SA200","SA400","SA50","SA600"),each=3))

gene_extend<-c(gene,env)


### average by gene
## calculate average 
velocity_average_by_gene<-do.call(rbind,lapply(unique(gene_extend),function(x){
	print(x)
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(velocity[geneInd,])
	}else{
		return(apply(velocity[geneInd,],2,mean,na.rm=T))	
	}
}))
rownames(velocity_average_by_gene)<-unique(gene_extend)
colnames(velocity_average_by_gene)<-colnames(velocity)

## calculate average difference
velocity_average_by_gene_relative<-do.call(rbind,lapply(1:dim(velocity_average_by_gene)[1],function(i){
	velocity_average_by_gene[i,]-velocity_average_by_gene[1,]
}))
rownames(velocity_average_by_gene_relative)<-unique(gene_extend)
colnames(velocity_average_by_gene_relative)<-colnames(velocity)

## calculate p-value  
Z_test<-function(value,vec){
	vec_mean<-mean(vec,na.rm=T)
	vec_sd<-sd(vec)
	accum_P<-pnorm(value,vec_mean,vec_sd)
	pvalue<-ifelse(accum_P<0.5,2*accum_P,2*(1-accum_P))
	return(pvalue)
}

velocity_average_by_gene_pvalue<-do.call(rbind,lapply(unique(gene_extend),function(x){
	print(x)
	print(match(x,unique(gene_extend)))
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(sapply(1:dim(velocity)[2],function(i){
			if(is.na(velocity[geneInd,i])){
				return(NA)
			}else{
				return(Z_test(velocity[geneInd,i],velocity[1:105,i]))
			}
		}))
	}else{
		return(sapply(1:dim(velocity)[2],function(i){
			if(length(velocity[geneInd,i])-length(which(is.na(velocity[geneInd,i])))>1){
				return(t.test(velocity[geneInd,i],velocity[1:105,i])$p.value)
			}
			else if(length(velocity[geneInd,i])-length(which(is.na(velocity[geneInd,i])))==1){
				return(Z_test(velocity[geneInd,i],velocity[1:105,i]))
			}else{
				return(NA)
			}
		}))
	}	
}))
rownames(velocity_average_by_gene_pvalue)<-unique(gene_extend)
colnames(velocity_average_by_gene_pvalue)<-colnames(velocity)

## calculate adjusted p-value 
velocity_average_by_gene_padj<-matrix(p.adjust(velocity_average_by_gene_pvalue,method="BH"),nrow=dim(velocity_average_by_gene_pvalue)[1])
velocity_average_by_gene_padj[which(velocity_average_by_gene_padj==0)]<-min(setdiff(velocity_average_by_gene_padj,0),na.rm=T)
rownames(velocity_average_by_gene_padj)<-unique(gene_extend)
colnames(velocity_average_by_gene_padj)<-colnames(velocity)


## calculate significant counts 
velocity_average_significant_count<-data.frame(
	gene=unique(gene_extend),
	sig_downreg=sapply(1:dim(velocity_average_by_gene)[1],function(i){length(which(velocity_average_by_gene_relative[i,]<0&velocity_average_by_gene_padj[i,]<0.05))}),
	sig_upreg=sapply(1:dim(velocity_average_by_gene)[1],function(i){length(which(velocity_average_by_gene_relative[i,]>0&velocity_average_by_gene_padj[i,]<0.05))}),
	sig_count=sapply(1:dim(velocity_average_by_gene)[1],function(i){length(which(velocity_average_by_gene_padj[i,]<0.05))})
)
velocity_average_significant_count$sig_diff_count=velocity_average_significant_count$sig_downreg-velocity_average_significant_count$sig_upreg

velocity_average_significant_count$sig_upreg_pvalue=sapply(1:dim(velocity_average_significant_count)[1],function(i){
	if(velocity_average_significant_count$sig_count[i]>0){
		return(binom.test(velocity_average_significant_count$sig_upreg[i],velocity_average_significant_count$sig_count[i],p=0.5,alternative="greater")$p.value)	
	}else{
		return(NA)
	}
})
velocity_average_significant_count$sig_downreg_pvalue=sapply(1:dim(velocity_average_significant_count)[1],function(i){
	if(velocity_average_significant_count$sig_count[i]>0){	
		return(binom.test(velocity_average_significant_count$sig_downreg[i],velocity_average_significant_count$sig_count[i],p=0.5,alternative="greater")$p.value)
	}else{
		return(NA)
	}
})


plot(velocity_average_significant_count[,2:3])
abline(0,1)

hist(velocity_average_significant_count[,5])

write.csv(velocity_average_significant_count,"velocity_average_significant_count.csv")



### variance by gene 
## calculate variance 
velocity_var_by_gene<-do.call(rbind,lapply(unique(gene_extend),function(x){
	print(x)
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(rep(NA,))
	}else{
		return(apply(velocity[geneInd,],2,var,na.rm=T))	
	}
}))
rownames(velocity_var_by_gene)<-unique(gene_extend)
colnames(velocity_var_by_gene)<-colnames(velocity)

## calculate log(var ratio)
velocity_var_by_gene_relative<-t(apply(velocity_var_by_gene,1,function(x){
	return(log((x+1e-12)/(velocity_var_by_gene[1,]+1e-12)))
}))
rownames(velocity_var_by_gene_relative)<-unique(gene_extend)
colnames(velocity_var_by_gene_relative)<-colnames(velocity)

## calculate p-value 
velocity_var_by_gene_pvalue<-do.call(rbind,lapply(unique(gene_extend),function(x){
	print(x)
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(rep(NA,dim(velocity)[2]))
	}else{
		return(sapply(1:dim(velocity)[2],function(i){
			if(length(velocity[geneInd,i])-length(which(is.na(velocity[geneInd,i])))>1){
				return(var.test(velocity[geneInd,i],velocity[1:105,i])$p.value)
			}else{
				return(NA)
			}			
		}))
	}	
}))
rownames(velocity_var_by_gene_pvalue)<-unique(gene_extend)
colnames(velocity_var_by_gene_pvalue)<-colnames(velocity)

## calculate adjusted p-value 
velocity_var_by_gene_padj<-matrix(p.adjust(velocity_var_by_gene_pvalue,method="BH"),nrow=dim(velocity_var_by_gene_pvalue)[1])
velocity_var_by_gene_padj[which(velocity_var_by_gene_padj==0)]<-min(setdiff(velocity_var_by_gene_padj,0),na.rm=T)
rownames(velocity_var_by_gene_padj)<-unique(gene_extend)
colnames(velocity_var_by_gene_padj)<-colnames(velocity)

## calculate significant counts 
velocity_var_significant_count<-data.frame(
	gene=unique(gene_extend),
	sig_downreg=sapply(1:dim(velocity_var_by_gene_relative)[1],function(i){length(which(velocity_var_by_gene_relative[i,]<0&velocity_var_by_gene_padj[i,]<0.05))}),
	sig_upreg=sapply(1:dim(velocity_var_by_gene_relative)[1],function(i){length(which(velocity_var_by_gene_relative[i,]>0&velocity_var_by_gene_padj[i,]<0.05))}),
	sig_count=sapply(1:dim(velocity_var_by_gene_relative)[1],function(i){length(which(velocity_var_by_gene_padj[i,]<0.05))})
)

velocity_var_significant_count$sig_diff_count=velocity_var_significant_count$sig_downreg-velocity_var_significant_count$sig_upreg

velocity_var_significant_count$sig_upreg_pvalue=sapply(1:dim(velocity_var_significant_count)[1],function(i){
	if(velocity_var_significant_count$sig_count[i]>0){
		return(binom.test(velocity_var_significant_count$sig_upreg[i],velocity_var_significant_count$sig_count[i],p=0.5,alternative="greater")$p.value)	
	}else{
		return(NA)
	}
})
velocity_var_significant_count$sig_downreg_pvalue=sapply(1:dim(velocity_var_significant_count)[1],function(i){
	if(velocity_var_significant_count$sig_count[i]>0){	
		return(binom.test(velocity_var_significant_count$sig_downreg[i],velocity_var_significant_count$sig_count[i],p=0.5,alternative="greater")$p.value)
	}else{
		return(NA)
	}
})

velocity_var_significant_count[velocity_var_significant_count$sig_count>0,]
plot(velocity_var_significant_count[,2:3])
abline(0,1)

hist(velocity_var_significant_count[,5])

write.csv(velocity_var_significant_count,"velocity_var_significant_count.csv")


### gene identification 
## sig genes 
sig_gene<-velocity_var_significant_count$gene[which(velocity_var_significant_count$sig_count>2)][1:51]
sig_gene_index<-match(sig_gene,rownames(velocity_var_by_gene))

sig_env<-velocity_var_significant_count$gene[which(velocity_var_significant_count$sig_count>2)][-c(1:51)]
sig_env_index<-match(sig_env,rownames(velocity_var_by_gene))

# insig genes 
#sig_gene_extend<-velocity_var_significant_count$gene[which(velocity_var_significant_count$sig_count>0)]
#set.seed(1)
#insig_gene_index<-sample(setdiff(2:641,
#								 c(match(sig_gene_extend,rownames(velocity_var_by_gene)),
#								   which(apply(velocity_var_by_gene_relative,1,function(x){length(which(!is.na(x)))})==0))),
#						 10)

# sample from velocity_var_by_gene_relative
#velocity_logRatio_pvalue<-data.frame(
#	gene=factor(tolower(rep(rownames(velocity_var_by_gene)[c(sig_gene_index,sig_env_index)],dim(velocity_var_by_gene)[2])),
#		level=tolower(rownames(velocity_var_by_gene)[c(sig_gene_index,sig_env_index)])),
#	cell_feature=factor(rep(colnames(velocity_var_by_gene),each=length(c(sig_gene_index,sig_env_index))),level=colnames(velocity_var_by_gene)),
#	values=as.numeric(velocity_var_by_gene_relative[c(sig_gene_index,sig_env_index),]),
#	pvalue=as.numeric(velocity_var_by_gene_pvalue[c(sig_gene_index,sig_env_index),]),
#	padj=as.numeric(velocity_var_by_gene_padj[c(sig_gene_index,sig_env_index),])
#)

velocity_logRatio_pvalue<-data.frame(
	gene=factor(tolower(rep(rownames(velocity_var_by_gene)[c(sig_gene_index)],dim(velocity_var_by_gene)[2])),
		level=tolower(rownames(velocity_var_by_gene)[c(sig_gene_index)])),
	cell_feature=factor(rep(colnames(velocity_var_by_gene),each=length(c(sig_gene_index))),level=colnames(velocity_var_by_gene)),
	values=as.numeric(velocity_var_by_gene_relative[c(sig_gene_index),]),
	pvalue=as.numeric(velocity_var_by_gene_pvalue[c(sig_gene_index),]),
	padj=as.numeric(velocity_var_by_gene_padj[c(sig_gene_index),])
)
velocity_logRatio_pvalue$values_adjust<-ifelse(velocity_logRatio_pvalue$values>0,velocity_logRatio_pvalue$values*3,velocity_logRatio_pvalue$values)

# plot 
ggplot(data = na.omit(velocity_logRatio_pvalue), 
       aes(x = gene, y = cell_feature, fill = values_adjust)) +
  geom_tile(color = NA) + 
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "log(variance ratio)"
  ) +
  geom_vline(xintercept = length(sig_gene) + 0.5, linetype = 2) +
  labs(
    x = "Gene (n = 51) & Env (n = 8)",
    y = "Cell",
    title = "log(variance ratio) of velocity"
  ) +
  theme(
    # remove y-axis text, ticks, and extra spacing
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    # remove panel spacing to eliminate left gap
    panel.spacing = unit(0, "mm"),
    # x-axis at top
    axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 4.8),
    axis.title.x.top = element_text(),
    axis.text.x.bottom = element_blank(),
    axis.title.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    # add border — now uniform
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
  ) +
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "left")


ggsave("E:/canalizationENV/figs/velocity_variance_sigGene.pdf",width=unit(4.2,"in"),height=unit(2.5,"in"))




### GO enrichment analysis 

## GO analysis 
sig_gene_all<-velocity_var_significant_count$gene[which(velocity_var_significant_count$sig_count>0)][1:169]

gene_mapped <- mapIds(
  org.Ce.eg.db,
  keys = tolower(sig_gene_all),
  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
  column = "ENTREZID"
)
gene_mapped

gene_mapped<-apply(data.frame(lower=mapIds(
			  org.Ce.eg.db,
			  keys = tolower(sig_gene_all),
			  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
			  column = "ENTREZID"
			),
			upper=mapIds(
			  org.Ce.eg.db,
			  keys = sig_gene_all,
			  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
			  column = "ENTREZID"
			)),1,function(x){ifelse(is.na(x[1]),x[2],x[1])})

names(which(is.na(gene_mapped)))


temp<-list(
c("b0261.1","vep-1","B0261.1"),
c("c30h7.2","erp-44.1","C30H7.2"),
c("cco-1","cox-5B","lpd-4","F26E4.9"),
c("eif-3.c","eif-3.C","T23D8.4"),
c("eif-3.e","eif-3.E","B0511.10"),
c("f26b1.2","hrpk-1","F26B1.2"),
c("f31c3.3","bltp-2","F31C3.3","fmp-27"),
c("f43g9.12","F43G9.12","mog-7"),
c("f52a8.6","kbrl-1","F52A8.6"),
c("f57b10.5","tmed-3","F57B10.5"),
c("inst-1","ints-1","let-384", "C06A5.1"),
c("m01a10.1","fubl-4","M01A10.1", "fubp-3.4"), 
c("pqn-44","tent-5","F55A12.9"),
c("t01h8.2","serr-1","T01H8.2"),
c("tag-179","algn-10","T24D1.4"),
c("w04a4.5","ints-4","W04A4.5"),
c("w09c5.8","cox-4","W09C5.8"),
c("y105e8a.3","slc-30A5","Y105E8A.d", "Y105E8A.3"),
c("y23h5b.6","ddx-10","Y74C10A_153.a", "Y23H5B.6", "Y23H5B.c"),
c("y44e3a.6","edc-4","Y44E3A.6"),
c("y48g1c.10","mtm-10","Y48G1C.f", "Y48G1C.10"),
c("y95b8a.8","zfr-1","Y95B8A.8")
)

temp_mapped<-sapply(1:length(temp),function(i){
	print(i)
	x=temp[[i]]

	return(na.omit(mapIds(org.Ce.eg.db,
	  keys = x,
	  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
	  column = "ENTREZID"
	)))	

})

gene_mapped[which(is.na(gene_mapped))]<-temp_mapped

genes_mapping=data.frame(gene=sig_gene_all,gene_mapped=gene_mapped)

selected_genes<-velocity_var_significant_count$gene[which(velocity_var_significant_count$sig_count>2)][1:51]

gene_mapped_selected <- na.omit(genes_mapping$gene_mapped[match(selected_genes,genes_mapping$gene)])

# The results of "ALL" mode are equal to the combination of "BP", "CC" and "MF" in this analysis
ego <- enrichGO(
  gene = gene_mapped_selected,
  OrgDb = org.Ce.eg.db,
  ont = "ALL",  # or BP, CC, MF
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego)
ego[ego$qvalue<0.05,]


ego_BP <- enrichGO(
  gene = gene_mapped_selected,
  OrgDb = org.Ce.eg.db,
  ont = "BP",  # or BP, CC, MF
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)
head(ego_BP)
ego_BP[ego_BP$qvalue<0.05,]


ego_MF <- enrichGO(
  gene = gene_mapped_selected,
  OrgDb = org.Ce.eg.db,
  ont = "MF",  # or BP, CC, MF
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego_MF)
ego_MF[ego_MF$qvalue<0.05,]


ego_CC <- enrichGO(
  gene = gene_mapped_selected,
  OrgDb = org.Ce.eg.db,
  ont = "CC",  # or BP, CC, MF
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego_CC)
ego_CC[ego_CC$qvalue<0.05,]

### KEGG enrichment analysis 

# kegg pathway 

library(clusterProfiler)
library(org.Ce.eg.db)
							
gene_mapped_selected_kegg<-data.frame(
	bitr_kegg(
		gene_mapped_selected,
		fromType = "ncbi-geneid",
		toType = "kegg",
		organism = "cel"
	))
gene_mapped_selected_kegg<-data.frame(
	gene_mapped_selected_kegg,
	ENTREZID=gene_mapped_selected,
	SYMBOL=selected_genes
)


ekegg <- enrichKEGG(
  gene = gene_mapped_selected_kegg[,2],
  organism = "cel",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

head(ekegg)
gene_mapped_selected_kegg[match(unlist(strsplit(ekegg$geneID,"/")),gene_mapped_selected_kegg[,2]),]

names(gene_mapped[match(gene_mapped_left,gene_mapped)])


### plot GO & KEGG results
## all combined results by dotplot 
df_go <- data.frame(ego)
df_go$Source <- "GO"

df_kegg <- data.frame(ekegg)
df_kegg$Source <- "KEGG"

df_kegg$ONTOLOGY <- "KEGG"

common_cols <- c("ONTOLOGY","ID","Description","GeneRatio","BgRatio",
                 "pvalue","p.adjust","qvalue","geneID","Count","Source")

df_go  <- df_go[,  common_cols]
df_kegg <- df_kegg[, common_cols]

dataTemp <- rbind(df_go, df_kegg)

calc_ratio <- function(x){
  v <- as.numeric(unlist(strsplit(x, "/")))
  v[1] / v[2]
}

dataTemp$Description <- factor(dataTemp$Description,
                               levels = rev(dataTemp$Description))
dataTemp$GeneRatio <- sapply(dataTemp$GeneRatio, calc_ratio)

dataTemp$geneID[19]<-paste(tolower(gene_mapped_selected_kegg[match(unlist(strsplit(dataTemp$geneID[19],"/")),gene_mapped_selected_kegg[,2]),4]),collapse="/")

# finally used 
ggplot(dataTemp, aes(x = GeneRatio, y = Description, col = p.adjust, size = Count)) +
  geom_point() +
  scale_colour_gradient(low = "#F8766D", high = "#C77CFF") +
  labs(x = "Gene ratio", y = "", title = "") +
  theme_get() +
  theme(
    # --- y-axis text ---
    axis.text.y = element_text(angle = 0, hjust = 1, size = rel(8/11)),
    # --- legend ---
    legend.title = element_text(size = rel(8/11)),
    # --- add ALL FOUR borders ---
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    # --- add horizontal grid lines only ---
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.2),
    # turn off vertical grid lines
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
	legend.position="bottom"
  ) +
  scale_x_continuous(breaks = c(0.05, 0.2)) +
  scale_size(range = c(0.5, 2))

ggsave("E:/canalizationENV/figs/velocity_variance_sigGene_GO_kegg_dotplot_all.pdf",width=unit(3.2,"in"),height=unit(2.7,"in"))

## plot gene-GO matrix 
# 1. create gene-GO dataframe
gene_GO_sparse<-data.frame(
	GO=unlist(lapply(1:dim(dataTemp)[1],function(i){rep(dataTemp$Description[i],length(unlist(strsplit(dataTemp$geneID[i],"/"))))})),
	Gene=tolower(unlist(lapply(1:dim(dataTemp)[1],function(i){unlist(strsplit(dataTemp$geneID[i],"/"))})))
)
gene_GO_sparse$value<-rep(1,dim(gene_GO_sparse)[1])
gene_GO_sparse$Gene[which(is.na(match(gene_GO_sparse$Gene,tolower(sig_gene))))]<-c("f26b1.2","inst-1","w04a4.5")

gene_GO_sparse$Gene<-factor(gene_GO_sparse$Gene,level=tolower(sig_gene))

gene_GO_sparse$GO<-factor(gene_GO_sparse$GO,level=rev(dataTemp$Description))

# 2. Expand the data to include ALL Gene × GO combinations
gene_GO_full <- merge(
  expand.grid(
    Gene = levels(gene_GO_sparse$Gene),
    GO   = levels(gene_GO_sparse$GO)
  ),
  gene_GO_sparse,
  by = c("Gene", "GO"),
  all.x = TRUE
)

# 3. Plot (empty cells will be NA → shown as blank)
ggplot(gene_GO_full, aes(x = Gene, y = GO, color = value)) +
  geom_point(show.legend = FALSE, na.rm = FALSE) +
  scale_color_gradient(low = "blue", high = "red", na.value = NA) +
  
  # Move x-axis to the top
  scale_x_discrete(position = "top") +
  
  labs(x = "Gene") +   # ← your x-axis label
  
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x.bottom = element_blank(),   # hide bottom axis text
    axis.ticks.x.bottom = element_blank(),  # hide bottom ticks
    axis.title.x.bottom = element_blank(),  # hide bottom axis label
    
    # show top ticks (optional; remove if not needed)
    axis.ticks.x.top = element_line(),

    # --- Add borders ---
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),

    # --- Add grid lines ---
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor.y = element_blank()
  )

ggsave("E:/canalizationENV/figs/velocity_variance_Gene_GO_relationship.pdf",width=unit(3.55,"in"),height=unit(1.8,"in"))

## total variance  comparison 
velocity_var_by_gene_varAverage<-apply(velocity_var_by_gene,1,mean,na.rm=T)

dataTemp<-data.frame(
	varAverage=velocity_var_by_gene_varAverage,
	type=rep("BgGene",dim(velocity_var_by_gene)[1])
)
dataTemp$type[1]<-"WT"
dataTemp$type[642:680]<-"BGEnv"
dataTemp$type[sig_gene_index]<-"Gene"
dataTemp$type[sig_env_index]<-"Env"

ggplot(dataTemp,aes(x=type,y=varAverage,col=type))+
geom_jitter(height=0,width=0.15)+
geom_hline(yintercept=dataTemp$varAverage[1],lty=2)

GO_enriched_genes<-names(table(as.character(gene_GO_sparse$Gene)))
GO_enriched_genes

dataTemp$GO_type<-rep("not",length(unique(gene_extend)))
dataTemp$GO_type[match(GO_enriched_genes,tolower(unique(gene_extend)))]<-"enriched"

ggplot(dataTemp[dataTemp$type=="Gene",],aes(x=varAverage,col=GO_type))+
geom_density()


## plot network by cnetplot 

df_go <- as.data.frame(ego)
df_go$Category <- "GO"

df_kegg <- as.data.frame(ekegg)
df_kegg$geneID<-paste(tolower(gene_mapped_selected_kegg[match(unlist(strsplit(ekegg$geneID,"/")),gene_mapped_selected_kegg[,2]),4]),collapse="/")

df_kegg$Category <- "KEGG"

df_kegg2 <- df_kegg[, c("ID","Description","GeneRatio","BgRatio","pvalue",
                        "p.adjust","qvalue","geneID","Count","Category")]

df_go2 <- df_go[, c("ID","Description","GeneRatio","BgRatio","pvalue",
                    "p.adjust","qvalue","geneID","Count","Category")]

df_all <- rbind(df_go2, df_kegg2)

library(clusterProfiler)

combined_er <- new(
  "enrichResult",
  result = df_all,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  organism = ego@organism,        # or "cel"
  ontology = "MIXED",             # GO + KEGG mixed
  gene = unique(unlist(strsplit(df_all$geneID, "/"))),
  keytype = ego@keytype
)

library(enrichplot)

# 
pdf("E:/canalizationENV/figs/velocity_variance_sigGene_GO_kegg_cnetplot_all.pdf",width=10,height=10)

set.seed(1) 
cnetplot(combined_er,cex_category=0.5,
         showCategory = 20,   # or any number you want
         circular = FALSE,
         colorEdge = FALSE)

dev.off()


### lineage variation of velcity for example velocity-regulation genes 





### normal distribution 
pdf("normal_distribution.pdf")
plot(density(rnorm(10000)))
dev.off()





# k forms the reverse force, increased k reduce force, therefore, causing increased variance 

# sisters inherite correlated force, followed by correlated velocity, phase deviation, and noise changes, therefore, similar k  


# gene expression 
geneExp<-readRDS("E:/canalization_and_divergence/raw/cel_cbr_embryo_single_cell-data/exp_vis copy/data/expr_data.rds")

unmatch_index<-which(sapply(sapply(sig_gene,function(x){grep(tolower(x),tolower(rownames(geneExp)[grep("C. elegans",rownames(geneExp))]))}),length)==0)
sig_gene_adj<-sig_gene 
sig_gene_adj[unmatch_index]<-c("hrpk-1","kbrl-1","ints-1","ints-4","cox-4","edc-4")


sigGene_withExp<-lapply(sig_gene_adj[match(tolower(gene_enriched),tolower(sig_gene))],function(x){grep(tolower(x),tolower(rownames(geneExp)[grep("C. elegans",rownames(geneExp))]))})
sigGene_geneExp<-t(geneExp[unlist(sigGene_withExp),1:303])
sigGene_geneExp<-sigGene_geneExp[,which(apply(sigGene_geneExp,2,var)>0)]

library(pheatmap)

pheatmap(cor(sigGene_geneExp,use="pairwise.complete.obs"))


#######################################################
# END 
#######################################################



pdf("E:/canalizationENV/figs/velocity_variance_sigGene_GO_analysis_cnetplot_BP.pdf",width=10,height=10)
set.seed(10)
cnetplot(ego,showCategory=ego$Description[1:9])

dev.off()



## for BP network  
pdf("E:/canalizationENV/figs/velocity_variance_sigGene_GO_analysis_cnetplot_BP.pdf",width=10,height=10)
set.seed(10)
cnetplot(ego,showCategory=ego$Description[1:9])

dev.off()



pdf("E:/canalizationENV/figs/velocity_variance_sigGene_GO_analysis_cnetplot_CC.pdf",width=10,height=10)
set.seed(10)

cnetplot(ego,showCategory=ego$Description[-c(1:9)])

dev.off()

















# various plots 

dotplot(ego, showCategory = 20)
barplot(ego, showCategory = 20)






pdf("E:/canalizationENV/figs/velocity_variance_sigGene_GO_analysis_dotplot_all.pdf",width=6.4,height=6)

dotplot(ego, showCategory = 20, label_format=100)+
facet_wrap(vars(ONTOLOGY),strip.position="right",nrow=2,scales="free_y")+
theme_get()+
theme(
	axis.text.y=element_text(angle=0,hjust=1,vjust=0.5),
	legend.title=element_text(size=rel(8/11))
)

dev.off()

calc_ratio <- function(x) {
  parts <- as.numeric(strsplit(x, "/")[[1]])
  parts[1] / parts[2]
}

dataTemp<-data.frame(ego)
#dataTemp$FC<-sapply(dataTemp$GeneRatio,calc_ratio)/sapply(dataTemp$BgRatio,calc_ratio)
dataTemp$Description<-factor(dataTemp$Description,level=rev(dataTemp$Description))
dataTemp$GeneRatio<-sapply(dataTemp$GeneRatio,calc_ratio)
ggplot(dataTemp, aes(x=GeneRatio,y=Description,col=p.adjust,size=Count))+
geom_point()+
scale_colour_gradient(
  low = "#F8766D",
  high = "#C77CFF",
)+
labs(x="Gene ratio",y="",title="")+
facet_wrap(vars(ONTOLOGY),strip.position="right",nrow=2,scales="free_y")+
theme_get()+
theme(
	axis.text.y=element_text(angle=0,hjust=1,vjust=0.5),
	legend.title=element_text(size=rel(8/11))
)


### ---- 3. Plot with wrapped y-axis text ----
# not used 
ggplot(dataTemp,
       aes(x = GeneRatio,
           y = Description,
           col = p.adjust,
           size = Count)) +
  geom_point() +
  scale_colour_gradient(low = "#F8766D", high = "#C77CFF") +
  labs(x = "Gene ratio", y = "", title = "") +
  facet_wrap(vars(ONTOLOGY),
             strip.position = "right",
             nrow = 3,
             scales = "free_y") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +   # ← auto line feed
  theme_get() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
    legend.title = element_text(size = rel(8/11))
  )+


pdf("E:/canalizationENV/figs/velocity_variance_sigGene_GO_analysis_cnetplot_all.pdf",width=10,height=10)
set.seed(10)

cnetplot(ego,showCategory=20)

dev.off()



pdf("E:/canalizationENV/figs/velocity_variance_sigGene_GO_analysis_dotplot_all.pdf",width=10,height=6, colormodel="cmyk")
dotplot(ego, showCategory = 20,font.size=15,label_format=60)
dev.off()



### left genes 

gene_enriched<-unique(unlist(lapply(ego$geneID,function(x){strsplit(x,"/")})))
gene_mapped_enriched<-unlist(mapIds(
			  org.Ce.eg.db,
			  keys = gene_enriched,
			  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
			  column = "ENTREZID"
			))

gene_mapped_left<-setdiff(gene_mapped_selected,gene_mapped_enriched)

ego_left <- enrichGO(
	gene = gene_mapped_left,
	OrgDb = org.Ce.eg.db,
	ont = "ALL",  # or BP, CC, MF
	pAdjustMethod = "BH",
	qvalueCutoff = 0.05,
	readable = TRUE
)


head(ego_left)
ego_left[ego_left$qvalue<0.05,]
cnetplot(ego_left)

##

match(gene_enriched,tolower(gene_extend))


###############################################################################
# pathway 
###############################################################################

dotplot(ekegg, showCategory = 20) +
  ggtitle("KEGG Pathway Enrichment - C. elegans")


# cnet plot 

# for kegg 
cnetplot(ekegg, showCategory = 10)

cnetplot(rbind(ego,ekegg), showCategory = 20)

#








#upsetplot(ego)


#write.csv(as.data.frame(ego),   "GO_results.csv")
#write.csv(as.data.frame(ekegg), "KEGG_results.csv")
#write.csv(as.data.frame(ereact),"Reactome_results.csv")


###############################################################################
# Phenotype enrichment — WormBase / WormCat / WormExp
###############################################################################

# WormBase enrichment: WormBase has phenotype annotations—there are APIs or downloads.

# WormCat: a tool tailored for C. elegans categorization (web / R wrappers).
# Use phenotype enrichment to connect genes to observable worm phenotypes (e.g., egg-laying, dauer).

















#######################################################
# END 
#######################################################

setwd("E:/canalizationENV/data/velocity")

source("E:/canalizationENV/functions/set_theme.R")

load("E:\\canalizationENV\\data\\basic\\q1_cell350_alignment_xyz.RData")

load("E:\\canalizationENV\\data\\basic\\q11_cell350_alignment_xyz.RData")
load("E:\\canalizationENV\\data\\basic\\cell_mo_da1_da2.RData")
load("E:\\canalizationENV\\data\\prediction\\embryo_lethal.RData")



{
bac_food<-c("BW25113","DA1877","HB101","HT115","OP50")
col1<-rep("other embryos",dim(q1_cell350_alignment_xyz)[1])
bac_food_index<-unlist(lapply(bac_food,function(x){which(gene_extend==x)}))
col1[bac_food_index]<-"bac_food"

dataTemp<-data.frame(XM=q1_cell350_alignment_xyz[,"Eala_y"],
					 XD=q1_cell350_alignment_xyz[,"Ealaa_y"],
					 delta_X=q1_cell350_alignment_xyz[,"Ealaa_y"]-
							 q1_cell350_alignment_xyz[,"Eala_y"],
					 col1=col1)
ggplot(dataTemp,aes(x=XM,y=XD,col=col1))+geom_point()+geom_smooth(method="lm")
ggplot(dataTemp,aes(x=XM,y=delta_X,col=col1))+geom_point()+geom_smooth(method="lm")
var.test(velocity[bac_food_index,"Eala_y"],velocity[1:105,"Eala_y"],alternative="less")
}



###

data.frame(gene=c("skr-1","skr-2","cul-1","chn-1"),sapply(c("skr-1","skr-2","cul-1","chn-1"),function(x){
	print(x)
	geneInd<-which(tolower(gene_extend)==x)
	return(length(geneInd))
}))





## refine data to remove NA lines and adjust those pvalue=0
velocity_var_by_gene_refine_RNAi<-velocity_var_by_gene_refine[2:615,]
velocity_var_by_gene_relative_refine_RNAi<-velocity_var_by_gene_relative_refine[2:615,]
velocity_var_by_gene_pvalue_refine_RNAi<-velocity_var_by_gene_pvalue_refine[2:615,]

velocity_var_by_gene_padj_RNAi<-matrix(p.adjust(velocity_var_by_gene_pvalue_refine_RNAi,method="BH"),nrow=dim(velocity_var_by_gene_pvalue_refine_RNAi)[1])
velocity_var_by_gene_padj_RNAi[which(velocity_var_by_gene_padj_RNAi==0)]<-min(setdiff(velocity_var_by_gene_padj_RNAi,0))
rownames(velocity_var_by_gene_padj_RNAi)<-rownames(velocity_var_by_gene_refine_RNAi)
colnames(velocity_var_by_gene_padj_RNAi)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]

velocity_var_by_gene_pvalue_refine_RNAi[which(velocity_var_by_gene_pvalue_refine_RNAi==0)]<-min(setdiff(velocity_var_by_gene_pvalue_refine_RNAi,0))
velocity_var_by_gene_pvalue_norm_RNAi<-qnorm(velocity_var_by_gene_pvalue_refine_RNAi)


velocity_var_significant_count_RNAi<-data.frame(gene=rownames(velocity_var_by_gene_refine_RNAi),
	sig_downreg=sapply(1:dim(velocity_var_by_gene_refine_RNAi)[1],function(i){length(which(velocity_var_by_gene_relative_refine_RNAi[i,]<0&velocity_var_by_gene_padj_RNAi[i,]<0.05))}),
	sig_upreg=sapply(1:dim(velocity_var_by_gene_refine_RNAi)[1],function(i){length(which(velocity_var_by_gene_relative_refine_RNAi[i,]>0&velocity_var_by_gene_padj_RNAi[i,]<0.05))}),
	sig_count=sapply(1:dim(velocity_var_by_gene_refine_RNAi)[1],function(i){length(which(velocity_var_by_gene_padj_RNAi[i,]<0.05))}))

write.csv(velocity_var_by_gene_pvalue_refine_RNAi,"velocity_var_by_gene_pvalue_refine_RNAi.csv")
write.csv(velocity_var_by_gene_refine_RNAi,"velocity_var_by_gene_refine_RNAi.csv")
write.csv(velocity_var_by_gene_relative_refine_RNAi,"velocity_var_by_gene_relative_refine_RNAi.csv")
write.csv(velocity_var_by_gene_padj_RNAi,"velocity_var_by_gene_padj_RNAi.csv")
write.csv(velocity_var_significant_count_RNAi,"velocity_var_significant_count_RNAi.csv")


## 
sig_gene<-rownames(velocity_var_by_gene_padj)[apply(velocity_var_by_gene_padj,1,function(x){length(which(x<0.05))})>0]
sig_gene_RNAi<-rownames(velocity_var_by_gene_padj_RNAi)[apply(velocity_var_by_gene_padj_RNAi,1,function(x){length(which(x<0.05))})>0]

#velocity_var_significant_count<-read.csv("E:/canalizationENV/data/velocity/velocity_var_significant_count.csv")
velocity_var_significant_count_sig<-velocity_var_significant_count[velocity_var_significant_count$sig_count>0,]
velocity_var_significant_count_sig[order(velocity_var_significant_count_sig$sig_count),]
#velocity_var_significant_count_RNAi<-read.csv("E:/canalizationENV/data/velocity/velocity_var_significant_count_RNAi.csv")
velocity_var_significant_count_RNAi_sig<-velocity_var_significant_count_RNAi[velocity_var_significant_count_RNAi$sig_count>0,]
velocity_var_significant_count_RNAi_sig[order(velocity_var_significant_count_RNAi_sig$sig_count),]



# distribution

velo_Var_statistics<-data.frame(gene=velocity_average_significant_count$gene[which(!is.na(velocity_var_by_gene[,1]))],
	do.call(rbind,lapply(1:dim(velocity_average_significant_count)[1],function(i){
		if (length(which(is.na(velocity_var_by_gene[i,])))==0){
			# average test 
			varRatio<-mean(log((velocity_var_by_gene[i,]+1e-8)/(velocity_var_by_gene[1,]+1e-8)))
			varRatio_pval<-wilcox.test(velocity_var_by_gene[i,],velocity_var_by_gene[1,],paired=T)$p.value
			less_pval<-wilcox.test(velocity_var_by_gene[i,],velocity_var_by_gene[1,],paired=T,alternative="less")$p.value
			greater_pval<-wilcox.test(velocity_var_by_gene[i,],velocity_var_by_gene[1,],paired=T,alternative="greater")$p.value
			return(c(varRatio,varRatio_pval,less_pval,greater_pval))		
		}

	})))

colnames(velo_Var_statistics)[2:5]<-c("varRatio","varRatio_pval","less_pval","greater_pval")

velo_Var_statistics$varRatio_padj<-p.adjust(velo_Var_statistics$varRatio_pval,"BH")
velo_Var_statistics$less_padj<-p.adjust(velo_Var_statistics$less_pval,"BH")
velo_Var_statistics$greater_padj<-p.adjust(velo_Var_statistics$greater_pval,"BH")
head(velo_Var_statistics)

sigGene_varRatio<-velo_Var_statistics$gene[which(velo_Var_statistics$varRatio_padj<0.05)]
sigGene_varRatio_less<-velo_Var_statistics$gene[which(velo_Var_statistics$less_padj<0.05)]
sigGene_varRatio_greater<-velo_Var_statistics$gene[which(velo_Var_statistics$greater_padj<0.05)]

dim(velo_Var_statistics[which(velo_Var_statistics$varRatio_padj<0.05),])
dim(velo_Var_statistics[which(velo_Var_statistics$less_padj<0.05),])
dim(velo_Var_statistics[which(velo_Var_statistics$greater_padj<0.05),])

par(mfrow=c(2,3))
hist(velo_Var_statistics[,"varRatio"],100)
hist(velo_Var_statistics[which(velo_Var_statistics$varRatio_padj<0.05),"varRatio"],100)
hist(velo_Var_statistics[which(velo_Var_statistics$varRatio_padj>0.05),"varRatio"],100)
hist(log(velocity_var_by_gene[1,]))
hist(log(velocity_var_by_gene[2,]/velocity_var_by_gene[1,]),100)
hist(log(velocity_var_by_gene[23,]/velocity_var_by_gene[1,]),100)

write.csv(velo_Var_statistics, "velo_Var_statistics.csv")

hist(velo_Var_statistics$varRatio[match(sig_gene,velo_Var_statistics$gene)])




library(foreach)
library(parallel)
library(doParallel)
cl <- makeCluster(28)
registerDoParallel(cl)

velocity_var_by_gene_pvalue_new_profile<-foreach(x=unique(gene_extend),.combine=c)%dopar%{
	print(x)
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(NA)
	}
	if(length(geneInd)>1){
		set.seed(1)
		rand_ctr<-sample(1:105,length(geneInd))	
		
		obs_stat<-sapply(1:dim(velocity)[2],function(i){var.test(velocity[geneInd,i],velocity[-rand_ctr,i])$p.value})
		
		null_stat<-sapply(1:dim(velocity)[2],function(i){var.test(velocity[rand_ctr,i],velocity[-rand_ctr,i])$p.value})
		
		return(wilcox.test(obs_stat, null_stat,alternative="less")$p.value)
	}
}

stopCluster(cl)

###############################################################################

# --------------------------------------------------------------------------
# robust_test: per-cell F-tests combined with Fisher and Stouffer; Wilcoxon
# Notes:
#  - F test per cell comparing s2_gi to s2_ctrl,i uses F = s2_gi / s2_ctrl,i
#  - For alternative="greater" we take upper-tail p = 1 - Fcdf(F, df1, df2)
#  - For two.sided we compute two-sided p by doubling the smaller tail
#  - Fisher's method: -2 sum log p_i ~ chi-square (2m df) under independence
#  - Stouffer: convert p->Z and sum Z / sqrt(m) ~ N(0,1)
# Limitations: these methods assume independence across cells (or weak dependence).
# --------------------------------------------------------------------------

# sampling can cause instability; profile-based test can be misled by minor effects; examine whether overal elevation of variance

robust_test <- function(S2, DF, ctrl_idx = 1, alternative = c("two.sided","greater","less")) {
	alternative <- match.arg(alternative)
	S2c <- as.matrix(S2)
	DFc <- as.matrix(DF)
	G <- nrow(S2c); m <- ncol(S2c)
	ctrl_row <- S2c[ctrl_idx, ]
	
	res <- data.frame(gene = rownames(S2c),
					fisher_p = NA_real_, stouffer_p = NA_real_, wilcox_p = NA_real_,
					stringsAsFactors = FALSE)
	
	for (g in 1:G) {
		print(g)
		if (g == ctrl_idx | any(is.na(S2c[g, ]))) next
		Fstat <- (S2c[g, ] + EPS) / (ctrl_row + EPS)
		
		if (alternative == "greater") {
			pvec <- 1 - pf(Fstat, df1 = DFc[g, ], df2 = DFc[ctrl_idx, ])
		} else if (alternative == "less") {
			pvec <- pf(Fstat, df1 = DFc[g, ], df2 = DFc[ctrl_idx, ])
		} else {
			# two-sided: double the smaller tail
			lower <- pf(Fstat, df1 = DFc[g, ], df2 = DFc[ctrl_idx, ])
			pvec <- 2 * pmin(lower, 1 - lower)
		}
		
		pvec[is.na(pvec)] <- 1
		pvec <- pmin(pmax(pvec, 1e-12), 1 - 1e-12)
		# Fisher
		chisq <- -2 * sum(log(pvec))
		fisher_p <- 1 - pchisq(chisq, df = 2 * m)
		
		# Stouffer (unweighted)
		zvec <- qnorm(1 - pvec)
		zst <- sum(zvec) / sqrt(m)
		stouffer_p <- 1 - pnorm(zst)
		
		# Wilcoxon signed-rank on log(Fstat) (paired across cells). This is robust
		logVarRatio <- log(Fstat)
		if (length(logVarRatio) < 2 || all(abs(logVarRatio) < 1e-12)) {
			w_p <- 1
		} else {
			alt_w <- if (alternative == "two.sided") "two.sided" else if (alternative == "greater") "greater" else "less"
			w_p <- wilcox.test(logVarRatio, mu = 0, alternative = alt_w, exact = FALSE)$p.value
		}
		
		res[g, c("fisher_p","stouffer_p","wilcox_p")] <- list(fisher_p, stouffer_p, w_p)
	}
	
	res$fisher_p_adj <- p.adjust(res$fisher_p, method = "BH")
	res$stouffer_p_adj <- p.adjust(res$stouffer_p, method = "BH")
	res$wilcox_p_adj <- p.adjust(res$wilcox_p, method = "BH")
	return(res)
}

# --------------------------------------------------------------------------
# run real data 
# --------------------------------------------------------------------------
EPS=1e-12

S2_velocity = velocity_var_by_gene
gene_extend_DF<-sapply(unique(gene_extend),function(x){
	print(x)
	geneInd<-which(gene_extend==x)	
	return(length(geneInd)-1)
})

DF_velocity = matrix(gene_extend_DF,nrow=dim(velocity_var_by_gene)[1],ncol=dim(velocity_var_by_gene)[2])

# for two.sided 
rob_res <- robust_test(S2=S2_velocity, DF=DF_velocity, ctrl_idx = 1, alternative = "two.sided")

sig_gene_fisher<-rob_res$gene[(rob_res$fisher_p_adj<0.05)&!is.na(rob_res$fisher_p_adj)]
sig_gene_stouffer<-rob_res$gene[(rob_res$stouffer_p_adj<0.05)&!is.na(rob_res$stouffer_p_adj)]
sig_gene_byCound<-velocity_var_significant_count$gene[which(velocity_var_significant_count[,4]>0)]

length(sig_gene_byCound)
length(sig_gene_fisher)
length(sig_gene_stouffer)
length(intersect(sig_gene_byCound,sig_gene_fisher))
length(intersect(sig_gene_byCound,sig_gene_stouffer))
length(intersect(sig_gene_fisher,sig_gene_stouffer))

write.csv(rob_res,"rob_res.csv")

# for less 
rob_res_less <- robust_test(S2=S2_velocity, DF=DF_velocity, ctrl_idx = 1, alternative = "less")

sig_gene_fisher<-rob_res_less$gene[(rob_res_less$fisher_p_adj<0.05)&!is.na(rob_res_less$fisher_p_adj)]
sig_gene_stouffer<-rob_res_less$gene[(rob_res_less$stouffer_p_adj<0.05)&!is.na(rob_res_less$stouffer_p_adj)]
sig_gene_byCound<-velocity_var_significant_count$gene[which(velocity_var_significant_count[,2]>velocity_var_significant_count[,3])]

length(sig_gene_byCound)
length(sig_gene_fisher)
length(sig_gene_stouffer)
length(intersect(sig_gene_byCound,sig_gene_fisher))
length(intersect(sig_gene_byCound,sig_gene_stouffer))
length(intersect(sig_gene_fisher,sig_gene_stouffer))

write.csv(rob_res_less,"rob_res_less.csv")


# for greater 
rob_res_greater <- robust_test(S2=S2_velocity, DF=DF_velocity, ctrl_idx = 1, alternative = "greater")

sig_gene_fisher<-rob_res_greater$gene[(rob_res_greater$fisher_p_adj<0.05)&!is.na(rob_res_greater$fisher_p_adj)]
sig_gene_stouffer<-rob_res_greater$gene[(rob_res_greater$stouffer_p_adj<0.05)&!is.na(rob_res_greater$stouffer_p_adj)]
sig_gene_byCound<-velocity_var_significant_count$gene[which(velocity_var_significant_count[,2]<velocity_var_significant_count[,3])]

length(sig_gene_byCound)
length(sig_gene_fisher)
length(sig_gene_stouffer)
length(intersect(sig_gene_byCound,sig_gene_fisher))
length(intersect(sig_gene_byCound,sig_gene_stouffer))
length(intersect(sig_gene_fisher,sig_gene_stouffer))

write.csv(rob_res_greater,"rob_res_greater.csv")

###############################################################################
# PPI analysis
###############################################################################

# ─────────────────────────────────────────────────────────
# Required libraries
# ─────────────────────────────────────────────────────────
library(igraph)
library(dplyr)
library(clusterProfiler)
library(org.Ce.eg.db)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggraph)


# ─────────────────────────────────────────────────────────
# Load PPI dataset
# ─────────────────────────────────────────────────────────
ppi <- read.csv("GGI_df.csv", stringsAsFactors = FALSE)
colnames(ppi) <- c("Gene_1", "Gene_2", "Interaction_type")


sig_gene_allias<-list(
c("ACS-13","Y65B4BL.5"),
c("ASPM-1","C45G3.1","tag-255"),
c("C16C2.4"),
c("CHD-7","tag-192","T04D1.4","scd-3"),
c("CSNK-1","Y106G6E.6"),
c("DHS-4","T05F1.10"),
c("DID-2","F23C8.6","phi-24"),
c("DRSH-1","F26E4.10"),
c("F26B1.2","hrpk-1"),
c("F31C3.3","fmp-27"),
c("F31C3.4","zipt-1"),
c("FHOD-1","C46H11.11","C46H11.a"),
c("FKB-8","Y18D10A.25"),
c("HLH-2","M05B5.5"),
c("HMP-2","mad-1","K05C4.6"),
c("HRPF-1","W02D3.11"),
c("IMP-1","C36B1.12"),
c("KAL-1","K03D10.1"),
c("KLP-15","M01E11.6"),
c("KLP-16","C41G7.2"),
c("LPD-3","Y47G6A.23"),
c("MAG-1","R09B3.5"),
c("MIG-1","cfz-1", "Y34D9B.1", "CNT1", "frizzled 1"),
c("PRP-4","C36B1.5"),
c("SEC-12","phi-15", "K02B12.3"),
c("SKR-2","F46A9.4"),
c("SUP-17","DY3.7"),
c("USP-5","T27A3.2"),
c("W04A4.5","ints-4"),
c("W09G3.6"),
c("WWP-1","CeWWP1", "Y65B4B_11.a", "Y65B4BR.4"),
c("Y34D9A.7","Y34D9A.j"),
c("Y44E3A.6","edc-4"),
c("ZHP-3","K02B12.8")
)


sig_gene_for_ppi<-na.omit(unlist(lapply(sig_gene_allias,function(x){union(ppi[match(tolower(x),tolower(ppi[,1])),1],ppi[match(tolower(x),tolower(ppi[,2])),2])})))

sig_gene_ppiIndex<-unique(unlist(lapply(sig_gene_for_ppi,function(x){union(which(ppi[,1]==x),which(ppi[,2]==x))})))


sig_gene_for_ppi_extend<-union(ppi[sig_gene_ppiIndex,1],ppi[sig_gene_ppiIndex,2])
sig_gene_ppiIndex_extend<-which(!is.na(match(ppi[,1],sig_gene_for_ppi_extend))&!is.na(match(ppi[,2],sig_gene_for_ppi_extend)))


ppi_related<-ppi[sig_gene_ppiIndex_extend,]

ppi_related<-ppi_related[which(ppi_related$Interaction_type=="Physical"),]


##
genes_mapping_extended<-c()
genes_mapping_extended<-data.frame(symbol=union(ppi_related[,1],ppi_related[,2]),
	entrez=apply(data.frame(mapIds(
	  org.Ce.eg.db,
	  keys = union(ppi_related[,1],ppi_related[,2]),
	  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
	  column = "ENTREZID"
	),mapIds(
	  org.Ce.eg.db,
	  keys = toupper(union(ppi_related[,1],ppi_related[,2])),
	  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
	  column = "ENTREZID"
	)),1,function(x){
		ifelse(is.na(x[1]),x[2],x[1])
	}))


genes_mapping_extended[which(is.na(genes_mapping_extended$entrez)),2]<-mapIds(
	  org.Ce.eg.db,
	  #keys = c("sec-16A.1","chd-7","M03A1.13","Y18D10A.32","lin-15B","F46A9.15"),
	  keys = c("sec-16A.1","lin-15B"),
	  keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
	  column = "ENTREZID"
	)


# ─────────────────────────────────────────────────────────
# Build graph & annotate nodes
# ─────────────────────────────────────────────────────────
ppi_filtered<-ppi_related
g <- graph_from_data_frame(ppi_filtered, directed = FALSE)

# Add interaction type as edge attribute
E(g)$type <- ppi_filtered$Interaction_type

# Map Wormbase ID back to Entrez for annotation
V(g)$entrez <- genes_mapping_extended[match(V(g)$name, genes_mapping_extended[,1]),2]

# ─────────────────────────────────────────────────────────
# Network topology metrics
# ─────────────────────────────────────────────────────────
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g)
V(g)$clustcoef <- transitivity(g, type = "local", isolates = "zero")

summary(V(g)$degree)
summary(V(g)$betweenness)

# Identify hubs (top 5%)
hub_threshold <- quantile(V(g)$degree, 0.95)
V(g)$hub <- V(g)$degree >= hub_threshold

# ─────────────────────────────────────────────────────────
# Community Detection: Leiden
# ─────────────────────────────────────────────────────────

set.seed(123)
lc <- cluster_leiden(g,objective_function = "modularity",resolution_parameter = 2,n_iterations = 10)
modules <- communities(lc)
mod_sizes <- sizes(lc)
print(mod_sizes)

# Assign module ID to nodes
V(g)$module <- lc$membership
V(g)$module_size<-mod_sizes[V(g)$module]
# ─────────────────────────────────────────────────────────
# Automate enrichment for each module
# ─────────────────────────────────────────────────────────
results_list <- list()

for (m in seq_along(mod_sizes)) {
  
  module_nodes <- V(g)$entrez[V(g)$module == m]
  module_nodes <- module_nodes[!is.na(module_nodes)]
  
  ego_focal <- enrichGO(
	gene = module_nodes,
	OrgDb = org.Ce.eg.db,
	ont = "ALL",
	pAdjustMethod = "BH",
	readable = TRUE
  )  

  #if (!is.null(ego_focal) && nrow(ego_focal@result) > 0) {
    results_list[[paste0("Module_", m)]] <- ego_focal

}


########################

##
V(g)$module_name <- sapply(results_list,function(x){
	return(x[1,"Description"])
})[lc$membership]

toHighlight_module = 11
module_name_display <- "cytoskeleton organization"

toHighlight_experimentalGene_in_highlightModule =
  intersect(modules[[toHighlight_module]], sig_gene_for_ppi)

toHighlight_otherGene_in_highlightModule =
  setdiff(modules[[toHighlight_module]], sig_gene_for_ppi)


# Default coloring
V(g)$node_group <- "Other_nodes"
# Highlight module nodes differently
V(g)$node_group[V(g)$name %in% toHighlight_otherGene_in_highlightModule] <- "Module_other_genes"
V(g)$node_group[V(g)$name %in% toHighlight_experimentalGene_in_highlightModule] <- "Module_exp_genes"

## node feature  
# Label only module nodes, blank others
V(g)$label <- ifelse(V(g)$name %in% toHighlight_experimentalGene_in_highlightModule,V(g)$name, "")
# Stroke widths: strong in focal module, small otherwise
V(g)$node_size <- scales::rescale(V(g)$degree, to = c(min(V(g)$degree),max(V(g)$degree)))
V(g)$node_stroke <- ifelse(V(g)$node_group == "Other_nodes",0.1,0.3)
V(g)$node_alpha <- ifelse(V(g)$node_group == "Other_nodes",0.3,0.5)

## edge features 
# Identify edges with both endpoints inside the focal module
edge_nodes <- ends(g, E(g))  # matrix: each row = (from,to)
edge_inside_module <- edge_nodes[,1] %in% modules[[toHighlight_module]] & edge_nodes[,2] %in% modules[[toHighlight_module]]

E(g)$edge_group <- ifelse(edge_inside_module, "Module_edges", "Other_edges")
E(g)$edge_width <- ifelse(edge_inside_module, 0.5, 0.2)
E(g)$edge_alpha <- ifelse(edge_inside_module, 0.5, 0.15)


# Module name for display

set.seed(16)

ggraph(g, layout = "fr") +
  geom_edge_link(aes(edge_colour = edge_group,edge_alpha = edge_alpha, edge_width = edge_width)) +
  scale_edge_color_manual(
    values = c(
      "Other_edges" = "grey80",
      "Module_edges" = "grey80"
    )
  )+
  scale_edge_width(range = c(0.2, 0.2), guide = "none") +
  scale_edge_alpha(range = c(0.5, 0.5), guide = "none") +  
  geom_node_point(aes(colour = node_group, fill = node_group, size = node_size, alpha=node_alpha),
	shape = 21, stroke=0.1) +
  scale_color_manual(
    values = c(
      "Other_nodes" = "grey50",
      "Module_other_genes" = "#00BFC4",
      "Module_exp_genes" = "#d73027"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Other_nodes" = "grey50",
      "Module_other_genes" = "#00BFC4",
      "Module_exp_genes" = "#d73027"
    )
  ) + 
  scale_alpha(range = c(0.2, 0.8))+
  scale_size(range = c(0.2, 6))+
  geom_node_text(aes(label = label, color = node_group),
                 fontface = "bold",
                 size = 2.5,
                 vjust = -0.8) +
  #theme_void() +
  annotate(
    "text",
    x = 3.5, y = -1,
    label = module_name_display,
    color = "#00BFC4",
    fontface = "bold",
    size = 3,
    vjust = 1.5
  )

ggsave("E:/canalizationENV/figs/velocity_sigGene_PPI_network.pdf",width=unit(8,"in"),height=unit(4,"in"))

####################



##########################


# Keep only nodes in the target set AND only edges among them
subset_gene<-modules[[11]]

sub_nodes <- setdiff(intersect(subset_gene, V(g)$name),c(""))
g_sub <- induced_subgraph(g, vids = sub_nodes)

# Label only module nodes, blank others
V(g_sub)$label <- V(g_sub)$name

set.seed(8)
ggraph(g_sub, layout = "fr") +
  # Weak edges globally; remove borders
  geom_edge_link(alpha = 0.35, width = 0.5, color = "#00BFC4") +
  # Node points: stroke and alpha fully controlled
  geom_node_point(aes(color = node_group, fill = node_group, size = node_size),
                  show.legend = TRUE, alpha = 0.35, shape = 21,stroke=0.5) +
  # Label only module nodes; no repel or segments
  geom_node_text(aes(label = label, color = node_group),
                 fontface = "bold",
                 size = 2.5,
                 vjust = -0.8) +
  # Consistent manual colors for nodes and text
  scale_color_manual(
    values = c(
      "Other_nodes" = "grey80",
      "Module_other_genes" = "#00BFC4",
      "Module_exp_genes" = "#d73027"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Other_nodes" = "grey80",
      "Module_other_genes" = "#00BFC4",
      "Module_exp_genes" = "#d73027"
    )
  ) +
  scale_size(range(0.3,1))+
  guides(color = "none", fill = "none") +
  # Add module name text near center of highlighted module
  annotate(
    "text",
    x = 1, 
    y = -1, 
    label = module_name_display,
    color = "#00BFC4",   # match non-exp gene highlight color
    fontface = "bold",
    size = 4,
    vjust = 1.5
  )

ggsave("E:/canalizationENV/figs/velocity_sigGene_PPI_network-subgraph.pdf",width=unit(4,"in"),height=unit(3,"in"))

############### individual embryo example 

sigCellIndex<-apply(velocity_var_by_gene_padj[match(c("IMP-1","KLP-15","KLP-16"),rownames(velocity_var_by_gene_padj)),],1,function(x){which(x<0.05)})

sigCellPadj<-apply(velocity_var_by_gene_padj[match(c("IMP-1","KLP-15","KLP-16"),rownames(velocity_var_by_gene_padj)),],1,function(x){x[which(x<0.05)]})

sigCell<-apply(velocity_var_by_gene_padj[match(c("IMP-1","KLP-15","KLP-16"),rownames(velocity_var_by_gene_padj)),],1,function(x){names(x)[which(x<0.05)]})

sigCellVariance<-lapply(1:3,function(i){
	temp<-velocity_var_by_gene[c(1,match(c("IMP-1","KLP-15","KLP-16")[i],rownames(velocity_var_by_gene))),sigCellIndex[[i]]]
	observed_fold<-temp[2,]/temp[1,]
	expected_fold<-apply(velocity_var_by_gene[-1,sigCellIndex[[i]]],2,mean,na.rm=T)/temp[1,]
	temp<-rbind(temp,observed_fold,expected_fold)
	rownames(temp)[3:4]<-c("FC_observed","FC_expected")
	return(temp)
})



## plot 

i=1
focalIndex<-match(do.call(rbind,strsplit(sigCell[[i]],"_"))[,1],cell714)
col1<-rep("grey40",714+length(focalIndex))
col1[c(focalIndex,c(1:length(focalIndex))+714)]<-c(rep("#F8766D",length(focalIndex)),rep("#00BFC4",length(focalIndex)))
fill1<-rep("grey40",714+length(focalIndex))
fill1[c(focalIndex,c(1:length(focalIndex))+714)]<-c(rep("#F8766D",length(focalIndex)),rep("#00BFC4",length(focalIndex)))
size1<-rep(0.6,714+length(focalIndex))
size1[c(focalIndex,c(1:length(focalIndex))+714)]<-size1[focalIndex]*sqrt(c(as.numeric(sigCellVariance[[i]][3,]),as.numeric(sigCellVariance[[i]][4,])))
stroke1<-rep(0.3,714)
stroke1[c(focalIndex,c(1:length(focalIndex))+714)]<-1
alpha1<-rep(0.5,714+length(focalIndex))
alpha1[c(focalIndex,c(1:length(focalIndex))+714)]<-1

dataTemp<-data.frame(x=c(apply(q1_cell350_alignment_xyz[,1:714],2,mean),apply(q1_cell350_alignment_xyz[,1:714],2,mean)[focalIndex]),
	y=c(apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean),apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)[focalIndex]))

ggplot(dataTemp,aes(x=x,y=y))+
geom_point(size=size1,stroke=0.1,col=col1,fill=fill1,alpha=alpha1)+
labs(x="Coordinate x",y="Coordinate y",title=c("IMP-1","KLP-15","KLP-16")[i])+
scale_y_continuous(breaks=c(-40,0,40))+
scale_x_continuous(breaks=c(-100,-50,0,50,100))


ggsave("E:/canalizationENV/figs/velocity_variance_foldChange_embryo_IMP-1.pdf",width=unit(2.1,"in"),height=unit(1.6,"in"))


i=2
focalIndex<-match(do.call(rbind,strsplit(sigCell[[i]],"_"))[,1],cell714)
col1<-rep("grey40",714+length(focalIndex))
col1[c(focalIndex,c(1:length(focalIndex))+714)]<-c(rep("#F8766D",length(focalIndex)),rep("#00BFC4",length(focalIndex)))
fill1<-rep("grey40",714+length(focalIndex))
fill1[c(focalIndex,c(1:length(focalIndex))+714)]<-c(rep("#F8766D",length(focalIndex)),rep("#00BFC4",length(focalIndex)))
size1<-rep(0.6,714+length(focalIndex))
size1[c(focalIndex,c(1:length(focalIndex))+714)]<-size1[focalIndex]*sqrt(c(as.numeric(sigCellVariance[[i]][3,]),as.numeric(sigCellVariance[[i]][4,])))
stroke1<-rep(0.3,714)
stroke1[c(focalIndex,c(1:length(focalIndex))+714)]<-1
alpha1<-rep(0.5,714+length(focalIndex))
alpha1[c(focalIndex,c(1:length(focalIndex))+714)]<-1

dataTemp<-data.frame(x=c(apply(q1_cell350_alignment_xyz[,1:714],2,mean),apply(q1_cell350_alignment_xyz[,1:714],2,mean)[focalIndex]),
	y=c(apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean),apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)[focalIndex]))

ggplot(dataTemp,aes(x=x,y=y))+
geom_point(size=size1,stroke=0.1,col=col1,fill=fill1,alpha=alpha1)+
labs(x="Coordinate x",y="Coordinate y",title=c("IMP-1","KLP-15","KLP-16")[i])+
scale_y_continuous(breaks=c(-40,0,40))+
scale_x_continuous(breaks=c(-100,-50,0,50,100))


ggsave("E:/canalizationENV/figs/velocity_variance_foldChange_embryo_KLP-15.pdf",width=unit(2.1,"in"),height=unit(1.6,"in"))


i=3
focalIndex<-match(do.call(rbind,strsplit(sigCell[[i]],"_"))[,1],cell714)
col1<-rep("grey40",714+length(focalIndex))
col1[c(focalIndex,c(1:length(focalIndex))+714)]<-c(rep("#F8766D",length(focalIndex)),rep("#00BFC4",length(focalIndex)))
fill1<-rep("grey40",714+length(focalIndex))
fill1[c(focalIndex,c(1:length(focalIndex))+714)]<-c(rep("#F8766D",length(focalIndex)),rep("#00BFC4",length(focalIndex)))
size1<-rep(0.6,714+length(focalIndex))
size1[c(focalIndex,c(1:length(focalIndex))+714)]<-size1[focalIndex]*sqrt(c(as.numeric(sigCellVariance[[i]][3,]),as.numeric(sigCellVariance[[i]][4,])))
stroke1<-rep(0.3,714)
stroke1[c(focalIndex,c(1:length(focalIndex))+714)]<-1
alpha1<-rep(0.5,714+length(focalIndex))
alpha1[c(focalIndex,c(1:length(focalIndex))+714)]<-1

dataTemp<-data.frame(x=c(apply(q1_cell350_alignment_xyz[,1:714],2,mean),apply(q1_cell350_alignment_xyz[,1:714],2,mean)[focalIndex]),
	y=c(apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean),apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)[focalIndex]))

ggplot(dataTemp,aes(x=x,y=y))+
geom_point(size=size1,stroke=0.1,col=col1,fill=fill1,alpha=alpha1)+
labs(x="Coordinate x",y="Coordinate y",title=c("IMP-1","KLP-15","KLP-16")[i])+
scale_y_continuous(breaks=c(-40,0,40))+
scale_x_continuous(breaks=c(-100,-50,0,50,100))


ggsave("E:/canalizationENV/figs/velocity_variance_foldChange_embryo_KLP-16.pdf",width=unit(2.1,"in"),height=unit(1.6,"in"))



# show legend 

dataTemp<-data.frame(y=rep(0,6),
	x=c(1:6))

ggplot(dataTemp,aes(x=x,y=y))+
geom_point(size=c(1:6)*0.6,stroke=0.1,col="black",fill=NA,alpha=1)+
scale_y_continuous(breaks=c(-40,0,40))+
scale_x_continuous(breaks=c(-100,-50,0,50,100))

ggsave("E:/canalizationENV/figs/velocity_variance_foldChange_embryo_legend.pdf",width=unit(2.1,"in"),height=unit(1.6,"in"))

###
pdf("E:/canalizationENV/figs/velocity_ppi_method.pdf",width=10,height=5)
set.seed(1)
par(mfrow=c(1,2))
plot(density(rnorm(100,0,1)),xlim=c(-10,10))
plot(density(rnorm(100,0,sqrt(sigCellVariance[[1]][3,1]))),xlim=c(-10,10))
dev.off()


###############

##
V(g)$module_name <- sapply(results_list,function(x){
	return(x[1,"Description"])
})[lc$membership]

## plot 
V(g)$size <- scales::rescale(sqrt(V(g)$degree), to = c(3, 8))

toHighlight_module = 11

toHighlight_experimentalGene_in_highlightModule =
  intersect(modules[[toHighlight_module]], sig_gene_for_ppi)

toHighlight_otherGene_in_highlightModule =
  setdiff(modules[[toHighlight_module]], sig_gene_for_ppi)

# Default coloring
V(g)$node_group <- "Other_nodes"

# Highlight module nodes differently
V(g)$node_group[V(g)$name %in% toHighlight_otherGene_in_highlightModule] <- "Module_other_genes"
V(g)$node_group[V(g)$name %in% toHighlight_experimentalGene_in_highlightModule] <- "Module_exp_genes"

# Label only module nodes, blank others
V(g)$label <- ifelse(V(g)$name %in% modules[[toHighlight_module]],
                     V(g)$name, "")

ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.2, stroke=0.1) +  # all weak edges
  geom_node_point(aes(color = node_group, size = size),
                  show.legend = TRUE,alpha=0.4, stroke=0.1) +
  geom_node_text(aes(label = label),
                 repel = TRUE,
                 size = 3) +
  scale_color_manual(
    values = c(
      "Other_nodes" = "grey80",       # weak background
      "Module_other_genes" = "#91bfdb",  # highlight non-exp
      "Module_exp_genes" = "#d73027"     # highlight experimental
    )
  ) +
  guides(size = "none") +
  theme_void() +
  ggtitle(paste0("Highlighted Module ", toHighlight_module))




ggraph(g, layout = "auto") +
  geom_edge_link(aes(color = type), alpha = 0.4) +
  geom_node_point(aes(size = size), color = "#91bfdb") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_discrete(name = "Interaction Type") +
  guides(size = "none") +
  theme_void() +
  ggtitle("Global PPI Network With Gene Names")




#wnt_ligands<-c("cwn-1", "cwn-2", "andmom-2")
#wnt_ligands_related<-unique(unlist(as.vector(ppi[na.omit(c(match(wnt_ligands,ppi[,1]),match(wnt_ligands,ppi[,2]))),1:2])))
#sort(unique(unlist(ppi[na.omit(c(match(wnt_ligands_related,ppi[,1]),match(wnt_ligands_related,ppi[,2]))),1:2])))
sig_gene_allias2<-list(
c("ACS-13","Y65B4BL.5"),
c("ASPM-1","C45G3.1","tag-255"),
c("ATM-1"),
c("C16C2.4"),
c("C27A12.7"),
c("C30H7.2"),
c("CHD-7","tag-192","T04D1.4","scd-3"),
c("CSNK-1","Y106G6E.6"),
c("DHS-4","T05F1.10"),
c("DID-2","F23C8.6","phi-24"),
c("DRSH-1","F26E4.10"),
c("DXBP-1"),
c("EARS-1"),
c("F26B1.2","hrpk-1"),
c("F31C3.3","fmp-27"),
c("F31C3.4","zipt-1"),
c("F43G9.12"),
c("F52A8.6"),
c("F57B10.5"),
c("FHOD-1","C46H11.11","C46H11.a"),
c("FKB-8","Y18D10A.25"),
c("HLH-2","M05B5.5"),
c("HMP-2","mad-1","K05C4.6"),
c("HRPF-1","W02D3.11"),
c("IMP-1","C36B1.12"),
c("INST-1"),
c("K12C11.1"),
c("KAL-1","K03D10.1"),
c("KLP-15","M01E11.6"),
c("KLP-16","C41G7.2"),
c("LPD-3","Y47G6A.23"),
c("MAG-1","R09B3.5"),
c("MIG-1","cfz-1", "Y34D9B.1", "CNT1", "frizzled 1"),
c("NCBP-2"),
c("NDX-7"),
c("PRP-4","C36B1.5"),
c("RPL-31"),
c("SEC-12","phi-15", "K02B12.3"),
c("SKR-1"),
c("SKR-2","F46A9.4"),
c("SNX-3"),
c("SUP-17","DY3.7"),
c("TOL-1"),
c("UBC-25"),
c("USP-5","T27A3.2"),
c("W04A4.5","ints-4"),
c("W09G3.6"),
c("WVE-1"),
c("WWP-1","CeWWP1", "Y65B4B_11.a", "Y65B4BR.4"),
c("Y34D9A.7","Y34D9A.j"),
c("Y44E3A.6","edc-4"),
c("Y47D9A.1"),
c("Y92H12A.2"),
c("ZHP-3","K02B12.8")  
)

#sig_gene_for_ppi<-na.omit(unlist(lapply(sig_gene_allias,function(x){union(ppi[match(tolower(x),tolower(ppi[,1])),1],ppi[match(tolower(x),tolower(ppi[,2])),2])})))

# focalIndex<-sort(sapply(unique(unlist(lapply(ego$geneID,function(x){unlist(strsplit(x,"/"))}))),function(y){
	# index<-c()
	# for(i in 1:length(sig_gene_allias)){
		# if(length(intersect(y,tolower(sig_gene_allias[[i]])))>0){
			# index<-c(index,i)
		# }
	# }
	# return(index)
# }))

#sig_gene_for_ppi<-na.omit(unlist(lapply(sig_gene_allias[focalIndex],function(x){union(ppi[match(tolower(x),tolower(ppi[,1])),1],ppi[match(tolower(x),tolower(ppi[,2])),2])})))
# direct_extended_genes<-union(ppi[,1],ppi[,2])[sapply(union(ppi[,1],ppi[,2]),function(x){
	# tempIndex<-union(which(ppi[,1]==x),which(ppi[,2]==x))
	# return(length(na.omit(match(union(ppi[tempIndex,1],ppi[tempIndex,2]),sig_gene_for_ppi)))>1)
# })]
# direct_extended_genes<-union(sig_gene_for_ppi,direct_extended_genes)

# sig_gene_ppiIndex<-which(!is.na(match(ppi[,1],direct_extended_genes))&!is.na(match(ppi[,2],direct_extended_genes)))

matrix(sig_gene_for_ppi[c(as.matrix(na.omit(data.frame(match(ppi_related[,1],sig_gene_for_ppi),match(ppi_related[,2],sig_gene_for_ppi)))))],nrow=3)
# genes_mapping_extended[which(is.na(genes_mapping_extended$entrez)),2]<-mapIds(
	  # org.Ce.eg.db,
	  # keys = c("sec-16A.1"),
	  # keytype = "SYMBOL",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
	  # column = "ENTREZID"
	# )
# sec-16a.1   sec-16a.1   <NA> # sec-16A.1
# scd-3           scd-3   <NA> # chd-7,tag-192
# m03a1.13     m03a1.13   <NA> # M03A1.13,cel_circ_008517
# y18d10a.32 y18d10a.32   <NA> # Y18D10A.32, cel_circ_008503 
# lin-15b       lin-15b   <NA> # lin-15, ZK662.4
# f46a9.15     f46a9.15   <NA> # F46A9.15, cel_circ_008445



# # Get genes for a specific GO ID
# library(org.Ce.eg.db)
# library(AnnotationDbi)

# go_id <- "GO:1902850"  # biological_process
# xx <- as.list(org.Ce.egGO2ALLEGS)
# if(length(xx) > 0){
    # # Gets the entrez gene ids for the top 2nd and 3nd GO identifiers
    # goids <- xx[2:3]
    # # Gets the entrez gene ids for the first element of goids
    # goids[[1]]
    # # Evidence code for the mappings
    # names(goids[[1]])
# }
# match(go_id,names(xx))

# print(xx[[go_id]])

# print(xx[["GO:1902850"]])
# focal_symbol<-to.lower(mapIds(
	  # org.Ce.eg.db,
	  # keys = xx[["GO:1902850"]],
	  # keytype = "ENTREZID",   # change to your input type: SYMBOL/WORMBASE/ENTREZID
	  # column = "SYMBOL"
	# ))
# focal_pairIndex<-which(apply(data.frame(match(ppi[,1],focal_symbol),match(ppi[,2],focal_symbol)),1,function(x){length(which(is.na(x)))})==0)
# ppi_related<-ppi[focal_pairIndex,]
#cat("Nodes:", gorder(), "Edges:", gsize(), "\n")

# ─────────────────────────────────────────────────────────
# Visualize top enriched terms across modules
# ─────────────────────────────────────────────────────────
plot_list <- list()
i <- 1
for (nm in names(results_list)) {
  plot_list[[i]] <- dotplot(results_list[[nm]], showCategory = 10) +
    ggtitle(paste0(nm, " GO BP"))
  i <- i + 1
}

# combine them
wrapped_plot <- wrap_plots(plot_list, ncol = 2)
wrapped_plot

set.seed(5)
cnetplot(results_list[[1]], showCategory = 20,layout = "gem") # 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'

genes_mapping_extended[match(V(g)$entrez[V(g)$module == 1],genes_mapping_extended[,2]),]

#
network_module_GO<-do.call(rbind,lapply(results_list,function(x){data.frame(x)}))
network_module_GO$overlapped_genes<-sapply(network_module_GO$geneID,function(x){paste0(sig_gene_for_ppi[na.omit(match(unlist(strsplit(x,"/")),sig_gene_for_ppi))],collapse="/")})
network_module_GO$overlapped_genes<-sapply(lapply(network_module_GO$overlapped_genes,function(x){
	if(x==""){
		return(NA)
	}else{
		focal_genes<-unlist(strsplit(x,"/"))
		focal_overlap<-c()
		for(gene in focal_genes){
			for(focal_gene_allias in sig_gene_allias){
				if(length(intersect(gene,tolower(focal_gene_allias)))==1){
					focal_overlap<-c(focal_overlap,tolower(focal_gene_allias[1]))
				}
			}
		}
		return(focal_overlap)	
	}

}),paste0,collapse="/")

network_module_GO$overlapped_geneCount<-sapply(network_module_GO$overlapped_genes,function(x){
	if(x=="NA"){
		return(0)
	}else{
		return(length(unlist(strsplit(x,"/"))))
	}
})
network_module_GO_overlapped<-network_module_GO[network_module_GO$overlapped_geneCount>1,]
network_module_GO_overlapped[order(network_module_GO_overlapped$overlapped_genes),]
network_module_GO_overlapped[order(network_module_GO_overlapped$qvalue),]

network_module_GO_overlapped$GeneRatio<-paste("#",network_module_GO_overlapped$GeneRatio)
network_module_GO_overlapped$BgRatio<-paste("#",network_module_GO_overlapped$BgRatio)

write.csv(network_module_GO_overlapped,"E:/canalizationENV/data/velocity/network_module_GO_overlapped.csv")

##


##

ego_global <- enrichGO(
	gene = V(g)$entrez,
	OrgDb = org.Ce.eg.db,
	ont = "ALL",
	pAdjustMethod = "BH",
	readable = TRUE
)

head(ego_global)

network_global_GO<-data.frame(ego_global)
network_global_GO$overlapped_genes_1<-sapply(network_global_GO$geneID,function(x){paste0(sig_gene_for_ppi[na.omit(match(unlist(strsplit(x,"/")),sig_gene_for_ppi))],collapse="/")})

network_global_GO$overlapped_genes_2<-network_global_GO$overlapped_genes_1
network_global_GO$overlapped_genes_2<-sapply(network_global_GO$overlapped_genes_1,function(x){
	if(x==""){
		return("")
	}else{
		focal_genes<-unlist(strsplit(x,"/"))
		focal_ppi_index<-na.omit(data.frame(match(ppi[,1],focal_genes),match(ppi[,2],focal_genes)))
		return(paste0(focal_genes[union(focal_ppi_index[,1],focal_ppi_index[,2])],collapse="/"))
	}
})

network_global_GO$overlapped_genes_3<-network_global_GO$overlapped_genes_2
network_global_GO$overlapped_genes_3<-sapply(network_global_GO$overlapped_genes_2,function(x){
	if(x==""){
		return(NA)
	}else{
		focal_genes<-unlist(strsplit(x,"/"))
		focal_overlap<-c()
		for(gene in focal_genes){
			for(focal_gene_allias in sig_gene_allias){
				if(length(intersect(gene,tolower(focal_gene_allias)))==1){
					focal_overlap<-c(focal_overlap,tolower(focal_gene_allias[1]))
				}
			}
		}
		return(paste0(focal_overlap,collapse="/"))	
	}

})

network_global_GO$overlapped_geneCount_1<-sapply(network_global_GO$overlapped_genes_1,function(x){
	if(is.na(x)){
		return(0)
	}else{
		return(length(unlist(strsplit(x,"/"))))
	}
})

network_global_GO$overlapped_geneCount_2<-sapply(network_global_GO$overlapped_genes_2,function(x){
	if(is.na(x)){
		return(0)
	}else{
		return(length(unlist(strsplit(x,"/"))))
	}
})

network_global_GO$overlapped_geneCount_3<-sapply(network_global_GO$overlapped_genes_3,function(x){
	if(is.na(x)){
		return(0)
	}else{
		return(length(unlist(strsplit(x,"/"))))
	}
})

network_global_GO_overlapped<-network_global_GO[network_global_GO$overlapped_geneCount_1>5,]
network_global_GO_overlapped[order(network_global_GO_overlapped$overlapped_genes_1),]

network_global_GO_overlapped<-network_global_GO[network_global_GO$overlapped_geneCount_2>1,]
network_global_GO_overlapped[order(network_global_GO_overlapped$overlapped_genes_2),]

network_global_GO_overlapped[order(network_global_GO_overlapped$qvalue),]


# Keep only nodes in the target set AND only edges among them
subset_gene<-setdiff(unique(unlist(strsplit(network_global_GO_overlapped[order(network_global_GO_overlapped$overlapped_genes_1)[2],"geneID"],"/"))),c("bir-1","rga-8","sax-7"))

sub_nodes <- setdiff(intersect(subset_gene, V(g)$name),c(""))
g_sub <- induced_subgraph(g, vids = sub_nodes)

##
V(g_sub)$size <- scales::rescale(V(g_sub)$degree, to = c(3, 8))

ggraph(g_sub, layout = "auto") +
  geom_edge_link(aes(color = type), alpha = 0.4) +
  geom_node_point(aes(size = size), color = "#91bfdb") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_discrete(name = "Interaction Type") +
  guides(size = "none") +
  theme_void() +
  ggtitle("Global PPI Network With Gene Names")


# Your selected target genes
target_gene_set <- c("hmp-2","csnk-1","skr-2","aspm-1","klp-15","fhod-1","klp-16")

# Create node color attribute: target genes vs. other genes
V(g_sub)$group <- ifelse(V(g_sub)$name %in% target_gene_set, "Experimented", "Other")

# Scale node size from degree
V(g_sub)$size <- scales::rescale(degree(g_sub), to = c(3, 8))

# Plot
ggraph(g_sub, layout = "fr") +
  geom_edge_link(aes(color = type), alpha = 0.4) +
  geom_node_point(aes(size = size, fill = group), shape = 21) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_fill_manual(values = c("Experimented" = "#d73027",   # red
                               "Other"  = "#91bfdb"))+  # blue
  scale_edge_color_discrete(name = "Interaction Type") +
  guides(size = "none") +
  theme_void() +
  ggtitle("cytoskeleton organization") +
  theme(plot.title = element_text(hjust = 0.5))


# Ensure type is a factor for consistent ordering
E(g_sub)$type <- factor(E(g_sub)$type, levels = c("Physical", "Regulatory", "Genetic"))

# Plot with clearer colors + thicker edges
ggraph(g_sub, layout = "fr") +
  geom_edge_link(aes(color = type), alpha = 0.5) +
  geom_node_point(aes(size = size, fill = group), shape = 21, color = "black", stroke = 0.6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_manual(
    name = "Interaction Type",
    values = c(
      "Physical"   = "blue",  # strong blue
      "Regulatory" = "red",  # vivid purple
      "Genetic"    = "#00BFC4"   # publication red
    )
  ) +
  scale_fill_manual(
    name = "Node Group",
    values = c(
      "Experimented" = "#ff7f00",   # orange pops out
      "Other"  = "#a6cee3"
    )
  ) +
  guides(size = "none") +
  theme_void() +
  ggtitle("cytoskeleton organization") +
  theme(plot.title = element_text(hjust = 0.5))



set.seed(123)
# Filter edges by type
phys_edges <- E(g_sub)[E(g_sub)$type == "Physical"]

# Keep nodes connected by these edges
nodes_phys <- unique(as.vector(ends(g_sub, phys_edges)))

# Subgraph containing those nodes + only physical edges
g_sub2 <- induced_subgraph(g = g_sub, vids = nodes_phys)
g_sub2 <- delete_edges(g_sub2, E(g_sub2)[E(g_sub2)$type != "Physical"])

ggraph(g_sub2, layout = "fr") +
  geom_edge_link(aes(color = type), alpha = 0.5) +
  geom_node_point(aes(size = size, fill = group), shape = 21, color = "black", stroke = 0.6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_manual(
    name = "Interaction Type",
    values = c(
      "Physical"   = "blue",  # strong blue
      "Regulatory" = "red",  # vivid purple
      "Genetic"    = "#00BFC4"   # publication red
    )
  ) +
  scale_fill_manual(
    name = "Node Group",
    values = c(
      "Experimented" = "#ff7f00",   # orange pops out
      "Other"  = "#a6cee3"
    )
  ) +
  guides(size = "none") +
  theme_void() +
  ggtitle("cytoskeleton organization") +
  theme(plot.title = element_text(hjust = 0.5))



lc_sub <- cluster_leiden(g_sub2,objective_function = "modularity",resolution_parameter = 3,n_iterations = 10)
modules <- communities(lc_sub)
modules





network_global_GO_overlapped$GeneRatio<-paste("#",network_global_GO_overlapped$GeneRatio)
network_global_GO_overlapped$BgRatio<-paste("#",network_global_GO_overlapped$BgRatio)

write.csv(network_global_GO_overlapped,"E:/canalizationENV/data/velocity/network_global_GO_overlapped.csv")





# Identify module 77 nodes
# module_id <- 77
# module77_nodes <- V(g)$name[V(g)$module == module_id]

# # Your five enriched genes among module 77
# # (from your printed mapping result)
# sig5_genes <- unique(unlist(lapply(results_list[[1]]$geneID,function(x){unlist(strsplit(x,"/"))})))


# # Create node attributes for plotting
# V(g)$highlight <- ifelse(V(g)$name %in% sig5_genes, "Enriched5",
                         # ifelse(V(g)$name %in% module77_nodes, "Module77", "Other"))

# # Color scale:
# color_map <- c(
  # "Other" = "grey80",
  # "Module77" = "#91bfdb",
  # "Enriched5" = "#d73027"
# )

# # Node size by degree
# V(g)$size <- scales::rescale(V(g)$degree, to = c(3, 8))

# # Plot network
# set.seed(77)
# ggraph(g, layout = "fr") +
  # geom_edge_link(alpha = 0.25) +
  # geom_node_point(aes(color = highlight, size = size)) +
  # scale_color_manual(values = color_map) +
  # guides(size = "none") +
  # theme_void() +
  # ggtitle("Global PPI network highlighting Module 77 and its enriched genes")



ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.2) +
  geom_node_point(aes(color = highlight, size = size)) +
  geom_node_text(data = . %>% filter(highlight == "Enriched5"),
                 aes(label = name),
                 size = 3, repel = TRUE) +
  scale_color_manual(values = color_map) +
  theme_void() +
  ggtitle("Highlighted enriched genes in Module 77")

library(concaveman)

ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.2) +
  geom_mark_hull(aes(x, y, filter = highlight == "Module77"),
                 concavity = 5, expand = unit(3, "mm"),
                 fill = "#91bfdb40", color = NA) +
  geom_node_point(aes(color = highlight, size = size)) +
  scale_color_manual(values = color_map) +
  theme_void()



###############################################################################
# Pathway visualization: Pathview for KEGG maps
###############################################################################

# BiocManager::install("pathview")
library(pathview)

# Suppose your gene list is Entrez IDs; create a named vector if you have fold-changes, else use 1
gene_data <- rep(1, length(gene_mapped_selected))
names(gene_data) <- gene_mapped_selected

# Example pathway id from KEGG results (choose a significant pathway)
kegg_id <- "cel04141" # replace with an actual KEGG pathway id from ekegg
pv <- pathview(gene.data = gene_data, pathway.id = kegg_id, species = "cel", out.suffix = "pathview")



###############################################################################
# human orthologs : Map C. elegans → human orthologs (biomaRt) and run Reactome 
###############################################################################

# Install if needed
# BiocManager::install(c("biomaRt","ReactomePA","clusterProfiler","org.Hs.eg.db"))

library(biomaRt)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

library(dplyr)
library(readr)

# update path to where you saved the file
orth_file <- "c_elegans.orthology.txt.gz"

orth_map <- read_tsv(orth_file, comment = "#") %>%
  filter(Species == "Homo sapiens") %>%
  select(
    wormbase_gene = GeneID,
    human_symbol = Ortholog,
    human_source = Source
  ) %>%
  distinct()

head(orth_map)


# orths will contain CEENTREZ / HS symbol mapping; extract human Entrez IDs:
human_entrez <- unique(na.omit(as.character(orths[, "EntrezGene"])) ) # column name may vary

# Reactome enrichment on mapped human genes
rreact <- enrichPathway(gene = human_entrez, organism = "human", pAdjustMethod = "BH")
head(rreact)
dotplot(rreact, showCategory = 20) + ggtitle("Reactome (via human orthologs)")



###############################################################################
# End of file
###############################################################################


focal_TF_allias<-list(
c("atf-7","C07G2.2"), 
c("ztf-6","dopy-1","W06H12.1"), 
c("sea-1","Ce-tbx-18", "tbx-18", "F19B10.9"), 
c("dpl-1","mex-4", "T23G7.1", "szy-10"), 
c("daf-8","R05D11.1", "smad-1"), 
c("crh-2","C27D6.4"), 
c("daf-19","F33H1.1", "daf-24"), 
c("syd-9","ztf-10", "ZK867.1", "tag-239"), 
c("K09A11.1"), 
c("svh-5","C33A11.4")
)

ppi_TF_index<-unique(unlist(lapply(focal_TF_allias[c(1,6,9)],function(x){
	return(unique(unlist(lapply(tolower(x),function(y){
		which(ppi[,1]==y|ppi[,2]==y)
	}))))
})))

ppi_TF<-ppi[ppi_TF_index,]

ppi_TF_index_extend<-unique(unlist(lapply(union(ppi_TF[,1],ppi_TF[,2]),function(x){union(which(ppi[,1]==x),which(ppi[,2]==x))})))

module11_gene_ppiIndex<-unique(unlist(lapply(modules[[11]],function(x){union(which(ppi[,1]==x),which(ppi[,2]==x))})))


ppi_TF_extend<-ppi[intersect(ppi_TF_index_extend,module11_gene_ppiIndex),] #2
ppi_TF_extend_index<-unique(unlist(lapply(union(ppi_TF_extend[,1],ppi_TF_extend[,2]),function(x){union(which(ppi[,1]==x),which(ppi[,2]==x))})))

ppi_TF_update<-ppi[intersect(ppi_TF_index,ppi_TF_extend_index),] #1


ppi_sig_gene_extend<-ppi[sig_gene_ppiIndex_extend,]
ppi_sig_gene_extend_update<-ppi_sig_gene_extend[which(!is.na(match(ppi_sig_gene_extend[,1],modules[[11]]))&!is.na(match(ppi_sig_gene_extend[,2],modules[[11]]))),]

ppi_TF_update #1
ppi_TF_extend #2
ppi_sig_gene_extend_update #3

rbind(ppi_TF_update,ppi_TF_extend,ppi_sig_gene_extend_update)
rbind(ppi_TF_update,ppi_TF_extend)
rbind(ppi_sig_gene_extend_update)



library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)

plot_ppi_network <- function(ppi_df,
                             seed = 1,
                             layout = "fr",
                             node_label_size = 3,
                             edge_alpha = 0.6) {
  
  # ---- 1. Build graph ----
  g <- graph_from_data_frame(ppi_df, directed = FALSE)

  # Convert to tbl_graph for tidygraph use if needed
  tg <- as_tbl_graph(g)

  # ---- 2. Node attributes ----
  tg <- tg %>%
    mutate(
      degree = centrality_degree(),
      size = scales::rescale(sqrt(degree), to = c(3, 10)),
      label = name
    )

  # ---- 3. Edge attributes ----
  # Edge color by interaction type
  interaction_colors <- c(
    "Physical"    = "#1f78b4",
    "Genetic"     = "#e31a1c",
    "Regulatory"  = "#33a02c"
  )

  # ---- 4. Plot ----
  set.seed(seed)

  ggraph(tg, layout = layout) +
    geom_edge_link(aes(color = Interaction_type),
                   width = 0.4,
                   alpha = edge_alpha,
                   show.legend = TRUE) +
    
    geom_node_point(aes(size = size),
                    shape = 21,
                    stroke = 0,
                    fill = "grey20",
                    color = "grey20",
                    alpha = 0.4) +
    
    geom_node_text(aes(label = label),
                   size = node_label_size,
                   repel = TRUE,
                   vjust = -0.8,
                   color = "blue") +
    
    scale_edge_color_manual(values = interaction_colors) +
    scale_size(range = c(3, 10), guide = "none") +
    theme_void() +
    theme(legend.position = "right")
}

ppi_df <- rbind(ppi_TF_update, ppi_TF_extend)
ppi_df_extend_index<-which(!is.na(match(ppi[,1],union(ppi_df[,1],ppi_df[,2])))&!is.na(match(ppi[,2],union(ppi_df[,1],ppi_df[,2]))))
ppi_df_extend<-ppi[ppi_df_extend_index,]


TF_gene<-c("atf-7","crh-2")
module11_gene<-c("klp-15","par-3","act-4","npp-19","mbk-1","plk-1","ced-3","pac-1","ebp-1","act-5","klp-16","lin-15","gdi-1","xpa-1","imp-1","acdh-2")
bridge_gene<-c("smo-1","lit-1","eor-1","unc-62","daf-16","pmk-1","hda-1","eya-1")
#bridge_gene<-setdiff(union(ppi_df_extend[,1],ppi_df_extend[,2]),c(TF_gene,modules[[11]]))
consider_gene<-c(TF_gene,bridge_gene,module11_gene)
gene<- readRDS("F:/downloads/15091632/cel_cbr_embryo_single_cell/exp_vis copy/data/expr_data.rds")
TF_gene_index<-grep(TF_gene[2],rownames(gene))
module11_gene_index<-sapply(module11_gene,function(x){
	focal_index<-grep(x,rownames(gene))
	ifelse(length(focal_index)>0,focal_index[1],NA)
})
bridge_gene_index<-sapply(bridge_gene,function(x){
	focal_index<-grep(x,rownames(gene))
	ifelse(length(focal_index)>0,focal_index[1],NA)
})



Celegans_idnex<-grep("C. elegans",rownames(gene))
par(mfrow=c(1,2))
hist(abs(cor(t(gene[TF_gene_index,]),t(gene[Celegans_idnex,]))[1,]))
hist(abs(cor(t(gene[TF_gene_index,]),t(gene[module11_gene_index,]))[1,]))

dataTemp<-data.frame(R_abs=c(abs(cor(t(gene[TF_gene_index,]),t(gene[Celegans_idnex,]))[1,]),abs(cor(t(gene[TF_gene_index,]),t(gene[module11_gene_index,]))[1,])),
	type=c(rep("all",length(cor(t(gene[TF_gene_index,]),t(gene[Celegans_idnex,]))[1,])),
		rep("module11",length(cor(t(gene[TF_gene_index,]),t(gene[module11_gene_index,]))[1,]))))

library(ggplot2)
ggplot(dataTemp,aes(x=R_abs,col=type,fill=type))+
geom_density(alpha=0.5)+
labs(title="correlation between TF(atf-7) and other genes")



dataTemp2<-data.frame(R_abs=c(abs(cor(t(gene[TF_gene_index,]),t(gene[Celegans_idnex,]))[1,]),abs(cor(t(gene[TF_gene_index,]),t(gene[bridge_gene_index,]))[1,])),
	type=c(rep("all",length(cor(t(gene[TF_gene_index,]),t(gene[Celegans_idnex,]))[1,])),
		rep("bridge",length(cor(t(gene[TF_gene_index,]),t(gene[bridge_gene_index,]))[1,]))))

library(ggplot2)
ggplot(dataTemp2,aes(x=R_abs,col=type,fill=type))+
geom_density(alpha=0.5)+
labs(title="correlation between TF(crh-2) and bridge genes")

wilcox.test(abs(cor(t(gene[TF_gene_index,]),t(gene[Celegans_idnex,]))[1,]),abs(cor(t(gene[TF_gene_index,]),t(gene[bridge_gene_index,]))[1,]))




library(pheatmap)
pheatmap(abs(cor(t(gene[unique(unlist(lapply(consider_gene,function(x){intersect(grep(x,rownames(gene)),grep("C. elegans",rownames(gene)))}))),]))),cluster_rows = FALSE,cluster_cols = FALSE)


module_gene<-intersect(union(ppi_df_extend[,1],ppi_df_extend[,2]),modules[[11]])

plot_ppi_network(ppi_df_extend,seed=2,TF_gene=TF_gene,module_gene=module_gene)


###

library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)

plot_ppi_network <- function(ppi_df,
                             seed = 1,
                             layout = "fr",
                             node_label_size = 5,
                             edge_alpha = 0.6,
                             TF_gene = NULL,
                             module_gene = NULL) {
  
  # ---- 1. Build graph ----
  g <- graph_from_data_frame(ppi_df, directed = FALSE)
  tg <- as_tbl_graph(g)

  # ---- 2. Node attributes ----
  tg <- tg %>%
    mutate(
      degree = centrality_degree(),
      size = scales::rescale(sqrt(degree), to = c(3, 10)),
      label = name,
      
      # Assign groups for node colors
      node_group = case_when(
        name %in% TF_gene ~ "TF_gene",
        name %in% module_gene ~ "module_gene",
        TRUE ~ "other"
      )
    )
  
  # ---- 3. Edge attributes ----
  interaction_colors <- c(
    "Physical"    = "#00BFC4",
    "Genetic"     = "#e31a1c",
    "Regulatory"  = "orange"
  )

  node_colors <- c(
    "TF_gene"     = "#d73027",  # red
    "module_gene" = "#00BFC4",  # green
    "other"       = "grey70"
  )

  # ---- 4. Plot ----
  set.seed(seed)

  ggraph(tg, layout = layout) +
    
    # Edges
    geom_edge_link(aes(color = Interaction_type),
                   width = 0.4,
                   alpha = edge_alpha,
                   show.legend = TRUE) +
    
    # Nodes
    geom_node_point(aes(size = size, fill = node_group),
                    shape = 21,
                    stroke = 0,
                    alpha = 0.4,
                    show.legend = TRUE) +
    
    # Node labels (colored same as nodes)
    geom_node_text(aes(label = label, color = node_group),
                   size = node_label_size,
                   repel = TRUE,
                   vjust = -0.8,
                   fontface = "bold") +
    
    scale_edge_color_manual(values = interaction_colors) +
    scale_fill_manual(values = node_colors, name = "Node class") +
    scale_color_manual(values = node_colors, guide = "none") +
    scale_size(range = c(3, 10), guide = "none") +
    theme_void() +
    theme(legend.position = "right")
}

ppi_df <- rbind(ppi_TF_update, ppi_TF_extend)

ppi_df_extend_index <- which(
  !is.na(match(ppi[,1], union(ppi_df[,1], ppi_df[,2]))) &
  !is.na(match(ppi[,2], union(ppi_df[,1], ppi_df[,2])))
)

ppi_df_extend <- ppi[ppi_df_extend_index,]

TF_gene <- c("atf-7", "crh-2")
module_gene <- intersect(
  union(ppi_df_extend[,1], ppi_df_extend[,2]),
  modules[[11]]
)

plot_ppi_network(ppi_df_extend,
                 seed = 2,
                 TF_gene = TF_gene,
                 module_gene = module_gene)









############################################################################################################







library(pheatmap)

custom_palette <- colorRampPalette(c("blue", "white", "red"))(50)
max_abs <- max(abs(velocity_average_by_gene_relative[-1,]))  # Maximum absolute value
breaks <- seq(-max_abs, max_abs, length.out = 51)
# Heatmap
pheatmap(
  velocity_average_by_gene_relative[-1,],
  color = custom_palette,
  breaks = breaks,
  main = "velocity average of gene or env - that of wild types)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = F, show_colnames = F
)
grid.text("gene/env", x = 0.02, y = 0.5, rot = 90, gp = gpar(fontsize = 12))  # Row axis title (left side)
grid.text("cell features", x = 0.5, y = 0.02, gp = gpar(fontsize = 12))          # Column axis title (bottom)

pheatmap(
  velocity_average_by_gene_relative_zero[-1,],
  color = custom_palette,
  breaks = breaks,
  main = "velocity average of gene or env - that of wild types)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = F, show_colnames = F
)
grid.text("gene/env", x = 0.02, y = 0.5, rot = 90, gp = gpar(fontsize = 12))  # Row axis title (left side)
grid.text("cell features", x = 0.5, y = 0.02, gp = gpar(fontsize = 12))          # Column axis title (bottom)


pheatmap(velocity_average_by_gene_relative_zero[-1,]^(1/3))

par(mfrow=c(2,1))
hist(apply(velocity_average_by_gene_relative_zero[-1,],1,function(x){length(which(x!=0))}),xlab="Number of cell features per gene/env",main="significant velocity change",100)
hist(apply(velocity_average_by_gene_relative_zero[-1,],2,function(x){length(which(x!=0))}),xlab="Number of gene/env per cell feature",main="significant velocity change",100)

hist(velocity_average_by_gene_relative[which(velocity_average_by_gene_padj<0.05)])

plot(velocity_average_by_gene_relative[,1:353],velocity_average_by_gene_relative[,c(1:353)+353])

## test variance ratio reative to initial_pos

velocity_var_by_gene<-do.call(rbind,lapply(unique(gene_extend),function(x){
	print(x)
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(rep(NA,))
	}else{
		return(apply(velocity[geneInd,],2,var))	
	}
}))
rownames(velocity_var_by_gene)<-unique(gene_extend)
colnames(velocity_var_by_gene)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]

velocity_var_by_gene_relative<-do.call(rbind,lapply(unique(gene_extend),function(x){
	print(x)
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(rep(NA,dim(velocity)[2]))
	}else{
		return(apply(velocity[geneInd,],2,var)/apply(velocity[1:105,],2,var))
	}	
}))
rownames(velocity_var_by_gene_relative)<-unique(gene_extend)
colnames(velocity_var_by_gene_relative)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]


hist(sapply(1:103,function(i){var.test(velocity[sample(1:105,3),1],velocity[1:105,1])$p.value}),xlab=NA)
hist(sapply(1:103,function(i){sum((velocity[c(i,(i+2)),1]-mean(velocity[c(i,(i+2)),1]))^2)/3}),xlab=NA)


velocity_var_by_gene_pvalue<-do.call(rbind,lapply(unique(gene_extend),function(x){
	print(x)
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(rep(NA,dim(velocity)[2]))
	}else{
		return(sapply(1:dim(velocity)[2],function(i){var.test(velocity[geneInd,i],velocity[1:105,i])$p.value}))
	}	
}))
rownames(velocity_var_by_gene_pvalue)<-unique(gene_extend)
colnames(velocity_var_by_gene_pvalue)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]


velocity_var_by_gene_refine<-velocity_var_by_gene[apply(velocity_var_by_gene_pvalue,1,function(x){length(which(is.na(x)))})==0,]
velocity_var_by_gene_relative_refine<-velocity_var_by_gene_relative[apply(velocity_var_by_gene_pvalue,1,function(x){length(which(is.na(x)))})==0,]
velocity_var_by_gene_pvalue_refine<-velocity_var_by_gene_pvalue[apply(velocity_var_by_gene_pvalue,1,function(x){length(which(is.na(x)))})==0,]

velocity_var_by_gene_padj<-matrix(p.adjust(velocity_var_by_gene_pvalue_refine,method="BH"),nrow=dim(velocity_var_by_gene_pvalue_refine)[1])
velocity_var_by_gene_padj[which(velocity_var_by_gene_padj==0)]<-min(setdiff(velocity_var_by_gene_padj,0))
rownames(velocity_var_by_gene_padj)<-rownames(velocity_var_by_gene_refine)
colnames(velocity_var_by_gene_padj)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]

write.csv(data.frame(gene=names(sort(apply(velocity_var_by_gene_padj,1,function(x){length(which(x<0.05))}))),
					 value=sort(apply(velocity_var_by_gene_padj,1,function(x){length(which(x<0.05))}))),
		  "E:/canalizationENV/ppt/velocity_var_by_gene_padj.csv")

## test the slope between velocity and initial_pos
k_wt<-sapply(1:dim(velocity)[2],function(i){coefficients(summary(lm(end_pos[1:105,i]~initial_pos[1:105,i])))[2,1:2]})
#pt(k_diff/k_diff_se,length(c(refInd,perturbInd))-2)
k_diff_matrix<-do.call(rbind,lapply(unique(gene_extend),function(x){
	print(x)
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(rep(NA,dim(velocity)[2]))
	}else{
		k_gene<-sapply(1:dim(velocity)[2],function(i){coefficients(summary(lm(end_pos[geneInd,i]~initial_pos[geneInd,i])))[2,1:2]})
		k_diff<-k_gene[1,]-k_wt[1,]
		return(k_diff)	
	}
}))
rownames(k_diff_matrix)<-unique(gene_extend)
colnames(k_diff_matrix)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]

k_diff_matrix_refine<-k_diff_matrix[apply(velocity_var_by_gene_pvalue,1,function(x){length(which(is.na(x)))})==0,]










velocity_slope_by_gene_pvalue<-do.call(rbind,lapply(unique(gene_extend),function(x){
	print(x)
	geneInd<-which(gene_extend==x)
	if(length(geneInd)==1){
		return(rep(NA,dim(velocity)[2]))
	}else{
		k_gene<-sapply(1:dim(velocity)[2],function(i){coefficients(summary(lm(velocity[geneInd,i]~initial_pos[geneInd,i])))[2,1:2]})
		k_diff<-k_gene[1,]-k_wt[1,]
		k_diff_se<-sqrt(k_gene[2,]^2+k_wt[2,]^2)
		k_diff_z<-k_diff/k_diff_se 
		k_diff_p<-1-pnorm(k_diff_z)
		#k_diff_p<-1-pt(k_diff/k_diff_se,105+length(geneInd)-2)
		return(k_diff_p)	
	}
}))
rownames(velocity_slope_by_gene_pvalue)<-unique(gene_extend)
colnames(velocity_slope_by_gene_pvalue)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]

velocity_slope_by_gene_pvalue_refine<-velocity_slope_by_gene_pvalue[apply(velocity_slope_by_gene_pvalue,1,function(x){length(which(is.na(x)))})==0,]

velocity_slope_by_gene_padj<-matrix(p.adjust(velocity_slope_by_gene_pvalue_refine,method="BH"),nrow=dim(velocity_slope_by_gene_pvalue_refine)[1])
rownames(velocity_slope_by_gene_padj)<-rownames(velocity_slope_by_gene_pvalue_refine)
colnames(velocity_slope_by_gene_padj)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]


write.csv(data.frame(gene=names(apply(velocity_slope_by_gene_padj,1,function(x){length(which(x<1e-15))})),
					 value=apply(velocity_slope_by_gene_padj,1,function(x){length(which(x<1e-15))})),
		  "E:/canalizationENV/ppt/velocity_slope_by_gene_padj.csv")



# velocity_var_by_gene_pvalue_smaller<-do.call(rbind,lapply(unique(gene_extend),function(x){
	# print(x)
	# geneInd<-which(gene_extend==x)
	# if(length(geneInd)==1){
		# return(rep(NA,dim(velocity)[2]))
	# }else{
		# return(sapply(1:dim(velocity)[2],function(i){var.test(velocity[geneInd,i],velocity[1:105,i],alternative="less")$p.value}))
	# }	
# }))
# rownames(velocity_var_by_gene_pvalue_smaller)<-unique(gene_extend)
# colnames(velocity_var_by_gene_pvalue_smaller)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]


# velocity_var_by_gene_refine_smaller<-velocity_var_by_gene[apply(velocity_var_by_gene_pvalue_smaller,1,function(x){length(which(is.na(x)))})==0,]
# velocity_var_by_gene_relative_refine_smaller<-velocity_var_by_gene_relative[apply(velocity_var_by_gene_pvalue_smaller,1,function(x){length(which(is.na(x)))})==0,]
# velocity_var_by_gene_pvalue_refine_smaller<-velocity_var_by_gene_pvalue_smaller[apply(velocity_var_by_gene_pvalue_smaller,1,function(x){length(which(is.na(x)))})==0,]

# velocity_var_by_gene_padj_smaller<-matrix(p.adjust(velocity_var_by_gene_pvalue_refine_smaller,method="BH"),nrow=dim(velocity_var_by_gene_pvalue_refine_smaller)[1])
# velocity_var_by_gene_padj_smaller[which(velocity_var_by_gene_padj_smaller==0)]<-min(setdiff(velocity_var_by_gene_padj_smaller,0))
# rownames(velocity_var_by_gene_padj_smaller)<-rownames(velocity_var_by_gene_refine_smaller)
# colnames(velocity_var_by_gene_padj_smaller)<-colnames(q1_cell350_alignment_xyz)[c(cell353_ind,cell353_ind+714,cell353_ind+714*2)]


# dataTemp<-data.frame(gene=rownames(velocity_var_by_gene_padj_smaller),
					 # value=sapply(1:dim(velocity_var_by_gene_padj_smaller)[1],function(i){length(which(velocity_var_by_gene_padj_smaller[i,]<0.05&velocity_var_by_gene_relative_refine_smaller[i,]<1))}))
dataTemp<-data.frame(gene=rownames(velocity_var_by_gene_padj),
					 value_smaller=sapply(1:dim(velocity_var_by_gene_padj)[1],function(i){length(which(velocity_var_by_gene_padj[i,]<0.05&velocity_var_by_gene_relative_refine[i,]<1))}),
					 value_larger=sapply(1:dim(velocity_var_by_gene_padj)[1],function(i){length(which(velocity_var_by_gene_padj[i,]<0.05&velocity_var_by_gene_relative_refine[i,]>1))}),
					 value_total=sapply(1:dim(velocity_var_by_gene_padj)[1],function(i){length(which(velocity_var_by_gene_padj[i,]<0.05))}))

write.csv(dataTemp,"E:/canalizationENV/ppt/velocity_var_by_gene_padj_smaller_larger_total.csv")

#

var_sigCell<-lapply(1:dim(velocity_var_by_gene_padj)[1],function(i){
	ind<-which(velocity_var_by_gene_padj[i,]<0.05)
	if(length(ind)!=0){
		return(data.frame(ind=ind,var_rel=round(as.numeric(velocity_var_by_gene_relative_refine[i,ind]),3)))
	}
	
})
names(var_sigCell)<-rownames(velocity_var_by_gene_padj)
var_sigCell<-var_sigCell[which(sapply(var_sigCell,length)>0)]

sort(table(unlist(lapply(var_sigCell,function(x){rownames(x)}))))


var_sigCell_p<-lapply(1:dim(velocity_var_by_gene_padj)[1],function(i){
	ind<-which(velocity_var_by_gene_padj[i,]<0.05)
	if(length(ind)!=0){
		return(data.frame(ind=ind,var_rel=as.numeric(velocity_var_by_gene_padj[i,ind])))
	}
	
})
names(var_sigCell_p)<-rownames(velocity_var_by_gene_padj)
var_sigCell_p<-var_sigCell_p[which(sapply(var_sigCell_p,length)>0)]

var_sigCell_kdiff<-lapply(1:dim(velocity_var_by_gene_padj)[1],function(i){
	ind<-which(velocity_var_by_gene_padj[i,]<0.05&velocity_var_by_gene_relative_refine[i,]>1)
	if(length(ind)!=0){
		return(data.frame(ind=ind,k_diff=as.numeric(k_diff_matrix_refine[i,ind])))
	}
	
})
names(var_sigCell_kdiff)<-rownames(velocity_var_by_gene_padj)
var_sigCell_kdiff<-var_sigCell_kdiff[which(sapply(var_sigCell_kdiff,function(x){ifelse(is.null(x),0,dim(x)[1])})>0)]
var_sigCell_kdiff
hist(do.call(rbind,var_sigCell_kdiff)[,2])

col1<-rep("grey",2039)
col1[which(gene_extend=="ACS-13")]<-"red"
dataTemp<-data.frame(XM=q1_cell350_alignment_xyz[,"ABpraaaap_x"],
					 XD=q1_cell350_alignment_xyz[,"ABpraaaapa_x"],
					 delta_X=q1_cell350_alignment_xyz[,"ABpraaaapa_x"]-
							 q1_cell350_alignment_xyz[,"ABpraaaap_x"],
					 col1=col1)
ggplot(dataTemp,aes(x=XM,y=XD,col=col1))+geom_point()+geom_smooth(method="lm")
ggplot(dataTemp,aes(x=XM,y=delta_X,col=col1))+geom_point()+geom_smooth(method="lm")


dataTemp<-data.frame(Var_ratio=do.call(rbind,var_sigCell)[,2],Var_ratio_pvalue=do.call(rbind,var_sigCell_p)[,2],k_diff=do.call(rbind,var_sigCell_kdiff)[,2])

ggplot(dataTemp,aes(x=Var_ratio,y=-log10(Var_ratio_pvalue)))+geom_point()+
geom_vline(xintercept=1,col="red")+
scale_x_log10()+
labs(x="log10(variance ratio between gene/env Vs wild types)",y="-log10(p)")

hist(dataTemp$k_diff)







var_sigCell
hist(unlist(var_sigCell))


dataTemp<-velocity_var_by_gene_relative_refine[apply(velocity_var_by_gene_padj,1,function(x){length(which(x<0.05))})>=3,]
dataTemp2<-data.frame(variance=apply(dataTemp,1,median),gene=rep(rownames(dataTemp),each=1))
ggplot(dataTemp2,aes(x=gene,y=variance))+
geom_point()+
geom_smooth()+
geom_hline(yintercept=1)

dataTemp3<-data.frame(log_variance=log(c(t(dataTemp))),gene=rep(rownames(dataTemp),each=1059))
ggplot(dataTemp3,aes(x=gene,y=log_variance))+
geom_boxplot()+
geom_hline(yintercept=0)

hist(velocity_var_by_gene_relative_refine[velocity_var_by_gene_padj<0.05])

which(apply(velocity_var_by_gene_padj,2,function(x){length(which(x<0.05))})>0)


velocity_var_by_gene_relative_one<-velocity_var_by_gene_relative_refine
velocity_var_by_gene_relative_one[which(velocity_var_by_gene_padj>0.05)]=1


custom_palette <- colorRampPalette(c("blue", "white", "red"))(50)
max_abs <- max(abs(log(velocity_var_by_gene_relative_refine[-1,1:353])))  # Maximum absolute value
breaks <- seq(-max_abs, max_abs, length.out = 51)
# Heatmap
pheatmap(
  log(velocity_var_by_gene_relative_refine[-1,]),
  color = custom_palette,
  breaks = breaks,
  main = "log(velocity variance of gene or env / variance of wild types)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = F, show_colnames = F
)
grid.text("gene/env", x = 0.02, y = 0.5, rot = 90, gp = gpar(fontsize = 12))  # Row axis title (left side)
grid.text("cell features", x = 0.5, y = 0.02, gp = gpar(fontsize = 12))          # Column axis title (bottom)

pheatmap(
  log(velocity_var_by_gene_relative_one[-1,]),
  color = custom_palette,
  breaks = breaks,
  main = "log(velocity variance of gene or env / variance of wild types)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = F, show_colnames = F
)
grid.text("gene/env", x = 0.02, y = 0.5, rot = 90, gp = gpar(fontsize = 12))  # Row axis title (left side)
grid.text("cell features", x = 0.5, y = 0.02, gp = gpar(fontsize = 12))          # Column axis title (bottom)






par(mfrow=c(1,3))
hist(velocity_var_by_gene_relative_refine[which(velocity_var_by_gene_padj<0.05)],1000)
hist(velocity_var_by_gene_relative_refine[which(velocity_var_by_gene_padj>0.05)],1000)
hist(log(velocity_var_by_gene_relative_one[-1,]))

plot(velocity_var_by_gene_relative_refine[,1:353],velocity_var_by_gene_relative_refine[,c(1:353)+353])


plot(velocity_var_by_gene_refine[1,],velocity_var_by_gene_refine[20,])



######################################################################
# Plot 
######################################################################

# velocity variance heatmap
library(pheatmap)
sig_cell<-colnames(velocity_var_by_gene_padj)[which(apply(velocity_var_by_gene_padj,2,function(x){length(which(x<0.05))})>0)]
pheatmap(-log10(velocity_var_by_gene_padj[match(sig_gene,rownames(velocity_var_by_gene_padj)),match(sig_cell,colnames(velocity_var_by_gene_padj))]))

pheatmap(-log10(velocity_var_by_gene_padj[match(selected_genes,rownames(velocity_var_by_gene_padj)),match(sig_cell,colnames(velocity_var_by_gene_padj))]))

table(sapply(sig_cell,function(x){unlist(strsplit(x,"_"))[2]}))


pheatmap(-log10(velocity_var_by_gene_padj[c(match(selected_genes,rownames(velocity_var_by_gene_padj)),setdiff(1:dim(velocity_var_by_gene_padj)[1],match(selected_genes,rownames(velocity_var_by_gene_padj)))),
	c(match(sig_cell,colnames(velocity_var_by_gene_padj)),setdiff(1:dim(velocity_var_by_gene_padj)[2],match(sig_cell,colnames(velocity_var_by_gene_padj))))]),cluster_rows=FALSE,cluster_cols=FALSE)


















