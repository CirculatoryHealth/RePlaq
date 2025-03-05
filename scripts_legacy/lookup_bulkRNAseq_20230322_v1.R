# Lookup in bulkRNAseq

#packages
library("ggplot2")
library("plotly")
library("readxl")
library("dplyr")
library("purrr")
library("tidyverse")
library("knitr")
library("janitor")
require(tidyverse)
require(haven)
require(DESeq2)
require(tidyverse) #to clean and change your dataset
require(tableone) # for creating baseline tables
library(ggplot2) # to make nice graphs
library(knitr) #to visualize tables in knitted output
require(data.table)
require(GenomicRanges)
library(org.Hs.eg.db)
require(ggpubr)
require(cowplot)
require(fabricatr)
require(Seurat)



lm_eqn <- function(df, y, x){
  formula = as.formula(sprintf('%s ~ %s', y, x))
  m <- lm(formula, data=df);
  # formating the values into a summary string to print out
  # ~ give some space, but equal size and comma need to be quoted
  if(coef(m)[2]>=0){
    eq <- substitute(italic(target) == a + b %.% italic(input)*","~~p~"="~italic(pvalue), 
                     list(target = y,
                          input = x,
                          a = format(as.vector(coef(m)[1]), digits = 2), 
                          b = format(as.vector(coef(m)[2]), digits = 2), 
                          # getting the pvalue is painful
                          pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=1)
                     )
    )
  }else{
    eq <- substitute(italic(target) == a - b %.% italic(input)*","~~p~"="~italic(pvalue), 
                     list(target = y,
                          input = x,
                          a = format(as.vector(coef(m)[1]), digits = 2), 
                          b = format(as.vector(abs(coef(m)[2])), digits = 2), 
                          # getting the pvalue is painful
                          pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=1)
                     )
    )
  }
  as.character(as.expression(eq));                 
}

# Loading data
load(file = "C://Users/ediez/Documents/Students/Kristina Dimitrova/Bulk_deseq2_object.Rdata")
load(file = "C://Users/ediez/Documents/Students/Kristina Dimitrova/Norm_count_deseq2_AERNA_no_duplicates.Rdata")
load(file = "C://Users/ediez/Documents/Students/Kristina Dimitrova/Norm_count_deseq2_AERNA_mean.Rdata")

# Formatting variables
tail(colnames(norm_count_met_s2),49)
norm_count_met_s2$Symptoms.5G<-factor(norm_count_met_s2$Symptoms.5G,levels=c("Asymptomatic","Other","Ocular","Retinal infarction","TIA","Stroke"))
norm_count_met_s2$Sex<-as.character(norm_count_met_s2$Sex)
norm_count_met_s2$Sex[which(norm_count_met_s2$Sex=="male")]<-"Male"
norm_count_met_s2$Sex[which(norm_count_met_s2$Sex=="female")]<-"Female"
norm_count_met_s2$Sex<-factor(norm_count_met_s2$Sex,levels=c("Female","Male"))
norm_count_met_s2$log_hsCRP_plasma<-log(norm_count_met_s2$hsCRP_plasma)
norm_count_met_s2$SmokerStatus<-factor(norm_count_met_s2$SmokerStatus,levels=c("Never smoked","Ex-smoker","Current smoker"))
norm_count_met_s2$OverallPlaquePhenotype<-as.character(norm_count_met_s2$OverallPlaquePhenotype)
norm_count_met_s2$OverallPlaquePhenotype[which(norm_count_met_s2$OverallPlaquePhenotype=="atheromatous")]<-"Atheromatous"
norm_count_met_s2$OverallPlaquePhenotype[which(norm_count_met_s2$OverallPlaquePhenotype=="fibrous")]<-"Fibrous"
norm_count_met_s2$OverallPlaquePhenotype[which(norm_count_met_s2$OverallPlaquePhenotype=="fibroatheromatous")]<-"Fibroatheromatous"
norm_count_met_s2$OverallPlaquePhenotype<-factor(norm_count_met_s2$OverallPlaquePhenotype,levels=c("Fibrous","Fibroatheromatous","Atheromatous"))
norm_count_met_s2$Plaque_Vulnerability_Index_num<-as.numeric(norm_count_met_s2$Plaque_Vulnerability_Index)
norm_count_met_s2$cluster[which(norm_count_met_s2$cluster==0)]<-"Fibro-collagenous"
norm_count_met_s2$cluster[which(norm_count_met_s2$cluster==1)]<-"Intermediate"
norm_count_met_s2$cluster[which(norm_count_met_s2$cluster==2)]<-"Lipomatous"
norm_count_met_s2$cluster[which(norm_count_met_s2$cluster==3)]<-"Fibro-inflammatory"
norm_count_met_s2$cluster[which(norm_count_met_s2$cluster==4)]<-"Fibro-cellular"
norm_count_met_s2$cluster<-factor(norm_count_met_s2$cluster)
# mean_normalized<-apply(norm_count_met_s2[,1:(dim(norm_count_met_s2)[2]-51)],2,function(x){mean(x,na.rm=TRUE)})
# save(mean_normalized,file = "C://Users/ediez/Documents/Students/Kristina Dimitrova/Norm_count_deseq2_AERNA_mean.Rdata")

deciles<-split_quantile(x = mean_normalized, type = 10)
data_dec<-as.data.frame(cbind(mean_normalized,deciles))
row.names(data_dec)<-names(mean_normalized)

# Plot for given genes
my_comp_sex <- list( c("Female", "Male"))
my_comp_S5G <- list( c("Asymptomatic", "Other"),c("Asymptomatic", "Ocular"),
                     c("Asymptomatic", "Retinal infarction"), c("Asymptomatic", "TIA"),
                     c("Asymptomatic", "Stroke"))
my_comp_EPc<- list(c("No composite endpoints","Composite endpoints"))
my_comp_EPm<- list(c("No major events (endpoints)","Major events (endpoints)"))
my_comp_Smoker<- list(c("Never smoked","Ex-smoker"),c("Never smoked","Current smoker"))
my_comp_plaque<-list(c("no/minor","moderate/heavy"))
my_comp_fat40<-list(c("<40%",">40%"))
my_comp_iph<-list(c("no","yes"))
my_comp_phen<-list(c("Fibrous","Fibroatheromatous"),c("Fibrous","Atheromatous"))

load(file = "C://Documents and Settings/ediez/Documents/PLINK/_AE_ORIGINALS/scRNAseq/20210811.46.patients.Koen.RDS")
pan_data<-readRDS(file = "C://Users/ediez/Documents/PLINK/_AE_ORIGINALS/scRNAseq/Pan_et_al_human_EDB.RDS")
wirka_pv<-readRDS(file = "C://Users/ediez/Documents/PLINK/_AE_ORIGINALS/scRNAseq/Wirka_2019_plaqview.rds")
seuset_emt<-readRDS(file = "C://Documents and Settings/ediez/Documents/PLINK/_AE_ORIGINALS/scRNAseq/20210812.46.patients.Koen.cleaned.endoMT.subset.RDS")
seuset_emt@active.ident<-factor(seuset_emt@active.ident,levels=c("EC1","EC2","EC3","EC4","SMC1","SMC2","SMC3","SMC4"))
seuset_emt_ec<-subset(seuset_emt,idents=c("EC1","EC2","EC3","EC4"))
seuset_emt_smc<-subset(seuset_emt,idents=c("SMC1","SMC2","SMC3","SMC4"))
seuset_emt_smc@active.ident<-factor(seuset_emt_smc@active.ident,levels = c("SMC1","SMC2","SMC3","SMC4"))

seuset$Sex[which(seuset$Sex=="male")]<-"Male"
seuset$Sex[which(seuset$Sex=="female")]<-"Female"


seuset_emt$Sex[which(seuset_emt$Sex=="male")]<-"Male"
seuset_emt$Sex[which(seuset_emt$Sex=="female")]<-"Female"


seuset$clusters_ct<-seuset@active.ident
seuset$clusters_ct<-as.character(seuset$clusters_ct)
seuset$clusters_ct[grep('KIT',seuset$clusters_ct)]<-"Mast Cells"
seuset$clusters_ct[grep("CD68\\+CD4\\+",seuset$clusters_ct)]<-"Monocytes"
seuset$clusters_ct[grep("CD68\\+ABCA",seuset$clusters_ct)]<-"Foam Cells"
seuset$clusters_ct[grep("CD68\\+IL18",seuset$clusters_ct)]<-"Resident Macrophages"
seuset$clusters_ct[grep("CD68\\+CASP1",seuset$clusters_ct)]<-"Inflammatory Macrophages"
seuset$clusters_ct[grep("CD68\\+CD1C",seuset$clusters_ct)]<-"Dendritic Cells"
seuset$clusters_ct[grep("ACTA2\\+",seuset$clusters_ct)]<-"SMCs"
seuset$clusters_ct[grep("CD3\\+CD56\\+",seuset$clusters_ct)]<-"NK-Cells"
seuset$clusters_ct[grep("CD34\\+ Endothelial Cells II",seuset$clusters_ct)]<-"ECs II"
seuset$clusters_ct[grep("CD34\\+ Endothelial Cells I",seuset$clusters_ct)]<-"ECs I"
seuset$clusters_ct[grep("CD3\\+",seuset$clusters_ct)]<-"T-Cells"
seuset$clusters_ct[grep("FOXP3\\+",seuset$clusters_ct)]<-"T-Cells"
seuset$clusters_ct[grep("CD79A",seuset$clusters_ct)]<-"Memory B-Cells"
seuset$clusters_ct[grep("CD79",seuset$clusters_ct)]<-"Plasma B-Cells"
seuset$clusters_ct<-factor(seuset$clusters_ct,levels=c("SMCs","ECs I","ECs II",
                                                       "Foam Cells","Resident Macrophages","Inflammatory Macrophages","Dendritic Cells","Monocytes","T-Cells","NK-Cells","Memory B-Cells","Plasma B-Cells",
                                                       "Mast Cells"))



list_genes<-c("ACTA2","IGFBP7")



for(gene_i in list_genes){
  
  gene_int<-gene_i
  
  if(!(gene_int%in%colnames(norm_count_met_s2))){
    print(paste("Gene ",gene_int," is not present in the bulkRNAseq dataset.",sep=""))
    mainDir <- "C:/Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq"
    subDir <- gene_int
    setwd(file.path(mainDir))
    if (file.exists(subDir)){
      setwd(file.path(mainDir, subDir))
    } else {
      dir.create(file.path(mainDir, subDir))
      setwd(file.path(mainDir, subDir))
      
    }
  }else{
    mainDir <- "C:/Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq"
    subDir <- gene_int
    setwd(file.path(mainDir))
    if (file.exists(subDir)){
      setwd(file.path(mainDir, subDir))
    } else {
      dir.create(file.path(mainDir, subDir))
      setwd(file.path(mainDir, subDir))
      
    }
    
    color_barplot<-rep("darkgrey",10)
    color_barplot[data_dec[gene_int,"deciles"]]<-"black"
    fill_barplot<-rep("lightgrey",10)
    fill_barplot[data_dec[gene_int,"deciles"]]<-"darkred"
    mean_gene<-data_dec[gene_int,"mean_normalized"]
    
    p_dec<-ggboxplot(data_dec,x="deciles",y="mean_normalized",color = color_barplot,fill=fill_barplot) + labs(x="Expression deciles (~5520 genes per decile)",y="Log Norm Expression") +
      geom_hline(yintercept = mean_gene, linetype = 2,color="darkred") +
      annotate(geom="text", x=1, y=mean_gene+1, label=paste("Mean expression ",gene_int," = ", round(mean_gene,digits=2),sep=""), color="darkred",hjust=0)
    
    p_age<-ggscatter(x="Age",y=gene_int,alpha = 0.05,conf.int = TRUE, data=norm_count_met_s2,add="reg.line",add.params = list(color="red",fill = "lightgray"), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="Age (years)")+
      annotate(geom="text", x=45, y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2,gene_int,"Age"), color="red",hjust=0,parse=T)
    
    p_age_sex<-ggscatter(x="Age",y=gene_int,alpha = 0.05, data=norm_count_met_s2,add="reg.line",color="Sex", palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="Age (years)")+
      annotate(geom="text", x=45, y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Female"),],gene_int,"Age"), color=get_palette("lancet",2)[1],hjust=0,parse=T)+
      annotate(geom="text", x=45, y=(0.10*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Male"),],gene_int,"Age"), color=get_palette("lancet",2)[2],hjust=0,parse=T)
    
    p_sex<-ggboxplot(x="Sex",y=gene_int,data=norm_count_met_s2,fill="Sex",add="jitter",add.params=list(alpha=0.05), palette ="lancet")+labs(y=paste("Log norm ",gene_int,sep="")) + 
      stat_compare_means(comparisons = my_comp_sex,label = "p.signif")+NoLegend()
    
    p_sympt<-ggboxplot(x="Symptoms.5G",y=gene_int,
                       data=norm_count_met_s2[-which(is.na(norm_count_met_s2$Symptoms.5G)),],fill="Symptoms.5G",xlab = "Symptoms at presentation",add="jitter",add.params=list(alpha=0.05), palette ="npg")+
      labs(y=paste("Log norm ",gene_int,sep="")) + stat_compare_means(comparisons = my_comp_S5G,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_cluster<-ggboxplot(x="cluster",y=gene_int,
                         data=norm_count_met_s2[-which(is.na(norm_count_met_s2$cluster)),],fill="cluster",xlab = "Cluster",add="jitter",add.params=list(alpha=0.05), palette ="nejm")+
      labs(y=paste("Log norm ",gene_int,sep="")) +NoLegend()+
      stat_compare_means(method = "anova")+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_epc<-ggboxplot(x="EP_composite",y=gene_int,
                     data=norm_count_met_s2[-which(is.na(norm_count_met_s2$EP_composite)),],fill="EP_composite",add="jitter",add.params=list(alpha=0.05), palette ="npg")+
      labs(y=paste("Log norm ",gene_int,sep="")) + stat_compare_means(comparisons = my_comp_EPc,label = "p.signif")+NoLegend()+theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    p_epm<-ggboxplot(x="EP_major",y=gene_int,
                     data=norm_count_met_s2[-which(is.na(norm_count_met_s2$EP_major)),],fill="EP_major",add="jitter",add.params=list(alpha=0.05), palette ="npg")+
      labs(y=paste("Log norm ",gene_int,sep="")) + stat_compare_means(comparisons = my_comp_EPm,label = "p.signif")+NoLegend()+theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    plots_top<-plot_grid(p_dec,p_age,ncol=2,labels="AUTO",label_size = 14)
    plots_middle<-plot_grid(p_sex,p_age_sex,ncol=2,labels=c("C","D"),label_size = 14)
    plots_middle2<-plot_grid(p_sympt,p_cluster,ncol=2,labels=c("E","F"),label_size = 14)
    plots_bottom<-plot_grid(p_epc,p_epm,ncol=2,labels=c("G","H"),label_size = 14)
    plots_all<-plot_grid(plots_top,plots_middle,plots_middle2,plots_bottom,ncol=1,labels=c("","","",""),label_size = 14)
    
    tiff(paste("C://Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq/",gene_int,"/Plot_1_general_",gene_int,".tiff",sep=""),width=550,height=1000)
    print(plots_all)
    dev.off()
    
    p_ldl<-ggscatter(x="LDL_final",y=gene_int,alpha = 0.05,conf.int = TRUE, data=norm_count_met_s2,add="reg.line",add.params = list(color="red",fill = "lightgray"), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="LDL (mmol/L)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$LDL_final,na.rm=TRUE)*0.25), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2,gene_int,"LDL_final"), color="red",hjust=0,parse=T)
    
    p_ldl_sex<-ggscatter(x="LDL_final",y=gene_int,alpha = 0.05, data=norm_count_met_s2,add="reg.line",color="Sex", palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="LDL (mmol/L)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$LDL_final,na.rm=TRUE)*0.25), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Female"),],gene_int,"LDL_final"), color=get_palette("lancet",2)[1],hjust=0,parse=T)+
      annotate(geom="text", x=(mean(norm_count_met_s2$LDL_final,na.rm=TRUE)*0.25), y=(0.10*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Male"),],gene_int,"LDL_final"), color=get_palette("lancet",2)[2],hjust=0,parse=T)
    
    p_hdl<-ggscatter(x="LDL_final",y=gene_int,alpha = 0.05,conf.int = TRUE, data=norm_count_met_s2,add="reg.line",add.params = list(color="red",fill = "lightgray"), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="HDL (mmol/L)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$HDL_final,na.rm=TRUE)*0.25), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2,gene_int,"HDL_final"), color="red",hjust=0,parse=T)
    
    p_hdl_sex<-ggscatter(x="HDL_final",y=gene_int,alpha = 0.05, data=norm_count_met_s2,add="reg.line",color="Sex", palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="HDL (mmol/L)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$HDL_final,na.rm=TRUE)*0.25), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Female"),],gene_int,"HDL_final"), color=get_palette("lancet",2)[1],hjust=0,parse=T)+
      annotate(geom="text", x=(mean(norm_count_met_s2$HDL_final,na.rm=TRUE)*0.25), y=(0.10*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Male"),],gene_int,"HDL_final"), color=get_palette("lancet",2)[2],hjust=0,parse=T)
    
    p_sbp<-ggscatter(x="systolic",y=gene_int,alpha = 0.05,conf.int = TRUE, data=norm_count_met_s2,add="reg.line",add.params = list(color="red",fill = "lightgray"), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="Systolic BP (mmHg)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$systolic,na.rm=TRUE)*0.65), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2,gene_int,"systolic"), color="red",hjust=0,parse=T)
    
    p_sbp_sex<-ggscatter(x="systolic",y=gene_int,alpha = 0.05, data=norm_count_met_s2,add="reg.line",color="Sex", palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="Systolic BP (mmHg)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$systolic,na.rm=TRUE)*0.65), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Female"),],gene_int,"systolic"), color=get_palette("lancet",2)[1],hjust=0,parse=T)+
      annotate(geom="text", x=(mean(norm_count_met_s2$systolic,na.rm=TRUE)*0.65), y=(0.10*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Male"),],gene_int,"systolic"), color=get_palette("lancet",2)[2],hjust=0,parse=T)
    
    p_dbp<-ggscatter(x="diastoli",y=gene_int,alpha = 0.05,conf.int = TRUE, data=norm_count_met_s2,add="reg.line",add.params = list(color="red",fill = "lightgray"), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="Diastolic BP (mmHg)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$diastoli,na.rm=TRUE)*0.65), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2,gene_int,"diastoli"), color="red",hjust=0,parse=T)
    
    p_dbp_sex<-ggscatter(x="diastoli",y=gene_int,alpha = 0.05, data=norm_count_met_s2,add="reg.line",color="Sex", palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="Diastolic BP (mmHg)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$diastoli,na.rm=TRUE)*0.65), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Female"),],gene_int,"diastoli"), color=get_palette("lancet",2)[1],hjust=0,parse=T)+
      annotate(geom="text", x=(mean(norm_count_met_s2$diastoli,na.rm=TRUE)*0.65), y=(0.10*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Male"),],gene_int,"diastoli"), color=get_palette("lancet",2)[2],hjust=0,parse=T)
    
    
    plots_top<-plot_grid(p_ldl,p_ldl_sex,ncol=2,labels="AUTO",label_size = 14)
    plots_middle<-plot_grid(p_hdl,p_hdl_sex,ncol=2,labels=c("C","D"),label_size = 14)
    plots_bottom<-plot_grid(p_sbp,p_sbp_sex,ncol=2,labels=c("E","F"),label_size = 14)
    plots_bottom2<-plot_grid(p_dbp,p_dbp_sex,ncol=2,labels=c("G","H"),label_size = 14)
    plots_all2<-plot_grid(plots_top,plots_middle,plots_bottom,plots_bottom2,ncol=1,labels=c("","","",""),label_size = 14)
    
    tiff(paste("C://Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq/",gene_int,"/Plot_2_Risk_factors_",gene_int,".tiff",sep=""),width=550,height=1000)
    print(plots_all2)
    dev.off()
    
    p_crp<-ggscatter(x="log_hsCRP_plasma",y=gene_int,alpha = 0.05,conf.int = TRUE, data=norm_count_met_s2,add="reg.line",add.params = list(color="red",fill = "lightgray"), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="log High-sensitivity CRP (mg/L)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$log_hsCRP_plasma,na.rm=TRUE))*(-2), y=(0.25*(mean_gene)), label=gsub("log_","log ",gsub("_plasma","",lm_eqn(norm_count_met_s2,gene_int,"log_hsCRP_plasma"))), color="red",hjust=0,parse=T)
    
    p_crp_sex<-ggscatter(x="log_hsCRP_plasma",y=gene_int,alpha = 0.05, data=norm_count_met_s2,add="reg.line",color="Sex", palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x="log High-sensitivity CRP (mg/L)")+
      annotate(geom="text", x=(mean(norm_count_met_s2$log_hsCRP_plasma,na.rm=TRUE))*(-2), y=(0.25*(mean_gene)), label=gsub("log_","log ",gsub("_plasma","",lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Female"),],gene_int,"log_hsCRP_plasma"))), color=get_palette("lancet",2)[1],hjust=0,parse=T)+
      annotate(geom="text", x=(mean(norm_count_met_s2$log_hsCRP_plasma,na.rm=TRUE))*(-2), y=(0.10*(mean_gene)), label=gsub("log_","log ",gsub("_plasma","",lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Male"),],gene_int,"log_hsCRP_plasma"))), color=get_palette("lancet",2)[2],hjust=0,parse=T)
    
    p_BMI<-ggscatter(x="BMI",y=gene_int,alpha = 0.05,conf.int = TRUE, data=norm_count_met_s2,add="reg.line",add.params = list(color="red",fill = "lightgray"), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("BMI (kg/m"^2~")"))+
      annotate(geom="text", x=(mean(norm_count_met_s2$BMI,na.rm=TRUE)*0.75), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2,gene_int,"BMI"), color="red",hjust=0,parse=T)
    
    p_BMI_sex<-ggscatter(x="BMI",y=gene_int,alpha = 0.05, data=norm_count_met_s2,add="reg.line",color="Sex", palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("BMI (kg/m"^2~")"))+
      annotate(geom="text", x=(mean(norm_count_met_s2$BMI,na.rm=TRUE)*0.75), y=(0.25*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Female"),],gene_int,"BMI"), color=get_palette("lancet",2)[1],hjust=0,parse=T)+
      annotate(geom="text", x=(mean(norm_count_met_s2$BMI,na.rm=TRUE)*0.75), y=(0.10*(mean_gene)), label=lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Male"),],gene_int,"BMI"), color=get_palette("lancet",2)[2],hjust=0,parse=T)
    
    p_smo<-ggboxplot(x="SmokerStatus",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$SmokerStatus)),],add="jitter",fill="SmokerStatus",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Smoking status"))+
      stat_compare_means(comparisons = my_comp_Smoker,label = "p.signif")+NoLegend()+theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_smo_sex<-ggboxplot(x="SmokerStatus",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$SmokerStatus)),],facet.by = "Sex",add="jitter",fill="SmokerStatus",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Smoking status"))+
      stat_compare_means(comparisons = my_comp_Smoker,label = "p.signif")+NoLegend()+theme(axis.title.x=element_blank(),axis.text.x = element_text(size=9,angle = 20,vjust=0.5,hjust = 0.5))
    
    plots_top<-plot_grid(p_crp,p_crp_sex,ncol=2,labels="AUTO",label_size = 14)
    plots_middle<-plot_grid(p_BMI,p_BMI_sex,ncol=2,labels=c("C","D"),label_size = 14)
    plots_bottom<-plot_grid(p_smo,p_smo_sex,ncol=2,labels=c("E","F"),label_size = 14)
    plots_all3<-plot_grid(plots_top,plots_middle,plots_bottom,ncol=1,labels=c("","",""),label_size = 14)
    
    tiff(paste("C://Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq/",gene_int,"/Plot_3_Risk_factors_",gene_int,".tiff",sep=""),width=550,height=750)
    print(plots_all3)
    dev.off()
    
    p_phen<-ggboxplot(x="OverallPlaquePhenotype",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$OverallPlaquePhenotype)),],add="jitter",fill="OverallPlaquePhenotype",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Plaque phenotype"))+
      stat_compare_means(comparisons = my_comp_phen,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_phen_sex<-ggboxplot(x="OverallPlaquePhenotype",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$OverallPlaquePhenotype)),],facet.by = "Sex",add="jitter",fill="OverallPlaquePhenotype",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Plaque phenotype"))+
      stat_compare_means(comparisons = my_comp_phen,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(size=9,angle = 20,vjust=0.5,hjust = 0.5))
    
    p_pvi<-ggboxplot(x="Plaque_Vulnerability_Index",y=gene_int, data=norm_count_met_s2,add="jitter",fill="Plaque_Vulnerability_Index",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Plaque Vulnerability Index"))+
      annotate(geom="text", x=(mean(as.numeric(norm_count_met_s2$Plaque_Vulnerability_Index),na.rm=TRUE)*0.25), y=(0.25*(mean_gene)), label=gsub("Plaque_Vulnerability_Index_num","PVI",lm_eqn(norm_count_met_s2,gene_int,"Plaque_Vulnerability_Index_num")), color="red",hjust=0,parse=T)+
      NoLegend()+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_pvi_sex<-ggboxplot(x="Plaque_Vulnerability_Index",y=gene_int, data=norm_count_met_s2,facet.by = "Sex",add="jitter",fill="Plaque_Vulnerability_Index",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Plaque Vulnerability Index"))+
      NoLegend()+theme(axis.text.x = element_text(size=9,angle = 20,vjust=0.5,hjust = 0.5))
    
    ann_text<-data.frame(gene_int=c((mean(as.numeric(norm_count_met_s2$Plaque_Vulnerability_Index,na.rm=TRUE))*0.10),
                                    (mean(as.numeric(norm_count_met_s2$Plaque_Vulnerability_Index,na.rm=TRUE))*0.10)),
                         Plaque_Vulnerability_Index=c(1,1),Sex=c("Male","Female"),label=c(gsub("Plaque_Vulnerability_Index_num","PVI",lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Female"),],gene_int,"Plaque_Vulnerability_Index_num")),gsub("Plaque_Vulnerability_Index_num","PVI",lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Male"),],gene_int,"Plaque_Vulnerability_Index_num"))))
    colnames(ann_text)[1]<-gene_int
    p_pvi_sex<-p_pvi_sex + geom_text(data = ann_text,label=ann_text$label,hjust=0,parse = T,size=2)
    
    p_smc<-ggboxplot(x="SMC.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$SMC.bin)),],add="jitter",fill="SMC.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression(alpha~"-SMA content"))+
      stat_compare_means(comparisons = my_comp_plaque,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_smc_sex<-ggboxplot(x="SMC.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$SMC.bin)),],facet.by = "Sex",add="jitter",fill="SMC.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression(alpha~"-SMA content"))+
      stat_compare_means(comparisons = my_comp_plaque,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(size=9,angle = 20,vjust=0.5,hjust = 0.5))
    
    p_mac<-ggboxplot(x="Macrophages.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$Macrophages.bin)),],add="jitter",fill="Macrophages.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression(paste("CD68"^+ ~ ""," Macrophage content",sep="")))+
      stat_compare_means(comparisons = my_comp_plaque,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_mac_sex<-ggboxplot(x="Macrophages.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$Macrophages.bin)),],facet.by = "Sex",add="jitter",fill="Macrophages.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression(paste("CD68"^+ ~ ""," Macrophage content",sep="")))+
      stat_compare_means(comparisons = my_comp_plaque,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(size=9,angle = 20,vjust=0.5,hjust = 0.5))
    
    
    
    plots_top<-plot_grid(p_phen,p_phen_sex,ncol=2,labels="AUTO",label_size = 14)
    plots_middle<-plot_grid(p_pvi,p_pvi_sex,ncol=2,labels=c("C","D"),label_size = 14)
    plots_bottom<-plot_grid(p_smc,p_smc_sex,ncol=2,labels=c("E","F"),label_size = 14)
    plots_bottom2<-plot_grid(p_mac,p_mac_sex,ncol=2,labels=c("G","H"),label_size = 14)
    plots_all4<-plot_grid(plots_top,plots_middle,plots_bottom,plots_bottom2,ncol=1,labels=c("","","",""),label_size = 14)
    
    tiff(paste("C://Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq/",gene_int,"/Plot_4_Plaque_characteristics_",gene_int,".tiff",sep=""),width=550,height=1000)
    print(plots_all4)
    dev.off()
    
    p_fat<-ggboxplot(x="Fat.bin_40",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$Fat.bin_40)),],add="jitter",fill="Fat.bin_40",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Fat content"))+
      stat_compare_means(comparisons = my_comp_fat40,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_fat_sex<-ggboxplot(x="Fat.bin_40",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$Fat.bin_40)),],facet.by = "Sex",add="jitter",fill="Fat.bin_40",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Fat content"))+
      stat_compare_means(comparisons = my_comp_fat40,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(size=9,angle = 20,vjust=0.5,hjust = 0.5))
    
    p_iph<-ggboxplot(x="IPH.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$IPH.bin)),],add="jitter",fill="IPH.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Plaque hemmoraghe"))+
      stat_compare_means(comparisons = my_comp_iph,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_iph_sex<-ggboxplot(x="IPH.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$IPH.bin)),],facet.by = "Sex",add="jitter",fill="IPH.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Plaque hemmoraghe"))+
      stat_compare_means(comparisons = my_comp_iph,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(size=9,angle = 20,vjust=0.5,hjust = 0.5))
    
    p_calc<-ggboxplot(x="Calc.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$Calc.bin)),],add="jitter",fill="Calc.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Calcification"))+
      stat_compare_means(comparisons = my_comp_plaque,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_calc_sex<-ggboxplot(x="Calc.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$Calc.bin)),],facet.by = "Sex",add="jitter",fill="Calc.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Calcification"))+
      stat_compare_means(comparisons = my_comp_plaque,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(size=9,angle = 20,vjust=0.5,hjust = 0.5))
    
    p_coll<-ggboxplot(x="Collagen.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$Collagen.bin)),],add="jitter",fill="Collagen.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Collagen content"))+
      stat_compare_means(comparisons = my_comp_plaque,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(angle = 20,vjust=0.5,hjust = 0.5))
    
    p_coll_sex<-ggboxplot(x="Collagen.bin",y=gene_int, data=norm_count_met_s2[-which(is.na(norm_count_met_s2$Collagen.bin)),],facet.by = "Sex",add="jitter",fill="Collagen.bin",add.params = list(alpha=0.02), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Collagen content"))+
      stat_compare_means(comparisons = my_comp_plaque,label = "p.signif")+NoLegend()+theme(axis.text.x = element_text(size=9,angle = 20,vjust=0.5,hjust = 0.5))
    
    p_neo<-ggscatter(x="VesselDensity_rankNorm",y=gene_int,alpha = 0.05,conf.int = TRUE, data=norm_count_met_s2,add="reg.line",add.params = list(color="red",fill = "lightgray"), palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Vessel Density rank-norm"))+
      annotate(geom="text", x=(mean(norm_count_met_s2$VesselDensity_rankNorm,na.rm=TRUE)*-15), y=(0.25*(mean_gene)), label=gsub("VesselDensity_rankNorm","Vessel Dens",lm_eqn(norm_count_met_s2,gene_int,"VesselDensity_rankNorm")), color="red",hjust=0,parse=T)
    
    p_neo_sex<-ggscatter(x="VesselDensity_rankNorm",y=gene_int,alpha = 0.05, data=norm_count_met_s2,add="reg.line",color="Sex", palette ="lancet")+
      labs(y=paste("Log norm ",gene_int,sep=""),x=expression("Vessel Density rank-norm"))+
      annotate(geom="text", x=(mean(norm_count_met_s2$VesselDensity_rankNorm,na.rm=TRUE)*-15), y=(0.25*(mean_gene)), label=gsub("VesselDensity_rankNorm","Vessel Dens",lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Female"),],gene_int,"VesselDensity_rankNorm")), color=get_palette("lancet",2)[1],hjust=0,parse=T)+
      annotate(geom="text", x=(mean(norm_count_met_s2$VesselDensity_rankNorm,na.rm=TRUE)*-15), y=(0.10*(mean_gene)), label=gsub("VesselDensity_rankNorm","Vessel Dens",lm_eqn(norm_count_met_s2[which(norm_count_met_s2$Sex=="Male"),],gene_int,"VesselDensity_rankNorm")), color=get_palette("lancet",2)[2],hjust=0,parse=T)
    
    
    plots_top<-plot_grid(p_fat,p_fat_sex,ncol=2,labels="AUTO",label_size = 14)
    plots_middle<-plot_grid(p_iph,p_iph_sex,ncol=2,labels=c("C","D"),label_size = 14)
    plots_bottom<-plot_grid(p_calc,p_calc_sex,ncol=2,labels=c("E","F"),label_size = 14)
    plots_bottom2<-plot_grid(p_coll,p_coll_sex,ncol=2,labels=c("G","H"),label_size = 14)
    plots_bottom3<-plot_grid(p_neo,p_neo_sex,ncol=2,labels=c("I","J"),label_size = 14)
    plots_all5<-plot_grid(plots_top,plots_middle,plots_bottom,plots_bottom2,plots_bottom3,ncol=1,labels=c("","","","",""),label_size = 14)
    
    tiff(paste("C://Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq/",gene_int,"/Plot_5_Plaque_characteristics_2_",gene_int,".tiff",sep=""),width=550,height=1250)
    print(plots_all5)
    dev.off()
  }
  if(!(gene_int%in%row.names(seuset@assays$RNA))){
    print(paste("Gene ",gene_int," is not present in the AE scRNAseq dataset.",sep=""))
  }else{
    
    
    
    zoom_out_seuset<-DimPlot(seuset,label=TRUE,label.size = 4,group.by="clusters_ct")+NoLegend()+theme(axis.title = element_text(size=9),
                                                                                                       axis.text = element_text(size=9),
                                                                                                       plot.title = element_blank())
    zoom_out_seuset_emt<-DimPlot(seuset_emt,label=TRUE,label.size = 4,reduction="umap")+NoLegend()+theme(axis.title = element_text(size=9),
                                                                                                         axis.text = element_text(size=9),
                                                                                                         plot.title = element_blank())
    Featureplot_seuset<-FeaturePlot(seuset,features = gene_int)+theme(axis.title = element_text(size=9),
                                                                      axis.text = element_text(size=9))
    Featureplot_seuset_emt<-FeaturePlot(seuset_emt,features = gene_int,reduction="umap")+theme(axis.title = element_text(size=9),
                                                                                               axis.text = element_text(size=9))
    
    VlnPlot_seuset<-VlnPlot(seuset,features = gene_int,group.by = "clusters_ct")+NoLegend()+labs(x="Slenders et al.",y=paste("Expression Level",gene_int))+theme(axis.text.x = element_blank(),
                                                                                                                                                                 axis.title.x = element_blank(),
                                                                                                                                                                 axis.ticks.x=element_blank())
    VlnPlot_seuset_sex<-VlnPlot(seuset,features = gene_int,split.by="Sex",group.by = "clusters_ct")+
      labs(x="Slenders et al.",y=paste("Expression Level",gene_int))+theme(axis.text.x = element_text(angle = 20,hjust=0.5,vjust=0.5),
                                                                           plot.title = element_blank(),
                                                                           legend.position="top",legend.justification = "center")
    
    VlnPlot_seuset_emt<-VlnPlot(seuset_emt,features = gene_int)+NoLegend()+labs(x="Slenders et al.",y=paste("Expression Level",gene_int))+
      theme(axis.title.y=element_text(size=9),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x=element_blank())
    VlnPlot_seuset_emt_sex<-VlnPlot(seuset_emt,features = gene_int,split.by="Sex")+
      labs(x="Slenders et al.",y=paste("Expression Level",gene_int))+theme(axis.title.y=element_text(size=9),
                                                                           axis.text.x = element_text(angle = 20,hjust=0.5,vjust=0.5),
                                                                           plot.title = element_blank(),
                                                                           legend.position="top",legend.justification = "center")
    
    
    top_plots<-plot_grid(zoom_out_seuset,Featureplot_seuset,ncol=2,labels="AUTO")
    bottom_plots<-plot_grid(VlnPlot_seuset,VlnPlot_seuset_sex,ncol=1,labels=c("C","D"),rel_heights = c(1,1.5))
    all_plots_seuset<-plot_grid(top_plots,bottom_plots,ncol=1)
    tiff(paste("C://Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq/",gene_int,"/Plot_6_scRNAseq_AE_",gene_int,".tiff",sep=""),width=1200,height=1250)
    print(all_plots_seuset)
    dev.off()
    
    top_plots<-plot_grid(zoom_out_seuset_emt,Featureplot_seuset_emt,ncol=2,labels="AUTO")
    bottom_plots<-plot_grid(VlnPlot_seuset_emt,VlnPlot_seuset_emt_sex,ncol=1,labels=c("C","D"),rel_heights = c(1,1.5))
    all_plots_seuset_emt<-plot_grid(top_plots,bottom_plots,ncol=1)
    tiff(paste("C://Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq/",gene_int,"/Plot_7_scRNAseq_AE_SMC_EC_",gene_int,".tiff",sep=""),width=600,height=625)
    print(all_plots_seuset_emt)
    dev.off()
    
  }
  
  if(!(gene_int%in%row.names(wirka_pv@assays$RNA)) | !(gene_int%in%row.names(pan_data@assays$RNA))){
    print(paste("Gene ",gene_int," is not present in the Wirka or Pan scRNAseq dataset.",sep=""))
  }else{
    
    zoom_out_wirka<-DimPlot(wirka_pv,label=TRUE,label.size = 3,group.by="Author_Provided")+NoLegend()+theme(axis.title = element_text(size=9),
                                                                                                            axis.text = element_text(size=9),
                                                                                                            plot.title = element_blank())
    zoom_out_pan<-DimPlot(pan_data,label=TRUE,label.size = 3,group.by="Author_Provided")+NoLegend()+theme(axis.title = element_text(size=9),
                                                                                                          axis.text = element_text(size=9),
                                                                                                          plot.title = element_blank())
    
    Feature_pan<-FeaturePlot(pan_data,features = gene_int)+theme(axis.title = element_text(size=9),
                                                                 axis.text = element_text(size=9))
    Feature_wirka<-FeaturePlot(wirka_pv,features = gene_int)+theme(axis.title = element_text(size=9),
                                                                   axis.text = element_text(size=9))
    
    VlnPlot_wirka<-VlnPlot(wirka_pv,features = gene_int,group.by = "Author_Provided")+NoLegend()+labs(x="Wirka et al.",y=paste("Expression Level",gene_int))+theme(axis.title.x = element_blank())
    VlnPlot_pan<-VlnPlot(pan_data,features = gene_int,group.by = "Author_Provided")+NoLegend()+labs(x="Pan et al.",y=paste("Expression Level",gene_int))+theme(axis.title.x = element_blank())
    top_plots<-plot_grid(zoom_out_pan,Feature_pan,ncol=2,labels="AUTO")
    top_plots2<-plot_grid(zoom_out_wirka,Feature_wirka,ncol=2,labels=c("D","E"))
    all_plots_other<-plot_grid(top_plots,VlnPlot_pan,top_plots2,VlnPlot_wirka,ncol=1,rel_heights = c(1.5,1,1.5,1))
    tiff(paste("C://Users/ediez/Documents/AE_Infrastructure/Lookup_AE_bulkRNAseq/",gene_int,"/Plot_8_scRNAseq_Others_",gene_int,".tiff",sep=""),width=800,height=1250)
    print(all_plots_other)
    dev.off()
  }
  
}

