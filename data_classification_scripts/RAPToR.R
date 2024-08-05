#! /usr/bin/env/ Rscript
# Run RAPToR on Worm ASE data
# By Francisco Valencia Avalos

## The objective of this script was to determine biological age of the worms used
## for RNAseq in another project in the lab. The main package used to do this was
## RAPToR (https://github.com/LBMC/RAPToR)

#load packages
library(RAPToR)
library(DESeq2)
library(wormRef)

## DESeq2 was used since the data has already been analyzed and its in an R object
## created with this package
## wormRef is part of the RAPToR package and contains the reference used to calculate ages

#load the R object from DESeq that contains the analyzed RNA seq data. These files are saved automatically by script ase_de_annotategenes_deseq2_fromemaseout.R
load(file.choose())

## for this script the R object loaded was dds_grp located in the lab's dropbox

#vst transform rna seq data
vst_dds_grp<-vst(dds.all,blind=T)
data_dds_grp<-assay(vst_dds_grp)

## the package recommends to normalize the expression data although it is not necesessary just for staging/
## the transformation of the data was done here with the vst command from DESeq2.
## Other methods of normalization were tried (limma, normalizeBetweenArrays) and the final
## results were similar or identical.

#load reference
ref<-prepare_refdata("Cel_YA_2","wormRef",600)

## in this case because of age of worms, Cel_YA_2 was used
## Age estimation was performed with different references (Cel_Ya_1, Cel_Larv_YA) and results were
## the same or low quality compared to Cel_YA_2.

#age estimation
ae_ptx<-ae(data_dds_grp,ref)

#convert to data frames and merge
sample_names<-as.data.frame(colData(dds.all))
age_estimates<-as.data.frame(ae_ptx$age.estimates)
age_estimate_final<-merge(sample_names,age_estimates,by=0)

## this step was to add sample information to the calculations for plotting purposes

#plotting with ggplot
library(ggplot2)

p<-ggplot(age_estimate_final,aes(x=factor(Strain),y=age.estimate,color=Strain))+
  geom_point(na.rm=TRUE, position=position_dodge2(0.9))+
  geom_errorbar(aes(ymin=lb,ymax=ub),position=position_dodge2(0.9))+
  scale_x_discrete(guide=guide_axis(angle=45))+
  ylab(expression(atop("Age estimation from gene expression profile",paste("(hrs after hatching)"))))+
  xlab("Strains")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), # angle = 45, vjust = 0.55, hjust = 0
        axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14), legend.position = "right",
        legend.text = element_text(size = 13), legend.title = element_text(size = 14))+
  facet_wrap(~Generation)

pdf("plot.pdf",10,5)
print(p)
dev.off()


# with ggplot to change labels
facet_labeller <- c(
  "CTR" = "control", "PAR1" = "par-1", "POS1" = "pos-1"
)
p + facet_wrap(Treatment ~ ., labeller = as_labeller(facet_labeller))

### OTHERS (not part of main script):

# check for correlation value using different references
plot_refs("wormRef")
ref_YA_1<-prepare_refdata("Cel_YA_1","wormRef",600)
ref_larv_YA<-prepare_refdata("Cel_larv_YA","wormRef",600)

ae_ptx_YA1<-ae(data_dds_grp,ref_YA_1)
ae_ptx_larvYA<-ae(data_dds_grp,ref_larv_YA) # this one has higher correlation scores, many samples con edge of calculation

plot(ae_ptx_larvYA)
plot(ae_ptx_YA1)
plot(ae_ptx)

#Comparing normalization methods

library(limma)
quantile_normalized<-data_dds_grp
quantile_normalized<-limma::normalizeBetweenArrays(quantile_normalized,method="quantile")
quantile_normalized<-log1p(quantile_normalized)
ae_qnorm<-ae(quantile_normalized,ref)

quantile_normalized_novst<-"counts"(dds_grp) #no transformation vst
quantile_normalized_novst<-limma::normalizeBetweenArrays(quantile_normalized_novst,method="quantile")
quantile_normalized_novst<-log1p(quantile_normalized_novst)
ae_qnormnovst<-ae(quantile_normalized_novst,ref)


View(ae_ptx$age.estimates)
View(ae_qnorm$age.estimates) #no changes to calculated age, minimal changes to error
View(ae_qnormnovst$age.estimates) # minimal changes to age calculated and error
