#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####  ### #### 
#### Chris Nowak 21 Jan 2017 - genetic associations with C10 and C12-carnitines in KORA/TwinsUK (fasting plasma)    ### #### 
#### and associations with HOMA-IR in n~46,000 in MAGIC ### #### #### #### #### #### #### #### #### #### #### ####  ### #### 

dir <- "/Users/nowakchr/Desktop/2nd project scripts for publication/"
setwd(dir)

library(gtx)
library(survival)
library(ggplot2)
library(gridExtra)
library(grid)
library(MendelianRandomization)

## Download from KORA website the full plasma association data from Shin et al. - include all SNPs with p<1e05
  
  ## get the identifier for C10 and C12 carnitine in KORA/TwinsUK
kora_id <- read.delim("metaboliteMap.txt")
c10_id <- as.character(kora_id$metabolonID[which(kora_id$metabolonDescription == "decanoylcarnitine")])
c12_id <- as.character(kora_id$metabolonID[which(kora_id$metabolonDescription == "laurylcarnitine")])
  
  ## upload GWAS results that match metabolonID
c10_kora <- read.table(paste(c10_id,".metabolites.out",sep=""), header=T)
c12_kora <- read.table(paste(c12_id,".metabolites.out",sep=""), header=T)

  ## align to carnitine-increasing allele
c10_kora$risk_allele <- ifelse(c10_kora$Effect < 0, as.character(c10_kora$Allele2), as.character(c10_kora$Allele1))
c10_kora$kora_effect <- abs(c10_kora$Effect) 
c10_kora$kora_stderr <- c10_kora$StdErr

c12_kora$risk_allele <- ifelse(c12_kora$Effect < 0, as.character(c12_kora$Allele2), as.character(c12_kora$Allele1))
c12_kora$kora_effect <- abs(c12_kora$Effect) 
c12_kora$kora_stderr <- c12_kora$StdErr

 ## upload MAGIC lnHOMAIR GWAS results
magic_ir <- read.delim("MAGIC_ln_HOMA.txt")

 ## upload global lipid consortium and GIANT bmi gwas results
ldl <- read.delim("jointGwasMc_LDL")
hdl <- read.delim("jointGwasMc_HDL")
tg <- read.delim("jointGwasMc_TG")
tc <- read.delim("jointGwasMc_TC")
BMI <- read.delim("SNP_gwas_mc_merge_nogc.tbl.uniq.txt")

 ## how many markers are shared between KORA/TwinsUK and MAGIC-HOMA-IR?
paste("For C10, out of",length(c10_kora$MarkerName), "SNPs in KORA,", length(which(as.character(c10_kora$MarkerName) %in% magic_ir$snp)),"are in MAGIC.",
      "For C12, out of",length(c12_kora$MarkerName), "SNPs in KORA,", length(which(as.character(c12_kora$MarkerName) %in% magic_ir$snp)),"are in MAGIC.")

  ####################################################################################################################################
  ########### [1] "For C10, out of 760 SNPs in KORA, 744 are in MAGIC. For C12, out of 419 SNPs in KORA, 405 are in MAGIC."
  ####################################################################################################################################
  
 ## merge KORA and MAGI
c10_merge <- merge(c10_kora, magic_ir, by.x = "MarkerName", by.y = "snp")
c12_merge <- merge(c12_kora, magic_ir, by.x = "MarkerName", by.y = "snp")

  ## align alleles (risk allele = carnitine-increasing allele)
c10_merge$ir_effect <- ifelse(c10_merge$risk_allele == c10_merge$effect_allele, c10_merge$effect, -c10_merge$effect)
c12_merge$ir_effect <- ifelse(c12_merge$risk_allele == c12_merge$effect_allele, c12_merge$effect, -c12_merge$effect)

  ## Make a function for genetic scatter plots
plot_scatter <- function(data, title , p_cutoff){
  metabo <- data
  if(p_cutoff == F){ 
    df <- data.frame(x = metabo$kora_effect, xmin = metabo$kora_effect - (1.96*metabo$kora_stderr), xmax = metabo$kora_effect + (1.96*metabo$kora_stderr),
                     y = metabo$ir_effect, ymin = metabo$ir_effect - (1.96*metabo$stderr), ymax = metabo$ir_effect + (1.96*metabo$stderr),
                     snp_name = metabo$MarkerName)
    p <- ggplot(data = df,aes(x = x,y = y)) + geom_point() + 
      geom_errorbar(aes(ymin = ymin, ymax = ymax)) + 
      geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +
      geom_smooth(method = 'lm', formula = y~x, se = F,col = "red", fullrange = T) +
      geom_hline(aes(yintercept = 0)) +
      geom_vline(aes(xintercept = 0)) +
      ggtitle(title) + 
      theme(axis.text = element_text(size=10), axis.title = element_text(size = 12)) +
      labs(y = "Genetic association with log-HOMA-IR in MAGIC",
           x = "Genetic association with plasma acylcarnitine in KORA/TwinsUK") 
    } else {
    metabo <- metabo[which(metabo$P.value < p_cutoff),]
    df <- data.frame(x = metabo$kora_effect, xmin = metabo$kora_effect - (1.96*metabo$kora_stderr), xmax = metabo$kora_effect + (1.96*metabo$kora_stderr),
                     y = metabo$ir_effect, ymin = metabo$ir_effect - (1.96*metabo$stderr), ymax = metabo$ir_effect + (1.96*metabo$stderr),
                     snp_name = metabo$MarkerName)
    p <- ggplot(data = df,aes(x = x, y = y)) + geom_point() + 
      geom_errorbar(aes(ymin = ymin, ymax = ymax)) + 
      geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +
      geom_smooth(method = 'lm', formula = y~x, se = F, col = "red", fullrange = T) +
      geom_hline(aes(yintercept = 0)) +
      geom_vline(aes(xintercept = 0)) +
      ggtitle(title) + 
      theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12)) +
      labs(y = "Genetic association with log-HOMA-IR in MAGIC",
           x = "Genetic association with plasma acylcarnitine in KORA/TwinsUK") 
  }
  return(p)
}

pdf("plot_1.pdf",w=12,h=12)
p1<-plot_scatter(c10_merge,title="C10",p_cutoff=F)
p2<-plot_scatter(c12_merge,title="C12",p_cutoff=F)
p3<-plot_scatter(c10_merge,title="C10, p<5e-08",p_cutoff=5e-08)
p4<-plot_scatter(c12_merge,title="C12, p<5e-08",p_cutoff=5e-08)
grid.arrange(p1,p2,p3,p4, ncol = 2,nrow=2)
dev.off()

 ## Is a GRS weighted by the effect size for natural-log C10 or C12 associated with natural-log-HOMA-IR?
 ## scale in KORA log with base 10. To convert to natural log = log10(x)*2.303
wi <- c10_merge$kora_effect*2.303
betai <- c10_merge$ir_effect
si <- c10_merge$stderr
grs <- grs.summary(wi, betai, si, 46000)
paste("C10-GRS, p<1e-05", "y=",grs$ahat, "[", grs$ahat - (1.96*grs$aSE), grs$ahat + (1.96*grs$aSE),"]", "p=",grs$pval)

  ####################################################################################################################################  
  ########### [1] "C10-GRS, p<1e-05 y= -0.0215869154074706 [ -0.0267150555267339 -0.0164587752882073 ] p= 1.57570033439956e-16"
  ####################################################################################################################################

wi <- c12_merge$kora_effect*2.303
betai <- c12_merge$ir_effect
si <- c12_merge$stderr
grs <- grs.summary(wi, betai, si, 46000)
paste("C12-GRS, p<1e-05", "y=",grs$ahat, "[", grs$ahat - (1.96*grs$aSE), grs$ahat + (1.96*grs$aSE),"]", "p=",grs$pval)

  ####################################################################################################################################
  ########### [1] "C12-GRS, p<1e-05 y= -0.0216156217545099 [ -0.030069423160966 -0.0131618203480537 ] p= 5.39943390357556e-07"
  ####################################################################################################################################

  ## GWAS sign only 
c10_gwas <- c10_merge[which(c10_merge$P.value < 5e-08),]
c12_gwas <- c12_merge[which(c12_merge$P.value < 5e-08),]

wi <- c10_gwas$kora_effect*2.303
betai <- c10_gwas$ir_effect
si <- c10_gwas$stderr
grs <- grs.summary(wi, betai, si, 46000)
paste("C10-GRS, p<5e08 snps,", "y=",grs$ahat, "[", grs$ahat - (1.96*grs$aSE), grs$ahat + (1.96*grs$aSE),"]", "p=",grs$pval)

  ####################################################################################################################################
  ########### [1] "C10-GRS, p<5e08 snps, y= -0.039520384393036 [ -0.0458914807688753 -0.0331492880171966 ] p= 5.2000261038348e-34"
  ####################################################################################################################################


wi <- c12_gwas$kora_effect*2.303
betai <- c12_gwas$ir_effect
si <- c12_gwas$stderr
grs <- grs.summary(wi, betai, si, 46000)
paste("C12-GRS, p<5e08 snps", "y=",grs$ahat, "[", grs$ahat - (1.96*grs$aSE), grs$ahat + (1.96*grs$aSE),"]", "p=",grs$pval)

  ####################################################################################################################################
  ########### [1] "C12-GRS, p<5e08 snps y= -0.0311393088063747 [ -0.0897563837236399 0.0274777661108905 ] p= 0.297775240449898"
  ####################################################################################################################################

  ## MR analysis for GWAS-significant carnitine SNPs
mr_c10 <- mr_allmethods(mr_input(bx = (c10_gwas$kora_effect*2.303), bxse = (c10_gwas$kora_stderr*2.303),
                                 by = c10_gwas$ir_effect, byse = c10_gwas$stderr))$Values
mr_c12 <- mr_allmethods(mr_input(bx = (c12_gwas$kora_effect*2.303), bxse = (c12_gwas$kora_stderr*2.303),
                                 by = c12_gwas$ir_effect, byse = c12_gwas$stderr))$Values
mr_c10
##                  Method     Estimate    Std Error      95% CI                     P-value
##1              Simple median -0.033337374 0.0045982051 -0.042349690 -2.432506e-02 4.165130e-13
##2            Weighted median -0.032199720 0.0045708276 -0.041158378 -2.324106e-02 1.859751e-12
##3  Penalized weighted median -0.032199720 0.0045708276 -0.041158378 -2.324106e-02 1.859751e-12
##4                        IVW -0.039520384 0.0032505594 -0.045891364 -3.314941e-02 2.950224e-34
##5              Penalized IVW -0.039520384 0.0032505594 -0.045891364 -3.314941e-02 2.950224e-34
##6                 Robust IVW -0.039039988 0.0026422630 -0.044218728 -3.386125e-02 6.867725e-38
##7       Penalized robust IVW -0.039039988 0.0026422630 -0.044218728 -3.386125e-02 6.867725e-38
##8                   MR-Egger -0.022224402 0.0122966736 -0.046325439  1.876635e-03 7.070761e-02
##9                (intercept) -0.001734645 0.0009553723 -0.003607140  1.378501e-04 6.942031e-02
##10        Penalized MR-Egger -0.021170157 0.0123242571 -0.045325257  2.984943e-03 8.583979e-02
##11               (intercept) -0.001851175 0.0009517378 -0.003716547  1.419688e-05 5.176930e-02
##12           Robust MR-Egger  0.075894937 0.0159375171  0.044657977  1.071319e-01 1.916551e-06
##13               (intercept) -0.012872985 0.0008213550 -0.014482811 -1.126316e-02 0.000000e+00
##14 Penalized robust MR-Egger  0.075894656 0.0159375851  0.044657563  1.071317e-01 1.916911e-06
##15               (intercept) -0.012872951 0.0008213709 -0.014482808 -1.126309e-02 0.000000e+00

mr_c12
##                Method     Estimate   Std Error       95% CI                    P-value
##1              Simple median -0.044869011 0.037573832 -0.1185123683  0.028774345 2.324170e-01
##2            Weighted median -0.044931367 0.035349561 -0.1142152328  0.024352499 2.037078e-01
##3  Penalized weighted median -0.044931367 0.035349561 -0.1142152328  0.024352499 2.037078e-01
##4                        IVW -0.031139309 0.029906671 -0.0897553066  0.027476689 2.263881e-01
##5              Penalized IVW -0.031139309 0.029906671 -0.0897553066  0.027476689 2.263881e-01
##6                 Robust IVW -0.043725869 0.017931760 -0.0788714737 -0.008580265 5.266362e-05
##7       Penalized robust IVW -0.043725869 0.017931760 -0.0788714737 -0.008580265 5.266362e-05
##8                   MR-Egger -0.266346891 0.807912550 -1.8498263921  1.317132610 7.416471e-01
##9                (intercept)  0.016614217 0.047035641 -0.0755739456  0.108802380 7.239189e-01
##10        Penalized MR-Egger -0.266346891 0.807912550 -1.8498263921  1.317132610 7.416471e-01
##11               (intercept)  0.016614217 0.047035641 -0.0755739456  0.108802380 7.239189e-01
##12           Robust MR-Egger -0.107000380 0.151199046 -0.4033450646  0.189344304 4.791447e-01
##13               (intercept)  0.004481536 0.002026671  0.0005093332  0.008453739 2.701651e-02
##14 Penalized robust MR-Egger -0.107000380 0.151199046 -0.4033450646  0.189344304 4.791447e-01
##15               (intercept)  0.004481536 0.002026671  0.0005093332  0.008453739 2.701651e-02


  ## Are any SNPs GWAS-significant associated with C10 also associated with plasma lipids or BMI?
length(which(ldl$P.value[which(ldl$rsid %in%c10_gwas$MarkerName)] < 5e-08))
  ####### [1] 0
length(which(hdl$P.value[which(hdl$rsid %in%c10_gwas$MarkerName)] < 5e-08))
  ####### [1] 0
length(which(tg$P.value[which(tg$rsid %in%c10_gwas$MarkerName)] < 5e-08))
  ####### [1] 0
length(which(tc$P.value[which(tc$rsid %in%c10_gwas$MarkerName)] < 5e-08))
  ####### [1] 0
length(which(BMI$p[which(BMI$SNP %in%c10_gwas$MarkerName)] < 5e-08))
  ####### [1] 0
