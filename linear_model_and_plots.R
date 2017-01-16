######################################################################################################################### ############ ############ #############
## Association between insulin sensitivity clamp M/I and blood metabolites measured at 0, 30, 120min during a 75g oral glucose challenge ########## #############
## 1) clamp = linear combination of metabolite values at any time point, adjusted for age, (all male), sample quality ### ############ ############ #############
##            ANOVA-type likelihood ratio test of age-quality baseline model vs. added metabolites ->  take bonferroni-sign metabolites forward ### #############
## 2) clamp = metab(30 - 0 min) + metab(120 - 30 min) + age + quality -> i.e. association between change over time and insulin sensitiviy - bonferroni adjusted #
## 3) plots = barplot of quality-age-adjusted metabolite levels at 0/30/120 min split into clamp M/I quartiles ########## ############ ############ #############
######################################################################################################################### ############ ############ #############

library(foreign)
library(car)
library(psych)
library(fmsb)
library(reshape2)

dir <- ("/ ... ")
setwd(dir)

ulsam_LM <- data.frame(read.table("ULSAM_MAIN_forLM", header = TRUE, quote = "\""))     # complete ULSAM dataset
sample_char <- data.frame(read.table("~/ ... /THAWED_UNTHAWED_for_R.txt", header = T))  # quality dummies as per lab logbooks

ulsam_LM2 <- merge(ulsam_LM, sample_char, by = "pat")
ulsam_LM2$sample_EDTA <- ifelse(ulsam_LM2$sample_type == "EDTA", 1, 0)

  ## split into DM and nonDM according to diagnosis of diabetes, FG >= 7.0 and/or 2h OGTT glucose >=11.1
ULSAM_nonDM <- subset(ulsam_LM2, z378 != 1 & z319.x < 7 & z323 < 11.1)  # non-diabetics
ULSAM_DM <- subset(ulsam_LM2, z378 == 1 & z319.x >= 7 & z323 >= 11.1)   # diabetics

# write.table(ULSAM_nonDM, "ULSAM_nonDM.txt")

  ## rename/code variables
ULSAM_nonDM$Clamp <- ULSAM_nonDM$z620 
ULSAM_nonDM$Age <- ULSAM_nonDM$age70 
ULSAM_nonDM$Storage_time <- ULSAM_nonDM$dat70_numeric - min(ULSAM_nonDM$dat70_numeric)  # in days, earliest = 0
ULSAM_nonDM$Thawed <- ULSAM_nonDM$thawed        # 1 = tinad/har tinats, i.e. suspicion of having been out of -80°C before
ULSAM_nonDM$Sample_type <- as.numeric(ULSAM_nonDM$sample_EDTA)
ULSAM_nonDM$Quality <- ULSAM_nonDM$low_quality  # 1 = haemolysed/other comments related to sample quality as per lab logbook
ULSAM_nonDM$waist <- ULSAM_nonDM$z020
ULSAM_nonDM$bmi <- ULSAM_nonDM$z290
ULSAM_nonDM$sbp <- ULSAM_nonDM$z013
ULSAM_nonDM$dbp <- ULSAM_nonDM$z014
ULSAM_nonDM$crp <- ULSAM_nonDM$z802
ULSAM_nonDM$fasting_ins <- ULSAM_nonDM$z606.x
ULSAM_nonDM$hdlc <- ULSAM_nonDM$z302
ULSAM_nonDM$ldlc <- ULSAM_nonDM$z324
ULSAM_nonDM$tg <- ULSAM_nonDM$z971
ULSAM_nonDM$chol <- ULSAM_nonDM$z972
ULSAM_nonDM$M0 <- ULSAM_nonDM[,grep("^X0_", colnames(ULSAM_nonDM))]     # 0min metabolite levels
ULSAM_nonDM$M30 <- ULSAM_nonDM[,grep("^X30_", colnames(ULSAM_nonDM))]   # 30min
ULSAM_nonDM$M120 <- ULSAM_nonDM[,grep("^X120_", colnames(ULSAM_nonDM))] # 120min

#########################################################################
#### CHECKS & DECRIPTIVES ###############################################
  
  ## check correlations clampM/I x main covariates  
attach(ULSAM_nonDM)
cor.test(Clamp, Thawed, method = "pearson")
cor.test(Clamp, Storage_time, method = "pearson")
cor.test(Clamp, Quality, method = "pearson")
cor.test(Clamp, Sample_type, method = "pearson")
cor.test(Clamp, Age, method = "pearson")

  ## normal distribution for clamp M/I -> histograms/qqplots by diabetes status
par(mfrow = c(2,2))
hist(ulsam_LM2$z620, prob = T, breaks = 30, ylim = c(0,0.4), main = "ClampMI, all", col = "grey")
hist(ulsam_LM2$z620[ulsam_LM2$z378 == 0], prob = T, breaks = 30, add = F, ylim = c(0,0.4), col = rgb(0,0,1,1), main = "ClampMI, NonDM (blue) vs. DM (red)")
hist(ulsam_LM2$z620[ulsam_LM2$z378 == 1], prob = T, breaks = 20, add = T, col = rgb(1,0,0,1/2), main = "")
qqnorm((ULSAM_nonDM$Clamp), main = "nonDM")
qqline(ULSAM_nonDM$Clamp)

  ## comparison clamp M/I transformations, raw/log/sqrt/power^2, regarding normal distribution in qqplots
par(mfrow = c(2,2)) 
qqnorm((ULSAM_nonDM$Clamp), main = "NonDM ClampMI")
qqline(ULSAM_nonDM$Clamp)
qqnorm(log(ULSAM_nonDM$Clamp), main = "nonDM natural log ClampMI")
qqline(log(ULSAM_nonDM$Clamp))
qqnorm(sqrt(ULSAM_nonDM$Clamp), main = "nonDM Sqrt ClampMI")
qqline(sqrt(ULSAM_nonDM$Clamp))
qqnorm(ULSAM_nonDM$Clamp^2, main = "nonDM Squared ClampMI")
qqline(ULSAM_nonDM$Clamp^2)
  ## ==================== raw clamp M/I reasonably normally distributed ==================== ##

  ## diagnotics baseline linear model: clamp on age + quality
LM_baseline <- lm(Clamp ~ Age + Storage_time + Thawed + Sample_type + Quality, data = ULSAM_nonDM) 
summary(LM_baseline)        # 8 excluded due to clamp-NA
par(mfrow = c(3,1))         # plots residuals
hist(LM_baseline$res, breaks = 20, main = "nonDM - rawClampMI - Residuals", cex.main = 2, col = "grey")
plot(LM_baseline, which = c(1,2))

  ## outlier/influencer observations? -> leverage, discrepancy plots
par(mfrow = c(2,1))
plot(LM_baseline, which =  4, main = "nonDM - Cook's D (rawClamp)", cex.main = 2)   # Cook's distance
outlierTest(LM_baseline)    # studentised residuals: bonferroni-adj outlier test + plot
qqPlot(LM_baseline, simulate = T, labels = row.names(data.frame(Clamp)), main = "Studentised Residuals")    

  ## interactive influence plot (left click to see ID of point): stud res + Cook's + hat values (leverage)
par(mfrow = c(1,1))
plot(hatvalues(LM_baseline), rstudent(LM_baseline), type = 'n', ylim = c(-3,6),
     xlab = "hat values (leverage/extreme)", ylab = "student. residuals (aberration from regression line)", main = "Influence plot (diameter = Cook's Dist.")
cook <- sqrt(cooks.distance(LM_baseline))
points(hatvalues(LM_baseline), rstudent(LM_baseline), cex = 10*cook/max(cook))
abline(v = 3/25,lty = 2)
abline(h = c(-2,0,2), lty = 2)
Cllm <- data.frame(cbind(ULSAM_nonDM$Clamp, ULSAM_nonDM$Age, ULSAM_nonDM$Storage_time, ULSAM_nonDM$Thawed, ULSAM_nonDM$Sample_type, ULSAM_nonDM$Quality)) 
Cll <- Cllm[!is.na(Cllm[,1]), ]     # remove 8 clamp NA
identify(hatvalues(LM_baseline), rstudent(LM_baseline), row.names(Cll))
  ## ==================== DECISION: no gross violations, take all 470 with complete clamp forward ==================== ## 

  ## descriptives
NA_clamp <- which(is.na(ULSAM_nonDM$Clamp))
describe(ULSAM_nonDM$Age[-NA_clamp])        
ULSAM_nonDM$box <- 0
ULSAM_nonDM$box[which(is.na(ULSAM_nonDM$Clamp))] <- 1
boxplot(ULSAM_nonDM$Age ~ ULSAM_nonDM$box)  # compare the 8 wo. clamp to the 470 with clamp
  ## etc. ...

#####################################################################################
#### Linear regression analysis: 1) LRT v. baseline model 2) change over time #######

  ## linear model 1, sd-unti metabolites
METABO_coeff <- sapply(1:192, function(x) 
  coef(lm(Clamp ~ Age + Storage_time + Thawed + Sample_type + Quality + scale(M0[,x]) + scale(M30[,x]) + scale(M120[,x]), data = ULSAM_nonDM)))
METABO_fstat <- sapply(1:192, function(x) 
    summary(lm(Clamp ~ Age + Storage_time + Thawed + Sample_type + Quality + scale(M0[,x]) + scale(M30[,x]) + scale(M120[,x]), data = ULSAM_nonDM))$fstatistic[1])
METABO_fstat_p <- sapply(1:192, function(x) pf(
    summary(lm(Clamp ~ Age + Storage_time + Thawed + Sample_type + Quality + scale(M0[,x]) + scale(M30[,x]) + scale(M120[,x]), data = ULSAM_nonDM))$fstatistic[1],
    summary(lm(Clamp ~ Age + Storage_time + Thawed + Sample_type + Quality + scale(M0[,x]) + scale(M30[,x]) + scale(M120[,x]), data = ULSAM_nonDM))$fstatistic[2],
    summary(lm(Clamp ~ Age + Storage_time + Thawed + Sample_type + Quality + scale(M0[,x]) + scale(M30[,x]) + scale(M120[,x]), data = ULSAM_nonDM))$fstatistic[3],
    lower.tail=F)) # model F-statistic 
  
  ## likelihood ratio test against age-quality-only model and bonferroni-adj. 
Baseline <- lm(Clamp ~ Age + Storage_time + Thawed + Sample_type + Quality, data = ULSAM_nonDM)
LRTest <- sapply(1:192, function(x)
  anova(lm(Clamp ~ Age + Storage_time + Thawed + Sample_type + Quality + scale(M0[,x]) + scale(M30[,x]) + scale(M120[,x]), data = ULSAM_nonDM), Baseline))
LRTest_p <- sapply(1:192, function(x) 
  LRTest[,x]$`Pr(>F)`[2]) 
Metanames <- colnames(M0)
t1 <- data.frame(cbind(as.character(Metanames), as.numeric(METABO_fstat) , as.numeric(METABO_fstat_p), as.numeric(LRTest_p), as.numeric(p.adjust(LRTest_p, method = "bonferroni"))))
t1$Sign <- ifelse(as.numeric(as.character(t1[,5])) <= 0.05, "sign.", "n.s.")
colnames(t1)<-c("Metabolite", "F_stat", "F_stat_p_value", "LR_Test_p_value", "Bonf_adj_LR_Test_p_value", "Bonf_sign")
# View(t1)
  
  ## produce output table and match MxxxTxxx with actual annotated metabolite names 
M_names <- colnames(ULSAM_nonDM[, grep("^X0_", colnames(ULSAM_nonDM))])
M_names <- gsub("^X0_", "", M_names)
t <- data.frame(cbind(M_names, t1))
matchlist <- data.frame(read.delim("~/ ... /metabo_id_match_list.txt"))
matchid <- (match(M_names,matchlist[, 3]))
matchlist$matchid = 1:192
matchlist <- matchlist[, c(1,3,4,7)]
t$matchid <- matchid
t2 <- merge(matchlist, t, by = "matchid")
output_table <- with(t2, t2[order(t2[,10]), ])

  ## check for multicollinearity - Variable Inflation Factor VIF 
VIF_M0 <- sapply(1:192, function(x) 
  VIF(lm(scale(M0[,x]) ~ Age + Storage_time + Thawed + Sample_type + Quality  + scale(M30[,x]) + scale(M120[,x]), data = ULSAM_nonDM)))
VIF_M30 <- sapply(1:192, function(x) 
  VIF(lm(scale(M30[,x]) ~ Age + Storage_time + Thawed + Sample_type + Quality  + scale(M0[,x]) + scale(M120[,x]), data = ULSAM_nonDM)))
VIF_M120 <- sapply(1:192, function(x) 
  VIF(lm(scale(M120[,x]) ~ Age + Storage_time + Thawed + Sample_type + Quality  + scale(M0[,x]) + scale(M30[,x]), data = ULSAM_nonDM)))
VIF <- data.frame(cbind(VIF_M0, VIF_M30, VIF_M120))
VIF$p0 <- VIF$VIF_M0 > 5 # or sqrt(VIF) > 2
VIF$p30 <- VIF$VIF_M30 > 5 # sqrt(VIF) shows how much larger the SE is
VIF$p120 <- VIF$VIF_M120 > 5  # compared to 1=no collinearity
#  View(VIF)
sum(sqrt(VIF$VIF_M0) > 2) # convention: check if sqrt(VIF)>2
sum(sqrt(VIF$VIF_M30) > 2)
sum(sqrt(VIF$VIF_M120) > 2)
  
  ## for the associated meabolites: association clamp M/I and change 0-30 min + 30-120 min 
# u_temp1 <- read.table("ULSAM_nonDM.txt", header = T) # Masterfile with non-diabetics

u_temp1 <- ULSAM_nonDM.txt
u_temp2 <- read.csv("Results_Overall_LM.csv") # results final linear model

  ## rename news variables
u_temp1$BMI <- u_temp1$z290
u_temp1$FI <- u_temp1$z613
u_temp1$FG <- u_temp1$z319.x

 ## select only bonferroni-significant metabolite from linear model 1
M_id <- which(u_temp2$METABO_anova_p < (0.05/192)) 
u_temp1$M0_lim <- u_temp1$M0[, M_id] 
u_temp1$M30_lim <- u_temp1$M30[, M_id]
u_temp1$M120_lim <- u_temp1$M120[, M_id]

 ## tie together
attach(u_temp1)
MAIN <- cbind(pat, Clamp, Age, Storage_time, Thawed, Sample_type, Quality, FI, FG,  BMI, M0_lim, M30_lim, M120_lim)
clampNA_ids <- which(is.na(u_temp1$Clamp) | is.na(u_temp1$FI)) 
MAIN1 <- MAIN[-clampNA_ids, ] # retain only non-NA
Mnames0 <- names(MAIN1[,grep("^X0", names(MAIN1))])
Mnames30 <- names(MAIN1[,grep("^X30", names(MAIN1))])
Mnames120 <- names(MAIN1[,grep("^X120", names(MAIN1))])
# write.table(MAIN1,"MAIN1.csv")

  ## metabolite levels adjusted for age-quality
M0_res <- sapply(grep("^X0", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality, data = MAIN1)))
M30_res <- sapply(grep("^X30", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality, data = MAIN1)))
M120_res <- sapply(grep("^X120", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality, data = MAIN1)))
M0_res <- data.frame(M0_res)
M30_res <- data.frame(M30_res)
M120_res <- data.frame(M120_res)
MERGED1 <- cbind(MAIN1$Clamp, M0_res, M30_res, M120_res)
colnames(MERGED1) <- c("Clamp", Mnames0, Mnames30, Mnames120)
MERGED1 <- data.frame(MERGED1)
MERGED2 <- MERGED1[order(MERGED1$Clamp),]
# write.table(MERGED2,"RESIDUALS_adj.csv")
  
  ## variables for metabolite change 0-30 min and 30-120 min
Residuals <- MERGED2
DIFF0_30 <- (sapply (2:36, function (x)  Residuals[,x+35] - Residuals[,x])) ## get the difference between the time points
DIFF30_120 <- (sapply (2:36, function (x)  Residuals[,x+35+35] - Residuals[,x+35]))
DIFF0_120 <- (sapply (2:36, function (x)  Residuals[,x+35+35] - Residuals[,x]))
colnames(DIFF0_30) <- names(Residuals)[2:36] # put names on metabolites
colnames(DIFF30_120) <- names(Residuals)[2:36]
colnames(DIFF0_120) <- names(Residuals)[2:36]
DIFF0_30 <- data.frame(DIFF0_30)
DIFF30_120 <- data.frame(DIFF30_120)
DIFF0_120 <- data.frame(DIFF0_120)
DIFF0_30 <- cbind(Residuals$Clamp, DIFF0_30) # merge differences with clamp variable
DIFF30_120 <- cbind(Residuals$Clamp, DIFF30_120)
DIFF0_120 <- cbind(Residuals$Clamp, DIFF0_120)

  ## linear models for clamp M/I = metabolite change
Fstat <- sapply (2:36, function (x)
  summary(lm(DIFF0_30[,1] ~ DIFF0_30[,x] + DIFF30_120[,x]))$fstatistic)
F_p_value <- sapply (1:35, function(x) pf(Fstat[1,x], Fstat[2,x], Fstat[3,x], lower = FALSE))
b_030 <- sapply (2:36, function (x)
  summary(lm(DIFF0_30[,1] ~ DIFF0_30[,x] + DIFF30_120[,x]))$coefficients[2,1])
se_030 <- sapply (2:36, function (x)
  summary(lm(DIFF0_30[,1] ~ DIFF0_30[,x] + DIFF30_120[,x]))$coefficients[2,2])
p_030 <- sapply (2:36, function (x)
  summary(lm(DIFF0_30[,1] ~ DIFF0_30[,x] + DIFF30_120[,x]))$coefficients[2,4])
b_30120 <- sapply (2:36, function (x)
  summary(lm(DIFF0_30[,1] ~ DIFF0_30[,x] + DIFF30_120[,x]))$coefficients[3,1])
se_30120 <- sapply (2:36, function (x)
  summary(lm(DIFF0_30[,1] ~ DIFF0_30[,x] + DIFF30_120[,x]))$coefficients[3,2])
p_30120 <- sapply (2:36, function (x)
  summary(lm(DIFF0_30[,1] ~ DIFF0_30[,x] + DIFF30_120[,x]))$coefficients[3,4])
Mnames <- names(DIFF0_30[2:36])
clamp_delta_results <- cbind(Mnames, F_p_value, b_030, se_030, p_030, b_30120, se_30120, p_30120)
# View(data.frame(clamp_delta_results))
# write.table(clamp_delta_results, "Clamp_DELTA_RESULTS_FI_ADJUSTED.csv")

  ## produce output table with matched metabolite names
M_names <- gsub("^X0_", "", Mnames)
clamp_delta_results <- data.frame(clamp_delta_results)
clamp_delta_results$M_names <- M_names
t <- clamp_delta_results
matchlist <- data.frame(read.delim("~/ ... /metabo_id_match_list.txt"))
matchid <- (match(M_names,matchlist[, 3]))
matchlist$matchid = 1:192
matchlist <- matchlist[, c(1,3,4,7)]
t$matchid <- matchid
t2 <- merge(matchlist, t, by = "matchid")
# View(with(t2,t2[order(t2[,10]),]))
# write.table(t2, "results_difference.txt")

##############################################################################################
#### PLOTS - sample script for additional adj. for fasting insulin, fasting glucose etc. #####

  ## barplot function - clamp M/I quartiles at 0/30/120 min
METABOplot <- function(x, y, z, v, vv, s, t, r, main = s, ylab = t, xlab = r, length = 0.05){
  xx <- tapply(y, list(z, x), mean)     # bar height = Mean
  yy <- tapply(y, list(z, x), sd)       # sd
  zz <- tapply(y, list(z, x), length)   # nr. of replicates
  er <- yy / sqrt(zz)     # sd conversion to se
  w <- length(levels(z))  # no. of bar colours
  barx <- barplot(xx, col = c(1:length(x)), beside = T, main = main, ylab = ylab, cex.main = 0.9, xlab = xlab, space = c(0.4, 2),
                  ylim = c(v, vv), xpd = FALSE, cex.lab = 0.8, cex.axis = 0.8, font.lab = 4, cex.names = 0.8)
  box()                                 # box around plot
  arrows(barx, xx + 1.96*er, barx, xx, angle = 90, code = 1, length = length)   # error bars = MEAN +- 95%CI
  arrows(barx, xx - 1.96*er, barx, xx, angle = 90, code = 1, length = length)
}

# MAIN1 <- read.table("MAIN1.csv")
Mnames0 <- names(MAIN1[,grep("^X0", names(MAIN1))])
Mnames30 <- names(MAIN1[,grep("^X30", names(MAIN1))])
Mnames120 <- names(MAIN1[,grep("^X120", names(MAIN1))])

    ## FI adjusted residuals
M0_res_fi <- sapply(grep("^X0", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality + FI, data = MAIN1)))
M30_res_fi <- sapply(grep("^X30", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality+ FI, data = MAIN1)))
M120_res_fi <- sapply(grep("^X120", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality+ FI, data = MAIN1)))
M0_res_fi <- data.frame(M0_res_fi)
M30_res_fi <- data.frame(M30_res_fi)
M120_res_fi <- data.frame(M120_res_fi)
MERGED1 <- cbind(MAIN1$Clamp, M0_res_fi, M30_res_fi, M120_res_fi)
colnames(MERGED1) <- c("Clamp", Mnames0, Mnames30, Mnames120)
MERGED1 <- data.frame(MERGED1)
MERGED2 <- MERGED1[order(MERGED1$Clamp),]
# write.table(MERGED2,"RESIDUALS_FI_adj.csv")

  ## BMI adjusted
M0_res_bmi <- sapply(grep("^X0", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI, data = MAIN1)))
M30_res_bmi <- sapply(grep("^X30", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality+ BMI, data = MAIN1)))
M120_res_bmi <- sapply(grep("^X120", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality+ BMI, data = MAIN1)))
M0_res_bmi <- data.frame(M0_res_bmi)
M30_res_bmi <- data.frame(M30_res_bmi)
M120_res_bmi <- data.frame(M120_res_bmi)
MERGED1 <- cbind(MAIN1$Clamp, M0_res_bmi, M30_res_bmi, M120_res_bmi)
colnames(MERGED1) <- c("Clamp", Mnames0, Mnames30, Mnames120)
MERGED1 <- data.frame(MERGED1)
MERGED2 <- MERGED1[order(MERGED1$Clamp),]
# write.table(MERGED2,"RESIDUALS_BMI_adj.csv")

  ## FG adjusted
M0_res_fg <- sapply(grep("^X0", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality + FG, data = MAIN1)))
M30_res_fg <- sapply(grep("^X30", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality+ FG, data = MAIN1)))
M120_res_fg <- sapply(grep("^X120", names(MAIN1)), function(x) residuals(lm(MAIN1[,x] ~ Age + Storage_time + Thawed + Sample_type + Quality+ FG, data = MAIN1)))
M0_res_fg <- data.frame(M0_res_fg)
M30_res_fg <- data.frame(M30_res_fg)
M120_res_fg <- data.frame(M120_res_fg)
MERGED1 <- cbind(MAIN1$Clamp, M0_res_fg, M30_res_fg, M120_res_fg)
colnames(MERGED1) <- c("Clamp", Mnames0, Mnames30, Mnames120)
MERGED1 <- data.frame(MERGED1)
MERGED2 <- MERGED1[order(MERGED1$Clamp),]
# write.table(MERGED2,"RESIDUALS_FG_adj.csv")

  ## plot raw value per quartile
p <- palette(c("white", "grey90", "grey30")) 
m <- MAIN1[order(MAIN1$Clamp), ]
m$id <- 1:nrow(m) 
m$clampQuartile <- "Q1" 
m$clampQuartile[117:234] <- "Q2"
m$clampQuartile[235:351] <- "Q3"
m$clampQuartile[352:nrow(m)] <- "Q4"
# pdf("C12carnitine.pdf", width = 8, height = 8)
par(mfrow = c(1,1))
gpe2_a <- m[, c("id", "clampQuartile", "X0_M344.280T218.202", "X30_M344.280T218.202", "X120_M344.280T218.202")]
gpe2_b <- melt(gpe2_a, id.vars = c("id", "clampQuartile"))
gpe2 <- gpe2_b[order(gpe2_b$clampQuartile, gpe2_b$variable), ]
METABOplot(gpe2$clampQuartile, gpe2$value, gpe2$variable, -0, 20, "C12-carnitine adj.: age","Residuals","Clamp M/I quartile")
abline(0, 0)
# dev.off()

  ## Extract the metabolite values for all-combined and per clamp M/I quartile - mean, sd, se
# C10 ##  gpe2_a<-m[,c("id","clampQuartile","X0_M316.249T186.936","X30_M316.249T186.936","X120_M316.249T186.936")]
# C12 ##  gpe2_a<-m[,c("id","clampQuartile","X0_M344.280T218.202","X30_M344.280T218.202","X120_M344.280T218.202")]
gpe2_b <- melt(gpe2_a, id.vars = c("id", "clampQuartile"))
gpe2 <- gpe2_b[order(gpe2_b$clampQuartile, gpe2_b$variable), ]
gpe2$ALL <- 1
n <- nrow(m)
gpe2$var1<-rep(1,nrow(gpe2))
gpe2$var1[1:116]<-1; gpe2$var1[117:234]<-2; gpe2$var1[235:351]<-3; gpe2$var1[352:n]<-4
gpe2$var1[(n+1):(n+116)]<-5; gpe2$var1[(n+117):(n+234)]<-6; gpe2$var1[(n+235):(n+351)]<-7; gpe2$var1[(n+352):(n+n)]<-8
gpe2$var1[(936+1):(936+116)]<-9; gpe2$var1[(936+117):(936+234)]<-10; gpe2$var1[(936+235):(936+351)]<-11; gpe2$var1[(936+352):(936+n)]<-12
temp1 <- cbind(aggregate(gpe2$value, by = list(gpe2$var1), FUN = mean), (aggregate(gpe2$value, by = list(gpe2$var1), FUN = sd)))
names(temp1) <- c("","mean","","sd")
rownames(temp1) <- c("1q0","1q30","1q120","2q0","2q30","2q120","3q0","3q30","3q120","4q0","4q30","4q120")
temp1$se <- temp1$sd / sqrt(117)
temp1 <- temp1[, -c(1,3)] 
temp2 <- cbind(aggregate(gpe2$value, by = list(gpe2$variable), FUN = mean),(aggregate(gpe2$value, by = list(gpe2$variable), FUN = sd)))
names(temp2) <- c("","mean","","sd")
rownames(temp2) <- c("0min","30min","120min")
temp2$se <- temp2$sd / sqrt(468)
temp2 <- temp2[, -c(1,3)]
# write.table(temp1, "carnitineXX_quartile.txt")
# write.table(temp2, "carnitineXX_all.txt")

    ## per-quartile plots, example
fiResiduals <- read.table("RESIDUALS_confoundercleared_FIadj.csv", header = TRUE, quote = "\"") 
fiResiduals$id <- 1:nrow(fiResiduals) 
fiResiduals$clampQuartile <- "Q1"
fiResiduals$clampQuartile[117:234] <- "Q2"
fiResiduals$clampQuartile[235:351] <- "Q3"
fiResiduals$clampQuartile[352:nrow(fiResiduals)] <- "Q4" 
p <- palette(c("dodgerblue3", "brown2", "goldenrod1"))
# pdf("LysoPE(181).pdf", width = 8, height = 8)
par(mfrow=c(1,1))
gpe2_a <- fiResiduals[,c("id", "clampQuartile", "X0_M480.309T397.138", "X30_M480.309T397.138", "X120_M480.309T397.138")]
gpe2_b <- melt(gpe2_a, id.vars = c("id", "clampQuartile"))
gpe2 <- gpe2_b[order(gpe2_b$clampQuartile, gpe2_b$variable), ]
METABOplot(gpe2$clampQuartile, gpe2$value, gpe2$variable, -0.4, 0.4, "LysoPE(18:1) adj.: FI", "Residuals", "Clamp M/I quartile")
abline(0,0)
# legend("topleft", c(levels(x)), fill = p, bty! = "n", cex = 1.1, legend = c("0min","30min","120min"), bg = "white", title = "LysoPE(18:1) (FI adjusted)") 
# dev.off()

######################################################################################################################### ############ ############ #############
######################################################################################################################### ############ ############ #############
######################################################################################################################### ############ ############ #############
 