####Packages ####
library(semTools)
library(reshape)
library(s20x)
library(scales)
library(dplyr)
library(survey)
library(srvyr)
library(lavaan)
library(semPlot)
library(sem)
library(car)
library(psych)
library(naniar)
library(haven)
library(GPArotation)
library(lavaan.survey)
library(survey)
library(reshape2)
library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(foreign)
library(poLCA)
library(polycor)
library(Amelia)
library(mice)

##### Read in datasets ##### 
combined_8Y <- fullbaseline

DEP_8Y <- combined_8Y[c("FAMID", "DP39_Y8M1", "DP2_Y8M1" , "DP5_Y8M1","DP32_Y8M1","DP18_Y8M1")] #Note: not using DP3_M54M as will use this for sensitivity
psych::describe(DEP_8Y)

DEP_8Y <- DEP_8Y[!apply(DEP_8Y == 99, 1, any), ]

DEP_8Y$sum <- rowSums(DEP_8Y[c("DP39_Y8M1", "DP2_Y8M1" , "DP5_Y8M1","DP32_Y8M1","DP18_Y8M1")])
table(DEP_8Y$sum)

DEP_8Y <- na.omit(DEP_8Y)
DEP_8Y$Y8 <- 1

View(miss_case_summary(DEP_8Y[2:6]))
summary(miss_var_summary(DEP_8Y[2:6]))

#### Cronbach for 8Y ####
DEP_8Y <- apply(DEP17_missing, 2, as.numeric)

DEP_8Y <- as.data.frame(DEP_8Y)
psych::omega(DEP_8Y[,2:6])

psych::describe(DEP_8Y[,2:6])
table(DEP_8Y$DP39_Y8M1)
table(DEP_8Y$DP2_Y8M1)
table(DEP_8Y$DP32_Y8M1)
table(DEP_8Y$DP5_Y8M1)
table(DEP_8Y$DP18_Y8M1)

DEP_8Y$sum <- rowSums(DEP_8Y[, 2:6])

corr_matrix <- tetrachoric(DEP_8Y[,2:6], y = NULL, na.rm = TRUE)
print(corr_matrix)

######### Factorability #########
#KMO
KMO(corr_matrix$rho)
KMO(DEP_8Y[,2:6])

#Correlations (Again)
df <- reshape2::melt(corr_matrix$rho)

ggplot(df, aes(Var1,Var2, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", midpoint = 0)+
  coord_fixed()+
  theme_minimal()+
  labs(x = "Variables", y = "Variables", fill = "Correlation")+
  ggtitle("Tetrachoric Correlation Matrix Heatmap")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

#### EFA ####

par(mfrow=c(1,2), mar=c(5,5,2,2))
factnumber <- psych::fa.parallel(DEP_8Y[,2:6], fm = 'pa', fa = 'fa', cor = "tet")

print(factnumber)
set.seed(123)
Loadings <- fa(DEP_8Y[,2:6], nfactors = 1, rotate = 'oblimin', fm = 'pa', cor = "tet")
print(Loadings$loadings,cutoff = 0.3)


DEP_8Y <- as.data.frame(DEP_8Y)

EFA <- DEP_8Y[c("DP32_Y8M1","DP33_Y8M1","DP34_Y8M1","DP35_Y8M1", "DP37_Y8M1",
                  "DP5_Y8M1", "DP39_Y8M1","DP2_Y8M1", "DP8_Y8M1",
                  "DP10_Y8M1", "DP51_Y8M1","DP42_Y8M1",
                  "DP11_Y8M1","DP45_Y8M1","DP13_Y8M1","DP46_Y8M1","DP52_Y8M1")]

model <- 'DEP =~ DP32_Y8M1 + DP33_Y8M1 + DP34_Y8M1 + DP35_Y8M1 + DP37_Y8M1 +
                 DP5_Y8M1 + DP39_Y8M1 + DP2_Y8M1 + DP8_Y8M1 +
                 DP10_Y8M1 + DP51_Y8M1 + DP42_Y8M1 +
                 DP11_Y8M1 + DP45_Y8M1 + DP13_Y8M1  + DP46_Y8M1 + DP52_Y8M1
DP32_Y8M1 ~~ DP33_Y8M1
DP13_Y8M1 ~~ DP52_Y8M1
DP51_Y8M1 ~~ DP42_Y8M1
DP13_Y8M1 ~~ DP46_Y8M1
DP35_Y8M1 ~~ DP37_Y8M1
DP46_Y8M1 ~~ DP52_Y8M1
DP5_Y8M1 ~~ DP42_Y8M1
DP35_Y8M1 ~~  DP5_Y8M1'


fit <- lavaan::cfa(model, data = EFA) #Confirmatory Factor Analysis
summary(fit ,fit.measures=TRUE, standardized=TRUE) #Summary stats for CFA.


modificationIndices(fit, sort = T)

#### Measurement invariance for 8 ####
DEP_8Y <- left_join(DEP_8Y , combined_8Y, by = "FAMID", suffix = c("", ".x"))
DEP_8Y <- DEP_8Y %>% dplyr::select(-contains(".x"))

DEP_8Y <- left_join(DEP_8Y , mother_AN, by = "MID", suffix = c("", ".x"))
DEP_8Y <- DEP_8Y %>% dplyr::select(-contains(".x"))

DEP_8Y_mi <- DEP_8Y[c("FAMID", "FAMID", "DP32_Y8M1","DP33_Y8M1","DP34_Y8M1","DP35_Y8M1", "DP37_Y8M1",
                          "DP5_Y8M1", "DP39_Y8M1","DP2_Y8M1", "DP8_Y8M1",
                          "DP10_Y8M1", "DP51_Y8M1","DP42_Y8M1",
                          "DP11_Y8M1","DP45_Y8M1","DP13_Y8M1","DP46_Y8M1","DP52_Y8M1", 
                          "FIN56_Y8M", "EDALL_AM", "AGE_GROUP_AM", "NZDEP2013_Y8M", 
                          "ETH5_NZDER_AM",
                          "ETH5_E_AM", "ETH5_O_AM","ETH5_MELA_AM", "ETH5_P_AM", "ETH5_A_AM", "ETH5_M_AM", 
                          "HH6_Y8M", "PQ5_Y8M")]

###Ethnicity
library(tidyverse)
DEP_8Y_mi$Ethnicity_AN[DEP_8Y_mi$ETH5_NZDER_AM ==1] <- "European"
DEP_8Y_mi$Ethnicity_AN[DEP_8Y_mi$ETH5_E_AM==1] <- "European"
DEP_8Y_mi$Ethnicity_AN[DEP_8Y_mi$ETH5_O_AM==1] <- "Other"
DEP_8Y_mi$Ethnicity_AN[DEP_8Y_mi$ETH5_MELA_AM==1] <- "Other"
DEP_8Y_mi$Ethnicity_AN[DEP_8Y_mi$ETH5_P_AM==1] <- "Pacific"
DEP_8Y_mi$Ethnicity_AN[DEP_8Y_mi$ETH5_A_AM==1] <- "Asian"
DEP_8Y_mi$Ethnicity_AN[DEP_8Y_mi$ETH5_M_AM==1] <- "M\u101ori"

measurementInvariance(model = model, data=DEP_8Y_mi, 
                      group = "Ethnicity_AN")

###Household tenure
table(DEP_8Y_mi$HH6_Y8M, exclude = NULL)
DEP_8Y_mi$HH6_Y8M <- recode(DEP_8Y_mi$HH6_Y8M, "c(99) = NA; c(98) = NA")
DEP_8Y_mi$HH6_Y8M <- as.factor(DEP_8Y_mi$HH6_Y8M)

# model comparison tests
measurementInvariance(model = model, data=DEP_8Y_mi, 
                      group = "HH6_Y8M") #Configural = Configural; Loadings = Weak; Intercepts = Strong; Means = Strict

### Relationship status
table(DEP_8Y_mi$PQ5_Y8M) 
DEP_8Y_mi$PQ5_Y8M <- recode(DEP_8Y_mi$PQ5_Y8M, "c(98) = NA")
DEP_8Y_mi$PQ5_Y8M <- as.factor(DEP_8Y_mi$PQ5_Y8M)

measurementInvariance(model = model, data=DEP_8Y_mi, 
                      group = "PQ5_Y8M")

##Factor scores
DEP_8Y_mi$refined.scores <- lavPredict(fit, newdata=DEP_8Y_mi, type= "lv", method = "regression" )
DEP_8Y_mi$sum.scores <- rowSums(DEP_8Y_mi[c("DP32_Y8M1","DP33_Y8M1","DP34_Y8M1","DP35_Y8M1", "DP37_Y8M1",
                                                "DP5_Y8M1", "DP39_Y8M1","DP2_Y8M1", "DP8_Y8M1",
                                                "DP10_Y8M1", "DP51_Y8M1","DP42_Y8M1",
                                                "DP11_Y8M1","DP45_Y8M1","DP13_Y8M1","DP46_Y8M1","DP52_Y8M1")])
sum.scores <- DEP_8Y_mi$sum.scores

table(sum.scores, exclude = NULL)
View(table(DEP_8Y_mi$refined.scores, exclude = NULL))

summary(lm(refined.scores ~ sum.scores, drop.unused.levels = TRUE))
plot(refined.scores, sum.scores)
cor.test(refined.scores, sum.scores, method = c("pearson"))


#### Distributions of Sum Scores at each time point ####

#8-year


DEP_8Y_mi %>%
  ggplot( aes(x=refined.scores)) +
  theme(panel.background = element_rect(fill='transparent')) +
  geom_histogram(fill="#fdb913", color="#fdb913", alpha=1, binwidth = 0.05)+
  theme_ipsum()+
  theme_classic()+
  labs(x = "Hardship Score", y = "Count")+
  ggtitle("Hardship Scores at 8-year")

####Create cut points each time point ####
DEP_8Y_mi$sum <- rowSums((DEP_8Y_mi[,3:19]))
table(DEP_8Y_mi$sum, exclude = NULL)
DEP_8Y_mi$categorical1[DEP_8Y_mi$sum >= 6] <- "material hardship"
DEP_8Y_mi$categorical1[DEP_8Y_mi$sum >= 9] <- "severe material hardship"
DEP_8Y_mi$categorical1[DEP_8Y_mi$sum < 6] <- "no/little material hardship"
prop.table(table(DEP_8Y_mi$categorical1))
table(DEP_8Y_mi$categorical1, exclude = NULL) 

DEP_8Y_mi$zscore <- scale(DEP_8Y_mi$sum)
View(DEP_8Y_mi$zscore)

options(digits=10)

####Latent Class####
LCA_8Y <- as.data.frame(DEP_8Y)

View(miss_case_summary(LCA_8Y))

LCA_8Y[,2:6][LCA_8Y[,2:6] >= 1] <- 2
LCA_8Y[,2:6][LCA_8Y[,2:6] < 1] <- 1

f <- cbind(DP39_Y8M1,DP5_Y8M1,DP2_Y8M1,DP32_Y8M1, DP18_Y8M1)~1

# Set up storage for AIC, BIC , adjusted BIC values
AIC <- rep(NA, 5)
BIC <- rep(NA, 5)
ADJBIC <- rep(NA, 5)
AWE <- rep(NA, 5)

# Fit LCA models with 1 to 10 classes and extract AIC and BIC values
set.seed(123)
for (i in 1:5) {
  fit <- poLCA(f, data = LCA_8Y[,2:6], nclass = i, maxiter = 1000)
  AIC[i] <- fit$aic
  BIC[i] <- fit$bic
  ADJBIC[i] <- {bic <- fit$bic
  sample_size <- fit$N
  df <- fit$npar
  adjusted_bic <- bic + (log(sample_size)/2) * df}
  AWE[i] <- {awe <- (-2 * fit$llik) + (fit$npar *(log(fit$N) + 1.5))}
}

AIC
BIC
ADJBIC
AWE

# Plot AIC and BIC against number of classes
plot(1:5, AIC, type = "b", pch = 16, col = "blue", xlab = "Number of Classes",
     ylab = "Value", main = "AIC, BIC, ADJBIC, and AWE")
lines(1:5, BIC, type = "b", pch = 16, col = "red")
lines(1:5, ADJBIC, type = "b", pch = 16, col = "green")
lines(1:5, AWE, type = "b", pch = 16, col = "orange")
legend("topright", legend = c("AIC", "BIC", "ADJBIC", "AWE"), col = c("blue", "red", "green", "orange"), pch = 16)

# Identify the number of classes with the lowest AIC and BIC values
min_AIC <- which.min(AIC)
min_BIC <- which.min(BIC)
min_ADJBIC <- which.min(ADJBIC)
min_AWE <- which.min(AWE)
cat("The optimal number of classes based on AIC is", min_AIC, "and based on BIC is", min_BIC, "and based on ADJBIC is", min_ADJBIC, 
    "and based on AWE is", min_AWE, "\n")

# Plot AIC, BIC, ADJBIC, and AWE
plot(1:5, AIC, type = "b", pch = 16, col = "blue", xlab = "Number of Classes",
     ylab = "Value", main = "AIC, BIC, ADJBIC, and AWE")
lines(1:5, BIC, type = "b", pch = 16, col = "red")
lines(1:5, ADJBIC, type = "b", pch = 16, col = "green")
lines(1:5, AWE, type = "b", pch = 16, col = "orange")
legend("topright", legend = c("AIC", "BIC", "ADJBIC", "AWE"), col = c("blue", "red", "green", "orange"), pch = 16)


set.seed(123)
fit_2 <- poLCA(f, data = LCA_8Y, nclass = 2, maxiter = 10000)
table(fit_2$predclass)

LCA_8Y <- cbind(LCA_8Y, "LCA_2" = fit_2$predclass)
table(LCA_8Y$LCA_2, exclude = NULL)

LCA_8Y$DP2_Y8M1 <- recode(LCA_8Y$DP2_Y8M1, "c(1) = 0; c(2) = 1")
LCA_8Y$DP5_Y8M1 <- recode(LCA_8Y$DP5_Y8M1, "c(1) = 0; c(2) = 1")
LCA_8Y$DP39_Y8M1 <- recode(LCA_8Y$DP39_Y8M1, "c(1) = 0; c(2) = 1")
LCA_8Y$DP32_Y8M1 <- recode(LCA_8Y$DP32_Y8M1, "c(1) = 0; c(2) = 1")
LCA_8Y$DP18_Y8M1 <- recode(LCA_8Y$DP18_Y8M1, "c(1) = 0; c(2) = 1")

LCA_8Y$sum <- rowSums(LCA_8Y[,2:6])
table(LCA_8Y$sum, LCA_8Y$LCA_2)
table(LCA_8Y$LCA_2)
prop.table(table(LCA_8Y$sum))*100
plot <- melt(fit_2$probs)
plot

plot$test[plot$Var2 == "Pr(1)"] <- "no hardship"
plot$test[plot$Var2 == "Pr(2)"] <- "hardship"
plot$L1[plot$L1 == "DP39_Y8M1"] <- "Buy cheaper or less meat "
plot$L1[plot$L1 == "DP5_Y8M1"] <- "Gone without fresh fruit and vegetables"
plot$L1[plot$L1 == "DP2_Y8M1"] <- "Put up with feeling cold"
plot$L1[plot$L1 == "DP32_Y8M1"] <- "Do not have two pairs of good shoes "
plot$L1[plot$L1 == "DP18_Y8M1"] <- "Received help"

zp2 <- ggplot(plot,
aes(x = Var1, y = value, fill = test))+ 
  geom_bar(stat = "identity", position = "stack")+ 
 facet_wrap(~ L1)+ 
 scale_x_discrete("Class", expand = c(0, 0)) + 
scale_y_continuous("Proportion", expand = c(0, 0))+ 
scale_fill_discrete(name = "Factor Level")+ 
theme_bw()
zp2

##### Any hardship #####

DEP_8Y <- as.data.frame(DEP_8Y)

DEP_8Y$anyhardship[DEP_8Y$DP39_Y8M1 < 0.5] <- 0
DEP_8Y$anyhardship[DEP_8Y$DP2_Y8M1 < 0.5] <- 0
DEP_8Y$anyhardship[DEP_8Y$DP32_Y8M1 < 0.5] <- 0
DEP_8Y$anyhardship[DEP_8Y$DP5_Y8M1 < 0.5] <- 0
DEP_8Y$anyhardship[DEP_8Y$DP18_Y8M1 < 0.5] <- 0
DEP_8Y$anyhardship[DEP_8Y$DP39_Y8M1 >= 0.5] <- 1
DEP_8Y$anyhardship[DEP_8Y$DP2_Y8M1 >= 0.5] <- 1
DEP_8Y$anyhardship[DEP_8Y$DP32_Y8M1 >= 0.5] <- 1
DEP_8Y$anyhardship[DEP_8Y$DP5_Y8M1 >= 0.5] <- 1
DEP_8Y$anyhardship[DEP_8Y$DP18_Y8M1 >= 0.5] <- 1

table(DEP_8Y$anyhardship)
prop.table(table(DEP_8Y$anyhardship))

### Cut points ####
DEP_8Y <- as.data.frame(DEP_8Y)
DEP_8Y$sum <- rowSums(DEP_8Y[,2:6])

combined_8Y <- dplyr::left_join(DEP_8Y, mother_8y, by = "FAMID", suffix = c("", ".a"))
combined_8Y <- combined_8Y %>% dplyr::select(-contains(".a"))

combined_8Y <- combined_8Y[c("DP31_Y8M",#food bank
                             "DP52_Y8M", #borrow money == 3
                             "FIN57_5_Y8M",#job seeker
                             "FIN57_6_Y8M" ,#domestic purposes/sole parent support 
                             "FIN57_7_Y8M" ,#sickness
                             "FIN57_17_Y8M" ,#dsiability
                             "FIN57_11_Y8M" , #student allowance
                             "FIN57_15_Y8M",#family tax credits
                             "FIN57_18_Y8M", #OSCAr (low or middle income)
                             "FIN57_10_Y8M", #accomodation supplement
                             "FIN57_19_Y8M", #training incentive
                             "FIN57_20_Y8M", #public housing/rent subsiy
                             "sum" )]
table(combined_8Y$DP52_Y8M, exclude = NULL)

combined_8Y$foodbank[combined_8Y$DP31_Y8M == 0] <- 0
combined_8Y$foodbank[combined_8Y$DP31_Y8M >= 1] <- 1
combined_8Y$foodbank[combined_8Y$DP31_Y8M == 99] <- NA
combined_8Y$borrow[combined_8Y$DP52_Y8M <= 2] <- 0
combined_8Y$borrow[combined_8Y$DP52_Y8M == 3 ] <- 1
combined_8Y$borrow[combined_8Y$DP52_Y8M == 99] <- NA

combined_8Y$benefit[combined_8Y$FIN57_5_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_6_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_7_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_17_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_11_Y8M == 0] <- 0

combined_8Y$benefit[combined_8Y$FIN57_18_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_10_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_19_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_20_Y8M == 0] <- 0

combined_8Y$benefit[combined_8Y$FIN57_5_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_6_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_7_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_17_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_11_Y8M == 1] <- 1

combined_8Y$benefit[combined_8Y$FIN57_18_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_10_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_19_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_20_Y8M == 1] <- 1


combined_8Y$famtax[combined_8Y$FIN57_15_Y8M == 0] <- 0
combined_8Y$famtax[combined_8Y$FIN57_15_Y8M == 1] <- 1


test <- combined_8Y[c("sum", "benefit", "foodbank", "borrow", "famtax")]
test <- na.omit(test)

test$economic[test$foodbank == 1] <- 1
test$economic[test$foodbank == 0] <- 0
test$economic[test$borrow == 1] <- 1
test$economic[test$borrow == 0] <- 0

test$validation[test$benefit == 1 & test$economic == 1 & test$famtax == 1 ] <- "Benefit/economic"
test$validation[test$benefit == 1 & test$economic == 1 & test$famtax == 0 ] <- "Benefit/economic"

test$validation[test$benefit == 0 & test$economic == 1 & test$famtax == 1 ] <- "Famtax/economic"

test$validation[test$benefit == 0 & test$economic == 1 & test$famtax == 0] <- "Economic only"

test$validation[test$benefit == 1 & test$economic == 0 & test$famtax == 1] <- "Benefit only"
test$validation[test$benefit == 1 & test$economic == 0 & test$famtax == 0] <- "Benefit only"

test$validation[test$benefit == 0 & test$economic == 0 & test$famtax == 1] <- "Famtax only"
test$validation[test$benefit == 0 & test$economic == 0 & test$famtax == 0] <- "None"

test$validation <- factor(test$validation, levels = c("Benefit/economic", "Famtax/economic", "Economic only",  "Benefit only", "Famtax only", "None"))

table(test$validation, exclude = NULL)
table(test$sum)
table(test$validation, test$sum, exclude = NULL)

### Robustness check 1 - lower threshold for borrowing ####
DEP_8Y <- as.data.frame(DEP_8Y)
DEP_8Y$sum <- rowSums(DEP_8Y[,2:6])

combined_8Y <- dplyr::left_join(DEP_8Y, mother_8y, by = "FAMID", suffix = c("", ".a"))
combined_8Y <- combined_8Y %>% dplyr::select(-contains(".a"))

combined_8Y <- combined_8Y[c("DP31_Y8M",#food bank
                             "DP52_Y8M", #borrow money == 3
                             "FIN57_5_Y8M",#job seeker
                             "FIN57_6_Y8M" ,#domestic purposes/sole parent support 
                             "FIN57_7_Y8M" ,#sickness
                             "FIN57_17_Y8M" ,#dsiability
                             "FIN57_11_Y8M" , #student allowance
                             "FIN57_15_Y8M",#family tax credits
                             "FIN57_18_Y8M", #OSCAr (low or middle income)
                             "FIN57_10_Y8M", #accomodation supplement
                             "FIN57_19_Y8M", #training incentive
                             "FIN57_20_Y8M", #public housing/rent subsiy
                             "sum" )]
table(combined_8Y$DP52_Y8M, exclude = NULL)

combined_8Y$foodbank[combined_8Y$DP31_Y8M == 0] <- 0
combined_8Y$foodbank[combined_8Y$DP31_Y8M >= 1] <- 1
combined_8Y$foodbank[combined_8Y$DP31_Y8M == 99] <- NA
combined_8Y$borrow[combined_8Y$DP52_Y8M == 1] <- 0 #Lower
combined_8Y$borrow[combined_8Y$DP52_Y8M >= 2 ] <- 1
combined_8Y$borrow[combined_8Y$DP52_Y8M == 99] <- NA

combined_8Y$benefit[combined_8Y$FIN57_5_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_6_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_7_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_17_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_11_Y8M == 0] <- 0

combined_8Y$benefit[combined_8Y$FIN57_18_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_10_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_19_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_20_Y8M == 0] <- 0

combined_8Y$benefit[combined_8Y$FIN57_5_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_6_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_7_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_17_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_11_Y8M == 1] <- 1

combined_8Y$benefit[combined_8Y$FIN57_18_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_10_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_19_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_20_Y8M == 1] <- 1


combined_8Y$famtax[combined_8Y$FIN57_15_Y8M == 0] <- 0
combined_8Y$famtax[combined_8Y$FIN57_15_Y8M == 1] <- 1


test <- combined_8Y[c("sum", "benefit", "foodbank", "borrow", "famtax")]
test <- na.omit(test)

test$economic[test$foodbank == 1] <- 1
test$economic[test$foodbank == 0] <- 0
test$economic[test$borrow == 1] <- 1
test$economic[test$borrow == 0] <- 0

test$validation[test$benefit == 1 & test$economic == 1 & test$famtax == 1 ] <- "Benefit/economic"
test$validation[test$benefit == 1 & test$economic == 1 & test$famtax == 0 ] <- "Benefit/economic"

test$validation[test$benefit == 0 & test$economic == 1 & test$famtax == 1 ] <- "Famtax/economic"

test$validation[test$benefit == 0 & test$economic == 1 & test$famtax == 0] <- "Economic only"

test$validation[test$benefit == 1 & test$economic == 0 & test$famtax == 1] <- "Benefit only"
test$validation[test$benefit == 1 & test$economic == 0 & test$famtax == 0] <- "Benefit only"

test$validation[test$benefit == 0 & test$economic == 0 & test$famtax == 1] <- "Famtax only"
test$validation[test$benefit == 0 & test$economic == 0 & test$famtax == 0] <- "None"

test$validation <- factor(test$validation, levels = c("Benefit/economic", "Famtax/economic", "Economic only",  "Benefit only", "Famtax only", "None"))

table(test$validation, exclude = NULL)
table(test$sum)
table(test$validation, test$sum, exclude = NULL)

### Robustness check 2 - separating ####
DEP_8Y <- as.data.frame(DEP_8Y)
DEP_8Y$sum <- rowSums(DEP_8Y[,2:6])

combined_8Y <- dplyr::left_join(DEP_8Y, mother_8y, by = "FAMID", suffix = c("", ".a"))
combined_8Y <- combined_8Y %>% dplyr::select(-contains(".a"))

combined_8Y <- combined_8Y[c("DP31_Y8M",#food bank
                             "DP52_Y8M", #borrow money == 3
                             "FIN57_5_Y8M",#job seeker
                             "FIN57_6_Y8M" ,#domestic purposes/sole parent support 
                             "FIN57_7_Y8M" ,#sickness
                             "FIN57_17_Y8M" ,#dsiability
                             "FIN57_11_Y8M" , #student allowance
                             "FIN57_15_Y8M",#family tax credits
                             "FIN57_18_Y8M", #OSCAr (low or middle income)
                             "FIN57_10_Y8M", #accomodation supplement
                             "FIN57_19_Y8M", #training incentive
                             "FIN57_20_Y8M", #public housing/rent subsiy
                             "sum" )]
table(combined_8Y$DP52_Y8M, exclude = NULL)

combined_8Y$foodbank[combined_8Y$DP31_Y8M == 0] <- 0
combined_8Y$foodbank[combined_8Y$DP31_Y8M >= 1] <- 1
combined_8Y$foodbank[combined_8Y$DP31_Y8M == 99] <- NA
combined_8Y$borrow[combined_8Y$DP52_Y8M <= 2] <- 0
combined_8Y$borrow[combined_8Y$DP52_Y8M == 3 ] <- 1
combined_8Y$borrow[combined_8Y$DP52_Y8M == 99] <- NA

combined_8Y$benefit[combined_8Y$FIN57_5_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_6_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_7_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_17_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_11_Y8M == 0] <- 0

combined_8Y$benefit[combined_8Y$FIN57_18_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_10_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_19_Y8M == 0] <- 0
combined_8Y$benefit[combined_8Y$FIN57_20_Y8M == 0] <- 0

combined_8Y$benefit[combined_8Y$FIN57_5_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_6_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_7_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_17_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_11_Y8M == 1] <- 1

combined_8Y$benefit[combined_8Y$FIN57_18_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_10_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_19_Y8M == 1] <- 1
combined_8Y$benefit[combined_8Y$FIN57_20_Y8M == 1] <- 1


combined_8Y$famtax[combined_8Y$FIN57_15_Y8M == 0] <- 0
combined_8Y$famtax[combined_8Y$FIN57_15_Y8M == 1] <- 1


test <- combined_8Y[c("sum", "benefit", "foodbank", "borrow", "famtax")]
test <- na.omit(test)

test$validation[test$benefit == 1 & test$borrow == 1 & test$foodbank == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 0 & test$borrow == 1 & test$foodbank == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 1 & test$borrow == 1 & test$foodbank == 1 & test$famtax == 0] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"

test$validation[test$benefit == 1 & test$borrow == 0 & test$foodbank == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 0 & test$borrow == 0 & test$foodbank == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 1 & test$borrow == 0 & test$foodbank == 1 & test$famtax == 0] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"

test$validation[test$benefit == 0 & test$borrow == 1 & test$foodbank == 1 & test$famtax == 0] <- "Either economic strain experiences & no benefit receipt"
test$validation[test$benefit == 0 & test$borrow == 0 & test$foodbank == 1 & test$famtax == 0] <- "Either economic strain experiences & no benefit receipt"
test$validation[test$benefit == 0 & test$borrow == 1 & test$foodbank == 0 & test$famtax == 0] <- "Either economic strain experiences & no benefit receipt"

test$validation[test$benefit == 1 & test$borrow == 1 & test$foodbank == 0 & test$famtax == 0] <- "borrow and Benefit receipt"
test$validation[test$benefit == 0 & test$borrow == 1 & test$foodbank == 0 & test$famtax == 1] <- "borrow and Benefit receipt"
test$validation[test$benefit == 1 & test$borrow == 1 & test$foodbank == 0 & test$famtax == 1] <- "borrow and Benefit receipt"

test$validation[test$benefit == 1 & test$borrow == 0 & test$foodbank == 0 & test$famtax == 1] <- "Benefit receipt only"
test$validation[test$benefit == 1 & test$borrow == 0 & test$foodbank == 0 & test$famtax == 0] <- "Benefit receipt only"
test$validation[test$benefit == 1 & test$borrow == 0 & test$foodbank == 0 & test$famtax == 1] <- "Benefit receipt only"

test$validation[test$benefit == 0 & test$borrow == 0 & test$foodbank == 0 & test$famtax == 1] <- "Fam tax only"
test$validation[test$benefit == 0 & test$borrow == 0 & test$foodbank == 0 & test$famtax == 0] <- "None"


test$validation <- factor(test$validation, levels = c("Food bank & benefit receipt / Both economic strain experiences and benefit receipt", 
                                                      "Either economic strain experiences & no benefit receipt",
                                                      "borrow and Benefit receipt", "Benefit receipt only",
                                                      "Fam tax only", "None"))

table(test$validation, exclude = NULL)
table(test$sum)
table(test$validation, test$sum, exclude = NULL)

### Applying the cut point ####
DEP_8Y$sum <- rowSums(DEP_8Y[c("DP39_Y8M1", "DP2_Y8M1" , "DP5_Y8M1","DP32_Y8M1","DP18_Y8M1")])
table(DEP_8Y$sum)
DEP_8Y$threshold8[DEP_8Y$sum >=2] <- 1
DEP_8Y$threshold8[DEP_8Y$sum < 2] <- 0
table(DEP_8Y$threshold8, exclude = NULL)

#### Comparison of all three methods ####
DEP_8Y <- dplyr::left_join(DEP_8Y, LCA_8Y, by = "FAMID", suffix = c("", ".a"))
DEP_8Y <- DEP_8Y %>% dplyr::select(-contains(".a"))

table(DEP_8Y$LCA_2)
table(DEP_8Y$anyhardship)
table(DEP_8Y$threshold8)
DEP_8Y$LCA_8Y_revised[DEP_8Y$LCA_2 == 1] <- 0
DEP_8Y$LCA_8Y_revised[DEP_8Y$LCA_2 == 2] <- 1

table(DEP_8Y$threshold8 , DEP_8Y$LCA_8Y_revised)
table(DEP_8Y$threshold8 , DEP_8Y$anyhardship)
table(DEP_8Y$LCA_8Y_revised , DEP_8Y$anyhardship)

DEP_8Y$congruence <- ifelse(DEP_8Y$anyhardship == 1 & DEP_8Y$threshold8 == 1 & DEP_8Y$LCA_8Y_revised == 1, 1, 
                             ifelse(DEP_8Y$anyhardship == 0 & DEP_8Y$threshold8 == 0 & DEP_8Y$LCA_8Y_revised == 0, 0, 2))
table(DEP_8Y$congruence)
