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
combined_12Y <- fullbaseline

DEP_12Y <- combined_12Y[c("FAMID", "DP39_Y12M1", "DP2_Y12M1" , "DP5_Y12M1","DP32_Y12M1","DP18_Y12M1")]
psych::describe(DEP_12Y)

DEP_12Y <- DEP_12Y[!apply(DEP_12Y == 99, 1, any), ]

DEP_12Y$sum <- rowSums(DEP_12Y[c("DP39_Y12M1", "DP2_Y12M1" , "DP5_Y12M1","DP32_Y12M1","DP18_Y12M1")])
table(DEP_12Y$sum)

DEP_12Y <- na.omit(DEP_12Y)
DEP_12Y$Y12 <- 1

View(miss_case_summary(DEP_12Y[2:6]))
summary(miss_var_summary(DEP_12Y[2:6]))

DEP17_missing12 <- DEP_12Y[which(rowMeans(!is.na(DEP_12Y[2:6])) > 0.6), ] #With >40% missing removed.
data <- DEP17_missing

describe(DEP_12Y[2:6])

mcar_test(DEP17_8Y[2:6])


#### Cronbach for 12Y ####
DEP_12Y <- apply(DEP17_missing12, 2, as.numeric)
DEP_12Y <- as.data.frame(DEP_12Y)

psych::describe(DEP_12Y[,2:6])
table(DEP_12Y$DP39_Y12M1)
table(DEP_12Y$DP2_Y12M1)
table(DEP_12Y$DP32_Y12M1)
table(DEP_12Y$DP5_Y12M1)
table(DEP_12Y$DP18_Y12M1)

DEP_12Y$sum <- rowSums(DEP_12Y[, 2:6])

corr_matrix <- tetrachoric(DEP_12Y[,2:6], y = NULL, na.rm = TRUE)
print(corr_matrix)

######### Factorability #########
#KMO
KMO(corr_matrix$rho)


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
factnumber <- psych::fa.parallel(DEP_12Y[,2:6], fm = 'pa', fa = 'fa', cor = "tet")

print(factnumber)
set.seed(123)
Loadings <- fa(DEP_12Y[,2:6], nfactors = 1, rotate = 'oblimin', fm = 'pa', cor = "tet")
print(Loadings$loadings,cutoff = 0.3)


#### Measurement invariance for 12 ####

####Latent Class####
LCA_12Y <- as.data.frame(DEP_12Y)
View(miss_case_summary(LCA_12Y))

LCA_12Y[,2:6][LCA_12Y[,2:6] >= 1] <- 2
LCA_12Y[,2:6][LCA_12Y[,2:6] < 1] <- 1
#write.csv(LCA_12Y, "LCA_12Y.csv")

f <- cbind(DP39_Y12M1,DP5_Y12M1,DP2_Y12M1,DP32_Y12M1, DP18_Y12M1)~1

# Set up storage for AIC, BIC , adjusted BIC values
AIC <- rep(NA, 5)
BIC <- rep(NA, 5)
ADJBIC <- rep(NA, 5)
AWE <- rep(NA, 5)

# Fit LCA models with 1 to 10 classes and extract AIC and BIC values
set.seed(123)
for (i in 1:5) {
  fit <- poLCA(f, data = LCA_12Y[,2:6], nclass = i, maxiter = 1000)
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
fit_2 <- poLCA(f, data = LCA_12Y, nclass = 2, maxiter = 10000)
table(fit_2$predclass)


LCA_12Y <- cbind(LCA_12Y, "LCA_2" = fit_2$predclass)
table(LCA_12Y$LCA_2, exclude = NULL)

LCA_12Y$DP2_Y12M1 <- recode(LCA_12Y$DP2_Y12M1, "c(1) = 0; c(2) = 1")
LCA_12Y$DP5_Y12M1 <- recode(LCA_12Y$DP5_Y12M1, "c(1) = 0; c(2) = 1")
LCA_12Y$DP39_Y12M1 <- recode(LCA_12Y$DP39_Y12M1, "c(1) = 0; c(2) = 1")
LCA_12Y$DP32_Y12M1 <- recode(LCA_12Y$DP32_Y12M1, "c(1) = 0; c(2) = 1")
LCA_12Y$DP18_Y12M1 <- recode(LCA_12Y$DP18_Y12M1, "c(1) = 0; c(2) = 1")

LCA_12Y$sum <- rowSums(LCA_12Y[,2:6])
table(LCA_12Y$sum, LCA_12Y$LCA_2)
table(LCA_12Y$LCA_2)
prop.table(table(LCA_12Y$sum))*100

##### Any hardship #####

DEP_12Y <- as.data.frame(DEP_12Y)

DEP_12Y$anyhardship[DEP_12Y$DP39_Y12M1 < 0.5] <- 0
DEP_12Y$anyhardship[DEP_12Y$DP2_Y12M1 < 0.5] <- 0
DEP_12Y$anyhardship[DEP_12Y$DP32_Y12M1 < 0.5] <- 0
DEP_12Y$anyhardship[DEP_12Y$DP5_Y12M1 < 0.5] <- 0
DEP_12Y$anyhardship[DEP_12Y$DP18_Y12M1 < 0.5] <- 0
DEP_12Y$anyhardship[DEP_12Y$DP39_Y12M1 >= 0.5] <- 1
DEP_12Y$anyhardship[DEP_12Y$DP2_Y12M1 >= 0.5] <- 1
DEP_12Y$anyhardship[DEP_12Y$DP32_Y12M1 >= 0.5] <- 1
DEP_12Y$anyhardship[DEP_12Y$DP5_Y12M1 >= 0.5] <- 1
DEP_12Y$anyhardship[DEP_12Y$DP18_Y12M1 >= 0.5] <- 1

table(DEP_12Y$anyhardship)
prop.table(table(DEP_12Y$anyhardship))*100
DEP_12Y$sum <- rowSums(DEP_12Y[c("DP39_Y12M1", "DP2_Y12M1", "DP5_Y12M1", "DP18_Y12M1", "DP32_Y12M1")])
table(DEP_12Y$sum)


### Cut points ####

DEP_12Y$sum <- rowSums(DEP_12Y[,2:6])

combined_12Y <- as.data.frame(DEP_12Y)


combined_12Y <- dplyr::left_join(DEP_12Y, dcw12_mother, by = "FAMID", suffix = c("", ".a"))
combined_12Y <- combined_12Y %>% dplyr::select(-contains(".a"))

combined_12Y <- combined_12Y[c("DP31_Y12M",#food bank
                               "DP52_Y12M", #borrow money == 3
                               "FIN57_5_Y12M",#job seeker 
                               "FIN57_6_Y12M",#sole parent
                               "FIN57_7_Y12M" ,#sickness
                               "FIN57_10_Y12M", #accomodation supplement
                               "FIN57_11_Y12M" , #student allowance
                               "FIN57_17_Y12M" ,#dsiability
                               "FIN57_15_Y12M",#family tax credits
                               "FIN57_18_Y12M", #OSCAr (low or middle income)
                               "FIN57_19_Y12M", #training incentive
                               "FIN57_20_Y12M", #public housing/rent subsiy
                               "sum" )]
table(combined_12Y$DP52_Y12M, exclude = NULL)

combined_12Y$foodbank[combined_12Y$DP31_Y12M == 0] <- 0
combined_12Y$foodbank[combined_12Y$DP31_Y12M >= 1] <- 1
combined_12Y$foodbank[combined_12Y$DP31_Y12M == 99] <- NA
combined_12Y$borrow[combined_12Y$DP52_Y12M <= 2] <- 0
combined_12Y$borrow[combined_12Y$DP52_Y12M == 3 ] <- 1
combined_12Y$borrow[combined_12Y$DP52_Y12M == 99] <- NA

combined_12Y$benefit[combined_12Y$FIN57_5_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_6_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_7_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_17_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_11_Y12M == 0] <- 0

combined_12Y$benefit[combined_12Y$FIN57_18_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_10_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_19_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_20_Y12M == 0] <- 0

combined_12Y$benefit[combined_12Y$FIN57_5_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_6_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_7_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_17_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_11_Y12M == 1] <- 1

combined_12Y$benefit[combined_12Y$FIN57_18_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_10_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_19_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_20_Y12M == 1] <- 1


combined_12Y$famtax[combined_12Y$FIN57_15_Y12M == 0] <- 0
combined_12Y$famtax[combined_12Y$FIN57_15_Y12M == 1] <- 1


test <- combined_12Y[c("sum", "benefit", "foodbank", "borrow", "famtax")]
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

combined_12Y <- dplyr::left_join(DEP_12Y, dcw12_mother, by = "FAMID", suffix = c("", ".a"))
combined_12Y <- combined_12Y %>% dplyr::select(-contains(".a"))
combined_12Y <- combined_12Y[c("DP31_Y12M",#food bank
                               "DP52_Y12M", #borrow money == 3
                               "FIN57_5_Y12M",#job seeker
                               "FIN57_6_Y12M" ,#domestic purposes/sole parent support 
                               "FIN57_7_Y12M" ,#sickness
                               "FIN57_17_Y12M" ,#dsiability
                               "FIN57_11_Y12M" , #student allowance
                               "FIN57_15_Y12M",#family tax credits
                               "FIN57_18_Y12M", #OSCAr (low or middle income)
                               "FIN57_10_Y12M", #accomodation supplement
                               "FIN57_19_Y12M", #training incentive
                               "FIN57_20_Y12M", #public housing/rent subsiy
                               "sum" )]
table(combined_12Y$DP52_Y12M, exclude = NULL)

combined_12Y$foodbank[combined_12Y$DP31_Y12M == 0] <- 0
combined_12Y$foodbank[combined_12Y$DP31_Y12M >= 1] <- 1
combined_12Y$foodbank[combined_12Y$DP31_Y12M == 99] <- NA
combined_12Y$borrow[combined_12Y$DP52_Y12M == 1] <- 0 #Lower
combined_12Y$borrow[combined_12Y$DP52_Y12M >= 2 ] <- 1
combined_12Y$borrow[combined_12Y$DP52_Y12M == 99] <- NA

combined_12Y$benefit[combined_12Y$FIN57_5_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_6_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_7_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_17_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_11_Y12M == 0] <- 0

combined_12Y$benefit[combined_12Y$FIN57_18_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_10_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_19_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_20_Y12M == 0] <- 0

combined_12Y$benefit[combined_12Y$FIN57_5_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_6_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_7_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_17_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_11_Y12M == 1] <- 1

combined_12Y$benefit[combined_12Y$FIN57_18_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_10_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_19_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_20_Y12M == 1] <- 1


combined_12Y$famtax[combined_12Y$FIN57_15_Y12M == 0] <- 0
combined_12Y$famtax[combined_12Y$FIN57_15_Y12M == 1] <- 1


test <- combined_12Y[c("sum", "benefit", "foodbank", "borrow", "famtax")]
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
combined_12Y <- dplyr::left_join(DEP_12Y, dcw12_mother, by = "FAMID", suffix = c("", ".a"))
combined_12Y <- combined_12Y %>% dplyr::select(-contains(".a"))

combined_12Y <- combined_12Y[c("DP31_Y12M",#food bank
                               "DP52_Y12M", #borrow money == 3
                               "FIN57_5_Y12M",#job seeker
                               "FIN57_6_Y12M" ,#domestic purposes/sole parent support 
                               "FIN57_7_Y12M" ,#sickness
                               "FIN57_17_Y12M" ,#dsiability
                               "FIN57_11_Y12M" , #student allowance
                               "FIN57_15_Y12M",#family tax credits
                               "FIN57_18_Y12M", #OSCAr (low or middle income)
                               "FIN57_10_Y12M", #accomodation supplement
                               "FIN57_19_Y12M", #training incentive
                               "FIN57_20_Y12M", #public housing/rent subsiy
                               "sum" )]
table(combined_12Y$DP52_Y12M, exclude = NULL)

combined_12Y$foodbank[combined_12Y$DP31_Y12M == 0] <- 0
combined_12Y$foodbank[combined_12Y$DP31_Y12M >= 1] <- 1
combined_12Y$foodbank[combined_12Y$DP31_Y12M == 99] <- NA
combined_12Y$borrow[combined_12Y$DP52_Y12M <= 2] <- 0
combined_12Y$borrow[combined_12Y$DP52_Y12M == 3 ] <- 1
combined_12Y$borrow[combined_12Y$DP52_Y12M == 99] <- NA

combined_12Y$benefit[combined_12Y$FIN57_5_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_6_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_7_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_17_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_11_Y12M == 0] <- 0

combined_12Y$benefit[combined_12Y$FIN57_18_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_10_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_19_Y12M == 0] <- 0
combined_12Y$benefit[combined_12Y$FIN57_20_Y12M == 0] <- 0

combined_12Y$benefit[combined_12Y$FIN57_5_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_6_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_7_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_17_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_11_Y12M == 1] <- 1

combined_12Y$benefit[combined_12Y$FIN57_18_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_10_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_19_Y12M == 1] <- 1
combined_12Y$benefit[combined_12Y$FIN57_20_Y12M == 1] <- 1


combined_12Y$famtax[combined_12Y$FIN57_15_Y12M == 0] <- 0
combined_12Y$famtax[combined_12Y$FIN57_15_Y12M == 1] <- 1


test <- combined_12Y[c("sum", "benefit", "foodbank", "borrow", "famtax")]
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
table(DEP_12Y$sum)
DEP_12Y$threshold12[DEP_12Y$sum >=2] <- 1
DEP_12Y$threshold12[DEP_12Y$sum < 2] <- 0
table(DEP_12Y$threshold12, exclude = NULL)
prop.table(table(DEP_12Y$threshold12))*100


#### Comparison of all three methods ####
DEP_12Y <- dplyr::left_join(DEP_12Y, LCA_12Y, by = "FAMID", suffix = c("", ".a"))
DEP_12Y <- DEP_12Y %>% dplyr::select(-contains(".a"))

table(DEP_12Y$LCA_2)
table(DEP_12Y$anyhardship)
table(DEP_12Y$threshold12)
DEP_12Y$LCA_12Y_revised[DEP_12Y$LCA_2 == 1] <- 0
DEP_12Y$LCA_12Y_revised[DEP_12Y$LCA_2 == 2] <- 1
table(DEP_12Y$LCA_12Y_revised)

table(DEP_12Y$threshold12 , DEP_12Y$LCA_12Y_revised)
table(DEP_12Y$threshold12 , DEP_12Y$anyhardship)
table(DEP_12Y$LCA_12Y_revised , DEP_12Y$anyhardship)

DEP_12Y$congruence <- ifelse(DEP_12Y$anyhardship == 1 & DEP_12Y$threshold12 == 1 & DEP_12Y$LCA_12Y_revised == 1, 1, 
                            ifelse(DEP_12Y$anyhardship == 0 & DEP_12Y$threshold12 == 0 & DEP_12Y$LCA_12Y_revised == 0, 0, 2))
table(DEP_12Y$congruence)

