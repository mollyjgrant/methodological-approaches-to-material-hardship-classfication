#### Libraries ####
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
library(rstatix)

##### Read in datasets ##### 
month9 <- fullbaseline

DEP_9M <- month9[c("FAMID", "dp1_m9m", "dp2_m9m" , "dp4_m9m","dp5_m9m","dp6_m9m")] #Note: not using dp3_m9m as will use this for sensitivity

DEP_9M <- DEP_9M[!apply(DEP_9M == 99, 1, any), ]
DEP_9M <- DEP_9M[!apply(DEP_9M == 98, 1, any), ]

DEP_9M$sum <- rowSums(DEP_9M[c( "dp1_m9m", "dp2_m9m" , "dp4_m9m","dp5_m9m","dp6_m9m")])
table(DEP_9M$sum)

DEP_9M <- na.omit(DEP_9M)

DEP_9M$M9 <- 1

#### Missing ####
View(miss_case_summary(DEP_9M[2:6]))
View(miss_var_summary(DEP_9M[2:6])) 
#No missing data

#### #### Item descriptives for 9m ####
DEP_9M <- DEP_9M[c("FAMID", "dp1_m9m", "dp2_m9m" , "dp4_m9m","dp5_m9m","dp6_m9m")]

psych::describe(DEP_9M[,2:6])
DEP_9M <- apply(DEP_9M, 2, as.numeric) # convert to numeric
DEP_9M <- as.data.frame(DEP_9M)

table(DEP_9M$dp1_m9m, exclude = NULL)
table(DEP_9M$dp2_m9m, exclude = NULL)
table(DEP_9M$dp4_m9m, exclude = NULL)
table(DEP_9M$dp5_m9m, exclude = NULL)
table(DEP_9M$dp6_m9m, exclude = NULL)

DEP_9M$sum <- rowSums(DEP_9M[, 2:6])

# Calculate tetrachoric correlations
corr_matrix <- tetrachoric(DEP_9M[,2:6], y = NULL, na.rm = TRUE)

print(corr_matrix)
heatmap(as.matrix(corr_matrix$rho), main = "Tetrachoric Correlation Heatmap", xlab = "", ylab = "", 
        Colv = NA, Rowv = NA, col = rev(heat.colors(12)), margins = c(5, 10))
legend("right", title = "Correlation", legend = round(seq(-1, 1, length.out = 12), 1), 
       fill = rev(heat.colors(12)), cex = 0.8, y.intersp = 1.2, bty = "n")

#### Cronbach for 9m ####
DEP_9M <- apply(DEP_9M, 2, as.numeric)
psych::omega(DEP_9M[,2:6])

######### Factorability #########
#KMO
KMO(corr_matrix$rho)
KMO(DEP_9M)

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
fact <- DEP_9M[,2:6]
par(mfrow=c(1,2), mar=c(5,5,2,2))
factnumber <- psych::fa.parallel(DEP_9M[,2:6], fm = 'pa', fa = 'fa', cor = "tet")
factnumber <- psych::fa.parallel(corr_matrix$rho, fm = 'pa', fa = 'fa', cor = "tet", n.obs = 6368)

print(factnumber)
set.seed(123)
Loadings <- fa(DEP_9M[,2:6], nfactors = 1, rotate = 'oblimin', fm = 'pa', cor = "tet")
print(Loadings$loadings,cutoff = 0.3)

#### Measurement invariance for 8 ####
names(month9)[names(month9) == "id_m9m"] <- "CID" 
names(DEP_9M)[names(DEP_9M) == "id_m9m"] <- "CID" 

DEP_9M <- left_join(DEP_9M , month9, by = "CID", suffix = c("", ".x"))
DEP_9M <- DEP_9M %>% dplyr::select(-contains(".x"))

DEP_9M <- left_join(DEP_9M , mother_AN, by = "MID", suffix = c("", ".x"))
DEP_9M <- DEP_9M %>% dplyr::select(-contains(".x"))

DEP_9M_mi <- DEP_9M[c("FAMID", "MID", "dp1_m9m", "dp2_m9m" ,"dp4_m9m","dp5_m9m","dp6_m9m",
                      "FIN_9M", "hh6_m9m", "pq5_m9m", "EDALL_AM", "AGE_GROUP_AM", "nzdep2006_m9m", "ETH5_NZDER_AM",
                      "ETH5_E_AM", "ETH5_O_AM","ETH5_MELA_AM", "ETH5_P_AM", "ETH5_A_AM", "ETH5_M_AM")]

###Ethnicity
library(tidyverse)
library(semTools)
DEP_9M_mi$Ethnicity_AN[DEP_9M_mi$ETH5_NZDER_AM ==1] <- "European"
DEP_9M_mi$Ethnicity_AN[DEP_9M_mi$ETH5_E_AM==1] <- "European"
DEP_9M_mi$Ethnicity_AN[DEP_9M_mi$ETH5_O_AM==1] <- "Other"
DEP_9M_mi$Ethnicity_AN[DEP_9M_mi$ETH5_MELA_AM==1] <- "Other"
DEP_9M_mi$Ethnicity_AN[DEP_9M_mi$ETH5_P_AM==1] <- "Pacific"
DEP_9M_mi$Ethnicity_AN[DEP_9M_mi$ETH5_A_AM==1] <- "Asian"
DEP_9M_mi$Ethnicity_AN[DEP_9M_mi$ETH5_M_AM==1] <- "M\u101ori"

measurementInvariance(model = model, data=DEP_9M_mi, 
                      group = "Ethnicity_AN")

###Household tenure
table(DEP_9M_mi$hh6_m9m, exclude = NULL)
DEP_9M_mi$hh6_m9m <- recode(DEP_9M_mi$hh6_m9m, "c(99) = NA; c(98) = NA")
DEP_9M_mi$hh6_m9m <- as.factor(DEP_9M_mi$hh6_m9m)

# model comparison tests
measurementInvariance(model = model, data=DEP_9M_mi, 
                      group = "hh6_m9m") #Configural = Configural; Loadings = Weak; Intercepts = Strong; Means = Strict

### Relationship status
table(DEP_9M_mi$pq5_m9m) 
DEP_9M_mi$pq5_m9m <- recode(DEP_9M_mi$pq5_m9m, "c(98) = NA")
DEP_9M_mi$pq5_m9m <- as.factor(DEP_9M_mi$pq5_m9m)

measurementInvariance(model = model, data=DEP_9M_mi, 
                      group = "pq5_m9m")

##Factor scores
DEP_9M_mi$refined.scores <- lavPredict(fit, newdata=DEP_9M_mi, type= "lv", method = "regression" )
DEP_9M_mi$sum.scores <- rowSums(DEP_9M_mi[c("dp1_m9m", "dp2_m9m" ,"dp4_m9m","dp5_m9m","dp6_m9m")])
sum.scores <- DEP_9M_mi$sum.scores

table(sum.scores, exclude = NULL)
View(table(DEP_9M_mi$refined.scores, exclude = NULL))

summary(lm(refined.scores ~ sum.scores, drop.unused.levels = TRUE))
plot(refined.scores, sum.scores)
cor.test(refined.scores, sum.scores, method = c("pearson"))



####Latent Class####
library(poLCA)
LCA_9M <- DEP_9M[c("dp1_m9m", "dp2_m9m" ,"dp4_m9m","dp5_m9m","dp6_m9m","FAMID")]
LCA_9M <- na.omit(LCA_9M)
LCA_9M[,1:5][LCA_9M[,1:5] >= 1] <- 2
LCA_9M[,1:5][LCA_9M[,1:5] < 1] <- 1
#write.csv(LCA_9M, "LCA_9M.csv")
View(miss_case_summary(DEP_9M[,2:6]))

f <- cbind(dp1_m9m, dp2_m9m ,dp4_m9m,dp5_m9m,dp6_m9m)~1

# Set up storage for AIC, BIC , adjusted BIC values
AIC <- rep(NA, 5)
BIC <- rep(NA, 5)
ADJBIC <- rep(NA, 5)
AWE <- rep(NA, 5)

# Fit LCA models with 1 to 10 classes and extract AIC and BIC values
set.seed(123)
for (i in 1:5) {
  fit <- poLCA(f, data = LCA_9M[,1:5], nclass = i, maxiter = 1000)
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
set.seed(123)
fit_3 <- poLCA(f, data = LCA_9M[,1:5], nclass = 3, maxiter = 10000, graph = T)
table(fit_3$predclass)
fit_3$posterior

LCA_9M <- cbind(LCA_9M, "LCA_9M" = fit_3$predclass)
table(LCA_9M$LCA_9M, exclude = NULL)

LCA_9M$dp1_m9m <- recode(LCA_9M$dp1_m9m, "c(1) = 0; c(2) = 1")
LCA_9M$dp2_m9m <- recode(LCA_9M$dp2_m9m, "c(1) = 0; c(2) = 1")
LCA_9M$dp4_m9m <- recode(LCA_9M$dp4_m9m, "c(1) = 0; c(2) = 1")
LCA_9M$dp5_m9m <- recode(LCA_9M$dp5_m9m, "c(1) = 0; c(2) = 1")
LCA_9M$dp6_m9m <- recode(LCA_9M$dp6_m9m, "c(1) = 0; c(2) = 1")

LCA_9M$sum <- rowSums((LCA_9M[,1:5]))

table(LCA_9M$sum, LCA_9M$LCA_3, exclude = NULL)
table(LCA_9M$LCA_3, exclude = NULL)

prop.table(table(LCA_9M$sum))*100

plot <- melt(fit_3$probs)
plot

plot$test[plot$Var2 == "Pr(1)"] <- "no hardship"
plot$test[plot$Var2 == "Pr(2)"] <- "hardship"
plot$L1[plot$L1 == "dp5_m9m"] <- "Gone without fresh fruit and vegetables"
plot$L1[plot$L1 == "dp1_m9m"] <- "Forced to buy cheaper food "
plot$L1[plot$L1 == "dp2_m9m"] <- "Put up with feeling cold"
plot$L1[plot$L1 == "dp4_m9m"] <- "Continued wearing shoes with holes"
plot$L1[plot$L1 == "dp6_m9m"] <- "Received help in the form of food, clothes, or money"

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
DEP_9M$anyhardship[DEP_9M$dp1_m9m < 0.5] <- 0
DEP_9M$anyhardship[DEP_9M$dp2_m9m < 0.5] <- 0
DEP_9M$anyhardship[DEP_9M$dp4_m9m < 0.5] <- 0
DEP_9M$anyhardship[DEP_9M$dp5_m9m < 0.5] <- 0
DEP_9M$anyhardship[DEP_9M$dp6_m9m < 0.5] <- 0
DEP_9M$anyhardship[DEP_9M$dp1_m9m >= 0.5] <- 1
DEP_9M$anyhardship[DEP_9M$dp2_m9m >= 0.5] <- 1
DEP_9M$anyhardship[DEP_9M$dp4_m9m >= 0.5] <- 1
DEP_9M$anyhardship[DEP_9M$dp5_m9m >= 0.5] <- 1
DEP_9M$anyhardship[DEP_9M$dp6_m9m >= 0.5] <- 1
table(DEP_9M$anyhardship)


################################################ Cut points based on cross-validation ####

DEP_9M <- as.data.frame(DEP_9M)
combined_9M <- dplyr::left_join(DEP_9M, month9, by = "FAMID", suffix = c("", ".a"))
combined_9M <- combined_9M %>% dplyr::select(-contains(".a"))

combined_9M$sum <- rowSums(combined_9M[, 2:6])
table(combined_9M$sum, exclude = NULL)
combined_9M <- combined_9M[c("dp3_m9m", #food bank
                             "fs4_m9m", #money worry
                             "nfin17_7_m9m", #unemployment benefit
                             "nfin17_9_m9m" , #domestic purposes/sole parent support 
                             "nfin17_8_m9m" , #sickness
                             "nfin17_10_m9m" , #invalid
                             "nfin17_11_m9m" , #student allowance
                             "nfin17_12_m9m", #Other govt support
                             "nfin17_16_m9m", #family tax credit/WFF
                             "sum" #sum score
                             )]
combined_9M <- na.omit(combined_9M)
combined_9M <- combined_9M[!apply(combined_9M == 99, 1, any), ]
combined_9M <- combined_9M[!apply(combined_9M == 98, 1, any), ]

combined_9M$benefit[combined_9M$nfin17_7_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_7_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_7_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_9_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_9_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_9_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_8_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_8_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_8_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_10_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_10_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_10_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_11_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_11_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_11_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_12_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_12_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_12_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_7_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_9_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_8_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_10_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_11_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_12_m9m == 1] <- 1

combined_9M$famtax[combined_9M$nfin17_16_m9m == 98] <- NA
combined_9M$famtax[combined_9M$nfin17_16_m9m == 99] <- NA
combined_9M$famtax[combined_9M$nfin17_16_m9m == 0] <- 0
combined_9M$famtax[combined_9M$nfin17_16_m9m == 1] <- 1

table(combined_9M$nfin17_9_m9m, combined_9M$famtax, exclude = NULL)

table(combined_9M$dp3_m9m, exclude = NULL) #Food bank usage
table(combined_9M$fs4_m9m, exclude = NULL) 

combined_9M$economic[combined_9M$fs4_m9m == 99] <- NA
combined_9M$economic[combined_9M$dp3_m9m == 99 ] <- NA
combined_9M$economic[combined_9M$fs4_m9m == 5] <- 0
combined_9M$economic[combined_9M$fs4_m9m  <= 3] <- 0
combined_9M$economic[combined_9M$dp3_m9m == 0 ] <- 0

combined_9M$economic[combined_9M$fs4_m9m  > 3] <- 1
combined_9M$economic[combined_9M$dp3_m9m == 1 ] <- 1

table(combined_9M$economic, exclude = NULL)

test <- combined_9M[c("sum", "benefit", "economic", "famtax")]
test <- na.omit(test)

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

### 9-month
DEP_9M$sum <- rowSums((DEP_9M[,2:6]))
prop.table(table(DEP_9M$sum, exclude = NULL))*100
table(DEP_9M$sum, exclude = NULL)

DEP_9M$threshold[DEP_9M$sum >= 1] <- 1
DEP_9M$threshold[DEP_9M$sum < 1] <- 0

prop.table(table(DEP_9M$threshold, exclude = NULL))*100
table(DEP_9M$threshold, exclude = NULL)

#Number in hardship at each threshold: 1 = 3525; 2 =	1563;	3 = 702; 4=302; 5 =	75
#prop in hardship at each threshold: 1 = 55.3549; 2=	24.544598; 3=	11.023869; 4=	4.742462; 5=	1.177764

#### Graph of six groups in hardship - 9 month ####
test$threshold[test$sum >= 1] <- 1
test$threshold[test$sum < 1] <- 0
table(test$threshold)
table(test$validation, test$threshold, exclude = NULL)

# Create a data frame with the values (divide prevalnce at each threshold by total sample size)
df <- data.frame(
  Category = c(">1", ">2", ">3", ">4", "5"),
  `Highest likelihood of hardship` = c(10.96105528,	7.851758794	,4.726758794,	2.653894472,	0.753768844	),
  `Higher likelihood of hardship` = c(5.873115578	,3.753140704,	2.182788945,	1.114949749,	0.314070352),
  `High likelihood of hardship` = c(4.302763819	,2.119974874,	0.769472362	,0.376884422	,0.06281407),
  `Low likelihood of hardship` = c(7.349246231	,3.046482412,	1.067839196	,0.141331658,	0.015703518),
  `Lower likelihood of hardship` = c(12.54711055	,4.616834171,	1.460427136,	0.282663317,	0.031407035),
  `Lowest likelihood of hardship` = c(14.32160804,	3.156407035,	0.816582915,	0.172738693,	0 )
)

names(df)[names(df) == "Lowest.likelihood.of.hardship"] <- "Lowest likelihood of hardship"
names(df)[names(df) == "Lower.likelihood.of.hardship"] <- "Lower likelihood of hardship"
names(df)[names(df) == "Low.likelihood.of.hardship"] <-  "Low likelihood of hardship"
names(df)[names(df) == "High.likelihood.of.hardship"] <-  "High likelihood of hardship"
names(df)[names(df) == "Higher.likelihood.of.hardship"] <-  "Higher likelihood of hardship"
names(df)[names(df) ==  "Highest.likelihood.of.hardship" ] <- "Highest likelihood of hardship"

# Reshape the data from wide to long format
df_long <- tidyr::pivot_longer(df, cols = -Category, names_to = "Group", values_to = "Value")

df_long$Group <- factor(df_long$Group, levels = c( "Lowest likelihood of hardship",
                                                   "Lower likelihood of hardship",
                                                   "Low likelihood of hardship",
                                                   "High likelihood of hardship",
                                                   "Higher likelihood of hardship",
                                                   "Highest likelihood of hardship"  
))


df <- data.frame(
  Category = c(">1", ">2", ">3", ">4", "5"),
  `High likelihood of hardship` = c(21.14	,13.72,	7.68,	4.15,	1.13	),
  `Low likelihood of hardship` = c(34.22,	10.82	,3.34,	0.60,	0.05)
)

names(df)[names(df) == "High.likelihood.of.hardship"] <- "High/higher/highest likelihood of hardship"
names(df)[names(df) == "Low.likelihood.of.hardship"] <- "Low/lower/lowest likelihood of hardship"


df_long <- tidyr::pivot_longer(df, cols = -Category, names_to = "Group", values_to = "Value")

df_long$Group <- factor(df_long$Group, levels = c( "Low/lower/lowest likelihood of hardship","High/higher/highest likelihood of hardship"
                                              
))

# Create the stacked bar graph
ggplot(df_long, aes(x = Category, y = Value , fill = Group)) +
  geom_bar(stat = "identity") +
  labs(x = "Material Hardship Threshold", y = "Percent of sample in hardship") +
  scale_fill_grey(start = 0.85, end = 0.45)  +
  geom_text(
    data = subset(df_long, Value > 1.2),  # Show labels only for values greater than 0.2
    aes(label = round(Value, 1)),
    position = position_stack(vjust = 0.70),
    family = "Calibri",  # Specify the desired font family
    size = 4,  # Specify the font size
    color = "black" , # Specify the font color
    fontface = "bold"
  )  +
  guides(fill = guide_legend(title = ""))+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(family = "Calibri", 
                            size = 12,
                            color = "black"))

### Second graph
# Create a data frame with the given values
df <- data.frame(
  Group = c("Highest likelihood of in hardship", "Higher likelihood of in hardship","High likelihood of in hardship",
            "Lower likelihood of in hardship", "Lowest likelihood of in hardship", "Low likelihood of in hardship"),
  `1 item to 2` = c(-3.11, -2.12, -2.18, 4.3, 7.93, 11.17),
  `2 items to 3` = c(-3.13	, -1.57, 	-1.35, 1.98, 3.16, 2.34),
  `3 items to 4` = c(-2.07	, -1.07	, -0/39, 0.93, 1.18, 0.64),
  `4 items to 5` = c(-1.90, -0.80, -0.31, 0.13, 0.25, 0.17))

names(df)[names(df) == "X1.item.to.2"] <- "1 item to 2"
names(df)[names(df) == "X2.items.to.3"] <- "2 items to 3"
names(df)[names(df) == "X3.items.to.4"] <- "3 items to 4"
names(df)[names(df) == "X4.items.to.5"] <- "4 items to 5"

df <- data.frame(
  Group = c("High/higher/highest likelihood of in hardship",
             "Low/lower/lowest likelihood of in hardship"),
  `1 item to 2` = c(-7.4, 23.4),
  `2 items to 3` = c(-6.1, 7.5),
  `3 items to 4` = c(-3.5, 2.8 ),
  `4 items to 5` = c(-3.0, 0.6))

names(df)[names(df) == "X1.item.to.2"] <- "1 item to 2"
names(df)[names(df) == "X2.items.to.3"] <- "2 items to 3"
names(df)[names(df) == "X3.items.to.4"] <- "3 items to 4"
names(df)[names(df) == "X4.items.to.5"] <- "4 items to 5"

# Convert the data frame to long format
df_long <- df %>% tidyr::pivot_longer(-Group, names_to = "Category", values_to = "Value")

df_long$Group <- factor(df_long$Group, levels = c(
                                                   "Low/lower/lowest likelihood of in hardship",
                                                   "High/higher/highest likelihood of in hardship"
                                                   
                                                   ))

# Create the side-by-side bar chart with grayscale color scheme
ggplot(df_long, aes(x = Category, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Threshold Change", y = "Percent of sample in hardship") +
  scale_fill_grey(start = 0.85, end = 0.45) +
  geom_text(data = subset(df_long, Value > -12),  # Show labels only for values greater than 0.2
            aes(label = round(Value, 1)),
            position = position_dodge(width = 0.9), vjust = -1.25,
            family = "Calibri",  # Specify the desired font family
            size = 3.25,  # Specify the font size
            color = "black" , # Specify the font color
            fontface = "bold"
  )  +
  theme_classic() +
  theme(legend.position = "bottom",
        text = element_text(family = "Calibri", 
                            size = 12,
                            color = "black"))+
  guides(fill = guide_legend(title = ""))

####Robustness test 1: Reduce the threshold of worrying about money to "Moderate Stress" and "High Stress" ####

DEP_9M <- as.data.frame(DEP_9M)
combined_9M <- dplyr::left_join(DEP_9M, month9, by = "FAMID", suffix = c("", ".a"))
combined_9M <- combined_9M %>% dplyr::select(-contains(".a"))

combined_9M$sum <- rowSums(combined_9M[, 2:6])

combined_9M <- combined_9M[c("dp3_m9m", #food bank
                             "fs4_m9m", #money worry
                             "nfin17_7_m9m", #unemployment benefit
                             "nfin17_9_m9m" , #domestic purposes/sole parent support 
                             "nfin17_8_m9m" , #sickness
                             "nfin17_10_m9m" , #invalid
                             "nfin17_11_m9m" , #student allowance
                             "nfin17_12_m9m", #Other govt support
                             "nfin17_16_m9m", #family tax credit/WFF
                             "sum" #sum score
)]
combined_9M <- na.omit(combined_9M)
combined_9M <- combined_9M[!apply(combined_9M == 99, 1, any), ]
combined_9M <- combined_9M[!apply(combined_9M == 98, 1, any), ]

combined_9M$benefit[combined_9M$nfin17_7_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_7_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_7_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_9_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_9_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_9_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_8_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_8_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_8_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_10_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_10_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_10_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_11_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_11_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_11_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_12_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_12_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_12_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_7_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_9_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_8_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_10_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_11_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_12_m9m == 1] <- 1

combined_9M$famtax[combined_9M$nfin17_16_m9m == 98] <- NA
combined_9M$famtax[combined_9M$nfin17_16_m9m == 99] <- NA
combined_9M$famtax[combined_9M$nfin17_16_m9m == 0] <- 0
combined_9M$famtax[combined_9M$nfin17_16_m9m == 1] <- 1

table(combined_9M$nfin17_9_m9m, combined_9M$famtax, exclude = NULL)

table(combined_9M$dp3_m9m, exclude = NULL) #Food bank usage
table(combined_9M$fs4_m9m, exclude = NULL) 

combined_9M$economic[combined_9M$fs4_m9m == 99] <- NA
combined_9M$economic[combined_9M$fs4_m9m == 5] <- 0
combined_9M$economic[combined_9M$fs4_m9m  <= 2] <- 0
combined_9M$economic[combined_9M$dp3_m9m == 0 ] <- 0
combined_9M$economic[combined_9M$dp3_m9m == 99 ] <- NA

combined_9M$economic[combined_9M$fs4_m9m  > 2] <- 1
combined_9M$economic[combined_9M$dp3_m9m == 1 ] <- 1

table(combined_9M$economic, exclude = NULL)

test <- combined_9M[c("sum", "benefit", "economic", "famtax")]
test <- na.omit(test)

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


####Robustness test 2:Separate all three aspects ####
DEP_9M <- as.data.frame(DEP_9M)
combined_9M <- dplyr::left_join(DEP_9M, month9, by = "FAMID", suffix = c("", ".a"))
combined_9M <- combined_9M %>% dplyr::select(-contains(".a"))

combined_9M$sum <- rowSums(combined_9M[, 2:6])

combined_9M <- combined_9M[c("dp3_m9m", #food bank
                             "fs4_m9m", #money worry
                             "nfin17_7_m9m", #unemployment benefit
                             "nfin17_9_m9m" , #domestic purposes/sole parent support 
                             "nfin17_8_m9m" , #sickness
                             "nfin17_10_m9m" , #invalid
                             "nfin17_11_m9m" , #student allowance
                             "nfin17_12_m9m", #Other govt support
                             "nfin17_16_m9m", #family tax credit/WFF
                             "sum" #sum score
)]
combined_9M <- na.omit(combined_9M)
combined_9M <- combined_9M[!apply(combined_9M == 99, 1, any), ]
combined_9M <- combined_9M[!apply(combined_9M == 98, 1, any), ]

combined_9M$benefit[combined_9M$nfin17_7_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_7_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_7_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_9_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_9_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_9_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_8_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_8_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_8_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_10_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_10_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_10_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_11_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_11_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_11_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_12_m9m == 98] <- NA
combined_9M$benefit[combined_9M$nfin17_12_m9m == 99] <- NA
combined_9M$benefit[combined_9M$nfin17_12_m9m == 0] <- 0

combined_9M$benefit[combined_9M$nfin17_7_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_9_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_8_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_10_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_11_m9m == 1] <- 1
combined_9M$benefit[combined_9M$nfin17_12_m9m == 1] <- 1

combined_9M$famtax[combined_9M$nfin17_16_m9m == 98] <- NA
combined_9M$famtax[combined_9M$nfin17_16_m9m == 99] <- NA
combined_9M$famtax[combined_9M$nfin17_16_m9m == 0] <- 0
combined_9M$famtax[combined_9M$nfin17_16_m9m == 1] <- 1

table(combined_9M$nfin17_9_m9m, combined_9M$famtax, exclude = NULL)

table(combined_9M$dp3_m9m, exclude = NULL) #Food bank usage
table(combined_9M$fs4_m9m, exclude = NULL) 

combined_9M$worry[combined_9M$fs4_m9m == 99] <- NA
combined_9M$worry[combined_9M$fs4_m9m == 5] <- 0
combined_9M$worry[combined_9M$fs4_m9m  <= 3] <- 0
combined_9M$worry[combined_9M$fs4_m9m  > 3] <- 1

combined_9M$food[combined_9M$dp3_m9m == 98 ] <- NA
combined_9M$food[combined_9M$dp3_m9m == 99 ] <- NA
combined_9M$food[combined_9M$dp3_m9m == 0 ] <- 0
combined_9M$food[combined_9M$dp3_m9m == 1 ] <- 1

table(combined_9M$food, exclude = NULL)

test <- combined_9M[c("sum", "benefit", "worry", "food", "famtax")]
test <- na.omit(test)

test$validation[test$benefit == 1 & test$worry == 1 & test$food == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 0 & test$worry == 1 & test$food == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 1 & test$worry == 1 & test$food == 1 & test$famtax == 0] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"

test$validation[test$benefit == 1 & test$worry == 0 & test$food == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 0 & test$worry == 0 & test$food == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 1 & test$worry == 0 & test$food == 1 & test$famtax == 0] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"

test$validation[test$benefit == 0 & test$worry == 1 & test$food == 1 & test$famtax == 0] <- "Either economic strain experiences & no benefit receipt"
test$validation[test$benefit == 0 & test$worry == 0 & test$food == 1 & test$famtax == 0] <- "Either economic strain experiences & no benefit receipt"
test$validation[test$benefit == 0 & test$worry == 1 & test$food == 0 & test$famtax == 0] <- "Either economic strain experiences & no benefit receipt"

test$validation[test$benefit == 1 & test$worry == 1 & test$food == 0 & test$famtax == 0] <- "Worry and Benefit receipt"
test$validation[test$benefit == 0 & test$worry == 1 & test$food == 0 & test$famtax == 1] <- "Worry and Benefit receipt"
test$validation[test$benefit == 1 & test$worry == 1 & test$food == 0 & test$famtax == 1] <- "Worry and Benefit receipt"

test$validation[test$benefit == 1 & test$worry == 0 & test$food == 0 & test$famtax == 1] <- "Benefit receipt only"
test$validation[test$benefit == 1 & test$worry == 0 & test$food == 0 & test$famtax == 0] <- "Benefit receipt only"
test$validation[test$benefit == 1 & test$worry == 0 & test$food == 0 & test$famtax == 1] <- "Benefit receipt only"

test$validation[test$benefit == 0 & test$worry == 0 & test$food == 0 & test$famtax == 1] <- "Fam tax only"
test$validation[test$benefit == 0 & test$worry == 0 & test$food == 0 & test$famtax == 0] <- "None"


test$validation <- factor(test$validation, levels = c("Food bank & benefit receipt / Both economic strain experiences and benefit receipt", 
                                                    "Either economic strain experiences & no benefit receipt",
                                                       "Worry and Benefit receipt", "Benefit receipt only",
                                                      "Fam tax only", "None"))

table(test$validation, exclude = NULL)
table(test$sum)
table(test$validation, test$sum, exclude = NULL)


### 9-month
DEP_9M$sum <- rowSums((DEP_9M[,2:6]))
prop.table(table(DEP_9M$sum, exclude = NULL))*100

#### Graph of four groups in hardship - 9 month ###
# Create a data frame with the values
df <- data.frame(
  Category = c(">1", ">2", ">3", ">4", "5"),
  `Highest likelihood of hardship` = c(3.93	,3.11,	2.12,	1.42,	0.42),
  `Higher likelihood of hardship` = c(13.23	,8.78	,4.95	,2.45	,0.66),
  `Lower likelihood of hardship` = c(23.91	,9.53	,3.15	,0.71	,0.09),
  `Lowest likelihood of hardship` = c(14.25	,3.15	,0.82	,0.17	,0.00)
)

names(df)[names(df) == "Lowest.likelihood.of.hardship"] <- "Lowest likelihood of hardship"
names(df)[names(df) == "Lower.likelihood.of.hardship"] <- "Lower likelihood of hardship"
names(df)[names(df) == "Higher.likelihood.of.hardship"] <-  "Higher likelihood of hardship"
names(df)[names(df) ==  "Highest.likelihood.of.hardship" ] <- "Highest likelihood of hardship"

# Reshape the data from wide to long format
df_long <- tidyr::pivot_longer(df, cols = -Category, names_to = "Group", values_to = "Value")

df_long$Group <- factor(df_long$Group, levels = c( "Lowest likelihood of hardship",
                                                   "Lower likelihood of hardship",
                                                   "Higher likelihood of hardship",
                                                   "Highest likelihood of hardship"  
))

# Create the stacked bar graph
ggplot(df_long, aes(x = Category, y = Value , fill = Group)) +
  geom_bar(stat = "identity") +
  labs(x = "Material Hardship Threshold", y = "% of sample in hardship") +
  scale_fill_grey(start = 0.85, end = 0.35)  +
  geom_text(
    data = subset(df_long, Value > 0.99),  # Show labels only for values greater than 0.2
    aes(label = round(Value, 1)),
    position = position_stack(vjust = 0.70),
    family = "Calibri",  # Specify the desired font family
    size = 3.25,  # Specify the font size
    color = "black" , # Specify the font color
    fontface = "bold"
  )  +
  guides(fill = guide_legend(title = ""))+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(family = "Calibri", 
                            size = 12,
                            color = "black"))

### Second graph
# Create a data frame with the given values
df <- data.frame(
  Group = c("Highest likelihood of in hardship", "Higher likelihood of in hardship",
            "Lower likelihood of in hardship", "Lowest likelihood of in hardship"),
  `1 item to 2` = c(-0.82, -4.45, 14.38, 11.10),
  `2 items to 3` = c(-0.99	, -3.82, 	6.39, 2.33),
  `3 items to 4` = c(-0.71	, -2.50	, 2.44, 0.64),
  `4 items to 5` = c(-0.99, -1.79, 0.61, 0.17))

names(df)[names(df) == "X1.item.to.2"] <- "1 item to 2"
names(df)[names(df) == "X2.items.to.3"] <- "2 items to 3"
names(df)[names(df) == "X3.items.to.4"] <- "3 items to 4"
names(df)[names(df) == "X4.items.to.5"] <- "4 items to 5"

# Convert the data frame to long format
df_long <- df %>% tidyr::pivot_longer(-Group, names_to = "Category", values_to = "Value")

df_long$Group <- factor(df_long$Group, levels = c( "Lowest likelihood of in hardship" ,
                                                   "Lower likelihood of in hardship", 
                                                   "Higher likelihood of in hardship",
                                                   "Highest likelihood of in hardship"
))

# Create the side-by-side bar chart with grayscale color scheme
ggplot(df_long, aes(x = Category, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Threshold Change", y = "% of sample in hardship") +
  scale_fill_grey(start = 0.85, end = 0.35) +
  geom_text(data = subset(df_long, Value > -12),  # Show labels only for values greater than 0.2
            aes(label = round(Value, 1)),
            position = position_dodge(width = 0.9), vjust = -2.25,
            family = "Calibri",  # Specify the desired font family
            size = 3.25,  # Specify the font size
            color = "black" , # Specify the font color
            fontface = "bold"
  )  +
  theme_classic() +
  theme(legend.position = "bottom",
        text = element_text(family = "Calibri", 
                            size = 12,
                            color = "black"))+
  guides(fill = guide_legend(title = ""))


### Applying the cut point ####
DEP_9M <- as.data.frame(DEP_9M)
DEP_9M$sum <- rowSums(DEP_9M[c("dp1_m9m", "dp2_m9m" , "dp4_m9m","dp5_m9m","dp6_m9m")])
table(DEP_9M$sum)
DEP_9M$threshold9[DEP_9M$sum >=3] <- 1
DEP_9M$threshold9[DEP_9M$sum < 3] <- 0
table(DEP_9M$threshold9, exclude = NULL)


#### Comparison of all three methods ####
DEP_9M <- dplyr::left_join(DEP_9M, LCA_9M, by = "FAMID", suffix = c("", ".a"))
DEP_9M <- DEP_9M %>% dplyr::select(-contains(".a"))

DEP_9M$LCA_9M_revised[DEP_9M$LCA_9M == 3] <- 0
DEP_9M$LCA_9M_revised[DEP_9M$LCA_9M == 1] <- 1
DEP_9M$LCA_9M_revised[DEP_9M$LCA_9M == 2] <- 1

table(DEP_9M$threshold9 , DEP_9M$LCA_9M_revised)
table(DEP_9M$threshold9 , DEP_9M$anyhardship)
table(DEP_9M$LCA_9M_revised , DEP_9M$anyhardship)

DEP_9M$congruence <- ifelse(DEP_9M$anyhardship == 1 & DEP_9M$threshold9 == 1 & DEP_9M$LCA_9M_revised == 1, 1, 
                            ifelse(DEP_9M$anyhardship == 0 & DEP_9M$threshold9 == 0 & DEP_9M$LCA_9M_revised == 0, 0, 2))
table(DEP_9M$congruence)
