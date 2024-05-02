
##### Read in datasets ##### 
month54 <- fullbaseline
DEP_54M <- month54[c("FAMID", "DP1_M54M", "DP2_M54M" ,"DP4_M54M","DP5_M54M","DP6_M54M")] #6074
DEP_54M <- DEP_54M[!apply(DEP_54M == 99, 1, any), ]

DEP_54M$sum <- rowSums(DEP_54M[c("DP1_M54M", "DP2_M54M" ,"DP4_M54M","DP5_M54M","DP6_M54M")])
table(DEP_54M$sum)

DEP_54M <- na.omit(DEP_54M)

DEP_54M$M54 <- 1

#### Missing ####
DEP_54M <- month54[c("FAMID", "DP1_M54M", "DP2_M54M" ,"DP4_M54M","DP5_M54M","DP6_M54M")] #6074

View(miss_case_summary(DEP_54M[2:6]))
View(miss_var_summary(DEP_54M[2:6]))
summary(miss_case_summary(DEP_54M[2:6]))
summary(miss_var_summary(DEP_54M[2:6]))

mcar_test(DEP_54M[2:6])

DEP_missing54 <- DEP_54M[which(rowMeans(!is.na(DEP_54M[2:6])) > 0.6), ] #With >40% missing removed.
DEP_missing54 <- DEP_missing54[!apply(DEP_missing54 == 99, 1, any), ]
DEP_missing54 <- DEP_missing54[!apply(DEP_missing54 == 98, 1, any), ]


data <- DEP_missing54

describe(DEP_54M[2:6])


#### Item descriptives for 9M #######
DEP_54M <- DEP_missing54[c("FAMID", "DP1_M54M", "DP2_M54M" , "DP4_M54M","DP5_M54M","DP6_M54M")] #Note: not using DP3_M54M as will use this for sensitivity
psych::describe(DEP_54M)

DEP_54M$sum <- rowSums(DEP_54M[c( "DP1_M54M", "DP2_M54M" , "DP4_M54M","DP5_M54M","DP6_M54M")])
table(DEP_54M$sum)


#### #### Item descriptives for 54m ####
psych::describe(DEP_54M[,2:6])
table(DEP_54M$DP1_M54M, exclude = NULL)
table(DEP_54M$DP6_M54M, exclude = NULL)

DEP_54M <- as.data.frame(DEP_54M) # convert to numeric

table(DEP_54M$DP6_M54M, exclude = NULL)

# Calculate tetrachoric correlations
corr_matrix <- tetrachoric(DEP_54M[,2:6], y = NULL,  na.rm = TRUE)
print(corr_matrix)

heatmap(as.matrix(corr_matrix$rho), main = "Tetrachoric Correlation Heatmap", xlab = "", ylab = "", 
        Colv = NA, Rowv = NA, col = rev(heat.colors(12)), margins = c(5, 10))
legend("right", title = "Correlation", legend = round(seq(-1, 1, length.out = 12), 1), 
       fill = rev(heat.colors(12)), cex = 0.8, y.intersp = 1.2, bty = "n")

#### Cronbach for 8Y ####
DEP_54M <- apply(DEP_54M, 2, as.numeric)
psych::alpha(DEP_54M[,2:6])

######### Factorability #########
#KMO
KMO(corr_matrix$rho)
KMO(DEP_54M)

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
factnumber <- psych::fa.parallel(DEP_54M[,2:6], fm = 'pa', fa = 'fa', cor = "tet")

print(factnumber)
set.seed(123)
Loadings <- fa(DEP_54M[,2:6], nfactors = 1, rotate = 'oblimin', fm = 'pa', cor = "tet")
print(Loadings$loadings,cutoff = 0.3)

###Ethnicity
library(tidyverse)
DEP_54M_mi$Ethnicity_AN[DEP_54M_mi$ETH5_NZDER_AM ==1] <- "European"
DEP_54M_mi$Ethnicity_AN[DEP_54M_mi$ETH5_E_AM==1] <- "European"
DEP_54M_mi$Ethnicity_AN[DEP_54M_mi$ETH5_O_AM==1] <- "Other"
DEP_54M_mi$Ethnicity_AN[DEP_54M_mi$ETH5_MELA_AM==1] <- "Other"
DEP_54M_mi$Ethnicity_AN[DEP_54M_mi$ETH5_P_AM==1] <- "Pacific"
DEP_54M_mi$Ethnicity_AN[DEP_54M_mi$ETH5_A_AM==1] <- "Asian"
DEP_54M_mi$Ethnicity_AN[DEP_54M_mi$ETH5_M_AM==1] <- "M\u101ori"

measurementInvariance(model = model, data=DEP_54M_mi, 
                      group = "Ethnicity_AN")

###Household tenure
table(DEP_54M_mi$HH6_M54M, exclude = NULL)
DEP_54M_mi$HH6_M54M <- recode(DEP_54M_mi$HH6_M54M, "c(99) = NA; c(98) = NA")
DEP_54M_mi$HH6_M54M <- as.factor(DEP_54M_mi$HH6_M54M)

# model comparison tests
measurementInvariance(model = model, data=DEP_54M_mi, 
                      group = "HH6_M54M") #Configural = Configural; Loadings = Weak; Intercepts = Strong; Means = Strict

### Relationship status
table(DEP_54M_mi$PQ5_M54M) 
DEP_54M_mi$PQ5_M54M <- recode(DEP_54M_mi$PQ5_M54M, "c(98) = NA")
DEP_54M_mi$PQ5_M54M <- as.factor(DEP_54M_mi$PQ5_M54M)

measurementInvariance(model = model, data=DEP_54M_mi, 
                      group = "PQ5_M54M")

##Factor scores
DEP_54M_mi$refined.scores <- lavPredict(fit, newdata=DEP_54M_mi, type= "lv", method = "regression" )
DEP_54M_mi$sum.scores <- rowSums(DEP_54M_mi[c("DP1_M54M", "DP2_M54M" ,"DP4_M54M","DP5_M54M","DP6_M54M")])
sum.scores <- DEP_54M_mi$sum.scores

table(sum.scores, exclude = NULL)
View(table(DEP_54M_mi$refined.scores, exclude = NULL))

summary(lm(refined.scores ~ sum.scores, drop.unused.levels = TRUE))
plot(refined.scores, sum.scores)
cor.test(refined.scores, sum.scores, method = c("pearson"))


#### Distributions of Sum Scores at each time point ####

DEP_54M_mi %>%
  ggplot( aes(x=sum.scores)) +
  theme(panel.background = element_rect(fill='transparent')) +
  geom_histogram(fill="#fdb913", color="#fdb913", alpha=1)+
  theme_ipsum()+
  theme_classic()+
  labs(x = "Hardship Score", y = "Count")+
  ggtitle("Hardship Scores at 54-month")

DEP_54M_mi %>%
  ggplot(aes(x=sum.scores)) +
  theme(panel.background = element_rect(fill='transparent')) +
  geom_density(fill="#fdb913", color="#fdb913", alpha=1, bw=0.4) + # Adjust the bandwidth value to control the smoothness of the curve
  theme_ipsum() +
  theme_classic() +
  labs(x = "Hardship Score", y = "Density") +
  ggtitle("Hardship Scores at 54-month")


####Create cut points each time point ####
DEP_54M_mi$sum <- rowSums((DEP_54M_mi[,2:6]))
table(DEP_54M_mi$sum, exclude = NULL)

DEP_54M_mi$zscore <- scale(DEP_54M_mi$sum)
table(DEP_54M_mi$zscore)

#Severe = 2.461894140
#Material = 1.42672144

DEP_54M_mi$categorical[DEP_54M_mi$zscore >= 1.42672144] <- "material hardship"
DEP_54M_mi$categorical[DEP_54M_mi$zscore >= 2.461894140] <- "severe material hardship"
DEP_54M_mi$categorical[DEP_54M_mi$zscore < 1.42672144] <- "no/little material hardship"
prop.table(table(DEP_54M_mi$categorical))
table(DEP_54M_mi$categorical, exclude = NULL) 

####Latent Class####
LCA_54M <- DEP_54M[c("DP1_M54M", "DP2_M54M" ,"DP4_M54M","DP5_M54M","DP6_M54M", "FAMID")]
View(miss_case_summary(LCA_54M[,1:5]))

LCA_54M[,1:5][LCA_54M[,1:5] >= 1] <- 2
LCA_54M[,1:5][LCA_54M[,1:5] < 1] <- 1

LCA_54M <- na.omit(LCA_54M)
DEP_54M <- na.omit(DEP_54M)
f <- cbind(DP1_M54M, DP2_M54M ,DP4_M54M,DP5_M54M,DP6_M54M)~1

# Set up storage for AIC, BIC , adjusted BIC values
AIC <- rep(NA, 5)
BIC <- rep(NA, 5)
ADJBIC <- rep(NA, 5)
AWE <- rep(NA, 5)

# Fit LCA models with 1 to 10 classes and extract AIC and BIC values
set.seed(123)
for (i in 1:5) {
  fit <- poLCA(f, data = LCA_54M, nclass = i, maxiter = 1000)
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
fit_3 <- poLCA(f, data = LCA_54M[,1:5], nclass = 3, maxiter = 10000)
table(fit_3$predclass)

LCA_54M <- cbind(LCA_54M, "LCA_3" = fit_3$predclass)
table(LCA_54M$LCA_3, exclude = NULL)

LCA_54M$DP1_M54M <- recode(LCA_54M$DP1_M54M, "c(1) = 0; c(2) = 1")
LCA_54M$DP2_M54M <- recode(LCA_54M$DP2_M54M, "c(1) = 0; c(2) = 1")
LCA_54M$DP4_M54M <- recode(LCA_54M$DP4_M54M, "c(1) = 0; c(2) = 1")
LCA_54M$DP5_M54M <- recode(LCA_54M$DP5_M54M, "c(1) = 0; c(2) = 1")
LCA_54M$DP6_M54M <- recode(LCA_54M$DP6_M54M, "c(1) = 0; c(2) = 1")

LCA_54M$sum <- rowSums(LCA_54M[,1:5])

table(LCA_54M$sum, LCA_54M$LCA_3, exclude = NULL)
table(LCA_54M$LCA_3)


prop.table(table(LCA_54M$sum))*100

test <- LCA_54M[c("")]

plot <- melt(fit_3$probs)
plot

plot$test[plot$Var2 == "Pr(1)"] <- "no hardship"
plot$test[plot$Var2 == "Pr(2)"] <- "hardship"
plot$L1[plot$L1 == "DP5_M54M"] <- "Gone without fresh fruit and vegetables"
plot$L1[plot$L1 == "DP1_M54M"] <- "Forced to buy cheaper food "
plot$L1[plot$L1 == "DP2_M54M"] <- "Put up with feeling cold"
plot$L1[plot$L1 == "DP4_M54M"] <- "Continued wearing shoes with holes"
plot$L1[plot$L1 == "DP6_M54M"] <- "Received help in the form of food, clothes, or money"

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
any_54M <- DEP_54M[c("DP1_M54M", "DP2_M54M" ,"DP4_M54M","DP5_M54M","DP6_M54M", "FAMID")]
any_54M <- na.omit(any_54M)

DEP_54M$anyhardship[DEP_54M$DP1_M54M < 0.5] <- 0
DEP_54M$anyhardship[DEP_54M$DP2_M54M < 0.5] <- 0
DEP_54M$anyhardship[DEP_54M$DP4_M54M < 0.5] <- 0
DEP_54M$anyhardship[DEP_54M$DP5_M54M < 0.5] <- 0
DEP_54M$anyhardship[DEP_54M$DP6_M54M < 0.5] <- 0
DEP_54M$anyhardship[DEP_54M$DP1_M54M >= 0.5] <- 1
DEP_54M$anyhardship[DEP_54M$DP2_M54M >= 0.5] <- 1
DEP_54M$anyhardship[DEP_54M$DP4_M54M >= 0.5] <- 1
DEP_54M$anyhardship[DEP_54M$DP5_M54M >= 0.5] <- 1
DEP_54M$anyhardship[DEP_54M$DP6_M54M >= 0.5] <- 1

table(DEP_54M$anyhardship)
prop.table(table(DEP_54M$anyhardship))

### Stack ###

LCA_54M$LCA_54M
LCA_54M$LCA_54M[LCA_54M$LCA_54M == 1] <- "Low hardship"
LCA_54M$LCA_54M[LCA_54M$LCA_54M == 2] <- "Moderate hardship"
LCA_54M$LCA_54M[LCA_54M$LCA_54M == 3] <- "High hardship"
names(LCA_54M)[names(LCA_54M) == "FAMID"] <- "MID" 

LCA_2Y <- month24[c("MID", "ls3_y2m1")]
LCA_2Y$ls3_y2m1[LCA_2Y$ls3_y2m1 == 0] <- "No hardship"
LCA_2Y$ls3_y2m1[LCA_2Y$ls3_y2m1 == 1] <- "Low hardship"
LCA_2Y$ls3_y2m1[LCA_2Y$ls3_y2m1 == 2] <- "Moderate hardship"
LCA_2Y$ls3_y2m1[LCA_2Y$ls3_y2m1 == 3] <- "High hardship"


LCA_54M$LCA_54
LCA_54M$LCA_54[LCA_54M$LCA_54 == 1] <- "No—Low hardship"
LCA_54M$LCA_54[LCA_54M$LCA_54 == 2] <- "Moderate hardship"
LCA_54M$LCA_54[LCA_54M$LCA_54 == 3] <- "Moderate—High hardship"


LCA_8Y$LCA_8
LCA_8Y$LCA_8[LCA_8Y$LCA_8 == 1] <- "Low—moderate hardship"
LCA_8Y$LCA_8[LCA_8Y$LCA_8 == 2] <- "Moderate—high hardship"
LCA_8Y$LCA_8[LCA_8Y$LCA_8 == 3] <- "No—Low hardship"
LCA_8Y$LCA_8[LCA_8Y$LCA_8 == 4] <- "Low—moderate hardship"
LCA_8Y$LCA_8[LCA_8Y$LCA_8 == 5] <- "High hardship"


combined_LCA <- dplyr::left_join(LCA_54M, LCA_2Y, by = "MID", suffix = c("", ".x"))
combined_LCA <- combined_LCA %>% dplyr::select(-contains(".x"))

combined_LCA <- dplyr::left_join(combined_LCA, LCA_54M, by = "MID", suffix = c("", ".x"))
combined_LCA <- combined_LCA %>% dplyr::select(-contains(".x"))

combined_LCA <- dplyr::left_join(combined_LCA, LCA_8Y, by = "MID", suffix = c("", ".x"))
combined_LCA <- combined_LCA %>% dplyr::select(-contains(".x"))

View(combined_LCA)

LCA <- combined_LCA[c("MID", "LCA_54M", "ls3_y2m1", "LCA_54", "LCA_8")]

LCA$LCA_54M <- factor(LCA$LCA_54M, levels = c("Low hardship", "Moderate hardship", "High hardship"))
LCA$ls3_y2m1 <- factor(LCA$ls3_y2m1, levels = c("No hardship", "Low hardship", "Moderate hardship", "High hardship"))
LCA$LCA_54 <- factor(LCA$LCA_54, levels = c("No—Low hardship", "Moderate hardship", "Moderate—High hardship"))
LCA$LCA_8 <- factor(LCA$LCA_8, levels = c("No—Low hardship", "Low—moderate hardship", "Moderate—high hardship", "High hardship"))


table(LCA$LCA_54M, exclude = NULL)
table(LCA$ls3_y2m1, exclude = NULL)
table(LCA$LCA_54, exclude = NULL)
table(LCA$LCA_8, exclude = NULL)

longitudinal_sequence <- LCA[c("MID","LCA_54M", "ls3_y2m1", "LCA_54", "LCA_8")]
View(longitudinal_sequence)
seqstatl(longitudinal_sequence[,2:5]) #This gives you the state names from your dataframe

dta.alphab <- c("High hardship"      ,    "Low—moderate hardship" , "Low hardship"  ,
                "Moderate—high hardship", "Moderate—High hardship" ,"Moderate hardship"  ,   
                "No—Low hardship" ,       "No hardship" )

seq <- seqdef(longitudinal_sequence, var = 2:5, xtstep = 1, alphabet = dta.alphab, cnames = c("LCA_54M", "ls3_y2m1", "LCA_54", "LCA_8"),
              tick.last = TRUE)


summary(seq) #200 unique sequences

dist.mostfreq<- seqdist(seq, method = "OM", sm = "TRATE", with.missing = T) #These two lines show all possible sequences (a little easier to interpret)
seqIplot(seq, border = NA,  sortv = dist.mostfreq)

### Cut points ####
month45 <- read_dta("DCW4/DCW4_m45M/STATA/DCW4M.DTA")
month24 <-  read_dta("DCW2/DCW2_y2M/STATA/DCW2Y2M.DTA")
names(month45)[names(month45) == "famid"] <- "FAMID" 
names(month24)[names(month24) == "famid"] <- "FAMID" 

DEP_54M$sum <- rowSums((DEP_54M[,2:6]))

combined_54 <- dplyr::left_join(DEP_54M, month54, by = "FAMID", suffix = c("", ".a"))
combined_54 <- combined_54 %>% dplyr::select(-contains(".a"))
combined_54 <- dplyr::left_join(combined_54, month45, by = "FAMID", suffix = c("", ".a"))
combined_54 <- combined_54 %>% dplyr::select(-contains(".a"))
combined_54 <- dplyr::left_join(combined_54, month24, by = "FAMID", suffix = c("", ".a"))
combined_54 <- combined_54 %>% dplyr::select(-contains(".a"))

combined_54 <- combined_54[c("DP3_M54M",#food bank
                             "ls3_y2m", #income meet daily needs
                             "finhhp_7_m45m",#unemployment benefit
                             "finhhp_9_m45m" ,#domestic purposes/sole parent support 
                             "finhhp_8_m45m" ,#sickness
                             "finhhp_10_m45m" ,#invalid
                             "finhhp_11_m45m" , #student allowance
                            
                             "finhhp_18_m45m",#family tax credits
                             "finhhp_16_m45m", #accomodation supplement
                             "finhhp_12_m45m",#Other govt support
                             "finhhp_s4_m45m", #Other government benefits
                             "sum" )]

#"finhhp_s3_m45m", #income tested benefit
table(combined_54$DP3_M54M,combined_54$ls3_y2m, exclude = NULL )
combined_54$foodbank[combined_54$DP3_M54M == 0] <- 0
combined_54$notenough[combined_54$ls3_y2m == 99 ] <- NA
combined_54$notenough[combined_54$ls3_y2m == 98 ] <- NA
combined_54$notenough[combined_54$ls3_y2m >= 2 ] <- 0
combined_54$foodbank[combined_54$DP3_M54M == 1] <- 1
combined_54$notenough[combined_54$ls3_y2m == 1 ] <- 1

combined_54$benefit[combined_54$finhhp_7_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_7_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_7_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_9_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_9_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_9_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_8_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_8_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_8_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_10_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_10_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_10_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_11_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_11_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_11_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_12_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_12_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_12_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_16_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_16_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_16_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_s4_m45m == 0] <- 0

combined_54$added[combined_54$finhhp_7_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_9_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_8_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_10_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_11_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_16_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_18_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_7_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_9_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_8_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_10_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_11_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_16_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_18_m45m == 1] <- 2

combined_54$benefit[combined_54$finhhp_s4_m45m == 1 & combined_54$added == 0] <- 1

combined_54$benefit[combined_54$finhhp_7_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_9_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_8_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_10_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_11_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_16_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_12_m45m == 1] <- 1

combined_54$famtax[combined_54$finhhp_18_m45m == 98] <- NA
combined_54$famtax[combined_54$finhhp_18_m45m == 99] <- NA
combined_54$famtax[combined_54$finhhp_18_m45m == 0] <- 0
combined_54$famtax[combined_54$finhhp_18_m45m == 1] <- 1

test <- combined_54[c("sum", "benefit", "foodbank","notenough", "famtax")]

test <- na.omit(test)
test$economic[test$foodbank == 0] <- 0
test$economic[test$notenough == 0] <- 0
test$economic[test$foodbank == 1] <- 1
test$economic[test$notenough == 1] <- 1


table(test$sum)
test$sum[test$sum < 0.5] <- 0
test$sum[test$sum >= 0.5 & test$sum < 1.5] <- 1
test$sum[test$sum >= 1.5 & test$sum < 2.5] <- 2
test$sum[test$sum >= 2.5 & test$sum < 3.5] <- 3

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


### Robustness check (all four) ####
combined_54 <- dplyr::left_join(DEP_54M, month54, by = "FAMID", suffix = c("", ".a"))
combined_54 <- combined_54 %>% dplyr::select(-contains(".a"))
combined_54 <- dplyr::left_join(combined_54, month45, by = "FAMID", suffix = c("", ".a"))
combined_54 <- combined_54 %>% dplyr::select(-contains(".a"))
combined_54 <- dplyr::left_join(combined_54, month24, by = "FAMID", suffix = c("", ".a"))
combined_54 <- combined_54 %>% dplyr::select(-contains(".a"))

combined_54 <- combined_54[c("DP3_M54M",#food bank
                             "ls3_y2m", #income meet daily needs
                             "finhhp_7_m45m",#unemployment benefit
                             "finhhp_9_m45m" ,#domestic purposes/sole parent support 
                             "finhhp_8_m45m" ,#sickness
                             "finhhp_10_m45m" ,#invalid
                             "finhhp_11_m45m" , #student allowance
                             
                             "finhhp_18_m45m",#family tax credits
                             "finhhp_16_m45m", #accomodation supplement
                             "finhhp_12_m45m",#Other govt support
                             "finhhp_s4_m45m", #Other government benefits
                             "sum" )]

#"finhhp_s3_m45m", #income tested benefit

table(combined_54$DP3_M54M,combined_54$ls3_y2m, exclude = NULL )
combined_54$foodbank[combined_54$DP3_M54M == 0] <- 0
combined_54$notenough[combined_54$ls3_y2m == 99 ] <- NA
combined_54$notenough[combined_54$ls3_y2m == 98 ] <- NA
combined_54$notenough[combined_54$ls3_y2m >= 2 ] <- 0
combined_54$foodbank[combined_54$DP3_M54M == 1] <- 1
combined_54$notenough[combined_54$ls3_y2m == 1 ] <- 1
table(combined_54$notenough, combined_54$ls3_y2m, exclude = NULL)

combined_54$benefit[combined_54$finhhp_7_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_7_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_7_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_9_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_9_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_9_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_8_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_8_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_8_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_10_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_10_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_10_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_11_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_11_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_11_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_12_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_12_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_12_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_16_m45m == 98] <- NA
combined_54$benefit[combined_54$finhhp_16_m45m == 99] <- NA
combined_54$benefit[combined_54$finhhp_16_m45m == 0] <- 0

combined_54$benefit[combined_54$finhhp_s4_m45m == 0] <- 0

combined_54$added[combined_54$finhhp_7_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_9_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_8_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_10_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_11_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_16_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_18_m45m == 0] <- 0
combined_54$added[combined_54$finhhp_7_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_9_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_8_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_10_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_11_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_16_m45m == 1] <- 1
combined_54$added[combined_54$finhhp_18_m45m == 1] <- 2

combined_54$benefit[combined_54$finhhp_s4_m45m == 1 & combined_54$added == 0] <- 1

combined_54$benefit[combined_54$finhhp_7_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_9_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_8_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_10_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_11_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_16_m45m == 1] <- 1
combined_54$benefit[combined_54$finhhp_12_m45m == 1] <- 1

combined_54$famtax[combined_54$finhhp_18_m45m == 98] <- NA
combined_54$famtax[combined_54$finhhp_18_m45m == 99] <- NA
combined_54$famtax[combined_54$finhhp_18_m45m == 0] <- 0
combined_54$famtax[combined_54$finhhp_18_m45m == 1] <- 1

test <- combined_54[c("sum", "benefit", "foodbank", "notenough", "famtax")]
test <- na.omit(test)

test$validation[test$benefit == 1 & test$notenough == 1 & test$foodbank == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 0 & test$notenough == 1 & test$foodbank == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 1 & test$notenough == 1 & test$foodbank == 1 & test$famtax == 0] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"

test$validation[test$benefit == 1 & test$notenough == 0 & test$foodbank == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 0 & test$notenough == 0 & test$foodbank == 1 & test$famtax == 1] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"
test$validation[test$benefit == 1 & test$notenough == 0 & test$foodbank == 1 & test$famtax == 0] <- "Food bank & benefit receipt / Both economic strain experiences and benefit receipt"

test$validation[test$benefit == 0 & test$notenough == 1 & test$foodbank == 1 & test$famtax == 0] <- "Either economic strain experiences & no benefit receipt"
test$validation[test$benefit == 0 & test$notenough == 0 & test$foodbank == 1 & test$famtax == 0] <- "Either economic strain experiences & no benefit receipt"
test$validation[test$benefit == 0 & test$notenough == 1 & test$foodbank == 0 & test$famtax == 0] <- "Either economic strain experiences & no benefit receipt"

test$validation[test$benefit == 1 & test$notenough == 1 & test$foodbank == 0 & test$famtax == 0] <- "Not enough and Benefit receipt"
test$validation[test$benefit == 0 & test$notenough == 1 & test$foodbank == 0 & test$famtax == 1] <- "Not enough and Benefit receipt"
test$validation[test$benefit == 1 & test$notenough == 1 & test$foodbank == 0 & test$famtax == 1] <- "Not enough and Benefit receipt"

test$validation[test$benefit == 1 & test$notenough == 0 & test$foodbank == 0 & test$famtax == 1] <- "Benefit receipt only"
test$validation[test$benefit == 1 & test$notenough == 0 & test$foodbank == 0 & test$famtax == 0] <- "Benefit receipt only"
test$validation[test$benefit == 1 & test$notenough == 0 & test$foodbank == 0 & test$famtax == 1] <- "Benefit receipt only"

test$validation[test$benefit == 0 & test$notenough == 0 & test$foodbank == 0 & test$famtax == 1] <- "Fam tax only"
test$validation[test$benefit == 0 & test$notenough == 0 & test$foodbank == 0 & test$famtax == 0] <- "None"


test$validation <- factor(test$validation, levels = c("Food bank & benefit receipt / Both economic strain experiences and benefit receipt", 
                                                      "Either economic strain experiences & no benefit receipt",
                                                      "Not enough and Benefit receipt", "Benefit receipt only",
                                                      "Fam tax only", "None"))

table(test$validation, exclude = NULL)
table(test$sum)
table(test$validation, test$sum, exclude = NULL)

### Applying the cut point ####
DEP_54M$sum <- rowSums((DEP_54M[,2:6]))
table(DEP_54M$sum, exclude = NULL)
DEP_54M$threshold54[DEP_54M$sum >=3] <- 1
DEP_54M$threshold54[DEP_54M$sum < 3] <- 0
table(DEP_54M$threshold54, exclude = NULL)

#### Comparison of all three methods ####

DEP_54M <- dplyr::left_join(DEP_54M, LCA_54M, by = "FAMID", suffix = c("", ".a"))
DEP_54M <- DEP_54M %>% dplyr::select(-contains(".a"))

table(DEP_54M$LCA_3)
table(DEP_54M$anyhardship)
DEP_54M$LCA_54M_revised[DEP_54M$LCA_3 == 3] <- 0
DEP_54M$LCA_54M_revised[DEP_54M$LCA_3 == 1] <- 1
DEP_54M$LCA_54M_revised[DEP_54M$LCA_3 == 2] <- 1

table(DEP_54M$threshold54 , DEP_54M$LCA_54M_revised)
table(DEP_54M$threshold54 , DEP_54M$anyhardship)
table(DEP_54M$LCA_54M_revised , DEP_54M$anyhardship)

DEP_54M$congruence <- ifelse(DEP_54M$anyhardship == 1 & DEP_54M$threshold54 == 1 & DEP_54M$LCA_54M_revised == 1, 1, 
                            ifelse(DEP_54M$anyhardship == 0 & DEP_54M$threshold54 == 0 & DEP_54M$LCA_54M_revised == 0, 0, 2))
table(DEP_54M$congruence)
table(DEP_54M$congruence,DEP_54M$anyhardship )
