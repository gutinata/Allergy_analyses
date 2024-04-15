install.packages("pROC", dependencies = T)
library("pROC")

############
#Tropomyosin
############

trop <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tropomyosin_identity.csv", header = T, row.names = 1)
trop <- as.matrix(trop)

trop_groups[1:10,]

trop_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tropomyosin_groups.csv", header = T)

#delete fish data from main matrices
trop_nf <- trop[-(64:67),-(64:67)]
trop_groups_nf <- trop_groups[-(64:67),]

tropomyosin_nf <- convert_matrix(trop_nf, trop_groups_nf, group_column = "group_type", type_column = "allergenic")
head(tropomyosin_nf)
levels(tropomyosin_nf$comparison_type)
summary(tropomyosin_nf$comparison)

tropomyosin_nf$allergenic <- NA
for (i in 1:nrow(tropomyosin_nf)){
  if (tropomyosin_nf$comparison_type[i] == "ALAL"){tropomyosin_nf$allergenic[i] <- 1}
  else if (tropomyosin_nf$comparison_type[i] == "ALNAL"){tropomyosin_nf$allergenic[i] <- 0}
  else {tropomyosin_nf$allergenic[i] <- NA}
}

tropomyosin_nf <- na.omit(tropomyosin_nf)

#order data frame by identity values
tropomyosin_nf <- tropomyosin_nf[order(tropomyosin_nf$identity), ]

#fit binomial
glm.fit=glm(tropomyosin_nf$allergenic~tropomyosin_nf$identity, family = binomial)
summary(glm.fit)
  #Coefficients:
  #             Estimate  Std. Error  z value   Pr(>|z|)    
  #(Intercept)  -35.3914      1.2871  -27.50    <2e-16 ***
  #  x$identity   0.6010      0.0223   26.95    <2e-16 ***
  #AIC: 2579.3

#plot data and model fit results
plot(tropomyosin_nf$identity, tropomyosin_nf$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(tropomyosin_nf$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(tropomyosin_nf$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(tropomyosin_nf$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
  #2340 controls (allergenic 0) < 2080 cases (allergenic 1)
  #Area under the curve: 93.53%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$tpp > 70 & roc.df$tpp < 90,]
  #        tpp        fpp thresholds
  #78 77.45192  0.4273504  0.6810743 
  #79 76.53846  0.0000000  0.7002112 *******

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #86.15385 14.52991  0.3559796

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(5.3, 4, 3.3, 2))

plot(tropomyosin_nf$identity, tropomyosin_nf$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(tropomyosin_nf$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(tropomyosin_nf$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(100, 76.53846, pch = 19)
text(98, 76.53846, paste("Threshold %ID value: 70% \nTrue positives: 76.5% \nFalse positives: 0%"), adj = c(0,1))
text(30, 20, "Area under the curve: 93.53%")




##########################
#Triosephosphate Isomerase
##########################

triose <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tim_identity.csv", header = T, row.names = 1)
triose <- as.matrix(triose)

triose_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tim_groups.csv", header = T)

tim <- convert_matrix(triose, triose_groups, group_column = "group_type", type_column = "allergenic")
head(tim)
levels(tim$comparison)
summary(tim$comparison)

tim$allergenic <- NA
for (i in 1:nrow(tim)){
  if (tim$comparison_type[i] == "ALAL"){tim$allergenic[i] <- 1}
  else if (tim$comparison_type[i] == "ALNAL"){tim$allergenic[i] <- 0}
  else {tim$allergenic[i] <- NA}
}

tim <- na.omit(tim)

#order data frame by identity values
tim <- tim[order(tim$identity), ]


#fit binomial
glm.fit=glm(tim$allergenic~tim$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#             Estimate  Std. Error  z value   Pr(>|z|)    
#(Intercept)  -0.44945    1.41228  -0.318    0.750
#tim$identity -0.01389    0.02088  -0.665    0.506
#AIC: 228.73

#plot data and model fit results
plot(tim$identity, tim$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(tim$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(tim$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(tim$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#180 controls (allergenic 0) > 45 cases (allergenic 1)
#Area under the curve: 47.39%%

#no further analyses since data cannot be discriminated. 


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(5.3, 4, 3.3, 2))

plot(tim$identity, tim$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(tim$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(tim$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
text(30, 20, "Area under the curve: 47.39%")



#################
#Beta Parvalbumin
#################

bparv <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/bparvalbumin_identity.csv", header = T, row.names = 1)
bparv <- as.matrix(bparv)

bparv_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/bparvalbumin_groups.csv", header = T)

bparvalbumin <- convert_matrix(bparv, bparv_groups, group_column = "group_type", type_column = "allergenic")
head(bparvalbumin)
levels(bparvalbumin$comparison_type)
summary(bparvalbumin$comparison)

bparvalbumin$allergenic <- NA
for (i in 1:nrow(bparvalbumin)){
  if (bparvalbumin$comparison_type[i] == "ALAL"){bparvalbumin$allergenic[i] <- 1}
  else if (bparvalbumin$comparison_type[i] == "ALNAL"){bparvalbumin$allergenic[i] <- 0}
  else {bparvalbumin$allergenic[i] <- NA}
}

bparvalbumin <- na.omit(bparvalbumin)

#order data frame by identity values
bparvalbumin <- bparvalbumin[order(bparvalbumin$identity), ]

#fit binomial
glm.fit=glm(bparvalbumin$allergenic~bparvalbumin$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#             Estimate  Std. Error  z value   Pr(>|z|)    
#(Intercept)           -36.90262    1.82642  -20.20   <2e-16 ***
#bparvalbumin$identity   0.64437    0.03237   19.91   <2e-16 ***
#AIC: 609.72

#plot data and model fit results
plot(bparvalbumin$identity, bparvalbumin$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(bparvalbumin$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(bparvalbumin$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(bparvalbumin$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#1344 controls (allergenic 0) < 1540 cases (allergenic 1)
#Area under the curve: 99.15%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$tpp > 90 & roc.df$tpp < 100,]
#        tpp        fpp   thresholds
#84.15584    0.07440476  0.956061208
#84.09091    0.00000000  0.961024617 
#94.09091     1.4136905  0.708771493 *****

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #96.75325 4.836310  0.3561168

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.6, 4, 2.6, 2))

plot(bparvalbumin$identity, bparvalbumin$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(bparvalbumin$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(bparvalbumin$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(100, 84.15584, pch = 19)
text(98, 84.15584, paste("Threshold %ID value: 70% \nTrue positives: 94.5% \nFalse positives: 1.4%"), adj = c(0,1))
text(30, 20, "Area under the curve: 99.15%")


###################
#Heat Shock Protein
###################

hsp <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/HSP_identity.csv", header = T, row.names = 1)
hsp <- as.matrix(hsp)

hsp_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/HSP_groups.csv", header = T)

heatshockp <- convert_matrix(hsp, hsp_groups, group_column = "group_type", type_column = "allergenic")
head(heatshockp)
levels(heatshockp$comparison_type)
summary(heatshockp$comparison)

heatshockp$allergenic <- NA
for (i in 1:nrow(heatshockp)){
  if (heatshockp$comparison_type[i] == "ALAL"){heatshockp$allergenic[i] <- 1}
  else if (heatshockp$comparison_type[i] == "ALNAL"){heatshockp$allergenic[i] <- 0}
  else {heatshockp$allergenic[i] <- NA}
}

heatshockp <- na.omit(heatshockp)

#order data frame by identity values
heatshockp <- heatshockp[order(heatshockp$identity), ]

#fit binomial
glm.fit=glm(heatshockp$allergenic~heatshockp$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#             Estimate  Std. Error  z value   Pr(>|z|)    
#(Intercept)         -2.43488    1.01037  -2.410    0.016 *
#heatshockp$identity  0.02461    0.01595   1.543    0.123  
#AIC: 120.07

#plot data and model fit results
plot(heatshockp$identity, heatshockp$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(heatshockp$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(heatshockp$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(heatshockp$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#72 controls (allergenic 0) > 28 cases (allergenic 1)
#Area under the curve: 58.13%

#no further analyses since data cannot be discriminated. 

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.6, 4, 2.6, 2))

plot(heatshockp$identity, heatshockp$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(heatshockp$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(heatshockp$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
text(30, 20, "Area under the curve: 58.13%")




################
#Arginine kinase
################

ak <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/argkin_identity.csv", header = T, row.names = 1)
ak <- as.matrix(ak)

ak_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/argkin_groups.csv", header = T)
levels(as.factor(ak_groups$group_type))

argkin <- convert_matrix(ak, ak_groups, group_column = "group_type", type_column = "allergenic")
head(argkin)
levels(argkin$comparison_type)
summary(argkin$comparison)

summary(argkin)

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(argkin)
argkin[1379,] <- c("InsectMollusk","ALAL", min(argkin$identity))
argkin[1380,] <- c("MammalMollusk","ALNAL", max(argkin$identity))

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
#min(argkin$identity[argkin$comparison_type == "ALAL"])
#max(argkin$identity[argkin$comparison_type == "ALNAL"])
#(50.5+52.9)/2
#argkin[1379,] <- c("InsectMollusk","ALAL", 50.5)
#argkin[1380,] <- c("MammalMollusk","ALNAL", 52.9)

tail(argkin)

argkin$identity <- as.numeric(argkin$identity)

argkin$allergenic <- NA
for (i in 1:nrow(argkin)){
  if (argkin$comparison_type[i] == "ALAL"){argkin$allergenic[i] <- 1}
  else if (argkin$comparison_type[i] == "ALNAL"){argkin$allergenic[i] <- 0}
  else {argkin$allergenic[i] <- NA}
}

argkin <- na.omit(argkin)

#order data frame by identity values
argkin <- argkin[order(argkin$identity), ]

#binomial fit
glm.fit=glm(argkin$allergenic~argkin$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#             Estimate  Std. Error  z value   Pr(>|z|)    
#(Intercept)     -27.36364    2.70981 -10.098   <2e-16 ***
#argkin$identity   0.52822    0.05399   9.783   <2e-16 ***
#AIC: 99.231

#Try with a different model/library
#library("logistf")
#glm.fit=logistf(argkin$allergenic~argkin$identity)
#str(glm.fit)
#glm.fit$predict

#plot data and model fit results
plot(argkin$identity, argkin$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(argkin$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(argkin$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(argkin$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#702(+1) controls (allergenic 0) > 325(+1) cases (allergenic 1)
#Area under the curve: 99.55%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0 & roc.df$thresholds < 1,]
#     tpp       fpp thresholds
#99.69325 0.1422475  0.4876606
#99.38650 0.1422475  0.6469555 *****
#98.77301 0.1422475  0.6869639
#98.15951 0.1422475  0.7261914
#97.85276 0.1422475  0.7464258

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #99.693252  0.2844950 0.281467874

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.6, 4, 2.6, 2))

plot(argkin$identity, argkin$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(argkin$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(argkin$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(99.38650, 100, pch = 19)
text(98, 99.38650, paste("Threshold %ID value: 65% \nTrue positives: 99.4% \nFalse positives: 0.1%"), adj = c(0,1))
text(30, 20, "Area under the curve: 99.55%")



############
#Casein alpha
############

c_a <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/casein_alpha_identity.csv", header = T, row.names = 1)
c_a <- as.matrix(c_a)

c_a_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/casein_alpha_groups.csv", header = T)

caseina <- convert_matrix(c_a, c_a_groups, group_column = "group_type", type_column = "allergenic")
head(caseina)
levels(caseina$comparison_type)
summary(caseina$comparison)

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(caseina)
caseina[277,] <- c("BovinaeBovinae","ALAL", min(caseina$identity))
caseina[278,] <- c("BovinaeCaprinae","ALNAL", max(caseina$identity))

caseina$identity <- as.numeric(caseina$identity)

caseina$allergenic <- NA
for (i in 1:nrow(caseina)){
  if (caseina$comparison_type[i] == "ALAL"){caseina$allergenic[i] <- 1}
  else if (caseina$comparison_type[i] == "ALNAL"){caseina$allergenic[i] <- 0}
  else {caseina$allergenic[i] <- NA}
}
tail(caseina)
caseina <- na.omit(caseina)

#order data frame by identity values
caseina <- caseina[order(caseina$identity), ]

#fit binomial
glm.fit=glm(caseina$allergenic~caseina$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)      -12.6239       1.9434   -6.496      8.26e-11 ***
#caseina$identity   0.1763       0.0275    6.413      1.43e-10 ***
#AIC: 32.614

#plot data and model fit results
plot(caseina$identity, caseina$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(caseina$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(caseina$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(caseina$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#141 controls (allergenic 0) < 46 cases (allergenic 1)
#Area under the curve: 97.14%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.50 & roc.df$thresholds < 0.99,]
  #   tpp   fpp       thresholds
  #90.68100 1.1389522  0.6404182
  #89.96416 1.1389522  0.6538169****

roc.df[roc.df$thresholds >  0 & roc.df$thresholds < 0.99,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #97.82609  0.7092199 0.470824675
  #min identity is ~40% so no point estimating 35%


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(caseina$identity, caseina$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(caseina$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(caseina$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(100, (100-1.1389522), pch = 19)
text(89.96416, 89.96416, paste("Threshold %ID value:65% \nTrue positives: 90% \nFalse positives: 1.1%"), adj = c(0,1))
text(30, 20, "Area under the curve: 97.14%")


############
#Casein beta
############

c_b <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/casein_beta_identity.csv", header = T, row.names = 1)
c_b <- as.matrix(c_b)

c_b_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/casein_beta_groups.csv", header = T)

caseinb <- convert_matrix(c_b, c_b_groups, group_column = "group_type", type_column = "allergenic")
head(caseinb)
levels(caseinb$comparison_type)
summary(caseinb$comparison)

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(caseinb)
caseinb[277,] <- c("BovinaeBovinae","ALAL", min(caseinb$identity))
caseinb[278,] <- c("BovinaeCaprinae","ALNAL", max(caseinb$identity))

caseinb$identity <- as.numeric(caseinb$identity)

caseinb$allergenic <- NA
for (i in 1:nrow(caseinb)){
  if (caseinb$comparison_type[i] == "ALAL"){caseinb$allergenic[i] <- 1}
  else if (caseinb$comparison_type[i] == "ALNAL"){caseinb$allergenic[i] <- 0}
  else {caseinb$allergenic[i] <- NA}
}
tail(caseinb)
caseinb <- na.omit(caseinb)

#order data frame by identity values
caseinb <- caseinb[order(caseinb$identity), ]

#fit binomial
glm.fit=glm(caseinb$allergenic~caseinb$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)      -21.71604    3.43605   -6.32 2.61e-10 ***
#caseinb$identity   0.26500    0.04281    6.19 6.01e-10 ***
#AIC: 35.736

#plot data and model fit results
plot(caseinb$identity, caseinb$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(caseinb$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(caseinb$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(caseinb$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#136 controls (allergenic 0) < 37 cases (allergenic 1)
#Area under the curve: 96.58%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.40 & roc.df$thresholds < 0.99,]
  #   tpp   fpp       thresholds
  #97.29730 0.7352941  0.4569704
  #91.89189 0.7352941  0.8783868****

roc.df[roc.df$thresholds >  0 & roc.df$thresholds < 0.9,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #97.29730   0.7352941 0.4569704128
  #min identity is ~60% so no point estimating 35%

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(caseinb$identity, caseinb$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(caseinb$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(caseinb$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(91.89189, (100 - 0.7352941), pch = 19)
text(98, 91.89189, paste("Threshold %ID value: 87% \nTrue positives: 92% \nFalse positives: 0.7%"), adj = c(0,1))
text(30, 20, "Area under the curve: 96.58%")


#############
#Vicilin 7SC
#############

vic_7SC <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/vic_7SC_all_identity.csv", header = T, row.names = 1)

vic_7SC <- as.matrix(vic_7SC)

vic_7SC_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/vic_7SC_groups.csv", header = T)

vicilin7SC <- convert_matrix(vic_7SC, vic_7SC_groups, group_column = "group_type", type_column = "allergenic")
head(vicilin7SC)

within <- c("FabalesFabales","FagalesFagales","RosalesRosales","SapindalesSapindales","LamialesLamiales","ProtealesProteales","PinopsidaPinopsida")

vicilin7SC$comparison_type <- NA
vicilin7SC$allergenic <- NA
for (i in 1:nrow(vicilin7SC)){
  value = vicilin7SC$comparison[i]
  if (value %in% within) {
    vicilin7SC$comparison_type[i] <- "Within"
    vicilin7SC$allergenic[i] <- 1
    }
  else {
    vicilin7SC$comparison_type[i] <- "Between"
    vicilin7SC$allergenic[i] <- 0
    }
}
vicilin7SC$comparison_type <- factor(vicilin7SC$comparison_type, levels = c("Within","Between"))

vicilin7SC <- vicilin7SC[order(vicilin7SC$identity), ]

#fit binomial
glm.fit=glm(vicilin7SC$allergenic~vicilin7SC$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)         -38.09698    3.94546  -9.656   <2e-16 ***
#vicilin7SC$identity   0.71189    0.07174   9.923   <2e-16 ***
#AIC:132.7

#plot data and model fit results
plot(vicilin7SC$identity, vicilin7SC$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(vicilin7SC$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(vicilin7SC$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(vicilin7SC$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#506 controls (between 0) < 484 cases (withhin 1)
#Area under the curve: 99.57%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.40 & roc.df$thresholds < 0.70,]
  #   tpp   fpp       thresholds
  #99.17355 1.581028  0.5461722
  #98.96694 1.581028  0.6466486

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #99.17355 2.766798  0.3780110

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(vicilin7SC$identity, vicilin7SC$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(vicilin7SC$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(vicilin7SC$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(99.17355, (100 - 1.581028),  pch = 19)
points(98.96694, (100 - 1.581028),  pch = 19)
text(96, 94.26523, paste("Threshold %ID values: 55% / 65% \nTrue positives: 99% / 98% \nFalse positives: 1.6%"), adj = c(0,1))
text(40, 20, "Area under the curve: 99.57%")



#############
#Vicilin 7SN
#############

vic_7SN <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/vic_7SN_all_identity.csv", header = T, row.names = 1)

vic_7SN <- as.matrix(vic_7SN)

vic_7SN_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/vic_7SN_groups.csv", header = T)

vicilin7SN <- convert_matrix(vic_7SN, vic_7SN_groups, group_column = "group_type", type_column = "allergenic")
head(vicilin7SN)

within <- c("FabalesFabales","FagalesFagales","RosalesRosales","SapindalesSapindales","LamialesLamiales","ProtealesProteales","PinopsidaPinopsida")

vicilin7SN$comparison_type <- NA
vicilin7SN$allergenic <- NA
for (i in 1:nrow(vicilin7SN)){
  value = vicilin7SN$comparison[i]
  if (value %in% within) {
    vicilin7SN$comparison_type[i] <- "Within"
    vicilin7SN$allergenic[i] <- 1
  }
  else {
    vicilin7SN$comparison_type[i] <- "Between"
    vicilin7SN$allergenic[i] <- 0
  }
}
vicilin7SN$comparison_type <- factor(vicilin7SN$comparison_type, levels = c("Within","Between"))

vicilin7SN <- vicilin7SN[order(vicilin7SN$identity), ]

#fit binomial
glm.fit=glm(vicilin7SN$allergenic~vicilin7SN$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)         -17.84136    1.18813  -15.02   <2e-16 ***
#vicilin7SN$identity   0.37519    0.02478   15.14   <2e-16 ***
#AIC:418.82

#plot data and model fit results
plot(vicilin7SN$identity, vicilin7SN$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(vicilin7SN$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(vicilin7SN$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(vicilin7SN$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#478 controls (between 0) < 425 cases (within 1)
#Area under the curve: 96.66%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.40 & roc.df$thresholds < 0.70,]
  #   tpp   fpp       thresholds
  #92.70588 7.531381  0.5618208******

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #96.70588 9.205021  0.3511251

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(vicilin7SN$identity, vicilin7SN$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(vicilin7SN$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(vicilin7SN$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(92.70588, (100 - 7.531381), pch = 19)
text(91, 92.70588, paste("Threshold %ID value: 55% \nTrue positives: 92% \nFalse positives: 7.5%"), adj = c(0,1))
text(40, 20, "Area under the curve: 96%")


###########
#2S albumin
###########

twoSalbumin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/2Salb_identity.csv", header = T, row.names = 1)

twoSalbumin <- as.matrix(twoSalbumin)

twoSalb_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/2Salb_groups.csv", header = T)


twoSalbumin_all <- convert_matrix(twoSalbumin, twoSalb_groups, group_column = "group_type", type_column = "allergenic")
head(twoSalbumin_all)

within <- c("FabalesFabales","FagalesFagales","MalpighialesMalpighiales","CucurbitalesCucurbitales","SapindalesSapindales","BrassicalesBrassicales","LamialesLamiales","EricalesEricales")

twoSalbumin_all$comparison_type <- NA
twoSalbumin_all$allergenic <- NA
for (i in 1:nrow(twoSalbumin_all)){
  value = twoSalbumin_all$comparison[i]
  if (value %in% within) {
    twoSalbumin_all$comparison_type[i] <- "Within"
    twoSalbumin_all$allergenic[i] <- 1
  }
  else {
    twoSalbumin_all$comparison_type[i] <- "Between"
    twoSalbumin_all$allergenic[i] <- 0
  }
}
twoSalbumin_all$comparison_type <- factor(twoSalbumin_all$comparison_type, levels = c("Within","Between"))

twoSalbumin_all <- twoSalbumin_all[order(twoSalbumin_all$identity), ]

#fit binomial
glm.fit=glm(twoSalbumin_all$allergenic~twoSalbumin_all$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)              -6.21154    0.39083  -15.89   <2e-16 ***
#twoSalbumin_all$identity  0.10789    0.00914   11.80   <2e-16 ***
#AIC:264.81

#plot data and model fit results
plot(twoSalbumin_all$identity, twoSalbumin_all$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(twoSalbumin_all$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(twoSalbumin_all$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(twoSalbumin_all$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#848 controls (between 0) < 98 cases (within 1)
#Area under the curve: 97.56%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.70 & roc.df$thresholds < 0.80,]
  #   tpp   fpp       thresholds
  #70.40816 1.179245  0.4001057****
  #66.32653 1.179245  0.4672287

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #70.40816 1.29717  0.3364538


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(twoSalbumin_all$identity, twoSalbumin_all$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(twoSalbumin_all$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(twoSalbumin_all$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(99, 70.40816, pch = 19)
text(97, 70.40816, paste("Threshold %ID value: 40% \nTrue positives: 70% \nFalse positives: 1.1%"), adj = c(0,1))
text(30, 20, "Area under the curve: 97.56%")

#############
#11S globulin
#############

eleven <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/11S_identity.csv", header = T, row.names = 1)

eleven <- as.matrix(eleven)

eleven_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/11S_groups.csv", header = T)

eleven_all <- convert_matrix(eleven, eleven_groups, group_column = "group_type", type_column = "allergenic")
head(eleven_all)

within <- c("FabalesFabales","FagalesFagales","RosalesRosales","CucurbitalesCucurbitales","SapindalesSapindales","BrassicalesBrassicales","LamialesLamiales","EricalesEricales","CaryophyllalesCaryophyllales")

eleven_all$comparison_type <- NA
eleven_all$allergenic <- NA
for (i in 1:nrow(eleven_all)){
  value = eleven_all$comparison[i]
  if (value %in% within) {
    eleven_all$comparison_type[i] <- "Within"
    eleven_all$allergenic[i] <- 1
  }
  else {
    eleven_all$comparison_type[i] <- "Between"
    eleven_all$allergenic[i] <- 0
  }
}

eleven_all$comparison_type <- factor(eleven_all$comparison_type, levels = c("Within","Between"))

eleven_all <- eleven_all[order(eleven_all$identity), ]

#fit binomial
glm.fit=glm(eleven_all$allergenic~eleven_all$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)         -12.82917    1.71002  -7.502 6.27e-14 ***
#eleven_all$identity   0.20673    0.03269   6.324 2.56e-10 ***
#AIC:131.54

#plot data and model fit results
plot(eleven_all$identity, eleven_all$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(eleven_all$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(eleven_all$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(eleven_all$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#360 controls (between 0) < 46 cases (within 1)
#Area under the curve: 88.48%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$tpp > 60 & roc.df$tpp < 80,]
  #   tpp   fpp       thresholds
  #71.73913   0  0.4324693
  #69.56522   0  0.4601106*******

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #73.91304   0  0.3616588


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(eleven_all$identity, eleven_all$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(eleven_all$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(eleven_all$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(100, 69.56522, pch = 19)
text(98, 69.56522, paste("Threshold %ID value: 46% \nTrue positives: 70% \nFalse positives: 0%"), adj = c(0,1))
text(30, 20, "Area under the curve: 88.48%")

##########
#thaumatin
##########

thaumatin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/thaumatin_identity.csv", header = T, row.names = 1)

thaumatin <- as.matrix(thaumatin)

thau_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/thaumatin_groups.csv", header = T)

thau <- convert_matrix(thaumatin, thau_groups, group_column = "group_type", type_column = "allergenic")
summary(thau)

within <- c("LamialesLamiales","SolanalesSolanales","EricalesEricales","RosalesRosales","PoalesPoales","ZingiberalesZingiberales","CupressalesCupressales")

thau$comparison_type <- NA
thau$allergenic <- NA
for (i in 1:nrow(thau)){
  value = thau$comparison[i]
  if (value %in% within) {
    thau$comparison_type[i] <- "Within"
    thau$allergenic[i] <- 1
  }
  else {
    thau$comparison_type[i] <- "Between"
    thau$allergenic[i] <- 0
  }
}

thau$comparison_type <- factor(thau$comparison_type, levels = c("Within","Between"))

thau <- thau[order(thau$identity), ]

#fit binomial
glm.fit=glm(thau$allergenic~thau$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)   -19.05230    2.78586  -6.839 7.98e-12 ***
#thau$identity   0.27767    0.04116   6.746 1.52e-11 ***
#AIC:74.803

#plot data and model fit results
plot(thau$identity, thau$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(thau$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(thau$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(thau$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#298 controls (between 0) > 108 cases (within 1)
#Area under the curve: 99.43%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.60 & roc.df$thresholds < 0.90,]
  #   tpp   fpp       thresholds
  #86.11111 2.013423  0.6083040*******
  #84.25926 2.013423  0.6452237

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #98.14815 4.026846  0.3685744


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(thau$identity, thau$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(thau$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(thau$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(98, 86.11111, pch = 19)
text(96, 86.11111, paste("Threshold %ID value: 60% \nTrue positives: 86% \nFalse positives: 2%"), adj = c(0,1))
text(30, 20, "Area under the curve: 99.43%")

#################
#gibberellin-like
#################

gibberellin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/gibberellin_identity.csv", header = T, row.names = 1)
gibberellin <- as.matrix(gibberellin)

gib_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/gibberellin_groups.csv", header = T)

gib <- convert_matrix(gibberellin, gib_groups, group_column = "group_type", type_column = "allergenic")
head(gib)

within <- c("RosalesRosales","SapindalesSapindales","SolanalesSolanales","CupressalesCupressales")

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(gib)
gib[37,] <- c("RosalesRosales","ALAL", min(gib$identity))
gib[38,] <- c("RosalesSapindales","ALAL", max(gib$identity))

tail(gib)

gib$identity <- as.numeric(gib$identity)

gib$comparison_type <- NA
gib$allergenic <- NA
for (i in 1:nrow(gib)){
  value = gib$comparison[i]
  if (value %in% within) {
    gib$comparison_type[i] <- "Within"
    gib$allergenic[i] <- 1
  }
  else {
    gib$comparison_type[i] <- "Between"
    gib$allergenic[i] <- 0
  }
}

gib$comparison_type <- factor(gib$comparison_type, levels = c("Within","Between"))

gib <- gib[order(gib$identity), ]

#fit binomial
glm.fit=glm(gib$allergenic~gib$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)  -15.41583    4.90711  -3.142  0.00168 **
#gib$identity   0.17218    0.05536   3.110  0.00187 **
#AIC:27.037

#plot data and model fit results
plot(gib$identity, gib$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(gib$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(gib$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(gib$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#28 controls (between 0) > 10 cases (within 1)
#Area under the curve: 88.21%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.60 & roc.df$thresholds < 0.90,]
  #   tpp       fpp       thresholds
  #     70 3.571429  0.6385076****
  #     60 3.571429  0.7241790

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #90 10.71429  0.3729836


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(gib$identity, gib$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(gib$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(gib$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(97, 70, pch = 19)
text(95, 70, paste("Threshold %ID value: 63% \nTrue positives: 70% \nFalse positives: 3.5%"), adj = c(0,1))
text(30, 20, "Area under the curve: 88.21%")


##########
#Lipocalin
##########

lipoc <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/lipocalin_identity.csv", header = T, row.names = 1)
lipoc <- as.matrix(trop)

lipoc_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/lipocalin_groups.csv", header = T)

lipocalin <- convert_matrix(lipoc, lipoc_groups, group_column = "group_type", type_column = "allergenic")
head(lipocalin)
levels(lipocalin$comparison_type)
summary(lipocalin$comparison_type)

for (i in 1:nrow(lipocalin)) {
  if (is.na(lipocalin$comparison_type[i])) {lipocalin$allergenic[i] <- NA} 
  else if (lipocalin$comparison_type[i] == "ALAL") {lipocalin$allergenic[i] <- 1} 
  else if (lipocalin$comparison_type[i] == "ALNAL") {lipocalin$allergenic[i] <- 0} 
  else {lipocalin$allergenic[i] <- NA}
}
lipocalin <- na.omit(lipocalin)

#order data frame by identity values
lipocalin <- lipocalin[order(lipocalin$identity), ]

#fit binomial
glm.fit=glm(lipocalin$allergenic~lipocalin$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#             Estimate  Std. Error  z value   Pr(>|z|)    
#(Intercept)        -0.89852    1.87295  -0.480    0.631
#lipocalin$identity  0.04505    0.03328   1.354    0.176
#AIC: 68.334

#plot data and model fit results
plot(lipocalin$identity, lipocalin$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(lipocalin$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(lipocalin$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(lipocalin$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#12 controls (allergenic 0) < 66 cases (allergenic 1)
#Area under the curve: 56.31%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$tpp > 30 & roc.df$tpp < 50,]
#        tpp        fpp thresholds
#37.87879  8.333333  0.8653238
#37.87879  0.000000  0.8724022

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(5.3, 4, 3.3, 2))

plot(lipocalin$identity, lipocalin$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(lipocalin$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(lipocalin$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
text(30, 20, "Area under the curve: 56.31%")



#############
#LTP
#############

LTP <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/LTP_identity.csv", header = T, row.names = 1)

LTP <- as.matrix(LTP)

LTP_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/LTP_groups.csv", header = T)


LTPs <- convert_matrix(LTP, LTP_groups, group_column = "group_type", type_column = "allergenic")
head(LTPs)

within <- c("PoalesPoales","ProtealesProteales","ApialesApiales","AsteralesAsterales","BrassicalesBrassicales","EricalesEricales","FabalesFabales","FagalesFagales","MalpighialesMalpighiales","MyrtalesMyrtales","RosalesRosales","SapindalesSapindales","SolanalesSolanales","VitalesVitales")

LTPs$comparison_type <- NA
LTPs$allergenic <- NA
for (i in 1:nrow(LTPs)){
  value = LTPs$comparison[i]
  if (value %in% within) {
    LTPs$comparison_type[i] <- "Within"
    LTPs$allergenic[i] <- 1
  }
  else {
    LTPs$comparison_type[i] <- "Between"
    LTPs$allergenic[i] <- 0
  }
}
LTPs$comparison_type <- factor(LTPs$comparison_type, levels = c("Within","Between"))

LTPs <- LTPs[order(LTPs$identity), ]

#fit binomial
glm.fit=glm(LTPs$allergenic~LTPs$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)   -22.14055    1.15251  -19.21   <2e-16 ***
#LTPs$identity   0.32715    0.01766   18.52   <2e-16 ***
#AIC:794.25

#plot data and model fit results
plot(LTPs$identity, LTPs$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(LTPs$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(LTPs$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(LTPs$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#1679 controls (between 0) < 466 cases (withhin 1)
#Area under the curve: 95.62%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.40 & roc.df$thresholds < 0.90,]
  #tpp   fpp       thresholds
  #83.90558 2.08457415  0.4762183
  #79.82833 1.13162597  0.5067899****

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #86.69528 4.526504  0.3482661


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(LTPs$identity, LTPs$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(LTPs$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(LTPs$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(99, 79.82833, pch = 19)
text(97, 79.82833, paste("Threshold %ID values: 50% \nTrue positives: 80% \nFalse positives: 1%"), adj = c(0,1))
text(40, 20, "Area under the curve: 95.62%")



##########
#Ovalbumin
##########

ovalbumin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/ovalbumin_identity.csv", header = T, row.names = 1)
ovalbumin <- as.matrix(ovalbumin)

ovalb_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/ovalbumin_groups.csv", header = T)

ovalb <- convert_matrix(ovalbumin, ovalb_groups, group_column = "group_type", type_column = "allergenic")
head(ovalb)

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(ovalb)
ovalb[137,] <- c("GalliformesGalliformes","ALAL", min(ovalb$identity))
ovalb[138,] <- c("GalliformesHuman","ALNAL", max(ovalb$identity))

tail(ovalb)

ovalb$identity <- as.numeric(ovalb$identity)

ovalb$allergenic <- NA
for (i in 1:nrow(ovalb)){
  if (ovalb$comparison_type[i] == "ALAL"){ovalb$allergenic[i] <- 1}
  else if (ovalb$comparison_type[i] == "ALNAL"){ovalb$allergenic[i] <- 0}
  else {ovalb$allergenic[i] <- NA}
}

ovalb <- na.omit(ovalb)

#order data frame by identity values
ovalb <- ovalb[order(ovalb$identity), ]

#binomial fit
glm.fit=glm(ovalb$allergenic~ovalb$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#             Estimate  Std. Error  z value   Pr(>|z|)    
#(Intercept)    -11.67334    1.96925  -5.928 3.07e-09 ***
#ovalb$identity   0.18767    0.03112   6.030 1.64e-09 ***
#AIC: 32.886

#plot data and model fit results
plot(ovalb$identity, ovalb$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(ovalb$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(ovalb$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(ovalb$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#61(+1) controls (allergenic 0) > 67(+1) cases (allergenic 1)
#Area under the curve: 96.99%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 45 & roc.df$thresholds < 75,]
  #     tpp       fpp thresholds
  #98.507463   1.639344 0.46191293
  #95.522388   1.639344 0.88047081*****

roc.df[roc.df$thresholds >  0 & roc.df$thresholds < 0.9,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #98.50746  1.639344 0.46191293

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.6, 4, 2.6, 2))

plot(ovalb$identity, ovalb$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(ovalb$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(ovalb$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(99, 95.522388, pch = 19)
text(97, 95.522388, paste("Threshold %ID value: 88% \nTrue positives: 95% \nFalse positives: 1.6%"), adj = c(0,1))
text(30, 20, "Area under the curve: 96.99%")


#############
#Profilin
#############

profilin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/profilin_identity.csv", header = T, row.names = 1)

profilin <- as.matrix(profilin)

prof_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/profilin_groups.csv", header = T)

prof <- convert_matrix(profilin, prof_groups, group_column = "group_type", type_column = "allergenic")
head(prof)

within <- c("ApialesApiales","AsteralesAsterales","SolanalesSolanales","LamialesLamiales","EricalesEricales","CucurbitalesCucurbitales","FabalesFabales","FagalesFagales","MalpighialesMalpighiales","RosalesRosales","BrassicalesBrassicales","SapindalesSapindales","CaryophyllalesCaryophyllales","AsparagalesAsparagales","ArecalesArecales","PoalesPoales","ZingiberalesZingiberales","ProtealesProteales")

prof$comparison_type <- NA
prof$allergenic <- NA
for (i in 1:nrow(prof)){
  value = prof$comparison[i]
  if (value %in% within) {
    prof$comparison_type[i] <- "Within"
    prof$allergenic[i] <- 1
  }
  else {
    prof$comparison_type[i] <- "Between"
    prof$allergenic[i] <- 0
  }
}
prof$comparison_type <- factor(prof$comparison_type, levels = c("Within","Between"))

prof <- prof[order(prof$identity), ]

#fit binomial
glm.fit=glm(prof$allergenic~prof$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)   -17.429741   0.526108  -33.13   <2e-16 ***
#prof$identity   0.191186   0.006411   29.82   <2e-16 ***
#AIC:5817

#plot data and model fit results
plot(prof$identity, prof$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(prof$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(prof$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(prof$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#8329 controls (between 0) < 1124 cases (within 1)
#Area under the curve: 73.9%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.40 & roc.df$thresholds < 0.90,]
  #tpp   fpp       thresholds
  #24.3772242 1.34469924  0.4029718****
  #24.3772242 1.30868051  0.4103415

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #25.26690 1.788930  0.3539941


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(prof$identity, prof$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(prof$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(prof$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(99, 24.3772242, pch = 19)
text(97, 24.3772242, paste("Threshold %ID values: 40% \nTrue positives: 250% \nFalse positives: 1%"), adj = c(0,1))
text(40, 50, "Area under the curve: 73.9%")

######################
#PR-10 -- Bet v 1 like
######################

pr10_m <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/PR10_identity.csv", header = T, row.names = 1)

pr10_m <- as.matrix(pr10_m)

pr10_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/PR10_groups.csv", header = T)

pr10 <- convert_matrix(pr10_m, pr10_groups, group_column = "group_type", type_column = "allergenic")
head(pr10)

within <- c("ApialesApiales","EricalesEricales","SolanalesSolanales","FabalesFabales","FagalesFagales","RosalesRosales","SapindalesSapindales")

pr10$comparison_type <- NA
pr10$allergenic <- NA
for (i in 1:nrow(pr10)){
  value = pr10$comparison[i]
  if (value %in% within) {
    pr10$comparison_type[i] <- "Within"
    pr10$allergenic[i] <- 1
  }
  else {
    pr10$comparison_type[i] <- "Between"
    pr10$allergenic[i] <- 0
  }
}
pr10$comparison_type <- factor(pr10$comparison_type, levels = c("Within","Between"))

pr10 <- pr10[order(pr10$identity), ]

#fit binomial
glm.fit=glm(pr10$allergenic~pr10$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)   -16.949458   0.281831  -60.14   <2e-16 ***
#pr10$identity   0.269272   0.004587   58.71   <2e-16 ***
#AIC:8107.1

#plot data and model fit results
plot(pr10$identity, pr10$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(pr10$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(pr10$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(pr10$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#10261 controls (between 0) < 7130 cases (within 1)
#Area under the curve: 96.22%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.40 & roc.df$thresholds < 0.90,]
  #tpp   fpp       thresholds
  #80.78541 2.8944547  0.6576530***

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #86.71809 11.110028  0.3499102


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(pr10$identity, pr10$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(pr10$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(pr10$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(97, 80.78541, pch = 19)
text(95, 80.78541, paste("Threshold %ID values: 65% \nTrue positives: 80% \nFalse positives: 3%"), adj = c(0,1))
text(40, 50, "Area under the curve: 96.22%")


############
#Cyclophilin
############

cyclophilin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/cyclophilin_identity.csv", header = T, row.names = 1)
cyclophilin <- as.matrix(cyclophilin)

cyc_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/cyclophilin_groups.csv", header = T)
levels(as.factor(cyc_groups$group_type))

cyc <- convert_matrix(cyclophilin, cyc_groups, group_column = "group_type", type_column = "allergenic")
head(cyc)
levels(cyc$comparison_type)
summary(cyc$comparison)

cyc$identity <- as.numeric(cyc$identity)

cyc$allergenic <- NA
for (i in 1:nrow(cyc)){
  if (cyc$comparison_type[i] == "ALAL"){cyc$allergenic[i] <- 1}
  else if (cyc$comparison_type[i] == "ALNAL"){cyc$allergenic[i] <- 0}
  else {cyc$allergenic[i] <- NA}
}

cyc <- na.omit(cyc)

#order data frame by identity values
cyc <- cyc[order(cyc$identity), ]

#binomial fit
glm.fit=glm(cyc$allergenic~cyc$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#             Estimate  Std. Error  z value   Pr(>|z|)    
#(Intercept)  -4.32024    2.87144  -1.505    0.132
#cyc$identity  0.05998    0.03881   1.546    0.122
#AIC: 118.96

#plot data and model fit results
plot(cyc$identity, cyc$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(cyc$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(cyc$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(cyc$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#40 controls (allergenic 0) > 45 cases (allergenic 1)
#Area under the curve: 50.42%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0 & roc.df$thresholds < 1,]
#     tpp       fpp thresholds

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.6, 4, 2.6, 2))

plot(cyc$identity, cyc$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(cyc$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(cyc$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
#points(100, 99.38650, pch = 19)
#text(98, 99.38650, paste("Threshold %ID value: 65% \nTrue positives: 99.4% \nFalse positives: 0.1%"), adj = c(0,1))
text(30, 20, "Area under the curve: 50.42%")



##########
#Polcalcin
##########

polcalcin <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/polcalcin_identity.csv", header = T, row.names = 1)

polcalcin <- as.matrix(polcalcin)

polc_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/polcalcin_groups.csv", header = T)

polc <- convert_matrix(polcalcin, polc_groups, group_column = "group_type", type_column = "allergenic")
head(polc)

within <- c("AsteralesAsterales","LamialesLamiales","FagalesFagales","RosalesRosales","BrassicalesBrassicales","CaryophyllalesCaryophyllales","PoalesPoales")

polc$comparison_type <- NA
polc$allergenic <- NA
for (i in 1:nrow(polc)){
  value = polc$comparison[i]
  if (value %in% within) {
    polc$comparison_type[i] <- "Within"
    polc$allergenic[i] <- 1
  }
  else {
    polc$comparison_type[i] <- "Between"
    polc$allergenic[i] <- 0
  }
}
polc$comparison_type <- factor(polc$comparison_type, levels = c("Within","Between"))

polc <- polc[order(polc$identity), ]

#fit binomial
glm.fit=glm(polc$allergenic~polc$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)   -52.6467    20.0049  -2.632  0.00850 **
#polc$identity   0.6446     0.2485   2.594  0.00949 **
#AIC:22.981

#plot data and model fit results
plot(polc$identity, polc$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(polc$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(polc$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(polc$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#104 controls (between 0) < 16 cases (within 1)
#Area under the curve: 98.98%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0.60 & roc.df$thresholds < 1.00,]
  #tpp   fpp       thresholds
  # 75.00 0.9615385  0.6156309****
  # 68.75 0.9615385  0.7123951

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #75 3.846154  0.3510188


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(polc$identity, polc$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Binomial Fit")
lines(polc$identity, glm.fit$fitted.values)
text(90, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(polc$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(99, 75, pch = 19)
text(97, 75, paste("Threshold %ID values: 60% \nTrue positives: 75% \nFalse positives: 1%"), adj = c(0,1))
text(40, 20, "Area under the curve: 98.98")


########
#Enolase
########

enolase <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/enolase_identity.csv", header = T, row.names = 1)
enolase <- as.matrix(enolase)

eno_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/enolase_groups.csv", header = T)
levels(as.factor(eno_groups$group_type))

eno <- convert_matrix(enolase, eno_groups, group_column = "group_type", type_column = "allergenic")
head(eno)
levels(eno$comparison_type)
summary(eno$comparison)

eno$identity <- as.numeric(eno$identity)

eno$allergenic <- NA
for (i in 1:nrow(eno)){
  if (eno$comparison_type[i] == "ALAL"){eno$allergenic[i] <- 1}
  else if (eno$comparison_type[i] == "ALNAL"){eno$allergenic[i] <- 0}
  else {eno$allergenic[i] <- NA}
}

eno <- na.omit(eno)

#order data frame by identity values
eno <- eno[order(eno$identity), ]

#binomial fit
glm.fit=glm(eno$allergenic~eno$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#             Estimate  Std. Error  z value   Pr(>|z|)    
#(Intercept)   0.7330604  0.8350214   0.878    0.380
#eno$identity -0.0005686  0.0117715  -0.048    0.961
#AIC: 405

#plot data and model fit results
plot(eno$identity, eno$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(eno$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(eno$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(eno$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#105 controls (allergenic 0) > 210 cases (allergenic 1)
#Area under the curve: 52.93%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$thresholds > 0 & roc.df$thresholds < 1,]
#     tpp       fpp thresholds

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.6, 4, 2.6, 2))

plot(eno$identity, eno$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(eno$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(eno$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
#points(100, 99.38650, pch = 19)
#text(98, 99.38650, paste("Threshold %ID value: 65% \nTrue positives: 99.4% \nFalse positives: 0.1%"), adj = c(0,1))
text(30, 20, "Area under the curve: 52.93%")

##################
#Alpha Lactalbumin
##################

alphalac <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/alpha_lactalbumin_identity.csv", header = T, row.names = 1)
alphalac <- as.matrix(alphalac)

alac_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/alpha_lactalbumin_groups.csv", header = T)

alac <- convert_matrix(alphalac, alac_groups, group_column = "group_type", type_column = "allergenic")
head(alac)
levels(alac$comparison_type)
summary(alac$comparison)

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(alac)
alac[106,] <- c("BovinaeBovinae","ALAL", min(alac$identity))
alac[107,] <- c("BovinaeCaprinae","ALNAL", max(alac$identity))


alac$identity <- as.numeric(alac$identity)

alac$allergenic <- NA
for (i in 1:nrow(alac)){
  if (alac$comparison_type[i] == "ALAL"){alac$allergenic[i] <- 1}
  else if (alac$comparison_type[i] == "ALNAL"){alac$allergenic[i] <- 0}
  else {alac$allergenic[i] <- NA}
}
tail(alac)
alac <- na.omit(alac)

#order data frame by identity values
alac <- alac[order(alac$identity), ]

#fit binomial
glm.fit=glm(alac$allergenic~alac$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)      -12.6239       1.9434   -6.496      8.26e-11 ***
#alac$identity   0.1763       0.0275    6.413      1.43e-10 ***
#AIC: 26.326

#plot data and model fit results
plot(alac$identity, alac$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(alac$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(alac$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(alac$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#55 controls (allergenic 0) < 16 cases (allergenic 1)
#Area under the curve: 92.22%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$tpp > 80 & roc.df$tpp < 100,]
  #   tpp   fpp       thresholds
  #87.50 1.818182  0.8391042***
  #81.25 1.818182  0.8529859

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #93.75   1.818182 0.443971064



#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(alac$identity, alac$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(alac$identity, glm.fit$fitted.values)
text(80, 0.5, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(alac$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(98, 87.50, pch = 19)
text(96, 87.50, paste("Threshold %ID value:83%% \nTrue positives: 87% \nFalse positives: 1.8%"), adj = c(0,1))
text(30, 20, "Area under the curve: 92.22%")


###################
#Beta Lactoglobulin
###################

betalac <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/beta_lactoglobulin_identity.csv", header = T, row.names = 1)
betalac <- as.matrix(betalac)

blac_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/beta_lactoglobulin_groups.csv", header = T)

blac <- convert_matrix(betalac, blac_groups, group_column = "group_type", type_column = "allergenic")
head(blac)
levels(blac$comparison_type)
summary(blac$comparison)

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(blac)
blac[56,] <- c("BovinaeBovinae","ALAL", min(blac$identity))
blac[57,] <- c("BovinaeCaprinae","ALNAL", max(blac$identity))


blac$identity <- as.numeric(blac$identity)

blac$allergenic <- NA
for (i in 1:nrow(blac)){
  if (blac$comparison_type[i] == "ALAL"){blac$allergenic[i] <- 1}
  else if (blac$comparison_type[i] == "ALNAL"){blac$allergenic[i] <- 0}
  else {blac$allergenic[i] <- NA}
}
tail(blac)
blac <- na.omit(blac)

#order data frame by identity values
blac <- blac[order(blac$identity), ]

#fit binomial
glm.fit=glm(blac$allergenic~blac$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)   -12.40521    2.81104  -4.413 1.02e-05 ***
#blac$identity   0.15744    0.03515   4.479 7.49e-06 ***
#AIC: 24.319

#plot data and model fit results
plot(blac$identity, blac$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(blac$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(blac$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(blac$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#29 controls (allergenic 0) < 22 cases (allergenic 1)
#Area under the curve: 92.63%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$tpp > 80 & roc.df$tpp < 100,]
  #   tpp   fpp       thresholds
  #95.45455  3.448276 0.49365062***
  #90.90909  3.448276 0.91056454

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
#     tpp        fpp thresholds
#95.454545   3.448276 0.49365062

#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(blac$identity, blac$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(blac$identity, glm.fit$fitted.values)
text(80, 0.2, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(blac$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(96, 95, pch = 19)
text(94, 95, paste("Threshold %ID value:50% \nTrue positives: 95% \nFalse positives: 3.5%"), adj = c(0,1))
text(30, 20, "Area under the curve: 92.63%")



##############
#Serum albumin
##############

serumalb <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/serum_albumin_identity.csv", header = T, row.names = 1)
serumalb <- as.matrix(serumalb)

sealb_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/serum_albumin_groups.csv", header = T)

sealb <- convert_matrix(serumalb, sealb_groups, group_column = "group_type", type_column = "allergenic")
head(sealb)
levels(sealb$comparison_type)
summary(sealb$comparison)

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(sealb)
sealb[79,] <- c("BovinaeBovinae","ALAL", min(sealb$identity))
sealb[80,] <- c("BovinaeCaprinae","ALNAL", max(sealb$identity))


sealb$identity <- as.numeric(sealb$identity)

sealb$allergenic <- NA
for (i in 1:nrow(sealb)){
  if (sealb$comparison_type[i] == "ALAL"){sealb$allergenic[i] <- 1}
  else if (sealb$comparison_type[i] == "ALNAL"){sealb$allergenic[i] <- 0}
  else {sealb$allergenic[i] <- NA}
}
tail(sealb)
sealb <- na.omit(sealb)

#order data frame by identity values
sealb <- sealb[order(sealb$identity), ]

#fit binomial
glm.fit=glm(sealb$allergenic~sealb$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)    -26.24695    5.90053  -4.448 8.66e-06 ***
#sealb$identity   0.30119    0.06939   4.341 1.42e-05 ***
#AIC: 28.725

#plot data and model fit results
plot(sealb$identity, sealb$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(sealb$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(sealb$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(sealb$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#43 controls (allergenic 0) < 16 cases (allergenic 1)
#Area under the curve: 91.79%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$tpp > 80 & roc.df$tpp < 100,]
  #   tpp   fpp       thresholds
  #93.75   2.325581 0.42598430
  #87.50   2.325581 0.72480923****

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #93.75   2.325581 0.42598430



#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(sealb$identity, sealb$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(sealb$identity, glm.fit$fitted.values)
text(95, 0.2, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(sealb$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(98, 87.50, pch = 19)
text(96, 87.50, paste("Threshold %ID value:72% \nTrue positives: 87.5% \nFalse positives: 2.3%"), adj = c(0,1))
text(30, 20, "Area under the curve: 91.79%")

#################
#Lactotransferrin
#################

lactran <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/lactotransferrin_identity.csv", header = T, row.names = 1)
lactran <- as.matrix(lactran)

lactr_groups <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/lactotransferrin_groups.csv", header = T)

lactr <- convert_matrix(lactran, lactr_groups, group_column = "group_type", type_column = "allergenic")
head(lactr)
levels(lactr$comparison_type)
summary(lactr$comparison)

#add dummy variables at 100,1 and 38,0 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(lactr)
lactr[92,] <- c("BovinaeHuman","ALAL", min(lactr$identity))
lactr[93,] <- c("BovinaeCaprinae","ALNAL", max(lactr$identity))


lactr$identity <- as.numeric(lactr$identity)

lactr$allergenic <- NA
for (i in 1:nrow(lactr)){
  if (lactr$comparison_type[i] == "ALAL"){lactr$allergenic[i] <- 1}
  else if (lactr$comparison_type[i] == "ALNAL"){lactr$allergenic[i] <- 0}
  else {lactr$allergenic[i] <- NA}
}
tail(lactr)
lactr <- na.omit(lactr)

#order data frame by identity values
lactr <- lactr[order(lactr$identity), ]

#fit binomial
glm.fit=glm(lactr$allergenic~lactr$identity, family = binomial)
summary(glm.fit)
#Coefficients:
#                 Estimate    Std. Error  z value    Pr(>|z|)    
#(Intercept)    -30.0723     9.0115  -3.337 0.000847 ***
#lactr$identity   0.3375     0.1004   3.361 0.000775 ***
#AIC: 14.989

#plot data and model fit results
plot(lactr$identity, lactr$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(lactr$identity, glm.fit$fitted.values)

#find ROC and plot it
roc(lactr$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")

#store ROC in a variable to find the optimal threshold value
roc.info <- roc(lactr$allergenic, glm.fit$fitted.values, legacy.axes=TRUE)
#46 controls (allergenic 0) < 11 cases (allergenic 1)
#Area under the curve: 98.02%

#extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df)
tail(roc.df)

roc.df[roc.df$tpp > 80 & roc.df$tpp < 100,]
  #   tpp   fpp       thresholds
  #90.90909 2.173913  0.7281103***
  #81.81818 2.173913  0.7373619

roc.df[roc.df$thresholds >  0.30 & roc.df$thresholds < 0.40,] #tpp and fpp at 35%
  #     tpp        fpp thresholds
  #90.90909 2.173913  0.3900608


#Final plot

par(mfcol=c(1,2), pty = "s", mar = c(4.1, 4, 2.1, 2))

plot(lactr$identity, lactr$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Binomial Fit")
lines(lactr$identity, glm.fit$fitted.values)
text(95, 0.2, paste("P-value: ", round(summary(glm.fit)$coefficients[2,4],4)))

roc(lactr$allergenic, glm.fit$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "ROC curve")
points(98, 90.90909, pch = 19)
text(96, 90.90909, paste("Threshold %ID value:72% \nTrue positives: 91% \nFalse positives: 2.1%"), adj = c(0,1))
text(30, 20, "Area under the curve: 98.02%")