"
List of variable names containing the tables with which I built the ROC curves

Animal proteins (divided in ALAL v. ALNAL)
 - tropomyosin_nf
 - tim **
 - bparvalbumin
 - heatshockp **
 - argkin #make sure it doesn't have dummy variables added to the results table
 - caseina #make sure it doesn't have dummy variables added to the results table
 - caseinb #make sure it doesn't have dummy variables added to the results table
 - lipocalin **
 - ovalb #make sure it doesn't have dummy variables added to the results table
 - cyc **
 - eno **
 - alac #make sure it doesn't have dummy variables added to the results table
 - blac #make sure it doesn't have dummy variables added to the results table
 - sealb #make sure it doesn't have dummy variables added to the results table
 - lactr #make sure it doesn't have dummy variables added to the results table

Plant proteins
 - vicilin7SC
 - vicilin7SN
 - twoSalbumin_all
 - eleven_all
 - thau
 - gib #make sure it doesn't have dummy variables added to the results table
 - LTPs
 - prof
 - pr10
 - polc
"

#make a list of all animal and all plant protein tables
allanimal <- list(tropomyosin_nf,tim,bparvalbumin,heatshockp,argkin,caseina,caseinb,lipocalin,ovalb,cyc,eno,alac,blac,sealb,lactr)

allplant <- list(vicilin7SC,vicilin7SN,twoSalbumin_all,eleven_all,thau,gib,LTPs,prof,pr10,polc)

#find the length of the all animal and all plant database
row_counts_animal <- lapply(allanimal, nrow)
sum(unlist(row_counts_animal))
row_counts_plant <- lapply(allplant, nrow)
sum(unlist(row_counts_plant))

#Concatenate all animal and all plant data
allanimal_data <- rbind(tropomyosin_nf,tim,bparvalbumin,heatshockp,argkin,caseina,caseinb,lipocalin,ovalb,cyc,eno,alac,blac,sealb,lactr)
nrow(allanimal_data)
allplant_data <- rbind(vicilin7SC,vicilin7SN,twoSalbumin_all,eleven_all,thau,gib,LTPs,prof,pr10,polc)
nrow(allplant_data)

#order by identity values
allanimal_data <- allanimal_data[order(allanimal_data$identity), ]
allplant_data <- allplant_data[order(allplant_data$identity), ]

#fit binomials
glm.fit.animal=glm(allanimal_data$allergenic~allanimal_data$identity, family = binomial)
summary(glm.fit.animal)
#Coefficients:
#                          Estimate Std. Error z value Pr(>|z|)    
#(Intercept)             -12.046061   0.241523  -49.88   <2e-16 ***
#allanimal_data$identity   0.194272   0.003994   48.64   <2e-16 ***
#AIC: 7538.1

glm.fit.plant=glm(allplant_data$allergenic~allplant_data$identity, family = binomial)
summary(glm.fit.plant)
#Coefficients:
#                         Estimate Std. Error z value Pr(>|z|)    
#(Intercept)            -5.7497188  0.0710917  -80.88   <2e-16 ***
#allplant_data$identity  0.0713987  0.0009666   73.87   <2e-16 ***
#AIC: 32414

#plot data and model fit results
plot(allanimal_data$identity, allanimal_data$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic")
lines(allanimal_data$identity, glm.fit.animal$fitted.values)

plot(allplant_data$identity, allplant_data$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups")
lines(allplant_data$identity, glm.fit.plant$fitted.values)

#find ROC and plot it
roc(allanimal_data$allergenic, glm.fit.animal$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")
#5299 controls (allergenic 0) < 4547 cases (allergenic 1)
#Area under the curve: 92%

roc(allplant_data$allergenic, glm.fit.plant$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %")
#22890 controls (allergenic 0) < 9906 cases (allergenic 1)
#Area under the curve: 77.32%


#store ROC in a variable to find the optimal threshold value
roc.info.animal <- roc(allanimal_data$allergenic, glm.fit.animal$fitted.values, legacy.axes=TRUE)
roc.info.plant <- roc(allplant_data$allergenic, glm.fit.plant$fitted.values, legacy.axes=TRUE)

#extract just the information that we want from that variable.
roc.df.animal <- data.frame(
  tpp=roc.info.animal$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info.animal$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info.animal$thresholds)

roc.df.animal[roc.df.animal$tpp > 80 & roc.df.animal$tpp < 99, ]
#   tpp   fpp       thresholds
# 78.13943 10.077373  0.4994589***
# 76.90785 10.058502  0.5006731

roc.df.plant <- data.frame(
  tpp=roc.info.plant$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info.plant$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info.plant$thresholds)

roc.df.plant [roc.df.plant$tpp > 90 & roc.df.plant$tpp < 95,]

#final plot

par(mfcol=c(2,2), pty = "s", mar = c(4.1, 4, 2.1, 2))


plot(allanimal_data$identity, allanimal_data$allergenic,
     xlab = "% identity", ylab = "0:Non-Allergenic; 1:Allergenic", main = "Animals, Binomial Fit")
lines(allanimal_data$identity, glm.fit.animal$fitted.values)
text(80, 0.2, paste("P-value: ", round(summary(glm.fit.animal)$coefficients[2,4],4)), cex= 0.8)

roc(allanimal_data$allergenic, glm.fit.animal$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "Animals, ROC curve")
points(78, 91, pch = 19)
text(76, 90, paste("Threshold %ID value:50% \nTrue positives: 78% \nFalse positives: 10%"), adj = c(0,1), cex= 0.8)
text(25, 20, "Area under the curve: 92%", cex= 0.8)


plot(allplant_data$identity, allplant_data$allergenic,
     xlab = "% identity", ylab = "0:Between groups; 1:Within groups", main = "Plants, Binomial Fit")
lines(allplant_data$identity, glm.fit.plant$fitted.values)
text(80, 0.2, paste("P-value: ", round(summary(glm.fit.plant)$coefficients[2,4],4)), cex= 0.8)

roc(allplant_data$allergenic, glm.fit.plant$fitted.values, 
    plot = T, legacy.axes = T, percent = T, xlab = "False positive %", ylab ="True positive %", main = "Plants, ROC curve")
#points(78, 90, pch = 19)
#text(76, 90, paste("Threshold %ID value:72% \nTrue positives: 91% \nFalse positives: 2.1%"), adj = c(0,1), cex= 0.8)
text(25, 20, "Area under the curve: 77.32%", cex= 0.8)

