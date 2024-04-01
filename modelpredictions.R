#Tropomyosin

trop <- read.csv("/Users/ngp/Documents/UNL/FARRP/FARRP_R_analyses/tropomyosin_identity.csv", header = T, row.names = 1)
trop <- as.matrix(trop)

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

#obtain id to human data
trop_to_human <- data.frame(type = trop_groups_nf$allergenic, id_to_human = as.numeric(trop_nf[,42]))

trop_to_human$allergenic <- NA
for (i in 1:nrow(trop_to_human)){
  if (trop_to_human$type[i] == "AL"){trop_to_human$allergenic[i] <- 1}
  else if (trop_to_human$type[i] == "NAL"){trop_to_human$allergenic[i] <- 0}
  else {trop_to_human$allergenic[i] <- NA}
}
trop_to_human <- na.omit(trop_to_human)

#add dummy variables at 48,0 and 100,1 to force convergence in the binomial fit (i.e. adding overlap between the identities of the different categories)
nrow(trop_to_human)
head(trop_to_human)
trop_to_human <- rbind(trop_to_human, data.frame(type = c("AL", "NAL"), id_to_human = c(min(trop_to_human$id_to_human), max(trop_to_human$id_to_human)), allergenic = c(0, 1)))

trop_to_human <- trop_to_human[order(as.numeric(trop_to_human$id_to_human)), ]

#fit binomial
glm.fit2=glm(trop_to_human$allergenic~trop_to_human$id_to_human, family = binomial)
summary(glm.fit2) 
#Coefficients:
#                           Estimate    Std. Error  z value   Pr(>|z|)    
#(Intercept)                13.62175    2.45588     5.547     2.91e-08 ***
#trop_to_human$id_to_human  -0.17552    0.03161    -5.552     2.82e-08 ***
#AIC:27.398

#generate probability data based on models 


#this function predicts the probability of being allergenic given by the previously fitted binomial models for values of identity between 0 and 100; glm1 should refer to the model fitted for probability depending on similarity to an allergen, glm2 should refer to the model fitted for probability depending on similarity to a human protein

predict_model <- function(glm1, glm2)  {
#create empty data frame
p = data.frame(rep(0:100, by = 1), rep(NA, 101), rep(NA, 101), rep(NA, 101), rep(NA, 101))
colnames(p) = c("identity", "log_prob_all", "log_prob_not_all", "prob_all", "prob_not_all")
for (i in rep(0:100, by = 1)){
  #use the coefficients given by the models
  p[i,2] = glm1$coefficients[1] + i*glm1$coefficients[2]
  p[i,3] = glm2$coefficients[1] + i*glm2$coefficients[2]
  #transform the log odds given by the previous step into probabilities between 0 and 1
  p[i,4] = exp(p[i,2]) / (1 + exp(p[i,2]))
  p[i,5] = exp(p[i,3]) / (1 + exp(p[i,3]))
}
return(p)
}

predict_model <- function(glm1, glm2)  {
  #create empty data frame
  p = data.frame(rep(0:100, by = 1), rep(NA, 101), rep(NA, 101))
  colnames(p) = c("identity", "prob_all", "prob_not_all")
  for (i in rep(0:100, by = 1)){
    #use the coefficients given by the models
    y1 = glm1$coefficients[1] + i*glm1$coefficients[2]
    y2 = glm2$coefficients[1] + i*glm2$coefficients[2]
    #transform the log odds given by the previous step into probabilities between 0 and 1
    p[i,2] = exp(y1) / (1 + exp(y1))
    p[i,3] = exp(y2) / (1 + exp(y2))
  }
  return(p)
}

trop_predictions = predict_model(glm.fit, glm.fit2)

trop_predictions <- na.omit(trop_predictions)

plot(trop_predictions$log_prob_all, trop_predictions$log_prob_not_all)
plot(trop_predictions$prob_all, trop_predictions$prob_not_all)

x= lm(trop_predictions$identity ~ trop_predictions$log_prob_all*trop_predictions$log_prob_not_all)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                     56.911      6.169   9.226 6.29e-15 ***
#  trop_predictions$prob_all       36.406      4.295   8.476 2.59e-13 ***
#  trop_predictions$prob_not_all  -29.296      5.873  -4.989 2.67e-06 ***
#Multiple R-squared:  0.8223,	Adjusted R-squared:  0.8186 
#F-statistic: 224.4 on 2 and 97 DF,  p-value: < 2.2e-16


x= glm(trop_predictions$identity ~ trop_predictions$prob_all*trop_predictions$prob_not_all, family= quasi)

normal AIC: 791.58
poisson AIC: 1182.2

summary(x)

plot(trop_predictions$prob_all~trop_predictions$identity)
lines(trop_predictions$prob_not_all~trop_predictions$identity)
