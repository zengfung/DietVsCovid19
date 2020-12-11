####################################
# load all data
fat = read.table("Fat_Supply_Quantity_Data.csv", header = TRUE, sep = ",")
kcal = read.table("Food_Supply_kcal_Data.csv", header = TRUE, sep = ",")
food = read.table("Food_Supply_Quantity_kg_Data.csv", header = TRUE, sep = ",")
protein = read.table("Protein_Supply_Quantity_Data.csv", header = TRUE, sep = ",")

####################################
# renaming variables
dfAsString = function(x){
  deparse(substitute(x))
}
colRename = function(x, name) {
  names(x) = c("Country", "Alc", "AnimalProd", "AnimalFats",
               "AqProd", "Cereal", "Eggs", "Seafood", "Fruits",
               "Meat", "Misc", "Milk", "Offals", "Oilcrops",
               "Pulses", "Spices", "StRt", "Stim", "SugarCrops",
               "Sweeteners", "Treenuts", "VegProd", "VegOil", "Veg", 
               "Obesity", "Undernourished", "Cfm", "Deaths", "Rec", 
               "Active", "Pop", "Unit")
  names(x)[2:24] = paste(names(x)[2:24], name, sep = "-")
  return(x)
}
fat = colRename(fat, dfAsString(fat))
kcal = colRename(kcal, dfAsString(kcal))
food = colRename(food, dfAsString(food))
protein = colRename(protein, dfAsString(protein))

####################################
## exploratory data analysis
# check correlation between measurements of same food categories
library(GGally)
for (i in 2:24){
  tempDf = as.data.frame(cbind(fat = fat[,i], kcal = kcal[,i], food = food[,i], protein = protein[,i]))
  title = strsplit(names(fat)[i], "-")[[1]]
  print(ggpairs(tempDf, title = title,
                upper = list(continuous = "points"), lower = list(continuous = "cor"), 
                progress = FALSE))
  print(colSums(cor(tempDf)))
}

# remove unnecessary columns
redDf = cbind(fat[c(1,2,4,6,9,11,12,18,19,20,23)],
              kcal[c(3,5,6,7,8,10,12,14:17,20:24)],
              food[c(2:4,13,18,19,20,22,24)],
              protein[c(4,9,11:31)])
dim(redDf)

# remove X variables with unproper distribution
# quick view of histograms for all X variables
for (i in 2:52){
  hist(redDf[,i], xlab = names(redDf)[i], main = names(redDf)[i])
}

# remove unnecessary columns
drops = c("Alc-fat", "SugarCrops-fat", "Sweeteners-fat", "AqProd-kcal", "Sweeteners-food", "SugarCrops-protein")
redDf = redDf[,!(names(redDf) %in% drops)]

# remove X variables with extremely high VIF
VIF = diag(solve(cor(redDf[,2:46])))
ord = order(VIF, decreasing = TRUE)
VIF[ord]

# remove unnecessary columns
drops = c("AnimalFats-food", "Veg-food", "AnimalProd-kcal", "VegProd-kcal")
redDf = redDf[,!(names(redDf) %in% drops)]
dim(redDf)

# checking correlation between Y variables
cor(na.omit(redDf)[,45:48])
dim(na.omit(redDf))

# best correlation between X and Y variables
bestCorDf = data.frame(Ord = 1:10)
for (i in 45:48){
  corrXY = cor(na.omit(redDf[,c(2:42,i)]))[,42]
  ord = order(abs(corrXY), decreasing = TRUE)
  bestCorDf = cbind(bestCorDf, 
                    names(corrXY[ord[2:11]]), 
                    corrXY[ord[2:11]])
}
names(bestCorDf) = c("Ord", "varCfm", "corCfm", "varD", "corD", "varRec", "corRec", "varAct", "corAct")
row.names(bestCorDf) = NULL
bestCorDf

# distribution of Y variables
hist(redDf$Cfm, xlab = "Cfm", main = "Confirmed Cases")
hist(redDf$Death, xlab = "Death", main = "Death")
hist(redDf$Rec, xlab = "Rec", main = "Recovered")
hist(redDf$Active, xlab = "Active", main = "Active")

# box cox transformation
library(MASS)
fit = lm(Cfm ~ ., data = redDf[c(2:42, 45)]); boxcox(fit)
fit = lm(I(Deaths+.00001) ~ ., data = redDf[c(2:42, 46)]); boxcox(fit)
fit = lm(I(Rec+.00001) ~ ., data = redDf[c(2:42, 47)]); boxcox(fit)
fit = lm(I(Active+.00001) ~ ., data = redDf[c(2:42, 48)]); boxcox(fit)

# log transformation of Y variables
hist(log(fat$Cfm), xlab = "log(Cfm)", main = "log(Confirmed Cases)")
hist(log(fat$Death), xlab = "log(Death)", main = "log(Death)")
hist(log(fat$Rec), xlab = "log(Rec)", main = "log(Recovered)")
hist(log(fat$Active), xlab = "log(Active)", main = "log(Active)")

##########################################
## model fitting
# omit NA's and form 4 separate data frames
cfmDf = na.omit(redDf[,c(1:42, 45)])
deathDf = na.omit(redDf[,c(1:42,46)])
recDf = na.omit(redDf[,c(1:42,47)])
actDf = na.omit(redDf[,c(1:42,48)])

# check dimensions of each dataframe
dim(cfmDf)
dim(deathDf)
dim(recDf)
dim(actDf)

# transformation of Y variables
cfmDf$Cfm = log(cfmDf$Cfm)
deathDf$Deaths = log(deathDf$Deaths + 1e-12)
recDf$Rec = log(recDf$Rec + 1e-12)
actDf$Active = log(actDf$Active + 1e-12)

# split data into train and validation
set.seed(10)
ind1 = sample(1:dim(cfmDf)[1], size = 0.7*dim(cfmDf)[1])
ind2 = sample(1:dim(actDf)[1], size = 0.7*dim(actDf)[1])
# train
cfmDf.train = cfmDf[ind1,]
deathDf.train = deathDf[ind1,]
recDf.train = recDf[ind1,]
actDf.train = actDf[ind2,]
#valid
cfmDf.valid = cfmDf[-ind1,]
deathDf.valid = deathDf[-ind1,]
recDf.valid = recDf[-ind1,]
actDf.valid = actDf[-ind2,]

# check distribution of response variables in train and validation sets
hist(cfmDf.train$Cfm, main = "Train")
hist(cfmDf.valid$Cfm, main = "Valid")
hist(deathDf.train$Deaths, main = "Train")
hist(deathDf.valid$Deaths, main = "Valid")
hist(recDf.train$Rec, main = "Train")
hist(recDf.valid$Rec, main = "Valid")
hist(actDf.train$Active, main = "Train")
hist(actDf.valid$Active, main = "Valid")

# model fitting for confirmed cases
cfm.fit1 = lm(Cfm ~ 1, data = cfmDf.train[,-1])

# fitting AIC model
cfm.fit2 = stepAIC(cfm.fit1, 
                   scope = list(upper = lm(Cfm ~ ., data = cfmDf.train[,-1]),
                                lower = cfm.fit1),
                   direction = "both", trace = 0, k = 2)
summary(cfm.fit2)
cfm.fit4 = stepAIC(cfm.fit1,
                   scope = list(upper = lm(Cfm ~ (`Misc-protein` + `Sweeteners-kcal` + `Meat-kcal` + 
                                                    `Pulses-kcal` + `Stim-protein` + `Oilcrops-kcal` + `Fruits-protein` + 
                                                    `Eggs-kcal` + `Alc-food` + `Treenuts-kcal` + `Spices-kcal` + 
                                                    `Sweeteners-protein`)^2, 
                                           data = cfmDf.train[, -1]),
                                lower= cfm.fit1),
                   direction = "both", trace = 0, k = 2)
summary(cfm.fit4)
# fitting BIC model
cfm.fit3 = stepAIC(cfm.fit1, 
                   scope = list(upper = lm(Cfm ~ .^2, data = cfmDf.train[,-1]),
                                lower = cfm.fit1),
                   direction = "both", trace = 0, k = log(nrow(cfmDf.train)))
summary(cfm.fit3)

# model fitting for death cases
death.fit1 = lm(Deaths ~ 1, data = deathDf.train[,-1])
# fitting AIC model
death.fit2 = stepAIC(death.fit1, 
                     scope = list(upper = lm(Deaths ~ ., data = deathDf.train[,-1]),
                                  lower = death.fit1),
                     direction = "both", trace = 0, k = 2)
summary(death.fit2)
death.fit4 = stepAIC(death.fit1, 
                     scope = list(upper = lm(formula = Deaths ~ (`Misc-protein` + `Treenuts-kcal` + `Meat-kcal` + 
                                                                   `Eggs-kcal` + `Spices-kcal` + `Pulses-protein` + `Stim-protein` + 
                                                                   `Oilcrops-kcal` + `Alc-food` + `Veg-protein` + `Veg-kcal` + 
                                                                   `Misc-fat`)^2,
                                             data = deathDf.train[, -1]),
                                  lower = death.fit1),
                     direction = "both", trace = 0, k = 2)
summary(death.fit4)
# fitting BIC model
death.fit3 = stepAIC(death.fit1, 
                     scope = list(upper = lm(Deaths ~ .^2, data = deathDf.train[,-1]),
                                  lower = death.fit1),
                     direction = "both", trace = 0, k = log(nrow(deathDf.train)))
summary(death.fit3)

# model fitting for recovery cases
rec.fit1 = lm(Rec ~ 1, data = recDf.train[,-1])
# fitting AIC model
rec.fit2 = stepAIC(rec.fit1, 
                   scope = list(upper = lm(Rec ~ .^2, data = recDf.train[,-1]),
                                lower = rec.fit1),
                   direction = "both", trace = 0, k = 2)
summary(rec.fit2)
# fitting BIC model
rec.fit3 = stepAIC(rec.fit1, 
                   scope = list(upper = lm(Rec ~ .^2, data = recDf.train[,-1]),
                                lower = rec.fit1),
                   direction = "both", trace = 0, k = log(nrow(recDf.train)))
summary(rec.fit3)

# model fitting for active cases
act.fit1 = lm(Active ~ 1, data = actDf.train[,-1])
# fitting AIC model
act.fit2 = stepAIC(act.fit1, 
                   scope = list(upper = lm(Active ~ ., data = actDf.train[,-1]),
                                lower = act.fit1),
                   direction = "both", trace = 0, k = 2)
summary(act.fit2)
act.fit4 = stepAIC(act.fit1,
                   scope = list(upper = lm(formula = Active ~ (`Oilcrops-kcal` + `Misc-protein` + `Meat-kcal` + 
                                                                 `Treenuts-protein` + `Pulses-kcal` + `Alc-food` + `Sweeteners-protein`)^2, 
                                           data = actDf.train[, -1]),
                                lower = act.fit1),
                   direction = "both", trace = 0, k = 2)
summary(act.fit4)
# fitting BIC model
act.fit3 = stepAIC(act.fit1, 
                   scope = list(upper = lm(Active ~ .^2, data = actDf.train[,-1]),
                                lower = act.fit1),
                   direction = "both", trace = 0, k = log(nrow(actDf.train)))
summary(act.fit3)

########################################################
## model diagnostics
# define functions
# compute sse, mse, R2, adjusted R2, Cp, Pressp statistics
modelDiag = function(obj, sigma2){
  n = length(obj$residuals)
  p = length(obj$coefficients)
  h = influence(obj)$hat
  sse = sum(obj$residuals^2)
  mse = sse/obj$df.residual
  r2 = summary(obj)$r.squared
  adjr2 = summary(obj)$adj.r.squared
  cp = sse/sigma2 - (n - 2*p)
  pressp = sum((obj$residuals/(1-h))^2)
  return (c(p = p, SSE = sse, MSE = mse, R2 = r2, Adj.R2 = adjr2, Cp = cp, PRESSp = pressp))
}

getMspe = function(fit, data){
  actual = data[,ncol(data)]
  pred = predict(fit, data)
  res = actual - pred
  m = nrow(data)
  return (c(MSPE = sum(res^2)/m))
}

modelValid = function(data, fit1, fit2){
  avg.sse1 = sum(fit1$residuals^2)/nrow(fit1$model)
  avg.sse2 = sum(fit2$residuals^2)/nrow(fit2$model)
  return (data.frame(AIC.valid = c(getMspe(fit1, data), `SSE/n` = avg.sse1), 
                     BIC.valid = c(getMspe(fit2, data), `SSE/n` = avg.sse2)))
}

# internal and external validation for confirmed cases
# confirmed cases
cfm.fitfull = lm(formula = Cfm ~ `Misc-protein` + `Oilcrops-kcal` + `Stim-protein` + 
                   `Fruits-protein` + `Pulses-kcal` + `Meat-kcal` + `Eggs-kcal` + 
                   `Sweeteners-kcal` + `Alc-food` + `Treenuts-kcal` + `Oilcrops-kcal`:`Eggs-kcal` + 
                   `Eggs-kcal`:`Sweeteners-kcal` + `Meat-kcal`:`Alc-food` + 
                   `Fruits-protein`:`Pulses-kcal` + `Misc-protein`:`Alc-food`, data = cfmDf.train[,-1])
cfm.fit4.diag = modelDiag(cfm.fit4, sum(cfm.fitfull$residuals^2)/cfm.fitfull$df.residual) # train model
cfm.fit3.diag = modelDiag(cfm.fit3, sum(cfm.fitfull$residuals^2)/cfm.fitfull$df.residual)
cfm.diagComp = data.frame(AIC.train = cfm.fit4.diag, BIC.train = cfm.fit3.diag)
cfm.diagComp

cfm.diagComp.v = modelValid(cfmDf.valid[,-1], cfm.fit4, cfm.fit3) # valid model
cfm.diagComp.v

# internal and external validation for death cases
# death cases
death.fitfull = lm(Deaths ~ `Misc-protein` + `Meat-kcal` + `Eggs-kcal` + 
                     `Spices-kcal` + `Pulses-protein` + `Stim-protein` + `Oilcrops-kcal` + 
                     `Treenuts-kcal` + `Veg-protein` + `Misc-fat` + `Alc-food` + 
                     `Eggs-kcal`:`Stim-protein` + `Eggs-kcal`:`Pulses-protein` + 
                     `Meat-kcal`:`Pulses-protein` + `Meat-kcal`:`Treenuts-kcal` + 
                     `Veg-protein`:`Misc-fat` + `Meat-kcal`:`Veg-protein` + `Misc-protein`:`Misc-fat` + 
                     `Eggs-kcal`:`Oilcrops-kcal` + `Spices-kcal`:`Treenuts-kcal` + 
                     `Meat-kcal`:`Misc-fat` + `Meat-kcal`:`Alc-food` + `Treenuts-kcal`:`Alc-food` + 
                     `Eggs-kcal`:`Veg-protein` + `Pulses-protein`:`Oilcrops-kcal` + 
                     `Pulses-protein`:`Treenuts-kcal` + `Misc-fat`:`Alc-food` + 
                     `Veg-protein`:`Alc-food` + `Misc-protein`:`Oilcrops-kcal` + 
                     `Stim-protein`:`Alc-food` + `Pulses-protein`:`Misc-fat` + 
                     `Stim-protein`:`Veg-protein` + `Eggs-kcal`:`Misc-fat` + `VegOil-fat`, deathDf.train[, -1])
death.fit4.diag = modelDiag(death.fit4, sum(death.fitfull$residuals^2)/death.fitfull$df.residual) # train model
death.fit3.diag = modelDiag(death.fit3, sum(death.fitfull$residuals^2)/death.fitfull$df.residual)
death.diagComp = data.frame(AIC.train = death.fit4.diag, BIC.train = death.fit3.diag)
death.diagComp

death.diagComp.v = modelValid(deathDf.valid[,-1], death.fit4, death.fit3) # valid model
death.diagComp.v

# internal and external validation for recovered cases
# recovery cases
rec.fitfull = rec.fit2
rec.fit2.diag = modelDiag(rec.fit2, sum(rec.fitfull$residuals^2)/rec.fitfull$df.residual) # train model
rec.fit3.diag = modelDiag(rec.fit3, sum(rec.fitfull$residuals^2)/rec.fitfull$df.residual)
rec.diagComp = data.frame(AIC.train = rec.fit2.diag, BIC.train = rec.fit3.diag)
rec.diagComp

rec.diagComp.v = modelValid(recDf.valid[,-1], rec.fit2, rec.fit3) # valid model
rec.diagComp.v

# internal and external validation for activ cases
# active cases
act.fitfull = lm(formula = Active ~ `Misc-protein` + `Offals-protein` + `Meat-kcal` + 
                   `Eggs-kcal` + `Spices-protein` + `Stim-fat` + `Meat-kcal`:`Spices-protein` + 
                   `Misc-protein`:`Stim-fat` + `Misc-protein`:`Meat-kcal` + `Offals-food`, 
                 data = actDf.train[, -1])
act.fit4.diag = modelDiag(act.fit4, sum(act.fitfull$residuals^2)/act.fitfull$df.residual) # train model
act.fit3.diag = modelDiag(act.fit3, sum(act.fitfull$residuals^2)/act.fitfull$df.residual)
act.diagComp = data.frame(AIC.train = act.fit4.diag, BIC.train = act.fit3.diag)
act.diagComp

act.diagComp.v = modelValid(actDf.valid[,-1], act.fit4, act.fit3) # valid model
act.diagComp.v

## checking for outlying cases
# function to calculate Cook's distance
getCooksDist = function(obj){
  e = obj$residuals
  p = length(obj$coefficients)
  mse = sum(e^2)/obj$df.residual
  h = influence(obj)$hat
  return (c(e^2*h/(p*mse*(1-h)^2)))
}

# model diagnostics for confirmed cases
# residual analysis
plot(cfm.fit4, which = 1:2, main = "AIC Model")
plot(cfm.fit3, which = 1:2, main = "BIC Model")
# outlying cases
plot(cfm.fit4, which = 4:6)
plot(cfm.fit3, which = 4:6)
# calculate cook's distance and determine outliers
cooks = getCooksDist(cfm.fit3)
cfmOutliers = cbind(cfmDf.train[which(cooks > 4/(cfm.fit3$df.residual)),c(1,43)], CooksDist = cooks[which(cooks > 4/(cfm.fit3$df.residual))])
cfmOutliers[order(cfmOutliers$CooksDist, decreasing = TRUE),]
# examine influence of possible outlier
cfm.fit3.rmOutlier = lm(formula = cfm.fit3$call$formula, data = cfmDf.train[-c(143,150,161,86),])
per.change = abs((cfm.fit3$fitted.values - predict.lm(cfm.fit3.rmOutlier, cfmDf.train))/cfm.fit3$fitted.values)*100
summary(per.change)
# remove outliers
cfmDf.train.rmOutlier = cfmDf.train[-c(143,150,161,86),]

# model diagnostics for death cases
# residual analysis
plot(death.fit4, which = 1:2, main = "AIC Model")
plot(death.fit3, which = 1:2, main = "BIC Model")
# outlying cases
plot(death.fit4, which = 4:6)
plot(death.fit3, which = 4:6)
# calculate cook's distance and determine outliers
cooks = getCooksDist(death.fit3)
deathOutliers = cbind(deathDf.train[which(cooks > 4/(death.fit3$df.residual)),c(1,43)], CooksDist = cooks[which(cooks > 4/(death.fit3$df.residual))])
deathOutliers[order(deathOutliers$CooksDist, decreasing = TRUE),]
# examine influence of possible outlier
death.fit3.rmOutlier = lm(formula = death.fit3$call$formula, data = deathDf.train[-c(102,165,86,25,42,151),])
per.change = abs((death.fit3$fitted.values - predict.lm(death.fit3.rmOutlier, deathDf.train))/death.fit3$fitted.values)*100
summary(per.change)
# remove outlier
deathDf.train.rmOutlier = deathDf.train[-c(102,165,86,25,42,151),]

# model diagnostics for recovered cases
# residual analysis
plot(rec.fit2, which = 1:2, main = "AIC Model")
plot(rec.fit3, which = 1:2, main = "BIC Model")
# outlying cases
plot(rec.fit2, which = 4:6)
plot(rec.fit3, which = 4:6)
# calculate cook's distance and determine outliers
cooks = getCooksDist(rec.fit3)
recOutliers = cbind(recDf.train[which(cooks > 4/(rec.fit3$df.residual)),c(1,43)], CooksDist = cooks[which(cooks > 4/(rec.fit3$df.residual))])
recOutliers[order(recOutliers$CooksDist, decreasing = TRUE),]
# examine influence of possible outlier
rec.fit3.rmOutlier = lm(formula = rec.fit3$call$formula, data = recDf.train[-c(165,15,146,143),])
per.change = abs((rec.fit3$fitted.values - predict.lm(rec.fit3.rmOutlier, recDf.train))/rec.fit3$fitted.values)*100
summary(per.change)
# remove outliers from dataset
recDf.train.rmOutlier = recDf.train[-c(165,15,146,143),]

# model diagnostics for active cases
# residual analysis
plot(act.fit4, which = 1:2, main = "AIC Model")
plot(act.fit3, which = 1:2, main = "BIC Model")
# outlying cases
plot(act.fit4, which = 4:6)
plot(act.fit3, which = 4:6)
# calculate cook's distance and determine outliers
cooks = getCooksDist(act.fit3)
actOutliers = cbind(actDf.train[which(cooks > 4/(act.fit3$df.residual)),c(1,43)], CooksDist = cooks[which(cooks > 4/(act.fit3$df.residual))])
actOutliers[order(actOutliers$CooksDist, decreasing = TRUE),]
# examine influence of possible outlier
act.fit3.rmOutlier = lm(formula = act.fit3$call$formula, data = actDf.train[-c(165,143),])
per.change = abs((act.fit3$fitted.values - predict.lm(act.fit3.rmOutlier, actDf.train))/act.fit3$fitted.values)*100
summary(per.change)

#########################################################################
## final models: refit models without outliers with BIC criterion
# confirmed cases
cfmDf.rmOutlier = rbind(cfmDf.train.rmOutlier, cfmDf.valid)
cfm.fit.final = lm(cfm.fit3$call$formula, data = cfmDf.rmOutlier)
summary(cfm.fit.final)
anova(cfm.fit.final)

# death cases
deathDf.rmOutlier = rbind(deathDf.train.rmOutlier, deathDf.valid)
death.fit.final = lm(death.fit3$call$formula, data = deathDf.rmOutlier)
summary(death.fit.final)
anova(death.fit.final)

# recovered cases
recDf.rmOutlier = rbind(recDf.train.rmOutlier, recDf.valid)
rec.fit.final = lm(rec.fit3$call$formula, data = recDf.rmOutlier)
summary(rec.fit.final)
anova(rec.fit.final)

# active cases
actDf.rmOutlier = rbind(actDf.train, actDf.valid)
act.fit.final = lm(act.fit4$call$formula, data = actDf.rmOutlier)
summary(act.fit.final)
anova(act.fit.final)
