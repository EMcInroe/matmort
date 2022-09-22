library(readr)
setwd("~/Statistics/STAT 581")
mortstat <- read_csv("mortstat.csv")
View(mortstat)

####################
### New Variables ##
####################
# Numeric Values to Percentages
mortstat$fempop<-((mortstat$Female)/100)*mortstat$Pop
mortstat$malepop<-mortstat$Pop-mortstat$fempop
mortstat$gradm_per<-((mortstat$GRadM)*1000)/mortstat$malepop*100
mortstat$bachm_per<-((mortstat$BachM)*1000)/mortstat$malepop*100
mortstat$AM_per<-((mortstat$AM)*1000)/mortstat$malepop*100
mortstat$HSM_per<-((mortstat$HSM)*1000)/mortstat$malepop*100
mortstat$GradF_per<-((mortstat$GradF)*1000)/mortstat$fempop*100
mortstat$BachF_per<-((mortstat$BachF)*1000)/mortstat$fempop*100
mortstat$AF_per<-((mortstat$AF)*1000)/mortstat$fempop*100
mortstat$HSF_per<-((mortstat$HSF)*1000)/mortstat$fempop*100

## New dataframe without numeric education data
mortclean <- mortstat[ , ! names(mortstat) %in% c("GRadM", "GradF", "BachM","BachF","AM",
                                                  "AF","HSM","HSF", "malepop","fempop",
                                                  "Graduate","bachelor","Associates","HighSchool")]
View(mortclean)

# Mortality as categorical data
mortclean$mortcat <- ifelse(mortclean$MatMort <  18, "Low",
                    ifelse(mortclean$MatMort > 26, "High", "Average"))

library(dplyr)
library(caret)
mortclean <- mortclean %>% relocate(mortcat, .before = MatMort)

mortclean.nonzero<-mortclean[-c(1:3),]
View(mortclean.nonzero)

#####################
###  EDA         ###
####################
summary(mortstat)
names(mortstat)

##### Correlation
cor(mortstat[,-1])
library(ggplot2)
library(GGally)
ggcorr(mortstat[,-1,-25])

### Visualizations
hist(mortstat$Pop)
hist(mortstat$Asian)
hist(mortstat$EdExpen)
plot(mortclean.nonzero$EdExpen, mortclean.nonzero$MatMort, xlab="Education Expenditure ($/student)",
     ylab="Maternal Mortality Rate", pch = 21, bg="darksalmon")

mortsum<-summary(mortstat$MatMort)

ggplot(mortstat, aes(x=MatMort, y=..density..))+geom_histogram(fill="lightblue", color="black", binwidth=4)+
  ggtitle("Distribution of Maternal Mortality") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_density()

ggplot(mortstat, aes(x=EdExpen, y=..density..))+geom_histogram(fill="lightblue", color="black", binwidth=1000)+
  ggtitle("Distribution of Education Expenditures (Amount Per Student)") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_density()


with(mortclean, plot(BachF_per, MatMort,pch = 21 , bg ="darkorchid3", xlab = "Percent Bachelor Degrees", 
        ylab = "Maternal Mortality (Per 100,000 births)"))
par(new = TRUE)
with(mortclean, plot(bachm_per, MatMort, xaxt = "n", yaxt = "n", ylab = "",
                 xlab = "", pch = 24, bg = "darkslategray4"))
legend("top", inset = c(0, -0.15), pch = c(21, 24), pt.bg = c("darkorchid3",
                                                              "darkslategray4"), 
      legend = c("Female", "Male"), ncol = 2, xpd = NA,
      cex = 0.9)

with(mortclean, plot(AF_per, MatMort,pch = 21 , bg ="darkorchid3", xlab = "Percent Associate Degrees", 
                     ylab = "Maternal Mortality (Per 100,000 births)"))
par(new = TRUE)
with(mortclean, plot(AM_per, MatMort, xaxt = "n", yaxt = "n", ylab = "",
                     xlab = "", pch = 24, bg = "darkslategray4"))
legend("top", inset = c(0, -0.15), pch = c(21, 24), pt.bg = c("darkorchid3",
                                                              "darkslategray4"), 
       legend = c("Female", "Male"), ncol = 2, xpd = NA,
       cex = 0.9)

with(mortclean, plot(HSF_per, MatMort,pch = 21 , bg ="darkorchid3", xlab = "Percent High School Diplomas", 
                     ylab = "Maternal Mortality (Per 100,000 births)"))
par(new = TRUE)
with(mortclean, plot(HSM_per, MatMort, xaxt = "n", yaxt = "n", ylab = "",
                     xlab = "", pch = 24, bg = "darkslategray4"))
legend("top", inset = c(0, -0.15), pch = c(21, 24), pt.bg = c("darkorchid3",
                                                              "darkslategray4"), 
       legend = c("Female", "Male"), ncol = 2, xpd = NA,
       cex = 0.9)

#### Transformations
#Transform Data
mort_numeric <- mortclean.nonzero[,-1]
pp_zscore <- preProcess(mort_numeric, method = c("center", "scale"))
mort_numeric_zscore <- predict(pp_zscore, mort_numeric)
head(mort_numeric_zscore)

############ Least Squares Regression
model1<-lm(MatMort~.-State-mortcat, data=mortclean)
summary(model1)

plot(mortclean$EdExpen, mortclean$MatMort)
model2<-lm(MatMort~EdExpen, data=mort_numeric_zscore)
summary(model2)
par(mfrow=c(2,2))
plot(model2)
par(mfrow=c(1,1))
plot(fitted(model2), model2$residuals)
abline(0,0)
hist(model2$residuals, main="Distribution of Residuals (Expenditure)")

mod.nonzero<-lm(MatMort~.-State-mortcat, data=mortclean.nonzero)
summary(mod.nonzero)

#### transformed LS
model3<-lm(MatMort~.-mortcat, data=mort_numeric_zscore)
summary(model3)
#produce residual vs. fitted plot
plot(fitted(model3), model3$residuals)

#add a horizontal line at 0 
abline(0,0)
hist(model3$residuals)

###########  Best Subset Selection
library(leaps)
reg.fit<-regsubsets(MatMort~.-State-mortcat, data=mortclean, nvmax=18)
reg.sum<-summary(reg.fit)
names(reg.sum)
reg.sum$rsq
plot(reg.sum$rss)
which.min(reg.sum$bic)
coef(reg.fit, 10)

model4<-lm(MatMort~AA+Graduate+Associates, data=mortstat)
summary(model4)

reg.fit.non<-regsubsets(MatMort~.-State-mortcat, data=mortclean.nonzero, nvmax=18)
reg.sum.non<-summary(reg.fit.non)
names(reg.sum.non)
reg.sum.non$rsq
plot(reg.sum.non$rss)
which.min(reg.sum.non$bic)
coef(reg.fit.non, 10)
coef(reg.fit.non, 3)

## transformed Best
reg.fit<-regsubsets(MatMort~.-mortcat, data=mort_numeric_zscore, nvmax=18)
reg.sum<-summary(reg.fit)
names(reg.sum)
reg.sum$rsq
plot(reg.sum$rss, ylab="Rooted Sum of Squares")
which.min(reg.sum$bic)
coef(reg.fit, 3)

model5<-lm(MatMort~AA+GradF_per+AF_per, data=mort_numeric_zscore)
summary(model5)

plot()

############# Ridge Regression
par(mfrow=c(1,1))
x<-model.matrix(MatMort~.-mortcat, mort_numeric_zscore)
y<-mort_numeric_zscore$MatMort

library(glmnet)
cv.out<-cv.glmnet(x,y,alpha=0)
bestlam<-cv.out$lambda.min
bestlam
ridge.mod<-glmnet(x,y,alpha=0, lambda=bestlam)
predict(ridge.mod,s=bestlam, type="coefficients")
plot(cv.out)

ridge.pred<-predict(ridge.mod, s=bestlam, newx=x[test,])
mean((ridge.pred-y[test,])^2)

########### Linear Discriminant Analysis
library(MASS)
x<-model.matrix(mortcat~AA+gradm_per+AF_per, mortclean.nonzero)
y<-mortclean.nonzero$mortcat
set.seed(208)
test<-sample(1:nrow(x),nrow(x)/5)
train<-(-test)
mortclean_train<-mortclean.nonzero[train,]
mortclean_test<-mortclean.nonzero[test,]
lda.fit<-lda(mortcat~AA+gradm_per+AF_per, data=mortclean, subset=train)
lda.fit
lda.pred<-predict(lda.fit, mortclean_test)
lda.mort<-lda.pred$class
table(lda.mort,mortclean_test$mortcat)
mean(lda.mort==mortclean_test$mortcat)

# with education expenditure
x<-model.matrix(mortcat~EdExpen, mortclean)
y<-mortclean$mortcat
set.seed(208)
test<-sample(1:nrow(x),nrow(x)/5)
train<-(-test)
mortclean_train<-mortclean[train,]
mortclean_test<-mortclean[test,]
lda.fit<-lda(mortcat~EdExpen, data=mortclean, subset=train)
lda.fit
lda.pred<-predict(lda.fit, mortclean_test)
lda.mort<-lda.pred$class
table(lda.mort,mortclean_test$mortcat)
mean(lda.mort==mortclean_test$mortcat)

##transformed
x<-model.matrix(mortcat~AA+gradm_per+AF_per, mort_numeric_zscore)
y<-mort_numeric_zscore$mortcat
set.seed(208)
test<-sample(1:nrow(x),nrow(x)/5)
train<-(-test)
mort_numeric_train<-mort_numeric_zscore[train,]
mort_numeric_test<-mort_numeric_zscore[test,]
lda.fit<-lda(mortcat~AA+gradm_per+AF_per, data=mort_numeric_zscore, subset=train)
lda.fit
lda.pred<-predict(lda.fit, mort_numeric_test)
lda.mort<-lda.pred$class
table(lda.mort,mort_numeric_test$mortcat)
mean(lda.mort==mort_numeric_test$mortcat)

###########################
### Clustering       ####
##########################
x1<-matrix(mort_numeric_zscore[,-1], ncol=19)
km.out<-kmeans(mort_numeric_zscore[,-1],3,nstart=1)
par(mfrow=c(1,2))
plot(x,col=(km.out$cluster+1),main="K-means Clustering Results with K=2",pch=20, cex=2)
km.out

pr.out<-prcomp(mort_numeric_zscore[,-1], scale=TRUE)
summary(pr.out)
plot(pr.out)
pve<-100*pr.out$sdev^2/sum(pr.out$sdev^2)
pve
par(mfrow=c(1,1))
plot(pve, type="o", ylab="PVE")

#### Profile plot

ggplot(mortclean.nonzero, aes(x=gradm_per, y=AF_per, color=mortcat, shape=mortcat)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)



