# ------------------- Relationship Survey -------------------
# Load Packages
library(Hmisc) #Describe Function
library(psych) #Multiple Functions for Statistics and Multivariate Analysis
library(GGally) #ggpairs Function
library(ggplot2) #ggplot2 Functions
library(vioplot) #Violin Plot Function
library(corrplot) #Plot Correlations
library(REdaS) #Bartlett's Test of Sphericity
library(psych) #PCA/FA functions
library(factoextra) #PCA Visualizations
library("FactoMineR") #PCA functions
library(ade4) #PCA Visualizations


# Read data
ECR <- read.csv('/Users/tsuerh/Documents/Depaul/DSC424AdvancedDataAnalysis/WK2/Assignment/ECR\ data.csv', header = TRUE, sep = ',')

# Check data
dim(ECR)
head(ECR)

# Check Missing Values
sum(is.na(ECR))

# Descriptive statistics
library(psych)
describe(ECR)

# Check multicollinearity
M <- cor(ECR, method = 'spearman')
M
corrplot(cor(M,method="spearman"), method = "number", type = "lower")


# PCA_Plot functions
#########################################################################

PCA_Plot = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = data.frame(pcaData$rotation, .names = row.names(pcaData$rotation))
  p + geom_text(data=loadings, mapping=aes(x = PC1, y = PC2, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")
}

PCA_Plot_Secondary = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = data.frame(pcaData$rotation, .names = row.names(pcaData$rotation))
  p + geom_text(data=loadings, mapping=aes(x = PC3, y = PC4, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC3", y = "PC4")
}

PCA_Plot_Psyc = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = as.data.frame(unclass(pcaData$loadings))
  s = rep(0, ncol(loadings))
  for (i in 1:ncol(loadings))
  {
    s[i] = 0
    for (j in 1:nrow(loadings))
      s[i] = s[i] + loadings[j, i]^2
    s[i] = sqrt(s[i])
  }
  
  for (i in 1:ncol(loadings))
    loadings[, i] = loadings[, i] / s[i]
  
  loadings$.names = row.names(loadings)
  
  p + geom_text(data=loadings, mapping=aes(x = PC1, y = PC2, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")
}

PCA_Plot_Psyc_Secondary = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = as.data.frame(unclass(pcaData$loadings))
  s = rep(0, ncol(loadings))
  for (i in 1:ncol(loadings))
  {
    s[i] = 0
    for (j in 1:nrow(loadings))
      s[i] = s[i] + loadings[j, i]^2
    s[i] = sqrt(s[i])
  }
  
  for (i in 1:ncol(loadings))
    loadings[, i] = loadings[, i] / s[i]
  
  loadings$.names = row.names(loadings)
  
  print(loadings)
  p + geom_text(data=loadings, mapping=aes(x = PC3, y = PC4, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC3", y = "PC4")
}

##################################################
#PCA/FA
##################################################

# KMO Sampling Adequacy
library(psych)
KMO(ECR)
#Overall MSA =  0.95

# Bartlett's Test of Sphericity
library(REdaS)
bart_spher(ECR)
#p-value < 2.22e-16

# Reliability Analysis test using Cronbach's Alpha
library(psych)
alpha(ECR,check.keys=TRUE)
#Alpha = 0.9

# Parallel Analysis (Horn's parallel analysis)
library(psych)
comp <- fa.parallel(ECR)
comp


##################################################
#Create PCA
p = prcomp(ECR, center=T, scale=T)
p

#Check Scree Plot
plot(p)
abline(1, 0)

#Check PCA Summary Information
summary(p)
print(p)

########################################################
#Check PCA visualizations
plot(p) #Scree Plot
PCA_Plot(p) #PCA_plot1
PCA_Plot_Secondary(p) #PCA_Plot2
biplot(p) #Biplot

#Calculating the Varimax Rotation Loadings manually
rawLoadings = p$rotation %*% diag(p$sdev, nrow(p$rotation), nrow(p$rotation))
print(rawLoadings)
v = varimax(rawLoadings)

p2 = psych::principal(ECR, rotate="varimax", nfactors=4, scores=TRUE)
p2
print(p2$loadings, cutoff=.5, sort=T)

#Calculating scores

scores <- p2$scores

min_score <- min(scores_1)
min_score


max_score <- max(scores_1)
max_score

summary(scores_1)
scores_1 <- scores[,1]
scores_1
mins1 <- min(scores_1)
mins1
max1 <- max(scores_1)
max1
match(mins1, scores_1)
match(max1, scores_1)

scores_2 <- scores[,2]
min2 <- min(scores_2)
min2
max2 <- max(scores_2)
max2

scores_3 <- scores[,3]
min3 <- min(scores_3)
min3
max3 <- max(scores_3)
max3

scores_4 <- scores[,4]
min4 <- min(scores_4)
min4
max4 <- max(scores_4)
max4


#Conducting Factor Analysis

fit = factanal(ECR, 4)
print(fit$loadings, cutoff=.5, sort=T)
summary(fit)
