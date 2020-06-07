#install.packages("cluster", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(cluster)
library(data.table)#Data. table is an extension of data. frame package in R. It is widely used for fast aggregation of large datasets,
library(Hmisc)#data analysis funs
library(dplyr)
library(tidyverse)
library(ggplot2)
library(plotly)
library(GGally)
library(ggthemes)
library(psych)
library(relaimpo)
library(e1071)
library(corrplot)
library(factoextra)
library(fpc)

##Loading Dataset into dataframe
Protein <- read.delim("C:/Alok/OneDrive/Rutgers_MITA/Semester2/MVA/MVA_Midterm_FinalFolder/Protein_Consumption.csv", header = TRUE,sep = ",")
Prot_DS <- Protein
#View(Prot_DS)
names(Prot_DS)
#Renaming 1st column to Country for simplicity
names(Prot_DS)[names(Prot_DS) == "ï..Country"] <- "Country"
attach(Prot_DS)

#*********************************
##Basic EDA
head(Prot_DS)
dim(Prot_DS)
#This dataset has 25 rows and 11 columns
#To check data types
str(Prot_DS)
glimpse(Prot_DS)
summary(Prot_DS)
grep('NA',Prot_DS)
#There are no NULL values in our dataset
#We will drop column 1 as it's categorical
Prot_DS.num <- Prot_DS[,-1]
Prot_DS.num
#I am not passing 'Total' Column for Techniques and correlation matrix etc as it is just an addition
#of all the values also, it will definitely have correlation with other columns and so may create multicollinearity issues

Prot_DS.num<- Prot_DS.num[,-10]
dim(Prot_DS.num)
head(Prot_DS.num)

# Computing the means of each variable in data frame 
colMeans(Prot_DS.num)
#avg found

# Covariance matrix without total column
cov(Prot_DS.num)
# Finding correlation -Correlation matrix takes units out and gives normalized values 
cor(Prot_DS.num)
#To check if variables are Normal individually for Milk and Fish just to get an idea
#If it is normal it shud show straight line
qqnorm(Prot_DS.num[,"Milk"], main = "Milk")
qqline(Prot_DS.num[,"Milk"]) #not very bad
#Milk appears almost normal 
qqnorm(Prot_DS.num[,"Fish"], main = "Fish Proteins")
qqline(Prot_DS.num[,"Fish"])
#Fish appears almost normal 

#now multivariate plot to check if the variables are multivariate normal 
names(Prot_DS.num)
x <- Prot_DS.num[, c("Red.Meat", "White.Meat", "Egg","Milk","Fish","Cereals",
                     "Starchy.Foods","Pulses.Nuts.and.Oilseeds","Fruits.and.Vegetables")]
cm <- colMeans(x)
#cm
S <- cov(x)
#S
d <- apply(x, MARGIN = 1, function(x)t(x - cm) %*% solve(S) %*% (x - cm))
#d
plot(qchisq((1:nrow(x) - 1/2) / nrow(x), df = 9), sort(d),
     xlab = expression(paste(chi[9]^2, " Quantile")),
     ylab = "Ordered distances")
abline(a = 0, b = 1)
#This graph shows that our data is a multivariate normal and can be passed for our models 
#without need of transformation.

#***************************
#Question1: PCA
##&&&&&&&&& PCA main start &&&&&&&&&&&&&&&&&&&&&&
##sparrows<-read.csv("C:/Alok/OneDrive/Rutgers_MITA/Semester2/MVA/Week4/Bumpus_sparrows.csv",stringsAsFactors = FALSE)

#Get the Correlations between the measurements
#Finding correlation
#View(Prot_DS.num)
cor.PT<-cor(Prot_DS.num)
cor.PT
#Plotting correlation
corrplot(cor.PT,method="number")
#As per above graph, most of the variables are correlated with each other.
#There seems to be Negative correlation between Eggs and Cereals
#There seems to be Negative correlation between Milk & Pulses.nuts.oilseeds
#There seems to be +ve correlation between Cereals and  Pulses.nuts.oilseeds
#There seems to be +ve correlation between White.Meat and  Pulses.nuts.oilseeds

# Using prcomp to compute the principal components (eigenvalues and eigenvectors). With scale=TRUE, variable means are set to zero, and variances set to one
Protein_pca <- prcomp(Prot_DS.num,scale=TRUE)
summary(Protein_pca)
Protein_pca
#9 Principal components are created as PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8 and PC9
head(Protein_pca) #op of this std deviation is in order PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8 and PC9
#Protein_pca
#Insights from Above PCA Output 
#Contents of Principal Components:
#PC1 is dominated by Negative effect of Cereals &Egg and Positive effect of Pulse.Nut.oilseeds 
#PC2 is dominated by Positive effect of Fish and Fruits.Veg 
#PC3 is dominated by Negative effect of Milk and Positive effect White.meat
#PC4 is dominated by Negative effect of Red.Meat and Fruits.Veg
#PC5 is dominated by Positive effect of Starchy.Foods

#From Summary of Pincipal components,
#Proportion of Variance, PC1, PC3 until PC5 explain 45%,  18%, 12%, 10% and 4% of variance respectively.
#'Cumulative Proportion' field, 90.5% of Cummulative variance is explained by PC1 until PC5
#So I will include PC1 until PC5 in my data input
#So My input variables will be reduced from 11 to 5

#Plotting Scree diagram
plot(eigen_Prot, xlab = "Component number", ylab = "Component variance", type = "l", main = "Scree diagram")
#Scree plot confirms that taking 5 Principals is enough without loosing much information.

#$x gives the new dataset #u need to rename these columns
head(Protein_pca$x)
(eigen_Prot <- Protein_pca$sdev^2) #singular values (square roots of eigenvalues) stored in sparrow_pca$sdev 
names(eigen_Prot) <- paste("PC",1:9,sep="") #Naming PC components
eigen_Prot
sumlambdas <- sum(eigen_Prot)
sumlambdas #sum of genvalues is total var of ur dataset
propvar <- eigen_Prot/sumlambdas
#Printing Proper variance per PC
propvar
#Percentage of total variance
percentvar <- (eigen_Prot/sumlambdas) *100
percentvar
#Bar plot of Percentage variance 
barplot(percentvar, main = "Bar Plot", xlab = "Principal Component", ylab = "Percentage Variance")
#As per above graph, PC1 holds 78% of ur total var, PC2 14% and so on
#Cummulative variance
cumvar_Prot <- cumsum(propvar)
cumvar_Prot
#Bar plot of Cummulative Percentage variance 
barplot(cumvar_Prot, main = "Bar Plot", xlab = "Principal Component", ylab = "Percentage Variance")

#Plotting log scree diagram
#plot(log(eigen_Prot), xlab = "Component number",ylab = "log(Component variance)", type="l",main = "Log(eigenvalue) diagram")
#Plotting Histogram
plot(Protein_pca)

#Printing our new Dataset after PCA 
#Binding with categorical columns from the original dataset
Prottyp_pca <- cbind(data.frame(Country),Protein_pca$x)
head(Prottyp_pca)
#Renaming Principal components PC1 to PC5
names(Prottyp_pca) <- c("Country", "NegCerlEgg_PostivePulse", "NegativeFish_FrutVeg", 
                           "NegateMilk_PostiveWhtMeat","NegateRedMet_FrutVeg","Postive_StarchFud",
                        "PC6","PC7","PC8","PC9")
#This is our new dataset
head(Prottyp_pca,5)

#PCA Conclusion:
#Principal Component analysis is a statistical technique that uses Orthogonal Transformation.
#It helps in reducing the number of input variables to be passed to a model.
#The principal componens are Non-correlated with each other.
#After performing PCA on Protein Consumption dataset, it can be concluded that:
#Based on per person protein consumption in Europian countries, 
#Total 9 Principal components are created.
#Contents of Principal Components:
#PC1 is dominated by Negative effect of Cereals &Egg and Positive effect of Pulse.Nut.oilseeds 
#PC2 is dominated by Positive effect of Fish and Fruits.Veg 
#PC3 is dominated by Negative effect of Milk and Positive effect White.meat
#PC4 is dominated by Negative effect of Red.Meat and Fruits.Veg
#PC5 is dominated by Positive effect of Starchy.Foods
#PC components renamed accordingly.
#From Summary of Pincipal components,
#Proportion of Variance, PC1, PC3 until PC5 explain 45%,  18%, 12%, 10% and 4% of variance respectively.
#'Cumulative Proportion' field, 90.5% of Cummulative variance is explained by PC1 until PC5
#So I will include PC1 until PC5 in my data input
#So My input variables will be reduced from 11 to 5.

##&&&&&&&&& PCA main end &&&&&&&&&&&&&&&&&&&&&&

#Question Number 3: Factor Analysis:
##&&&&&&&&& Fractal main start &&&&&&&&&&&&&&&&&&&&&&
#Creating new input dataframe for Factor analysis With ROW NAMES.
PT_Fact_1 <- read.csv("C:/Alok/OneDrive/Rutgers_MITA/Semester2/MVA/MVA_Midterm_FinalFolder/Protein_Consumption.csv",row.names=1, fill = TRUE)
head(PT_Fact_1)
names(PT_Fact_1)
#Removing Total Column
PT_Fact_1<- PT_Fact_1[,-10]
dim(PT_Fact_1)
attach(PT_Fact_1)
#Finding correlation
cor.PT<-cor(PT_Fact_1)
cor.PT
#Plotting correlation
corrplot(cor.PT,method="number")
#As per above graph, most of the variables are correlated with each other.
#There seems to be Negative correlation between Eggs and Cereals
#There seems to be Negative correlation between Milk & Pulses.nuts.oilseeds
#There seems to be +ve correlation between Cereals and  Pulses.nuts.oilseeds
#There seems to be +ve correlation between White.Meat and  Pulses.nuts.oilseeds

#To check how many factors needed, Plotting Scree diagram
#Same scree diagram used in PCA
plot(eigen_Prot, xlab = "Component number", ylab = "Component variance", type = "l", main = "Scree diagram")
#As per scree plot, there should be 5 factors, will see what parallel analysis suggests
fa.parallel(PT_Fact_1) 
#Parallel analysis Also suggests that the number of factors are 5 or 6.
vss(PT_Fact_1) # See Factor recommendations for a simple structure
nfactors(PT_Fact_1)
#nfactors suggests we can either go with 5 factors or 6 factors,
#I would prefer 5 factors as suggested in scree plot
#Factoring 
fit.PT2 <- principal(PT_Fact_1, nfactors=5, rotate="varimax") 
fit.PT2 #2 factors RC1, RC2, RC3, RC4 and RC5 are created
round(fit.PT2$values, 3)
#Above are factor values for all 9 Protein variables 
fit.PT2$loadings
# Above are the Loadings for all 9 Protein variables
for (i in c(1,2,3,4,5)) { print(fit.PT2$loadings[[1,i]])}
#Printing Communalities
fit.PT2$communality
#Rotated factor scores
head(fit.PT2$scores)
#Plotting
plot(fit.PT2)

# Plotting the relationship and mapping between variables and factors with weights
fa.diagram(fit.PT2)
#Above, output gives weigths going in RCs
#Red line indicates negative relation

#Now lets rename these factors as per their contributing variables as per above graph
colnames(fit.PT2$loadings) <- c("WhtMet_NegPulse","RedMet_Egg_Milk","Fish_NegCerl"
                                ,"FrutVeg","StrchFud")
fit.PT2$loadings

#Factor Analysis Conclusion:
#Factor analysis is a technique used to reduce number of columns.
#Factor analysis tries to find if there is any underlying latent variable in your input columns.
#After performing Factor analysis on Protein Consumption dataset, it can be concluded that:
#Based on per personProtein consumption in Europian countries,
#Total 5 factors have been formed with common variance of different Protein sources contributing to them.
#For example, Fish and Cereals contributing to RC2 positively and Negatively respectively.
#As per Above Factors, 
#For example, It did form the latent variable using the collinear variables 'RedMeat, Egg and milk;.
#The factors are Renamed accordingly.
#As per above diagram, almost all the factors have significant contribution and so 
#Though RC2 and RC3 covering lesser variance, by omitting them I would loose 4 variables 
#So, its better not to loose any of 5 factors
#So we will take All 5 Factors, RC1 to RC5 as inputs for our models
#Above factor analysis, we can conclude to reduce number of 9 Protien variables to 5 in our input dataset.

##&&&&&&&&& Fractal main end &&&&&&&&&&&&&&&&&&&&&&

#Question 2: Cluster Analysis
##&&&&&&&&& Cluster main start &&&&&&&&&&&&&&&&&&&&&&
#Cluster analysis is a technique that groups the observations into clusters based on similarities
#Clustering types: Hierarchical and nonhierarchical

#install.packages("cluster", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
#library(cluster)

#Creating new input dataframe for cluster analysis
PT_Clust_1 <- read.csv("C:/Alok/OneDrive/Rutgers_MITA/Semester2/MVA/MVA_Midterm_FinalFolder/Protein_Consumption.csv",row.names=1, fill = TRUE)
head(PT_Clust_1)
attach(PT_Clust_1)
#Removing Total Column
PT_Clust_1<- PT_Clust_1[,-10]
dim(PT_Clust_1)
# Standardizing the data with scale()
matstd.PT <- scale(PT_Clust_1)
head(matstd.PT)

######## Hierarchical Clustering ##########
# Creating a (Euclidean) distance matrix of the standardized data 
dist.PT_Clust_1 <- dist(matstd.PT, method="euclidean")
# Invoking hclust command (cluster analysis by single linkage method)      
clusPT.nn <- hclust(dist.PT_Clust_1, method = "single") 
# Plotting vertical dendrogram      
# create extra margin room in the dendrogram, on the bottom 
par(mar=c(6, 4, 4, 2) + 0.1)
plot(as.dendrogram(clusPT.nn),ylab="Distance between Countries-Single Linkage",ylim=c(0,2.5),main="Dendrogram of Protein consumption for Countries")

#Average
clusPT.avl <- hclust(dist.PT_Clust_1,method="average")
plot(clusPT.avl,hang=-1,xlab="Object",ylab="Distance",
     main="Dendrogram. Group average linkage")
#Dendogrm shows that countries have roughly 2 main group which are subdivided into smalled grps.

#Lazy option --> agnes is 1 liner command for clustering 
# We will use agnes function as it allows us to select option for data standardization, the distance measure and clustering algorithm in one single function
(agn.PT <- agnes(PT_Clust_1, metric="euclidean", stand=TRUE, method = "single"))

#  Description of cluster merging
agn.PT$merge

#Dendogram 
plot(as.dendrogram(agn.PT), xlab= "Distance between Countries",xlim=c(8,0),
     horiz = TRUE,main="Agnes Dendrogram  \n Protein Consumption for countries")

#Default - Complete Linkage
clusPT.fn <- hclust(dist.PT_Clust_1) 
plot(clusPT.fn,hang=-1,xlab="Object",ylab="Distance",
     main="Dendrogram. Farthest neighbor linkage")
#Dendogrm shows that countries have roughly 2 main group which are subdivided into smalled grps.
#If you cut at level 4 then you get around 6 Different Clusters of the countries.
#If you cut at level 5 then you get around 3 Different Clusters of the countries.

######## Non Hierarchical clustering -- K-Means Clustering##########
#Creating new input dataframe for cluster analysis
head(PT_Clust_1)
dim(PT_Clust_1)
# Standardizing the data with scale()
matstd.PT <- scale(PT_Clust_1[,1:9])
head(matstd.PT)
#matstd.PT

#Implementing K-Means Clustering with different values of k.
# K-means, k=2, 3, 4, 5, 6
# Centers (k's) are numbers thus, 10 random sets are chosen
#k=2
(kmeans2.PT <- kmeans(matstd.PT,2,nstart = 10))
# Computing the percentage of variation accounted for. Two clusters
perc.var.2 <- round(100*(1 - kmeans2.PT$betweenss/kmeans2.PT$totss),1)
names(perc.var.2) <- "Perc. 2 clus"
perc.var.2
#61% variance with k=2

# Computing the percentage of variation accounted for. Three clusters
(kmeans3.PT <- kmeans(matstd.PT,3,nstart = 10))
perc.var.3 <- round(100*(1 - kmeans3.PT$betweenss/kmeans3.PT$totss),1)
names(perc.var.3) <- "Perc. 3 clus"
perc.var.3
#49% variance with k=2

# Computing the percentage of variation accounted for. Four clusters
(kmeans4.PT <- kmeans(matstd.PT,4,nstart = 10))
perc.var.4 <- round(100*(1 - kmeans4.PT$betweenss/kmeans4.PT$totss),1)
names(perc.var.4) <- "Perc. 4 clus"
perc.var.4

# Computing the percentage of variation accounted for. Five clusters
(kmeans5.PT <- kmeans(matstd.PT,5,nstart = 10))
perc.var.5 <- round(100*(1 - kmeans5.PT$betweenss/kmeans5.PT$totss),1)
names(perc.var.5) <- "Perc. 5 clus"
perc.var.5

(kmeans6.PT <- kmeans(matstd.PT,6,nstart = 10))
# Computing the percentage of variation accounted for. Six clusters
perc.var.6 <- round(100*(1 - kmeans6.PT$betweenss/kmeans6.PT$totss),1)
names(perc.var.6) <- "Perc. 6 clus"
perc.var.6
#
(kmeans7.PT <- kmeans(matstd.PT,7,nstart = 10))
# Computing the percentage of variation accounted for. Six clusters
perc.var.7 <- round(100*(1 - kmeans7.PT$betweenss/kmeans7.PT$totss),1)
names(perc.var.7) <- "Perc. 7 clus"
perc.var.7

#It can be seen that Variance goes down as K increases...
#matstd.PT
#To Identify the Best number of K Clusters, plotting Elbow Plot
wss=c()########## empty vector to hold wss
for(i in 2:10)#### from 2 to 10 cluster
{
  km = kmeans(matstd.PT[,1:9],i)
  wss[i-1]=km$tot.withinss
}
wss
#Creating a 'elbowdt' data table with column names num and wss with the contents of wss
elbowdt = data.table(num=2:10,wss)
elbowdt
#Plotting
ggplot(elbowdt,aes(x=num,y=wss)) + geom_line()
#For k = 3 the between sum of square/total sum of square ratio tends to change slowly 
#and remain less changing as compared to others.
#Also this dataset has only 25 rows so more than 3/4 clusters would not make much sense to me.
#Therefore, k = 3 should be a good choice for the number of clusters.
#For 3 clusters, k-means = 3
# Centers (k's) are numbers thus, 10 random sets are chosen
(kmeans3.PT <- kmeans(matstd.PT,3,nstart = 10))
perc.var.3 <- round(100*(1 - kmeans3.PT$betweenss/kmeans3.PT$totss),1)
names(perc.var.3) <- "Perc. 3 clus"
perc.var.3
kmeans3.PT
kmeans3.PT$cluster
#plotting output of kmeans for 3 clusters
fviz_cluster(kmeans3.PT,data=matstd.PT)
#Clusters plotting in another way to see them more clearly
plotcluster(matstd.PT,kmeans3.PT$cluster)
#Creating separate matrices for clusters
clus1 <- matrix(names(kmeans3.PT$cluster[kmeans3.PT$cluster == 1]), 
                ncol=1, nrow=length(kmeans3.PT$cluster[kmeans3.PT$cluster == 1]))
colnames(clus1) <- "Cluster 1"
clus1

clus2 <- matrix(names(kmeans3.PT$cluster[kmeans3.PT$cluster == 2]), 
                ncol=1, nrow=length(kmeans3.PT$cluster[kmeans3.PT$cluster == 2]))
colnames(clus2) <- "Cluster 2"
clus2
clus3 <- matrix(names(kmeans3.PT$cluster[kmeans3.PT$cluster == 3]), 
                ncol=1, nrow=length(kmeans3.PT$cluster[kmeans3.PT$cluster == 3]))
colnames(clus3) <- "Cluster 3"
clus3


#clus4 <- matrix(names(kmeans4.PT$cluster[kmeans4.PT$cluster == 4]), 
#                ncol=1, nrow=length(kmeans4.PT$cluster[kmeans4.PT$cluster == 4]))
#colnames(clus4) <- "Cluster 4"
#clus4

#clus5 <- matrix(names(kmeans5.PT$cluster[kmeans5.PT$cluster == 5]), 
#                ncol=1, nrow=length(kmeans5.PT$cluster[kmeans5.PT$cluster == 5]))
#colnames(clus5) <- "Cluster 5"
#clus5

#clus6 <- matrix(names(kmeans6.PT$cluster[kmeans6.PT$cluster == 6]), 
#                ncol=1, nrow=length(kmeans6.PT$cluster[kmeans6.PT$cluster == 6]))
#colnames(clus6) <- "Cluster 6"
#clus6

#clus7 <- matrix(names(kmeans7.PT$cluster[kmeans7.PT$cluster == 7]), 
#                ncol=1, nrow=length(kmeans7.PT$cluster[kmeans7.PT$cluster == 7]))
#colnames(clus7) <- "Cluster 7"
#clus7

#Displaying all the Clolonies as in their respective clusters
list(clus1,clus2,clus3)
#Making Subsets for 3 clusters using Row filtering from the Original dataset
#(Not the scaled one)
#So below are the 3 cluster sets of Original entire dataset
#BF_Clust_1
#Using original dataframe to capture Country names as clusters have been made based on country
names
head(Prot_DS)
PT_Cl1_Dt<-subset(Prot_DS,Prot_DS$Country %in% clus1)
#PT_Cl1_Dt
PT_Cl2_Dt<-subset(Prot_DS,Prot_DS$Country %in% clus2)
#PT_Cl2_Dt
PT_Cl3_Dt<-subset(Prot_DS,Prot_DS$Country %in% clus3)
#PT_Cl3_Dt
#PT_Cl4_Dt<-subset(Prot_DS,Prot_DS$Country %in% clus4)
#PT_Cl4_Dt
#PT_Cl5_Dt<-subset(Prot_DS,Prot_DS$Country %in% clus5)
#PT_Cl5_Dt
#PT_Cl6_Dt<-subset(Prot_DS,Prot_DS$Country %in% clus6)
#PT_Cl6_Dt
#PT_Cl7_Dt<-subset(Prot_DS,Prot_DS$Country %in% clus7)
#PT_Cl7_Dt

#Printing all the columns of the Clusters formed 
#Original observations after clustering with all the variables
list(PT_Cl1_Dt,PT_Cl2_Dt,PT_Cl3_Dt)

#Trying to plot different variables to see if I can spot any groups within them.
#Plotting observations against Fruits and vegetables 
ggplot(Prot_DS,
       aes(x=Prot_DS$Country,y=Prot_DS$Fruits.and.Vegetables))+
  geom_point(size=1.8,color='dark blue')

#Plotting observations against Cereals 
ggplot(Prot_DS,
       aes(x=Prot_DS$Country,y=Prot_DS$Cereals))+
  geom_point(size=1.8,color='dark blue')

#Plotting Scatter plotss with different variables to check on what basis clusters may have formed.
#Plotting observations against Fish
ggplot(Prot_DS,
       aes(x=Prot_DS$Country,y=Prot_DS$Fish))+
  geom_point(size=1.8,color='dark blue')

#Plotting observations against Cereals 
ggplot(Prot_DS,
       aes(x=Prot_DS$Country,y=Prot_DS$Starchy.Food))+
  geom_point(size=1.8,color='dark blue')

#As per above 3 clusters formed based on similarities of per person Protien consumprion for Europian countries,
#I can roughly see that  
#Clustering might have formed based on variables 'Fruits and Veg' , Fish and upto some extent
#based on Cereals and starchy foods.
#Because There seems to be 3 different ranges of these variables in the allocated 3 clusters.
#So roughly, there are 3 groups formed manily based on the similarities of the per person Protein consumption 
#through intake of Fruits and Veg, Fish, Cereals and Starchy foods.
#Hierarchical and non hierarchical clustering has been performed on Protein consumptiondataset above.
##&&&&&&&&& Cluster main end &&&&&&&&&&&&&&&&&&&&&&


