# Big-Data-Assessment-1-

### Install Packages 
>library(corrplot) # for correlation matrix graph visualization.

>library(ggplot2) # for graph visualizations. 

>library(ggthemr)

>library(ggthemes)

>install.packages('devtools')

>devtools::install_github('bbc/bbplot')

>install.packages('naniar')

>library(naniar)

>library(performance)

>library(see)

>library(car) # for VIF test.

>library(caTools)

>library(dplyr)

>library(MASS) # for AIC test.

>library(tidyverse)

>library(Amelia)

>library(mlbench)

>library(flexmix)

>library(caret)

>library(effects)

>install.packages("janitor")

>library(janitor)


### Install Data Set 
>data <- read.csv("IntOrg_NCD_variables_2022_02_02.csv", header = TRUE) # imports the data set. The header=True command tells RStudio to use the first row of the data file as the names of each variable/column. 
>attach(data) # attaches the data to your environment so that you can directly refer to the variable by name.
>names(data) # shows the name of variables in the data set.
>head(data)
>summary(data)
>str(data) # shows the observations and variables of the data.

### Exploratory Data Analysis
>missmap(data, col=c("blue", "red"), legend = FALSE) # Checks for missing data. 

### Coding/Data Cleaning
#generate ID variable
>MyDataCleaned<-tibble::rowid_to_column(data, "Identification No.")
#renaming variables
>MyDataCleaned<-dplyr::rename(MyDataCleaned, ID = "Identification No.")
#re-coding
>MyDataCleaned[MyDataCleaned == "" | MyDataCleaned == " "] <- NA
>MyDataCleaned[MyDataCleaned == "N/A" | MyDataCleaned == "n/a" | MyDataCleaned == "N/a"| MyDataCleaned == "-"] <- NA

>missmap(MyDataCleaned, col=c("blue", "red"), legend = FALSE) # Checks for missing data. 
>pct_miss(MyDataCleaned) # Percent of ALL data frame values that are missing
>pct_miss_case(MyDataCleaned) # Percent of rows with any value missing
>pct_complete_case(MyDataCleaned) # Percent of rows that are complete (no values missing) 

>MyDataCleaned2 <-na.omit(MyDataCleaned)
>nrow(MyDataCleaned)
>nrow(MyDataCleaned2)
>na.fail(MyDataCleaned)
>na.fail(MyDataCleaned2)
#drop variables
>MyDataCleaned2<-dplyr::select(MyDataCleaned2, -ID)
#generate ID variable
>MyDataCleaned3<-tibble::rowid_to_column(MyDataCleaned2, "Identification No.")
#renaming variables
>MyDataCleaned3<-dplyr::rename(MyDataCleaned3, ID = "Identification No.")

>missmap(MyDataCleaned3, col=c("blue", "red"), legend = FALSE) # Checks for missing data.
>pct_miss(MyDataCleaned3) # Percent of ALL data frame values that are missing
>pct_miss_case(MyDataCleaned3) # Percent of rows with any value missing
>pct_complete_case(MyDataCleaned3)
>names(MyDataCleaned3) 
>head(MyDataCleaned3)
>str(MyDataCleaned3)
>summary(MyDataCleaned3)

#Correlation Analysis (correlation matrix and scatter plot matrix)
>datacorrelation <- subset(MyDataCleaned3, select =-c(1,2,3,4,5,17,18)) # subset data to remove x,y, month and day as they are factors. 
>nums <- unlist(lapply(datacorrelation, is.numeric)) # used to create correlation plot from numeric values in fire data.
>Datacorr <- cor(datacorrelation[, nums])
>corrplot(Datacorr, type = "upper", method = "number")# correlation matrix graph visualization. 
>pairs(datacorrelation) # graphical matrix visualization. 

>attach(MyDataCleaned3)
>get_dupes(MyDataCleaned3, Country, ISO, Sex, Year, Mean_BMI_children, Prevalence_obesity_children, Prevalence_overweight_children, Prevalence_underweight_children, Mean_BMI_adults, Prevalence_obesity_adults, Prevalence_underweight_adults, Prevalence_morbid_obesity_adults, Diabetes_prevalence, Systolic_blood_pressure, Prevalence_raised_blood_pressure, Region, Superregion, Years_of_education, Urbanisation, Western_diet_score, GDP_USD)

>typeof(MyDataCleaned3)
>names(MyDataCleaned3)
>unique(Country)
>unique(ISO)
>unique(Sex)
>unique(Region)
>unique(Superregion)

>Country <- factor(Country)
>ISO <- factor(ISO)
>Sex <- factor(Sex)
>Region <- factor(Region)
>Superregion <- (Superregion)




#histogram of education on its own 
>h <- hist(Years_of_education,ylim=c(0,1500), main = "histogram showing years of education", breaks = 20, xlab = "years of education", col = "Red")

>text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))

#histogram of diabetes prevalence on its own
>g <- hist(Diabetes_prevalence,ylim=c(0,5000), main = "histogram showing prevalance of diabetes", breaks = 20, xlab = "prevalance of diabetes", col = "Red")

>text(g$mids,g$counts,labels=g$counts, adj=c(0.5, -0.5))

#histogram of systolic blood pressure on its own
>d<- hist(Systolic_blood_pressure, ylim=c(0,3000), main = "histogram showing systolic blood pressure", breaks = 15, xlab = "systolic blood pressure", col = "Red")

>text(d$mids,d$counts,labels=d$counts, adj=c(0.5, -0.5))

#histogram of GDP_USD on its own
>hist(GDP_USD,ylim=c(0,3000), main = "histogram showing GDP per capita in USD", breaks = 100, xlab = "GDP_USD", col = "Red")

#histogram showing GDP per capita against super region
>ggplot(df, aes(x=GDP_USD, color=Superregion)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")




>df <- MyDataCleaned3
>df <- MyDataCleaned3$GDP_USD
>df <- na.omit(df)
>df <- scale(df)
>head(df)


>plot(df)
>distance <- get_dist(df)
>fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))



>k2 <- kmeans(df, centers = 2, nstart = 25)
>k3 <- kmeans(df, centers = 3, nstart = 25)
>k4 <- kmeans(df, centers = 4, nstart = 25)
>k5 <- kmeans(df, centers = 5, nstart = 25)
>fviz_cluster(k2, data = df)
>p1 <- fviz_cluster(k2, geom = "point", data = df) + ggtitle("k = 2")
>p2 <- fviz_cluster(k3, geom = "point",  data = df) + ggtitle("k = 3")
>p3 <- fviz_cluster(k4, geom = "point",  data = df) + ggtitle("k = 4")
>p4 <- fviz_cluster(k5, geom = "point",  data = df) + ggtitle("k = 5")

