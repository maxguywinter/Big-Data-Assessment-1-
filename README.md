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

