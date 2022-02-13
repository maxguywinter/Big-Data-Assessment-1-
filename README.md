# Big-Data-Assessment-1-

### Install Packages 
library(corrplot) # for correlation matrix graph visualization.
library(ggplot2) # for graph visualizations. 
library(ggthemr)
library(ggthemes)
install.packages('devtools')
devtools::install_github('bbc/bbplot')
install.packages('naniar')
library(naniar)
library(performance)
library(see)
library(car) # for VIF test.
library(caTools)
library(dplyr)
library(MASS) # for AIC test.
library(tidyverse)
library(Amelia)
library(mlbench)
library(flexmix)
library(caret)
library(effects)

### Install Data Set 
>data <- read.csv("IntOrg_NCD_variables_2022_02_02.csv", header = TRUE) # imports the data set. The header=True command tells RStudio to use the first row of the data file as the names of each variable/column. 

>attach(data) # attaches the data to your environment so that you can directly refer to the variable by name.

>names(data) # shows the name of variables in the data set.

>head(data)

>str(data) # shows the observations and variables of the data.

### Exploratory Data Analysis

>missmap(data, col=c("blue", "red"), legend = FALSE) # Checks for missing data. 

## Coding/Data Cleaning

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

nrow(MyDataCleaned)

>MyDataCleaned[complete.cases(MyDataCleaned),]
>na.omit(MyDataCleaned)
>na.exclude(MyDataCleaned)

>MyDataCleaned %>% drop_na(Country, ISO, Sex, Year, Mean_BMI_children, Prevalence_obesity_children, Prevalence_overweight_children, Prevalence_underweight_children, Mean_BMI_adults, Prevalence_obesity_adults,Prevalence_underweight_adults,Prevalence_morbid_obesity_adults,Diabetes_prevalence, Systolic_blood_pressure, Prevalence_raised_blood_pressure, Region, Superregion, Years_of_education, Urbanisation, Western_diet_score, GDP_USD)

>nrow(MyDataCleaned)
>na.fail(MyDataCleaned)
>missmap(MyDataCleaned, col=c("blue", "red"), legend = FALSE) # Checks for missing data.
