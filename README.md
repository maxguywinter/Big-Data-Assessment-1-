# Big-Data-Assessment-1-

## Install Packages 
library(corrplot) # for correlation matrix graph visualization.
library(ggplot2) # for graph visualizations. 
library(ggthemr)
library(ggthemes)
install.packages('devtools')
devtools::install_github('bbc/bbplot')
library(performance)
library(see)
library(car) # for VIF test.
library(caTools)
library(dplyr)
library(MASS) # for AIC test.
library(tidyverse)
library(Amelia) #this package will enable us to use the missmap function to find missing variables within the dataset.
library(mlbench)
library(flexmix)
library(caret)
library(effects)

## Install Data Set 
data <- read.csv("IntOrg_NCD_variables_2022_02_02.csv", header = TRUE) # imports the data set. The header=True command tells RStudio to use the first row of the data file as the names of each variable/column. 
attach(data) # attaches the data to your environment so that you can directly refer to the variable by name.
names(data) # shows the name of variables in the data set.
head(data)
str(data) # shows the observations and variables of the data.

# Exploratory Data Analysis

missmap(IntOrg_NCD_variables_2022_02_02, col=c("black", "purple"), legend = FALSE)
