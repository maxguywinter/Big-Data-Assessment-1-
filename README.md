################################################################################
# This code processes and analyses non-communicable disease (NCD) variables    #
# from the IntOrg_NCD_variables_2022_02_02 data set. The code was created and  #
# tested on RStudio version 4.1.1 (2021-08-10).                                #
#                                                                              #
# Code developed by Max Winter (mw636), Andre Faid (af525)                     #
# and Jessica Ndoci (jn339).                                                   #
#                                                                              #
# Firstly, the code checks the data quality and subsequently cleans and codes  # 
# the data set suitable for analysis. Next, the code performs an exploratory   #
# data analysis (EDA) to determine which analysis is appropriate for the       #
# data set. Consequently, clustering (K-means and hierarchical clustering) and #
# classification (Decision trees, Random Forest and Generalised boosted        #
# regression model) methods were utilised to analyse the data set.             # 
################################################################################

################################################################################
# We have aimed to limit the code to 80 characters per line in order for it to #
# fit comfortably on a printed page with a reasonably sized font.              #
# However, some of the code line widths are more than 80 characters so that    #
# the whole of particular sections of the code can be seen simultaneously.     #
# The code displays most clearly on our laptops and desktops in the RStudio    #
# script editor.                                                               #
################################################################################

##### SET WORKING DIRECTORY ####################################################
# User should set the relevant working directory.                              #
# Either Session --> Set Working Directory --> To Source File Location         #
# To Source File Location means that the data will be imported and saved       #
# IN THE SAME FOLDER where you saved the R Script you're working with!         #
## OR                                                                          #
## setwd("")                                                                   #
# This code used setwd("~/Desktop/Big data assessment 1") as the working       #
# directory that contained the IntOrg_NCD_variables_2022_02_02 data set.       # 
################################################################################

##### RELEVANT PACKAGES ########################################################
##### Users should install the relevant packages below                         #
#                                                                              #
# Use install.packages('') if unable install the required packages from        #
# library ()                                                                   #
################################################################################
library(corrplot) # for correlation matrix graph visualization.
library(ggplot2) # for graph visualizations. 
library(ggthemes) # for extra themes, geoms and scales from package "ggplot2"
library(tidyverse) # set of packages for data representation and design
library(Amelia) # for incomplete data
library(naniar) # for plotting missing data and imputations
library(gridExtra) # for drawing up tables and plots
library(devtools) # for development tools - functions
library(dplyr) # for handling large data
library(effects) # for displaying graphs and tables
library(janitor) # for cleaning data
library(ggpubr) # for customizing graph visualizations from package "ggplot2"

library(cluster) # for cluster analysis
library(factoextra) # to help extract and visualize our analyses
library(dendextend) # for creating visualizzing and comparing hierarchical clustering
library(gplots) # for plotting data

library(mlbench)
library(caTools) # for test and train split
library(caret) # for plotting classification and regression models
library(dismo) # for prediction of environmental similarity
library(gbm) # for boosted regression model in classification
library(rpart) # for building classificationa nd regression trees
library(rpart.plot) # to plot "rpart" function
library(randomForest) # for classification and regression of trees

##### DATA INPUT ###############################################################
##### The data are age- & sex-stratified summary statistics for each study.    #
##### The data file must be a csv file with the following columns:             #
##### id_study: alphanumeric string (we suggest using only A-Z, 0-9 and _).    #
##### We suggest "iso3 code"_"data collection mid-year"_"survey name"          #
##### DATA CHARACTERISTICS                                                     #
##### sex: must be either "male" or "female", as the two sexes are modelled    #
##### separately                                                               #
##### age: positive number >=18 for adult-only analyses; 'age' is the mid-age  #
##### of the age group for this datapoint                                      #
##### mean_bmi: positive number                                                #
##### se_bmi: positive number                                                  #
##### STUDY CHARACTERISTICS                                                    #
##### mid_year: integer year, as required by the non-linear component of the   #
##### model                                                                    #
##### survey_type: "National", "Subnational" or "Community"; used in the       #
##### analysis to account for possible bias in non-national studies            #
##### urban_rural: "urban", "rural" or "both"; used in the analysis to account #
##### for possible bias in urban-only or rural-only studies                    #
##### Country: country names, which must match covariate file                  #
##### Region: region names, which must match covariate file                    #
##### Superregion: superregion names, which must match covariate file          #
################################################################################
data <- read.csv("IntOrg_NCD_variables_2022_02_02.csv", header = TRUE) # imports the data set. The header=True command tells RStudio to use the first row of the data file as the names of each variable/column. 
attach(data) # attaches the data to your environment so that you can directly refer to the variable by name.
names(data) # shows the name of variables in the data set.
head(data)
summary(data)
str(data) # shows the observations and variables of the data.

##### Exploratory Data Analysis ################################################
##### The data are age- & sex-stratified summary statistics for each study.    #
##### The data file must be a csv file with the following columns:             #
##### id_study: alphanumeric string (we suggest using only A-Z, 0-9 and _).    #
##### We suggest "iso3 code"_"data collection mid-year"_"survey name"          #
##### DATA CHARACTERISTICS                                                     #
##### sex: must be either "male" or "female", as the two sexes are modelled    #
##### separately                                                               #
##### age: positive number >=18 for adult-only analyses; 'age' is the mid-age  #
##### of the age group for this datapoint                                      #
##### mean_bmi: positive number                                                #
##### se_bmi: positive number                                                  #
##### STUDY CHARACTERISTICS                                                    #
##### mid_year: integer year, as required by the non-linear component of the   #
##### model                                                                    #
##### survey_type: "National", "Subnational" or "Community"; used in the       #
##### analysis to account for possible bias in non-national studies            #
##### urban_rural: "urban", "rural" or "both"; used in the analysis to account #
##### for possible bias in urban-only or rural-only studies                    #
##### Country: country names, which must match covariate file                  #
##### Region: region names, which must match covariate file                    #
##### Superregion: superregion names, which must match covariate file          #
################################################################################

##### Data Coding/Cleaning #####################################################
data[data == "" | data == " "] <- NA 
data[data == "N/A" | data == "n/a" | data == "N/a"| data == "-"] <- NA
data$Country <- factor(data$Country) # Changing categorical variables to factors 
data$ISO <- factor(data$ISO)
data$Sex <- factor(data$Sex)
data$Region <- factor(data$Region)
data$Superregion<- factor(data$Superregion)
missmap(data, col=c("red", "#1380A1"), x.cex = 0.3,legend = TRUE) # Checks for missing data. 
pct_miss(data) # Percent of ALL data frame values that are missing
pct_miss_case(data) # Percent of rows with any value missing
pct_complete_case(data) # Percent of rows that are complete (no values missing) 
get_dupes(data) # Check for duplicates and wrong input 
unique(data$Country)
unique(data$ISO)
unique(data$Sex)
unique(data$Region)
unique(data$Superregion)

##### Cleaned Data #############################################################
MyDataCleaned<-na.omit(data) # removes NA's
nrow(data) # shows obs
nrow(MyDataCleaned) # shows obs of new data set
na.fail(data) # checks for missing values 
na.fail(MyDataCleaned)
MyDataCleaned<-tibble::rowid_to_column(MyDataCleaned, "Identification No.") # generate ID variable
MyDataCleaned<-dplyr::rename(MyDataCleaned, ID = "Identification No.") # renaming variables
missmap(MyDataCleaned,col=c("red", "#1380A1"), x.cex = 0.3, legend = TRUE) # Checks for missing data.
pct_miss(MyDataCleaned) # Percent of ALL data frame values that are missing
pct_miss_case(MyDataCleaned) # Percent of rows with any value missing
pct_complete_case(MyDataCleaned)
names(MyDataCleaned) 
summary(MyDataCleaned)
str(MyDataCleaned)

##### Graphical visualization ##################################################
#                                                                              #
#
#                                                                            
#
#
#
#
#
#
#
#                
#                                                                              #
################################################################################

##### Correlation Analysis #####################################################
datacorrelation <- subset(MyDataCleaned, select =-c(1,2,3,4,5,17,18)) # subset data to remove factors. 
nums <- unlist(lapply(datacorrelation, is.numeric)) # used to create correlation plot from numeric values
Datacorr <- cor(datacorrelation[, nums])
corrplot(Datacorr, type = "upper", method = "number", number.cex = 0.5, tl.cex = 0.5, title = "\n\n Correlation Plot \n") # correlation matrix graph visualization. 
pairs(datacorrelation) # graphical matrix visualization.

##### lollipop chart ###########################################################
Fig.5 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$ISO, y = MyDataCleaned$Year)) + # Countries (ISO) and Years
  geom_point(size = 2, color = "#1380A1" ) +
  ggtitle("Country (ISO)", subtitle = "Country data aviable for each year") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "Countries (ISO)", y = "Years", tag = "Fig.5") +
  theme_economist()
Fig.5 

##### Bar graphs ###############################################################
Fig.6 <- ggplot(MyDataCleaned, aes(MyDataCleaned$Sex)) + # sex 
  geom_bar(fill = "#1380A1") +
  labs(title = "Sex", subtitle = "Number of Males and Females", x= "Sex", tag = "Fig.6") +
  theme_economist()
Fig.6

Fig.7 <- ggplot(MyDataCleaned, aes(MyDataCleaned$Region)) + # Region
  geom_bar(fill = "#1380A1") +
  labs(title = "Region", subtitle = "Number of Regions", x= "Sex", tag = "Fig.7") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "Countries (ISO)", y = "Years", tag = "Fig.5") +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))
Fig.7

Fig.8 <- ggplot(MyDataCleaned, aes(MyDataCleaned$Superregion)) + # Super Region  
  geom_bar(fill = "#1380A1") +
  labs(title = "Super Region", subtitle = "Number of Super Regions", x= "Sex", tag = "Fig.8") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "Countries (ISO)", y = "Years", tag = "Fig.5") +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))
Fig.8

##### Box plots ################################################################
Fig.9 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_obesity_children)) + # Prevalence obesity children
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Child Obesity", y = "Child Obesity Prevalence") +
  theme_economist()
Fig.10 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_overweight_children)) + # Prevalence overweight children
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Child Overweight", y = "Child Overweight Prevalence") +
  theme_economist()
Fig.11 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_underweight_children)) + # Prevalence underweight children
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Child Underweight", y = "Child Underweight Prevalence") +
  theme_economist()
Fig.12 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Mean_BMI_children)) + # Mean BMI children
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Children Mean BMI", y = "Mean BMI") +
  theme_economist()

ggarrange(Fig.9, Fig.10, Fig.11, Fig.12 + rremove("x.text"), 
          labels = c("Fig.7", "Fig.8", "Fig.9", "Fig.10"),
          ncol = 2, nrow = 2)

Fig.13 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_morbid_obesity_adults)) + # Prevalence morbid_obesity adults
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Morbid Obesity", y = "Adult Morbid Obesity Prevalence") + 
  theme_economist()
Fig.14 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_underweight_adults)) + # Prevalence underweight adults
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Underweight", y = "Adult Underweight Prevalence") +
  theme_economist()
Fig.15 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_obesity_adults)) + # Prevalence obesity adults
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Obesity", y = "Adult Obesity Prevalence") +
  theme_economist()
Fig.16 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Mean_BMI_adults)) + # Mean BMI adults
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Mean BMI", y = "Mean BMI") +
  theme_economist()

ggarrange(Fig.13, Fig.14, Fig.15, Fig.16 + rremove("x.text"), 
          labels = c("Fig.11", "Fig.12", "Fig.13", "Fig.14"),
          ncol = 2, nrow = 2)

Fig.17 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Urbanisation)) + # urbanization 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Urbanisation", y = "Urbanisation Score (0-1)") +
  theme_economist()
Fig.18 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Western_diet_score)) + # Western diet score 
  geom_boxplot() + 
  labs(title = "Western Diet Score", y = "Western Diet Score (-2.5-4.5)") +
  theme_economist()

ggarrange(Fig.15, Fig.16 + rremove("x.text"), 
          labels = c("Fig.17", "Fig.18"),
          ncol = 2)

##### Histograms ###############################################################
Fig.19 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Diabetes_prevalence)) + # Diabetes Prevalence
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Diabetes Prevalence", x = "Diabtes Prevalence") +
  theme_economist()
Fig.19
Fig.20 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Systolic_blood_pressure)) + # Systolic blood Pressure
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Systolic blood Pressure", x = "Systolic blood Pressure") +
  theme_economist()
Fig.20
Fig.21 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_raised_blood_pressure)) + # Prevalence Raised Blood Pressure
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Prevalence Raised Blood Pressure", x = "Prevalence Raised Blood Pressure") +
  theme_economist()
Fig.21
Fig.22 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Years_of_education)) + # Years of Education
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Years of Education", x = "Years of Education") +
  theme_economist()
Fig.22
Fig.23 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$GDP_USD)) + # GDP
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "GDP($)", x = "GDP($)") +
  theme_economist()
Fig.23

##### Variable Comparisons ######################################################

##### Box plots ################################################################
Fig.24 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Years_of_education, x=  MyDataCleaned$Superregion)) + # Super Region and Years of education 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Years of Education for each Super Region", y = "Years of Education", x = "Super Region", tag = "Fig.24") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))

Fig.25<- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Diabetes_prevalence, x=  MyDataCleaned$Region)) + # Region and Diabetes prevalence 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Prevalence of Diabetes in each Region", y = "Disabetes Prevalence", x = "Region", tag = "Fig.25") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))

Fig.26 <-ggplot(MyDataCleaned, aes(y = MyDataCleaned$Urbanisation, x=  MyDataCleaned$Region)) + # Region and Urbanisation 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Urbanisation of each Region", y = "Urbanisation Score (0-1)", x = "Region", tag = "Fig.26") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))

Fig.27 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Western_diet_score, x=  MyDataCleaned$Superregion)) + # Super region and Western diet score
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Westen Diet Score each Super Region", y = "Western Diet Score (-2.5-4.5)", x = "Super Region", tag = "Fig.27") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))

ggarrange(Fig.24, Fig.25, Fig.26,Fig.27 + rremove("x.text"), 
          ncol = 2, nrow = 2)

Fig.28 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Diabetes_prevalence, x=  MyDataCleaned$Sex)) + # Sex and Diabetes (box plot)
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Prevalence of Diabetes for each Sex", y = "Disabetes Prevalence", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()

Fig.29 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Years_of_education, x=  MyDataCleaned$Sex)) + # Sex and education (box plot)
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Years of Education for each Sex", y = "Years of Education", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()

Fig.30 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_morbid_obesity_adults, x=  MyDataCleaned$Sex)) + # Sex and obesity (box plot)
  geom_boxplot() + 
  labs(title = "Prevalence of Adult Morbid Obesity for each Sex", y = "Morbid Obesity Prevalence", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()

Fig.31 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_underweight_adults, x=  MyDataCleaned$Sex)) + # Sex and obesity (box plot)
  geom_boxplot() + 
  labs(title = "Prevalence of Underweight Adults for each Sex", y = "Prevalence of Underweight Adults", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()

ggarrange(Fig.28, Fig.29, Fig.30, Fig.31 + rremove("x.text"), 
          labels = c("Fig.28", "Fig.29", "Fig.30", "Fig.31"),
          ncol = 2, nrow = 2)

##### Scatter plots ############################################################
Fig.32 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Western_diet_score, y = MyDataCleaned$Prevalence_obesity_children)) + # West diet score and Prevalence obesity children
  geom_point(colour = "#1380A1") +
  labs(title = "Western Diest and Prevalence of Child Obesity", y = "Prevalence of Child Obesity", x = "Western Diet Score (-2.5-4.5)", tag = "Fig.32") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.33 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Western_diet_score, y = MyDataCleaned$Prevalence_obesity_adults)) + # West diet score and Prevalence obesity adults
  geom_point(colour = "#1380A1") +
  labs(title = "Western Diest and Prevalence of Adult Obesity", y = "Disabetes Prevalence", x = "Western Diet Score (-2.5-4.5)", tag = "Fig.33") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

ggarrange(Fig.32, Fig.33+ rremove("x.text"), 
          ncol = 2)

Fig.34 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Years_of_education, y = MyDataCleaned$Prevalence_obesity_children)) + # education and Prevalence obesity children
  geom_point(colour = "#1380A1") +
  labs(title = "Years of Education and Prevalence of Child Obesity", y = "Prevalence of Child Obesity", x = "Years of Education", tag = "Fig.34") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.35 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Years_of_education, y = MyDataCleaned$Prevalence_obesity_adults)) + # education and Prevalence obesity adults
  geom_point(colour = "#1380A1") +
  labs(title = "Western Diest and Prevalence of Adult Obesity", y = "Prevalence of Adult Obesity", x = "Years of Education", tag = "Fig.35") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

ggarrange(Fig.34, Fig.35 + rremove("x.text"), 
          ncol = 2)

Fig.36 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_underweight_children , y = MyDataCleaned$GDP_USD)) + # GDP and Prevalence underweight children
  geom_point(colour = "#1380A1") +
  labs(title = "GDP($) and Prevalence of underweight Children", y = "GDP($)", x = "Prevalence of underweight Children", tag = "Fig.36") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.37 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_underweight_adults, y = MyDataCleaned$GDP_USD)) + # GDP and Prevalence underweight adults
  geom_point(colour = "#1380A1") +
  labs(title = "GDP($) and Prevalence of underweight Adults", y = "GDP($)", x = "Prevalence of underweight Adults", tag = "Fig.37") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

ggarrange(Fig.36, Fig.37+ rremove("x.text"), 
          ncol = 2)

Fig.38 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Urbanisation, y = MyDataCleaned$Prevalence_overweight_children)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Urbanisation and Prevalence of Overweight Children", y = "Prevalence of Oveweight Children", x = "Urbanisation Score (0-1)", tag = "Fig.38") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()
Fig.38 # Urbanization and Prevalence of Overweight Children

Fig.39 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Mean_BMI_adults, y = MyDataCleaned$Systolic_blood_pressure)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Sysolic Blood Pressure and Mean BMI (Adults)", y = "Systolic Blood Pressure", x = "Mean BMI", tag = "Fig.39") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()
Fig.39 # Systolic blood pressure and Mean BMI Adults

Fig.40 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Mean_BMI_adults, y = MyDataCleaned$Prevalence_obesity_adults)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Prevalence of Obesity and Mean BMI (Adults)", y = "Prevalence of Obesity", x = "Mean BMI", tag = "Fig.40") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()
Fig.40 # Adult Mean BMI and Prevalence obesity adults

Fig.41 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_overweight_children, y = MyDataCleaned$Prevalence_obesity_adults)) + # Prevalence overweight children and Prevalence obesity adults
  geom_point(colour = "#1380A1") +
  labs(title = "Prevalence of Adult Obesity and Overweight Children", y = "Prevalence of Adult Obesity", x = "Prevalence of Overweight Children", tag = "Fig.41") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.42 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_underweight_children, y = MyDataCleaned$Prevalence_underweight_adults)) + # Prevalence underweight children and Prevalence underweight adults
  geom_point(colour = "#1380A1") +
  labs(title = "Prevalence of Underweight Adults and Children", y = "Prevalence Underweight Adults", x = "Prevalence Underweight Children", tag = "Fig.42") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

ggarrange(Fig.41, Fig.42+ rremove("x.text"), 
          ncol = 2)

##### Clustering ################################################################


##### Data Coding/Cleaning for clustering ######################################
df <- subset(MyDataCleaned, select =-c(1,2,3,4,5,17,18))
df <- na.omit(df)
df2 <- scale(df)
head(df2)
str(df2)

##### Euclidean distances between the points and visualizations ################
distance <- get_dist(df2)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

##### K-means ##################################################################
k2 <- kmeans(df2, centers = 2, nstart = 25)
k3 <- kmeans(df2, centers = 3, nstart = 25)
k4 <- kmeans(df2, centers = 4, nstart = 25)
k5 <- kmeans(df2, centers = 5, nstart = 25)

p1 <- fviz_cluster(k2, geom = "point", data = df2) + ggtitle("k = 2") # plots to compare
p2 <- fviz_cluster(k3, geom = "point",  data = df2) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df2) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 5")
grid.arrange(p1, p2, p3, p4, nrow = 2)

##### Elbow Method #############################################################
set.seed(123)
fviz_nbclust(df2, kmeans, method = "wss") # 4? clusters k

##### Silhoutte Method #########################################################
fviz_nbclust(df2, kmeans, method = "silhouette") # 3 clusters k 

##### Gap Statistic ############################################################
gap_stat <- clusGap(df2,FUN = kmeans, nstart = 25,K.max = 10, B = 50)
fviz_gap_stat(gap_stat) # 10? clusters k
set.seed(123)
gap_stat <- clusGap(df2, FUN = kmeans, nstart = 25, # compute gap statistic
                    K.max = 10, B = 50)

##### Interpretation of the clusters ###########################################
final <- kmeans(df2, 4, nstart = 25)
fviz_cluster(final, data = df2)

df %>% 
  mutate(Cluster = final$cluster) %>% ##########################################
  group_by(Cluster) %>% 
  summarise_all("mean") 

##### Hierarchical clustering ##################################################

##### arcsin transformation #################################################### 
arcsin_transformation <- function(x) asin(x/100)

dend <- df2 %>% arcsin_transformation %>%
  dist %>% hclust(method = "com") %>% 
  as.dendrogram %>%
  set("branches_k_color", k = 3) %>% 
  ladderize

dend2 <- dend %>% 
  set("labels", c(1:14000)) # change labels

dend2 %>% labels

##### Heat Map #################################################################
gplots::heatmap.2(as.matrix(df2), 
                  main = "Heat Map",srtCol = 60, dendrogram = "row",
                  Rowv = dend2,Colv = "NA", trace="none", margins = c(8,8), cexCol = 0.5,     
                  key.xlab = "Health Factors", denscol = "grey", density.info = "density",col = colorspace::diverge_hcl)

##### Hierarchical clustering methods ##########################################
hclust_methods <- c("single", "complete", "average", "median")
df2_dendlist <- dendlist()

for(i in seq_along(hclust_methods)) {
  tmp_dend <- df2 %>% arcsin_transformation %>% dist %>%hclust(method=hclust_methods[i]) %>% as.dendrogram 
  df2_dendlist<-dendlist(df2_dendlist,tmp_dend)
}
names(df2_dendlist) <- hclust_methods

corrplot::corrplot(cor.dendlist(df2_dendlist),"pie","lower")

##### compare trees produced by different methods using a tanglegram ###########
dend1 <- df2 %>% arcsin_transformation %>% dist %>% hclust(method = "complete") %>% 
  as.dendrogram %>% color_branches(k=3) %>% ladderize

dend2 <- df2 %>% arcsin_transformation %>% dist %>% hclust(method = "average") %>% 
  as.dendrogram %>% color_branches(k=3) %>% ladderize

dends<-dendlist(tree1 = dend1, tree2 = dend2) 
dends %>% tanglegram(margin_inner = 7)

##### Classification ###########################################################

##### Decision trees ###########################################################

##### Decision tree Sex ########################################################
tree_default <- rpart(MyDataCleaned$Sex ~ ., data = MyDataCleaned)
rpart.plot(tree_default,extra=2, under = TRUE, varlen=0, faclen=0)

tree_full <- rpart(MyDataCleaned$Sex ~., data= MyDataCleaned, control=rpart.control(minsplit=2, cp=0))
rpart.plot(tree_full,extra=2,under=TRUE,varlen=0,faclen=0,cex=.7)

confusion_table<-table(MyDataCleaned$Sex, 	predict(tree_default,MyDataCleaned,type="class"))
confusion_table

correct <- sum(diag(confusion_table))
error <- sum(confusion_table)-correct
accuracy <- correct / (correct+error);accuracy
confusionMatrix(data= predict(tree_default,MyDataCleaned,type="class"), reference = MyDataCleaned$Sex)

##### Decision tree Super Region ###############################################
tree_default2 <- rpart(MyDataCleaned$Superregion ~ ., data = MyDataCleaned)
rpart.plot(tree_default2,extra=2, under = TRUE, varlen=0, faclen=0)

tree_full2 <- rpart(MyDataCleaned$Superregion ~., data= MyDataCleaned, control=rpart.control(minsplit=2, cp=0))
rpart.plot(tree_ful2l,extra=2,under=TRUE,varlen=0,faclen=0,cex=.7)

confusion_table2<-table(MyDataCleaned$Superregion, 	predict(tree_default2,MyDataCleaned,type="class"))
confusion_table2

correct2 <- sum(diag(confusion_table2))
error2 <- sum(confusion_table2)-correct2
accuracy2 <- correct2 / (correct2+error2);accuracy2

##### Random Forest ############################################################

##### Data Coding/Cleaning for classification ##################################
df3 <- MyDataCleaned[,c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,18)] 

##### Data Train and Test sets #################################################
set.seed(123) # Makes simulations random numbers the same to ensure all results, figures  are reproducible.

split<-sample.split(df3, SplitRatio = 0.7)  
training_set<-subset(df3,split==TRUE)
test_set<-subset(df3,split==FALSE) 
dim(training_set);dim(test_set)
topredict_set<-test_set[1:18]  
topredict_set2<-test_set[2:19]  

##### Random Forest Super Region ###############################################
model_rf<-randomForest(training_set$Superregion~.,data=training_set,importance=TRUE, ntree=1000) 
preds_rf <- predict(model_rf, topredict_set)              
(conf_matrix_forest <- table(preds_rf, test_set$Superregion))
confusionMatrix(conf_matrix_forest) 

##### Random Forest Sex ########################################################
model_rf2<-randomForest(training_set$Sex~.,data=training_set,importance=TRUE, ntree=1000) 
preds_rf2 <- predict(model_rf2, topredict_set2)              
(conf_matrix_forest2 <- table(preds_rf2, test_set$Sex))
confusionMatrix(conf_matrix_forest2) 

##### Generalised boosted regression model #####################################

##### Data Coding/Cleaning for gbm #############################################
df4 <- MyDataCleaned[,c(4,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,18)] 
df4$Sex<-ifelse(df4$Sex=="Female",1,0)

##### Data Train and Test sets #################################################
set.seed(123) # Makes simulations random numbers the same to ensure all results, figures  are reproducible.

split2<-sample.split(df4, SplitRatio = 0.7)  
training_set2<-subset(df4,split==TRUE)
test_set2<-subset(df4,split==FALSE) 
dim(training_set2);dim(test_set2)
 
##### Sex gbm ################################################################## (change learning rate in names to 05)
sex.tc5.lr01 <- gbm.step(data=training_set2, gbm.x = 2:21, 
                            gbm.y = 1, family = "bernoulli", tree.complexity = 5, 	learning.rate = 0.005, bag.fraction = 0.5) 

sex.simp <- gbm.simplify(sex.tc5.lr01, n.drops = 5)
sex.simp$pred.list[[4]]

sex.tc5.lr005.simp <- gbm.step(training_set2, 	gbm.x=sex.simp$pred.list[[4]], gbm.y=1, 
                                  tree.complexity=5, learning.rate=0.005)


preds <- predict.gbm(sex.tc5.lr005.simp, test_set2, 	n.trees=sex.tc5.lr005.simp$gbm.call$best.trees, type="response")
pred.limit<-0.25
confusionMatrix(table(as.numeric(preds>pred.limit),
                      test_set2[,1]),positive="1")




























