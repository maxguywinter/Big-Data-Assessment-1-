Brief: You work for a statistical consultancy that has been commissioned to carry out analysis for a major global organisation. This organisation has collected a large dataset but doesn’t have the expertise to analyse it. Given that they are working with big data, they are also not sure of the quality of the data. They have asked you to write R code for processing and analysing the data. The code needs to have extensive comments, so that when your client runs it, they will understand what each part does, and what types of analytical techniques you have used. You will also need to include comments to tell them which version of R to use, and the versions of any packages required for the analysis.  You should start by checking data quality, then you should carry out EDA (exploratory data analysis) to decide which specific types of analysis are appropriate. Finally, you should carry out the analysis, including appropriate plots. Both clustering and classification methods are generally appropriate for this dataset. You will need to decide how to split the work within your team. The client is paying you for about 20 to 30 hours work if you are in a group of 4. (Grade 78%)



```{r} 

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
library(caret) # for plotting classification models
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
datacorrelation <- subset(MyDataCleaned, select =-c(1,2,3,4,5,17,18)) # This code subsets the data in order to remove factors. 
nums <- unlist(lapply(datacorrelation, is.numeric)) # This code creates a correlation plot from numeric values
Datacorr <- cor(datacorrelation[, nums])
corrplot(Datacorr, type = "upper", method = "number", number.cex = 0.5, tl.cex = 0.5, title = "\n\n Correlation Plot \n") # This code creates a correlation matrix graph visualization in order to help us visualize the correlation between all variables. 
pairs(datacorrelation) # This code displays a graphical matrix visualization.

##### lollipop chart ###########################################################
Fig.5 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$ISO, y = MyDataCleaned$Year)) + 
  geom_point(size = 2, color = "#1380A1" ) +
  ggtitle("Country (ISO)", subtitle = "Country data aviable for each year") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "Countries (ISO)", y = "Years", tag = "Fig.5") +
  theme_economist()
Fig.5 # Figure 5 compares the variables countries (ISO) and years

##### Bar graphs ###############################################################
Fig.6 <- ggplot(MyDataCleaned, aes(MyDataCleaned$Sex)) + 
  geom_bar(fill = "#1380A1") +
  labs(title = "Sex", subtitle = "Number of Males and Females", x= "Sex", tag = "Fig.6") +
  theme_economist()
Fig.6 # Figure 6 displays a bar graph showing the number of males and females within the variable "sex"

Fig.7 <- ggplot(MyDataCleaned, aes(MyDataCleaned$Region)) + 
  geom_bar(fill = "#1380A1") +
  labs(title = "Region", subtitle = "Number of Regions", x= "Sex", tag = "Fig.7") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "Countries (ISO)", y = "Years", tag = "Fig.5") +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))
Fig.7 # Figure 7 is a bar graph displaying the number of regions within the variable "region"

Fig.8 <- ggplot(MyDataCleaned, aes(MyDataCleaned$Superregion)) +
  geom_bar(fill = "#1380A1") +
  labs(title = "Super Region", subtitle = "Number of Super Regions", x= "Sex", tag = "Fig.8") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "Countries (ISO)", y = "Years", tag = "Fig.5") +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))
Fig.8  # Figure 8 is a bar graph displaying the number of super regions within the variable "super region"

##### Box plots ################################################################
Fig.9 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_obesity_children)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Child Obesity", y = "Child Obesity Prevalence") +
  theme_economist() 
Fig.9 # Figure 9 displays a box plot showing the prevalence of obese children

Fig.10 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_overweight_children)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Child Overweight", y = "Child Overweight Prevalence") +
  theme_economist()
Fig.10 # Figure 10 is a box plot that shows the prevalence of overweight children
  
Fig.11 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_underweight_children)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Child Underweight", y = "Child Underweight Prevalence") +
  theme_economist()
Fig.11 # Figure 11 is a box plot that displays the prevelance of underweight children
  
Fig.12 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Mean_BMI_children)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Children Mean BMI", y = "Mean BMI") +
  theme_economist()
Fig.12 # Figure 12 shows a box plot illustrating the mean BMI within children

ggarrange(Fig.9, Fig.10, Fig.11, Fig.12 + rremove("x.text"), 
          labels = c("Fig.7", "Fig.8", "Fig.9", "Fig.10"),
          ncol = 2, nrow = 2) # This code arranges the above figures into one plot in order to see all figures better and easier

Fig.13 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_morbid_obesity_adults)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Morbid Obesity", y = "Adult Morbid Obesity Prevalence") + 
  theme_economist()
Fig.13 # Figure 13 shows a box plot illustrating the prevalence of morbid obesity within adults
  
Fig.14 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_underweight_adults)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Underweight", y = "Adult Underweight Prevalence") +
  theme_economist()
Fig.14 # Figure 14 is a box plot that displays the prevelance of underweight adults
  
Fig.15 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_obesity_adults)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Obesity", y = "Adult Obesity Prevalence") +
  theme_economist()
Fig.15 # Figure 15 is a box plot detailing the prevalence of obese adults
  
Fig.16 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Mean_BMI_adults)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Mean BMI", y = "Mean BMI") +
  theme_economist()
Fig.16 # Figure 16 con sists of a box plot displaying the mean BMI within adults

ggarrange(Fig.13, Fig.14, Fig.15, Fig.16 + rremove("x.text"), 
          labels = c("Fig.11", "Fig.12", "Fig.13", "Fig.14"),
          ncol = 2, nrow = 2) # 

Fig.17 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Urbanisation)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Urbanisation", y = "Urbanisation Score (0-1)") +
  theme_economist()
Fig.17 # Figure 17 is a box plot that shows the urbanisation score of each country within the data
  
Fig.18 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Western_diet_score)) + 
  geom_boxplot() + 
  labs(title = "Western Diet Score", y = "Western Diet Score (-2.5-4.5)") +
  theme_economist()
Fig.18 # Figure 18 is aa box plot that displays the western diet score of each country within the data

ggarrange(Fig.15, Fig.16 + rremove("x.text"), 
          labels = c("Fig.17", "Fig.18"),
          ncol = 2) # This code arranges the above figures into one plot in order to see all the figures in one plot

##### Histograms ###############################################################
Fig.19 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Diabetes_prevalence)) + 
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Diabetes Prevalence", x = "Diabtes Prevalence") +
  theme_economist()
Fig.19 # Figure 19 is a histogram that displays the prevelance of diebetes within countries

Fig.20 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Systolic_blood_pressure)) + 
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Systolic blood Pressure", x = "Systolic blood Pressure") +
  theme_economist()
Fig.20 # Figure 20 is a histogram that illustrates systolic blood pressure within each country from the data 

Fig.21 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_raised_blood_pressure)) + # Prevalence Raised Blood Pressure
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Prevalence Raised Blood Pressure", x = "Prevalence Raised Blood Pressure") +
  theme_economist()
Fig.21 # Figure 21 is a histogram showing Prevalence Raised Blood Pressure from the cleaned data set.
Fig.22 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Years_of_education)) + # Years of Education
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Years of Education", x = "Years of Education") +
  theme_economist()
Fig.22 # Figure 21 is a histogram showing Years of education using the cleaned data set.
Fig.23 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$GDP_USD)) + # GDP
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "GDP($)", x = "GDP($)") +
  theme_economist()
Fig.23 # Figure 23 histogram displaying GDP in USD from the cleaned data sets.

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
  
Fig.24 # Figure 24 displays a box plot comparing Super Region and Years of education from the cleaned data set.

Fig.25<- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Diabetes_prevalence, x=  MyDataCleaned$Region)) + # Region and Diabetes prevalence 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Prevalence of Diabetes in each Region", y = "Disabetes Prevalence", x = "Region", tag = "Fig.25") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))
  
Fig.25 # Figure 25 displays a box plot comparing region with prevalence of diabetes using the cleaned data set.

Fig.26 <-ggplot(MyDataCleaned, aes(y = MyDataCleaned$Urbanisation, x=  MyDataCleaned$Region)) + # Region and Urbanisation 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Urbanisation of each Region", y = "Urbanisation Score (0-1)", x = "Region", tag = "Fig.26") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))

Fig.26 # Figure 26 displays a box plot that looks to compare Region with urbanisation.

Fig.27 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Western_diet_score, x=  MyDataCleaned$Superregion)) + # Super region and Western diet score
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Western Diet Score in each Super Region", y = "Western Diet Score (-2.5-4.5)", x = "Super Region", tag = "Fig.27") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))

Fig.27 # Figure 27 displays a box plot that compares Super Region and Western diet score 

#ggarrange allows us to display the four figures 24,25,26 and 27 which show compare region with Prevelance of Diabetes and Urbanisation as well as Super region against Years of Education and Western Diet score.
ggarrange(Fig.24, Fig.25, Fig.26,Fig.27 + rremove("x.text"), 
          ncol = 2, nrow = 2)

Fig.28 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Diabetes_prevalence, x=  MyDataCleaned$Sex)) + # Sex and Diabetes (box plot)
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Prevalence of Diabetes for each Sex", y = "Disabetes Prevalence", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()
  
Fig.28 # Figure 28 displays a box plot that compares Sex and prevalence of Diabetes.

Fig.29 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Years_of_education, x=  MyDataCleaned$Sex)) + # Sex and education (box plot)
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Years of Education for each Sex", y = "Years of Education", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()

Fig.29 # Figure 2 displays a box plot comparing both genders with average years of education. 

Fig.30 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_morbid_obesity_adults, x=  MyDataCleaned$Sex)) + # Sex and obesity (box plot)
  geom_boxplot() + 
  labs(title = "Prevalence of Adult Morbid Obesity for each Sex", y = "Morbid Obesity Prevalence", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()
  
Fig.30 # Figure 30 displays a box plot comparing Sex with prevalence of morbid adult obesity. 

Fig.31 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_underweight_adults, x=  MyDataCleaned$Sex)) + # Sex and obesity (box plot)
  geom_boxplot() + 
  labs(title = "Prevalence of Underweight Adults for each Sex", y = "Prevalence of Underweight Adults", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()
  
Fig.31 # Figure 31 illustrates a box plot comparing Sex with prevalence of underweight adults.

# ggarange illustrates Figures 28,29,30 and 31 which show Prevalence of Morbid Adult Obesity, Underweight Adults, Diabetes and Years of Education.

ggarrange(Fig.28, Fig.29, Fig.30, Fig.31 + rremove("x.text"), 
          labels = c("Fig.28", "Fig.29", "Fig.30", "Fig.31"),
          ncol = 2, nrow = 2)

##### Scatter plots ############################################################
Fig.32 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Western_diet_score, y = MyDataCleaned$Prevalence_obesity_children)) + # West diet score and Prevalence obesity children
  geom_point(colour = "#1380A1") +
  labs(title = "Western Diet and Prevalence of Child Obesity", y = "Prevalence of Child Obesity", x = "Western Diet Score (-2.5-4.5)", tag = "Fig.32") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.32 # Figure 32 displays a scatter plot showing Western Diet score being plotted against prevalence of childhood obesity.

Fig.33 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Western_diet_score, y = MyDataCleaned$Prevalence_obesity_adults)) + # West diet score and Prevalence obesity adults
  geom_point(colour = "#1380A1") +
  labs(title = "Western Diet and Prevalence of Adult Obesity", y = "Diabetes Prevalence", x = "Western Diet Score (-2.5-4.5)", tag = "Fig.33") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.33 # Figure 33 displays a scatter plot displaying Western Diet score against prevalence of adult obesity.

# ggarange shows us Figures 32 and 33 to show Western Diet score against childhood and adult obesity. 
ggarrange(Fig.32, Fig.33+ rremove("x.text"), 
          ncol = 2)


Fig.34 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Years_of_education, y = MyDataCleaned$Prevalence_obesity_children)) + # education and Prevalence obesity children
  geom_point(colour = "#1380A1") +
  labs(title = "Years of Education and Prevalence of Child Obesity", y = "Prevalence of Child Obesity", x = "Years of Education", tag = "Fig.34") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.34 # Figure 34 displays a scatter plot showing Years of Education against Prevelance of Child Obesity. 

Fig.35 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Years_of_education, y = MyDataCleaned$Prevalence_obesity_adults)) + # education and Prevalence obesity adults
  geom_point(colour = "#1380A1") +
  labs(title = "Years of Education and prevelance of Adult Obesity", y = "Prevalence of Adult Obesity", x = "Years of Education", tag = "Fig.35") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.35 # Figure 35 displays a scatter plot showing Years of Education against Prevelance of Adult Obesity.

#ggarrange allows us to display both Figure 34 and 35 to show Years of Education against both prevelance of Adult and Child Obesity.
ggarrange(Fig.34, Fig.35 + rremove("x.text"), 
          ncol = 2)

Fig.36 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_underweight_children , y = MyDataCleaned$GDP_USD)) + # GDP and Prevalence underweight children
  geom_point(colour = "#1380A1") +
  labs(title = "GDP($) and Prevalence of underweight Children", y = "GDP($)", x = "Prevalence of underweight Children", tag = "Fig.36") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.36 # Figure 36 displays a scatter plot showing GDP in USD compared with Prevelance of underweight Children.

Fig.37 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_underweight_adults, y = MyDataCleaned$GDP_USD)) + # GDP and Prevalence underweight adults
  geom_point(colour = "#1380A1") +
  labs(title = "GDP($) and Prevalence of underweight Adults", y = "GDP($)", x = "Prevalence of underweight Adults", tag = "Fig.37") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.37 # Figure 37 shows a scatter plot displaying GDP in USD against Prevelance of Underweight Adults.

#ggarrange allows us to display the two figures 36 and 37 together which show GDP in USD compared with both prevelance of underweight Adults and Children.
ggarrange(Fig.36, Fig.37+ rremove("x.text"), 
          ncol = 2)

Fig.38 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Urbanisation, y = MyDataCleaned$Prevalence_overweight_children)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Urbanisation and Prevalence of Overweight Children", y = "Prevalence of Oveweight Children", x = "Urbanisation Score (0-1)", tag = "Fig.38") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()
Fig.38 # Figure 38 displays a scatter plot showing Urbanisation against Prevelance of Overweight Children.

Fig.39 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Mean_BMI_adults, y = MyDataCleaned$Systolic_blood_pressure)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Sysolic Blood Pressure and Mean BMI (Adults)", y = "Systolic Blood Pressure", x = "Mean BMI", tag = "Fig.39") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()
Fig.39 # Figure 39 displays a scatter plot showing Systolic blood pressure against Mean BMI Adults.

Fig.40 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Mean_BMI_adults, y = MyDataCleaned$Prevalence_obesity_adults)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Prevalence of Obesity and Mean BMI (Adults)", y = "Prevalence of Obesity", x = "Mean BMI", tag = "Fig.40") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()
Fig.40 # Figure 40 displays a scatter plot showing Mean BMI against Prevalence obesity adults.

Fig.41 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_overweight_children, y = MyDataCleaned$Prevalence_obesity_adults)) + # Prevalence overweight children and Prevalence obesity adults
  geom_point(colour = "#1380A1") +
  labs(title = "Prevalence of Adult Obesity and Overweight Children", y = "Prevalence of Adult Obesity", x = "Prevalence of Overweight Children", tag = "Fig.41") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.41 # Figure 41 displays a scatter plot displaying Prevelance of Overweight Children against Prevelance of Adult Obesity.

Fig.42 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_underweight_children, y = MyDataCleaned$Prevalence_underweight_adults)) + # Prevalence underweight children and Prevalence underweight adults
  geom_point(colour = "#1380A1") +
  labs(title = "Prevalence of Underweight Adults and Children", y = "Prevalence of Underweight Adults", x = "Prevalence of Underweight Children", tag = "Fig.42") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.42 # Figure 42 displays a scatter plot showing Prevelance of Underweight Children against Prevelance of Underweight Adults.

#ggarrange allows us to display both Figure 41 and 42 which display Prevelance of Underweight Children against Prevelance of Underweight Adults as well as Prevelance of Overweight Children and Prevelance of Adult Obesity.
ggarrange(Fig.41, Fig.42+ rremove("x.text"), 
          ncol = 2)

##### Clustering ################################################################


##### Data Coding/Cleaning for clustering ######################################
df <- subset(MyDataCleaned, select =-c(1,2,3,4,5,17,18)) # We use this code in order to clean our data
df <- na.omit(df) # This code removes any missing data i.e. where there are any NAs within the large data set
df2 <- scale(df) # We scale the the variable in case one variable has a larger range than another variable
head(df2) # This code retrieves the first numbers of rows within the data frame - this helps to quickly test whether the data frame has the right type of data in it
str(df2) # This code displays the variables along with the observations within the data set to ensure we are using the right data

##### Euclidean distances between the points and visualizations ################
distance <- get_dist(df2) # This code measures the straight line distance between the rows of data
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) # This code helps us to visualise the distance matrix 

##### K-means ##################################################################
k2 <- kmeans(df2, centers = 2, nstart = 25)
k3 <- kmeans(df2, centers = 3, nstart = 25)
k4 <- kmeans(df2, centers = 4, nstart = 25)
k5 <- kmeans(df2, centers = 5, nstart = 25) # K-means clustering groups similar items in the form of clusters - it is recommended to do between two and five clusters 

p1 <- fviz_cluster(k2, geom = "point", data = df2) + ggtitle("k = 2") 
p2 <- fviz_cluster(k3, geom = "point",  data = df2) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df2) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 5") # These codes allow us to visualise the clusters in plots
grid.arrange(p1, p2, p3, p4, nrow = 2) # This code allows us to see all four k-means clusters in one plot in order to compare them all

##### Elbow Method #############################################################
set.seed(123) 
fviz_nbclust(df2, kmeans, method = "wss") # This code displays the elbow method graph - from the graph we can see that 4 is the optimal number of clusters 

##### Silhouette Method #########################################################
fviz_nbclust(df2, kmeans, method = "silhouette") # The code displays the silhouette method graph - from the graph we can see that 3 is the optimal number of clusters 

##### Gap Statistic ############################################################
gap_stat <- clusGap(df2,FUN = kmeans, nstart = 25,K.max = 10, B = 50)
fviz_gap_stat(gap_stat) # The code displays the gap statistic graph - from the graph we can see that 10 is the optimal number of clusters 

set.seed(123) # This code allows us to compute the gap statistic
gap_stat <- clusGap(df2, FUN = kmeans, nstart = 25, 
                    K.max = 10, B = 50) # This code allows us to visualize the results and we can see that four is the optimal number of clusters

##### Interpretation of the clusters ###########################################
final <- kmeans(df2, 4, nstart = 25)
fviz_cluster(final, data = df2) # Most of the approaches suggest that four is the optimal number of clusters

df %>% 
  mutate(Cluster = final$cluster) %>% ##########################################
  group_by(Cluster) %>% 
  summarise_all("mean") # This code extracts and add descriptive statistics at the cluster level

##### Hierarchical clustering ##################################################

##### arcsin transformation #################################################### 
arcsin_transformation <- function(x) asin(x/100) # There is no package for the arcsin transformation therefore, here we are creating our own function to make it easier to use

dend <- df2 %>% arcsin_transformation %>%
  dist %>% hclust(method = "com") %>% 
  as.dendrogram %>%
  set("branches_k_color", k = 3) %>% 
  ladderize # This code applies the arcsin transformation to the model using a clustering algorithm

dend2 <- dend %>% 
  set("labels", c(1:14000)) # This code organises the labels within the vectors 1 to 14000

dend2 %>% labels

##### Heat Map #################################################################
gplots::heatmap.2(as.matrix(df2), 
                  main = "Heat Map",srtCol = 60, dendrogram = "row",
                  Rowv = dend2,Colv = "NA", trace="none", margins = c(8,8), cexCol = 0.5,     
                  key.xlab = "Health Factors", denscol = "grey", density.info = "density",col = colorspace::diverge_hcl) # This code allows us to examine our results through the use of a plot - we decided that a heat map was the best way to illustrate our clusters 

##### Hierarchical clustering methods ##########################################
hclust_methods <- c("single", "complete", "average", "median")
df2_dendlist <- dendlist() 
for(i in seq_along(hclust_methods)) {
  tmp_dend <- df2 %>% arcsin_transformation %>% dist %>%hclust(method=hclust_methods[i]) %>% as.dendrogram 
  df2_dendlist<-dendlist(df2_dendlist,tmp_dend)
}
names(df2_dendlist) <- hclust_methods # This entire code allows us to look at the cophenetic correlation which is a measure of how similar the trees are to one another 

corrplot::corrplot(cor.dendlist(df2_dendlist),"pie","lower") # This code allows us to visualise the data from above using a correlation matrix

##### compare trees produced by different methods using a tanglegram ###########
dend1 <- df2 %>% arcsin_transformation %>% dist %>% hclust(method = "complete") %>% 
  as.dendrogram %>% color_branches(k=3) %>% ladderize # This code displays the first dendogram  by using the complete linkage method

dend2 <- df2 %>% arcsin_transformation %>% dist %>% hclust(method = "average") %>% 
  as.dendrogram %>% color_branches(k=3) %>% ladderize # The code display the second dendogram by using the mean linkage method

dends<-dendlist(tree1 = dend1, tree2 = dend2) 
dends %>% tanglegram(margin_inner = 7) # This code displays both dendograms in the tanglegram method - this helps us visualise btoh dendograms against one another 

##### Classification ###########################################################

##### Decision trees ###########################################################

##### Decision tree Sex ########################################################
#creating an initial tree using Sex as a response variable.
tree_default <- rpart(MyDataCleaned$Sex ~ ., data = MyDataCleaned)
#plotting the intial tree
rpart.plot(tree_default,extra=2, under = TRUE, varlen=0, faclen=0)

#creating a full tree that displays all variables necessary. 
tree_full <- rpart(MyDataCleaned$Sex ~., data= MyDataCleaned, control=rpart.control(minsplit=2, cp=0))
#plotting the full tree,
rpart.plot(tree_full,extra=2,under=TRUE,varlen=0,faclen=0,cex=.7)

#creating a confusion table to display the tree based on Sex.
confusion_table<-table(MyDataCleaned$Sex, 	predict(tree_default,MyDataCleaned,type="class"))
#viewing the confusion_table.
confusion_table

#calculating the number of predictions that are correct.
correct <- sum(diag(confusion_table))
#calculating the number of predictions that are incorrect.
error <- sum(confusion_table)-correct
#using both correct and incorrect values to calculate accuracy of predictions.
accuracy <- correct / (correct+error);accuracy
#creating a confusion matrix to display the predictions.
confusionMatrix(data= predict(tree_default,MyDataCleaned,type="class"), reference = MyDataCleaned$Sex)

##### Decision tree Super Region ###############################################
#creating an second tree that looks at Super Region.
tree_default2 <- rpart(MyDataCleaned$Superregion ~ ., data = MyDataCleaned)
#plotting the second tree.
rpart.plot(tree_default2,extra=2, under = TRUE, varlen=0, faclen=0)

#creating a second full tree that displays all variables for Super Region.
tree_full2 <- rpart(MyDataCleaned$Superregion ~., data= MyDataCleaned, control=rpart.control(minsplit=2, cp=0))
#plotting the second full tree.
rpart.plot(tree_ful2l,extra=2,under=TRUE,varlen=0,faclen=0,cex=.7)

#creating another confusion table based off the second initial tree for Super Region.
confusion_table2<-table(MyDataCleaned$Superregion, 	predict(tree_default2,MyDataCleaned,type="class"))
#showing the values of the second confusion table.
confusion_table2

#calculating the number of predictions based off the second confusion table that are correct.
correct2 <- sum(diag(confusion_table2))
#calculating the number of predictions based off the second confusion table that are incorrect.
error2 <- sum(confusion_table2)-correct2
#using both correct and incorrect values to calculate the accuracy of second confusion table.
accuracy2 <- correct2 / (correct2+error2);accuracy2

##### Random Forest ############################################################

##### Data Coding/Cleaning for classification ##################################
df3 <- MyDataCleaned[,c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,18)] 

##### Data Train and Test sets #################################################
set.seed(123) # Makes simulations random numbers the same to ensure all results, figures  are reproducible.

#splitting the data so that the models are using 70% of the data. 
split<-sample.split(df3, SplitRatio = 0.7)  
training_set<-subset(df3,split==TRUE)
test_set<-subset(df3,split==FALSE) 

#the dim function shows us the number of variables and data entries for both training and test sets for df3.
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

#splitting the data so that the models are using 70% of the data. 
split2<-sample.split(df4, SplitRatio = 0.7)  
training_set2<-subset(df4,split==TRUE)
test_set2<-subset(df4,split==FALSE) 
#the dim function shows us the number of variables and data entries for both training and test sets for df4.
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





```






















