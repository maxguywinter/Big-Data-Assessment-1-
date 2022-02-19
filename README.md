### Install Packages 
library(corrplot) # for correlation matrix graph visualization.
library(ggplot2) # for graph visualizations. 
library(ggthemr)
library(ggthemes)
install.packages('hrbrthemes')
library(hrbrthemes)
library(devtools)
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
library(janitor)
library(ggpubr)

### Install Data Set 
data <- read.csv("IntOrg_NCD_variables_2022_02_02.csv", header = TRUE) # imports the data set. The header=True command tells RStudio to use the first row of the data file as the names of each variable/column. 
attach(data) # attaches the data to your environment so that you can directly refer to the variable by name.
names(data) # shows the name of variables in the data set.
head(data)
summary(data)
str(data) # shows the observations and variables of the data.

### Exploratory Data Analysis

### Coding/Data Cleaning
data[data == "" | data == " "] <- NA
data[data == "N/A" | data == "n/a" | data == "N/a"| data == "-"] <- NA

### Changing categorical variables to factors 
data$Country <- factor(data$Country)
data$ISO <- factor(data$ISO)
data$Sex <- factor(data$Sex)
data$Region <- factor(data$Region)
data$Superregion<- factor(data$Superregion)

### Missing Data 
missmap(data, col=c("red", "#1380A1"), x.cex = 0.3,legend = TRUE) # Checks for missing data. 
pct_miss(data) # Percent of ALL data frame values that are missing
pct_miss_case(data) # Percent of rows with any value missing
pct_complete_case(data) # Percent of rows that are complete (no values missing) 

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

### Check for duplicates and wrong input 
get_dupes(MyDataCleaned)
unique(Country)
unique(ISO)
unique(Sex)
unique(Region)
unique(Superregion)

# Graphical visualization 

## Correlation Analysis (correlation matrix and scatter plot matrix)
datacorrelation <- subset(MyDataCleaned, select =-c(1,2,3,4,5,17,18)) # subset data to remove factors. 
nums <- unlist(lapply(datacorrelation, is.numeric)) # used to create correlation plot from numeric values
Datacorr <- cor(datacorrelation[, nums])
corrplot(Datacorr, type = "upper", method = "number", number.cex = 0.5, tl.cex = 0.5, title = "\n\n Correlation Plot \n") # correlation matrix graph visualization. 
pairs(datacorrelation) # graphical matrix visualization.

## lollipop chart
Fig.5 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$ISO, y = MyDataCleaned$Year)) + # Countries (ISO) and Years
  geom_point(size = 2, color = "#1380A1" ) +
  ggtitle("Country (ISO)", subtitle = "Country data aviable for each year") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "Countries (ISO)", y = "Years", tag = "Fig.5") +
  theme_economist()
Fig.5 

## Bar graph 
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

## Box plots 

### Child weights 
Fig.9 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_obesity_children)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Child Obesity", y = "Child Obesity Prevalence") +
  theme_economist()

Fig.9 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_overweight_children)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Child Overweight", y = "Child Overweight Prevalence") +
  theme_economist()

Fig.10 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_underweight_children)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Child Underweight", y = "Child Underweight Prevalence") +
  theme_economist()

Fig.11 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Mean_BMI_children)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Children Mean BMI", y = "Mean BMI") +
  theme_economist()

ggarrange(Fig.8, Fig.9, Fig.10, Fig.11 + rremove("x.text"), 
          labels = c("Fig.7", "Fig.8", "Fig.9", "Fig.10"),
          ncol = 2, nrow = 2)

### Adult weights 
Fig.12 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_morbid_obesity_adults)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Morbid Obesity", y = "Adult Morbid Obesity Prevalence") +
  theme_economist()

Fig.13 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_underweight_adults)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Underweight", y = "Adult Underweight Prevalence") +
  theme_economist()

Fig.14 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_obesity_adults)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Obesity", y = "Adult Obesity Prevalence") +
  theme_economist()

Fig.15 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Mean_BMI_adults)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Adult Mean BMI", y = "Mean BMI") +
  theme_economist()

ggarrange(Fig.12, Fig.13, Fig.14, Fig.15 + rremove("x.text"), 
          labels = c("Fig.11", "Fig.12", "Fig.13", "Fig.14"),
          ncol = 2, nrow = 2)


Fig.16 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Urbanisation)) + # urbanization 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Urbanisation", y = "Urbanisation Score (0-1)") +
  theme_economist()
Fig.17 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Western_diet_score)) + # Western diet score 
  geom_boxplot() + 
  labs(title = "Western Diet Score", y = "Western Diet Score (-2.5-4.5)") +
  theme_economist()
ggarrange(Fig.15, Fig.16 + rremove("x.text"), 
          labels = c("Fig.16", "Fig.17"),
          ncol = 2)

# Hist
Fig.18 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Diabetes_prevalence)) + # Diabetes Prevalence
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Diabetes Prevalence", x = "Diabtes Prevalence") +
  theme_economist()
Fig.18

Fig.19 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Systolic_blood_pressure)) + # Systolic blood Pressure
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Systolic blood Pressure", x = "Systolic blood Pressure") +
  theme_economist()
Fig.19

Fig.20 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_raised_blood_pressure)) + # Prevalence Raised Blood Pressure
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Prevalence Raised Blood Pressure", x = "Prevalence Raised Blood Pressure") +
  theme_economist()
Fig.20

Fig.21 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Years_of_education)) + # Years of Education
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "Years of Education", x = "Years of Education") +
  theme_economist()
Fig.21

Fig.22 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$GDP_USD)) + # GDP
  geom_histogram(colour = 4, fill = "#1380A1", 
                 bins = 30) +
  labs(title = "GDP($)", x = "GDP($)") +
  theme_economist()
Fig.22

#Comparisons 

# Regions
Fig.23 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Years_of_education, x=  MyDataCleaned$Superregion)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Years of Education for each Super Region", y = "Years of Education", x = "Super Region", tag = "Fig.23") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))
Fig.23

Fig.24<- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Diabetes_prevalence, x=  MyDataCleaned$Region)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Prevalence of Diabetes in each Region", y = "Disabetes Prevalence", x = "Region", tag = "Fig.24") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))
Fig.24

Fig.25 <-ggplot(MyDataCleaned, aes(y = MyDataCleaned$Urbanisation, x=  MyDataCleaned$Region)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Urbanisation of each Region", y = "Urbanisation Score (0-1)", x = "Region", tag = "Fig.25") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))
Fig.25

Fig.26 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Western_diet_score, x=  MyDataCleaned$Superregion)) +
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Westen Diet Score each Super Region", y = "Western Diet Score (-2.5-4.5)", x = "Super Region", tag = "Fig.26") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist() + scale_colour_economist() +
  theme(axis.text = element_text(size = 5))

ggarrange(Fig.23, Fig.24, Fig.25,Fig.26 + rremove("x.text"), 
          ncol = 2, nrow = 2)

# Sex 

Fig.27 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Diabetes_prevalence, x=  MyDataCleaned$Sex)) + # Sex and Diabetes (box plot)
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Prevalence of Diabetes for each Sex", y = "Disabetes Prevalence", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()

Fig.28 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Years_of_education, x=  MyDataCleaned$Sex)) + # Sex and education (box plot)
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot() + 
  labs(title = "Years of Education for each Sex", y = "Years of Education", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()

Fig.29 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_morbid_obesity_adults, x=  MyDataCleaned$Sex)) + # Sex and obesity (box plot)
  geom_boxplot() + 
  labs(title = "Prevalence of Adult Morbid Obesity for each Sex", y = "Morbid Obesity Prevalence", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()

Fig.30 <- ggplot(MyDataCleaned, aes(y = MyDataCleaned$Prevalence_underweight_adults, x=  MyDataCleaned$Sex)) + # Sex and obesity (box plot)
  geom_boxplot() + 
  labs(title = "Prevalence of Underweight Adults for each Sex", y = "Prevalence of Underweight Adults", x = "Sex") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  theme_economist()

ggarrange(Fig.27, Fig.28, Fig.29, Fig.30 + rremove("x.text"), 
          labels = c("Fig.27", "Fig.28", "Fig.29", "Fig.30"),
          ncol = 2, nrow = 2)

# West diet score and Obese child/adult 
Fig.31 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Western_diet_score, y = MyDataCleaned$Prevalence_obesity_children)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Western Diest and Prevalence of Child Obesity", y = "Prevalence of Child Obesity", x = "Western Diet Score (-2.5-4.5)", tag = "Fig.31") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.32 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Western_diet_score, y = MyDataCleaned$Prevalence_obesity_adults)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Western Diest and Prevalence of Adult Obesity", y = "Disabetes Prevalence", x = "Western Diet Score (-2.5-4.5)", tag = "Fig.32") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

ggarrange(Fig.31, Fig.32+ rremove("x.text"), 
          ncol = 2)

# education and Obesity child/adult 
Fig.33 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Years_of_education, y = MyDataCleaned$Prevalence_obesity_children)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Years of Education and Prevalence of Child Obesity", y = "Prevalence of Child Obesity", x = "Years of Education", tag = "Fig.33") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.34 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Years_of_education, y = MyDataCleaned$Prevalence_obesity_adults)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Western Diest and Prevalence of Adult Obesity", y = "Prevalence of Adult Obesity", x = "Years of Education", tag = "Fig.34") +
  scale_x_discrete(guide = guide_axis(n.dodge=10)) +
  guides(x = guide_axis(angle = 90)) +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

ggarrange(Fig.33, Fig.34 + rremove("x.text"), 
          ncol = 2)

# GDP and underweight child/adult
Fig.35 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_underweight_children , y = MyDataCleaned$GDP_USD)) +
  geom_point(colour = "#1380A1") +
  labs(title = "GDP($) and Prevalence of underweight Children", y = "GDP($)", x = "Prevalence of underweight Children", tag = "Fig.35") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.36 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_underweight_adults, y = MyDataCleaned$GDP_USD)) +
  geom_point(colour = "#1380A1") +
  labs(title = "GDP($) and Prevalence of underweight Adults", y = "GDP($)", x = "Prevalence of underweight Adults", tag = "Fig.36") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

ggarrange(Fig.35, Fig.36+ rremove("x.text"), 
          ncol = 2)

# Urbanization and Prevalence of Overweight Children

Fig.37 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Urbanisation, y = MyDataCleaned$Prevalence_overweight_children)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Urbanisation and Prevalence of Overweight Children", y = "Prevalence of Oveweight Children", x = "Urbanisation Score (0-1)", tag = "Fig.37") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()
Fig.37

# Systolic blood pressure and Mean BMI Adults
Fig.38 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Mean_BMI_adults, y = MyDataCleaned$Systolic_blood_pressure)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Sysolic Blood Pressure and Mean BMI (Adults)", y = "Systolic Blood Pressure", x = "Mean BMI", tag = "Fig.38") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()
Fig.38

# Adult Mean BMI and Obesity 
Fig.39 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Mean_BMI_adults, y = MyDataCleaned$Prevalence_obesity_adults)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Prevalence of Obesity and Mean BMI (Adults)", y = "Prevalence of Obesity", x = "Mean BMI", tag = "Fig.39") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()
Fig.39

# Child and Adult weights 

Fig.40 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_overweight_children, y = MyDataCleaned$Prevalence_obesity_adults)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Prevalence of Adult Obesity and Overweight Children", y = "Prevalence of Adult Obesity", x = "Prevalence of Overweight Children", tag = "Fig.40") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

Fig.41 <- ggplot(MyDataCleaned, aes(x = MyDataCleaned$Prevalence_underweight_children, y = MyDataCleaned$Prevalence_underweight_adults)) +
  geom_point(colour = "#1380A1") +
  labs(title = "Prevalence of Underweight Adults and Children", y = "Prevalence Underweight Adults", x = "Prevalence Underweight Children", tag = "Fig.41") +
  geom_smooth(method='lm', formula= y~x, se= FALSE, colour = "red") +
  theme_economist()

ggarrange(Fig.40, Fig.41+ rremove("x.text"), 
          ncol = 2)

# Clustering 

### Install Packages 
library(gridExtra)
library(tidyverse)  
library(cluster)   
install.packages("factoextra")
library(factoextra) 

### Data set for clustering 
df <- subset(MyDataCleaned, select =-c(1,2,3,4,5,17,18))
df <- na.omit(df)
df <- scale(df)
head(df)
str(df)

distance <- get_dist(df)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# K-means 
k2 <- kmeans(df, centers = 2, nstart = 25)
k3 <- kmeans(df, centers = 3, nstart = 25)
k4 <- kmeans(df, centers = 4, nstart = 25)
k5 <- kmeans(df, centers = 5, nstart = 25)
# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = df) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = df) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df) + ggtitle("k = 5")
grid.arrange(p1, p2, p3, p4, nrow = 2)

# Elbow Method
set.seed(123)
fviz_nbclust(df, kmeans, method = "wss") # 4? clusters k

# Silhoutte Method
fviz_nbclust(df, kmeans, method = "silhouette") # 3 clusters k 

# Gap Statistic
gap_stat <- clusGap(df,FUN = kmeans, nstart = 25,K.max = 10, B = 50)
fviz_gap_stat(gap_stat) # 10? clusters k

# compute gap statistic
set.seed(123)
gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)

final <- kmeans(df, 4, nstart = 25)
fviz_cluster(final, data = df)

df %>% 
  mutate(Cluster = final$cluster) %>% 
  group_by(Cluster) %>% 
  summarise_all("mean") 

















