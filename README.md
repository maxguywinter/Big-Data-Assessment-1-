# Big-Data-Assessment-1-
#installing relevant packages that will be used to analyse the data provided.
install.packages("Amelia)
library(Amelia) #this package will enable us to use the missmap function to find missing variables within the dataset.



missmap(IntOrg_NCD_variables_2022_02_02, col=c("black", "purple"), legend = FALSE)
