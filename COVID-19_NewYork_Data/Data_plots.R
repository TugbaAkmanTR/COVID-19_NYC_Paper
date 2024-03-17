# Load packages
library(lava)
library(mice)
library(data.table)
library(ggplot2)

#Type getwd() to the Console to get the address
setwd("G:/Drive'ım/Preprint/InfectDiseases/COVID-19_NewYork_Data")
# To import this datafile 
datafull <- read.csv(file = 'COVID-19_Daily_Counts_of_Cases__Hospitalizations__and_Deaths.csv')
str(datafull)
summary(datafull)


# Extract the colums for case count, hospitalized count and death count
length_data_file = 290 # until 14 Dec 2021
data_Cases_Hosp_Death <- fread("COVID-19_Daily_Counts_of_Cases__Hospitalizations__and_Deaths.csv", 
                  select = c("date_of_interest","CASE_COUNT", "HOSPITALIZED_COUNT", "DEATH_COUNT")) 


#only import rows where date_of_interest > 90
data_Cases_Hosp_Death_2020 = head(data_Cases_Hosp_Death,length_data_file)

str(data_Cases_Hosp_Death_2020)
summary(data_Cases_Hosp_Death_2020)

# Extract the data and write it to a csv file

data_dates = data_Cases_Hosp_Death_2020[,1]

data_Cases <- fread("COVID-19_Daily_Counts_of_Cases__Hospitalizations__and_Deaths.csv", 
                               select = c("date_of_interest","CASE_COUNT")) 
data_Cases = head(data_Cases,length_data_file)

data_Hosp <- fread("COVID-19_Daily_Counts_of_Cases__Hospitalizations__and_Deaths.csv", 
                               select = c("date_of_interest", "HOSPITALIZED_COUNT"))
data_Hosp = head(data_Hosp,length_data_file)

data_Death <- fread("COVID-19_Daily_Counts_of_Cases__Hospitalizations__and_Deaths.csv", 
                               select = c("date_of_interest", "DEATH_COUNT"))
data_Death = head(data_Death,length_data_file)

write.csv(as.data.frame(data_dates),"G:/Drive'ım/Preprint/InfectDiseases/COVID-19_NewYork_Data/Dates.csv")
write.csv(as.data.frame(data_Cases),"G:/Drive'ım/Preprint/InfectDiseases/COVID-19_NewYork_Data/Cases.csv")
write.csv(as.data.frame(data_Hosp),"G:/Drive'ım/Preprint/InfectDiseases/COVID-19_NewYork_Data/Hosp.csv")
write.csv(as.data.frame(data_Death),"G:/Drive'ım/Preprint/InfectDiseases/COVID-19_NewYork_Data/Death.csv")


# Plots
ggplot(data_Cases, aes(x = date_of_interest, y = CASE_COUNT)) +
  geom_point()

ggplot(data_Cases_Hosp_Death_2020,
  aes(x = date_of_interest)) +
  geom_point(aes(y = CASE_COUNT), color = "green") +
  geom_point(aes(y = HOSPITALIZED_COUNT), color = "blue") +
  geom_point(aes(y = DEATH_COUNT), color = "red") +
  geom_density(alpha = 0.5)



  