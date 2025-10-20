#Answers to R6 Exercises in https://r4biochemists.github.io/
#Nov, 24th 2024

########################
###    Exercise 1    ###
########################

#Read data
library(data.table)
c19test <- fread("https://opendata.ecdc.europa.eu/covid19/testing/csv/data.csv")
#check it out
str(c19test)
head(c19test)
tabla1 <- xtabs(tests_done ~ as.factor(country) + as.factor(year_week),c19test)
#a flat table is more useful
ftable(tabla1,row.vars=1:2)

#second table

# install.packages("ISOweek") #install if needed
library(ISOweek)
#addapt the format to standard ISO 8601, adding the day
c19test$date<-paste(c19test$year_week,"1",sep="-")
#transform ISO to date
c19test$date<-as.Date(ISOweek2date(c19test$date))
#re-format to months
c19test$month <- format(as.Date(c19test$dat), "%b-%Y")
#xtabulate
tabla2<-xtabs(as.integer(tests_done)~as.factor(country)+as.factor(month),c19test)
#flat table
head(ftable(tabla2,row.vars=1:2))

########################
###    Exercise 2    ###
########################
c19test <- fread("https://opendata.ecdc.europa.eu/covid19/testing/csv/data.csv")

#Two options
#1. you can drop character columns...
minitest <- subset(c19test, select = !sapply(c19test, is.character) )
head(minitest)
#... or keep only integer ones
minitest2 <- subset(c19test, select = sapply(c19test, is.numeric) | sapply(c19test, is.integer))
head(minitest2)
#is the same?
all.equal(minitest2,minitest)


########################
###    Exercise 3    ###
########################
apply(minitest,2,summary)
#alternative
ex3 <- apply(minitest,2,
             function(x)c(min(na.omit(x)),max(na.omit(x)),median(na.omit(x)),mean(na.omit(x)))) 
row.names(ex3) <- c("Min","Max","Median","Mean")
ex3

########################
###    Exercise 4    ###
########################

c19test2 <- na.omit(c19test)
mean <- tapply(c19test2$new_cases,c19test2$country, mean) # mean per country
median <- tapply(c19test2$new_cases,c19test2$country, median) #median per country
min <- tapply(c19test2$new_cases,c19test2$country, min) #median per country
df <- data.frame(cbind(min,mean,median)) #create the dataframe
colnames(df) <- c("Min", "Mean","Median") 
rownames(df) <- names(mean) #names are the same in the three vectors
