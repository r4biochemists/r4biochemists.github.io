#Answers to R3 Problems in https://r4biochemists.github.io/
#By Blanca Lacruz, November 2024


setwd("C:/Users/Blanca.lacruz/HPBBM") #in Windows
setwd("~/Documents/GitHub/r4biochemists.github.io/")  #Mac

# 1. Load colis.csv and explore dataset structure ########

# load file
colis <- read.table("data/colis3.csv",
                    header = TRUE,
                    row.names = 1,
                    sep = ";")

# we could also use read.csv2:
# colis <- read.csv2("colis3.csv",
#                     row.names = 1)

# explore structure
dim(colis) # number of rows and columns
head(colis) # first 6 rows
tail(colis) # last 6 rows
str(colis)


# 2. Calculate means ########################################

mean(colis$Year) # there are missing values!


# Dealing with NA's:
mean(colis$Year, na.rm = TRUE) # the easiest way
mean(colis$Year[!is.na(colis$Year)]) # subsetting with is.na
mean(colis$Year[complete.cases(colis$Year)]) # subsetting with complete.cases

# OK, but what if there are NAs in other variables too?

coli_genomes_renamed <- na.omit(colis) # colis, but without NA's
dim(coli_genomes_renamed) # we lose 4 rows

mean(coli_genomes_renamed$Year)


## Other columns ##########
# A quick-and-dirty solution using Apply() - see Lesson R6:
apply(coli_genomes_renamed, 2, function(x) mean(as.numeric(x)))

# Only in numeric columns:
apply(coli_genomes_renamed[, sapply(coli_genomes_renamed, is.numeric)], 
      2, 
      function(x) mean(as.numeric(x)))




# Mean sequencing date:
mean(coli_genomes_renamed$seqs) # doesn't work
mean(as.Date(coli_genomes_renamed$seqs))

# 3. Save RData ##############################
save(colis, coli_genomes_renamed, file = "./coli_tables.Rdata")


# 4. Add means as last row #############
#Mean calculation without apply(), one by one
year <- mean(coli_genomes_renamed$Year)
crispr <- mean(coli_genomes_renamed$crispr)
amr <- mean(coli_genomes_renamed$AMR)
vf <- mean(coli_genomes_renamed$Vir) #na.rm not needed, but no problem
integron <- mean(coli_genomes_renamed$Integron)
seq <- mean(as.Date(coli_genomes_renamed$seqs))
#if we don't want to remove the NAs
year <- mean(colis$Year,na.rm=TRUE)
crispr <- mean(colis$crispr,na.rm=TRUE)
amr <- mean(colis$AMR,na.rm=TRUE)
vf <- mean(colis$Vir,na.rm=TRUE) #na.rm not needed here, but I copy/paste code
integron <- mean(colis$Integron,na.rm=TRUE)
seq <- mean(as.Date(colis$seqs))

#add the data in the last row without knowing how many rows
coli_genomes_renamed[nrow(coli_genomes_renamed)+1,] <- c("mean",year,crispr,NA,amr,vf,integron,seq)
#alternative with rbind()
coli_genomes_renamed <- rbind(coli_genomes_renamed, c("mean",year,crispr,NA,amr,vf,integron,seq))

coli_genomes_renamed <- rbind(coli_genomes_renamed, means)
tail(coli_genomes_renamed)


# we could save it again:
save(colis, coli_genomes_renamed, file = "./coli_tables.Rdata")
