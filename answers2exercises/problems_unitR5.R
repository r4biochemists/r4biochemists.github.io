#Answers to R5 Problems in https://r4biochemists.github.io/

########################
###    Exercise 1    ###
########################

#Data generation
(GeneA <- rnorm(50))
(GeneB <- c(rep(-1, 30), rep(2, 20)) + rnorm(50))
(tumor <- factor(c(rep("Colon", 30), rep("Lung", 20))))

#We generated the data to be normal, thus we can use t.test
(testA <- t.test(GeneA ~ tumor)) #p>0.01, we accept the Null Hypothesis
(testb <- t.test(GeneB ~ tumor)) #p>0.01, we reject the Null Hypothesis, i.e., there is significative diference
#Let's explore the "test objects"
str(testA)
str(testB)

#Plots
par(mfrow = c(2, 2)) ## a 2-by-2 layout 
# gene A 
boxplot(GeneA ~ tumor, col = c("orange", "red"), xlab = "Type of tumor", 
        ylab = "Gene A relative expression", main = "Expression of GeneA") 
plot(density(GeneA[tumor == "Lung"]), col = "red", xlab = "Gene A relative expression", ylab = "Density curve", main = "Expression of GeneA", ylim = c(0, 0.5), xlim = c(-4, 4), ) 
polygon(density(GeneA[tumor == "Lung"]), col = rgb(1, 0, 0, 0.5), border = rgb(1, 0, 0, 0.5)) 
lines(density(GeneA[tumor == "Colon"]), col = "orange") 
polygon(density(GeneA[tumor == "Colon"]), col = rgb(1, 0.6, 0, 0.5), border = rgb(1, 0.6, 0, 0.5)) 
legend("topright", fill = c("red", "orange"), legend = c("Lung", "Colon"))
# gene B 
boxplot(GeneB ~ tumor, col = c("blue", "green"), xlab = "Type of tumor", ylab = "Gene B relative expression", main = "Expression of GeneB") 
plot(density(GeneB[tumor == "Lung"]), col = "green", ylim = c(0, 0.5), xlim = c(-4.5, 5), xlab = "Gene A relative expression", ylab = "Density curve", main = "Expression of GeneB") 
polygon(density(GeneB[tumor == "Lung"]), col = rgb(0, 1, 0, 0.5), border = rgb(0, 1, 0, 0.5)) 
lines(density(GeneB[tumor == "Colon"]), col = "blue") 
polygon(density(GeneB[tumor == "Colon"]), col = rgb(0, 0, 1, 0.5), border = rgb(0, 0, 1, 0.5)) 
legend("topright", fill = c("green", "blue"), legend = c("Lung", "Colon"), horiz = TRUE)

#What do you see?
#The plots show that the expression level of geneA is similar in both cancer types, but that’s not the case for geneB
par(mfrow = c(1, 1))

########################
###    Exercise 2    ###
########################
setwd("~/MEGAsync/teaching_2024_25/2024-25_HPBBM/")
uniprot <-read.csv("data/aapc.csv",sep= ";", dec= ".") 
# order the data 
uniprot <-uniprot[order(uniprot$Percent,decreasing= TRUE), ] 
# color by category 
colorines <- c(Aliphatic= "Grey", Basic= "Red", Amide = "Yellow", Acidic= "Blue", Sulfur= "Green", Aromatic ="Orange", Hydroxy= "Brown")

# two alternatives for order the data
#1:using reorder() 
barplot(uniprot$Percent ~ reorder(uniprot$Name3, uniprot$Percent, decreasing= TRUE),ylab= "Protein aminoacid composition (%)", xlab= "", las= 2, ylim= c(0, 10),col= colorines) 
legend("topright",fill= unique(colorines), legend= unique(names(colorines)), ncol= 2)

#2: convert the variable to a factor, making sure that the order in the dataframe is kept
uniprot$Name3<-factor(uniprot$Name3, levels= uniprot$Name3) 
barplot(uniprot$Percent ~ uniprot$Name3, ylab= "Protein aminoacid composition (%)", xlab= "", las= 2, ylim= c(0, 10),col= colorines[uniprot$Type]) 
legend("topright",fill= colorines, legend= names(colorines), ncol= 2)

########################
###    Exercise 3    ###
########################
# read the data 
(bsa <- read.csv2("data/bsa_pattern.csv") )
(wt <- read.csv2("data/curve_wt.csv") )
(mut <- read.csv2("data/curve_mut.csv"))

# there are many columns but only two variables, so we stack the bsa 
bsa2 <- cbind(bsa$ng, stack(bsa[, 2:4])) 
names(bsa2) <- c("ng", "AUC", "exp")
bsa2
#wait! stack is in Lesson R6, let's think in an alternative
bsa2b <- data.frame(ng=rep(bsa$ng,3),AUC=c(bsa$AUC1,bsa$AUC2,bsa$AUC3),exp=c(rep("AUC1",8),rep("AUC2",8),rep("AUC3",8)))
bsa2b


# linear model & correlation test
bsa_model <- lm(ng ~ AUC,data=bsa2) 
summary(bsa_model)
(test1 <- cor.test(bsa2$ng, bsa2$AUC))

#plot
colorines2 <- c(AUC1 = rgb(1, 0, 0, 0.5), AUC2 = rgb(0, 1, 0, 0.5), AUC3 = rgb(0, 0, 1, 0.5))
plot(bsa2$ng ~ bsa2$AUC, pch = 19, col = colorines2[bsa2$exp], ylab = "BSA (ng)", xlab = "AUC (arbitrary units)") 
abline(bsa_model, col = "red") 
text(5500000, 700, paste("Pearson R²=", round(test1$estimate, 2))) 
legend("bottomright", fill = colorines2, legend = names(colorines2), ncol = 2)                

#now we calculate the concentration

##Option 1: use the model coefficients to calculate
# wt amount in each well 
wt$ng <- (bsa_model$coefficients[2] * wt$AUC) + bsa_model$coefficients[1] 
# concentration in ng/µL 
wt$conc <- wt$ng/wt$μl # average concentration 
paste("The concentration of the WT protein is", round(mean(wt$conc/99), 2), "µM.")
# same for mutant 
mut$ng <- bsa_model$coefficients[2] * mut$AUC + bsa_model$coefficients[1] 
mut$conc <- mut$ng/mut$μl
paste("The concentration of the mutant protein is", round(mean(mut$conc/99), 2), "µM.")

##Option 2: use the model to interpolate directly
wt$ng2<- predict(bsa_model,data.frame(AUC=wt$AUC))
wt$conc2 <- wt$ng2/wt$μl 
paste("The concentration of the WT protein is", round(mean(wt$conc2)/99, 2), "µM.")

mut$ng2<- predict(bsa_model,data.frame(AUC=wt$AUC))
mut$conc2 <- wt$ng2/wt$μl
paste("The concentration of the mutant protein is", round(mean(mut$conc2)/99, 2), "µM.")

