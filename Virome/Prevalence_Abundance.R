#Linear mixed modeling for contig prevalence and abundance preceding NEC onset (Supplementary Figure 3).
#Run model first with interaction term to check for interaction between time and case/control status. If there is a significant interaction, run model separately for cases and controls.

library(nlme)

#Prevalence
data <- read.delim("ContigPrevalence.txt")
nlme.test <- lme(Prevalence ~ TimeFromNEC + Case_or_Control + (TimeFromNEC*Case_or_Control), random =~1|ContigName, data = data)
summary(nlme.test)

data <- read.delim("ContigPrevalence_controls.txt")
nlme.test <- lme(Prevalence ~ TimeFromNEC, random =~1|ContigName, data = data)
summary(nlme.test)

#Average abundance
data <- read.delim("ContigAverageAbundance.txt")
nlme.test <- lme(Abundance ~ TimeFromNEC + Case_or_Control + (TimeFromNEC*Case_or_Control), random =~1|ContigName, data = data)
summary(nlme.test)

data <- read.delim("ContigAverageAbundance_controls.txt")
nlme.test <- lme(Abundance ~ TimeFromNEC, random =~1|ContigName, data = data)
summary(nlme.test)

#Summed abundance
data <- read.delim("ContigSummedAbundance.txt")
nlme.test <- lme(SummedAbundance ~ TimeFromNEC + Case_or_Control + (TimeFromNEC*Case_or_Control), random =~1|InfantID, data = data)
summary(nlme.test)

data <- read.delim("ContigSummedAbundance_controls.txt")
nlme.test <- lme(SummedAbundance ~ TimeFromNEC, random =~1|InfantID, data = data)
summary(nlme.test)
