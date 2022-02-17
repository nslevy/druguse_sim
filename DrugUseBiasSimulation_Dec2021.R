###############################################
# Simulation of Misclassification of Drug Use #
###############################################

# Load needed packages
library(tidyverse)
library(survey, quietly = T)
library(srvyr, quietly= T)
library(cowplot)

setwd("/PATH")
options(scipen=999)
rm(list = ls())

##########################################
# 1 - Create simulated demographics data #
##########################################

# Read in 2018 Census ACS Data Query, cleaned up a little bit in Excel because Census MDAT output is v. messy
# Create labeled factor versions of numeric demographic variables
raw.data <- 
  read.csv("/PATH/Census Custom Over 18_Clean.csv") %>%
  mutate(SEX_FACT = factor(SEX, 
                           labels = c("Male","Female")),
         AGE_FACT = factor(AGE, 
                           labels = c("18-25","26-34","35-49","50+"),
                           ordered = TRUE),
         EDUC_FACT = factor(EDUC, 
                            labels = c("Less than highschool","High school","Some college","Associates","College or more"))) %>%
  arrange(SEX, AGE, EDUC)

# Convert Census data into proportions
# Create rounded aggregated totals for simulated dataset of n=100,000
# Turn into long dataset
# Create demographic variables that will match with NSDUH categories
# Create demographic variable flag for later merging purposes

total <- sum(raw.data[,4:11])

census.data <- 
  raw.data %>%
  mutate_at(c("HISP","WHITE","BLACK","AI.AN","ASIAN","NH.PA","OTHER","MULT"), ~./total) %>%
  mutate_at(c("HISP","WHITE","BLACK","AI.AN","ASIAN","NH.PA","OTHER","MULT"), ~.*100000) %>%
  mutate_at(c("HISP","WHITE","BLACK","AI.AN","ASIAN","NH.PA","OTHER","MULT"), round) %>%
  pivot_longer(cols = c("HISP","WHITE","BLACK","AI.AN","ASIAN","NH.PA","OTHER","MULT"), names_to = "RACE_CHAR") %>%
  mutate(RACE = as.integer(case_when(
             RACE_CHAR == "WHITE" ~ 1,
             RACE_CHAR == "BLACK" ~ 2,
             RACE_CHAR == "AI.AN" ~ 3,
             RACE_CHAR == "NH.PA" ~ 4,
             RACE_CHAR == "ASIAN" ~ 5,
             RACE_CHAR == "MULT" ~ 6,
             RACE_CHAR == "HISP" ~ 7,
             RACE_CHAR == "OTHER" ~ 8)),
    RACE_FACT = as.factor(RACE_CHAR),
    RACE2 = ifelse(RACE==8, 6, RACE),
    DEMFLAG = paste(SEX, AGE, EDUC, RACE, sep = "_"),
    DEMFLAG2 = paste(SEX, AGE, EDUC, RACE2, sep = "_")) %>%
  select(SEX, SEX_FACT, AGE, AGE_FACT, EDUC, EDUC_FACT, RACE, RACE_FACT, DEMFLAG, DEMFLAG2, value) %>%
  arrange(SEX, AGE, EDUC, RACE)


###################################################################
# 2 - Create self-reported drug use variables based on 2018 NSDUH #
###################################################################

# Read in and recode 2018 NSDUH Data
nsduh2018 <-
  read.csv(file = "/PATH/NSDUH_2018_Tab.tsv", sep = "\t") %>%
  rename_all(toupper) %>%
  filter(CATAG3 > 1) %>%
  select(QUESTID2, IRSEX, IREDUHIGHST2, CATAG3, NEWRACE2, ILLYR, MRJYR, MRJMON, COCYR, COCMON, 
         METHAMYR, METHAMMON, ANALWT_C, VESTR, VEREP) %>%
  rename(SEX = IRSEX,
         AGE = CATAG3,
         RACE = NEWRACE2) %>%
  mutate(
    EDUC = as.integer(case_when(
      IREDUHIGHST2 < 8 ~ 1,
      IREDUHIGHST2 == 8 ~ 2,
      IREDUHIGHST2 == 9 ~ 3,
      IREDUHIGHST2 == 10 ~ 4,
      IREDUHIGHST2 == 11 ~ 5)),
    DEMFLAG2 = paste(SEX, AGE, EDUC, RACE, sep = "_"))

# Create survey design object
dnsduh <-
  svydesign(
    id = ~ VEREP, ### PSU - first stage cluster sampling
    strata = ~ VESTR, ### strata
    weights = ~ ANALWT_C, ### probability of being sampled
    data = nsduh2018, ### data
    nest = TRUE ### clusters nested within strata
  )

# Calculate prevalence of self-reported past-year marijuana, cocaine, and meth use by sex, age, race, and education
# Need to capture both "yes" and "no" use responses

# Totals by sex, age, race, and education
nsduh.pop <- 
  as.data.frame(svytable(~SEX + AGE + EDUC + RACE, dnsduh)) %>%
  mutate(DEMFLAG2 = paste(SEX, AGE, EDUC, RACE, sep = "_")) %>%
  arrange(SEX, AGE, EDUC, RACE)

# Totals by past-year drug use, sex, age, race, and education
drug.counts <- function(formula, x){
  formula <- make.formula(formula)
  df <- as.data.frame(svytable(formula, dnsduh)) %>%
    mutate(DEMFLAG2 = paste(SEX, AGE, EDUC, RACE, sep = "_")) %>%
    dplyr::select(DEMFLAG2, !!x, Freq) %>%
    dplyr::arrange(DEMFLAG2, !!x)
  return(df)
}

mj.pop <- drug.counts("MRJMON + SEX + AGE + EDUC + RACE", quo(MRJMON))
cocaine.pop <- drug.counts("COCMON + SEX + AGE + EDUC + RACE", quo(COCMON))
meth.pop <- drug.counts("METHAMMON + SEX + AGE + EDUC + RACE", quo(METHAMMON))

# Proportions of past-month drug use by sex, age, race, and education
drug.prev1 <- 
  left_join(mj.pop, nsduh.pop, by = "DEMFLAG2") %>%
  mutate(MJPREV = ifelse(is.na(Freq.x/Freq.y), 0, Freq.x/Freq.y)) %>%
  rename(MJPOP = Freq.x,
         TOTALPOP = Freq.y) %>%
  select(SEX, AGE, EDUC, RACE, DEMFLAG2, TOTALPOP, MRJMON, MJPOP, MJPREV)

drug.prev2 <- 
  left_join(cocaine.pop, nsduh.pop, by = "DEMFLAG2") %>%
  mutate(COCPREV = ifelse(is.na(Freq.x/Freq.y), 0, Freq.x/Freq.y)) %>%
  rename(COCPOP = Freq.x,
         TOTALPOP = Freq.y,
         DEMFLAG.COC = DEMFLAG2) %>%
  select(DEMFLAG.COC, COCMON, COCPOP, COCPREV)

drug.prev3 <- 
  left_join(meth.pop, nsduh.pop, by = "DEMFLAG2") %>%
  mutate(METHPREV = ifelse(is.na(Freq.x/Freq.y), 0, Freq.x/Freq.y)) %>%
  rename(METHPOP = Freq.x,
         TOTALPOP = Freq.y, 
         DEMFLAG.METH = DEMFLAG2) %>%
  select(DEMFLAG.METH, METHAMMON, METHPOP, METHPREV)

drug.prev <-
  cbind(drug.prev1, drug.prev2, drug.prev3)

###################################################################################
# 3 - Create combined dataset of demographic and self-reported drug use variables #
###################################################################################

# Merge drug use prevalence with total counts by demographic characteristics from Census data
# Generate population sizes by self-reported past-year drug use
census.drug <- 
  full_join(drug.prev, census.data, by = "DEMFLAG2") %>%
  rename(SEX = SEX.y, AGE = AGE.y, EDUC = EDUC.y, RACE = RACE.y) %>%
  select(SEX, SEX_FACT, AGE, AGE_FACT, EDUC, EDUC_FACT, RACE, RACE_FACT,
         DEMFLAG, DEMFLAG2, MRJMON, MJPREV, COCMON, COCPREV, METHAMMON, METHPREV, value) %>%
  arrange(DEMFLAG) %>%
  mutate(MJ.SELFREPORT = ifelse(is.na(MRJMON), value, round(value*MJPREV)),
         COC.SELFREPORT = ifelse(is.na(COCMON), value, round(value*COCPREV)),
         METH.SELFREPORT = ifelse(is.na(METHAMMON), value, round(value*METHPREV)))

####################################################################
# 4 - Generate true prevalence of drug use based on Fendrich paper #
####################################################################

# Add values of underreporting percentages per Fendrich papers by demographics
census.drug <-
  census.drug %>%
  mutate(
    AGE.BIAS.MJ = ifelse(SEX==1 & AGE==2, 0.0556,
                         ifelse(SEX==1 & AGE==3, 0.0449,
                                ifelse(SEX==1 & AGE>3, 0,
                                       ifelse(SEX==2 & AGE==2, 0.1180,
                                              ifelse(SEX==2 & AGE==3, 0.0965, 0))))),
                         
    AGE.BIAS.COC = ifelse(SEX==1 & AGE==2, 0.0222,
                          ifelse(SEX==1 & AGE==3, 0.0778,
                                 ifelse(SEX==1 & AGE>3, 0.0345,
                                        ifelse(SEX==2 & AGE==2, 0.0348,
                                               ifelse(SEX==2 & AGE==3, 0.1181, 0.0537))))),
                          
    RACE.BIAS.MJ = ifelse(SEX==1 & RACE==1, 0.0049,
                          ifelse(SEX==1 & RACE>1 & RACE!=7, 0.0667,
                                 ifelse(SEX==1 & RACE==7, 0.0316,
                                        ifelse(SEX==2 & RACE==1, 0.0111,
                                               ifelse(SEX==2 & RACE==7, 0.0691, 0.1397))))),
    
    RACE.BIAS.COC = ifelse(SEX==1 & RACE==1, 0.0141,
                           ifelse(SEX==1 & RACE>1 & RACE!=7, 0.0541,
                                  ifelse(SEX==1 & RACE==7, 0.0105,
                                         ifelse(SEX==2 & RACE==1, 0.0222,
                                                ifelse(SEX==2 & RACE==7, 0.0165, 0.0832))))),
    
    EDUC.BIAS.MJ = ifelse(SEX==1 & EDUC<=2, 0.1059,
                          ifelse(SEX==1 & EDUC>2, 0,
                                 ifelse(SEX==2 & EDUC<=2, 0.2121, 0))),
    
    EDUC.BIAS.COC = ifelse(SEX==1 & EDUC<=2, 0.1047,
                           ifelse(SEX==1 & EDUC==3, 0.0204,
                                  ifelse(SEX==1 & EDUC>3, 0,
                                         ifelse(SEX==2 & EDUC<=2, 0.1566,
                                                ifelse(SEX==2 & EDUC==3, 0.0320, 0))))))

# For each demographic grouping select the mean amount of bias
census.drug$MJ.BIAS <- ifelse(census.drug$AGE==0, 0, apply(census.drug[,c("AGE.BIAS.MJ","RACE.BIAS.MJ","EDUC.BIAS.MJ")], 1, mean))
census.drug$COC.BIAS <- ifelse(census.drug$AGE==0, 0, apply(census.drug[,c("AGE.BIAS.COC","RACE.BIAS.COC","EDUC.BIAS.COC")], 1, mean))

# Create true drug use variables
census.drug$MJ.TRUE <- rep(NA, nrow(census.drug))
census.drug$COC.TRUE <- rep(NA, nrow(census.drug))
census.drug$METH.TRUE <- rep(NA, nrow(census.drug))

for (i in 1:nrow(census.drug)){
  census.drug$MJ.TRUE[i] <-
    ifelse(census.drug$MRJMON[i]==1, round(census.drug$MJ.SELFREPORT[i]/(1-census.drug$MJ.BIAS[i])), 
           ifelse(census.drug$MRJMON[i]==0, census.drug$MJ.SELFREPORT[i] - round(census.drug$MJ.BIAS[i]*round(census.drug$MJ.SELFREPORT[i+1]/(1-census.drug$MJ.BIAS[i]))), NA))

  census.drug$COC.TRUE[i] <-
    ifelse(census.drug$COCMON[i]==1, round(census.drug$COC.SELFREPORT[i]/(1-census.drug$COC.BIAS[i])), 
           ifelse(census.drug$COCMON[i]==0, census.drug$COC.SELFREPORT[i] - round(census.drug$COC.BIAS[i]*round(census.drug$COC.SELFREPORT[i+1]/(1-census.drug$COC.BIAS[i]))), NA))
}

###################################################
# 5 - Generate individual-level simulated dataset #
###################################################

# Create individual observations based on TRUE drug use
true.data.mj <- 
  uncount(census.drug, MJ.TRUE) %>%
  rowid_to_column() %>%
  mutate(MJ.TRUE = ifelse(MRJMON == 0, 0, 1)) %>%
  select(rowid, SEX, SEX_FACT, AGE, AGE_FACT, EDUC, EDUC_FACT, RACE, RACE_FACT, DEMFLAG, DEMFLAG2, MJ.TRUE)

true.data.coc <- 
  uncount(census.drug, COC.TRUE) %>%
  rowid_to_column() %>%
  mutate(COC.TRUE = ifelse(COCMON == 0, 0, 1)) %>%
  select(rowid, COC.TRUE)

# Create individual observations based on BIASED drug use
bias.data.mj <- 
  uncount(census.drug, MJ.SELFREPORT) %>%
  rowid_to_column() %>%
  mutate(MJ.SELFREPORT = ifelse(MRJMON == 0, 0, 1)) %>%
  select(rowid, MJ.SELFREPORT)

bias.data.coc <- 
  uncount(census.drug, COC.SELFREPORT) %>%
  rowid_to_column() %>%
  mutate(COC.SELFREPORT = ifelse(COCMON == 0, 0, 1)) %>%
  select(rowid, COC.SELFREPORT)


# Join dataframes and limit to n = 100,000
sim.data.final <- 
  left_join(true.data.mj, true.data.coc, by = "rowid") %>%
  left_join(., bias.data.mj, by = "rowid") %>%
  left_join(., bias.data.coc, by = "rowid") %>%
  sample_n(100000)

# Save files
write.csv(sim.data.final, file = "simdatafinal.csv", row.names = FALSE)
save(sim.data.final, file = "/PATH/simdatafinal.RData")
save(census.data, file = "/PATH/censusdata.RData")
save(census.drug, file = "/PATH/censusdrug.RData")
save(nsduh.pop, file = "/PATH/nsduhpop.RData")
save(mj.pop, file = "/PATH/mj.pop.RData")
save(cocaine.pop, file = "/PATH/cocaine.pop.RData")

#######################################
# 6 - Examples of correcting the bias #
#######################################

load("DrugUseBias.RData")

##############################################################
# COMPARISON OF ORIGINAL DATA AND TOTAL SIMULATED POPULATION #
##############################################################

# Calculate demographic proportions in Census and simulated population
## Total
census.total <- raw.data %>%
  summarise(TOTAL = sum(HISP, WHITE, BLACK, AI.AN, ASIAN, NH.PA, OTHER, MULT))

## Sex
census.sex <- raw.data %>%
  group_by(SEX_FACT) %>%
  summarise(TOTAL = sum(HISP, WHITE, BLACK, AI.AN, ASIAN, NH.PA, OTHER, MULT))
census.sex$CENSUS <- census.sex$TOTAL/census.total$TOTAL
sim.sex <- data.frame(prop.table(table(sim.data.final$SEX_FACT))) %>% rename(SEX_FACT = Var1, SIMULATED = Freq)
sex.prop <- left_join(census.sex, sim.sex, by = "SEX_FACT") %>% rename(VAR = SEX_FACT)

## Age
census.age <- raw.data %>%
  group_by(AGE_FACT) %>%
  summarise(TOTAL = sum(HISP, WHITE, BLACK, AI.AN, ASIAN, NH.PA, OTHER, MULT))
class(census.age$AGE_FACT) <- "factor"
census.age$CENSUS <- census.age$TOTAL/census.total$TOTAL
sim.age <- data.frame(prop.table(table(sim.data.final$AGE_FACT))) %>% rename(AGE_FACT = Var1, SIMULATED = Freq)
age.prop <- left_join(census.age, sim.age, by = "AGE_FACT") %>% rename(VAR = AGE_FACT)

## Education
census.educ <- raw.data %>%
  group_by(EDUC_FACT) %>%
  summarise(TOTAL = sum(HISP, WHITE, BLACK, AI.AN, ASIAN, NH.PA, OTHER, MULT))
class(census.educ$EDUC_FACT) <- "factor"
census.educ$CENSUS <- census.educ$TOTAL/census.total$TOTAL
sim.educ <- data.frame(prop.table(table(sim.data.final$EDUC_FACT))) %>% rename(EDUC_FACT = Var1, SIMULATED = Freq)
educ.prop <- left_join(census.educ, sim.educ, by = "EDUC_FACT") %>% rename(VAR = EDUC_FACT)

## Race
census.race <- raw.data %>%
  summarise(HISP = sum(HISP)/census.total$TOTAL,
            WHITE = sum(WHITE)/census.total$TOTAL,
            BLACK = sum(BLACK)/census.total$TOTAL,
            AI.AN = sum(AI.AN)/census.total$TOTAL,
            ASIAN = sum(ASIAN)/census.total$TOTAL,
            NH.PA = sum(NH.PA)/census.total$TOTAL,
            OTHER = sum(OTHER)/census.total$TOTAL,
            MULT = sum(MULT)/census.total$TOTAL) %>%
  pivot_longer(cols = c(HISP, WHITE, BLACK, AI.AN, ASIAN, NH.PA, OTHER, MULT),
               names_to = "RACE_FACT") %>%
  mutate(RACE_FACT = factor(RACE_FACT)) %>%
  rename(CENSUS = value)

sim.race <- data.frame(prop.table(table(sim.data.final$RACE_FACT))) %>%
  rename(RACE_FACT = Var1,
         SIMULATED = Freq)

race.prop <- left_join(census.race, sim.race, by = "RACE_FACT") %>% rename(VAR = RACE_FACT) %>% mutate(TOTAL = NA)

table1_demo <- rbind(sex.prop, age.prop, educ.prop, race.prop)
write.csv(table1_demo, file = "Table 1 Demographics.csv", row.names = FALSE)

## Calculate Self-reported Drug Use in NSDUH and Simulated Population
mj.n <- data.frame(prop.table(svytable(~MRJMON, dnsduh))) %>% rename(MJ = MRJMON, NSDUH = Freq)
mj.sim <- data.frame(prop.table(table(sim.data.final$MJ.SELFREPORT))) %>% rename(MJ = Var1, SIMULATED = Freq)
coc.n <- data.frame(prop.table(svytable(~COCMON, dnsduh))) %>% rename(COC = COCMON, NSDUH = Freq)
coc.sim <- data.frame(prop.table(table(sim.data.final$COC.SELFREPORT))) %>% rename(COC = Var1, SIMULATED = Freq)

table1_drugs <- cbind(mj.n, mj.sim, coc.n, coc.sim)
write.csv(table1_drugs, file = "Table 1 Drug Prevalence.csv", row.names = FALSE)

###########################################
# RESULTS FROM TOTAL SIMULATED POPULATION #
###########################################

# Marijuana
mjsr <- 
  as.data.frame(prop.table(table(sim.data.final$MJ.SELFREPORT))) %>%
  rename(SelfReport = Freq,
         Drug = Var1)
mjtrue <- 
  as.data.frame(prop.table(table(sim.data.final$MJ.TRUE))) %>%
  rename(Truth = Freq,
         Drug = Var1)

mjresults <- 
  left_join(mjsr, mjtrue, by="Drug") %>%
  mutate(Drug = ifelse(Drug==1,"MJ",""))

mjresults <- mjresults[2,]

# Cocaine
cocsr <- 
  as.data.frame(prop.table(table(sim.data.final$COC.SELFREPORT))) %>%
  rename(SelfReport = Freq,
         Drug = Var1)
coctrue <- 
  as.data.frame(prop.table(table(sim.data.final$COC.TRUE))) %>%
  rename(Truth = Freq,
         Drug = Var1)

cocresults <- 
  left_join(cocsr, coctrue, by="Drug") %>%
  mutate(Drug = ifelse(Drug==1,"Cocaine",""))

cocresults <- cocresults[2,]

drugresults <- rbind(mjresults, cocresults) %>%
  mutate(Sens = SelfReport/Truth)

write.csv(drugresults, file = "Simulated Results Overall.csv", row.names = FALSE)

totalresults <- pivot_longer(drugresults,
                             c(SelfReport, Truth),
                             names_to = "Type") %>%
  mutate(Drug = factor(Drug, levels = c("MJ", "Cocaine")),
         Type = factor(Type, levels = c("SelfReport", "Truth")))

true.props <- 
  filter(totalresults, Type=="Truth") %>%
  select(-Sens) %>%
  rename(Prev = value)

(totalresults.plot <- 
    ggplot(data = totalresults,
           aes(x = Drug, y = value, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("#453781FF","#B8DE29FF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

(totalresults.plot.noleg <- 
    ggplot(data = totalresults,
           aes(x = Drug, y = value, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), show.legend = FALSE) +
    scale_fill_manual(values = c("#453781FF","#B8DE29FF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

save_plot("Total Sim Pop Results.pdf", totalresults.plot)
save_plot("Total Sim Pop Results - No Legend.pdf", totalresults.plot.noleg)

########################
# SIMPLE RANDOM SAMPLE #
########################

##############################
# Draw single SRS as example #
##############################
set.seed(060520)
samp <- data.frame(sample(sim.data.final$rowid, 1000, replace = FALSE))
samp$rowid <- samp[,1]
samp <- select(samp, rowid)
sim.study <- left_join(samp, sim.data.final, by = "rowid")

sim.study.totals <-
  sim.study %>%
  group_by(DEMFLAG) %>%
  summarise(TOTAL = n())

selfreport.mj <-
  sim.study %>%
  filter(MJ.SELFREPORT == 1) %>%
  group_by(DEMFLAG) %>%
  summarise(MJ.SELFREPORT = n())

selfreport.coc <-
  sim.study %>%
  filter(COC.SELFREPORT == 1) %>%
  group_by(DEMFLAG) %>%
  summarise(COC.SELFREPORT = n())

sim.selfreport <- left_join(sim.study.totals, selfreport.mj, by = "DEMFLAG")
sim.selfreport <- left_join(sim.selfreport, selfreport.coc, by = "DEMFLAG")

# True counts of drug use by sex, age, education, and race
true.mj <-
  sim.study %>%
  filter(MJ.TRUE == 1) %>%
  group_by(DEMFLAG) %>%
  summarise(MJ.TRUE = n())

true.coc <-
  sim.study %>%
  filter(COC.TRUE == 1) %>%
  group_by(DEMFLAG) %>%
  summarise(COC.TRUE = n())


sim.true <- left_join(sim.study.totals, true.mj, by = "DEMFLAG")
sim.true <- left_join(sim.true, true.coc, by = "DEMFLAG")


# Corrected counts of drug use by sex, age, education, and race
# Replace N/As (true zeros) with 0.1 to avoid undefined sensitivities
sim.corrected <- left_join(sim.selfreport, sim.true, by = c("DEMFLAG","TOTAL"))
sim.corrected <- mutate_if(sim.corrected, is.integer, as.numeric)
sim.corrected[is.na(sim.corrected)] <- 0.1

# Calculate sensitivity of self-report by sex, age, education, and race
sim.corrected <-
  sim.corrected %>%
  mutate(MJ.SENS = MJ.SELFREPORT/MJ.TRUE,
         COC.SENS = COC.SELFREPORT/COC.TRUE)

sim.corrected <-
  sim.corrected %>%
  mutate(MJ.CORRECTED = MJ.SELFREPORT/MJ.SENS,
         COC.CORRECTED = COC.SELFREPORT/COC.SENS)

sim.corrected[sim.corrected == 0.1] <- 0

corrected.data <- sim.corrected[ , c(-1, -5:-8)]
corrected.counts <- summarise_all(corrected.data, sum)

corrected.samp <- 
  corrected.counts %>%
  pivot_longer(MJ.SELFREPORT:COC.CORRECTED) %>%
  mutate(Drug = factor(rep(c("MJ","Cocaine"), 2), levels = c("MJ","Cocaine")),
         Type = factor(rep(c("SelfReport","Corrected"), each = 2), 
                       levels = c("SelfReport","Corrected")),
         Prev = value/TOTAL) %>%
  select(Drug, Type, Prev)

corrected.samp.final <- bind_rows(corrected.samp, true.props)

write.csv(corrected.samp.final, "First Simple Random Sample Results.csv", row.names=F)

(corrected.samp.plot <- 
    ggplot(data = corrected.samp.final,
           aes(x = Drug, y = Prev, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("#453781FF","#1F968BFF","#B8DE29FF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

save_plot(filename = "Simulated SRS Correction.pdf", plot = corrected.samp.plot)

#########################################################################################################
# Draw multiple simple random samples to estimate mean and IQR of self-reported and corrected estimates #
#########################################################################################################
corrected.res <- list()
set.seed(031120)
for(i in 1:500){
  # Create "study sample" from full simulated dataset as a simple random sample
  samp <- sample(sim.data.final$rowid, 1000, replace = FALSE)
  
  sim.study <- sim.data.final %>%
    filter(rowid %in% samp)
  
  # Self-reported counts and prevalence of drug use by sex, age, education, and race
  sim.study.totals <-
    sim.study %>%
    group_by(DEMFLAG) %>%
    summarise(TOTAL = n())
  
  selfreport.mj <-
    sim.study %>%
    filter(MJ.SELFREPORT == 1) %>%
    group_by(DEMFLAG) %>%
    summarise(MJ.SELFREPORT = n())
  
  selfreport.coc <-
    sim.study %>%
    filter(COC.SELFREPORT == 1) %>%
    group_by(DEMFLAG) %>%
    summarise(COC.SELFREPORT = n())

  
  sim.selfreport.mult <- left_join(sim.study.totals, selfreport.mj, by = "DEMFLAG")
  sim.selfreport.mult <- left_join(sim.selfreport.mult, selfreport.coc, by = "DEMFLAG")
  
  # True counts of drug use by sex, age, education, and race
  true.mj <-
    sim.study %>%
    filter(MJ.TRUE == 1) %>%
    group_by(DEMFLAG) %>%
    summarise(MJ.TRUE = n())
  
  true.coc <-
    sim.study %>%
    filter(COC.TRUE == 1) %>%
    group_by(DEMFLAG) %>%
    summarise(COC.TRUE = n())
  
  sim.true.mult <- left_join(sim.study.totals, true.mj, by = "DEMFLAG")
  sim.true.mult <- left_join(sim.true.mult, true.coc, by = "DEMFLAG")
  
  # Corrected counts of drug use by sex, age, education, and race
  # Replace N/As (true zeros) with 0.1 to avoid undefined sensitivities
  sim.corrected.mult <- left_join(sim.selfreport.mult, sim.true.mult, by = c("DEMFLAG","TOTAL"))
  sim.corrected.mult <- mutate_if(sim.corrected.mult, is.integer, as.numeric)
  sim.corrected.mult[is.na(sim.corrected.mult)] <- 0.1
  
  # Calculate sensitivity of self-report by sex, age, education, and race
  sim.corrected.mult <-
    sim.corrected.mult %>%
    mutate(MJ.SENS = MJ.SELFREPORT/MJ.TRUE,
           COC.SENS = COC.SELFREPORT/COC.TRUE)
  
  sim.corrected.mult <-
    sim.corrected.mult %>%
    mutate(MJ.CORRECTED = MJ.SELFREPORT/MJ.SENS,
           COC.CORRECTED = COC.SELFREPORT/COC.SENS)
  
  sim.corrected.mult[sim.corrected.mult == 0.1] <- 0
  
  corrected.data.mult <- sim.corrected.mult[ , c(-1, -5:-8)]
  corrected.counts.mult <- summarise_all(corrected.data.mult, sum)
  
  corrected.res[[i]] <- corrected.counts.mult 
}

corrected.results <- bind_rows(corrected.res) %>%
  mutate(across(-TOTAL, ~ ./ TOTAL)) %>%
  mutate(MJ.TRUE = true.props$Prev[1],
         COC.TRUE = true.props$Prev[2],
         MJ.SELFREPORT.BIAS = MJ.SELFREPORT - MJ.TRUE,
         MJ.CORRECTED.BIAS = MJ.CORRECTED - MJ.TRUE,
         COC.SELFREPORT.BIAS = COC.SELFREPORT - COC.TRUE,
         COC.CORRECTED.BIAS = COC.CORRECTED - COC.TRUE)

corrected.results.sum <-
  corrected.results %>%
  summarise(MJ.SELFREPORT.MEANBIAS = mean(MJ.SELFREPORT.BIAS),
            MJ.SELFREPORT.SEBIAS = sd(MJ.SELFREPORT.BIAS) / sqrt(500),
            MJ.SELFREPORT.MEAN = mean(MJ.SELFREPORT),
            MJ.SELFREPORT.SE = sd(MJ.SELFREPORT) / sqrt(500),
            MJ.SELFREPORT.LCL = MJ.SELFREPORT.MEAN - 1.96*MJ.SELFREPORT.SE,
            MJ.SELFREPORT.UCL = MJ.SELFREPORT.MEAN + 1.96*MJ.SELFREPORT.SE,
            
            MJ.CORRECTED.MEANBIAS = mean(MJ.CORRECTED.BIAS),
            MJ.CORRECTED.SEBIAS = sd(MJ.CORRECTED.BIAS) / sqrt(500),
            MJ.CORRECTED.MEAN = mean(MJ.CORRECTED),
            MJ.CORRECTED.SE = sd(MJ.CORRECTED) / sqrt(500),
            MJ.CORRECTED.LCL = MJ.CORRECTED.MEAN - 1.96*MJ.CORRECTED.SE,
            MJ.CORRECTED.UCL = MJ.CORRECTED.MEAN + 1.96*MJ.CORRECTED.SE,
            
            COC.SELFREPORT.MEANBIAS = mean(COC.SELFREPORT.BIAS),
            COC.SELFREPORT.SEBIAS = sd(COC.SELFREPORT.BIAS) / sqrt(500),
            COC.SELFREPORT.MEAN = mean(COC.SELFREPORT),
            COC.SELFREPORT.SE = sd(COC.SELFREPORT) / sqrt(500),
            COC.SELFREPORT.LCL = COC.SELFREPORT.MEAN - 1.96*COC.SELFREPORT.SE,
            COC.SELFREPORT.UCL = COC.SELFREPORT.MEAN + 1.96*COC.SELFREPORT.SE,
            
            COC.CORRECTED.MEANBIAS = mean(COC.CORRECTED.BIAS),
            COC.CORRECTED.SEBIAS = sd(COC.CORRECTED.BIAS) / sqrt(500),
            COC.CORRECTED.MEAN = mean(COC.CORRECTED),
            COC.CORRECTED.SE = sd(COC.CORRECTED) / sqrt(500),
            COC.CORRECTED.LCL = COC.CORRECTED.MEAN - 1.96*COC.CORRECTED.SE,
            COC.CORRECTED.UCL = COC.CORRECTED.MEAN + 1.96*COC.CORRECTED.SE) %>%
  
  mutate(MJ.SELFREPORT.MSE = (MJ.SELFREPORT.SE^2) - (MJ.SELFREPORT.MEAN^2),
         MJ.CORRECTED.MSE = (MJ.CORRECTED.SE^2) - (MJ.CORRECTED.MEAN^2),
         COC.SELFREPORT.MSE = (COC.SELFREPORT.SE^2) - (COC.SELFREPORT.MEAN^2),
         COC.CORRECTED.MSE = (COC.CORRECTED.SE^2) - (COC.CORRECTED.MEAN^2),
         MJ.TRUE.MEAN = true.props$Prev[1],
         COC.TRUE.MEAN = true.props$Prev[2])

corrected.results.final <- corrected.results.sum %>%
  pivot_longer(everything()) %>%
  separate(col = name,
           into = c("Drug", "Type", "Result")) %>%
  pivot_wider(id_cols = c(Drug, Type), names_from = Result, values_from = value) %>%
  mutate(Drug = factor(Drug, levels = c("MJ", "COC")),
         Type = factor(Type, levels = c("SELFREPORT", "CORRECTED", "TRUE")))

write.csv(corrected.results.final, file = "Simple Random Sample Results.csv", row.names = FALSE)

(correctedresults.plot <- 
    ggplot(data = corrected.results.final,
           aes(x = Drug, y = MEAN, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), position = position_dodge2(width = 0.5, padding = 0.5), size=0.3) +
    scale_fill_manual(values = c("#453781FF","#1F968BFF","#B8DE29FF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

(correctedresults.noleg <- 
    ggplot(data = corrected.results.final,
           aes(x = Drug, y = MEAN, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), show.legend = FALSE) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), position = position_dodge2(width = 0.5, padding = 0.5), size=0.3) +
    scale_fill_manual(values = c("#453781FF","#1F968BFF","#B8DE29FF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

save_plot(filename = "Simulated SRS Correction with SE.pdf", plot = correctedresults.plot)


############################
# STRATIFIED RANDOM SAMPLE #
############################

# Create "study sample" from full simulated dataset as a stratified random sample
set.seed(062020)

bydemflag <- 
  sim.data.final %>% 
  group_by(DEMFLAG) %>%
  select(rowid, DEMFLAG) %>%
  add_count(DEMFLAG)

raw1 <- filter(bydemflag, n < 100) %>% group_by(DEMFLAG)
raw2 <- filter(bydemflag, n >= 100) %>% group_by(DEMFLAG)

############################################
# Draw single stratified sample as example #
############################################
samp1 <- slice_sample(raw1, n = 5, replace = FALSE)
samp2 <- slice_sample(raw2, prop = 0.01, replace = FALSE)
samp3 <- 
  bind_rows(samp1, samp2) %>% 
  group_by(DEMFLAG) %>%
  arrange (DEMFLAG, rowid)

samp3 <- 
  add_count(samp3, DEMFLAG, name = "TOTAL") %>%
  mutate(selection.prob = TOTAL/n)

sim.study.stratified <- 
  left_join(samp3, sim.data.final, by = c("rowid","DEMFLAG")) %>%
  select(-n, -TOTAL)

# Self-reported counts and prevalence of drug use by sex, age, education, and race
sim.selection.prob <-
  sim.study.stratified %>%
  select(DEMFLAG, selection.prob) %>%
  distinct(.)

sim.study.totals <-
  sim.study.stratified %>%
  group_by(DEMFLAG) %>%
  summarise(TOTAL = n())

selfreport.mj <-
  sim.study.stratified %>%
  filter(MJ.SELFREPORT == 1) %>%
  group_by(DEMFLAG) %>%
  summarise(MJ.SELFREPORT = n())

selfreport.coc <-
  sim.study.stratified %>%
  filter(COC.SELFREPORT == 1) %>%
  group_by(DEMFLAG) %>%
  summarise(COC.SELFREPORT = n())

sim.selfreport.stratified <- left_join(sim.study.totals, selfreport.mj, by = "DEMFLAG")
sim.selfreport.stratified <- left_join(sim.selfreport.stratified, selfreport.coc, by = "DEMFLAG")

# True counts of drug use by sex, age, education, and race
true.mj <-
  sim.study.stratified %>%
  filter(MJ.TRUE == 1) %>%
  group_by(DEMFLAG) %>%
  summarise(MJ.TRUE = n())

true.coc <-
  sim.study.stratified %>%
  filter(COC.TRUE == 1) %>%
  group_by(DEMFLAG) %>%
  summarise(COC.TRUE = n())

sim.true.stratified <- left_join(sim.study.totals, true.mj, by = "DEMFLAG")
sim.true.stratified <- left_join(sim.true.stratified, true.coc, by = "DEMFLAG")

# Corrected counts of drug use by sex, age, education, and race
# Replace N/As (true zeros) with 0.1 to avoid undefined sensitivities
sim.corrected.stratified <- left_join(sim.selfreport.stratified, sim.true.stratified, by = c("DEMFLAG","TOTAL"))
sim.corrected.stratified <- mutate_if(sim.corrected.stratified, is.integer, as.numeric)
sim.corrected.stratified[is.na(sim.corrected.stratified)] <- 0.1

# Calculate sensitivity of self-report by sex, age, education, and race
sim.corrected.stratified <-
  sim.corrected.stratified %>%
  mutate(MJ.SENS = MJ.SELFREPORT/MJ.TRUE,
         COC.SENS = COC.SELFREPORT/COC.TRUE)

sim.corrected.stratified <-
  sim.corrected.stratified %>%
  mutate(MJ.CORRECTED = MJ.SELFREPORT/MJ.SENS,
         COC.CORRECTED = COC.SELFREPORT/COC.SENS)

sim.corrected.stratified[sim.corrected.stratified == 0.1] <- 0

# Multiply all values by inverse of selection probability
corrected.data <- sim.corrected.stratified[ , c(-1, -5:-8)] * (1/sim.selection.prob$selection.prob)

corrected.res <- summarise_all(corrected.data, sum)

corrected.results.strat <-
  corrected.res %>%
  pivot_longer(MJ.SELFREPORT:COC.CORRECTED) %>%
  mutate(Drug = factor(rep(c("MJ","Cocaine"), 2), levels = c("MJ","Cocaine")),
         Type = factor(rep(c("SelfReport","Corrected"), each = 2), 
                       levels = c("SelfReport","Corrected")),
         Prev = value/TOTAL) %>%
  select(Drug, Type, Prev)

corrected.strat.final <- bind_rows(corrected.results.strat, true.props)

write.csv(corrected.strat.final, "First Stratified Random Sample Results.csv", row.names = F)

######################################################################################################
# Draw multiple stratified samples to estimate mean and IQR of self-reported and corrected estimates #
######################################################################################################
corrected.res.strat <- list()
set.seed(031120)
for(i in 1:500){
  samp1 <- slice_sample(raw1, n = 5, replace = FALSE)
  samp2 <- slice_sample(raw2, prop = 0.01, replace = FALSE)
  samp3 <- 
    bind_rows(samp1, samp2) %>% 
    group_by(DEMFLAG) %>%
    arrange (DEMFLAG, rowid)
  
  samp3 <- 
    add_count(samp3, DEMFLAG, name = "TOTAL") %>%
    mutate(selection.prob = TOTAL/n)
  
  sim.study.stratified <- 
    left_join(samp3, sim.data.final, by = c("rowid","DEMFLAG")) %>%
    select(-n, -TOTAL)
  
  # Self-reported counts and prevalence of drug use by sex, age, education, and race
  sim.selection.prob <-
    sim.study.stratified %>%
    select(DEMFLAG, selection.prob) %>%
    distinct(.)
  
  sim.study.totals <-
    sim.study.stratified %>%
    group_by(DEMFLAG) %>%
    summarise(TOTAL = n())
  
  selfreport.mj <-
    sim.study.stratified %>%
    filter(MJ.SELFREPORT == 1) %>%
    group_by(DEMFLAG) %>%
    summarise(MJ.SELFREPORT = n())
  
  selfreport.coc <-
    sim.study.stratified %>%
    filter(COC.SELFREPORT == 1) %>%
    group_by(DEMFLAG) %>%
    summarise(COC.SELFREPORT = n())
  
  sim.selfreport.stratified.mult <- left_join(sim.study.totals, selfreport.mj, by = "DEMFLAG")
  sim.selfreport.stratified.mult <- left_join(sim.selfreport.stratified.mult, selfreport.coc, by = "DEMFLAG")
  
  # True counts of drug use by sex, age, education, and race
  true.mj <-
    sim.study.stratified %>%
    filter(MJ.TRUE == 1) %>%
    group_by(DEMFLAG) %>%
    summarise(MJ.TRUE = n())
  
  true.coc <-
    sim.study.stratified %>%
    filter(COC.TRUE == 1) %>%
    group_by(DEMFLAG) %>%
    summarise(COC.TRUE = n())
  
  sim.true.stratified.mult <- left_join(sim.study.totals, true.mj, by = "DEMFLAG")
  sim.true.stratified.mult <- left_join(sim.true.stratified.mult, true.coc, by = "DEMFLAG")
  
  # Corrected counts of drug use by sex, age, education, and race
  # Replace N/As (true zeros) with 0.1 to avoid undefined sensitivities
  sim.corrected.stratified.mult <- left_join(sim.selfreport.stratified.mult, sim.true.stratified.mult, by = c("DEMFLAG","TOTAL"))
  sim.corrected.stratified.mult <- mutate_if(sim.corrected.stratified.mult, is.integer, as.numeric)
  sim.corrected.stratified.mult[is.na(sim.corrected.stratified.mult)] <- 0.1
  
  # Calculate sensitivity of self-report by sex, age, education, and race
  sim.corrected.stratified.mult <-
    sim.corrected.stratified.mult %>%
    mutate(MJ.SENS = MJ.SELFREPORT/MJ.TRUE,
           COC.SENS = COC.SELFREPORT/COC.TRUE)
  
  sim.corrected.stratified.mult <-
    sim.corrected.stratified.mult %>%
    mutate(MJ.CORRECTED = MJ.SELFREPORT/MJ.SENS,
           COC.CORRECTED = COC.SELFREPORT/COC.SENS)
  
  sim.corrected.stratified.mult[sim.corrected.stratified.mult == 0.1] <- 0
  
  # Multiply all values by inverse of selection probability
  corrected.data.mult <- sim.corrected.stratified.mult[ , c(-1, -5:-8)] * (1/sim.selection.prob$selection.prob)
  
  corrected.res.mult <- summarise_all(corrected.data.mult, sum)
  
  corrected.res.strat[[i]] <- corrected.res.mult
}

corrected.results.strat.mult <- bind_rows(corrected.res.strat) %>%
  mutate(across(-TOTAL, ~ ./ TOTAL)) %>%
  mutate(MJ.TRUE = true.props$Prev[1],
         COC.TRUE = true.props$Prev[2],
         MJ.SELFREPORT.BIAS = MJ.SELFREPORT - MJ.TRUE,
         MJ.CORRECTED.BIAS = MJ.CORRECTED - MJ.TRUE,
         COC.SELFREPORT.BIAS = COC.SELFREPORT - COC.TRUE,
         COC.CORRECTED.BIAS = COC.CORRECTED - COC.TRUE)

corrected.results.strat.sum <-
  corrected.results.strat.mult %>%
  summarise(MJ.SELFREPORT.MEANBIAS = mean(MJ.SELFREPORT.BIAS),
            MJ.SELFREPORT.SEBIAS = sd(MJ.SELFREPORT.BIAS) / sqrt(500),
            MJ.SELFREPORT.MEAN = mean(MJ.SELFREPORT),
            MJ.SELFREPORT.SE = sd(MJ.SELFREPORT) / sqrt(500),
            MJ.SELFREPORT.LCL = MJ.SELFREPORT.MEAN - 1.96*MJ.SELFREPORT.SE,
            MJ.SELFREPORT.UCL = MJ.SELFREPORT.MEAN + 1.96*MJ.SELFREPORT.SE,
            
            MJ.CORRECTED.MEANBIAS = mean(MJ.CORRECTED.BIAS),
            MJ.CORRECTED.SEBIAS = sd(MJ.CORRECTED.BIAS) / sqrt(500),
            MJ.CORRECTED.MEAN = mean(MJ.CORRECTED),
            MJ.CORRECTED.SE = sd(MJ.CORRECTED) / sqrt(500),
            MJ.CORRECTED.LCL = MJ.CORRECTED.MEAN - 1.96*MJ.CORRECTED.SE,
            MJ.CORRECTED.UCL = MJ.CORRECTED.MEAN + 1.96*MJ.CORRECTED.SE,
            
            COC.SELFREPORT.MEANBIAS = mean(COC.SELFREPORT.BIAS),
            COC.SELFREPORT.SEBIAS = sd(COC.SELFREPORT.BIAS) / sqrt(500),
            COC.SELFREPORT.MEAN = mean(COC.SELFREPORT),
            COC.SELFREPORT.SE = sd(COC.SELFREPORT) / sqrt(500),
            COC.SELFREPORT.LCL = COC.SELFREPORT.MEAN - 1.96*COC.SELFREPORT.SE,
            COC.SELFREPORT.UCL = COC.SELFREPORT.MEAN + 1.96*COC.SELFREPORT.SE,
            
            COC.CORRECTED.MEANBIAS = mean(COC.CORRECTED.BIAS),
            COC.CORRECTED.SEBIAS = sd(COC.CORRECTED.BIAS) / sqrt(500),
            COC.CORRECTED.MEAN = mean(COC.CORRECTED),
            COC.CORRECTED.SE = sd(COC.CORRECTED) / sqrt(500),
            COC.CORRECTED.LCL = COC.CORRECTED.MEAN - 1.96*COC.CORRECTED.SE,
            COC.CORRECTED.UCL = COC.CORRECTED.MEAN + 1.96*COC.CORRECTED.SE) %>%
  
  mutate(MJ.SELFREPORT.MSE = (MJ.SELFREPORT.SE^2) - (MJ.SELFREPORT.MEAN^2),
         MJ.CORRECTED.MSE = (MJ.CORRECTED.SE^2) - (MJ.CORRECTED.MEAN^2),
         COC.SELFREPORT.MSE = (COC.SELFREPORT.SE^2) - (COC.SELFREPORT.MEAN^2),
         COC.CORRECTED.MSE = (COC.CORRECTED.SE^2) - (COC.CORRECTED.MEAN^2),
         MJ.TRUE.MEAN = true.props$Prev[1],
         COC.TRUE.MEAN = true.props$Prev[2])

corrected.results.strat.final <- corrected.results.strat.sum %>%
  pivot_longer(everything()) %>%
  separate(col = name,
           into = c("Drug", "Type", "Result")) %>%
  pivot_wider(id_cols = c(Drug, Type), names_from = Result, values_from = value) %>%
  mutate(Drug = factor(Drug, levels = c("MJ", "COC")),
         Type = factor(Type, levels = c("SELFREPORT", "CORRECTED", "TRUE")))

write.csv(corrected.results.strat.final, "Stratified Random Sample Results.csv")

(correctedresults.stratified.plot <- 
    ggplot(data = corrected.results.strat.final,
           aes(x = Drug, y = MEAN, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), position = position_dodge2(width = 0.5, padding = 0.5), size=0.3) +
    scale_fill_manual(values = c("#453781FF","#1F968BFF","#B8DE29FF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

(correctedresults.strat.noleg <- 
    ggplot(data = corrected.results.strat.final,
           aes(x = Drug, y = MEAN, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), show.legend = FALSE) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), position = position_dodge2(width = 0.5, padding = 0.5), size=0.3) +
    scale_fill_manual(values = c("#453781FF","#1F968BFF","#B8DE29FF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

save_plot(filename = "Stratified Sample Simulated Correction.pdf", plot = correctedresults.stratified.plot)

sim.corrected.plots <- plot_grid(correctedresults.plot, correctedresults.stratified.plot, align = "h")
save_plot(filename = "Sampled Simulation Correction.pdf", plot = sim.corrected.plots)

sim.corrected.plots.noleg <- plot_grid(correctedresults.noleg, correctedresults.strat.noleg, align = "h")
save_plot(filename = "Sampled Simulation Correction - No Legend.pdf", plot = sim.corrected.plots.noleg)

################################################
# POWER CURVE - REPEATED SIMPLE RANDOM SAMPLES #
################################################

# Create "study sample" from full simulated dataset simple random samples of different sizes

srs.func <- function(x){
  sampv <- sample(sim.data.final$rowid, x, replace = FALSE)
  sampdf <- data.frame(sampv)
  sampdf <- rename(sampdf, "rowid" = "sampv")
  df <- left_join(sampdf, sim.data.final, by = "rowid")

  estimates <-
    summarise(df, n = nrow(df), MJ.EST = sum(MJ.TRUE)/nrow(df), COC.EST = sum(COC.TRUE)/nrow(df))

  estimates$MJ.DIFF = abs(estimates$MJ.EST - as.numeric(true.props[1,3]))
  estimates$COC.DIFF = abs(estimates$COC.EST - as.numeric(true.props[2,3]))
  return(estimates)
}

sampsizes <- seq(100, 5000, by = 100)
power.analysis <- list()
power.data <- list()
set.seed(060520)

for(r in 1:500){
  for(i in seq_along(sampsizes)){
    power.analysis[[i]] <- srs.func(sampsizes[i])
  }
  power.data[[r]] <- bind_rows(power.analysis)
}

power.sim <-
  bind_rows(power.data) %>%
  group_by(n)

save(power.sim, file = "/PATH/power.sim.RData")

(mj.power.plot <-
    ggplot(data = power.sim,
           aes(x = n, y = MJ.DIFF)) +
    geom_point(colour = "#1F968BFF") +
    geom_smooth(colour = "black") +
    ylim(0, 0.15) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 300),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")) +
    scale_x_continuous(n.breaks = 11))

save_plot(filename = "MJ Power by Sample Size.pdf", plot = mj.power.plot)

(coc.power.plot <-
    ggplot(data = power.sim,
           aes(x = n, y = COC.DIFF)) +
    geom_point(colour = "#453781FF") +
    geom_smooth(colour = "black") +
    ylim(0, 0.15) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 300),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")) +
    scale_x_continuous(n.breaks = 11))

save_plot(filename = "Cocaine Power by Sample Size.pdf", plot = coc.power.plot)

power.plots <- plot_grid(mj.power.plot, coc.power.plot, align = "h", nrow = 1)
save_plot(filename = "Power Plots.pdf", plot = power.plots)

###############################
# 7 - Correct NSDUH 2018 Data #
###############################

# Get total counts by demographics and for each drug from NSDUH data
nsduh.pop2 <- 
  nsduh.pop %>%
  select(DEMFLAG2, Freq) %>%
  arrange(DEMFLAG2, Freq)

mj.pop2 <- filter(mj.pop, MRJMON == 1) %>%
  rename(MJ.SELFREPORT = Freq) %>%
  select(-MRJMON)

coc.pop2 <- filter(cocaine.pop, COCMON == 1) %>%
  rename(COC.SELFREPORT = Freq) %>%
  select(-COCMON)

# Sensitivities based on simple random sample
sensitivity <- 
  sim.corrected %>%
  select(DEMFLAG, MJ.SENS, COC.SENS)

nsduh.corrected <- left_join(nsduh.pop2, mj.pop2, by = "DEMFLAG2")
nsduh.corrected <- left_join(nsduh.corrected, coc.pop2, by = "DEMFLAG2")
nsduh.corrected <- left_join(nsduh.corrected, sensitivity, by = c("DEMFLAG2" = "DEMFLAG"))

# Sensitivities based on stratified random sample
sensitivity.strat <-
  sim.corrected.stratified %>%
  select(DEMFLAG, MJ.SENS, COC.SENS)

nsduh.corrected.str <- left_join(nsduh.pop2, mj.pop2, by = "DEMFLAG2")
nsduh.corrected.str <- left_join(nsduh.corrected.str, coc.pop2, by = "DEMFLAG2")
nsduh.corrected.str <- left_join(nsduh.corrected.str, sensitivity.strat, by = c("DEMFLAG2" = "DEMFLAG"))

# NOTE: In samples, we didn't sample every demographic combo in our study data so setting missing values to 1
# Effectively assuming perfect sensitivity since we don't know
# Also need to change zero sensitivity values to a small number to avoid undefined values
# Correct self-reported values by dividing by sensitivity

# SIMPLE RANDOM SAMPLE #
nsduh.corrected <-
  nsduh.corrected %>%
  mutate(MJ.SENS = if_else(is.na(MJ.SENS), 1, if_else(MJ.SENS == 0, 0.1, MJ.SENS)),
         COC.SENS = if_else(is.na(COC.SENS), 1, if_else(COC.SENS == 0, 0.1, COC.SENS)),
         MJ.CORRECTED = MJ.SELFREPORT/MJ.SENS,
         COC.CORRECTED = COC.SELFREPORT/COC.SENS)

nsduh.corrected2 <- nsduh.corrected[ , c(-1, -5:-6)]
nsduh.corrected2 <- summarise_all(nsduh.corrected2, sum)

nsduh.corrected2 <-
  nsduh.corrected2 %>%
  pivot_longer(MJ.SELFREPORT:COC.CORRECTED) %>%
  mutate(Drug = factor(rep(c("MJ","Cocaine"), 2), levels = c("MJ","Cocaine")),
         Type = factor(rep(c("SelfReport","Corrected"), each = 2), 
                       levels = c("SelfReport","Corrected")),
         Prev = value/Freq) %>%
  select(Drug, Type, Prev)

write.csv(nsduh.corrected2, "Corrected NSDUH - Simple Random Sample.csv")

(correctednsduh.plot <- 
    ggplot(data = nsduh.corrected2,
           aes(x = Drug, y = Prev, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    ylim(0, 0.20) +
    scale_fill_manual(values = c("#453781FF","#1F968BFF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

(correctednsduh.plot.noleg <- 
    ggplot(data = nsduh.corrected2,
           aes(x = Drug, y = Prev, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), show.legend = FALSE) +
    ylim(0, 0.20) +
    scale_fill_manual(values = c("#453781FF","#1F968BFF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

save_plot("Corrected NSDUH.pdf", correctednsduh.plot)

# STRATIFIED RANDOM SAMPLE #
# Change zero sensitivity values to a small number to avoid undefined values
# Correct self-report values by dividing by sensitivity
nsduh.corrected.str <-
  nsduh.corrected.str %>%
  mutate(MJ.SENS = if_else(is.na(MJ.SENS), 1, if_else(MJ.SENS == 0, 0.1, MJ.SENS)),
         COC.SENS = if_else(is.na(COC.SENS), 1, if_else(COC.SENS == 0, 0.1, COC.SENS)),
         MJ.CORRECTED = MJ.SELFREPORT/MJ.SENS,
         COC.CORRECTED = COC.SELFREPORT/COC.SENS)

nsduh.corrected.str2 <- nsduh.corrected.str[ , c(-1, -5:-6)]
nsduh.corrected.str2 <- summarise_all(nsduh.corrected.str2, sum)

nsduh.corrected.str2 <-
  nsduh.corrected.str2 %>%
  pivot_longer(MJ.SELFREPORT:COC.CORRECTED) %>%
  mutate(Drug = factor(rep(c("MJ","Cocaine"), 2), levels = c("MJ","Cocaine")),
         Type = factor(rep(c("SelfReport","Corrected"), each = 2), 
                       levels = c("SelfReport","Corrected")),
         Prev = value/Freq) %>%
  select(Drug, Type, Prev)

write.csv(nsduh.corrected.str2, "Corrected NSDUH - Stratified Random Sample.csv")

(correctednsduh.plot2 <- 
    ggplot(data = nsduh.corrected.str2,
           aes(x = Drug, y = Prev, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    ylim(0, 0.20) +
    scale_fill_manual(values = c("#453781FF","#1F968BFF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

(correctednsduh.plot2.noleg <- 
    ggplot(data = nsduh.corrected.str2,
           aes(x = Drug, y = Prev, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), show.legend = FALSE) +
    ylim(0, 0.20) +
    scale_fill_manual(values = c("#453781FF","#1F968BFF")) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_line(colour = "black")))

save_plot("Corrected NSDUH Stratified Sample.pdf", correctednsduh.plot2)

corrections.plots <- plot_grid(correctednsduh.plot.noleg, correctednsduh.plot2.noleg, align = "h")
save_plot("Corrected NSDUH Samples.pdf", corrections.plots)

#######################
# Save image of files #
#######################

# rm(list = ls())
# load("simdatafinal.RData")
# load("censusdata.RData")
# load("censusdrug.RData")
# load("nsduhpop.RData")
# load("mj.pop.RData")
# load("cocaine.pop.RData")
# load("meth.pop.RData")
# load("power.sim.RData")
# save.image("DrugUseBias.RData")
