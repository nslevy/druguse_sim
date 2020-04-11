###############################################
# Simulation of Misclassification of Drug Use #
###############################################

# Load needed packages
library(tidyverse)
library(survey, quietly = T)
library(srvyr, quietly= T)

setwd("/Users/nslevy/Documents/Columbia PhD/Misclassification of Drug Use")


##########################################
# 1 - Create simulated demographics data #
##########################################

# Read in 2018 Census ACS Data Query, cleaned up a little bit in Excel because Census MDAT output is v. messy
# Create labeled factor versions of numeric demographic variables
raw.data <- 
  read.csv("Census Custom LONG_Clean.csv") %>%
  mutate(SEX_FACT = factor(SEX, 
                           labels = c("Male","Female")),
         AGE_FACT = factor(AGE, 
                           labels = c("1-11","12-15","16-20","21-30","31-40","41-50","51-60","60+")),
         SCHL_FACT = factor(SCHL, 
                            labels = c("None/NA", "Less than highschool","High school","More than highschool"))) %>%
  arrange(SEX, AGE, SCHL)

# Convert Census data into proportions
total <- sum(raw.data[,4:10])

raw.prop <- 
  raw.data %>%
  mutate_at(c("AI.AN.NH","WHITE","BLACK","ASIAN","OTHER","MULT","HISP"), ~./total)

# Create rounded aggregated totals for simulated dataset of n=100000
raw.totals <-
  raw.prop %>%
    mutate_at(c("AI.AN.NH","WHITE","BLACK","ASIAN","OTHER","MULT","HISP"), ~.*100000) %>%
    mutate_at(c("AI.AN.NH","WHITE","BLACK","ASIAN","OTHER","MULT","HISP"), round)

# Turn into long dataset
# Create demographic variables that will match with NSDUH categories
# Create demographic variable flag for later merging purposes
census.data <-
  raw.totals %>%
  pivot_longer(cols = c("AI.AN.NH","WHITE","BLACK","ASIAN","OTHER","MULT","HISP"), names_to = "RACE_CHAR") %>%
  mutate(RACE = as.integer(case_when(
      RACE_CHAR == "AI.AN.NH" ~ 1,
      RACE_CHAR == "WHITE" ~ 2,
      RACE_CHAR == "BLACK" ~ 3,
      RACE_CHAR == "ASIAN" ~ 4,
      RACE_CHAR == "OTHER" ~ 5,
      RACE_CHAR == "MULT" ~ 6,
      RACE_CHAR == "HISP" ~ 7)),
      RACE_FACT = as.factor(RACE_CHAR),
      DEMFLAG = paste(SEX, AGE, SCHL, RACE, sep = "_"),
      AGE2 = ifelse(AGE==6, 5, AGE),
      RACE2 = ifelse(RACE==5, 6, RACE),
      SCHL2 = ifelse(SCHL==1, 2, SCHL),
      DEMFLAG2 = paste(SEX, AGE2, SCHL2, RACE2, sep = "_")) %>%
  select(SEX, SEX_FACT, AGE, AGE_FACT, SCHL, SCHL_FACT, RACE, RACE_FACT, DEMFLAG, DEMFLAG2, value) %>%
  arrange(SEX, AGE, SCHL, RACE)


#########################################################
# 2 - Create drug use variable based on 2018 NSDUH Data #
#########################################################

# Read in and recode 2018 NSDUH Data
load("/Users/nslevy/Documents/Columbia PhD/Silvia/NSDUH/Data/NSDUH_2018.RData")

names(PUF2018_100819) <- tolower(names(PUF2018_100819))

nsduh2018 <-
  PUF2018_100819 %>%
  select(questid2, irsex, ireduhighst2, catag6, catag7, newrace2, illyr, analwt_c, vestr, verep) %>%
  rename(SEX = irsex) %>%
  mutate(AGE = as.integer(case_when(
                  catag7 == 1 | catag7 == 2 ~ 2,
                  catag7 == 3 | catag7 == 4 ~ 3,
                  catag7 == 5 | catag7 == 6 ~ 4,
                  catag6 == 4 ~ 5,
                  catag6 == 5 ~ 7,
                  catag6 == 6 ~ 8)),
         SCHL = as.integer(case_when(
                  ireduhighst2 < 8 ~ 2,
                  ireduhighst2 == 8 ~ 3,
                  ireduhighst2 > 8 ~ 4)),
         RACE = as.integer(case_when(
                  newrace2 == 3 | newrace2 == 4 ~ 1,
                  newrace2 == 1 ~ 2,
                  newrace2 == 2 ~ 3,
                  newrace2 == 5 ~ 4,
                  newrace2 == 7 ~ 6,
                  newrace2 == 6 ~ 7)))

# Create survey design object
dnsduh <-
  svydesign(
    id = ~ verep, ### PSU - first stage cluster sampling
    strata = ~ vestr, ### strata
    weights = ~ analwt_c, ### probability of being sampled
    data = nsduh2018, ### data
    nest = TRUE ### clusters nested within strata
  )

# Calculate prevalence of any past-year ilicit drug use by sex, age, race, and education

# Totals by sex, age, race, and education
nsduh.pop <- 
  as.data.frame(svytable(~SEX + AGE + SCHL + RACE, dnsduh)) %>%
  mutate(DEMFLAG2 = paste(SEX, AGE, SCHL, RACE, sep = "_")) %>%
  arrange(SEX, AGE, SCHL, RACE)

# Totals by past-year drug use, sex, age, race, and education
drug.pop <- 
  as.data.frame(svytable(~illyr + SEX + AGE + SCHL + RACE, dnsduh)) %>%
  mutate(DEMFLAG2 = paste(SEX, AGE, SCHL, RACE, sep = "_")) %>%
  select(SEX, AGE, SCHL, RACE, DEMFLAG2, illyr, Freq) %>%
  arrange(SEX, AGE, SCHL, RACE, illyr)

# Proportions of past-year drug use by sex, age, race, and education
drug.prev <- 
  left_join(drug.pop, nsduh.pop, by = "DEMFLAG2") %>%
  mutate(DRUGPREV = ifelse(is.na(Freq.x/Freq.y), 0, Freq.x/Freq.y)) %>%
  rename(SEX = SEX.x,
         AGE = AGE.x,
         SCHL = SCHL.x,
         RACE = RACE.x,
         DRUGPOP = Freq.x,
         TOTALPOP = Freq.y) %>%
  select(SEX, AGE, SCHL, RACE, DEMFLAG2, illyr, DRUGPOP, TOTALPOP, DRUGPREV)


######################################################################
# 3 - Create combined dataset of demographic and drug use variables #
######################################################################

# Merge drug use prevalence with total counts by demographic characteristics from Census data
census.drug <- 
  full_join(drug.prev, census.data, by = "DEMFLAG2") %>%
  rename(SEX = SEX.y, AGE = AGE.y, SCHL = SCHL.y, RACE = RACE.y) %>%
  select(SEX, SEX_FACT, AGE, AGE_FACT, SCHL, SCHL_FACT, RACE, RACE_FACT,
         DEMFLAG, DEMFLAG2, illyr, DRUGPOP, TOTALPOP, DRUGPREV, value) %>%
  arrange(DEMFLAG, illyr)

# Generate simulation population sizes by past-year drug use
census.drug$DRUGCOUNT <- round(ifelse(is.na(census.drug$DRUGPREV), census.drug$value, 
                                census.drug$value*census.drug$DRUGPREV))


##############################################
# 4 - Generate biased prevalence of drug use #
##############################################

# Randomly generate proportions of bias for each category of demographics
# Assuming all those reporting no drug use are unbiased
# Assuming youngest age group does not misrepresent their drug use
set.seed(041020)
census.drug <-
  census.drug %>%
  mutate(BIAS = ifelse(illyr==0 | is.na(illyr), 1, ifelse(AGE==1, 0, 1 - runif(nrow(census.drug), 0, 0.5))))

# Create biased DRUGCOUNT variable
census.drug$DRUGCOUNT.BIAS <- rep(NA, nrow(census.drug))

for (i in 1:nrow(census.drug)){
  census.drug$DRUGCOUNT.BIAS[i] <-
    ifelse(is.na(census.drug$illyr[i]), census.drug$DRUGCOUNT[i],
        ifelse(census.drug$illyr[i]==0, 
               round((census.drug$DRUGCOUNT[i]*census.drug$BIAS[i])+
                       (census.drug$DRUGCOUNT[i+1]*(1-census.drug$BIAS[i+1]))),
               round(census.drug$DRUGCOUNT[i]*census.drug$BIAS[i])))
}


###################################################
# 5 - Generate individual-level simulated dataset #
###################################################

# Create individual observations based on TRUE drug use
sim.data <- 
  uncount(census.drug, DRUGCOUNT) %>%
  rowid_to_column() %>%
  mutate(TRUE.DRUGUSE = ifelse(illyr == 0 | is.na(illyr), 0, 1)) %>%
  select(rowid, SEX, SEX_FACT, AGE, AGE_FACT, SCHL, SCHL_FACT, RACE, RACE_FACT, DEMFLAG, DEMFLAG2, TRUE.DRUGUSE)

# Create individual observations based on BIASED drug use
sim.data2 <- 
  uncount(census.drug, DRUGCOUNT.BIAS) %>%
  rowid_to_column() %>%
  mutate(BIASED.DRUGUSE = ifelse(illyr == 0 | is.na(illyr), 0, 1)) %>%
  select(rowid, BIASED.DRUGUSE)

# Join dataframes
sim.data.final <- left_join(sim.data, sim.data2, by="rowid")

prop.table(table(sim.data.final$TRUE.DRUGUSE))
prop.table(table(sim.data.final$BIASED.DRUGUSE))


########################
# 6 - Correct the bias #
########################

# True counts of drug use by sex, age, education, and race
true.counts <- 
  sim.data.final %>%
  group_by(DEMFLAG, TRUE.DRUGUSE) %>%
  count(TRUE.DRUGUSE) %>%
  rename(TRUE.DRUGCOUNT = n)

# Biased counts of drug use by sex, age, education, and race
biased.counts <-
  sim.data.final %>%
  group_by(DEMFLAG, BIASED.DRUGUSE) %>%
  count(BIASED.DRUGUSE) %>%
  rename(BIASED.DRUGCOUNT = n)

# Corrected counts of drug use by sex, age, education, and race
corrected.data <- left_join(true.counts, biased.counts, by=c("DEMFLAG"="DEMFLAG", "TRUE.DRUGUSE" = "BIASED.DRUGUSE"))

corrected.data <-
  corrected.data %>%
  mutate(BIASPROP = BIASED.DRUGCOUNT/TRUE.DRUGCOUNT,
         CORRECTED.DRUGCOUNT = BIASED.DRUGCOUNT/BIASPROP) %>%
  rename (DRUGUSE = TRUE.DRUGUSE)

corrected.data %>%
  group_by(DRUGUSE) %>%
  summarise(TRUE.DRUGCOUNT = sum(TRUE.DRUGCOUNT),
            BIASED.DRUGCOUNT = sum(BIASED.DRUGCOUNT),
            CORRECTED.DRUGCOUNT = sum(CORRECTED.DRUGCOUNT))
