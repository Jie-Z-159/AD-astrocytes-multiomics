#NHANES analysis 

###Extract the required variables -----
library(haven)
library(plyr)
library(dplyr) 
library(arsenal) 
library(survey)
library(nhanesA)
demo.g <- read_xpt("DEMO_G.XPT")
demo.h <- read_xpt("DEMO_H .XPT")
smq.h<-nhanes('SMQ_H')
smq.g <- read_xpt("SMQ_G.XPT")
smq.data.file <- dplyr::bind_rows(list(smq.g, smq.h))
colnames(smq.data.file)
smq.data <- smq.data.file[,c('SEQN', 'SMQ020', 'SMQ040', 'SMQ050Q', 'SMQ050U')]
alq.g <- read_xpt("ALQ_G.XPT")
alq.h <- read_xpt("ALQ_H.XPT")
alq.data.file <- dplyr::bind_rows(list(alq.g, alq.h))
colnames(alq.data.file)
alq.data <- alq.data.file[,c('SEQN', 'ALQ101', 'ALQ110', 'ALQ120Q', 'ALQ120U')]
bmx.g <- read_xpt("BMX_G.XPT") 
bmx.h <- read_xpt("BMX_H.XPT")
bmx.data.file <- dplyr::bind_rows(list(bmx.g, bmx.h))
bmx.data <- bmx.data.file[,c('SEQN', 'BMXBMI', 'BMXWAIST')]
dr1tot.g <- read_xpt('DR1TOT_G.XPT')
dr1tot.h <- read_xpt('DR1TOT_H.XPT')
dr1tot.data.file <- dplyr::bind_rows(list(dr1tot.g, dr1tot.h))
dr1tot.data <- dr1tot.data.file[,c('SEQN', 'DR1TSFAT','DR1TMFAT','DR1TPFAT', 'DR1TKCAL')]###我改
dr2tot.g <- read_xpt('DR2TOT_G.XPT')
dr2tot.h <- read_xpt('DR2TOT_H.XPT')
dr2tot.data.file <- dplyr::bind_rows(list(dr2tot.g, dr2tot.h))
dr2tot.data <- dr2tot.data.file[,c('SEQN', 'DR2TSFAT','DR2TMFAT','DR2TPFAT', 'DR2TKCAL')]
dr.data <- merge(dr2tot.data, dr1tot.data)
ds1tot.g <- read_xpt('DS1TOT_G.XPT')
ds1tot.h <- read_xpt('DS1TOT_H.XPT')
ds1tot.data.file <- dplyr::bind_rows(list(ds1tot.g, ds1tot.h))
ds1tot.data <- ds1tot.data.file[,c('SEQN', 'DS1TSFAT','DS1TMFAT','DS1TPFAT')]
ds2tot.g <- read_xpt('DS2TOT_G.XPT')
ds2tot.h <- read_xpt('DS2TOT_H.XPT')
ds2tot.data.file <- dplyr::bind_rows(list(ds2tot.g, ds2tot.h))
ds2tot.data <- ds2tot.data.file[,c('SEQN', 'DS2TSFAT','DS2TMFAT','DS2TPFAT')]
ds.data <- merge(ds1tot.data, ds2tot.data)
cfq.g <- read_xpt('CFQ_G.XPT')
cfq.h <- read_xpt('CFQ_H.XPT')
cfq.data.file <- dplyr::bind_rows(list(cfq.g, cfq.h))
cfq.data <- cfq.data.file[,c('SEQN', 'CFDCST1', 'CFDCST2', 'CFDCST3',
                             'CFDCSR', 'CFDAST', 'CFDDS')]
weight.data <- dr1tot.data.file[,c('SEQN', 'WTDRD1')]
weight.data$WTDRD1 <- weight.data$WTDRD1/2 
survey.design.data <- demo.data.file[,c('SEQN', 'SDMVPSU', 'SDMVSTRA')]
output <- plyr::join_all(list(demo.data, smq.data, alq.data, bmx.data,
                        dr.data, ds.data, cfq.data, weight.data, survey.design.data),
                        by='SEQN', type='left')

#quality control
data.age.60 <- subset.data.frame(output, RIDAGEYR >= 60 )
dim(data.age.60) 
which(!is.na(data.age.60$WTDRD1))
weight.non.na.index <- which(!is.na(data.age.60$WTDRD1))
sum(data.age.60$WTDRD1[weight.non.na.index]) 
paper.data <- subset.data.frame(output, RIDAGEYR >= 60 &
                                  (!is.na(CFDCSR)) & 
                                  (!is.na(output$CFDAST)) &
                                  (!is.na(output$CFDDS))& 
                                  (!is.na(output$CFDCST1))&
                                  (!is.na(output$CFDCST2))&
                                  (!is.na(output$CFDCST3))&
                                  (!is.na(DR1TSFAT))&
                                  (!is.na(DR1TMFAT))&
                                  (!is.na(DR1TPFAT))& 
                                  (!is.na(DR1TKCAL))& 
                                  (!is.na(RIDAGEYR)) & 
                                  (!is.na(RIAGENDR)) & 
                                  (!is.na(RIDRETH3)) & 
                                  (!is.na(DMDEDUC2)) & 
                                  (!is.na(SMQ020)) 
                                 
)
dim(paper.data)
paper.exclude.data <- subset.data.frame(output, RIDAGEYR >= 60 &(
                                          (is.na(CFDCSR)) |
                                          (is.na(CFDAST)) |
                                          (is.na(CFDDS)) | 
                                            (is.na(CFDCST1))|
                                            (is.na(CFDCST2))|
                                            (is.na(CFDCST3))|
                                            (is.na(DR1TSFAT))|
                                            (is.na(DR1TMFAT))|
                                          (is.na(DR1TPFAT))| 
                                          (is.na(DR1TKCAL))| 
                                          (is.na(RIDAGEYR)) | 
                                          (is.na(RIAGENDR)) | 
                                          (is.na(RIDRETH3)) | 
                                          (is.na(DMDEDUC2)) | 
                                          (is.na(SMQ020))) 
)
dim(paper.exclude.data)

index.60.na <- which(output$RIDAGEYR <60)
length(index.60.na)#16299=19931-3632
index_60_plus <- which(output$RIDAGEYR >= 60)
index_DR_na <- which(!is.na(output$DR1TSFAT[index_60_plus]))
index_60_DR_na <- index_60_plus[index_DR_na]


index_CERAD_na <-  (is.na(output$CFDCSR)) | 
  (is.na(output$CFDAST)) |
  (is.na(output$CFDDS))| 
  (is.na(output$CFDCST1))|
  (is.na(output$CFDCST2))|
  (is.na(output$CFDCST3))
index.na <- which(index_CERAD_na[index_60_DR_na])
length(index.na) 

dim(paper.exclude.data)   
dim(paper.data)    
paper.exclude.data$Type <- 'exclude'
paper.data$Type <- 'include'
compare.data <- dplyr::bind_rows(paper.data, paper.exclude.data)

###Variable Derivation----
paper.data$Sex <- ifelse(paper.data$RIAGENDR == 1, 'male', 'female')
table(paper.data$Sex)
paper.data$age.group <- ifelse(paper.data$RIDAGEYR >= 60 & paper.data$RIDAGEYR < 69, '60-69 years',
                               ifelse(paper.data$RIDAGEYR >=70 & paper.data$RIDAGEYR < 79, '70-79 years',
                                      '80+ years'))
table(paper.data$age.group)
table(paper.data$RIDRETH3)
race <- recode_factor(paper.data$RIDRETH3, 
                      `1` = 'Mexican American',
                      `2` = 'Other Hispanic',
                      `3` = 'Non-Hispanic White',
                      `4` = 'Non-Hispanic Black',
                      `6` = 'Other/multiracial',
                      `7` = 'Other/multiracial'
)
paper.data$race <- race
table(paper.data$race)
ddply(paper.data, .(ALQ101, ALQ110), summarize, n = length(SEQN))
ddply(alq.g, .(ALQ101, ALQ110), summarize, n = length(SEQN))
ori.alq.unit <- paper.data$ALQ120U
table(ori.alq.unit)
trans.unit.month <- ifelse(ori.alq.unit == 1, 4, 
                           ifelse(ori.alq.unit == 3, 1/12, 
                                  ifelse(ori.alq.unit == 7|ori.alq.unit == 9, NA, 1)))
paper.data$trans.unit.month <- trans.unit.month
ori.alq.quantity <- paper.data$ALQ120Q
trans.quantity.month <- ifelse(ori.alq.quantity >= 0,  
                               ori.alq.quantity * trans.unit.month, NA)
paper.data$trans.quantity.month <- trans.quantity.month
alq101 <- paper.data$ALQ101
trans.quantity.month.factor <- ifelse(trans.quantity.month >=1 & trans.quantity.month <5, '1-5 drinks/month',
                                      ifelse(trans.quantity.month >=5 & trans.quantity.month <10, '5-10 drinks/month',
                                             ifelse(trans.quantity.month >=10, '10+ drinks/month', 'wait')))
paper.data$trans.quantity.month.factor <- trans.quantity.month.factor
ddply(paper.data, .(ALQ101, trans.quantity.month.factor), summarise, n = length(SEQN))
table(trans.quantity.month.factor) 
index.1 <- which((trans.quantity.month.factor=='wait' | is.na(trans.quantity.month.factor)) & alq101 == 1) #844
trans.quantity.month.factor[index.1] <- '1-5 drinks/month'
table(trans.quantity.month.factor)
paper.data$trans.quantity.month.factor <- trans.quantity.month.factor
ddply(paper.data, .(ALQ101, trans.quantity.month.factor), summarise, n = length(SEQN))
index.less.1 <- which(alq101 == 2)
trans.quantity.month.factor[index.less.1] <- 'Non-drinker'
table(trans.quantity.month.factor)
paper.data$alq.group <- trans.quantity.month.factor
table(paper.data$alq.group)
paper.data$total.calories <- paper.data$DR1TKCAL
CERAD.total <- apply(paper.data[,c('CFDCST1', 'CFDCST2', 'CFDCST3')], 1, sum)
paper.data$CERAD.total <- CERAD.total
paper.data$BMI.group <- ifelse(paper.data$BMXBMI <18.5, 'Underweight(<18.5)',
                               ifelse(paper.data$BMXBMI >=18.5 & paper.data$BMXBMI < 25, 'Normal(18.5 to <25)',
                                      ifelse(paper.data$BMXBMI >=25 & paper.data$BMXBMI < 30, 'Overweight(25 to <30)',
                                             'Obese(30 or greater)')))
paper.data <- mutate(paper.data, smoke.group = case_when(
  SMQ020 == 2 ~ 'Never smoker',
  SMQ020 == 1 & SMQ040 == 3 ~ 'Former smoker',
  SMQ020 == 1 & SMQ040 <= 2 ~ 'Current smoker'
))
education.attainment <- recode_factor(paper.data$DMDEDUC2, 
                                      `1` = 'Less Than 9th Grade',
                                      `2` = '9-11th Grade',
                                      `3`= 'High School Grad/GED',
                                      `4`= 'Some College or AA degree',
                                      `5`= 'College Graduate or above')

paper.data$education.attainment <- education.attainment

paper.data <- subset.data.frame(paper.data, 
                                (!is.na(SFAT)) &
                                  (!is.na(MFAT)) &
                                  (!is.na(PFAT)) &
                                  (!is.na(smoke.group))& 
                                  (!is.na(education.attainment))&
                                  (!is.na(CFDCSR)) & #CERAD 延迟回忆评分
                                  (!is.na(CFDAST)) & #语言流畅性评分Animal Fluency
                                  (!is.na(CFDDS))&
                                  (!is.na(alq.group))) #数字符号替代测试(DSST)) #是否吸烟至少100支

dim(paper.data) 
colnames(paper.data)
analyze.variable <- c("SEQN", "WTDRD1", "SDMVPSU", "SDMVSTRA",
                      "Sex", "RIDAGEYR", "age.group", "race", "education.attainment", "INDFMPIR",
                      "alq.group", "smoke.group","BMXBMI", "BMI.group", "BMXWAIST", 
                      "total.calories", "SFAT", "MFAT","PFAT",
                      "CFDCST1", "CFDCST2", "CFDCST3",  "CERAD.total", 
                      "CFDCSR", "CFDAST", "CFDDS")

paper.data <- paper.data[, analyze.variable]
colnames(paper.data) <- c("SEQN", "WTDRD1", "SDMVPSU", "SDMVSTRA",
                          "Sex", 'Age', "Age.group", "Race", "Education.attainment", "PIR",
                          "Alq.group", "Smoke.group", "BMI", "BMI.group", "Waist", 
                          "Total.calories",  "Total.SFAT",
                         "Total.MFAT",
                          "Total.PFAT",
                          "CERAD1", "CERAD2", "CERAD3", "CERAD.total",
                          "CERAD.delay.recall", "Animal.Fluency", "DSST")
#Generate variables by fatty acid intake ratio----
paper.datanew<-mutate(paper.data,totFA=Total.SFAT+Total.MFAT+Total.PFAT)
paper.datanew<-mutate(paper.datanew,raSFAT=Total.SFAT/totFA)
paper.datanew<-mutate(paper.datanew,raMFAT=Total.MFAT/totFA)
paper.datanew<-mutate(paper.datanew,raPFAT=Total.PFAT/totFA)
paper.data<-paper.datanew
paper.data<-read.csv(file="paperdata_new.csv",sep=',',header = T)

###Generate complex sampling objects-----
NHANES_design <- svydesign(data = paper.data, ids = ~SDMVPSU, strata = ~SDMVSTRA, nest = TRUE, weights = ~WTDRD1, survey.lonely.psu = "adjust") 
paper.data$Alq.group <- factor(paper.data$Alq.group, 
                               levels = c('1-5 drinks/month', 
                                          '5-10 drinks/month',
                                          '10+ drinks/month',
                                          'Non-drinker'))
paper.data$BMI.group <- factor(paper.data$BMI.group, 
                               levels = c('Underweight(<18.5)', 
                                          'Normal(18.5 to <25)',
                                          'Overweight(25 to <30)',
                                          'Obese(30 or greater)'))
levels(fct_infreq(paper.data$Race)) 
paper.data$Race <- factor(paper.data$Race, levels(fct_infreq(paper.data$Race)))
svytotal(~ Race, NHANES_design, na.rm=TRUE) 
NHANES_design <- svydesign(data = paper.data, ids = ~SDMVPSU, strata = ~SDMVSTRA, nest = TRUE, weights = ~WTDRD1, survey.lonely.psu = "adjust") 

#### 2. Table1 ####
tbl <- tbl_svysummary(NHANES_design,  
by = Sex,
include = c(Age.group, Sex,  Race, BMI.group, Alq.group, Smoke.group, Education.attainment,
            Age, PIR, BMI, Waist, Total.calories,  Total.SFAT,
            Total.MFAT,Total.PFAT,
            CERAD1, CERAD2, CERAD3, CERAD.total, CERAD.delay.recall, Animal.Fluency, DSST),
label = list(Age ~ 'Age (years)', PIR ~ 'PIR'),
#type = list(Age ~ "categorical"),
statistic = list(all_continuous()  ~ "{median} ({p25}, {p75})", 
                 all_categorical() ~ "{n_unweighted} ({p}%)"),
digits = list(Age ~ 1, PIR ~ 2),  
sort = list(Race ~ "frequency", BMI.group ~ "alphanumeric" ),
missing = 'no') %>% 
  add_overall() %>%
  modify_header(
    all_stat_cols() ~ "**{level}**, N = {n_unweighted} ({style_percent(p)}%)",
  ) %>% # update the column header
  modify_spanning_header(
    stat_0 ~ NA, 
    update = all_stat_cols() ~ "**Sex**") %>% 
  modify_footnote(
    update = all_stat_cols() ~ "median (IQR) for continuous; n (%) for categorical") %>%
  bold_labels() %>%
  italicize_levels() 
#tbl
tbl %>%
as_flex_table() %>% 
 flextable::save_as_docx(path = 'ratable1.docx')

##### 3 table2 Total ####

raSFAT.quantile.res <- svyquantile(~raSFAT, NHANES_design, quantiles = c(0, 0.25, 0.5, 0.75, 1))
paper.data$raSFAT.quantile.var <- cut(paper.data$raSFAT,
                                          breaks =raSFAT.quantile.res$raSFAT[,'quantile'],
                                          labels = c('Q1', 'Q2', 'Q3', 'Q4'))

paper.data$raSFAT.quantile.var[which(is.na(paper.data$raSFAT.quantile.var))] <- 'Q1'
NHANES_design <- svydesign(data = paper.data, ids = ~SDMVPSU, strata = ~SDMVSTRA, nest = TRUE, weights = ~WTDRD1, survey.lonely.psu = "adjust") 
raSFATTable2<-tbl_svysummary(NHANES_design, by =raSFAT.quantile.var, missing = 'no',
               include = c(raSFAT, Animal.Fluency, CERAD1, CERAD2, CERAD3, CERAD.total, CERAD.delay.recall,DSST),
               label = list(raSFAT ~ 'raSFAT', Animal.Fluency ~ 'Animal Fluency: Score Total')) %>%           
  modify_header(all_stat_cols() ~ "**{level}**, N = {n_unweighted} ({style_percent(p)}%)") %>%               
  add_p() %>%as.data.frame()

tbl_svysummary(NHANES_design, by =raSFAT.quantile.var, missing = 'no',
               include = c(raSFAT, Animal.Fluency, CERAD1, CERAD2, CERAD3, CERAD.total, CERAD.delay.recall,DSST),
               label = list(raSFAT ~ 'raSFAT intake', Animal.Fluency ~ 'Animal Fluency: Score Total')) %>%               
  modify_header(all_stat_cols() ~ "**{level}**, N = {n_unweighted} ({style_percent(p)}%)") %>%               
  add_p() %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = 'raSFATTable2_Result.docx')

#MFAT
raMFAT.quantile.res <- svyquantile(~raMFAT, NHANES_design, quantiles = c(0, 0.25, 0.5, 0.75, 1))
paper.data$raMFAT.quantile.var <- cut(paper.data$raMFAT,
                                          breaks =raMFAT.quantile.res$raMFAT[,'quantile'],
                                          labels = c('Q1', 'Q2', 'Q3', 'Q4'))

paper.data$raMFAT.quantile.var[which(is.na(paper.data$raMFAT.quantile.var))] <- 'Q1'
NHANES_design <- svydesign(data = paper.data, ids = ~SDMVPSU, strata = ~SDMVSTRA, nest = TRUE, weights = ~WTDRD1, survey.lonely.psu = "adjust") 
raMFATTable2<-tbl_svysummary(NHANES_design, by =raMFAT.quantile.var, missing = 'no',
               include = c(raMFAT, Animal.Fluency, CERAD1, CERAD2, CERAD3, CERAD.total, CERAD.delay.recall,DSST),
               label = list(raMFAT ~ 'raMFAT intake', Animal.Fluency ~ 'Animal Fluency: Score Total')) %>%              
  modify_header(all_stat_cols() ~ "**{level}**, N = {n_unweighted} ({style_percent(p)}%)") %>%               
  add_p()



#PFAT
raPFAT.quantile.res <- svyquantile(~raPFAT, NHANES_design, quantiles = c(0, 0.25, 0.5, 0.75, 1))
paper.data$raPFAT.quantile.var <- cut(paper.data$raPFAT,
                                          breaks =raPFAT.quantile.res$raPFAT[,'quantile'],
                                          labels = c('Q1', 'Q2', 'Q3', 'Q4'))

paper.data$raPFAT.quantile.var[which(is.na(paper.data$raPFAT.quantile.var))] <- 'Q1'
NHANES_design <- svydesign(data = paper.data, ids = ~SDMVPSU, strata = ~SDMVSTRA, nest = TRUE, weights = ~WTDRD1, survey.lonely.psu = "adjust") 

raPFATTable2<-tbl_svysummary(NHANES_design, by =raPFAT.quantile.var, missing = 'no',
               include = c(raPFAT, Animal.Fluency, CERAD1, CERAD2, CERAD3, CERAD.total, CERAD.delay.recall,DSST),
               label = list(raPFAT ~ 'raPFAT intake', Animal.Fluency ~ 'Animal Fluency: Score Total')) %>%            
  modify_header(all_stat_cols() ~ "**{level}**, N = {n_unweighted} ({style_percent(p)}%)") %>%               
  add_p()
raPFATTable2<-as.data.frame(raPFATTable2)

tbl_svysummary(NHANES_design, by =raPFAT.quantile.var, missing = 'no',
               include = c(raPFAT, Animal.Fluency, CERAD1, CERAD2, CERAD3, CERAD.total, CERAD.delay.recall,DSST),
               label = list(raPFAT ~ 'raPFAT intake', Animal.Fluency ~ 'Animal Fluency: Score Total')) %>%              
  modify_header(all_stat_cols() ~ "**{level}**, N = {n_unweighted} ({style_percent(p)}%)") %>%               
  add_p() %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = 'raPFATTable2_Result.docx')

#### 4. Table3  #####

# Fully adjusted
SFAT.full <- svyglm(CERAD.delay.recall ~ raSFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
SFAT.Q <- svyglm(CERAD.delay.recall ~ raSFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)
###MFAT
# Fully adjusted
MFAT.full <- svyglm(CERAD.delay.recall ~ raMFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
MFAT.Q <- svyglm(CERAD.delay.recall ~ raMFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)

###PFAT
# Fully adjusted
PFAT.full <- svyglm(CERAD.delay.recall ~ raPFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
PFAT.Q <- svyglm(CERAD.delay.recall ~ raPFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)
#DSST
tbl_SFAT_full_CERAD.delay.recall <- tbl_regression(SFAT.full, include = c(raSFAT))
tbl_SFAT_Q_CERAD.delay.recall <- tbl_regression(SFAT.Q, include = c(raSFAT.quantile.var))

tbl_MFAT_full_CERAD.delay.recall <- tbl_regression(MFAT.full, include = c(raMFAT))
tbl_MFAT_Q_CERAD.delay.recall <- tbl_regression(MFAT.Q, include = c(raMFAT.quantile.var))

tbl_PFAT_full_CERAD.delay.recall <- tbl_regression(PFAT.full, include = c(raPFAT))
tbl_PFAT_Q_CERAD.delay.recall <- tbl_regression(PFAT.Q, include = c(raPFAT.quantile.var))



###CERAD.total----
###SFAT
# Fully adjusted
SFAT.full <- svyglm(CERAD.total ~ raSFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
SFAT.Q <- svyglm(CERAD.total ~ raSFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)
###MFAT
# Fully adjusted
MFAT.full <- svyglm(CERAD.total ~ raMFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
MFAT.Q <- svyglm(CERAD.total ~ raMFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)

###PFAT
# Fully adjusted
PFAT.full <- svyglm(CERAD.total ~ raPFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
PFAT.Q <- svyglm(CERAD.total ~ raPFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)
#CERAD.total
tbl_SFAT_full_CERAD.total <- tbl_regression(SFAT.full, include = c(raSFAT))
tbl_SFAT_Q_CERAD.total <- tbl_regression(SFAT.Q, include = c(raSFAT.quantile.var))

tbl_MFAT_full_CERAD.total <- tbl_regression(MFAT.full, include = c(raMFAT))
tbl_MFAT_Q_CERAD.total <- tbl_regression(MFAT.Q, include = c(raMFAT.quantile.var))

tbl_PFAT_full_CERAD.total <- tbl_regression(PFAT.full, include = c(raPFAT))
tbl_PFAT_Q_CERAD.total <- tbl_regression(PFAT.Q, include = c(raPFAT.quantile.var))



###Animal.Fluency----
###SFAT
# Fully adjusted
SFAT.full <- svyglm(Animal.Fluency ~ raSFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
SFAT.Q <- svyglm(Animal.Fluency ~ raSFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)
###MFAT
# Fully adjusted
MFAT.full <- svyglm(Animal.Fluency ~ raMFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
MFAT.Q <- svyglm(Animal.Fluency ~ raMFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)

###PFAT
# Fully adjusted
PFAT.full <- svyglm(Animal.Fluency ~ raPFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
PFAT.Q <- svyglm(Animal.Fluency ~ raPFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)
#
tbl_SFAT_full_Animal.Fluency <- tbl_regression(SFAT.full, include = c(raSFAT))
tbl_SFAT_Q_Animal.Fluency <- tbl_regression(SFAT.Q, include = c(raSFAT.quantile.var))

tbl_MFAT_full_Animal.Fluency <- tbl_regression(MFAT.full, include = c(raMFAT))
tbl_MFAT_Q_Animal.Fluency <- tbl_regression(MFAT.Q, include = c(raMFAT.quantile.var))

tbl_PFAT_full_Animal.Fluency <- tbl_regression(PFAT.full, include = c(raPFAT))
tbl_PFAT_Q_Animal.Fluency <- tbl_regression(PFAT.Q, include = c(raPFAT.quantile.var))




###DSST----
###SFAT
# Fully adjusted
SFAT.full <- svyglm(DSST ~ raSFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
SFAT.Q <- svyglm(DSST ~ raSFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)
###MFAT
# Fully adjusted
MFAT.full <- svyglm(DSST ~ raMFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
MFAT.Q <- svyglm(DSST ~ raMFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)

###PFAT
# Fully adjusted
PFAT.full <- svyglm(DSST ~ raPFAT + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                    design = NHANES_design)
# Q4 vs Q1 
PFAT.Q <- svyglm(DSST ~ raPFAT.quantile.var + Age + Sex + BMI + Alq.group + Smoke.group + PIR + Education.attainment, 
                 design = NHANES_design)
#DSST
tbl_SFAT_full_DSST <- tbl_regression(SFAT.full, include = c(raSFAT))
tbl_SFAT_Q_DSST <- tbl_regression(SFAT.Q, include = c(raSFAT.quantile.var))

tbl_MFAT_full_DSST <- tbl_regression(MFAT.full, include = c(raMFAT))
tbl_MFAT_Q_DSST <- tbl_regression(MFAT.Q, include = c(raMFAT.quantile.var))

tbl_PFAT_full_DSST <- tbl_regression(PFAT.full, include = c(raPFAT))
tbl_PFAT_Q_DSST <- tbl_regression(PFAT.Q, include = c(raPFAT.quantile.var))
#table3----
SFAT_merged_table <- tbl_merge(
  list(tbl_SFAT_full_CERAD.total, tbl_SFAT_Q_CERAD.total, tbl_SFAT_full_CERAD.delay.recall, tbl_SFAT_Q_CERAD.delay.recall,
       tbl_SFAT_full_DSST, tbl_SFAT_Q_DSST, tbl_SFAT_full_Animal.Fluency, tbl_SFAT_Q_Animal.Fluency),
  tab_spanner = c("CERAD.total ", "CERAD.total Quantile", "CERAD.delay.recall ", "CERAD.delay.recall Quantile", 
                  "DSST ", "DSST Quantile","Animal.Fluency ", "Animal.Fluency Quantile")
)
SFAT_ft_table <- as_flex_table(SFAT_merged_table)
SFAT_ft_table <- flextable::set_caption(SFAT_ft_table, caption = "SFAT")

MFAT_merged_table <- tbl_merge(
  list(tbl_MFAT_full_CERAD.total, tbl_MFAT_Q_CERAD.total, tbl_MFAT_full_CERAD.delay.recall, tbl_MFAT_Q_CERAD.delay.recall,
       tbl_MFAT_full_DSST, tbl_MFAT_Q_DSST, tbl_MFAT_full_Animal.Fluency, tbl_MFAT_Q_Animal.Fluency),
  tab_spanner = c("CERAD.total ", "CERAD.total Quantile", "CERAD.delay.recall ", "CERAD.delay.recall Quantile", 
                  "DSST ", "DSST Quantile","Animal.Fluency ", "Animal.Fluency Quantile")
)
MFAT_ft_table <- as_flex_table(MFAT_merged_table)
MFAT_ft_table <- flextable::set_caption(MFAT_ft_table, caption = "MFAT")

PFAT_merged_table <- tbl_merge(
  list(tbl_PFAT_full_CERAD.total, tbl_PFAT_Q_CERAD.total, tbl_PFAT_full_CERAD.delay.recall, tbl_PFAT_Q_CERAD.delay.recall,
       tbl_PFAT_full_DSST, tbl_PFAT_Q_DSST, tbl_PFAT_full_Animal.Fluency, tbl_PFAT_Q_Animal.Fluency),
  tab_spanner = c("CERAD.total ", "CERAD.total Quantile", "CERAD.delay.recall ", "CERAD.delay.recall Quantile", 
                  "DSST ", "DSST Quantile","Animal.Fluency ", "Animal.Fluency Quantile")
)
PFAT_ft_table <- as_flex_table(PFAT_merged_table)
PFAT_ft_table <- flextable::set_caption(PFAT_ft_table, caption = "PFAT")
flextable::save_as_docx(SFAT_ft_table, path = 'raSFAT_Table3.docx')
flextable::save_as_docx(MFAT_ft_table, path = 'raMFAT_Table3.docx')
flextable::save_as_docx(PFAT_ft_table, path = 'raPFAT_Table3.docx')

#all_objects <- mget(ls())
#saveRDS(all_objects, file = "all_objects.rds")

###Visualization of the results of table 2----
library(ggplot2)
library(MoMAColors)
library(patchwork)
library(reshape2)
###Animal.Fluency

long_format_data <- melt(paper.data, id.vars = "Animal.Fluency", 
                         measure.vars = c("raSFAT.quantile.var", "raMFAT.quantile.var", "raPFAT.quantile.var"),
                         variable.name = "FattyAcidType", value.name = "Quantile")

long_format_data$FattyAcidType <- gsub("raSFAT.quantile.var", "SFAT", long_format_data$FattyAcidType)
long_format_data$FattyAcidType <- gsub("raMFAT.quantile.var", "MFAT", long_format_data$FattyAcidType)
long_format_data$FattyAcidType <- gsub("raPFAT.quantile.var", "PFAT", long_format_data$FattyAcidType)
long_format_data$Quantile <- factor(long_format_data$Quantile, levels = c("Q1", "Q2", "Q3", "Q4"))

P1=ggplot(long_format_data, aes(x = interaction(Quantile, FattyAcidType), y = Animal.Fluency, fill = FattyAcidType)) +
  geom_boxplot() +
  labs(x = "Fatty Acid Type and Ratio Quartile", y = "Animal Fluency Score") +
  scale_fill_manual(values = c('#aee0f9','#89d4ce','#fbd3d7')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = c(2, 6, 10), y = rep(max(long_format_data$Animal.Fluency) * 1.1, 3), 
           label = c("SFAT (p<0.001)", "MFAT (p=0.4)", "PFAT (p=0.006)"),
           size = 4, hjust = 0.5, vjust = 0)
###CERAD.total
long_format_data <- melt(paper.data, id.vars = "CERAD.total", 
                         measure.vars = c("raSFAT.quantile.var", "raMFAT.quantile.var", "raPFAT.quantile.var"),
                         variable.name = "FattyAcidType", value.name = "Quantile")

long_format_data$FattyAcidType <- gsub("raSFAT.quantile.var", "SFAT", long_format_data$FattyAcidType)
long_format_data$FattyAcidType <- gsub("raMFAT.quantile.var", "MFAT", long_format_data$FattyAcidType)
long_format_data$FattyAcidType <- gsub("raPFAT.quantile.var", "PFAT", long_format_data$FattyAcidType)
long_format_data$Quantile <- factor(long_format_data$Quantile, levels = c("Q1", "Q2", "Q3", "Q4"))
long_format_data$FattyAcidType<-factor(long_format_data$FattyAcidType,levels = c("SFAT","MFAT","PFAT"))

P2=ggplot(long_format_data, aes(x = interaction(Quantile, FattyAcidType), y = CERAD.total, fill = FattyAcidType)) +
  geom_boxplot() +
  labs(x = "Fatty Acid Type and Ratio Quartile", y = "CERAD Total Score") +
  scale_fill_manual(values = c('#aee0f9','#89d4ce','#fbd3d7')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = c(2, 6, 10), y = rep(max(long_format_data$CERAD.total) * 1.1, 3), 
           label = c("SFAT (p=0.001)", "MFAT (p>0.9)", "PFAT (p=0.001)"),
           size = 4, hjust = 0.5, vjust = 0)

###CERAD.delay.recall

long_format_data <- melt(paper.data, id.vars = "CERAD.delay.recall", 
                                      measure.vars = c("raSFAT.quantile.var", "raMFAT.quantile.var", "raPFAT.quantile.var"),
                                      variable.name = "FattyAcidType", value.name = "Quantile")

long_format_data$FattyAcidType <- gsub("raSFAT.quantile.var", "SFAT", long_format_data$FattyAcidType)
long_format_data$FattyAcidType <- gsub("raMFAT.quantile.var", "MFAT", long_format_data$FattyAcidType)
long_format_data$FattyAcidType <- gsub("raPFAT.quantile.var", "PFAT", long_format_data$FattyAcidType)
long_format_data$Quantile <- factor(long_format_data$Quantile, levels = c("Q1", "Q2", "Q3", "Q4"))
long_format_data$FattyAcidType<-factor(long_format_data$FattyAcidType,levels = c("SFAT","MFAT","PFAT"))

P3=ggplot(long_format_data, aes(x = interaction(Quantile, FattyAcidType), y = CERAD.delay.recall, fill = FattyAcidType)) +
  geom_boxplot() +
  labs(x = "Fatty Acid Type and Ratio Quartile", y =  "CERAD Delay Recall Score") +
  scale_fill_manual(values = c('#aee0f9','#89d4ce','#fbd3d7')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = c(2, 6, 10), y = rep(max(long_format_data$CERAD.delay.recall) * 1.1, 3), 
           label = c("SFAT (p=0.015)", "MFAT (p=0.9)", "PFAT (p<0.001)"),
           size = 4, hjust = 0.5, vjust = 0)

###DSST

long_format_data <- melt(paper.data, id.vars = "DSST", 
                              measure.vars = c("raSFAT.quantile.var", "raMFAT.quantile.var", "raPFAT.quantile.var"),
                              variable.name = "FattyAcidType", value.name = "Quantile")

long_format_data$FattyAcidType <- gsub("raSFAT.quantile.var", "SFAT", long_format_data$FattyAcidType)
long_format_data$FattyAcidType <- gsub("raMFAT.quantile.var", "MFAT", long_format_data$FattyAcidType)
long_format_data$FattyAcidType <- gsub("raPFAT.quantile.var", "PFAT", long_format_data$FattyAcidType)
long_format_data$Quantile <- factor(long_format_data$Quantile, levels = c("Q1", "Q2", "Q3", "Q4"))
long_format_data$FattyAcidType<-factor(long_format_data$FattyAcidType,levels = c("SFAT","MFAT","PFAT"))

P4=ggplot(long_format_data, aes(x = interaction(Quantile, FattyAcidType), y = DSST, fill = FattyAcidType)) +
  geom_boxplot() +
  labs(x = "Fatty Acid Type and Ratio Quartile", y =  "DSST Score") +
  scale_fill_manual(values = c('#aee0f9','#89d4ce','#fbd3d7')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = c(2, 6, 10), y = rep(max(long_format_data$DSST) * 1.1, 3), 
           label = c("SFAT (p=0.026)", "MFAT (p=0.4)", "PFAT (p=0.049)"),
           size = 4, hjust = 0.5, vjust = 0)

library(ggpubr)
ggarrange(P2,P3,P4,P1,common.legend = T,legend = 'top',ncol = 2,nrow = 2)




