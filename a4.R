#=========================================================================
# Econ 613 A3
# Yuqi Zhou
#=========================================================================

library(tidyverse)
library(readr)
library(esquisse)
library(dplyr)
library(AER)
library(panelr)
getwd()
# Exercise 1 Preparing the Data
dat_A4 <- read_csv("/Users/jennyzhou/Desktop/ECON 613/a4/dat_A4.csv")
dat_A4_panel <- read_csv("/Users/jennyzhou/Desktop/ECON 613/a4/dat_A4_panel.csv")
drop1 <- which(dat_A4$CV_HGC_BIO_DAD_1997 == 95)
drop2 <- which(dat_A4$CV_HGC_BIO_MOM_1997 == 95)
drop3 <- which(dat_A4$CV_HGC_RES_DAD_1997 == 95)
drop4 <- which(dat_A4$CV_HGC_RES_MOM_1997 == 95)
dat_A4$CV_HGC_BIO_DAD_1997[drop1] = NA
dat_A4$CV_HGC_BIO_MOM_1997[drop2] = NA
dat_A4$CV_HGC_RES_DAD_1997[drop3] = NA
dat_A4$CV_HGC_RES_MOM_1997[drop4] = NA

# (a)
dat_A4 <- dat_A4 %>%
  mutate(age = 2019 - KEY_BDATE_Y_1997,
         work_exp = rowSums(dat_A4[,18:28], na.rm=TRUE) / 48)
# (b)
dat_A4 <- dat_A4 %>%
  mutate(parent_edu = rowSums(dat_A4[,8:11], na.rm=TRUE) / 4,
         highest_degree = case_when(YSCH.3113_2019 == 1 ~ "<12",
                                    YSCH.3113_2019 == 2 ~ "~12",
                                    YSCH.3113_2019 == 3 ~ "12",
                                    YSCH.3113_2019 == 4 ~ "14",
                                    YSCH.3113_2019 == 5 ~ "16",
                                    YSCH.3113_2019 == 6 ~ "18",
                                    YSCH.3113_2019 == 7 ~ "21",
                                    YSCH.3113_2019 == 8 ~ "21"), na.rm=TRUE)
dat_A4$highest_degree <- as.numeric(dat_A4$highest_degree, na.rm=TRUE)

# (c)
income_data <- dat_A4 %>% 
  filter(YINC_1700_2019 != '0' & YINC_1700_2019 != 'NA')

# Age 
age_income <- income_data %>% 
  group_by(age) %>%
  summarise(average_wage = mean(YINC_1700_2019, na.rm = TRUE))

ggplot(age_income) +
  aes(x = age, y = average_wage) +
  geom_line(size = 0.5, colour = "#112446") +
  labs(x = "Age", y = "Average Income") +
  theme_minimal()

# Gender
male_data <- income_data %>% 
  filter(KEY_SEX_1997 == '1')
ggplot(male_data) +
  aes(x = YINC_1700_2019) +
  geom_histogram(bins = 30L, fill = "#112446") +
  labs(x = "Male Income Distribution") +
  theme_minimal()

female_data <- income_data %>% 
  filter(KEY_SEX_1997 == '1')
ggplot(female_data) +
  aes(x = YINC_1700_2019) +
  geom_histogram(bins = 30L, fill = "#112446") +
  labs(x = "Female Income Distribution") +
  theme_minimal()

# Number of children
children_income <- income_data %>% 
  group_by(CV_BIO_CHILD_HH_U18_2019) %>%
  summarise(average_wage = mean(YINC_1700_2019, na.rm = TRUE))

ggplot(children_income) +
  aes(x = CV_BIO_CHILD_HH_U18_2019, y = average_wage) +
  geom_line(size = 0.5, colour = "#112446") +
  labs(x = "Number of Children", y = "Average Income") +
  theme_minimal()

dat_A4$YINC_1700_2019[is.na(dat_A4$YINC_1700_2019)] <- 0
table <- dat_A4 %>% 
  group_by(age) %>%
  summarise(Number_of_People = n(),
            No_income_People = length(which(YINC_1700_2019 == 0)),
            Share = No_income_People / Number_of_People)

# Exercise 2 Heckman Selection Model
# (a)
OLS_reg <- lm(YINC_1700_2019 ~ age + work_exp + parent_edu + KEY_SEX_1997, CV_MARSTAT_COLLAPSED_2019, data = income_data)
summary(OLS_reg) 
# There might be a selection problem because we rule out all the 0s and NAs in the dataset.

# (b)
# (c)
dat_A4 <- dat_A4 %>% mutate(Intercept = 1, Income_or_not = 0) 
dat_A4$Income_or_not[which(dat_A4$YINC_1700_2019 > 0)] <- 1
x1 = dat_A4$Intercept
x2 = dat_A4$age
x3 = dat_A4$work_exp
x4 = dat_A4$parent_edu
x5 = dat_A4$KEY_SEX_1997
x6 = dat_A4$CV_MARSTAT_COLLAPSED_2019
# Probit Function
flike = function(par,x1,x2,x3,x4,x5,x6,Income_or_not){
  yhat = par[1]*x1 + par[2]*x2 + par[3]*x3 + par[4]*x4 + par[5]*x5
  prob = pnorm(yhat)
  prob[prob>0.999999] = 0.999999
  prob[prob<0.000001] = 0.000001
  like = Income_or_not*log(prob) + (1-Income_or_not)*log(1-prob)
  return(-sum(like))
}

predict <- function(par,x1,x2,x3,x4,x5,x6,Income_or_not){
  yhat = par[1]*x1 + par[2]*x2 + par[3]*x3 + par[4]*x4 + par[5]*x5
  return(yhat)
}
Income_or_not = dat_A4$Income_or_not
start <- runif(7, -1, 1)
result  <- optim(start, fn = flike, method="BFGS", control = list(trace = 6, REPORT = 1, maxit = 1000),
                 x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6, Income_or_not = Income_or_not, hessian = TRUE)
result$par

predictor <- predict(result$par,x1,x2,x3,x4,x5,x6,Income_or_not)
Inv_M_ratio <- dnorm(predictor) / pnorm(predictor)
Hecman_reg <- lm(dat_A4$YINC_1700_2019 ~ x2 + x3 + x4 + x5 + x6 + Inv_M_ratio)

# Exercise 3 Censoring
# (a)
ggplot(income_data) +
  aes(x = YINC_1700_2019) +
  geom_histogram(bins = 30L, fill = "#112446") +
  theme_minimal()
# From the histogram, we can get that the censored value is $100,000.

# (b)
# We can use the Tobit model to fix the censored problem.

# (c)
# Use the Tobit package to get a good initial par
reg_tobit <- tobit(dat_A4$YINC_1700_2019 ~ x2 + x3 + x4 + x5 + x6, left=-Inf,right = 10000)
summary(reg_tobit)
initial_par <- as.vector(c(reg_tobit$coefficients,10))
start <- initial_par #+ runif(7, -10, 10)

# Set up data
dat_A4 <- dat_A4 %>% 
  filter(dat_A4$YINC_1700_2019 != 0 & dat_A4$YINC_1700_2019 != 'NA') %>%
  mutate(censored = 1) 
dat_A4$censored[which(dat_A4$YINC_1700_2019 < 100000)] <- 0
censored <- dat_A4$censored
income <- dat_A4$YINC_1700_2019
x1 = dat_A4$Intercept
x2 = dat_A4$age
x3 = dat_A4$work_exp
x4 = dat_A4$YSCH.3113_2019
x5 = dat_A4$KEY_SEX_1997
x6 = dat_A4$CV_MARSTAT_COLLAPSED_2019

# This is the Tobit log-likelihood
tobit_likelihood <- function(par,x1,x2,x3,x4,x5,x6,censored,income){
  yhat = par[1]*x1 + par[2]*x2 + par[3]*x3 + par[4]*x4 + par[5]*x5 + par[6]*x6
  sigma = exp(par[7])
  residual = income - yhat
  standardization = (100000-yhat)/sigma
  like = censored*log(1 - pnorm(standardization)) + (1-censored)*log(dnorm(residual/sigma)/sigma) 
  return(-sum(like))
}

result2 <- optim(start,fn=tobit_likelihood,method="BFGS",control=list(trace=6,REPORT=1,maxit=1000),x1=NULL,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,censored=censored,income=income,hessian=TRUE)
result2$par
summary(reg_tobit)


# (d) 
# Exercise 4 Panel Data
# Convert the panel data from wide to long
dat_A4_panel <- read_csv("dat_A4_panel.csv")
dat_A4_panel <- dat_A4_panel %>% rename(CV_HIGHEST_DEGREE_EVER_EDT_1998 = CV_HIGHEST_DEGREE_9899_1998,
                                        CV_HIGHEST_DEGREE_EVER_EDT_1999 = CV_HIGHEST_DEGREE_9900_1999,
                                        CV_HIGHEST_DEGREE_EVER_EDT_2000 = CV_HIGHEST_DEGREE_0001_2000,
                                        CV_HIGHEST_DEGREE_EVER_EDT_2001 = CV_HIGHEST_DEGREE_0102_2001, 
                                        CV_HIGHEST_DEGREE_EVER_EDT_2002 = CV_HIGHEST_DEGREE_0203_2002,
                                        CV_HIGHEST_DEGREE_EVER_EDT_2003 = CV_HIGHEST_DEGREE_0304_2003,
                                        CV_HIGHEST_DEGREE_EVER_EDT_2004 = CV_HIGHEST_DEGREE_0405_2004,
                                        CV_HIGHEST_DEGREE_EVER_EDT_2005 = CV_HIGHEST_DEGREE_0506_2005, 
                                        CV_HIGHEST_DEGREE_EVER_EDT_2006 = CV_HIGHEST_DEGREE_0607_2006, 
                                        CV_HIGHEST_DEGREE_EVER_EDT_2007 = CV_HIGHEST_DEGREE_0708_2007, 
                                        CV_HIGHEST_DEGREE_EVER_EDT_2008 = CV_HIGHEST_DEGREE_0809_2008,
                                        CV_HIGHEST_DEGREE_EVER_EDT_2009 = CV_HIGHEST_DEGREE_0910_2009,
                                        CV_HIGHEST_DEGREE_EVER_EDT_2010 = CV_HIGHEST_DEGREE_1011_2010)
dat_A4_long_panel <- long_panel(dat_A4_panel, prefix='_', begin  = 1997, end = 2019, label_location = "end")

# Create education, marital status and working experience variables in the panel data.
dat_A4_long_panel <- as.data.frame(dat_A4_long_panel)
dat_A4_long_panel <- dat_A4_long_panel %>%
  mutate(work_exp = rowSums(dat_A4_long_panel[,c(10:16,35:37)], na.rm = TRUE) / 48,
         education = dat_A4_long_panel$CV_HIGHEST_DEGREE_EVER_EDT,
         marital_status = dat_A4_long_panel$CV_MARSTAT_COLLAPSED)

#Between Estimator: exper/edu/marital status
estimate <- dat_A4_long_panel %>% 
  group_by(id) %>% 
  mutate(m_income = mean(`YINC-1700`,na.rm = TRUE),
         m_work_exp = mean(work_exp,na.rm = TRUE),
         m_education = mean(education,na.rm = TRUE),
         m_marital_status = mean(marital_status,na.rm = TRUE))

between_estimator <- lm(m_income ~ m_work_exp + m_education + m_marital_status, data = estimate)
summary(between_estimator)

# Within Estimator: exper/edu/marital status
estimate <- estimate %>% 
  group_by(id) %>% 
  mutate(diff_income = `YINC-1700` - m_income,
         diff_education = education - m_education,
         diff_work_exp = work_exp - m_work_exp,
         diff_marital_status = marital_status - m_marital_status)

within_estimator <- lm(diff_income ~ diff_education + diff_work_exp + diff_marital_status, data = estimate)
summary(within_estimator)

# Difference Estimator: exper/edu/marital status(first difference)
estimate <- estimate %>% 
  group_by(id) %>% 
  mutate(first_diff_income = `YINC-1700` - lag(`YINC-1700`),
         first_diff_education = education - lag(education),
         first_diff_work_exp = work_exp - lag(work_exp),
         first_diff_marital_status = marital_status - lag(marital_status))
         
difference_estimator <- lm(first_diff_income ~ first_diff_education + first_diff_work_exp + first_diff_marital_status, data = estimate)
summary(difference_estimator)

