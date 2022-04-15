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

# Exercise 1 Preparing the Data
dat_A4 <- read_csv("Desktop/ECON 613/a4/dat_A4.csv")
dat_A4_panel <- read_csv("Desktop/ECON 613/a4/dat_A4_panel.csv")
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
  mutate(education = rowSums(dat_A4[,8:11], na.rm=TRUE) / 4)
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
OLS_reg <- lm(YINC_1700_2019 ~ age + work_exp + education + KEY_SEX_1997, CV_MARSTAT_COLLAPSED_2019, data = income_data)
summary(OLS_reg) 
# There might be a selection problem because we rule out all the 0s and NAs in the dataset.

# (b)

# (c)
dat_A4 <- dat_A4 %>% mutate(Intercept = 1, Income_or_not = 0) 
dat_A4$Income_or_not[which(dat_A4$YINC_1700_2019 > 0)] <- 1
x1 = dat_A4$Intercept
x2 = dat_A4$age
x3 = dat_A4$work_exp
x4 = dat_A4$education
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
start <- runif(6, -0.5, 0.5)
result  <- optim(start, fn = flike, method="BFGS", control = list(trace = 6, REPORT = 1, maxit = 1000),
              x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6, Income_or_not = dat_A4$Income_or_not, hessian = TRUE)
result$par
 
# Use the Probit package to check the result
reg_package <- glm(Income_or_not ~ x2 + x3 + x3 + x4 + x5 + x6, family = binomial(link = "probit"))
summary(reg_package)
reg_package$coefficients

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
reg_tobit <- tobit(income ~ x2 + x3 + x4 + x5 + x6, left=-Inf,right = 100000)
summary(reg_tobit)
initial_par <- as.vector(c(reg_tobit$coefficients))
start <- initial_par + runif(7, -10, 10)

# Set up data
dat_A4 <- dat_A4 %>% 
  filter(dat_A4$YINC_1700_2019 != 0 & dat_A4$YINC_1700_2019 != 'NA') %>%
  mutate(censored = 1) 
dat_A4$censored[which(dat_A4$YINC_1700_2019<100000)] <- 0
censored <- dat_A4$censored
income <- dat_A4$YINC_1700_2019
x1 = dat_A4$Intercept
x2 = dat_A4$age
x3 = dat_A4$work_exp
x4 = dat_A4$education
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

result <- optim(start,fn=tobit_likelihood,method="BFGS",control=list(trace=6,REPORT=1,maxit=1000),x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,censored=censored,income=income,hessian=TRUE)
result$par
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
         education = rowSums(dat_A4_long_panel[,19:21], na.rm=TRUE) / 4,
         marital_status = dat_A4_long_panel$CV_MARSTAT_COLLAPSED)

#Between Estimator: gender/exper/edu/marital status
between_estimate <- dat_A4_long_panel %>% 
  group_by(id) %>% 
  summarize(income = mean(`YINC-1700`,na.rm = TRUE),
            work_exp = mean(work_exp,na.rm = TRUE),
            education = mean(CV_HIGHEST_DEGREE_EVER_EDT,na.rm = TRUE),
            marital_status = mean(CV_MARSTAT_COLLAPSED,na.rm = TRUE))
between_estimator <- lm(income ~ work_exp + education + marital_status, data = between_estimate)
summary(between_estimator)

# Within Estimator: exper/edu/marital status
within_estimate <- dat_A4_long_panel %>% 
  mutate(difference_income = `YINC-1700` - ,
         difference_education = CV_HIGHEST_DEGREE_EVER_EDT - ,
         difference_work_exp = work_exp - ,
         difference_marital_status = CV_MARSTAT_COLLAPSED - ,)

dat_A4_long_panel$meanincome <- ave(dat_A4_long_panel$`YINC-1700`, dat_A4_long_panel$id, FUN=function(x)mean(x, na.rm=T))
data$meanedu <- ave(data$CV_HIGHEST_DEGREE_EVER_EDT, data$id, FUN=function(x)mean(x, na.rm=T)) 
data$meanexper <- ave(data$work_exp, data$id, FUN=function(x)mean(x, na.rm=T)) 
data$meanms<- ave(data$CV_MARSTAT_COLLAPSED, data$id, FUN=function(x)mean(x, na.rm=T))

data$d_income <- data$`YINC-1700` - data$meanincome
data$d_edu <- data$CV_HIGHEST_DEGREE_EVER_EDT - data$meanedu
data$d_exper <- data$work_exp - data$meanexper
data$d_ms <- data$CV_MARSTAT_COLLAPSED - data$meanms

panel_within_estimator <- lm(data$d_income~ data$d_edu + data$d_exper + data$d_ms)
summary(panel_within_estimator)
# Difference estimator

