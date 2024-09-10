
## Load packages
library(mglm4twin)
library(readxl)
library(dplyr)
library(tidyr)

## Load data set
url <- "/home/wagner/Downloads/Random test data.xlsx"
dados <- read_excel(url)
head(dados)
View(dados)

## Checking number of observations per pair
dados |>
  group_by(ID_pair) |>
  summarize("N" = n()) |>
  arrange(N)

# It's fine !

## Ordering zygosity
dados$Zygosity <- as.factor(dados$Zygosity)
dados <- dados |>
  arrange(Zygosity)
View(dados)

# Fill NA (just to have the complete data set)
m1 = mean(dados$Right_upper_meibography, na.rm = TRUE)
m2 = mean(dados$Right_lower_meibography, na.rm = TRUE)
m3 = mean(dados$Left_upper_meibography, na.rm = TRUE)
m4 = mean(dados$Left_lower_meibography, na.rm = TRUE)
dados <- dados |>
  replace_na(list(Right_upper_meibography = m1 ,
                  Right_lower_meibography = m2,
                  Left_upper_meibography = m3,
                  Left_lower_meibography =  m4))
View(dados)


## Data set is fine!

# Step 1: Linear predictor
form_rum <- Right_upper_meibography ~ 1
form_rlm <- Right_lower_meibography ~ 1
form_lum <- Left_upper_meibography ~ 1
form_llm <- Left_lower_meibography ~ 1

# Step 2: Matrix linear predictor
table(dados$Zygosity)/2
mult_E <- mt_twin(N_DZ = 21, N_MZ = 23, n_resp = 4, model = "E")
length(mult_E) # 10 dispersion parameters
lapply(mult_E, dim) # Size of each matrix should be dim(dados)[1]*n_resp - OK

# Step 3: Link and variance functions
link_fc <- rep("logit", 4)
variance_fc <- rep("binomialP", 4)

# Step 4: Fitting the model
fitE <- mglm4twin(linear_pred = c(form_rum, form_rlm, form_lum, form_llm),
                  matrix_pred = mult_E,
                  link = link_fc, variance = variance_fc,
                  data = dados)

summary(fitE, model = "E")


## Model ACE
mult_AE <- mt_twin(N_DZ = 21, N_MZ = 23, n_resp = 4, model = "AE")
fitAE <- mglm4twin(linear_pred = c(form_rum, form_rlm, form_lum, form_llm),
                   matrix_pred = mult_AE,
                   link = link_fc, variance = variance_fc,
                   control_algorithm = list(verbose = TRUE, tuning = 0.5,
                                            max_iter = 200, correct = TRUE),
                   data = dados)
summary(fitAE, model = "AE")
summary(fitAE, model = "AE", biometric = TRUE)
