pkgname <- "mglm4twin"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "mglm4twin-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('mglm4twin')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("anthro")
### * anthro

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: anthro
### Title: Anthropometric measures (weight and height)
### Aliases: anthro
### Keywords: datasets

### ** Examples

require(mglm4twin)
data(anthro, package="mglm4twin")
anthro$age <- (anthro$age - mean(anthro$age))/sd(anthro$age)
anthro$weight <- (anthro$weight - mean(anthro$weight))/sd(anthro$weight)
anthro$height <- (anthro$height - mean(anthro$height))/sd(anthro$height)
form_Wt <- weight ~ age + Group*Twin_pair
form_Ht <- height ~ age + Group*Twin_pair
biv0 <- list("formE1" = ~ age, "formE2" = ~ age, "formE12" = ~ age,
             "formA1" = ~ age, "formA2" = ~ age, "formA12" = ~ age,
             "formC1" = ~ age, "formC2" = ~ age, "formC12" = ~ age)
Z_biv0 <- mt_twin(N_DZ = 327, N_MZ = 534, n_resp = 2, model = "ACE",
                 formula = biv0, data = anthro)
control_initial <- list()
control_initial$regression <- list("R1" = c(0.13, 0.10, -0.20, -0.02, 0.037),
                                   "R2" = c(0.23, 0.01, -0.27, -0.11, 0.11))
control_initial$power <- list(c(0), c(0))
control_initial$tau <- c(0.15, 0, 0.12, rep(0,15))
fit_0 <- mglm4twin(linear_pred = c(form_Wt, form_Ht), matrix_pred = Z_biv0,
                   control_initial = control_initial,
                   control_algorithm = list(tuning = 0.5),
                   power_fixed = c(TRUE, TRUE), data = anthro)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("anthro", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bmi")
### * bmi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bmi
### Title: Body mass index
### Aliases: bmi
### Keywords: datasets

### ** Examples

require(mglm4twin)
data(bmi, package="mglm4twin")
form = bmi ~ Group*Twin_pair
ACE = mt_twin(N_DZ = 5576/2, N_MZ = 2966/2, n_resp = 1, model = "ACE")
fit_ACE <- mglm4twin(linear_pred = c(form), matrix_pred = ACE, data = bmi)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bmi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bpdrds")
### * bpdrds

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bpdrds
### Title: Bronchopulmonary dysplasia and respiratory distress syndrome on
###   preterm infants
### Aliases: bpdrds
### Keywords: datasets

### ** Examples

require(mglm4twin)
data(bpdrds, package="mglm4twin")
form_BPD <- BPD ~ BW + GA + gender + Group*Twin_pair
form_RDS <- RDS ~ BW + GA + gender + Group*Twin_pair
AE <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 2, model = "AE")
fitAE <- mglm4twin(linear_pred = c(form_BPD, form_RDS), matrix_pred = AE,
                   link = c("logit","logit"),
                   variance = c("binomialP","binomialP"), data = bpdrds)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bpdrds", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_link_function")
### * mt_link_function

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_link_function
### Title: Link Functions
### Aliases: mt_link_function mt_logit mt_probit mt_cauchit mt_cloglog
###   mt_loglog mt_identity mt_log mt_sqrt mt_invmu2 mt_inverse

### ** Examples

x1 <- seq(-1, 1, l = 5)
X <- model.matrix(~ x1)
mt_link_function(beta = c(1,0.5), X = X,
                 offset = NULL, link = 'log')
mt_link_function(beta = c(1,0.5), X = X,
                 offset = rep(10,5), link = 'identity')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_link_function", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_matrix_linear_predictor")
### * mt_matrix_linear_predictor

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_matrix_linear_predictor
### Title: Matrix Linear Predictor
### Aliases: mt_matrix_linear_predictor

### ** Examples

require(Matrix)
Z0 <- Diagonal(5, 1)
Z1 <- Matrix(rep(1,5)%*%t(rep(1,5)))
Z <- list(Z0, Z1)
mt_matrix_linear_predictor(tau = c(1,0.8), Z = Z)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_matrix_linear_predictor", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mt_variance_function")
### * mt_variance_function

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mt_variance_function
### Title: Variance Functions
### Aliases: mt_variance_function mt_tweedie mt_binomialP mt_binomialPQ
###   mt_constant

### ** Examples

x1 <- seq(-1, 1, l = 5)
X <- model.matrix(~x1)
mu <- mt_link_function(beta = c(1, 0.5), X = X, offset = NULL,
                       link = "logit")
mt_variance_function(mu = mu$mu, power = c(2, 1), Ntrial = 1,
                     variance = "binomialPQ",
                     derivative_power = TRUE, derivative_mu = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mt_variance_function", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("psydis")
### * psydis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: psydis
### Title: Psychiatric disorders
### Aliases: psydis
### Keywords: datasets

### ** Examples

require(mglm4twin)
data(psydis, package="mglm4twin")
ex1_form <- y ~ 1
ex1_AE <- mt_twin(N_DZ = 440, N_MZ = 590, n_resp = 1, model = "AE")
ex1_AE <- mglm4twin(c(ex1_form), matrix_pred = ex1_AE,
                    link = c("logit"), variance = c("binomialP"),
                   data = psydis)
summary(ex1_AE, model = "AE")
summary(ex1_AE, model = "AE", biometric = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("psydis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("t0psqi")
### * t0psqi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: t0psqi
### Title: Sleep’s quality
### Aliases: t0psqi
### Keywords: datasets

### ** Examples

require(mglm4twin)
form_T0 <- T0 ~ Age + Gender + Group + Type*Twin_pair
form_PSQI <- PSQI ~ Age + Gender + Group + Type*Twin_pair
AE <- mt_twin(N_DZ = 135, N_MZ = 116, n_resp = 2, model = "AE")
fit_AE <- mglm4twin(linear_pred = c(form_T0,  form_PSQI),
                   matrix_pred = AE,
                   link = c("log","logit"),
                   variance = c("tweedie","binomialP"),
                   control_algorithm = list(tuning = 0.25, max_iter = 100),
                   power_fixed = c(FALSE,FALSE), data = t0psqi)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("t0psqi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
