
gen.example <- function(){
  
  set.seed(1)
  
  n <- 600
  bmi <- stats::rchisq(n, df = 6)
  bmi <- round(16 + (40-16)/(max(bmi)-min(bmi)) * (bmi-min(bmi)), digits = 1)
  
  age <- stats::rnorm(n)
  age <- round(40 + (80-40)/(max(age)-min(age)) * (age-min(age)), digits = 0)
  
  breastfeed <- sample(c('yes', 'no'), n, TRUE, c(.85, .15))
  
  parity <- sample(0:5, n, TRUE, c(.05, 0.15, 0.45, 0.3, 0.03, 0.02))
  
  p53 <- sample(1:3, n, TRUE, c(.3, .3, .4))
  
  subtype <- sample(c('A', 'B', 'C'), n, TRUE, c(.65, .2, .15))
  
  dat <- data.frame(subj.id = paste0('SID-',1:n), density = 0, bmi, age, breastfeed, parity, p53, subtype, stringsAsFactors = FALSE)
  
  mm <- model.matrix(density~I(age >= 50 & age < 60) + I(age >= 60) + as.factor(breastfeed) + I(parity > 0) + I(p53 !=1)*I(subtype != 'A'), data=dat)
  
  para <- c(1, 0.1, 0.2, -0.1, -0.3, 0.05, 0., 0.15)
  
  dat$density <- exp(as.vector(mm %*% para) + stats::rnorm(n, sd = sqrt(0.05)))
  
  int <- dat[1:200, ]
  ext1 <- dat[201:400, c('subj.id', 'density', 'age', 'parity', 'subtype')]
  ext2 <- dat[201:300, c('subj.id', 'density', 'age', 'breastfeed', 'subtype')]
  ext3 <- dat[301:600, c('subj.id', 'density', 'p53', 'subtype')]
  
  int$age_50to59 <- I(int$age >= 50 & int$age < 60)+1-1
  int$age_gt60 <- I(int$age >= 60)+1-1
  formula <- as.formula(log(density) ~ age_50to59 + age_gt60 + breastfeed + parity + as.factor(p53)*I(subtype=='A'))
  
  fit0 <- glm(formula, data = int, family = 'gaussian')
  fit1 <- glm(log(density) ~ I(age < 60) + I(parity > 0) + subtype, data = ext1, family = 'gaussian')
  fit2 <- glm(log(density) ~ I(age >= 70) + breastfeed + subtype, data = ext2, family = 'gaussian')
  fit31 <- glm(log(density) ~ as.factor(p53), data = ext3, family = 'gaussian')
  fit32 <- glm(log(density) ~ subtype, data = ext3, family = 'gaussian')
  
  summary(fit0)
  summary(fit1)
  summary(fit2)
  
  cbind(colnames(mm), para)
  
  model <- list()
  #model[[1]] <- list(formula='log(density) ~ I(age < 60) + I(parity > 0) + subtype', 
  #                   info=data.frame(var = names(coef(fit1))[3:5], 
  #                              bet = coef(fit1)[3:5]))
  
  model[[1]] <- list(formula='log(density) ~ I(age >= 70) + breastfeed + subtype', 
                     info=data.frame(var = names(coef(fit2))[2:4], 
                                bet = coef(fit2)[2:4]))
  
  #model[[2]] <- list(formula='log(density) ~ as.factor(p53)', 
  #                   info=data.frame(var = names(coef(fit31))[-1], 
  #                              bet = coef(fit31)[-1]))
  
  model[[2]] <- list(formula='log(density) ~ subtype', 
                     info=data.frame(var = names(coef(fit32))[-1], 
                                bet = coef(fit32)[-1]))
  
  n1 <- nrow(ext1)
  n2 <- nrow(ext2)
  n3 <- nrow(ext3)
  n0 <- length(intersect(ext1$subj.id, ext3$subj.id))
  #nsample <- matrix(c(n1, n0, 0, 0, n0, n2, 0, 0, 0, 0, n3, n3, 0, 0, n3, n3), 4,4)
  nsample <- matrix(c(n1, n0, n0, n3), 2, 2)
  
  form <- formula
  dat <- int
  #fit <- gim(form, 'gaussian', dat, model, nsample)
  
  lt <- list(form = form, dat = dat, model = model, nsample = nsample)
  #lt <- list(fit = fit, form = form, dat = dat, model = model, nsample = nsample)
  #save(list = names(lt), file = 'data/example.rda')
  lt
  
}

