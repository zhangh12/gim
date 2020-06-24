\name{gim}
\alias{gim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Generalized Integration Models
}
\description{
\code{gim} is used to fit generalized integration models, which assume linear or logistic regression model on an (internal) individual-level data, while integrating auxiliary or summary information of relevant variables that are estimated from external data, on which different working models could be assumed. \code{gim} can work even if partial information from working models are available. Compared to conventional regression model, e.g., \code{\link[stats]{glm}}, that is based on internal data, the estimate of \code{gim} method gains additional power by making maximum use of all kinds of available data. 
}
\usage{
gim(formula, family, data, model, nsample = NULL, 
    ncase = NULL, nctrl = NULL, ref = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class "\code{\link[stats]{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted on the given dataset. More details of model specification are illustrated in 'Details' and 'Examples'.
}
  \item{family}{
a character. \code{"gaussian"} for linear regression. For binary outcome fitted by logistic regression, use \code{"binomial"} for random sample, or \code{"case-control"} for case-control data. \code{gim} employs different methods to make inference on random sample and case-control data. If your data are collected in case-control studies, do not use \code{"binomial"}, otherwise inference may be problematic. 
}
  \item{data}{
a data frame containing all variables that are specified in \code{formula} and \code{model}. Incomplete lines will be discarded. 
}
  \item{model}{
a list describing auxiliary information and working models that are used to generate such information. See 'Details' and 'Examples' for more details.
}
  \item{nsample}{
a matrix specifying the number of samples shared in datasets that are used to fit the working models given in \code{model}. Specify this argument when \code{family} is \code{"gaussian"} or \code{"binomial"}, otherwise \code{NULL}. See 'Details' and 'Examples' for more details.
}
  \item{ncase}{
a matrix specifying the number of cases shared in datasets that are used to fit the working models given in \code{model}. Specify this argument when \code{family} is \code{"case-control"}, otherwise \code{NULL}. \code{ncase} and \code{nctrl} should be specified simultaneously. See 'Details' and 'Examples' for more details.
}
  \item{nctrl}{
a matrix specifying the number of controls shared in datasets that are used to fit the working models given in \code{model}. Specify this argument when \code{family} is \code{case-control}, otherwise \code{NULL}. See 'Details' and 'Examples' for more details.
}
  \item{ref}{
a data frame containing the covariates specified in \code{formula} and \code{model}. It is the reference sample for modeling summary statistics in \code{model}. This assumes that \code{ref} is sampled from the external population. By default it is \code{NULL} which means that the internal and external populations are the same, therefore \code{gim} will use \code{data} as the reference. Outcome could be absent or missing in \code{ref} because \code{gim} will anyway ignore it. See 'Details' for more details. 
}
  \item{...}{
for test purpose, use its default value. 
}
}
\details{

\bold{\code{formula} }
\code{formula} is the model to be used to fit a conventional regression model if no additional information is available. It could be very general as long as it is acceptable to the \code{glm} or \code{lm} functions. It can eliminate the intercept, \code{y ~ .-1}, or involve arithmetic expressions, e.g., \code{log(x)}, or other operators like \code{*} for interactions \code{as.factor(x1)*I(x2 > 0)}.

\bold{\code{model} }
Summary information are calculated on data of external studies, but we do not have access to their raw data. Instead, estimates from working model fitted on external data are given (e.g., reported in literature). The argument \code{model} is a list, each component contains information of a working model. Specifically, a component is also a list of two entries \code{form} and \code{info}, where \code{form} is a formula representing the fitted working model, and \code{info} is a data frame with two columns \code{var} and \code{bet}, the names of variables and their estimates from the working model, respectively. Usually the estimate of intercept of a working model is unavailable as people fit but do not reporte it. If user is able to provide such an estimate, the name in column \code{var} must be \code{"(Intercept)"}. See below for an example. 

Note that multiple working models could be fitted on the same external data, in that case, the summary information of each working model should be given in \code{model} separately. For example, on an external dataset, if two models \code{y ~ x1} and \code{y ~ x2} are fitted, then the estimates of \code{x1} and \code{x2} should be given as two components in \code{model}. This happens as many research groups can study the same datasets from different angles. 

\bold{\code{data} }
\code{gim} requires an internal dataset \code{data} in which individual-level samples are available. Statistically, this data is critical to provide information of correlation between covariates. This data is also known as the reference data in the literatures. Since general formula is supported in \code{gim}, it is important to provide variables in \code{data} so that \code{R} can find columns of all variables parsed from formulas in \code{formula} and \code{model}. Read vignettes (upcoming) for more examples about how to create a proper \code{data} for \code{gim}. We will also release a function to help users with this. \code{gim} will discard incomplete lines in \code{data}. 

\bold{\code{nsample} }
Some of summary information can be calculated from datasets that share samples. Ignoring this will lead to underestimated standard error. For example, if a dataset is studied by two different models, the estimates from these two models are not independent but highly correlated. Therefore, this correlation must be properly handled when calculating the standard error of \code{gim} estimate, from which a hypothesis testing is conducted. \code{nsample} is a squared matrix of dimension \code{p}, which is equal to the length of \code{model}. Thus, the (i,i) entry in \code{nsample} is the number of samples used in fitting the working model specified in \code{model[[i]]$formula}, while the (i,k) entry is the number of samples that are involved in fitting working models \code{model[[i]]$formula} and \code{model[[k]]$formula}. For example, if two working models, e.g., \code{y ~ x1} and \code{y ~ x2} are fitted on the same dataset of 100 samples, then \code{nsample} is a matrix of all entries being 100. Read example below and vignettes (upcoming) for more examples. 

\bold{\code{ncase} } and \bold{\code{nctrl} }
Specify these two arguments when data are sampled from case-control studies. Refer to \code{nsample} for their formats. 

\bold{\code{ref} }
By default, \code{ref} is \code{NULL} if it is not specified explicitly. This assumes that the internal and external populations are the same, and \code{gim} will assign \code{data} to \code{ref} implicitly. If this assumption holds, and you have additional covariates data (no outcome), e.g. \code{add.ref}, that also comes from the internal population, you can specified \code{ref} as \code{rbind(data, add.ref)} where the column of missing outcome in \code{add.ref} is set as \code{NA}. You can also rbind \code{data} and \code{add.ref}, with outcome in \code{data} being deleted. If the external population is different from the internal population, you have to assign \code{add.ref} to \code{ref} as reference. 

}
\value{
\code{gim} returns an object of class "\code{gim}". The function \code{\link[base]{summary}} can be used to print a summary of the results. We will support the use of \code{\link[stats]{anova}} in later versions. 

The generic accessor functions \code{\link[stats]{coef}}, \code{\link[stats]{confint}}, and \code{\link[stats]{vcov}} can be used to extract coefficients, confidence intervals, and variance-covariance of estimates from the object returned by \code{gim}. 

An object of class "\code{gim}" is a list containing the following components:
  \item{coefficients }{a named vector of coefficients}
  \item{vcov }{the variance-covariance matrix of estimates, including the intercept}
  \item{sigma2 }{estimated variance of error term in a linear model. Only available for the \code{gaussian} family}
  \item{call }{the matched call}
  \item{V.bet }{the variance-covariance matrix of external estimate \code{bet} in \code{model}}
%% ...
}
\references{
Zhang, H., Deng, L., Schiffman, M., Qin, J., Yu, K. (2020) Generalized integration model for improved statistical inference by leveraging external summary data. Biometrika. asaa014, https://doi.org/10.1093/biomet/asaa014
}
\author{
Han Zhang
}
\examples{
## An artificial dataset is lazyloaded to illustrate the concept of GIM method
## It contains:
## A continuous outcome y. 
## Four covariates x1, x2, x3, x4 (character).
## A binary outcome d

head(dat)

## internal data of 500 samples
dat0 <- dat[1:500, ]

## three external datasets.
## dat2 and dat3 share some samples
dat1 <- dat[501:1500, c('y', 'x1', 'x2')]
dat2 <- dat[1501:2500, c('y', 'x1', 'x3', 'x4')]
dat3 <- dat[2001:3000, c('y', 'x3', 'x4')]

## four working models are fitted
form1 <- 'y ~ I(x1 < 0) + I(x2 > 0)'
form2 <- 'y ~ x3 + x4'
form3 <- 'y ~ I(x4 == "a")'
form4 <- 'y ~ sqrt(x3)'

## two working models are fitted on dat3
## thus nsample is a 4x4 matrix
nsample <- matrix(c(1000, 0, 0, 0, 
                    0, 1000, 500, 500, 
                    0, 500, 1000, 1000, 
                    0, 500, 1000, 1000), 
                  4, 4)

fit1 <- summary(lm(form1, dat1))$coef
fit2 <- summary(lm(form2, dat2))$coef
fit3 <- summary(lm(form3, dat3))$coef ## <-- dat3 is used twice
fit4 <- summary(lm(form4, dat3))$coef ## <-- dat3 is used twice

options(stringsAsFactors = FALSE)
model <- list()
## partial information is available
model[[1]] <- list(form = form1, 
                   info = data.frame(var = rownames(fit1)[2], 
                                     bet = fit1[2, 1]))

## intercept is provided, but miss estimate of a covariate
model[[2]] <- list(form = form2, 
                   info = data.frame(var = rownames(fit2)[1:2], 
                                     bet = fit2[1:2, 1]))

model[[3]] <- list(form = form3, 
                   info = data.frame(var = rownames(fit3)[2], 
                                     bet = fit3[2, 1]))

model[[4]] <- list(form = form4, 
                   info = data.frame(var = rownames(fit4)[2], 
                                     bet = fit4[2, 1]))

form <- 'y ~ I(x1 < 0) + I(x1 > 1) + x2 * x4 + log(x3) - 1'
fit <- gim(form, 'gaussian', dat0, model, nsample)

summary(fit)
coef(fit)
confint(fit)

# one can compare the gim estimates with those estimated from internal data
fit0 <- lm(form, dat0)
summary(fit0)

# by default, covariates in dat is used as reference in gim
# which assumes that the external and internal populations are the same
fit1 <- gim(form, 'gaussian', dat0, model, nsample, ref = dat0)
all(coef(fit) == coef(fit1)) # TRUE

# if additional reference is available, 
# and it comes from the internal population from which dat is sampled
# gim can use it
add.ref <- dat[3001:3500, ]
add.ref$y <- NA ## <-- outcome is unavailable in reference
ref <- rbind(dat0, add.ref)
fit2 <- gim(form, 'gaussian', dat0, model, nsample, ref = ref)

# if the external population is different from the internal population
# then reference for summary data specified in model needs to be provided
ext.ref <- dat[3501:4000, ] ## <-- as an example, assume ext.ref is different 
                            ##    from dat0
fit3 <- gim(form, 'gaussian', dat0, model, nsample, ref = ext.ref)

}
