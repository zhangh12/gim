\name{gim}
\alias{gim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Generalized Integration Models
}
\description{
\code{gim} is used to fit generalized integration models, which assume linear or logistic regression model on a given (internal) data, while integrating auxiliary or summary information of relevant variables that are estimated from external data, on which different working models could be assumed. Compared to conventional regression model, e.g., \code{\link[stats]{glm}}, based on internal data, the estimate of \code{gim} method gains additional power by making maximum use of all kinds of available data. 
}
\usage{
gim(formula, family, data, model, nsample, ref = NULL, niter = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class "\code{\link[stats]{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted on the given dataset. More details of model specification are illustrated in 'Details' and 'Examples'.
}
  \item{family}{
'gaussian' for linear regression or 'binomial' for logistic regression. 
}
  \item{data}{
a data frame containing the variables specified in \code{formula} and \code{model}. 
}
  \item{model}{
a list describing auxiliary information and working models that are used to generate such information. See 'Details' and 'Examples' for more details.
}
  \item{nsample}{
a matrix specifying the number of samples shared in datasets that are used to fit the working models given in \code{model}. See 'Details' and 'Examples' for more details.
}
  \item{ref}{
a data frame containing the covariates specified in \code{formula} and \code{model}. It is the reference sample for modeling summary statistics in \code{model}. This assumes that \code{ref} is sampled from the external population. 
}
  \item{niter}{
for test purpose, use its default value. 
}
}
\details{

\bold{\code{formula} }
\code{formula} is the model to be used to fit a conventional regression model if no additional information is available. It could be very general as long as it is acceptable to the \code{glm} or \code{lm} functions. It can eliminate the intercept, \code{y ~ .-1}, or involve arithmetic expressions, e.g., \code{log(x)}, or other operators like \code{*} for interactions \code{as.factor(x1)*I(x2 > 0)}.

\bold{\code{model} }
Summary information are calculated on data of external studies, but we do not have access to their raw data. Instead, estimates from working model fitted on external data are given (e.g., reported in literature). The argument \code{model} is a list, each component contains information of a working model. Specifically, a component is also a list of two entries \code{form} and \code{info}, where \code{form} is a formula representing the fitted working model, and \code{info} is a data frame with two columns \code{var} and \code{bet}, the names of variables and their estiamtes from the working model, respectively. Usually the estimate of intercept of a working model is unavailable as people fit but do not reporte it. If user is able to provide such an estimate, the name in column \code{var} must be \code{"(Intercept)"}. See below for an example. 

Note that multiple working models could be fitted on the same external data, in that case, the summary information of each working model should be given in \code{model} separately. For example, on an external dataset, if two models \code{y ~ x1} and \code{y ~ x2} are fitted, then the estimates of \code{x1} and \code{x2} should be given as two components in \code{model}. This happens as many research groups can study the same datasets from different angles. 

\bold{\code{data} }
\code{gim} requires an internal dataset \code{data} in which individual-level samples are available. Statistically, this data is critical to provide information of correlation between covariates. This data is also known as the reference data in the literatures. Since general formula is supported in \code{gim}, it is important to provide variables in \code{data} so that \code{R} can find columns of all variables parsed from formulas in \code{formula} and \code{model}. Read vignettes (upcoming) for more examples about how to create a proper \code{data} for \code{gim}. We will also release a function to help users with this. 

\bold{\code{nsample} }
Some of summary information can be calculated from datasets that share samples. Ignoring this will lead to underestimated standard error. For example, if a dataset is studied by two different models, the estimates from these two models are not independent but highly correlated. Therefore, this correlation must be properly handled when calculating the standard error of \code{gim} estimate, from which a hypothesis testing is conducted. \code{nsample} is a squared matrix of dimension \code{p}, which is equal to the length of \code{model}. Thus, the (i,i) entry in \code{nsample} is the number of samples used in fitting the working model specified in \code{model[[i]]$form}, while the (i,k) entry is the number of samples that are involved in fitting working models \code{model[[i]]$form} and \code{model[[k]]$form}. For example, if two working models, e.g., \code{y ~ x1} and \code{y ~ x2} are fitted on the same dataset of 100 samples, then \code{nsample} is a matrix of all entries being 100. Read example below and vignettes (upcoming) for more examples. 

}
\value{
\code{gim} returns an object of class "\code{gim}". The function \code{\link[base]{summary}} can be used to print a summary of the results. We will support the use of \code{\link[stats]{anova}} in later versions. 

The generic accessor functions \code{\link[stats]{coefficients}}, \code{\link[stats]{confint}}, and \code{\link[stats]{vcov}} can be used to extract coefficients, confidence intervals, and variance-covariance of estimates from the object returned by \code{gim}. 

An object of class "\code{gim}" is a list containing the following components:
  \item{coefficients }{a named vector of coefficients}
  \item{vcov }{the variance-covariance matrix of estimates, including the intercept}
  \item{sigma2 }{estimated variance of error term in a linear model. Only available for the \code{gaussian} family}
  \item{call }{the matched call}
  \item{V.bet }{the variance-covariance matrix of external estimate \code{bet} in \code{model}}
%% ...
}
\references{
Zhang, H., Deng, L., Schiffman, M., Qin, J., Yu, K. (2018) Generalized integration model for improved statistical inference by leveraging external summary data. Under review. 
}
\author{
Han Zhang
}
\examples{
# All data used in this example are lazyloaded
# An artificial dataset is attached to illustrate the concept of GIM method
# It contains:
# A data frame "dat"" contains samples from an internal study. 
# A formula specified in "form" is used to describe the full model.
# A list "model" specifying summary information from (four) external studies. 
# A matrix "nsample" specifying sample size of (four) external studies. 
head(dat)
form

# model is a list, each component contains summary data fitted from 
# a working model.
#
# In this example, summary data are calculated from two working models, 
# thus model is a list of length 2. 
#
# Estimates of a subset of variables in model[[1]]$form are given in 
# model[[1]]$info, except for the intercept and subtypeC (estimated but missed). 
#
# If you want to provide estimate of intercept, you can add it to $info by 
# giving its name '(Intercept)'. See vignettes (upcoming ...). 
model

# nsample is a 2x2 matrix, with its (i,k) entry means the number of samples 
# that are used in fitting models i and k. 
#
# two working models are fitted on two datasets (n=200, 300), 
# shared 100 subjects 
nsample

fit <- gim(form, 'gaussian', dat, model[1], nsample[1, 1])
summary(fit)
coef(fit)
confint(fit)

# one can compare the gim estimates with those estimated from internal data
fit1 <- glm(form, 'gaussian', dat)
summary(fit1)

# by default, covariates in dat is used as reference in gim
# which assumes that the external and internal populations are the same
fit3 <- gim(form, 'gaussian', dat, model[1], nsample[1, 1], ref = dat)
all(coef(fit) == coef(fit3)) # TRUE

# if additional reference is available, 
# and it comes from the internal population from which dat is sampled
# gim can use it
add.ref <- data.frame(p53 = sample(1:3, 200, TRUE, c(.3, .3, .4)), 
                      subtype = sample(c('A','B','C'), 200, TRUE, 
                                       c(0.6, 0.2, 0.2)), 
                      stringsAsFactors = FALSE)
ref <- rbind(dat[, c('p53', 'subtype')], add.ref)
fit2 <- gim(form, 'gaussian', dat, model[1], nsample[1, 1], ref = ref)

# if the external population is different from the internal population
# then reference for summary data specified in model needs to be provided
ext.ref <- data.frame(p53 = sample(1:3, 200, TRUE, c(.3, .3, .4)), 
                      subtype = sample(c('A','B','C'), 200, TRUE, 
                                       c(0.6, 0.2, 0.2)), 
                      stringsAsFactors = FALSE)
fit3 <- gim(form, 'gaussian', dat, model[1], nsample[1, 1], ref = ext.ref)

}
