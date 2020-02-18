\name{gim-package}
\alias{gim-package}
\docType{package}
\title{
Generalized Integration Model
}
\description{
\code{gim} implements the generalized integration model proposed in Zhang et al. (2018). \code{gim} integrates individual-level data and summary statistics under a generalized linear model framework. It supports continuous and binary outcomes to be modeled by the linear and logistic regression models. For binary outcome, data can be sampled in prospective cohort studies or case-control studies. \code{gim} uses \bold{\emph{different}} methods to model binary outcome under \bold{\emph{different}} designs. 
}
\details{
\tabular{ll}{
Package: \tab gim\cr
Type: \tab Package\cr
Version: \tab 0.17.0\cr
Date: \tab 2019-02-25\cr
License: \tab MIT + file LICENSE\cr
}

Meta-analysis has become a powerful tool for enhanced inference by gathering evidence from multiple sources. It pools summary-level data, i.e. estiamtes of model coefficients, from different studies to improve estimating efficiency under the assumption that all participating studies are analyzed under the same statistical model. This assumption, however, is usually not true in practice as studies may adjust for different covariates according to specific purpose, or are conduct on partial observed covariates they collect. Meta-analysis can lead to biased estimates when this assumption is violated. 

It is challenging to integrate external summary data calculated from different models with a newly conducted internal study in which individual-level data is collected. \code{gim} is a novel statistical inference framework based on the \bold{G}eneralized \bold{I}ntegration \bold{M}odel, which effectively synthesizes internal and external information according to their variations for multivariate analysis. This new framework is versatile to incorporate various types of summary data from multiple sources. It can be showed that the \code{gim} estimate is theoretically more efficient than the internal data based maximum likelihood estimate, and the recently developed constraint maximum likelihood estimate that incorporates the outside information. 

The \code{gim} function implemented in this package accounts for the sample sizes shared by different studies. Ignoring this sample overlap may lead to inflated false positive. \code{gim} requires estimates of coefficients in external working models, but do not rely on their standard errors for two reasons. (1) It is more convenient to request less information from users, especially when the unrequired information could be estimated with other given information. (2) The standard errors reported in literatures where users collect their external data maybe underestimated as a working model rather than a true underlying model was assumed. 

The \code{gim} always requests a set of raw data in which both outcome and independent variables are available. This dataset is called the reference set or internal data. This requirement seems unconvenient as outcome can sometimes be expensive, and some other approaches may be applicable with a reference where only independent variables are collected. Theoretically, an slightly extended \code{gim} can become workable under the same circumstances, however, we suggest to be more careful to do so. A reference is used to estimate correlation between variables (including outcome) in the population of interest. In practice, there could be a difference in the population of your own study and external studies (from which summary information are collected), therefore, the correlation differ among studies. Conducting an analysis soly relying on a set of independent variables may be biased if the difference between studies is significant. In comparison, a full set of reference consisting of outcome can make the analysis more robust to the potential population difference.

The main function in this package is \code{\link{gim}}. 

}
\author{
Han Zhang, Kai Yu

Maintainer: Han Zhang <han.zhang2@nih.gov>
}
\references{
Zhang, H., Deng, L., Schiffman, M., Qin, J., Yu, K. (2018) Generalized integration model for improved statistical inference by leveraging external summary data. Biometrika. Accepted. 
}

\keyword{Generalized Integration Model}

