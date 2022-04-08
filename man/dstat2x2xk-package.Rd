\name{dstat2x2xk-package}
\alias{dstat2x2xk-package}
\docType{package}
\title{
\packageTitle{dstat2x2xk}
}
\description{
\packageDescription{dstat2x2xk}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{dstat2x2xk}
\packageIndices{dstat2x2xk}
There is one function, dstat2x2xk, intended to increase
design sensitivity when analyzing a 2x2xk table,
treatment x binary outcome x strata, from an
observational study.  See Rosenbaum (2022).
}
\author{
\packageAuthor{dstat2x2xk}

Maintainer: \packageMaintainer{dstat2x2xk}
}
\references{
Rosenbaum, P. R. (2020) <doi:10.1093/biomet/asaa032> A conditional test with demonstrated insensitivity to unmeasured bias in matched observational studies. Biometrika, 107(4): 827–840.

Rosenbaum, P. R. (2022) A Mantel-Haenszel-Birch statistic with demonstrated insensitivity to unmeasured bias in observational studies.
Manuscript.
}

\keyword{ package }

\examples{
cataract<-c(33,8,356,65,23,6,222,35,
            139,19,477,82,114,30,451,45,
            76,33,210,23,99,26,282,27,
            172,121,133,25,364,165,237,32)

cataract<-array(cataract,c(2,2,8))
dimnames(cataract)<-list(light=c("less","more"),
  status=c("case","referent"),
  str=c("M20-44","F20-44","M45-64","F45-64","M65-74","F65-74","M75+","F75+")
)
cataract2<-cataract[,2:1,]

dstat2x2xk(cataract2,gamma=2.365,kappa=1)
}
