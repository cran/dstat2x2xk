\name{dstat2x2xk}
\alias{dstat2x2xk}
\concept{Sensitivity analysis}
\concept{Causal inference}
\concept{Adaptive inference}
\concept{d-statistic}
\concept{Mantel-Haenszel test}
\concept{Contingency table}
\concept{Observational study}
\title{
Sensitivity Analysis in 2x2xK Tables with Demonstrated Insensitivity to Unmeasured Bias
}
\description{
The function dstat2x2xk performs a sensitivity analysis for an observational study that produced a 2x2xk
contingency table recording binary treatment x binary outcome x k strata defined by covariates.  The d-statistic uses only those strata that demonstrate a degree of insensitivity to unmeasured biases in treatment assignment.
}
\usage{
dstat2x2xk(tab, gamma = 1, kappa = NULL, lambda = NULL, rnd = 2, warn = TRUE)
}

\arguments{
  \item{tab}{
A 2x2xK contingency table, with K greater than or equal to two.  That is, tab is a 3-dimensional array of nonnegative integers, with dimensions 2, 2 and K.
}
  \item{gamma}{
Sensitivity parameter.  A number greater than or equal to 1.
}
  \item{kappa}{
An optional parameter.  The d-statistic truncates the extended hypergeometric distribution at kappa times its expectation.  A sensible value of kappa is kappa = 1.  An error will result if both kappa and
lambda are not NULL; give a value to either kappa or lambda, not both.
If lambda and kappa are both NULL, then the result is a sensitivity analysis using the Mantel-Haenszel-Birch test.
}
  \item{lambda}{
An optional parameter.  The d-statistic truncates the extended hypergeometric distribution at its lambda quantile.  Sensible values of lambda are lambda = 0.25 for the lower quartile or lambda = .5 for the median.  A cautious value of lambda is lambda = 0.10. An error will result if both kappa and
lambda are not NULL; give a value to either kappa or lambda, not both.
If lambda and kappa are both NULL, then the result is a sensitivity analysis using the Mantel-Haenszel-Birch test.  Values of
kappa strictly greater than 1 are not recommended: they will generate a warning.
}
  \item{rnd}{
Certain output is rounded to rnd digits for attractive display.
}
  \item{warn}{
If one of the 2x2 subtables is degenerate -- if its conditional
distribution given its margins has only one support point -- then that subtable is removed before the analysis begins.  If warn=TRUE, then
a harmless but informative warning is displayed when one or more subtables is removed.
}
}
\details{
There are (2^K)-1 nontrivial subtables of a 2x2xK table, and the d-statistic considers them all, searching for the subtable least sensitive
to unmeasured bias, paying a limited price for multiple testing.

The method is described in Rosenbaum (2022).  A related method for
continuous responses in matched pairs is described in Rosenbaum (2020)
and is implemented in package \code{\link[dstat]{dstat}}.  Adaptation
between 2 Mantel-Haenszel statistics, rather than (2^K)-1 statistics, is discussed in Rosenbaum and Small (2017) and implemented in \code{\link[sensitivity2x2xk]{adaptmh}}.
}
\value{
\item{tstat }{Value of the test statistic.  See below.}
\item{pval }{Upper bound on the on-sided P-value testing the null
hypothesis of no treatment effect in the presence of a bias in
treatment assignment of magnitude at most gamma.}
\item{detail }{A dataframe showing detailed results, stratum by
stratum.  If the strata have dimnames, then the rows of detail
are labeled with this names.  Column n11 is the observed count in
the 11 upper left corner cell of a 2x2 subtable of the 2x2xK table.
Column E(n11) is the null expectation of n11 in the upper bounding
extended hypergeometric distribution with parameter gamma.  Column
cut is the lower truncation point for n11, and column use is 1
if n11 >= cut; otherwise use=0.  The test statistic, tstat,
is the sum of n11 for strata with use=1, and columns
E(n11|use=1) and var(n11|use=1) give moments of its conditional
distribution given that it is used.  Column OR is simply the sample odds
ratio for a 2x2 subtable.  In the notation of Rosenbaum (2022), the column labeled cut is \eqn{t_{s}} and the column labeled use is \eqn{A_{s}}.}
}

\references{

Birch, M. W. (1964) <doi:10.1111/j.2517-6161.1964.tb00564.x> The detection of partial association, I: the 2 x 2 case. Journal of the Royal Statistical Society: Series B (Methodological), 26(2), 313-324.

Hiller, R., Giacometti, L., and Yuen, K. (1977)
<doi:10.1093/oxfordjournals.aje.a112404> Sunlight and cataract: an epidemiologic investigation. American Journal of Epidemiology, 105(5), 450-459.

Mantel, N., and Haenszel, W. (1959) <doi:10.1093/jnci/22.4.719> Statistical aspects of the analysis of data from retrospective studies of disease." Journal of the National Cancer Institute 22(4): 719-748.

Rosenbaum, P. R. and Small, D. S. (2015) <doi:10.1111/biom.12591> An adaptive Mantel-Haenszel test for sensitivity analysis in observational studies.  Biometrics 73, no. 2 (2017): 422-430.

Rosenbaum, P. R. (2020) <doi:10.1093/biomet/asaa032> A conditional test with demonstrated insensitivity to unmeasured bias in matched observational studies. Biometrika, 107(4): 827–840.

Rosenbaum, P. R. (2022) A Mantel-Haenszel-Birch statistic with demonstrated insensitivity to unmeasured bias in observational studies.
Manuscript.
}
\author{
Paul R. Rosenbaum
}
\note{
The concept of a d-statistic is introduced in Rosenbaum (2020) for matched pairs with continuous responses, and this version of a d-statistic is implemented in the package "dstat".  Rosenbaum (2022) extends the idea
to contingency tables, and dstat2x2xk implements this extension.
}




\examples{
# Data from Hiller et al. (1977)

cataract<-c(33,8,356,65,23,6,222,35,
            139,19,477,82,114,30,451,45,
            76,33,210,23,99,26,282,27,
            172,121,133,25,364,165,237,32)

cataract<-array(cataract,c(2,2,8))
dimnames(cataract)<-list(light=c("less","more"),
  status=c("case","referent"),
  str=c("M20-44","F20-44","M45-64","F45-64",
  "M65-74","F65-74","M75+","F75+")
)
cataract2<-cataract[,2:1,]

# The following example compares the d-statistic and the usual
# Mantel-Haenszel-Birch exact statistic.
# With kappa=1, the d-statistic truncates at the expectation
# of the null distribution, eliminating three strata.
# With kappa=0, the d-statistic does no truncation and
# equals the Mantel-Haenszel-Birch exact statistic.
# The d-statistic reports greater insensitivity to unmeasured bias in
# this example.

# Truncation at the expectation
dstat2x2xk(cataract2,gamma=2.365,kappa=1)


# Truncation at the lower quartile and median
dstat2x2xk(cataract2,gamma=2.49,lambda=.25)
dstat2x2xk(cataract2,gamma=2.50,lambda=.25)

# No truncation = exact, usual Mantel-Haenszel-Birch (MHB) test
dstat2x2xk(cataract2,gamma=2.365,kappa=0)
# The same usual untruncated MHB answer is
# produced by the mh() function in package sensitivity2x2xk
# mh(cataract2,Gamma=2.365)
}

\keyword{ htest }

