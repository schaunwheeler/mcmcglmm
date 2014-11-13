# mcmcglmm

*PLEASE NOTE: THIS PACKAGE WAS ONLY EVER EXPERIMENTAL. IT ONLY EVER WORKED
WITH THE "GAUSSIAN" FAMILY. IT WAS NEVER EXTENSIVELY TESTED. AND I NO
LONGER ACTIVELY MAINTAIN THIS CODE. I DON'T EVEN WRITE R CODE ANYMORE.*

The aim of `mcmcglmm` is to provide some added functionality to the MCMCglmm
package by facilitating cross validation through both data preparation and 
new-data prediction, by implementing default priors, and by providing a few 
tools to evaluate parameters and model fit. You can track and contribute to 
development of `mcmcglmm` at https://github.com/schaunwheeler/mcmcglmm.

## Package (little) mcmcglmm

Data Preparation:

* `SplitData` takes a data frame and splits it into a large subset, to be used
  for model training, and a small subset, to be used for cross validation. The
  function checks the small subset to make sure it does not contain variable 
  options not included in the large subset, thus ensuring that cross-valiation
  checks will be possible (since it's not very easy to predict based on
  variables that weren't included in the original model).

Bayesian Modeling:
  
* `mcmcglmm` is a wrapper for the MCMCglmm() function in the `MCMCglmm` package
  developed by Jerrod Hadfield (http://cran.r-project.org/web/packages/MCMCglmm/
  index.html). The wrapper function allows for two variants of two defualt
  priors on the covariance matrices. The two defaults are InvW for an inverse-
  Wishart prior, which sets the degrees of freedom parameter equal to the 
  dimension of each covariance matrix, and InvG for an inverse-Gamma prior, 
  which sets the degrees of freedom parameter to 0.002 more than one less than 
  the dimensions of the covariance matrix. "-pe" can be added to the call for 
  either of these priors to use parameter-expanded priors. The function also
  saves the levels for each variables stored as a character or factor, which
  faclitates the PredictNew() function. Unlike the MCMCglmm function, mcmcglmm
  saves the random effects values as a default.

Model Evaluation:

* `QuickSummary` provides some slightly-more-than-basic measures for evaluating
  an 'mcmcglmm' output. Given the output, the function calculates the posterior
  mean, the highest posterior density intervals for a given probability (set
  through the "prob" option), the "type S" error (probability that the estimate
  actually is of the opposite sign of the posterior mean), and the "type M" error 
  (probability that the estimate is the same sign but substantially smaller than
  the posterior mean - defaults to measuring the probability that the estimate
  is less than one half the size of the mean). The function also allows for 
  rounding of the output for convenience - defaults to four decimal places.
  
* `PredictNew` is a modified version of predict.MCMCglmm() that allows for 
  prediction based on new data. Additionally, PredictNew() differs from 
  predict.MCMCglmm() in that by default it marginalizes none of the random
  effects, whereas predict.MCMCglmm() marginalizes all random effects by 
  default, and defaults to predicting values in post-link-function (Gaussian)
  scale, rather than on their original scale. This function also allows a 
  vector to be passed to an "index" option, which will append that index to
  the output.
