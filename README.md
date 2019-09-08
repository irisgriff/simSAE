simSAE
================

Introduction
============

The U.S. Forest Service's Forest Inventory and Analysis (FIA) program conducts its forest inventory using a sample of ground plots which are surveyed every ten years. Researchers with the FIA are often required to estimate forest attributes for small geographic domains of interest, and in some cases, these small domains include very few to zero sample plots. The field of small area estimation (SAE) is concerned with the task of constructing reliable estimates for responses within these so called "small areas". As a result, the challenges associated with SAE calls for unique estimation methods that perform better in the small area setting.

There are three classes of estimators available in the SAE literature: Direct, indirect, and small area model-based estimators (Rao 2005). Direct estimators only use the observed data for both domain estimation and estimation of respective MSEs. In the context of FIA's work, a direct estimator uses only sample plot level data to generate an estimate. In contrast to direct estimators, indirect estimators are a class of model assisted estimators that use auxiliary population level data to generate estimates, while borrowing strength from geographically similar observations outside of the domain of interest. In particular, small area model-based estimators are indirect estimators that incorporate a random effects model to bolster the process of gathering strength. More specifically, these notes aim to compare the performances of the direct simple random sampling estimator (SRS), the indirect synthetic regression estimator (SRE), the generalized regression estimator (GREG), and the Empirical Best Linear Unbiased Prediction (EBLUP), a small area model-based estimator. Much of the work outlined in these notes draw from the 2012 Breidenbach and Astrup paper, which describes these four estimators as the following:

![SRS\_i = \\bar{y\_i}](https://latex.codecogs.com/png.latex?SRS_i%20%3D%20%5Cbar%7By_i%7D "SRS_i = \bar{y_i}")

![SRE\_i = \\bar X\_i^T\\hat \\beta\_{OLS}](https://latex.codecogs.com/png.latex?SRE_i%20%3D%20%5Cbar%20X_i%5ET%5Chat%20%5Cbeta_%7BOLS%7D "SRE_i = \bar X_i^T\hat \beta_{OLS}")

![GREG\_i = \\bar X\_i^T\\hat \\beta\_{OLS} + \\frac{1}{n\_i}\\sum\_{j \\in S\_i} \\varepsilon\_j](https://latex.codecogs.com/png.latex?GREG_i%20%3D%20%5Cbar%20X_i%5ET%5Chat%20%5Cbeta_%7BOLS%7D%20%2B%20%5Cfrac%7B1%7D%7Bn_i%7D%5Csum_%7Bj%20%5Cin%20S_i%7D%20%5Cvarepsilon_j "GREG_i = \bar X_i^T\hat \beta_{OLS} + \frac{1}{n_i}\sum_{j \in S_i} \varepsilon_j")

![EBLUP\_i = \\bar X\_i^T\\hat \\beta\_{LMM} + \\frac{\\gamma\_i}{n\_i}\\sum\_{j \\in S\_i} \\epsilon\_j](https://latex.codecogs.com/png.latex?EBLUP_i%20%3D%20%5Cbar%20X_i%5ET%5Chat%20%5Cbeta_%7BLMM%7D%20%2B%20%5Cfrac%7B%5Cgamma_i%7D%7Bn_i%7D%5Csum_%7Bj%20%5Cin%20S_i%7D%20%5Cepsilon_j "EBLUP_i = \bar X_i^T\hat \beta_{LMM} + \frac{\gamma_i}{n_i}\sum_{j \in S_i} \epsilon_j")

Where ![\\bar X^T](https://latex.codecogs.com/png.latex?%5Cbar%20X%5ET "\bar X^T") is a design matrix of population level domain means of the predictor variables, and ![\\gamma\_i = \\frac{\\sigma\_\\nu ^2}{\\sigma\_\\nu ^2 + \\sigma\_\\epsilon ^2/n\_i}](https://latex.codecogs.com/png.latex?%5Cgamma_i%20%3D%20%5Cfrac%7B%5Csigma_%5Cnu%20%5E2%7D%7B%5Csigma_%5Cnu%20%5E2%20%2B%20%5Csigma_%5Cepsilon%20%5E2%2Fn_i%7D "\gamma_i = \frac{\sigma_\nu ^2}{\sigma_\nu ^2 + \sigma_\epsilon ^2/n_i}") is a mixture parameter that determines how much bias correction is introduced, and ![\\sigma\_\\nu^2](https://latex.codecogs.com/png.latex?%5Csigma_%5Cnu%5E2 "\sigma_\nu^2") is the variance of the random intercept of the EBLUP's linear mixed model.

Our work aims to streamline the small area estimation process by determining a set of best practices for deciding when to use which estimator. In order to acomplish outlining this set of guideline, we devised a simulation study to test the empirical MSEs of these four estimators. The simulation is constructed using the plot level fire data in the states of Montana and Idaho. In order to obtain "truth" or population values to compare with the results of the various estimation techniques, we considered all plots within a given domain with size 100 or greater as our "population". Then, our simulation would select a small random sample of plots within each domain (of a random size proportional to the size of the domain) from this "population" which would then act as our sample for that iteration. From there, estimates of our response variable (tree basal area) were produced using the four different estimators. For the synthetic and GREG estimators, an OLS model was fit to predict tree basal area using the predictors `tcc` and `bio`, which represent tree canopy cover and biomass, respectively. And for the EBLUP estimator, a linear mixed effects model was fit using the same variables, but introduced random intercept terms to control for domain specific effects.

References
==========

Breidenbach, Johannes, and Rasmus Astrup. 2012. “Small Area Estimation of Forest Attributes in the Norwegian National Forest Inventory.” *European Journal of Forest Research* 131 (4). Springer Science; Business Media LLC: 1255–67. doi:[10.1007/s10342-012-0596-7](https://doi.org/10.1007/s10342-012-0596-7).

Rao, J.N.K. 2005. *Small Area Estimation*. Wiley Series in Survey Methodology. Wiley. <https://books.google.com/books?id=f8NY6M-5EEwC>.
