# Change log of the R package 'layeranalyzer'

## layeranalyzer 0.4.1 (2026-07-09)

Merging of code by Trond Reitan and Adam Kocsis. Adam has been working on making the code CRAN-compliant. Trond has updated layeranalyzer with internal IO logic. 

## layeranalyzer 0.4.0 (2026-06-19)

### 11/6-2026

Extra sanity checks of covariance matrices (positive definiteness, checked on eigenvalues) for laxness levels "moderate" and "low". Sanity check for laxness level "high" (default) has been sped up. Also, sanity check is now performed on `S_k`, the measurement covariance matrix. Changed from using `mkl_set_num_threads` to `flexiblas_set_num_threads` on Linux, because the mkl version crashed on Fedora 44.

### 11/3-2026

Added an alternative way of specifyng connections (used instead of 'causal' and 'corr') if one wants a compact (but less user-friendly) way to specify connections. This alternative way of specifying a connection was introduced in order to make it easy to produce different connection models by simply traversing the possible connection between each pair of processes. If one use a script to traverse all possible connections models, reads a script and corrects the connection specification, one has a tool for outputing one script for each connection models. These scripts can then possibly be run in parallel on a computer with many cores or a computing cluster.

In addition, support for distance-based correlation (exponential decrease of correlation with distance) has been introduced. One needs to specify a distance matrix between the different sites. This for using layeranalyzer to analyze spatiotemporal processes. Note however that there is much more freedom in specifying the temporal correlation structure than the spatial correlation structure.

### 11/02-2026

Added support for external time series. That is series that are not treated as stochastic processes to be modelled, but only considered as possible causal drivers of the observational time series that are modelled. In order for this to work properly, they should have many times more sample points than the observational time series (and also going further back in time), as their influence on the observed time series are calculated using numerical integration. Consider interpolation algorithms if that is not the case. This is mainly meant for time series that do not represent directly observed processes, but instead come about by modelling, interpolation and possible smoothing (such as paleontological climate series).


## layeranalyzer 0.3.4 

### 17/10-2025

Memory leak issue in log-likelihood calculations due to early return in error situations fixed and tested. This was the version, 0.3.4, found in `archive_/layeranalyzer_0.3.4.tar.gz`.


## layeranalyzer 0.3.3

### 02/10-2025

Combined new code (for using stationary standard deviation instead of stochastic contribution for stationary processes) with CRAN compliance corrections, to make Version 0.3.3,.


## layeranalyzer 0.3.1

### 17/06-2025

Big functionality upgrade. Possible to start ML optimization from user-specified parameter values in 'layer.analyzer.timeseries.list' or to simply ask for number of parameters and parameter names. This has to be considered "expert mode" susage. Simpler version avaiable in 'layer.analyzer', where one can start optimization for a connection model using parameter values from a non-connected model (run on the fly). Also added option for consistency checks ("laxness" options). Some general debugging, including changing from one optimization routine ("L-BFGS-B") to another (Nelder-Mead). (May later consider to let choice of optimization strategy to be user-specified.) Version 0.3.1 is now default, found in `layeranalyzer_0.3.1.tar.gz`.

## layeranalyzer 0.3.0

### 07/04-2025

Massive debug effort.

## layeranalyzer 0.2.1

### 01/4-2024

Reverting to old lapack with new call structure.. Version 0.2.1 is now default, found in layeranalyzer_0.2.1.tar.gz.

## layeranalyzer 0.2.0

### 06/10-2023

Support for new Lapack structure added. Found in layeranalyzer_0.2.0.tar.gz.

## layeranalyzer 0.1.1

### 10/06-2022

New version, 0.1.1 is now default. Found in layeranalyzer_0.1.1.tar.gz.

### 15/06-2019

Vignettes added.

## layeranalyzer 0.1.0

### 05/01-2019

R package added: layeranalyzer_0.1.0.tar.gz.

### 23/11-2017

New version of the analysis program simply called "layeranalyzer" added. This 
  version does not depend on anything but the Lapack library. All parts of the Hydrasub library 
  that were needed, has been included in the program code itself. This should make compilation 
  easier and allow for compilation on more platforms. The change was made in order to prepare 
  the groun for an R package version.
  
### 25/2-2015

Debug on layer_analyzer_complex has been performed on the option for linear time 
  trends and the option has been tested on simulated files. The input for multiple files has 
  been shortened down to avoid using having dual ways of specifying the same options. Run 
  layer_analyzer_complex without input arguments to see the changes.
  
### 21/3-2013

The program (layer_analyzer_complex) can now subtract a periodic pattern 
  (sine/cosine) from the observations, thus allowing seasonality or other strict periodic 
  phenomena to be removed. The length of the period is specified. Multiple periodicities 
  can be specified, making ti possible to handle arbitrarily complex period patterns. The 
  coefficients in fron of the sine/cosine terms enter in as parameters. In addition, some 
  extra debug info was added as options. Also the parallel tempering was given an extra 
  option, where the multiplicative steps of the 'temperature' now can be specified. Also, 
  starting parameters for the MCMC chain can now be set. A short intro to using the program 
  was made, here. The into does not show all the options. Call the program without input in 
  order to see that.
  
### 21/5-2012

The program layer_analyzer_complex is now able to deal with multiple correlated 
  series, using an extra file specifying the correlation in the observations. Thus allometry 
  in biology is dealt with in the model. Correlations as between series as well as causal 
  pathways are possible. This machinery is currently being tested on microfossil data. The 
  results so far look encouraging.
  
### 23/3-2012

A large expansion that will later replace the current layer_analyzer is now made 
  available as a separate program named layer_analyzer_complex. The program deals with 
  separately specified feedback loops and also with the resulting problem of complex eigenvalues.
  
### 7/2-2012

New code for dealing with pair-wise inter-regional correlations in a specified layer. 
  (Note that there should be no more than 6 or at worst 7 sites for the methodology to be efficient. 
  If there's more, the resources to the initial Monte Carlo estimation of the probability of 
  drawing positive covariance structures needs to be adjusted in the code. Either that or a 
  completely different prior for the correlation structure, like the inverse-Wishart must be 
  implemented.) Also, an indicator on/off option for allowing the analysis to switch off 
  inter-regional correlations for some sites.
  
### 10/2-2012

Bug fix for 'layer_analyzer' concerning the prior identification strategies in 
  combination with itnegrated layers. Also, test code for importance sampling added.
  
### 17/1-2012

Bug fix on 'layer_analyzer' on the priors for hierarchical models with identification 
  restriction. Possibility of switching off identification restrictions. Realizations from 
  smoothing estimates as well as simulations based on given parameters also implemented as options. 
  Possibility of plugging in secondary data sources that will also be modelled using linear SDEs. 
  The point is to make 'layer_analyzer' a single tool for all practical linear SDE problems. 

### Older code:

- layer_analyzer_complex.C, layer_analyzer_complex.H: As layer_analyzer but is able to deal with 
  feedback loops and complex eigenvalues (thus making this a tractable way of continuous time 
  ecological modelling). May run slower on such problems, but should run just as fast as usual 
  on models without feedback loops. Also able to deal with non-stationary (exploding) time series, 
  as long as an initial value treatment is specified. Also, multiple measurements at the same time 
  possble (as long as they don't belong to the same site and the same series). Thus phenotypical 
  traits with allometry can be handled. Replaced the old "layer_analyzer" as the default analysis 
  program, but has now been replaced by the program "layeranalyzer" (see above, and sorry for the 
  lack of imagination), which doesn not use the Hydrasub library and has replaced the use of the 
  GSL library for the more universally available "lapack" library.
  
- layer_analyzer.C, layer_analyzer.H: A generalized layered analyzer, though has been replaced by an 
  even more generalized program (now called layer_analyzer_complex), which allows for feedback 
  loops (something this version does not). It does Bayesian analysis (and also a subsequent ML 
  analysis, if that is wanted) on a layered stochastic differential equation set. More than one 
  time series can be analyzed at a time, with options for letting a layer belonging to one series 
  affect a layer belonging to another.
  
- many_mcmc_layer: A Perl-script for traversing a host of models (710 in number) in order to find 
  the posterior model probabilities and find the most probable model. Uses 'layer_analyzer'.
  
- supplementary.zip: A bundle of the supplementary text, the dataset and the main source codes 
  (layer_analyzer). 
