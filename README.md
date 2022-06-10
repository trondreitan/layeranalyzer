# layeranalyzer
R package for causal and correlative time series analysis with hidden layers, using linear stochastic differential equations

In addition to the R package itself, this page presents the C++ code (used in a standalone program as well as the R package), documentation and example datasets and code for multi-layered analysis using vectorial linear stochastic differential equations. The reason for using stochastic differential equation is to facilitate studying irregularly spaced observations like those we find in our study dataset of micro-fossils. This work was performed in conjuction with doing the research for and writing the manuscript "Phenotypic Evolution studied by Layered Stochastic Differential Equations" by T. Reitan, T. Schweder and J. Henderiks.

Originally, this was used for modelling phenotypic evolution towards a optima that itself could be a stochastic process depending responding to layers below it (climate or primary optimum). However, the method have be used for other applications as well.


Code:

* The R package, layeranalyzer_0.1.1.gz, uses the same underlying C++ code (layeranalyzer.cpp). It can be installed with the R code: install.packages("https://github.com/trondreitan/layeranalyzer/raw/master/layeranalyzer_0.1.1.tar.gz",type="source") or install_github(repo="trondreitan/layeranalyzer",dependencies=FALSE,build_vignettes=TRUE) if devtools is installed. Note that if vignettes are activated, you need the rmarkdown and markdown packages also! On Linux at least, you need the program 'pandoc' also (not in R but on the Linux machine). You might need sudo rights for that. If this is too troublesome, use the option 'vignettes=FALSE' instead.  See "troubleshooting" for technical issues. 

* Note that a new install file for upcoming version has been added. The old one can be installed with install.packages("https://github.com/trondreitan/layeranalyzer/raw/master/layeranalyzer_0.1.0.tar.gz",type="source"). The new version already has a few debug changes due to changes in the "class" function and one due to maximum likelihood calculations that sometimes failed in a certain context. I have also added code for doing prediction as a separate task from parameter inference ("calibration"). However, the R interface hasn't been finished and then there is the testing, so more will come.

* layeranalyzer.cpp. This is the current source version of the program. This version does not depend on any library but the Lapack library (for linear algebra operations), which typically can be found in an R installation. The program can then typically be compiled on a Linux machine like this: "g++ -DMAIN -DENGLISH_LANGUAGE -I/usr/include/R -I/usr/include/R/R_ext -o layeranalyzer layeranalyzer.cpp -lm -llapack". If the Lapack library files or header files are somewhere else on your computer, you should change the "-I" (which tells the compiler where to look for header files) and "-L" (which tells the compiler where to search for library files other than in the /usr/lib directory). This has already been tested on multiple platforms. See "help" texts for more on usage. See below in this text for some examples. 



Documentation:

Presentations:

* BISP_small.pdf/pptx is a presentation made for the Seventh Conference on Bayesian Inference in Stochastic Processes in Madrid, Spain, BISP7.
* stoch_layers7.pdf/pptx is a simpler presentation made for CEES, UiO of the basic modelling framework and how it connects to the application can be found here.
* layer_analyzer_complex.pptx is a how-to concerning the program usage (an older version called 'layer_analyzer_complex' is used, but it has the same input structure as the current 'layeranalyzer' program), with a short intro to continuous time modelling. 
* linkoping2021.pptx - A presentation I did over Zoom for the statistics department at the university of Linköping. Focuses slightly more on the mathematics than the other presentations.
    
Papers on theory and application:

    Reitan, T., Schweder, T., Henderiks, J. (2012)
    Phenotypic Evolution studied by Layered Stochastic Differential Equations
    Annals of Applied Statistics, Volume 6 (4): 1531-1551.
    DOI: 10.1214/12-AOAS559
    Link: https://projecteuclid.org/euclid.aoas/1356629050.
    supplementary.pdf: A supplementary to the article 
      "Phenotypic Evolution studied by Layered Stochastic Differential Equations", 
      which goes more into detail concerning the algorithmical aspects of the inference.

    Liow, L.H., Reitan, T., Harnik, P.G. (2015)
    Ecological interactions on macroevolutionary time scales; clams and brachiopods are more than ships that pass in the night
    Ecology Letters, Volume 18(10): 1030-1039.
    DOI: 10.1111/ele.12485
    Link: http://onlinelibrary.wiley.com/doi/10.1111/ele.12485/full.

    Reitan, T., Liow, L. H. (2017)
    An unknown Phanerozoic driver of brachiopod extinction rates unveiled by multivariate linear stochastic differential equations.
    Paleobiology, Published online.
    DOI: 10.1017/pab.2017.11
    Link: https://www.cambridge.org/core/journals/paleobiology/article/an-unknown-phanerozoic-driver-of-brachiopod-extinction-rates-unveiled-by-multivariate-linear-stochastic-differential-equations/0EAB1A42D58CCD17012D9CF3020D57EA.    

Examples
R example code and vignettes as html files:

    analyze_malta.R/analyze_malta_data.html: Analyzis of reed warbler body size, provided by 
      Camilla Lo Cascio Sætre, published in Nature Communications.
    
    normalize_harelynx.R/analyze_normalized_harelynx.R/analyze_normalized_harelynx.html: 
      Analyzis for hare/lynx data from the Hudson Bay Company. The data is first normalized (using 
      normalize_harelynx.R) before analysis is performed. This script, normalize_harelynx.R, is not 
      strictly necessary to run, for the analyze_normalized_harelynx.R script, since the normalized 
      datasets ('hare.norm' and 'lynx.norm') are provided in the package. However, the transformation 
      function is defined in the first sccript and is used for plotting in the second script.The code 
      does a brief foray into the process structure of each time series separately, before doing a 
      "simple" connection analysis with the simplest process structure available. Some residual 
      analysis is then performed, which reveals that changes in the process structure may be 
      warranted. The process structure for lynx is then expanded and a new connection analysis is 
      performed.
      
    brachbivalve_2layer_local_id4.R/brachbivalve_2layer_local_id4.html: Analysis of bivalve+brachiopod 
      diversifications rates. This analysis seeks to check the model deemed best in Reitan & Liow (2017), 
      using the R package. Some extra residual analysis is also done. 

PS: All example datasets are provided in the package, so that the file reading is not necessary.

Example datasets (see in addition the example code):

    full_sd_smoothed.txt: This is the aggregated dataset with smoothed standard deviations analyzed 
      in the paper "Phenotypic Evolution studied by Layered Stochastic Differential Equations", and 
      was collected by Jorijntje Henderiks. full_sd_smoothed_prior2.txt contains the prior 
      specification, used in layered_analyzer.
      
    default_prior.txt. Defines the default prior used in the R package.
    
    test_1layer.txt, test_2layer.txt, test_3layer.txt, test_cause.txt+test_effect.txt, 
      test_corr1.txt+test_corr2.txt, test_twoway1.txt+test_twoway2.txt. Various similated test files.
      
    temp.txt: Three series of water temperature from the Norwegian Water Resources Administration, NVE. 
      Some data has been removed to check the state inference. The full data is found in the file temp.txt 
      which can be used for comparison in layered_analyzer. Prior file: temp_prior.txt . 

Simulations:
Simulations can be used for testing the behaviour of the analysis.

    numeric_stoch_diff.R. Contains a general purpose stochastic differential equation 
      simulations tool called explicit.stoch.
    
    generate_1_2_3_layers.R. Uses the simulations tool above to generate 1-, 2- and 3-layered 
      linear SDEs and were used for generating the test files test_1layer.txt, test_2layer.txt and 
      test_3layer.txt above.
    
    generate_causal.R. Generates two processes, one OU process (the cause) and one process 
      causally tracking the OU process (the effect). Used for generating the test files t
      est_cause.txt and test_effect.txt.
    
    generate_corr.R. Generates two OU processes having a correlative connection. Used for 
      generating the test files test_corr1.txt and test_corr2.txt.
    
    generate_twoway_causal.R. Generates two processes having a positive-negative feedback loop. 
      Used for generating the test files test_twoway1.txt and test_twoway2.txt. 

Updates and older versions
Update history:

    10/06-2022: New version, 1.1 is now default. Found in layeranalyzer_0.1.1.tar.gz.

    15/06-2019: Vignettes added.

    05/01-2019: R package added: layeranalyzer_0.1.0.tar.gz.
    
    23/11-2017: New version of the analysis program simply called "layeranalyzer" added. This 
      version does not depend on anything but the Lapack library. All parts of the Hydrasub library 
      that were needed, has been included in the program code itself. This should make compilation 
      easier and allow for compilation on more platforms. The change was made in order to prepare 
      the groun for an R package version.
      
    25/2-2015: Debug on layer_analyzer_complex has been performed on the option for linear time 
      trends and the option has been tested on simulated files. The input for multiple files has 
      been shortened down to avoid using having dual ways of specifying the same options. Run 
      layer_analyzer_complex without input arguments to see the changes.
      
    21/3-2013: The program (layer_analyzer_complex) can now subtract a periodic pattern 
      (sine/cosine) from the observations, thus allowing seasonality or other strict periodic 
      phenomena to be removed. The length of the period is specified. Multiple periodicities 
      can be specified, making ti possible to handle arbitrarily complex period patterns. The 
      coefficients in fron of the sine/cosine terms enter in as parameters. In addition, some 
      extra debug info was added as options. Also the parallel tempering was given an extra 
      option, where the multiplicative steps of the 'temperature' now can be specified. Also, 
      starting parameters for the MCMC chain can now be set. A short intro to using the program 
      was made, here. The into does not show all the options. Call the program without input in 
      order to see that.
      
    21/5-2012: The program layer_analyzer_complex is now able to deal with multiple correlated 
      series, using an extra file specifying the correlation in the observations. Thus allometry 
      in biology is dealt with in the model. Correlations as between series as well as causal 
      pathways are possible. This machinery is currently being tested on microfossil data. The 
      results so far look encouraging.
      
    23/3-2012: A large expansion that will later replace the current layer_analyzer is now made 
      available as a separate program named layer_analyzer_complex. The program deals with 
      separately specified feedback loops and also with the resulting problem of complex eigenvalues.
      
    7/2-2012: New code for dealing with pair-wise inter-regional correlations in a specified layer. 
      (Note that there should be no more than 6 or at worst 7 sites for the methodology to be efficient. 
      If there's more, the resources to the initial Monte Carlo estimation of the probability of 
      drawing positive covariance structures needs to be adjusted in the code. Either that or a 
      completely different prior for the correlation structure, like the inverse-Wishart must be 
      implemented.) Also, an indicator on/off option for allowing the analysis to switch off 
      inter-regional correlations for some sites.
      
    10/2-2012: Bug fix for 'layer_analyzer' concerning the prior identification strategies in 
      combination with itnegrated layers. Also, test code for importance sampling added.
      
    17/1-2012: Bug fix on 'layer_analyzer' on the priors for hierarchical models with identification 
      restriction. Possibility of switching off identification restrictions. Realizations from 
      smoothing estimates as well as simulations based on given parameters also implemented as options. 
      Possibility of plugging in secondary data sources that will also be modelled using linear SDEs. 
      The point is to make 'layer_analyzer' a single tool for all practical linear SDE problems. 

Older code:

    layer_analyzer_complex.C, layer_analyzer_complex.H: As layer_analyzer but is able to deal with 
      feedback loops and complex eigenvalues (thus making this a tractable way of continuous time 
      ecological modelling). May run slower on such problems, but should run just as fast as usual 
      on models without feedback loops. Also able to deal with non-stationary (exploding) time series, 
      as long as an initial value treatment is specified. Also, multiple measurements at the same time 
      possble (as long as they don't belong to the same site and the same series). Thus phenotypical 
      traits with allometry can be handled. Replaced the old "layer_analyzer" as the default analysis 
      program, but has now been replaced by the program "layeranalyzer" (see above, and sorry for the 
      lack of imagination), which doesn not use the Hydrasub library and has replaced the use of the 
      GSL library for the more universally available "lapack" library.
      
    layer_analyzer.C, layer_analyzer.H: A generalized layered analyzer, though has been replaced by an 
      even more generalized program (now called layer_analyzer_complex), which allows for feedback 
      loops (something this version does not). It does Bayesian analysis (and also a subsequent ML 
      analysis, if that is wanted) on a layered stochastic differential equation set. More than one 
      time series can be analyzed at a time, with options for letting a layer belonging to one series 
      affect a layer belonging to another.
      
    many_mcmc_layer: A Perl-script for traversing a host of models (710 in number) in order to find 
      the posterior model probabilities and find the most probable model. Uses 'layer_analyzer'.
      
    supplementary.zip: A bundle of the supplementary text, the dataset and the main source codes 
      (layer_analyzer). 

The latest version of the analysis program and R package uses only the Lapack library. Older versions use the Hydrasub library for support. This library in turn uses GSL (Gnu Scientific Library). Plotting programs, such as 'vvgraph' and 'histogramme' are used for displaying the results. PS: If plotting in the program (rather than the R package) is to be used, then even the latest version requires the 'vvgraph' and 'histogramme' programs from the Hydrasub library. However, there are options for not doing plots and for sending plots to file instead.

The software (as well as the previously used underlying program library Hydrasub) is licenced as LGPL-based, meaning that anyone can modify and use it for other purposes than the original intent. 

This information is also available as a web page on http://folk.uio.no/trondr/layered. The content of this page can sometimes be a bit newer as it represents "bleeding edge" rather than "leading edge" development of the package.


Troubleshooting:

* The layeranalyzer package requires the Rcpp, coda, and for installations with vignettes also knitr, rmarkdown and markdown. When you start installing layeranalyzer, these packages should be automatically included. However, there may be contexts where it pays off to install them in advance.
    
* On Windows computers, you need the "rtools" program in addition to a standard R installation. Rtools can be found here: https://cran.r-project.org/bin/windows/Rtools/ .
    
* On Windows machines where rtools is installed via a company software center and you do not have direct admin rights, rtools may not function as intended. This will cause an error message when installing layeranalyzer saying "make not found". The path needs to be updated to where rtools store its program library, typically "c:\rtools\usr\bin". Search for "make.exe" on your computer and you should be able to find it. You then either need temporary admin rights to add this directory to the enviromental variable called "Path" or you need your IT support to do it for you.
    
* On Mac computers, the fortran compiler may be missing or outdated, see for instance https://stackoverflow.com/questions/58610155/problem-installing-r-package-ld-warning-directory-not-found-for-option and https://cran.r-project.org/bin/macosx/. This will cause a complaint about "gfortran not found" or similar messages. You then need to (re-)install fortran. A fortan compiler for Mac can be found at https://mac.r-project.org/tools/.
