# Stochastic_treatment_paper
 
**Abstract**

Motivated by an application to pre-exposure prophylaxis (PrEP) initiation studies, we propose a new treatment intervention dependent on the observed treatment process. We show there exist 1) estimators that are doubly and multiply robust to model misspecification, and 2) estimators that when used with machine learning algorithms can attain fast convergence rates for our proposed intervention.  

**Availability**

As per data confidentiality agreements, the data cannot be made publicly available because of the sensitivity of the information. Because of this, we have also provided numerous simulations in the paper using data we have generated ourselves. 

**Simulated data**

We conducted two simulation studies for our proposed multiplicative shift intervention distribution. Our first simulation findings show that the weighted ICE is more robust to model misspecification than IPW and ICE when the nuisance functions are estimated using parametric models. Our second simulation findings show that the TMLE with sample-splitting and cross-fitting is consistent as long as the nuisance functions are estimated consistently at fast enough rates using machine learning methods, which may not necessarily be n^−1/2. 

**Code**

Description:

This code is delivered via the files described below. There is no licensing information for this code (and it will be developed into a complete package after the submission of this paper). 
This code includes all of the functions necessary to run the results found in Section 8 of the paper (the results run on simulated, non-proprietary data). The main two folders are:
- parametric_models: containing codes to run various estimators described in the paper with parametric models at different sample sizes and at different delta values
- ML_sims: containing codes to run various estimators described in the paper with machine learning methods at different sample sizes and at delta=0.5

Instructions for use:

The contents of the parametric_models folder are as follows:
- n*/delta_*/datagen.R: code to generate the simulated datasets
- n*/delta_*/datagen_wide_true.R: code to generate true values 
- n*/delta_*/*method*.R: code to run the methods described in the paper (IPW, ICE, DR (weighted ICE); note that dr_wrongOR.R is an R script where outcome regression models are misspecified, and dr_wrongPS.R is an R script where propensity score models are misspecified).
- sim_analysis: code to compare bias, SE and MSE found in the paper
- jplusone/n*/delta*/...: contains the files to run the methods described in the paper to show J+1 multiple robustness of the weighted ICE estimator

The contents of the ML_sims folder are as follows:
- n*/delta_*/datagen.R: code to generate the simulated datasets
- n*/delta_*/datagen_wide_true.R: code to generate true values 
- n*/delta_*/*method*.R: code to run the methods described in the paper (IPW, ICE, TMLE)
- n*/delta_*/*dr_ss_asymvar*.R: code to calculate asymptotic variance of the TMLE estimator with sample-splitting and cross-fitting
- sim_analysis: code to compare bias, SE and MSE found in the paper

Software requirements:
R packages to run necessary functions include:
- library(SuperLearner): Version 2.0-26
- library(dplyr): Version 	1.0.5
- library(glm2): Version 1.2.1
- library(data.table): Version 1.14.0
In addition, parallel processing is done using the doParallel (Version 1.0.16) and foreach (Version 1.5.1) library. 

The following workflow should be followed if running the weighted ICE procedure in the parametric_models folder for sample size of 500 and delta of 0.5:
1. Go to parametric_models/n500/delta_0.5/dr.R
2. Install and load necessary packages
3. Set the number of clusters for parallel processing of code
4. Source the function code: myfunc()
5. Run the foreach(m=1:1000) %dopar% myfunc(m) and save results in CSV file
