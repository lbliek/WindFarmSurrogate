# PIBO: Optimal Transport based Permutation Invariant Bayesian Optimization

This repository refers to the paper **_"Optimal Transport based Permutation Invariant Bayesian Optimization of offshore wind farm layouts"_**

(add link to paper...)

# Repository' s organization:
 * **run_PIBO.R** - R scprit to run the PIBO (Permutation Invariant Bayesian Optimization) agorithm
 * **run_vanilla_BO_on_flows.R** - R script to run _vanilla_ BO within the flows space (i.e., the same space of PIBO, but without dealing with permutation invariance)
 * **run_vanilla_BO_on_pointclouds.R** - R script to run _vanilla_ BO directly in the space of point-clouds (i.e., the physical space where the m=5 wid turbines must be optimally placed).
 * **utils.R** - R script consisting on a single functions to evaluate the permutation-invariant objectiv function
 * **ensemble_model_predict.py** - Python script to load and use an ensemble of 5 xgboost models predicting the generated power for the input wind farm layout. The model has been trained to be consistent with permutation invariance but its input it is not (that is why we created the function in the 'utils.R' script).
 * **permutation_invariant_ensemble_model.pkl** - the trained ensemble of 5 xgboost models predicting the generated power for a given wind far layout
 * ........   

# Why Permutation-Invariant?
(add the picture of the Bird function?)

# Important conclusions
(Convergene Plot)
(Solutions in the Physical Space - relevance of OT)
