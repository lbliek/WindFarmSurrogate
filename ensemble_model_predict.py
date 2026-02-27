# Python code for running the ensemble model that predicts the amount of
# generated power for a given (vectorized) point-clouds of 5 wind turbines,
# (x1,y1,x2,y2,x3,y3,x4,y4,x5,y5).
# The objective function in this Python script is not permutation-invariant: this
# aspect is managed in the 'utils.R' file containing the permutation-invariant
# objective function


import pickle
import numpy as np
import xgboost as xgb


# Define median-based ensemble model
class medClassifier:
    def __init__(self, classifiers=None):
        self.classifiers = classifiers

    def predict(self, X):
        self.predictions_ = list()
        for classifier in self.classifiers:
            try:
                self.predictions_.append(classifier.predict(X)) # used for the random forest that is part of the ensemble
            except:
                X = xgb.DMatrix(X)
                self.predictions_.append(classifier.predict(X)) # used for the XGBoost models that are part of the ensemble
        med1 = np.median(self.predictions_, axis=0) # median of predictions
        mean1 = np.mean(self.predictions_, axis=0) # mean of predictions
        out = med1 + np.random.rand()*np.abs(med1-mean1) # add noise if median is far from mean, indicating more uncertainty, also all noise is positive to focus on minimizing parts with more certainty
        return out

      
# Load Ensemble model
Ensemblefile = "permutation_invariant_ensemble_model.pkl"
with open(Ensemblefile, 'rb') as file:
    Ensemble = pickle.load(file)


# This is the first objective function
def objective(x):
    x = np.array([x])
    pred = Ensemble.predict(x)
    return pred[0]
