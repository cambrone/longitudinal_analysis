#%%
import pandas as pd
import numpy as np 
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

from datetime import datetime
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss, average_precision_score,roc_auc_score
from sklearn.preprocessing import StandardScaler

from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant
from scipy.interpolate import make_interp_spline, BSpline
from scipy.special import expit
import random 

np.random.seed(222)
    
pd.set_option('display.max_columns', None)