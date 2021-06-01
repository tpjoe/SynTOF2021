"""
This script perform machine learning on the clustered data for AD and LBD
"""

#### import libraries
import numpy as np
import pandas as pd
# import cudf
# import cupy as cp
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn import preprocessing
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon, mannwhitneyu, spearmanr
# from utils_ML import cuLogisticRegressionCV, customLogisticRegressionCV
# from cuml.linear_model import LogisticRegression as cuLogisticRegression
from scipy.stats import norm
from xgboost import XGBClassifier, XGBRegressor
from sklearn.svm import SVC, SVR
from sklearn.linear_model import Ridge, Lasso, ElasticNetCV, ElasticNet
from joblib import parallel_backend
from dask.distributed import Client, LocalCluster
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier


#### Importing data ------------------------------------------------------------------
regions = ['BA9', 'DLCau', 'Hipp']
for region in regions:
    df_ = pd.read_csv(''.join(['R_py_exchange/df_meanAllMarkers_', region, '_noStd_exp2mc_5_13.csv'])).iloc[:, 1:]
    if region == 'BA9':
        df_.columns = df_.columns[0:2].tolist() + [region + '_' + i for i in df_.columns[2:df_.shape[1]].tolist()]
        df = df_
    else:
        #check that samples are in the same order
        print(all(df.iloc[:, 1] == df_.iloc[:, 1]))
        df_ = df_.drop(['group', 'sample'], axis=1)
        df_.columns = [region + '_' + i for i in df_.columns.tolist()]
        df = pd.concat([df.reset_index(drop=True), df_.reset_index(drop=True)], axis=1)


# filter data
y_pred = []
pair = ['LowNo', 'LBD']     #####<<<<<<<<<<<<< This is where you select either ['LowNo', 'LBD'] or ['LowNo', 'PHAD']
df_ = df.loc[df.group.isin(pair), :]
X = df_.drop(['group', 'sample'], axis=1)
y = df_.group.apply(lambda x: 0 if x=='LowNo' else 1).astype('float64').to_numpy()

# remove super highly correlated features
if pair[1] == 'PHAD':
    # X = X.loc[:, ~(np.array(['Hipp' in i for i in X.columns]) & np.array(['p-Tau' in i for i in X.columns]))]
    logP = X.apply(lambda x: -(np.log10(2*norm.cdf(-np.abs(np.sqrt((X.shape[0]-3)/1.06)*np.arctanh(spearmanr(x, y)[0]))))))
    exclude = logP[logP>-(np.log10(0.05/len(logP)))].index
    X = X.loc[:, ~X.columns.isin(exclude)]

if pair[1] == 'LBD':
    logP = X.apply(lambda x: -(np.log10(2*norm.cdf(-np.abs(np.sqrt((X.shape[0]-3)/1.06)*np.arctanh(spearmanr(x, y)[0]))))))
    exclude = logP[logP>3.5].index
    X = X.loc[:, ~X.columns.isin(exclude)]


X = X.fillna(X.mean())
X = X.drop(X.columns[X.std()==0], axis=1)
X = X.loc[:, X.sum(axis=0) != 0]

# preprocessing and group filter
colnames = X.columns
X = X.to_numpy()
X = preprocessing.scale(X)


#### with concatenated features ---------------------------------------------------------------------------------------------------------
models = ['LASSO', 'Ridge', 'EN', 'Random Forest', 'SVM']
res = pd.DataFrame(np.zeros((len(models), 2)), index=models, columns=['AUC', 'P-value'])
y_preds = {k:pd.DataFrame(np.zeros((len(y), 1)), columns=['y_pred']) for k in models}

# with parallel_backend(backend="dask"):
for algo in models:
    loo = LeaveOneOut()
    for train_index, test_index in loo.split(X):
        models = {
                    'EN': LogisticRegression(penalty='elasticnet', l1_ratio=0.5, fit_intercept=False, solver='saga', max_iter=10000), #'EN': customLogisticRegressionCV(cv=5, l1_ratios=[0, 0.25, 0.5, 0.75, 1], fit_intercept=False),
                    'LASSO': LogisticRegression(penalty='l1', fit_intercept=False, solver='saga', max_iter=10000), #customLogisticRegressionCV(cv=5, l1_ratios=[1], fit_intercept=False), #
                    'Ridge': LogisticRegression(penalty='l2', fit_intercept=False), #customLogisticRegressionCV(cv=5, l1_ratios=[0], fit_intercept=False), #
                    'Random Forest': RandomForestClassifier(),
                    'KNN': KNeighborsClassifier(),
                    'SVM': SVC(probability=True),  #rbf for AD, linear for LBD
                    'XGBoost': XGBClassifier(n_jobs=-1, booster='gblinear', objective='binary:logistic')
                }
        clf = models[algo]
        X_train, X_test = X[train_index, :], X[test_index, :]
        y_train, y_test = y[train_index], y[test_index]
        clf.fit(np.array(X_train), np.array(y_train))
        # if algo in ['LASSO', 'Ridge']:
            # y_pred = np.array(clf.best_estimator_.predict_proba(np.array(X_test))[:, 1])
        # else:
        y_pred = np.array(clf.predict_proba(np.array(X_test))[:, 1])
        y_preds[algo].loc[test_index, :] = y_pred
    # evaluation
    res.loc[algo, 'P-value'] = mannwhitneyu(y_preds[algo][y==0], y_preds[algo][y==1]).pvalue
    res.loc[algo, 'AUC'] = roc_auc_score(y, y_preds[algo])
    print(res.loc[algo, 'AUC'])



best_mod = res.loc[res.AUC==res.AUC.max(), :].index[0]
res.to_csv('R_py_exchange_afterCluster/res_' + str(pair[1]) + '_noPTau.csv')
pd.concat([y_preds[best_mod], pd.Series(y)], axis=1).to_csv('R_py_exchange_afterCluster/preds_' + str(pair[1]) + '_'+ best_mod + '_noPTau.csv')





# Compute ROC curve and ROC area for each class
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib.figure import figaspect
import matplotlib as mpl

fpr = dict()
tpr = dict()
roc_auc = dict()
gr_truth = y


font = "Arial"
hfont = {'fontname':'Helvetica'}
for i, algo in enumerate(y_preds.keys()):
    fpr[i], tpr[i], _ = roc_curve(gr_truth, np.array(y_preds[algo]))
    roc_auc[i] = auc(fpr[i], tpr[i])

plt.rc('grid', linestyle="dotted", color='lightgray')
color = ['#22162B', '#724E91', '#EE4266', '#f8a530', '#93B5C6']
w, h = figaspect(1.)
plt.figure(figsize=(w, h))
mpl.rcParams['axes.linewidth'] = 0.3
lw = 1.5
fontsize = 14
plt.plot([0, 1], [0, 1], color='lightgrey', lw=lw, linestyle='--')
for i, algo in enumerate(y_preds.keys()):
    plt.plot(fpr[i], tpr[i], color=color[i], 
            lw=lw, label='{} = {:.2f}'.format(algo, roc_auc[i]))

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xticks(fontsize=fontsize, fontname=font)
plt.yticks(fontsize=fontsize, fontname=font)
plt.xlabel('False Positive Rate', size=fontsize, **hfont)
plt.ylabel('True Positive Rate', size=fontsize, **hfont)
plt.title('Control/LBD ROC and AUC', size=fontsize)
plt.text(0.815, 0.4, '$\it{n=}$' + str(len(y)), fontsize=fontsize)
plt.legend(loc="lower right", fontsize=fontsize-4)
plt.tight_layout()
plt.grid(True)
plt.savefig('figures/ROC/LBD_all_noPerfect2.pdf')





#### rerun the best model again for getting prediction values and weights --------------------------------------------------------
y_pred = np.zeros((len(y)))
wt = np.zeros((len(y), X.shape[1]))
loo = LeaveOneOut()

for i, (train_index, test_index) in enumerate(loo.split(X)):
    clf = LogisticRegression(fit_intercept=False, penalty='elasticnet', solver='saga', l1_ratio=0.5, max_iter=10000) # <--- for AD
    # clf = LogisticRegression(penalty='l2', fit_intercept=False)   <--- for LBD
    X_train, X_test = X[train_index, :], X[test_index, :]
    y_train, y_test = y[train_index], y[test_index]
    clf.fit(np.array(X_train), np.array(y_train))
    y_pred[i] = np.array(clf.predict_proba(np.array(X_test))[:, 1])
    wt[i, ] = clf.coef_[0]

# evaluation
p_value = mannwhitneyu(y_pred[y==0], y_pred[y==1]).pvalue
auc = roc_auc_score(y, y_pred)

# save model weights
if pair[1]=='PHAD':
    wt_df = pd.DataFrame(wt, columns=df.columns[~df.columns.isin(exclude)][2:])
else:
    wt_df = pd.DataFrame(wt, columns=df.columns[~df.columns.isin(exclude)][2:])
wt_df.to_csv('R_py_exchange_afterCluster/wt_' + pair[0] + '_' + pair[1] + '_Ridge.csv')

# save predictions
aa = pd.concat([pd.Series(y_pred), pd.Series(y), df.loc[df.group.isin(pair), ['group', 'sample']].reset_index(drop=True)], axis=1)
aa.sort_values([1, 0])

