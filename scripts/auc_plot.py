import pandas as pd
import numpy as np
import plotly.express as px
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from sklearn.datasets import make_classification
import plotly
import sklearn.metrics as metrics
import os


dc_mayo_ac = pd.read_excel('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\7.PASNet\\Activity\\Output\\Mayo_test\\Mayo_TEST.xlsx')
fl_mayo_ac = pd.read_csv('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\7.PASNet\\Activity\\Output\\Mayo_test\\PASNet_pred_0_0.txt', sep=',', names=['first', 'second'])

dc_stanford_ac = pd.read_excel('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\7.PASNet\\Activity\\Output\\Stanford_test\\STANFORD_TEST.xlsx')
fl_stanford_ac = pd.read_csv('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\7.PASNet\\Activity\\Output\\Stanford_test\\PASNet_pred_0_0.txt', sep=',', names=['first', 'second'])

dc_scmgh_ac = pd.read_excel('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_2\\7.PASNet\\Output\\scMGH_TEST.xlsx')
fl_scmgh_ac = pd.read_csv('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_2\\7.PASNet\\Output\\PASNet_pred_0_0.txt', sep=',', names=['first', 'second'])

dc_mayo_ex = pd.read_excel('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\7.PASNet\\Expression\\Output\\Mayo_test\\MAYO_TEST.xlsx')
fl_mayo_ex = pd.read_csv('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\7.PASNet\\Expression\\Output\\Mayo_test\\PASNet_pred_0_0.txt', sep=',', names=['first', 'second'])

dc_stanford_ex = pd.read_excel('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\7.PASNet\\Expression\\Output\\Stanford_test\\STANFORD_TEST.xlsx')
fl_stanford_ex = pd.read_csv('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\7.PASNet\\Expression\\Output\\Stanford_test\\PASNet_pred_0_0.txt', sep=',', names=['first', 'second'])



fl_mayo_ac['combined']= fl_mayo_ac.values.tolist()
fl_stanford_ac['combined']= fl_stanford_ac.values.tolist()
fl_scmgh_ac['combined']= fl_scmgh_ac.values.tolist()
fl_mayo_ex['combined']= fl_mayo_ex.values.tolist()
fl_stanford_ex['combined']= fl_stanford_ex.values.tolist()


y_pred_mayo_ac = np.array(list(fl_mayo_ac.combined))
y_pred_stanford_ac = np.array(list(fl_stanford_ac.combined))
y_pred_scmgh_ac = np.array(list(fl_scmgh_ac.combined))
y_pred_mayo_ex = np.array(list(fl_mayo_ex.combined))
y_pred_stanford_ex = np.array(list(fl_stanford_ex.combined))


def vectorized_label(target, n_class):

	TARGET = np.array(target).reshape(-1)

	return np.eye(n_class)[TARGET]

y_mayo_ac = dc_mayo_ac['Condition'].astype(int)
y_true_mayo_ac = vectorized_label(y_mayo_ac, 2)

y_stanford_ac = dc_stanford_ac['Condition'].astype(int)
y_true_stanford_ac = vectorized_label(y_stanford_ac, 2)

y_scmgh_ac = dc_scmgh_ac['Condition'].astype(int)
y_true_scmgh_ac = vectorized_label(y_scmgh_ac, 2)

y_mayo_ex = dc_mayo_ex['Condition'].astype(int)
y_true_mayo_ex = vectorized_label(y_mayo_ex, 2)

y_stanford_ex = dc_stanford_ex['Condition'].astype(int)
y_true_stanford_ex = vectorized_label(y_stanford_ex, 2)

fpr_list = []
tpr_list = []
threshold_list = []
roc_auc_list = []

for i in range(1): # you can make this more general
    fpr_mayo_ac, tpr_mayo_ac, threshold_mayo_ac = metrics.roc_curve(y_true_mayo_ac[:, i], y_pred_mayo_ac[:, i])
    roc_auc_mayo_ac = metrics.auc(fpr_mayo_ac, tpr_mayo_ac)
    fpr_list.append(fpr_mayo_ac)
    tpr_list.append(tpr_mayo_ac)
    threshold_list.append(threshold_mayo_ac)
    roc_auc_list.append(roc_auc_mayo_ac)
    
for i in range(1): # you can make this more general
    fpr_stanford_ac, tpr_stanford_ac, threshold_stanford_ac = metrics.roc_curve(y_true_stanford_ac[:, i], y_pred_stanford_ac[:, i])
    roc_auc_stanford_ac = metrics.auc(fpr_stanford_ac, tpr_stanford_ac)
    fpr_list.append(fpr_stanford_ac)
    tpr_list.append(tpr_stanford_ac)
    threshold_list.append(threshold_stanford_ac)
    roc_auc_list.append(roc_auc_stanford_ac)    
  
for i in range(1): # you can make this more general
    fpr_scmgh_ac, tpr_scmgh_ac, threshold_scmgh_ac = metrics.roc_curve(y_true_scmgh_ac[:, i], y_pred_scmgh_ac[:, i])
    roc_auc_scmgh_ac = metrics.auc(fpr_scmgh_ac, tpr_scmgh_ac)
    fpr_list.append(fpr_scmgh_ac)
    tpr_list.append(tpr_scmgh_ac)
    threshold_list.append(threshold_scmgh_ac)
    roc_auc_list.append(roc_auc_scmgh_ac)  

for i in range(1): # you can make this more general
    fpr_mayo_ex, tpr_mayo_ex, threshold_mayo_ex = metrics.roc_curve(y_true_mayo_ex[:, i], y_pred_mayo_ex[:, i])
    roc_auc_mayo_ex = metrics.auc(fpr_mayo_ex, tpr_mayo_ex)
    fpr_list.append(fpr_mayo_ex)
    tpr_list.append(tpr_mayo_ex)
    threshold_list.append(threshold_mayo_ex)
    roc_auc_list.append(roc_auc_mayo_ex)  

for i in range(1): # you can make this more general
    fpr_stanford_ex, tpr_stanford_ex, threshold_stanford_ex = metrics.roc_curve(y_true_stanford_ex[:, i], y_pred_stanford_ex[:, i])
    roc_auc_stanford_ex = metrics.auc(fpr_stanford_ex, tpr_stanford_ex)
    fpr_list.append(fpr_stanford_ex)
    tpr_list.append(tpr_stanford_ex)
    threshold_list.append(threshold_stanford_ex)
    roc_auc_list.append(roc_auc_stanford_ex)   
  
d_mayo_ac = {'FPR':fpr_mayo_ac,'TPR':tpr_mayo_ac}
d_stanford_ac = {'FPR':fpr_stanford_ac,'TPR':tpr_stanford_ac}
d_scmgh_ac = {'FPR':fpr_scmgh_ac,'TPR':tpr_scmgh_ac}
d_mayo_ex = {'FPR':fpr_mayo_ex,'TPR':tpr_mayo_ex}
d_stanford_ex = {'FPR':fpr_stanford_ex,'TPR':tpr_stanford_ex}

df_mayo_ac = pd.DataFrame(d_mayo_ac)
df_stanford_ac = pd.DataFrame(d_stanford_ac)
df_scmgh_ac = pd.DataFrame(d_scmgh_ac)
df_mayo_ex = pd.DataFrame(d_mayo_ex)
df_stanford_ex = pd.DataFrame(d_stanford_ex)

df_mayo_ac['Model'] = 'APNET_Mayo (AUC: 0.96)'
df_stanford_ac['Model'] = 'APNet_Stanford (AUC: 0.91)'
df_scmgh_ac['Model'] = 'APNet_scMGH (AUC: 0.99)'
df_mayo_ex['Model'] ='PASNet_exp_Mayo (AUC: 0.64)'
df_stanford_ex['Model'] = 'PASNet_exp_Stanford (AUC: 0.89)'

## Values loaded from Random Forest
d_mayo_rf = {'FPR':[0,0.01104972,1],
              'TPR':[0,0.31343284,1],
              'Model':['RF_Mayo (AUC: 0.65)', 'RF_Mayo (AUC: 0.65)', 'RF_Mayo (AUC: 0.65)']}
df_mayo_rf = pd.DataFrame(d_mayo_rf)

d_stanford_rf = {'FPR':[0,0.025,1],
              'TPR':[0,0.5,1],
              'Model':['RF_Stanford (AUC: 0.73)', 'RF_Stanford (AUC: 0.73)', 'RF_Stanford (AUC: 0.73)']}
df_stanford_rf = pd.DataFrame(d_stanford_rf)


d_scmgh_rf = {'FPR':[0,0.17029549,1],
              'TPR':[0,0.90401052,1],
              'Model':['PASNet_exp_scMGH (AUC: 0.87)', 'PASNet_exp_scMGH (AUC: 0.87)', 'PASNet_exp_scMGH (AUC: 0.87)']}
df_scmgh_rf = pd.DataFrame(d_scmgh_rf)

df = pd.concat([df_mayo_ac, df_stanford_ac, df_scmgh_ac, df_mayo_ex, 
               df_stanford_ex, df_mayo_rf, df_stanford_rf,
               df_scmgh_rf])

os.mkdir('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Figures\\tables_auc')
df.to_csv('c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\Figures\\tables_auc\\df_auc.csv', index=False)

