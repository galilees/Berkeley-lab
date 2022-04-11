
import numpy as np
import pandas as pd
import json
import os

def open_file():
    path_to_json = '/net/cci-filer3/home/galilees/json/json_files/'
    data_list = []
    for file_name in [file for file in os.listdir(path_to_json) if file.endswith('.json')]:
        with open(path_to_json + file_name) as json_file:
            data = json.load(json_file)
            data_list.append(data)
        for data in data_list:
            success = data['success']
            if success == False:
                print(file_name)
    return data_list
   
def data_fail_count(data_list):
    count = 0 
    for data in data_list:
        success = data['success']
        if success == False:
            print(data)
            count += 1
            print()
    return count 
    

def data_check(data_list):
    checked_data = []

    for data in data_list:
    
        success = data['success']
        if success == True:
            checked_data.append(data)
            
    return checked_data

def get_data(checked_data, resname):
    data_list = []

    for dict in checked_data:
        data = dict[resname]
        data_list.append(data)
   
    return  data_list

def get_nearby_res(data_list):
    
    nearby_dict_list = []
    for i in data_list:
        for k, v in i.items():
             for key, val in v.items():
            
                 if key == "nearby_res":
                      nearby_dict_list.append(val)
            
    return (nearby_dict_list)
          

def get_hbonds(data_list):
    h_bond_dict_list = []

    for i in data_list:
        for k, v in i.items():
             for ke, va in v.items():
            
                 if ke == "n_hbonds":
                      h_bond_dict_list.append(va)
    return h_bond_dict_list


def data_frame(nearby_dict_list):
    df = pd.DataFrame(nearby_dict_list)
    return df

data_list = open_file()
#print("all data from json files: ", data_list)
print("Total number of files read:", len(data_list))

data_fail = data_fail_count(data_list)
print("The number of files that failed check: ", data_fail)

list_checked_data = data_check(data_list)
print("The number of files that passed the fail check: ",  len(list_checked_data))
#print("Data from files that that passed the fail check: ", list_checked_data )

gol_data = get_data(list_checked_data,'GOL')
#print("GOL data from files" , gol_data )
nearby_res_gol= get_nearby_res(gol_data)
#print("nearby residues:", nearby_res_gol)
h_bonds = get_hbonds(gol_data)
#print(h_bonds)

hoh_data = get_data(data_list, 'HOH')
#print("HOH data from files" , hoh_data )
nearby_res_hoh= get_nearby_res(hoh_data)
#print("nearby residues hoh:", nearby_res_hoh)

# original data with counts for glycerol and hoh extracted from the json files is added to frames called df_gol and df_hoh. hbonds data are added to df_hbonds  

df_gol = data_frame(nearby_res_gol)
print("Total number of glycerol: " , len(df_gol))
#print("gol data" , df_gol)
df_hbonds = data_frame(h_bonds)
#print(" hbonds", df_hbonds.columns)
df_hoh = data_frame(nearby_res_hoh)
#print(df_hoh)

# All residue counts besides amino acids and HOH are combined into a single column called Other the dataframes are then called df_gol_aa and df_hoh_aa 
final_table_columns =['HOH','ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
df_gol_sample_aa = df_gol[df_gol.columns.intersection(final_table_columns)]

df_hoh_sample_aa = df_hoh[df_hoh.columns.intersection(final_table_columns)]
#print(df_gol_sample_aa)
#print(df_gol_sample_aa.columns)

df_gol_sample_other = df_gol.drop(columns=[col for col in df_gol if col in final_table_columns])
df_hoh_sample_other = df_hoh.drop(columns=[col for col in df_hoh if col in final_table_columns])

df_gol_other = df_gol_sample_other.sum(axis=1)
df_hoh_other = df_hoh_sample_other.sum(axis=1)
#print(df_gol_other )

df_gol_aa = df_gol_sample_aa.assign(Other = df_gol_other)
#print(df_gol_aa)
df_hoh_aa = df_hoh_sample_aa.assign(Other = df_hoh_other) #use to assign new random batch
#print(df_hoh_aa)



df_gol_0 = df_gol_aa.fillna(0)
df_hoh_0 = df_hoh_aa.fillna(0)
#print(df_gol_0)
#print(df_hoh_0)
#print("columns" , df_gol_0.columns)

electroneg = [{'HOH':7,'ARG':9,'HIS':6,'LYS':3,'ASP':7,'GLU':7,'SER':3.5,'THR':3.5,'ASN':6.5,'GLN':6.5,'CYS':2.5, 'GLY':0,'PRO':0,'ALA':0,'VAL':0,'ILE':0,'LEU':0,'MET':2.5,'PHE':0,'TYR':3.5,'TRP':3}]
df_electroneg = pd.DataFrame(electroneg)
df_electroneg
df_repeated_gol = pd.concat([df_electroneg]*len(df_gol_0.index), ignore_index=True)

df_gol_aa_eneg = df_gol_0.mul(df_repeated_gol, axis='columns', level=None, fill_value=None)


df_repeated_hoh = pd.concat([df_electroneg]*len(df_hoh_0.index), ignore_index=True)

df_hoh_aa_eneg = df_hoh_0.mul(df_repeated_hoh, axis='columns', level=None, fill_value=None)
df_gol_aa_eneg = df_gol_aa_eneg.fillna(0)
df_hoh_aa_eneg = df_hoh_aa_eneg.fillna(0) 
#print(df_gol_aa_eneg)


eneg_gol = df_gol_aa_eneg.sum(axis=1)
eneg_hoh = df_hoh_aa_eneg.sum(axis=1)


#print(df_gol_0)

df_gol_0_1 = df_gol_0.mask(df_gol_0 > 0, 1)
df_hoh_0_1 = df_hoh_0.mask(df_hoh_0 > 0, 1)
#print(df_gol_0_1)


df_gol_sample_extended = df_gol_0_1.join(df_gol_0 ,how='inner',rsuffix="_count")
df_hoh_sample_extended = df_hoh_0_1.join(df_hoh_0 ,how='inner',rsuffix="_count")
df_gol_sample_extended['eneg_gol'] = eneg_gol
df_hoh_sample_extended['eneg_hoh'] = eneg_hoh
df_gol_sample_extended['h_bonds'] = df_hbonds

#print(df_gol_sample_extended)
#df_hoh_sample_extended.columns


df_gol_sample = df_gol_sample_extended
df_hoh_sample = df_hoh_sample_extended

"Choose a random selection of negative samples. The sample indexes can be saved for repeated testing "
df_hoh_sample = df_hoh_sample.sample(n=384)

df_hoh_sample_labels = df_hoh_sample.copy(deep = True) 
df_hoh_sample_labels['label'] = False

df_gol_sample_labels = df_gol_sample.copy(deep = True) 
df_gol_sample_labels['label'] = True
#print(df_gol_sample_labels)


features_df = pd.concat([df_gol_sample,df_hoh_sample])
features_df = features_df.fillna(0)
#features_df = features_df.drop(columns=[], axis=1)

#The features are amino acid one hot encoded values, amino acid counts, electronegativity measures, and H_bonds to GOL
print("Features for machine learning:", features_df.columns)
print(features_df)

labels_df = pd.concat([df_gol_sample_labels['label'],df_hoh_sample_labels['label']])


import sklearn
from sklearn.utils import shuffle
labels_df = shuffle(labels_df)

X = np.array(features_df)

y = np.array(labels_df)

#print(X,y)

from sklearn.preprocessing import StandardScaler
object = StandardScaler()
object.fit_transform(X)

# split into train test sets
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33)

# Use if not splitting the data sets
#X_train, X_test, y_train, y_test = X, X, y, y


# Another thing to consider is imbalance, if you have too few positive samples then the model might not work well. 
# Ideal would be 0.5, but having a larger number of negative is more realistic
imbalance = y.sum()/y.shape[0] # ratio of positive samples over total samples
print("Data imbalence:", imbalance)

# Logistic regression model
# from sklearn.datasets import make_classification
# from sklearn.linear_model import LogisticRegression
# from sklearn.model_selection import train_test_split
# from sklearn.pipeline import make_pipeline
# from sklearn.preprocessing import StandardScaler

# # # X, y = make_classification(random_state=42)
# # # X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
# # # pipe = make_pipeline(StandardScaler(), LogisticRegression())
# # # pipe.fit(X_train, y_train)  # apply scaling on training data
# # # Pipeline(steps=[('standardscaler', StandardScaler()),
# # #                 ('logisticregression', LogisticRegression())])

# # # pipe.score(X_test, y_test)  # apply scaling on testing data, without leaking training data.
# # # y_pred = pipe.predict(X_test)

# # # from sklearn.linear_model import LogisticRegression 
# # # from sklearn.metrics import accuracy_score 
# # # model = LogisticRegression(max_iter=1000)
# # # model.fit(X_train, y_train)
# # # y_pred = model.predict(X_test)


# # # #Import scikit-learn metrics module for accuracy calculation
# # # from sklearn import metrics
# # # # Model Accuracy, how often is the classifier correct?
# # # print("Accuracy:",metrics.accuracy_score(y_test, y_pred))

# # # # Model Precision: What proportion of positive identifications was actually correct?
# # # print("Precision:",metrics.precision_score(y_test, y_pred))

# # # # Model Recall:What proportion of true positives was identified correctly?
# # # #A model that produces no false negatives has a recall of 1.0.
# # # print("Recall:",metrics.recall_score(y_test, y_pred))


# # # from sklearn.metrics import confusion_matrix
# # # confusion_matrix(y_test,y_pred)
# # # pd.crosstab(y_test, y_pred, rownames = ['Actual'], colnames =['Predicted'], margins = True)



# Random Forest Model
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
#Create a Gaussian Classifier
clf=RandomForestClassifier(n_estimators=100)

#Train the model using the training sets y_pred=clf.predict(X_test)
clf.fit(X_train,y_train)

y_pred=clf.predict(X_test)

#Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics
# Model Accuracy, how often is the classifier correct?
print("Accuracy:",metrics.accuracy_score(y_test, y_pred))

# Model Precision: What proportion of positive identifications was actually correct?
print("Precision:",metrics.precision_score(y_test, y_pred))

# Model Recall:What proportion of true positives was identified correctly?
#A model that produces no false negatives has a recall of 1.0.
print("Recall:",metrics.recall_score(y_test, y_pred))

from sklearn.ensemble import GradientBoostingClassifier


import matplotlib.pyplot as plt
import seaborn as sns
def plot_cm(labels, predictions, p=0.5):
  cm = confusion_matrix(labels, predictions > p)
# #   plt.figure(figsize=(5,5))
# #   sns.heatmap(cm, annot=True, fmt="d")
# #   plt.title('Confusion matrix @{:.2f}'.format(p))
# #   plt.ylabel('Actual label')
# #   plt.xlabel('Predicted label')

  print('(True Negatives): ', cm[0][0])
  print('(False Positives): ', cm[0][1])
  print('(False Negatives): ', cm[1][0])
  print('(True Positives): ', cm[1][1])
  print('Total GOL: ', np.sum(cm[1]))

print(plot_cm(y_test,y_pred))
