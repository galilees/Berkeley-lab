import numpy as np
import pandas as pd
import json
import os, time
import random

import tensorflow as tf
from tensorflow import keras
import tempfile
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

path_to_json = '/net/anaconda/raid1/dorothee/14_frontiers_QR_restraints/galilee/pdb_survey/'
#path_to_json = '/net/anaconda/raid1/dorothee/14_frontiers_QR_restraints/galilee/pdb_survey/'



def data_frame(list_):
    df = pd.DataFrame(list_)
    return df

def run():
  ''' Parse json files create dataframes and neural network 
  '''
  os.chdir(path_to_json)
  n_json = 0
  n_success = 0
  list_fail = []
  nearby_dict_list_gol = []
  nearby_dict_list_hoh = []
  n_hbonds_list = []
  n_ss_helix_list = []
  n_ss_beta_list = []
  for file_name in [file for file in os.listdir(path_to_json) if file.endswith('.json')]:
  #for root, dirs, files in os.walk(path_to_json):
    #if os.path.basename(root).startswith('queue_logs'): continue
    #for file in files:
    #  if not file.endswith('.json'): continue
      with open(file_name, 'r') as fp:
        n_json += 1
        if n_json > 8000: break
      
        #print(file_name)
        data = json.load(fp)
        success = data['success']
        if not success:
          list_fail.append(data['pdb_code'])
          continue
        n_success += 1
        # GOL
        all_gols_dict = data['GOL']
        #print(all_gols_dict)

        for sel_str, gol_data in all_gols_dict.items():
          if not gol_data: continue
          if 'n_hbonds' not in gol_data: continue
          if 'nearby_res' not in gol_data: continue
          if 'n_ss_helix'not in gol_data: continue
          if 'n_ss_beta' not in gol_data: continue
          if gol_data['n_hbonds'] is None: continue
          if gol_data['n_ss_helix'] is None: continue
          if gol_data['n_ss_beta'] is None: continue

          n_ss_beta_list_gol = []
          n_hbonds_list.append(gol_data['n_hbonds'])
          n_ss_helix_list.append(gol_data['n_ss_helix'])
          n_ss_beta_list.append(gol_data['n_ss_beta'])
          nearby_dict_list_gol.append(gol_data['nearby_res'])
      
        
        all_hohs_dict = data['HOH']

        for sel_str, hoh_data in all_hohs_dict.items():
          if 'nearby_res' not in hoh_data: continue
          if 'n_ss_helix'not in hoh_data: continue
          if hoh_data['n_ss_helix'] is None: continue
          if hoh_data['n_ss_beta'] is None: continue
          if 'nearby_res' not in hoh_data: continue
  
          nearby_dict_list_hoh.append(hoh_data['nearby_res'])
          n_ss_helix_list.append(hoh_data['n_ss_helix'])
          n_ss_beta_list.append(hoh_data['n_ss_beta'])
      
  # print(len(n_hbonds_list))
  # print(len(nearby_dict_list_gol))
  assert(len(n_hbonds_list) == len(nearby_dict_list_gol))
 
  print('Number of json files read: ', n_json)
  print('Number of success: ', n_success)
  print('Number of failues', len(list_fail))
  print("neural network")
  #print('failures')
  # for p in list_fail:
  #   print(p)

#print("Total number of glycerol: " , len(gol_data))

#gol_data_sample = gol_data
  df_gol = data_frame(list_ = nearby_dict_list_gol)
  #print(df_gol.head(10))
  df_hbonds = data_frame(n_hbonds_list)
  df_helix = data_frame(n_ss_helix_list)
  df_beta = data_frame(n_ss_beta_list)
  #print(" hbonds", df_hbonds)
#1550
  print("Total number of negative samples", len(nearby_dict_list_hoh))
  print("Total number of positive samples",  len(nearby_dict_list_gol))

  hoh_data_sample = random.sample(nearby_dict_list_hoh,len(nearby_dict_list_hoh))
  df_hoh = data_frame(hoh_data_sample)
  #print(df_hoh.head(10))
#

  ## All residue counts besides amino acids and HOH are combined into a single column called Other the dataframes are then called df_gol_aa and df_hoh_aa
  final_table_columns =['HOH','ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
  df_gol_sample_aa = df_gol[df_gol.columns.intersection(final_table_columns)]
  #
  df_hoh_sample_aa = df_hoh[df_hoh.columns.intersection(final_table_columns)]
  # print(df_gol_sample_aa)
  # print(df_gol_sample_aa.columns)
  #
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
  #plot the maximum count of amino acids nearby gol
  df_gol_sample_aa_count = df_gol_0.drop(columns= ['HOH','Other'])
  df_gol_sample_aa_count.rename(columns ={'ARG':'R','HIS':'H','LYS':'K','ASP':'D','GLU':'E','SER':'S','THR':'T','ASN':'N','GLN':'Q','CYS':'C','GLY':'G','PRO':'P','ALA':'A','VAL':'V','ILE':'I','LEU':'L','MET':'M','PHE':'F','TYR':'Y','TRP':'W'}, inplace = True)
  y = df_gol_sample_aa_count.max().to_numpy().flatten()
  x = df_gol_sample_aa_count.columns
  print('amino acids', x)
  print('counts', y)

  fig, ax = plt.subplots()
  plt.xlabel('Amino Acids')
  plt.ylabel('Count')
  ax.bar(x,y)
  plt.show()
  
 
  #print("columns" , df_gol_0.columns)

  electroneg = [{'HOH':3.5,'ARG':9,'HIS':6,'LYS':3,'ASP':7,'GLU':7,'SER':3.5,'THR':3.5,'ASN':6.5,'GLN':6.5,'CYS':2.5,'GLY':0,'PRO':0,'ALA':0,'VAL':0,'ILE':0,'LEU':0,'MET':2.5,'PHE':0,'TYR':3.5,'TRP':3}]
  df_electroneg = pd.DataFrame(electroneg)

  df_repeated_gol = pd.concat([df_electroneg]*len(df_gol_0.index), ignore_index=True)

  df_gol_aa_eneg = df_gol_0.mul(df_repeated_gol, axis='columns', level=None, fill_value=None)


  df_repeated_hoh = pd.concat([df_electroneg]*len(df_hoh_0.index), ignore_index=True)

  df_hoh_aa_eneg = df_hoh_0.mul(df_repeated_hoh, axis='columns', level=None, fill_value=None)
  df_gol_aa_eneg = df_gol_aa_eneg.fillna(0)
  df_hoh_aa_eneg = df_hoh_aa_eneg.fillna(0)
  #print(df_gol_aa_eneg)


  eneg_gol = df_gol_aa_eneg.sum(axis=1)
  eneg_hoh = df_hoh_aa_eneg.sum(axis=1)
  eneg = eneg_gol.append(eneg_hoh).reset_index(drop=True)
 
  # print(eneg_gol)
  # print(eneg_hoh)
  #print(eneg)

  carbons = [{'HOH':0,'ARG':10.2,'HIS':10.2,'LYS':10.2,'ASP':5.1,'GLU':7.65,'SER':2.55,'THR':5.1,'ASN':5.1,'GLN':7.65,'CYS':5.1, 'GLY':0,'PRO':7.65,'ALA':2.55,'VAL':7.65,'ILE':7.65,'LEU':7.65,'MET':10.2,'PHE':17.85,'TYR':17.85,'TRP':22.95}]
  df_carbons = pd.DataFrame(carbons)

  df_repeated_gol_ = pd.concat([df_carbons]*len(df_gol_0.index), ignore_index=True)

  df_gol_aa_carbons = df_gol_0.mul(df_repeated_gol_, axis='columns', level=None, fill_value=None)


  df_repeated_hoh_ = pd.concat([df_carbons]*len(df_hoh_0.index), ignore_index=True)

  df_hoh_aa_carbons = df_hoh_0.mul(df_repeated_hoh_, axis='columns', level=None, fill_value=None)
  df_gol_aa_carbons = df_gol_aa_carbons.fillna(0)
  df_hoh_aa_carbons = df_hoh_aa_carbons.fillna(0)
  carbons_gol = df_gol_aa_carbons.sum(axis=1)
  carbons_hoh = df_hoh_aa_carbons.sum(axis=1)
  carbons = carbons_gol.append(carbons_hoh).reset_index(drop=True)


  # print(df_gol_0)

  df_gol_0_1 = df_gol_0.mask(df_gol_0 > 0, 1)
  df_hoh_0_1 = df_hoh_0.mask(df_hoh_0 > 0, 1)
  #print(df_gol_0_1)

  # print(df_gol_0.head(10)
  # print(df_hoh_0.head(10))
  df_gol_sample_extended = df_gol_0_1.join(df_gol_0 ,how='inner',rsuffix="_count")
  df_hoh_sample_extended = df_hoh_0_1.join(df_hoh_0 ,how='inner',rsuffix="_count")
  # df_gol_sample_extended['h_bonds'] = df_hbonds
 

  #print(df_gol_sample_extended)
  #df_hoh_sample_extended.columns


  df_gol_sample = df_gol_sample_extended
  df_hoh_sample = df_hoh_sample_extended
  

  df_hoh_sample_labeled = df_hoh_sample.copy(deep = True)
  df_hoh_sample_labeled['label'] = False

  df_gol_sample_labeled = df_gol_sample.copy(deep = True)
  df_gol_sample_labeled['label'] = True
  #print(df_gol_sample_labels)
 
  features_df = pd.concat([df_gol_sample,df_hoh_sample])
  features_df['eneg'] = eneg
  features_df['n_helix'] = df_helix
  features_df['n_beta'] = df_beta
  features_df['carbons'] = carbons
  
  #print(features_df)
 

  labels_df = pd.concat([df_gol_sample_labeled['label'],df_hoh_sample_labeled['label']])
  #print(X,y)
  
  X = np.array(features_df)
  
  y = np.array(labels_df)
  
  # scaler = StandardScaler()
  # # transform data
  # scaled = scaler.fit_transform(X)
  # X = scaled
  
  
  # split into train test sets
  from sklearn.model_selection import train_test_split
  X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33)
  X_train, X_val_and_test, y_train, y_val_and_test = train_test_split(X, y, test_size=0.33)
  X_val, X_test, y_val, y_test = train_test_split(X_val_and_test, y_val_and_test, test_size=0.5)
  
  
  # Another thing to consider is imbalance, if you have too few positive samples then the model might not work well.
  # Ideal would be 0.5, but having a larger number of negative is more realistic
  imbalance = y.sum()/y.shape[0] # ratio of positive samples over total samples
  print("Data imbalence:", imbalance)

  
  mpl.rcParams['figure.figsize'] = (12, 10)
  colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
  
  METRICS = [
       keras.metrics.TruePositives(name='tp'),
       keras.metrics.FalsePositives(name='fp'),
       keras.metrics.TrueNegatives(name='tn'),
       keras.metrics.FalseNegatives(name='fn'),
       keras.metrics.BinaryAccuracy(name='accuracy'),
       keras.metrics.Precision(name='precision'),
       keras.metrics.Recall(name='recall'),
       
  ]
  
  def make_model(metrics=METRICS, output_bias=None):
   if output_bias is not None:
     output_bias = tf.keras.initializers.Constant(output_bias)
   model = keras.Sequential([
       keras.layers.Dense(
           16, activation='relu',
           input_shape=(X_train.shape[-1],)),
       keras.layers.Dropout(0.5),
       keras.layers.Dense(1, activation='sigmoid',
                          bias_initializer=output_bias),
   ])
  
   model.compile(
       optimizer=keras.optimizers.Adam(learning_rate=1e-3),
       loss=keras.losses.BinaryCrossentropy(),
       metrics=metrics)
  
   return model
  
  EPOCHS = 100
  BATCH_SIZE = 1000
  
  early_stopping = tf.keras.callbacks.EarlyStopping(
     monitor='val_prc',
     verbose=1,
     patience=10,
     mode='max',
     restore_best_weights=True)
  
  
  neg, pos = np.bincount(labels_df )
  total = neg + pos
  print('Examples:\n    Total: {}\n    Positive: {} ({:.2f}% of total)\n'.format(
     total, pos, 100 * pos / total))
  print("Features for neural network:", features_df.columns)
  initial_bias = np.log([pos/neg])


  # Scaling by total/2 helps keep the loss to a similar magnitude.
# The sum of the weights of all examples stays the same.
  weight_for_0 = (1 / neg) * (total / 2)
  weight_for_1 = (1 / pos) * (total / 2)

  class_weight = {0: weight_for_0, 1: weight_for_1}
  print('Weight for class 0: {:.2f}'.format(weight_for_0))
  print('Weight for class 1: {:.2f}'.format(weight_for_1))

  
  model = make_model(output_bias=initial_bias)
 
  model.fit(
     X_train,
     y_train,
     batch_size=BATCH_SIZE,
     epochs=EPOCHS,
     callbacks=[early_stopping],
     validation_data=(X_val, y_val),
     class_weight=class_weight)

  test_predictions_baseline = model.predict(X_test, batch_size=BATCH_SIZE)
  def plot_cm(labels, predictions, p=0.5):
   cm = confusion_matrix(labels, predictions > p)
   plt.figure(figsize=(5,5))
   sns.heatmap(cm, annot=True, fmt="d")
   plt.title('Confusion matrix @{:.2f}'.format(p))
   plt.ylabel('Actual label')
   plt.xlabel('Predicted label')
  
   print('(True Negatives): ', cm[0][0])
   print('(False Positives): ', cm[0][1])
   print('(False Negatives): ', cm[1][0])
   print('(True Positives): ', cm[1][1])
   print('Total GOL: ', np.sum(cm[1]))
   print('Examples:\n    Total: {}\n    Positive: {} ({:.2f}% of total)\n'.format(
   total, pos, 100 * pos / total))
   print('Number of json files read: ', n_json)
   print('Number of success: ', n_success)
   print('Number of failues', len(list_fail))
   print("neural network")
  baseline_results = model.evaluate(X_test, y_test,
                                   batch_size=BATCH_SIZE, verbose=0)
  for name, value in zip(model.metrics_names, baseline_results):
   print(name, ': ', value)
 
  print()
 
  
  plot_cm(y_test,test_predictions_baseline)#

  return


if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print ('\nFinished. Time:', round(time.time()-t0, 2))
