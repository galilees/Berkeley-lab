import numpy as np
import pandas as pd
import json
import os, time
import random

path_to_json = '/net/anaconda/raid1/dorothee/14_frontiers_QR_restraints/galilee/pdb_survey/'

def data_frame(list_):
    df = pd.DataFrame(list_)
    return df

def run():
  '''
  '''
  os.chdir(path_to_json)
  n_json = 0
  n_success = 0
  list_fail = []
  nearby_dict_list_gol = []
  nearby_dict_list_hoh = []
  n_ss_helix_list      = []
  n_hbonds_list        = []
  for file_name in [file for file in os.listdir(path_to_json) if file.endswith('.json')]:
  #for root, dirs, files in os.walk(path_to_json):
    #if os.path.basename(root).startswith('queue_logs'): continue
    #for file in files:
    #  if not file.endswith('.json'): continue
      with open(file_name, 'r') as fp:
        n_json += 1
        #if n_json > 500: break
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
          if gol_data['n_hbonds'] is None: continue
          n_hbonds_list.append(gol_data['n_hbonds'])
          #n_ss_helix_list.append(gol_data['n_ss_helix'])
          nearby_dict_list_gol.append(gol_data['nearby_res'])
        #
        all_hohs_dict = data['HOH']
        for sel_str, hoh_data in all_hohs_dict.items():
          if 'nearby_res' not in hoh_data: continue
          nearby_dict_list_hoh.append(hoh_data['nearby_res'])

  print(len(n_hbonds_list))
  print(len(nearby_dict_list_gol))
  assert(len(n_hbonds_list) == len(nearby_dict_list_gol))

  print('Number of json files read: ', n_json)
  print('Number of success: ', n_success)
  print('Number of failues', len(list_fail))
  #print('failures')
  #for p in list_fail:
  #  print(p)

#print("Total number of glycerol: " , len(gol_data))

#gol_data_sample = gol_data
  df_gol = data_frame(list_ = nearby_dict_list_gol)
  #print(df_gol.head(10))
  df_hbonds = data_frame(n_hbonds_list)
  #print(" hbonds", df_hbonds)
#
  hoh_data_sample = random.sample(nearby_dict_list_hoh, 70000)
  df_hoh = data_frame(hoh_data_sample)
  #print(df_hoh.head(10))
#
  return
## All residue counts besides amino acids and HOH are combined into a single column called Other the dataframes are then called df_gol_aa and df_hoh_aa
#final_table_columns =['HOH','ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
#df_gol_sample_aa = df_gol[df_gol.columns.intersection(final_table_columns)]
#
#df_hoh_sample_aa = df_hoh[df_hoh.columns.intersection(final_table_columns)]
##print(df_gol_sample_aa)
##print(df_gol_sample_aa.columns)
#
#df_gol_sample_other = df_gol.drop(columns=[col for col in df_gol if col in final_table_columns])
#df_hoh_sample_other = df_hoh.drop(columns=[col for col in df_hoh if col in final_table_columns])
#
#df_gol_other = df_gol_sample_other.sum(axis=1)
#df_hoh_other = df_hoh_sample_other.sum(axis=1)
##print(df_gol_other )
#
#df_gol_aa = df_gol_sample_aa.assign(Other = df_gol_other)
##print(df_gol_aa)
#df_hoh_aa = df_hoh_sample_aa.assign(Other = df_hoh_other) #use to assign new random batch
##print(df_hoh_aa)
#
#
#df_gol_0 = df_gol_aa.fillna(0)
#df_hoh_0 = df_hoh_aa.fillna(0)
#
##print("columns" , df_gol_0.columns)
#
#electroneg = [{'HOH':7,'ARG':9,'HIS':6,'LYS':3,'ASP':7,'GLU':7,'SER':3.5,'THR':3.5,'ASN':6.5,'GLN':6.5,'CYS':2.5, 'GLY':0,'PRO':0,'ALA':0,'VAL':0,'ILE':0,'LEU':0,'MET':2.5,'PHE':0,'TYR':3.5,'TRP':3}]
#df_electroneg = pd.DataFrame(electroneg)
#df_electroneg
#df_repeated_gol = pd.concat([df_electroneg]*len(df_gol_0.index), ignore_index=True)
#
#df_gol_aa_eneg = df_gol_0.mul(df_repeated_gol, axis='columns', level=None, fill_value=None)
#
#
#df_repeated_hoh = pd.concat([df_electroneg]*len(df_hoh_0.index), ignore_index=True)
#
#df_hoh_aa_eneg = df_hoh_0.mul(df_repeated_hoh, axis='columns', level=None, fill_value=None)
#df_gol_aa_eneg = df_gol_aa_eneg.fillna(0)
#df_hoh_aa_eneg = df_hoh_aa_eneg.fillna(0)
##print(df_gol_aa_eneg)
#
#
#eneg_gol = df_gol_aa_eneg.sum(axis=1)
#eneg_hoh = df_hoh_aa_eneg.sum(axis=1)
#
#
##print(df_gol_0)
#
#df_gol_0_1 = df_gol_0.mask(df_gol_0 > 0, 1)
#df_hoh_0_1 = df_hoh_0.mask(df_hoh_0 > 0, 1)
##print(df_gol_0_1)
#
#print(df_gol_0.head(10))
#print(df_hoh_0.head(10))
#df_gol_sample_extended = df_gol_0_1.join(df_gol_0 ,how='inner',rsuffix="_count")
#df_hoh_sample_extended = df_hoh_0_1.join(df_hoh_0 ,how='inner',rsuffix="_count")
#df_gol_sample_extended['eneg_gol'] = eneg_gol
#df_hoh_sample_extended['eneg_hoh'] = eneg_hoh
##df_gol_sample_extended['h_bonds'] = df_hbonds
#
##print(df_gol_sample_extended)
##df_hoh_sample_extended.columns
#
#
#df_gol_sample = df_gol_sample_extended
#df_hoh_sample = df_hoh_sample_extended
#
#"Choose a random selection of negative samples. The sample indexes can be saved for repeated testing "
##df_hoh_sample = df_hoh_sample.sample(n=70000)
#
#df_hoh_sample_labels = df_hoh_sample.copy(deep = True)
#df_hoh_sample_labels['label'] = False
#
#df_gol_sample_labels = df_gol_sample.copy(deep = True)
#df_gol_sample_labels['label'] = True
##print(df_gol_sample_labels)
#
#
#features_df = pd.concat([df_gol_sample,df_hoh_sample])
#features_df = features_df.fillna(0)
##features_df = features_df.drop(columns=[], axis=1)
#
##The features are amino acid one hot encoded values, amino acid counts, electronegativity measures, and H_bonds to GOL
#print("Features for machine learning:", features_df.columns)
#print(features_df)
#
#Stop()
#
#labels_df = pd.concat([df_gol_sample_labels['label'],df_hoh_sample_labels['label']])
##print(X,y)
#
#X = np.array(features_df)
#
#y = np.array(labels_df)
#
##print(X,y)
#
#
## split into train test sets
#from sklearn.model_selection import train_test_split
#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33)
#
## Use if not splitting the data sets
##X_train, X_test, y_train, y_test = X, X, y, y
#
#
## Another thing to consider is imbalance, if you have too few positive samples then the model might not work well.
## Ideal would be 0.5, but having a larger number of negative is more realistic
#imbalance = y.sum()/y.shape[0] # ratio of positive samples over total samples
#print("Data imbalence:", imbalance)

#import tensorflow as tf
#from tensorflow import keras
#
#import os
#import tempfile
#
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import numpy as np
#import pandas as pd
#import seaborn as sns
#
#import sklearn
#from sklearn.metrics import confusion_matrix
#from sklearn.model_selection import train_test_split
#from sklearn.preprocessing import StandardScaler
#
#mpl.rcParams['figure.figsize'] = (12, 10)
#colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#
#METRICS = [
#      keras.metrics.TruePositives(name='tp'),
#      keras.metrics.FalsePositives(name='fp'),
#      keras.metrics.TrueNegatives(name='tn'),
#      keras.metrics.FalseNegatives(name='fn'),
#      keras.metrics.BinaryAccuracy(name='accuracy'),
#      keras.metrics.Precision(name='precision'),
#      keras.metrics.Recall(name='recall'),
#      keras.metrics.AUC(name='auc'),
#      keras.metrics.AUC(name='prc', curve='PR'), # precision-recall curve
#]
#
#def make_model(metrics=METRICS, output_bias=None):
#  if output_bias is not None:
#    output_bias = tf.keras.initializers.Constant(output_bias)
#  model = keras.Sequential([
#      keras.layers.Dense(
#          16, activation='relu',
#          input_shape=(X_train.shape[-1],)),
#      keras.layers.Dropout(0.5),
#      keras.layers.Dense(1, activation='sigmoid',
#                         bias_initializer=output_bias),
#  ])
#
#  model.compile(
#      optimizer=keras.optimizers.Adam(learning_rate=1e-3),
#      loss=keras.losses.BinaryCrossentropy(),
#      metrics=metrics)
#
#  return model
#
#EPOCHS = 100
#BATCH_SIZE = 1000
#
#early_stopping = tf.keras.callbacks.EarlyStopping(
#    monitor='val_prc',
#    verbose=1,
#    patience=10,
#    mode='max',
#    restore_best_weights=True)
#
#model = make_model()
#model.summary()
#
#model.predict(X_train[:10])
#
#results = model.evaluate(X_train, y_train, batch_size=BATCH_SIZE, verbose=0)
#print("Loss: {:0.4f}".format(results[0]))
#
#neg, pos = np.bincount(labels_df )
#total = neg + pos
#print('Examples:\n    Total: {}\n    Positive: {} ({:.2f}% of total)\n'.format(
#    total, pos, 100 * pos / total))
#
#initial_bias = np.log([pos/neg])
#initial_bias
#
#model = make_model(output_bias=initial_bias)
#model.predict(X_train[:10])
#
#results = model.evaluate(X_train, y_train, batch_size=BATCH_SIZE, verbose=0)
#print("Loss: {:0.4f}".format(results[0]))
#
#initial_weights = os.path.join(tempfile.mkdtemp(), 'initial_weights')
#model.save_weights(initial_weights)
#
#model = make_model()
#model.load_weights(initial_weights)
#model.layers[-1].bias.assign([0.0])
#zero_bias_history = model.fit(
#    X_train,
#    y_train,
#    batch_size=BATCH_SIZE,
#    epochs=20,
#    #validation_data=(val_features, val_labels),
#    verbose=0)
#
#model = make_model()
#model.load_weights(initial_weights)
#careful_bias_history = model.fit(
#    X_train,
#    y_train,
#    batch_size=BATCH_SIZE,
#    epochs=20,
#    #validation_data=(val_features, val_labels),
#    verbose=0)
#
#
#def plot_loss(history, label, n):
#  # Use a log scale on y-axis to show the wide range of values.
#  plt.semilogy(history.epoch, history.history['loss'],
#               color=colors[n], label='Train ' + label)
##   plt.semilogy(history.epoch, history.history['val_loss'],
##                color=colors[n], label='Val ' + label,
##                linestyle="--")
#  plt.xlabel('Epoch')
#  plt.ylabel('Loss')
#
#plot_loss(zero_bias_history, "Zero Bias", 0)
#plot_loss(careful_bias_history, "Careful Bias", 1)
#
#model = make_model()
#model.load_weights(initial_weights)
#baseline_history = model.fit(
#    X_train,
#    y_train,
#    batch_size=BATCH_SIZE,
#    epochs=EPOCHS,
#    callbacks=[early_stopping])
#    #, validation_data=(val_features, val_labels))
#def plot_metrics(history):
#  metrics = ['loss', 'prc', 'precision', 'recall']
#  for n, metric in enumerate(metrics):
#    name = metric.replace("_"," ").capitalize()
#    plt.subplot(2,2,n+1)
#    plt.plot(history.epoch, history.history[metric], color=colors[0], label='Train')
##     plt.plot(history.epoch, history.history['val_'+metric],
##              color=colors[0], linestyle="--", label='Val')
#    plt.xlabel('Epoch')
#    plt.ylabel(name)
#    if metric == 'loss':
#      plt.ylim([0, plt.ylim()[1]])
#    elif metric == 'auc':
#      plt.ylim([0.8,1])
#    else:
#      plt.ylim([0,1])
#
#    plt.legend();
#
#plot_metrics(baseline_history)
#
#train_predictions_baseline = model.predict(X_train, batch_size=BATCH_SIZE)
#test_predictions_baseline = model.predict(X_test, batch_size=BATCH_SIZE)
#
#
#def plot_cm(labels, predictions, p=0.5):
#  cm = confusion_matrix(labels, predictions > p)
#  plt.figure(figsize=(5,5))
#  sns.heatmap(cm, annot=True, fmt="d")
#  plt.title('Confusion matrix @{:.2f}'.format(p))
#  plt.ylabel('Actual label')
#  plt.xlabel('Predicted label')
#
#  print('(True Negatives): ', cm[0][0])
#  print('(False Positives): ', cm[0][1])
#  print('(False Negatives): ', cm[1][0])
#  print('(True Positives): ', cm[1][1])
#  print('Total GOL: ', np.sum(cm[1]))
#
#baseline_results = model.evaluate(X_test, y_test,
#                                  batch_size=BATCH_SIZE, verbose=0)
#for name, value in zip(model.metrics_names, baseline_results):
#  print(name, ': ', value)
#print()
#
#plot_cm(y_test,test_predictions_baseline)#


if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print ('\nFinished. Time:', round(time.time()-t0, 2))
