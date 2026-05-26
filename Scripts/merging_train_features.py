import os
import pandas as pd
import numpy as np
import pickle
import torch


with open("ccmpred_rpair.pickle", 'rb') as handle:
    train_ccmpred_rpair = pickle.load(handle)

with open ("row_att_rpair.pickle", 'rb') as handle:
    train_row_att_rpair = pickle.load(handle)

with open ("sxxxs_rpair.pickle", 'rb') as handle:
    train_sxxxs_rpair = pickle.load(handle)

with open ("pxxxp_rpair.pickle", 'rb') as handle:
    train_pxxxp_rpair = pickle.load(handle)

with open ("gxxxg_rpair.pickle", 'rb') as handle:
    train_gxxxg_rpair = pickle.load(handle)

with open ("monomers_distance_map_from_dimer.pickle", 'rb') as handle:
    train_monomer_intra_dist_map_rpair = pickle.load(handle)

with open ("pssm_rpair.pickle", 'rb') as handle:
    train_pssm_rpair = pickle.load(handle)

with open ("lips_rpair.pickle", 'rb') as handle:
    train_lips_rpair = pickle.load(handle)

with open ("new_train_dimers_distance_map.pkl", 'rb') as handle:
    train_labels = pickle.load(handle)

print(train_ccmpred_rpair["TACD2_HUMAN"].shape)
print(train_row_att_rpair["TACD2_HUMAN"].shape)
print(train_sxxxs_rpair["TACD2_HUMAN"].shape)
print(train_pxxxp_rpair["TACD2_HUMAN"].shape)
print(train_gxxxg_rpair["TACD2_HUMAN"].shape)
# print(train_monomer_intra_dist_map_rpair["TACD2_HUMAN"].shape)
print(train_pssm_rpair["TACD2_HUMAN"].shape)
print(train_lips_rpair["TACD2_HUMAN"].shape)



train_set_df=pd.DataFrame()
for key in train_labels.keys():
  df=pd.DataFrame()
  #add a column containing the key
  k=[key]*train_labels[key].shape[0]
  df=pd.concat([df,pd.DataFrame(k,columns=['Dimer'])],axis=1)
  df=pd.concat([df,train_lips_rpair[key].reset_index(drop=True)],axis=1)

  gxxx=pd.DataFrame(train_gxxxg_rpair[key])
  gxxx=gxxx.reset_index(drop=True,inplace=False)
  df=pd.concat([df,gxxx],axis=1)

  pxxxp=pd.DataFrame(train_pxxxp_rpair[key])
  pxxxp=pxxxp.reset_index(drop=True,inplace=False)
  df=pd.concat([df,pxxxp],axis=1)

  sxxxs=pd.DataFrame(train_sxxxs_rpair[key])
  sxxxs=sxxxs.reset_index(drop=True,inplace=False)
  df=pd.concat([df,sxxxs],axis=1)

  m_distance=pd.DataFrame(train_monomer_intra_dist_map_rpair[key],columns=['res1_num','res2_num','res1_count','res2_count','intra_distance'])
  m_distance=m_distance['intra_distance']
  df=pd.concat([df,m_distance],axis=1)

  df=pd.concat([df,pd.DataFrame(train_pssm_rpair[key].reset_index(drop=True))],axis=1)
  df=pd.concat([df,pd.DataFrame(train_ccmpred_rpair[key])],axis=1)
  df=pd.concat([df,pd.DataFrame(train_row_att_rpair[key])],axis=1)

  #concatenate the train label
  df=pd.concat([df,pd.DataFrame(train_labels[key]['min_distance'])],axis=1)
  train_set_df=pd.concat([train_set_df,df])
  print(f"{key} processed")

print(train_set_df.shape)

train_set_df.to_csv("train_set.csv",index=False)

