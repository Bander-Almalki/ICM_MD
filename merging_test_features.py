import pickle
import numpy as np
import pandas as pd

test_dir="Features/test/features_pickle/"

with open (test_dir+"test_ccmpred_rpair.pickle", "rb") as f:
  test_ccmpred = pickle.load(f)
with open (test_dir+"test_lips_rpair.pickle", "rb") as f:
  test_lips = pickle.load(f)
with open (test_dir+"test_gxxxg_rpair.pickle", "rb") as f:
  test_gxxxg = pickle.load(f)
with open(test_dir+"test_row_att_rpair.pickle", "rb") as f:
  test_msa_row_att = pickle.load(f)
with open(test_dir+"test_pxxxp_rpair.pickle", "rb") as f:
  test_pxxxp = pickle.load(f)
with open(test_dir+"test_sxxxs_rpair.pickle", "rb") as f:
  test_sxxxs = pickle.load(f)
with open(test_dir+"test_pssm_rpair.pickle", "rb") as f:
  test_PSSM = pickle.load(f)
with open(test_dir+"test_monomers_distance_map_from_dimer.pickle", "rb") as f:
  test_monomer_distances = pickle.load(f)
# importing test lables
with open (test_dir+"test_dimers_distance_map.pkl", "rb") as f:
  test_labels = pickle.load(f)



test_tm_domain_start_end={
    "2hac":[4,24],

    "2k1k":[13,35],
    "2l34":[5,26],
    "2loh":[18,43],
    "2l6w":[10,33],
    "2lzl":[18,41],
    "2mk9":[4,27],
    "2ka1":[9,34],
    "2l2t":[10,34],
    "2n90":[12,35]
}
tmdf=pd.DataFrame(test_monomer_distances['2ka1'],columns=['res1_num','res2_num','res1_count','res2_count','intra_distance'])

#get the dictionary keys as a list
test_keys = list(test_lips.keys())

new_test_monomer_distances={}
for key,value in test_monomer_distances.items():
  #choosing only specific test proteins
  if key not in test_tm_domain_start_end.keys():
    continue
  tmdf=pd.DataFrame(value,columns=['res1_num','res2_num','res1_count','res2_count','intra_distance'])
  dr=tmdf[~((tmdf['res1_count'] <test_tm_domain_start_end[key][0]-1) | (tmdf['res2_count'] <test_tm_domain_start_end[key][0]-1))]
  dr.reset_index(drop=True,inplace=True)
  dr2=dr[~((dr['res1_count'] >test_tm_domain_start_end[key][1]-2) | (dr['res2_count'] >test_tm_domain_start_end[key][1]-2))]
  dr2.reset_index(drop=True,inplace=True)
  new_test_monomer_distances[key]=dr2

  def capitalize_keys(d):
    return {k.upper(): v for k, v in d.items()}

# my_dict = {'name': 'Alice', 'age': 30}
# new_dict = capitalize_keys(my_dict)
# print(new_dict)  # Output: {'NAME': 'Alice', 'AGE': 30}
test_ccmpred=capitalize_keys(test_ccmpred)
test_lips=capitalize_keys(test_lips)
test_gxxxg=capitalize_keys(test_gxxxg)
test_msa_row_att=capitalize_keys(test_msa_row_att)
test_pxxxp=capitalize_keys(test_pxxxp)
test_sxxxs=capitalize_keys(test_sxxxs)
test_PSSM=capitalize_keys(test_PSSM)
new_test_monomer_distances=capitalize_keys(new_test_monomer_distances)
test_monomer_distances=capitalize_keys(test_monomer_distances)
test_labels=capitalize_keys(test_labels)

def capitalize_list(lst):
    
    return [x.upper() for x in lst]
two_chains_dimers=["2hac","2k1k","2l34","2loh","2l6w","2lzl","2ka1","2l2t"]
two_chains_dimers=capitalize_list(two_chains_dimers)
test_keys=capitalize_list(test_keys)

test_set_df=pd.DataFrame()
for key in test_keys:
  #print(f"{key} exists")
  if key not in two_chains_dimers:
    continue
  df=pd.DataFrame()
  #add a column containing the key
  k=[key]*test_ccmpred[key].shape[0]
  df=pd.concat([df,pd.DataFrame(k,columns=['Dimer'])],axis=1)
  df=pd.concat([df,test_lips[key].reset_index(drop=True)],axis=1)

  gxxx=pd.DataFrame(test_gxxxg[key])
  gxxx=gxxx.reset_index(drop=True,inplace=False)
  df=pd.concat([df,gxxx],axis=1)

  pxxxp=pd.DataFrame(test_pxxxp[key])
  pxxxp=pxxxp.reset_index(drop=True,inplace=False)
  df=pd.concat([df,pxxxp],axis=1)

  sxxxs=pd.DataFrame(test_sxxxs[key])
  sxxxs=sxxxs.reset_index(drop=True,inplace=False)
  df=pd.concat([df,sxxxs],axis=1)

  #concatenate monomer distance
  m_distance=new_test_monomer_distances[key]['intra_distance']
  df=pd.concat([df,m_distance],axis=1)

  df=pd.concat([df,pd.DataFrame(test_PSSM[key].reset_index(drop=True))],axis=1)
  df=pd.concat([df,pd.DataFrame(test_ccmpred[key])],axis=1)
  df=pd.concat([df,pd.DataFrame(test_msa_row_att[key])],axis=1)


  #concatenate the test label
  df=pd.concat([df,test_labels[key]['min_distance']],axis=1)
  test_set_df=pd.concat([test_set_df,df])

print("saving test_set_df")
test_set_df.to_csv("test_set_df.csv",index=False)