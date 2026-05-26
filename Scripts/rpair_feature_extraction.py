import os
import pandas as pd
import numpy as np
import pickle
import torch

#import train features from pickle file

# Extracting residue pair values of CCMpred
'''
ccmpred_rpair={}
CCMpred_dir="Features/train/CCMpred_output/"
for filename in os.listdir(CCMpred_dir):
    print(f"Processing {filename}")
    df=pd.read_csv(CCMpred_dir+filename,sep='\s+',header=None)
    ccm_rp_feat=df.melt(ignore_index=True)
    ccm_rp_feat.drop(columns=['variable'],inplace=True)
    ccmpred_rpair[filename.split(".fasta")[0]]=ccm_rp_feat
    print(f"Processing {filename}==> orig_shape: {df.shape}, new shape: {ccm_rp_feat.shape}")

print("Processing CCMpred files Done!")

#save ccmpred_rpair as a picle file 
with open('ccmpred_rpair.pickle', 'wb') as handle:
    pickle.dump(ccmpred_rpair, handle, protocol=pickle.HIGHEST_PROTOCOL)
'''

#----------------------------------------------------------------------------------------
'''

#Extracting residue pair values of row att
MSATransformer_file="Features/train/MSAT_output/msat_row_att.pickle"

with open(MSATransformer_file, 'rb') as handle:
    new_train_msa_row_att = pickle.load(handle)
row_att_rpair={}
for key,value in new_train_msa_row_att.items():
  row_att_rp_feat=value
  row_att_rp_feat=row_att_rp_feat.cpu()
  row_att_rp_feat=row_att_rp_feat.numpy().reshape(-1, 144)
  row_att_rp_feat=pd.DataFrame(row_att_rp_feat)
  row_att_rpair[key.split(".fasta")[0]]=row_att_rp_feat
print("Done!")
print(row_att_rpair.keys())
#print(row_att_rpair['CD72_HUMAN'].shape)

with open('row_att_rpair.pickle', 'wb') as handle:
    pickle.dump(row_att_rpair, handle, protocol=pickle.HIGHEST_PROTOCOL)
'''

#----------------------------------------------------------------------------------------
#lips score rpair extraction

'''
lips_rpair={}
lips_dir="Features/train/lips_outputs/"

for file in os.listdir(lips_dir):
    if file.endswith(".csv"):
        df2=pd.DataFrame()
        print(f"Processing {file}")
        df=pd.read_csv(lips_dir+file)
        for ind,row in df.iterrows():
            for ind2,row2 in df.iterrows():
                conc=np.concatenate([row.values,row2.values])
                df1=pd.DataFrame(conc.reshape(1,-1))
                df2=pd.concat([df2,df1],axis=0)
        df2.columns=['residue1_num', 'residue1_name', 'res1_LIPS_polarity', 'res1_LIPS_entropy',
                     'res1_LIPS_surface','res2_residue_num', 'res2_residue_name', 'res2_LIPS_polarity', 'res2_LIPS_entropy',
                     'res2_LIPS_surface']
        lips_rpair[file.split(".csv")[0]]=df2

with open('lips_rpair.pickle', 'wb') as handle:
    pickle.dump(lips_rpair, handle, protocol=pickle.HIGHEST_PROTOCOL)

'''



#----------------------------------------------------------------------------------------

#Extracting residue pair values of gxxxg

'''
gxxxg_rpair={}

gxxxg_dir="Features/train/motifs_outputs/gxxxg/"

for file in os.listdir(gxxxg_dir):
    if file.endswith(".csv"):
        df1=pd.DataFrame()
        print(f"Processing {file}")
        df=pd.read_csv(gxxxg_dir+file,index_col=0)
        for i in range (df.shape[0]):
            for j in range (df.shape[0]):
                hstack=np.hstack([df.iloc[i].values,df.iloc[j].values])
                df1=pd.concat([df1,pd.DataFrame(hstack.reshape(1,-1))],axis=0)
                
            #df1.columns=['residue1_motif','residue2_motif']
            gxxxg_rpair[file.split(".csv")[0]]=df1
with open('gxxxg_rpair.pickle', 'wb') as handle:
    pickle.dump(gxxxg_rpair, handle, protocol=pickle.HIGHEST_PROTOCOL)

'''
#----------------------------------------------------------------------------------------

#Extracting residue pair values of sxxxs

'''
sxxxs_rpair={}

sxxxs_dir="Features/train/motifs_outputs/sxxxs/"

for file in os.listdir(sxxxs_dir):
    if file.endswith(".csv"):
        df1=pd.DataFrame()
        print(f"Processing {file}")
        df=pd.read_csv(sxxxs_dir+file,index_col=0)
        for i in range (df.shape[0]):
            for j in range (df.shape[0]):
                hstack=np.hstack([df.iloc[i].values,df.iloc[j].values])
                df1=pd.concat([df1,pd.DataFrame(hstack.reshape(1,-1))],axis=0)
                
            #df1.columns=['residue1_motif','residue2_motif']
            sxxxs_rpair[file.split(".csv")[0]]=df1
with open('sxxxs_rpair.pickle', 'wb') as handle:
    pickle.dump(sxxxs_rpair, handle, protocol=pickle.HIGHEST_PROTOCOL)
'''





#----------------------------------------------------------------------------------------

#Extracting residue pair values of pxxxp
 


# pxxxp_rpair={}

# pxxxp_dir="Features/train/motifs_outputs/pxxxp/"

# for file in os.listdir(pxxxp_dir):
#     if file.endswith(".csv"):
#         df1=pd.DataFrame()
#         print(f"Processing {file}")
#         df=pd.read_csv(pxxxp_dir+file,index_col=0)
#         for i in range (df.shape[0]):
#             for j in range (df.shape[0]):
#                 hstack=np.hstack([df.iloc[i].values,df.iloc[j].values])
#                 df1=pd.concat([df1,pd.DataFrame(hstack.reshape(1,-1))],axis=0)
                
#             #df1.columns=['residue1_motif','residue2_motif']
#             pxxxp_rpair[file.split(".csv")[0]]=df1
# with open('pxxxp_rpair.pickle', 'wb') as handle:
#     pickle.dump(pxxxp_rpair, handle, protocol=pickle.HIGHEST_PROTOCOL)


#----------------------------------------------------------------------------------------

#Extracting residue pair values of pssm 

'''
pssm_rpair={}

pssm_dir="Features/train/pssm_outputs/"

for file in os.listdir(pssm_dir):
    if file.endswith(".csv"):
        df=pd.read_csv(pssm_dir+file,index_col=0)
        df1=pd.DataFrame()
        for i in range (df.shape[0]):
            for j in range (df.shape[0]):
                hstack=np.hstack([df.iloc[i].values,df.iloc[j].values])
                df1=pd.concat([df1,pd.DataFrame(hstack.reshape(1,-1))],axis=0)
    print(f"Processing {file}==> orig_shape: {df.shape}, new shape: {df1.shape}")
    pssm_rpair[file.split(".csv")[0]]=df1

with open('pssm_rpair.pickle', 'wb') as handle:
    pickle.dump(pssm_rpair, handle, protocol=pickle.HIGHEST_PROTOCOL)

'''

#----------------------------------------------------------------------------------------

#Extracting residue pair values of monomer_intra_dist_map

 # Extracting residue pair values of CCMpred

'''
monomer_intra_dist_map_rpair={}
monomer_intra_dist_map_dir="Features/train/monomer_features_outputs/"

for file in os.listdir(monomer_intra_dist_map_dir):
    if file.endswith(".csv"):
        print(f"processing{file}")
        df=pd.read_csv(monomer_intra_dist_map_dir+file,header=None)
        monom_rp_feat=df.melt(ignore_index=True)
        monom_rp_feat.drop(columns=['variable'],inplace=True)
        monomer_intra_dist_map_rpair[file.split(".csv")[0]]=monom_rp_feat
        print(f"Processing {file}==> orig_shape: {df.shape}, new shape: {monom_rp_feat.shape}")

with open('monomer_intra_dist_map_rpair.pickle', 'wb') as handle:
    pickle.dump(monomer_intra_dist_map_rpair, handle, protocol=pickle.HIGHEST_PROTOCOL)

'''





