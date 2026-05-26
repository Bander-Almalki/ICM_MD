import pandas as pd
from tensorflow.keras.models import load_model
import tensorflow as tf
from tensorflow.keras.models import load_model
from tensorflow.keras import backend as K
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score
import math



def focal_loss_fixed(y_true, y_pred):
    gamma = 2.0
    alpha = 0.25

    y_true = tf.cast(y_true, tf.float32)
    y_pred = K.clip(y_pred, K.epsilon(), 1.0 - K.epsilon())

    loss = -alpha * y_true * K.pow(1.0 - y_pred, gamma) * K.log(y_pred) \
           - (1.0 - alpha) * (1.0 - y_true) * K.pow(y_pred, gamma) * K.log(1.0 - y_pred)

    return K.mean(loss)


# Define the focal loss function
def focal_loss(gamma=2.0, alpha=0.25):
    def focal_loss_fixed(y_true, y_pred):
        epsilon = tf.keras.backend.epsilon()
        y_pred = tf.clip_by_value(y_pred, epsilon, 1 - epsilon)
        y_true = tf.cast(y_true, tf.float32)

        alpha_t = y_true * alpha + (1 - y_true) * (1 - alpha)
        p_t = y_true * y_pred + (1 - y_true) * (1 - y_pred)
        focal_loss = -alpha_t * tf.pow((1 - p_t), gamma) * tf.math.log(p_t)

        return tf.reduce_mean(focal_loss)
    return focal_loss_fixed


def precision_at(y_true, y_scores, L):
    # Sort predictions by their scores, descending
    sorted_indices = np.argsort(y_scores)[::-1]
    top_L_indices = sorted_indices[:L]

    # Select top L true labels and count relevant items
    top_L_true_labels = y_true[top_L_indices]
    top_L_predicted_labels = y_scores[top_L_indices]
    top_L_predicted_labels=np.ones_like(top_L_predicted_labels)
    precision=precision_score(top_L_true_labels,top_L_predicted_labels)
    return precision

test_dir="Features/test/"
test_set_df=pd.read_csv(test_dir+"test_set_df.csv")
test_set_df2=test_set_df.copy()
test_set_df2.drop(columns=['Dimer','residue1_num','residue1_name','res2_residue_num','res2_residue_name'],inplace=True)

x_test=test_set_df2.drop(columns=['min_distance'])
x_test.columns=x_test.columns.astype(str)
y_test=test_set_df2['min_distance']
y_test=y_test.apply(lambda x: 1 if x<=8 else 0)

model = load_model("model.keras", custom_objects={"focal_loss_fixed": focal_loss_fixed})
predictions = model.predict(x_test)



#prediction
pred=model.predict(x_test)

pred_binary=np.where(pred>0.5,1,0)


test_grouper=test_set_df.groupby('Dimer')
test_grouper.groups.keys()


#itterate over test grouper
test_pred_proba={}
true_labels={}
top_pred={}

for key,value in test_grouper.groups.items():
  L=int(math.sqrt(value.shape[0]))
  print(key, L)  
  v=test_grouper.get_group(key)
  res_nums=v[['residue1_num','res2_residue_num']]
  v=v.drop(columns=['Dimer','residue1_num','residue1_name','res2_residue_num','res2_residue_name'])
  x_test=v.drop(columns=['min_distance'])
  x_test = x_test.apply(pd.to_numeric, errors='coerce')
  pred=model.predict(x_test)
  pred_binary=np.where(pred>0.5,1,0)
  df=pd.concat([res_nums,pd.DataFrame(pred)],axis=1)
  df.columns=['residue1_num','res2_residue_num','pred']
  test_pred_proba[key]=df
  y_test=v['min_distance']
  y_test=y_test.apply(lambda x: 1 if x<=8 else 0)
  true_labels[key]=y_test
  top_10=precision_at(y_test.values,pred[:,0],10)
  top_25=precision_at(y_test.values,pred[:,0],25)
  top_l_5=precision_at(y_test.values,pred[:,0],round(L/5))
  #top_l_10=precision_at(y_test.values,pred[:,0],round(L/10))
  top_l=precision_at(y_test.values,pred[:,0],round(L))
  top_pred[key]=[top_10,top_25,top_l_5,top_l]
  #print(key,f1_score(y_test,pred_binary),precision_score(y_test,pred_binary),recall_score(y_test,pred_binary))
  print(key,top_10,top_25,top_l_5,top_l)

