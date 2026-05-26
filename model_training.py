import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
data_folder="Features/train/"
train_set_df=pd.read_csv(data_folder+"train_set.csv")
train_set_df.drop(columns=['Dimer','residue1_num','residue1_name','res2_residue_num','res2_residue_name'],inplace=True)

x_train=train_set_df.drop(columns=['min_distance'])
x_train.columns=x_train.columns.astype(str)

y_train=train_set_df['min_distance']
# a thrashold on y_train
y_train=y_train.apply(lambda x: 1 if x<=8 else 0)


#------------------------------------*********--------------------------------------------

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

# Build a simple neural network classifier
model = Sequential()
model.add(Dense(100, input_dim=198, activation='relu'))  # Replace 20 with the number of features in your dataset
# add a dropout layer
model.add(tf.keras.layers.Dropout(0.2))
model.add(Dense(50, activation='relu'))
#add a dropout layer
model.add(tf.keras.layers.Dropout(0.2))
model.add(Dense(10, activation='relu'))
#add a dropout layer
# model.add(tf.keras.layers.Dropout(0.5))
# model.add(Dense(25, activation='relu'))

# model.add(tf.keras.layers.Dropout(0.5))
# model.add(Dense(10, activation='relu'))

model.add(Dense(1, activation='sigmoid'))

# Compile the model using focal loss
model.compile(optimizer=Adam(learning_rate=0.001), loss=focal_loss(gamma=1.5, alpha=0.25), metrics=['accuracy'])

# Train the model (replace X_train, y_train, X_val, y_val with your data)
# model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=20, batch_size=32)



model.fit(x_train, y_train, epochs=20, batch_size=32)
model.save_weights("model_weights.h5")

weights = model.get_weights()

# Save model weights to a file
with open("model_weights.txt", "w") as f:
    for layer_weights in weights:
        f.write(str(layer_weights) + "\n")

#save model 
model.save("model.keras")
