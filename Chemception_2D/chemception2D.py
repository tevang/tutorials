import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')    # temporary solution to avoid "ImportError: No module named PyQt5" which is mainly for Python 3
import matplotlib.pyplot as plt
print "RDKit: %s"%rdkit.__version__


import keras
from keras.models import Sequential, Model
from keras.layers import Conv2D, MaxPooling2D, Input, GlobalMaxPooling2D
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.optimizers import Adam
from keras.preprocessing.image import ImageDataGenerator
from keras.callbacks import ReduceLROnPlateau
print("Keras: %s"%keras.__version__)


data = pd.read_hdf("Sutherland.h5","table")
data["mol"] = data["smiles"].apply(Chem.MolFromSmiles)


def chemcepterize_mol(mol, embed=20.0, res=0.5):
    dims = int(embed*2/res)
    cmol = Chem.Mol(mol.ToBinary())
    cmol.ComputeGasteigerCharges()
    AllChem.Compute2DCoords(cmol)
    coords = cmol.GetConformer(0).GetPositions()
    vect = np.zeros((dims,dims,4))
    #Bonds first
    for i,bond in enumerate(mol.GetBonds()):
        bondorder = bond.GetBondTypeAsDouble()
        bidx = bond.GetBeginAtomIdx()
        eidx = bond.GetEndAtomIdx()
        bcoords = coords[bidx]
        ecoords = coords[eidx]
        frac = np.linspace(0,1,int(1/res*2)) #
        for f in frac:
            c = (f*bcoords + (1-f)*ecoords)
            idx = int(round((c[0] + embed)/res))
            idy = int(round((c[1]+ embed)/res))
            #Save in the vector first channel
            vect[ idx , idy ,0] = bondorder
    #Atom Layers
    for i,atom in enumerate(cmol.GetAtoms()):
            idx = int(round((coords[i][0] + embed)/res))
            idy = int(round((coords[i][1]+ embed)/res))
            #Atomic number
            vect[ idx , idy, 1] = atom.GetAtomicNum()
            #Gasteiger Charges
            charge = atom.GetProp("_GasteigerCharge")
            vect[ idx , idy, 3] = charge
            #Hybridization
            hyptype = atom.GetHybridization().real
            vect[ idx , idy, 2] = hyptype
    return vect


# To better understand what the code has done, lets try to “chemcepterize” a molecule and show it as an image.
# The embedding and the resolution are set lower than they will be for the final dataset. Matplotlib only supports
# RGB, so only the first three channels are used.
mol = data["mol"][0]
v = chemcepterize_mol(mol, embed=10, res=0.2)
print(v.shape)
plt.imshow(v[:,:,:3])


# Next step is to “chemcepterize” the entire collection of RDKit molecules and add a new column with the “images” to the dataframe
def vectorize(mol):
    return chemcepterize_mol(mol, embed=12)
data["molimage"] = data["mol"].apply(vectorize)

# The dataset already had a split value indicating if it should be train or test set. The shape of the final numpy arrays are
# (samples, height, width, channels)
X_train = np.array(list(data["molimage"][data["split"]==1]))
X_test = np.array(list(data["molimage"][data["split"]==0]))
print(X_train.shape)
print(X_test.shape)


# We also need to the prepare the values to predict. Here it is the IC50 for some DHFR inhibitors. The data is converted to log space and
# the robust scaler from scikit-learn is used to scale the data to somewhat between -1 and 1 (neural networks like this range and it makes
# training somewhat easier).
assay = "PC_uM_value"
y_train = data[assay][data["split"]==1].values.reshape(-1,1)
y_test = data[assay][data["split"]==0].values.reshape(-1,1)
from sklearn.preprocessing import RobustScaler
rbs = RobustScaler(with_centering=True, with_scaling=True, quantile_range=(5.0, 95.0), copy=True)
y_train_s = rbs.fit_transform(np.log(y_train))
y_test_s = rbs.transform(np.log(y_test))
h = plt.hist(y_train_s, bins=20)


input_shape = X_train.shape[1:]
print input_shape


def Inception0(input):
    tower_1 = Conv2D(16, (1, 1), padding='same', activation='relu')(input)
    tower_1 = Conv2D(16, (3, 3), padding='same', activation='relu')(tower_1)
    tower_2 = Conv2D(16, (1, 1), padding='same', activation='relu')(input)
    tower_2 = Conv2D(16, (5, 5), padding='same', activation='relu')(tower_2)
    tower_3 = Conv2D(16, (1, 1), padding='same', activation='relu')(input)
    output = keras.layers.concatenate([tower_1, tower_2, tower_3], axis=-1)
    return output


def Inception(input):
    tower_1 = Conv2D(16, (1, 1), padding='same', activation='relu')(input)
    tower_1 = Conv2D(16, (3, 3), padding='same', activation='relu')(tower_1)
    tower_2 = Conv2D(16, (1, 1), padding='same', activation='relu')(input)
    tower_2 = Conv2D(16, (5, 5), padding='same', activation='relu')(tower_2)
    tower_3 = MaxPooling2D((3, 3), strides=(1, 1), padding='same')(input)
    tower_3 = Conv2D(16, (1, 1), padding='same', activation='relu')(tower_3)
    output = keras.layers.concatenate([tower_1, tower_2, tower_3], axis=-1)
    return output


input_img = Input(shape=input_shape)
x = Inception0(input_img)
x = Inception(x)
x = Inception(x)
od=int(x.shape[1])
x = MaxPooling2D(pool_size=(od,od), strides=(1,1))(x)
x = Flatten()(x)
x = Dense(100, activation='relu')(x)
output = Dense(1, activation='linear')(x)
model = Model(inputs=input_img, outputs=output)
print model.summary()


# For the optimization I use the Adam optimizer and the mean absolute error as a loss function.
optimizer = Adam(lr=0.00025)
model.compile(loss="mae", optimizer=optimizer)


# The next part is crucial to avoid overfitting. Here the ImageDataGenerator object is used to perform random rotations and flips
# of the images before the training as a way of augmenting the training dataset. By doing this, the network will learn how to handle
# rotations and seeing the features in different orientations will help the model generalize better. Not including this will lead to
# completely overfit models. We have not encoded stereochemical information in the images, otherwise the flipping should be done by
# other means. The training set is concatenated to 50 times the length to have some sensible size epochs.

from image import ImageDataGenerator
generator = ImageDataGenerator(rotation_range=180,
                               width_shift_range=0.1,height_shift_range=0.1,
                               fill_mode="constant",cval = 0,
                               horizontal_flip=True, vertical_flip=True,data_format='channels_last',
                               )
#Concatenate for longer epochs
Xt = np.concatenate([X_train]*50, axis=0)
yt = np.concatenate([y_train_s]*50, axis=0)
batch_size=128
g = generator.flow(Xt, yt, batch_size=batch_size, shuffle=True)
steps_per_epoch = 10000/batch_size



# Now for the interesting part: Training. To lower the learning rate once the validation loss starts to plateau off I use
# the ReduceLROnPlateau callback avaible as part of Keras.
reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.5,patience=10, min_lr=1e-6, verbose=1)
history = model.fit_generator(g,
                              steps_per_epoch=len(Xt)//batch_size,
                              epochs=150,
                              validation_data=(X_test,y_test_s),
                              callbacks=[reduce_lr])

# Models can be saved and loaded. The history objects history dictionary is pickled.
name = "Chemception_std_notebook_demo"
model.save("%s.h5"%name)
hist = history.history
import pickle
pickle.dump(hist, file("%s_history.pickle"%name,"w"))
#from keras.model import load_model
#model = load_model("%s.h5"%name)


# The convergence of the training can be judged from a plot of the learning process. Somewhat unusual, when there's
# no regularization: The validation loss drops before the loss. The validation set is not augmented and thus consists of
# some “perfect” pictures, whereas maybe it may take the network some longer to deal with all the rotations, which also
# introduces some pixel artifacts due to the low resolution.
for label in ['val_loss','loss']:
    plt.plot(hist[label], label = label)
plt.legend()
plt.yscale("log")
plt.xlabel("Epochs")
plt.ylabel("Loss/lr")

# Plotting and Evaluating the Performance
y_pred_t = rbs.inverse_transform(model.predict(X_train))
y_pred = rbs.inverse_transform(model.predict(X_test))
plt.scatter(np.log(y_train), y_pred_t, label="Train")
plt.scatter(np.log(y_test), y_pred, label="Test")
plt.xlabel("log(PC_uM)")
plt.ylabel("predicted")
plt.plot([-10,6],[-10,6])
plt.legend()

corr2 = np.corrcoef(np.log(y_test).reshape(1,-1), y_pred.reshape(1,-1))[0][1]**2
rmse = np.mean((np.log(y_test) - y_pred)**2)**0.5
print("R2 : %0.2F"%corr2)
print("RMSE : %0.2F"%rmse)


# Visualizing the Layers
# It can be interesting to try and understand how the model "sees" the molecules. For this I’ll take an example molecule
# and plot some of the outputs from the different layers. I’ve taken the compound with the lowest IC50, number 143 in the dataset.
molnum = 143
molimage = np.array(list(data["molimage"][molnum:molnum+1]))
mol = data["mol"][molnum]

# The molecule looks like this
from rdkit.Chem import Draw
Draw.MolToImage(mol)

# And has this “chemcepterized” image as shown below
plt.imshow(molimage[0,:,:,:3])

# The first example is the third layer, which is the 1,1 convolution which feeds the 3,3 convolutional layer in tower 2.
layer1_model = Model(inputs=model.input,
                    outputs=model.layers[2].output)
kernels1 = layer1_model.predict(molimage)[0]
def plot_kernels(kernels):
    fig, axes = plt.subplots(2,3, figsize=(12,8))
    for i,ax in enumerate(axes.flatten()):
        ax.matshow(kernels[:,:,i])
        ax.set_title("Kernel %s"%i)
plot_kernels(kernels1)

# Lets go deeper...
for layer in [7,13,15,19,20]:
    print("Layer %i"%layer)
    plot_kernels(Model(inputs=model.input,outputs=model.layers[layer].output).predict(molimage)[0])
    plt.show()







