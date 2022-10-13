import pandas as pd 
import numpy as np
import pickle
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from umap import UMAP
import matplotlib.patches as mpatches


# ### Load Training Data

def load_training_data():
    bm=pd.read_csv('TRAINING_DATA/BM_Microarray_Selected_INTENSITY.csv')
    bm.drop(['SYMBOL'],axis=1,inplace=True)

    kg1a=pd.read_csv('TRAINING_DATA/KG1A_Microarray_Selected_INTENSITY.csv')
    kg1a.drop(['SYMBOL'],axis=1,inplace=True)

    bm_kg1a=pd.concat([bm,kg1a],axis=1)
    train_df=bm_kg1a
    
    X_train=train_df.values
    X_train=np.transpose(X_train)
    
    return X_train


# ### Load Validation Data

def load_validation_data():
    val_df=pd.read_csv('INPUT_FILES/Input_INTENSITIES.txt',sep='\t')

    if val_df.shape[0]>=25000:
        file_name='SELECTED_GENES_PROBES/Selected_Probes_HG-U133_Plus_2.csv'
    else:
        file_name='SELECTED_GENES_PROBES/Selected_Probes_HG-U133A.csv'

    selec_probes=pd.read_csv(file_name)
    val_df.rename({'Unnamed: 0':'PROBEID'},axis=1,inplace=True)
    val_df=pd.merge(val_df,selec_probes,on='PROBEID',how='inner')
    val_df.drop(['PROBEID'],axis=1,inplace=True)

    gene_sym=pd.read_csv('SELECTED_GENES_PROBES/Microarray_Selected_Genes.csv')
    val_df=val_df[val_df.SYMBOL.isin(gene_sym.SYMBOL.values)]
    val_df=pd.DataFrame(val_df.groupby(['SYMBOL']).mean()).reset_index()
    val_df.drop(['SYMBOL'],axis=1,inplace=True)
    
    X_val=val_df.values
    X_val=np.transpose(X_val)
    
    return X_val,val_df


# ### Preprocess Dataset

def preprocess_data(X_train,X_val):
    scaler1=StandardScaler()
    X_train_norm=scaler1.fit_transform(X_train)

    scaler2=pickle.load(open('MODELS/Microarray_StandardScaler.sav', 'rb'))

    
    X_val_norm=scaler2.transform(X_val)
    X_val_norm=np.concatenate([X_train_norm,X_val_norm])
    return X_val_norm


# ### Load Model

def load_model(X_val_norm):

    umap1=pickle.load(open('MODELS/Microarray_UMAP_Models.sav', 'rb'))

    X_val_umap=umap1.transform(X_val_norm)
    return X_val_umap


# ### Plot Specific Points

def plot_specifics(X_val_umap,val_df):
    words=['BM_c','BM_c','BM_c','BM_c','BM_c','BM_A','BM_A','BM_A',
           'KG1A_c','KG1A_c','KG1A_c','KG1A_c','KG1A_c','KG1A_A','KG1A_A']

    val_words=list(val_df.columns)
    words=words+val_words


    plt.figure(figsize=(10,10))
    fsize=9
    for i,w in enumerate(words):
        if i<=4:
            plt.scatter(X_val_umap[i,0], X_val_umap[i,1], marker='x', color='red')
        if (i>=5)&(i<=7):
            plt.scatter(X_val_umap[i,0], X_val_umap[i,1], marker='o', color='blue')
        if(i>=8)&(i<=12):
            plt.scatter(X_val_umap[i,0], X_val_umap[i,1], marker='x', color='red')
        if(i>=13)&(i<=14):
            plt.scatter(X_val_umap[i,0], X_val_umap[i,1], marker='o', color='blue')
        if (i>=15):
            plt.scatter(X_val_umap[i,0], X_val_umap[i,1], marker='x', color='green')
            fsize=15

        plt.text(X_val_umap[i,0]+.01, X_val_umap[i,1]+.01, w, fontsize=fsize)
    plt.title('UMAP visualization with sample specification of Microarray Data')
    plt.show()


# ### Plot as a Whole

def plot_a_whole(X_val_umap,val_df):
    label_val=['red','red','red','red','red','blue','blue','blue',
              'red','red','red','red','red','blue','blue'
              ]

    val_words=['green']*val_df.shape[1]
    label_val=label_val+val_words

    classes=['CD34+ Control','AML','New_Samples']
    fig, ax = plt.subplots(1, figsize=(10, 10))
    plt.scatter(*X_val_umap.T, s=70,c=np.array(label_val), alpha=1.0)

    class_colours = ['red','blue','green']
    recs = []
    for i in range(0,len(class_colours)):
        recs.append(mpatches.Rectangle((0,0),1,1,fc=class_colours[i]))
    plt.legend(recs,classes,loc=1,fontsize=15)
    plt.title('UMAP visualization of Microarray Data')
    plt.show()


if __name__ == "__main__":
    
    print('Load Training Data .....')
    X_train=load_training_data()
    
    print('Load Validation Data .....')
    X_val,val_df=load_validation_data()
    
    print('Preprocess Data .....')
    X_val_norm=preprocess_data(X_train,X_val)
    
    print('Load UMAP Model .....')
    X_val_umap=load_model(X_val_norm)
    
    plot_specifics(X_val_umap,val_df)
    
    plot_a_whole(X_val_umap,val_df)
