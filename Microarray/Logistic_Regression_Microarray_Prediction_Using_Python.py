import pandas as pd 
import numpy as np
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


# ### Load Data

def load_data():
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

def preprocess_data(X_val):

    scaler2=pickle.load(open('MODELS/Microarray_StandardScaler.sav', 'rb'))

    
    X_val_norm=scaler2.transform(X_val)
    return X_val_norm


# ### Load Model

def load_model():
    model=pickle.load(open('MODELS/Microarray_LR.sav', 'rb'))

    return model


# ### Predict

def predict_data(model,X_val_norm,val_df):
    df=pd.DataFrame(columns=['Sample','Predicted_Result','Prob@AML'])
    df['Sample']=val_df.columns
    df['Predicted_Result']=model.predict(X_val_norm)
    df['Prob@AML']=model.predict_proba(X_val_norm)[:,1]
    df.to_csv('OUTPUT_FILES/Microarray_LR_Prediction.csv',index=False,sep='\t')
    return df



if __name__ == "__main__":
    
    print('Load Data .....')
    X_val,val_df=load_data()
    
    print('Preprocess Data .....')
    X_val_norm=preprocess_data(X_val)
    
    print('Load Logistic Regression Model .....')
    model=load_model()
    
    print('Predict Data .....')
    df=predict_data(model,X_val_norm,val_df)

