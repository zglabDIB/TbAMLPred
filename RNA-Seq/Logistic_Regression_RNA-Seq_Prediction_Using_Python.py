import pandas as pd 
import numpy as np
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


# ### Load Data

def load_data():
    bm=pd.read_csv('TRAINING_DATA/BM_RNA_Seq_Selected_COUNT.csv')
    ensg_gene=bm.Gene.values
    
    val_df=pd.read_csv('INPUT_FILES/Input_Normalized_Count.txt',sep='\t')
    val_df=val_df[val_df.Gene.isin(ensg_gene)]
    val_df=val_df.set_index('Gene')
    val_df=val_df.loc[ensg_gene]
    val_df=val_df.reset_index()
    val_df.drop(['Gene'],axis=1,inplace=True)


    X_val=val_df.values
    X_val=np.transpose(X_val)
    
    return X_val,val_df


# ### Preprocess Dataset

def preprocess_data(X_val):
    scaler2=pickle.load(open('MODELS/RNA_Seq_StandardScaler.sav', 'rb'))

    
    X_val_norm=scaler2.transform(X_val)
    return X_val_norm


# ### Load Model

def load_model():
    model=pickle.load(open('MODELS/RNA-Seq_LR.sav', 'rb'))

    return model


# ### Predict

def predict_data(model,X_val_norm,val_df):
    df=pd.DataFrame(columns=['Sample','Predicted_Result','Prob@AML'])
    df['Sample']=val_df.columns
    df['Predicted_Result']=model.predict(X_val_norm)
    df['Prob@AML']=model.predict_proba(X_val_norm)[:,1]
    df.to_csv('OUTPUT_FILES/RNA-Seq_LR_Prediction.csv',index=False,sep='\t')
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
