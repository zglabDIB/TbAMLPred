import pandas as pd
import numpy as np

check_df=pd.read_csv('ASSOCIATED_FILES/cd34_PB-AML_matrix.txt',sep='\t')

df=pd.read_csv('INPUT_FILES/ALL_Input_Count.txt',sep='\t')
df=df[df.Geneid.isin(check_df.Geneid.values)]

df.to_csv('INPUT_FILES/Input_Count.txt',sep='\t',index=False)
print(df.shape)
