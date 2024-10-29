import os
import pandas as pd
data = pd.read_table('a/compound_new.txt')
print(data.columns)
data[['inchikey','entry','name']].to_csv('a/mini_compounds_new.csv',index=False)