import  pandas as pd

data = pd.read_table('mnx_reactions.tsv')
print(data.columns)
filter_data = data[['id','mnx_equation']]
print(filter_data.shape)
print(filter_data.head())
filter_data.to_csv('mnx_equation_filter.csv',index=False)