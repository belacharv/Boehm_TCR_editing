import pandas as pd
#sol = []
df_sol = pd.DataFrame([])
for i in range(5):
    sol3 = [{
        'n4': 1,
        'value': 2,
        'depth':i,
        'equation':3}]
    df_new = pd.DataFrame(sol3)
    df_sol = pd.concat([df_sol,df_new],ignore_index=True)

print(df_sol)
for i in range(5):
    n4,value,depth = df_sol.iloc[0,0],df_sol.iloc[0,1],df_sol.iloc[0,2]
    print(n4,value,depth)

row = df_sol.iloc[0]
print(row.iloc[0])