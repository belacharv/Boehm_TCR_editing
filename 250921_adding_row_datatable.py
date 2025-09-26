import pandas as pd
# Test with known parameters
#true_Pc = 0.9151860764850289 #0.8311920832252657 
#i = 1
pc = 0.91
#i = 1
true_n4 = 372 #107
true_n8 = 14 #7
total_seq = true_n4 + true_n8
ci = []
df_ci = pd.DataFrame(ci)

#model = TcellLineageModel(pc, true_n4, total_seq,"a")
#print_results(pc,true_n4,i)
def functino(pc):
    CIs = []
    CIs.append({
            'mean_pc':pc,
            'ci_pc_1':1,
            'ci_pc_2':2,
            'mean_n4':3,
            'ci_n4_1':4,
            'ci_n4_2':5,
            'mean_cd4_cd4':6,
            'ci_44_1':7,
            'ci_44_2':8,

        })
    df_CIs = pd.DataFrame(CIs)
    return df_CIs


for i in range(1,12):
    print(f"Simulations of Pc = {pc}")
    print("="*60)
    df_ci_new_row = functino(pc)
    df_ci = pd.concat([df_ci,df_ci_new_row],ignore_index=True)
    pc = pc - 0.01
    
    
#df_ci = pd.DataFrame(ci)
print(df_ci)
df_ci.to_csv("250921_ci_simulations_table_new.csv", index=False)