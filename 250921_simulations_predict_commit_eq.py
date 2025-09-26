import numpy as np
#import matplotlib.pyplot as plt
import sympy as sp
import pandas as pd
from scipy import stats
from b_predict_commit_eq import TcellLineageModel
from b_predict_commit_eq import fit_model_to_data
import os
# Parameter recovery
os.chdir("changing_pc_value_simulations")
def simulate(pc, true_n4):
# Generate multiple datasets and see if we can recover parameters
    recovery_results = []
    n_simulations = 100

    for i in range(n_simulations):
    # Simulate data
        cd4_cd4, cd4_cd8, cd8_cd8 = model.simulate_dataset(seed=i)
        print(cd4_cd4, cd4_cd8, cd8_cd8)
    
    # Fit model to simulated data
        solutions = fit_model_to_data(cd4_cd4, cd8_cd8, total_seq,"a")
    
        if len(solutions) >= 1:
        # Take the first valid solution
            fitted_Pc, fitted_n4 = solutions[1]
            fitted_n8 = model.total_sequences - fitted_n4
        
            #print(f"n4: {fitted_n4}, n8: {fitted_n8}")
            #print(fitted_Pc, fitted_n4)
        
            recovery_results.append({
                'sim_id': i,
                'true_Pc': pc,
                'true_n4': true_n4,
                'fitted_Pc': fitted_Pc,
                'fitted_n4': fitted_n4,
                'cd4_cd4': cd4_cd4,
                'cd8_cd8': cd8_cd8,
                'error_Pc': abs(fitted_Pc - pc),
                'error_n4': abs(fitted_n4 - true_n4),
                'error_cd4_cd4': abs(cd4_cd4 - 337)
            })
        else:
            print("*"*50)
            print(f"No valid solutions found for simulation {i}")
            print("*"*50)
            print(cd4_cd4, cd4_cd8, cd8_cd8)
            print("*"*50)

    return recovery_results

def print_results(pc,true_n4,i):
    recovery_results = simulate(pc,true_n4)
    df_recovery = pd.DataFrame(recovery_results)
    n_simulations = 100
    print(df_recovery)
    print(f"Successful recoveries: {len(df_recovery)}/{n_simulations}")
    CIs = []

    if len(df_recovery) > 0:
        print(f"Mean Pc error: {df_recovery['error_Pc'].mean():.4f} ± {df_recovery['error_Pc'].std():.4f}")
        print(f"Mean n4 error: {df_recovery['error_n4'].mean():.4f} ± {df_recovery['error_n4'].std():.4f}")
        print(f"Mean cd4_cd4 error: {df_recovery['error_cd4_cd4'].mean():.4f} ± {df_recovery['error_cd4_cd4'].std():.4f}")
        name = str(i)+"_pc_"+str(pc)+"_simulations.csv"
        df_recovery.to_csv(name, index=False)

        fitted_Pc = df_recovery['fitted_Pc']

        confidence_level = 0.95
        mean = np.mean(fitted_Pc)
        sem = stats.sem(fitted_Pc)  # Standard error of the mean
        ci = stats.t.interval(confidence_level, len(fitted_Pc)-1, loc=mean, scale=sem)


        print(f"Mean: {mean}")
        print(f"95% Confidence Interval fitted_Pc: [{ci[0]:.8f}, {ci[1]:.8f}]")

        fitted_n4 = df_recovery['fitted_n4']

        confidence_level = 0.95
        mean4 = np.mean(fitted_n4)
        sem4 = stats.sem(fitted_n4)  # Standard error of the mean
        ci4 = stats.t.interval(confidence_level, len(fitted_n4)-1, loc=mean4, scale=sem4)

        print(f"Mean: {mean4}")
        print(f"95% Confidence Interval fitted_n4: [{ci4[0]:.8f}, {ci4[1]:.8f}]")

        cd4_cd4 = df_recovery['cd4_cd4']

        mean44 = np.mean(cd4_cd4)
        sem44 = stats.sem(cd4_cd4)  # Standard error of the mean
        ci44 = stats.t.interval(confidence_level, len(cd4_cd4)-1, loc=mean44, scale=sem44)

        print(f"Mean: {mean44}")
        print(f"95% Confidence Interval cd4_cd4: [{ci44[0]:.8f}, {ci44[1]:.8f}]")
        CIs.append({
            'true_pc':pc,
            'mean_pc':mean,
            'ci_pc_1':ci[0],
            'ci_pc_2':ci[1],
            'mean_n4':mean4,
            'ci_n4_1':ci4[0],
            'ci_n4_2':ci4[1],
            'mean_cd4_cd4':mean44,
            'ci_44_1':ci44[0],
            'ci_44_2':ci44[1],

        })

    else:
        print("No successful parameter recoveries found!")
    
    df_CIs = pd.DataFrame(CIs)
    return df_CIs

#print_results()

# Test with known parameters
true_Pc = 0.9151860764850289 #0.8311920832252657 
#i = 1
pc = 0.91
#i = 1
true_n4 = 372 #107
true_n8 = 14 #7
total_seq = true_n4 + true_n8
ci = []
df_ci = pd.DataFrame(ci)

print(f"Simulations of Pc = {true_Pc}")
print("="*60)
model = TcellLineageModel(true_Pc, true_n4, total_seq, "a")
df_new = print_results(true_Pc, true_n4,0) 
df_ci = pd.concat([df_ci,df_new],ignore_index=True)

for i in range(12,22):
    print(f"Simulations of Pc = {pc}")
    print("="*60)
    model = TcellLineageModel(pc, true_n4, total_seq,"a")
    print(f"True parameters: Pc = {pc}, n4 = {true_n4}")
    df_ci_new_row = print_results(pc,true_n4,i)
    df_ci = pd.concat([df_ci,df_ci_new_row],ignore_index=True)
    pc = pc + 0.01
   
#model = TcellLineageModel(true_Pc, true_n4, total_seq, "a")
#df_new = print_results(true_Pc, true_n4,12) 
#df_ci = pd.DataFrame(ci)
print(df_ci)
df_ci.to_csv("250923_ci_simulations_table_new.csv", index=False)