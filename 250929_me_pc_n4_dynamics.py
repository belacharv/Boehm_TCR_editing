# SUBSTITUTING N4 AND LOOKING AT PC
import sympy as sp
import pandas as pd


def solvingEq(n4,chain,obs_counts,all_seq):
    Pc = sp.symbols('Pc')
    # constants for each combination
    # c44 = constant CD4-CD4
    if chain == "a":
        c44 = sp.Rational(10*9, 15*14)   # 3/7
        c48 = sp.Rational(10*5, 15*14)   # 5/21
        c88 = sp.Rational(5*4, 15*14)
    elif chain == "b":
        c44 = sp.Rational(21*20, 31*30)
        c48 = sp.Rational(21*10, 31*30)
        c88 = sp.Rational(10*9, 31*30)
    else:
        raise Exception("Wrong input, uncorrect TCR chain selection")

    # For CD4 biased sequences  
    p4_44 = Pc**2 * c44 # both of the T cells chose the right fate (according to the bias of the sequence)
    p4_48 = 2*Pc*(1-Pc) * c48 # one was correct, one was wrong
    p4_88 = (1-Pc)**2 * c88 # both were wrong
    p4_total = p4_44 + p4_48 + p4_88

    # For CD8 biased sequences
    p8_44 = (1-Pc)**2 * c44
    p8_48 = 2*Pc*(1-Pc) * c48
    p8_88 = Pc**2 * c88
    p8_total = p8_44 + p8_48 + p8_88

    # Method 1: Using CD4-CD4 eq
    eq1 = sp.Eq(n4 * (p4_44 / p4_total) + (all_seq - n4) * (p8_44 / p8_total), obs_counts[0])
    print(eq1)
    solutions = []
    df_sol = pd.DataFrame(solutions)
    #print(f"Equation: {eq1}")
    for guess in [0.8,0.95]:
        sol = sp.nsolve(eq1, Pc,guess)
        sol_float = float(sol)
        if 0 <= sol_float <= 1:
            print(f"  Pc = {sol_float:.10f}")
            if sol in df_sol:
                pass
            else:
                sol1 = [{
                    'n4': n4,
                    'value': sol_float,
                    'method':"cd4-cd4"}]
                df_new = pd.DataFrame(sol1)
                df_sol = pd.concat([df_sol,df_new],ignore_index=True)
    #print(f"\nSolutions from CD4-CD4 constraint:")
    
    # Method 2: Using CD8-CD8 eq
    expected_cd8_cd8 = n4 * (p4_88/p4_total) + (all_seq-n4) * (p8_88/p8_total)
    eq2 = sp.Eq(expected_cd8_cd8, obs_counts[2])
    for guess in [0.8,0.95]:
        try:
            sol = sp.nsolve(eq2, Pc,guess)
            sol_float = float(sol)
            if 0 <= sol_float <= 1:
                print(f"  Pc = {sol_float:.10f}")
                if sol in solutions:
                    pass
                else:
                    sol2 = [{
                        'n4': n4,
                        'value': sol_float,
                        'method': "cd8-cd8"}]
                    df_new = pd.DataFrame(sol2)
                    df_sol = pd.concat([df_sol,df_new],ignore_index=True)
        except:
            pass
            
    
    # Method 3: Using CD4-CD8 count
    #print("\n" + "-"*60)
    #print("Method 3: Using CD4-CD8 equation")
    expected_cd4_cd8 = n4 * (p4_48/p4_total) + (all_seq-n4) * (p8_48/p8_total)
    eq3 = sp.Eq(expected_cd4_cd8, obs_counts[1])
    #print(f"Equation: {eq3}")
    
    #solutions_3 = sp.solve(eq3, Pc)
    #print(f"\nSolutions from CD4-CD8 constraint:")
    for guess in [0.8,0.95]:
        sol = sp.nsolve(eq3, Pc,guess)
        sol_float = float(sol)
        if 0 <= sol_float <= 1:
            print(f"  Pc = {sol_float:.10f}")
            if sol in solutions:
                pass
            else:
                sol3 = [{
                    'n4': n4,
                    'value': sol_float,
                    'method':"cd4-cd8"}]
                df_new = pd.DataFrame(sol3)
                df_sol = pd.concat([df_sol,df_new],ignore_index=True)
    
    #print(solutions)
    #df_sol = pd.DataFrame(solutions)
    print(df_sol)
    return df_sol

def verify_solution(Pc_val, n4_val, n8_val, c44, c48, c88,obs_counts):
    # Verify the solution
    print("\n" + "="*60)
    print("VERIFICATION")
    print("="*60)
    """Verify a solution by calculating expected counts"""
    # For CD4 biased sequences
    p4_44_val = Pc_val**2 * float(c44)
    p4_48_val = 2*Pc_val*(1-Pc_val) * float(c48)
    p4_88_val = (1-Pc_val)**2 * float(c88)
    p4_total_val = p4_44_val + p4_48_val + p4_88_val
    
    # For CD8 biased sequences
    p8_44_val = (1-Pc_val)**2 * float(c44)
    p8_48_val = 2*Pc_val*(1-Pc_val) * float(c48)
    p8_88_val = Pc_val**2 * float(c88)
    p8_total_val = p8_44_val + p8_48_val + p8_88_val
    
    exp_44 = n4_val * (p4_44_val/p4_total_val) + n8_val * (p8_44_val/p8_total_val)
    exp_48 = n4_val * (p4_48_val/p4_total_val) + n8_val * (p8_48_val/p8_total_val)
    exp_88 = n4_val * (p4_88_val/p4_total_val) + n8_val * (p8_88_val/p8_total_val)
    
    print(f"\nFor Pc={Pc_val:.10f}, n4={n4_val}, n8={n8_val}:")
    print(f"  Expected CD4-CD4: {exp_44:.4f} (observed: {obs_counts[0]})")
    print(f"  Expected CD4-CD8: {exp_48:.4f} (observed: {obs_counts[1]})")
    print(f"  Expected CD8-CD8: {exp_88:.4f} (observed: {obs_counts[2]})")
    print(f"  Total: {exp_44+exp_48+exp_88:.4f}")
    
    return exp_44, exp_48, exp_88

def generate_table(all_seq,chain,counts):
    sol = []
    df_sol_all = pd.DataFrame(sol)
    for n4 in range(counts[0],all_seq+1):
        #n8 = all_seq - n4
        df_sol_new = solvingEq(n4,chain,counts,all_seq)
    # Verify valid solutions from each method
    # try:
    #     for sol in df_sol_new:
    #             sol_float = float(sol)
    #             if 0 <= sol_float <= 1:
    #                 verify_solution(sol_float, n4, all_seq-n4)
    #     except:
    #         pass
        df_sol_all = pd.concat([df_sol_all,df_sol_new],ignore_index=True)

    print("\n" + "="*60)
    print(df_sol_all)
    filename = "250930_"+ chain +"_pc_estim_n4_given.csv"
    df_sol_all.to_csv(filename, index=False)

    return df_sol_all
#soloo = solvingEq(372,"a")
#df_ci = pd.DataFrame(ci)

import matplotlib.pyplot as plt
def draw_plots(df):
# Assuming your dataframe is called 'df'
    plt.figure(figsize=(10, 6))

    ## all 
    for method in df['method'].unique():
        mask = df['method'] == method
        method_data = df[mask].sort_values('value')  # Sort by x-axis for proper line plotting
        plt.plot(method_data['value'], method_data['n4'], marker='o', label=method, linewidth=2)

    #plt.scatter(df_sol_all['value'], df_sol_all['n4'])
    plt.xlabel('Pc')
    plt.ylabel('true_cd4')
    plt.legend()
    plt.title('Estimating Pc while true_cd4 fixed')
    plt.grid(True, alpha=0.3)
    plt.show()

    colors = ['blue', 'red', 'green']
    i=0
    ## each method separately
    for method in df['method'].unique():
        df_filtered = df[df['method'] == method]
        plt.scatter(df_filtered['value'], df_filtered['n4'], color=colors[i])
        plt.xlabel('Pc')
        plt.ylabel('true_cd4')
        plt.title(f"Estimating Pc while true_cd4 fixed from {method} eq")
        plt.grid(True, alpha=0.3)
        plt.show()
        i += 1

##### ALPHA ############################################################################################
#n4 = 337
# Observed counts
obs_alpha = [337,39,10]
all_seq_a = 386

df_alpha = generate_table(all_seq_a,"a",obs_alpha)
print(df_alpha)
draw_plots(df_alpha)


##### BETA #############################################################################################
obs_beta = [89,21,4]
all_seq_b = 114

df_beta = generate_table(all_seq_b,"b",obs_beta)
print(df_beta)
draw_plots(df_beta)


#plt.show()