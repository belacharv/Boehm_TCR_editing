import sympy as sp
import pandas as pd
#import pandas as pd

###################################################################################################################
#### 0. Input
# input_row = [maxdepth, s4, s8, all_seq] + cd8_values

# input function that for each depth return information about counts specific for that depth(=layer)
    # input_row = [maxdepth, s4, s8, all_seq] + cd8_values
def input_function(depth,chain):  
    
    # read the counts of each group in each layer
    if chain == "alpha":
        s4 = 10 # cd4 samples
        s8 = 5 # cd8 samples
        df = pd.read_csv('251016_a_seq_counts_distribution_groups_layers.csv')
    elif chain == "beta":
        s4 = 21
        s8 = 10
        df = pd.read_csv('251016_b_seq_counts_distribution_groups_layers.csv')
    else:
        raise Exception("wrong input")
    # we start from 2nd layer (seq present in 2 samples), we need to adjust for beginning of the list
    i = depth - 2
    maxdepth = df.iloc[i, 0]
    all_seq = df.iloc[i, 1]
    cd8_values = df.iloc[i, 2:].tolist()

    # vital row storing all the info specific for each layer
        # maxdepth refers to current layer
        # all_seq is Sequence count in that particular layer
        # cd8_values are group divided on the amount of cd8s in that group
            # cd8_0 (cd4_only), cd8_1, cd8_2 ... cd8_max (all the samples of cd8 present)
                # cd8_max means that in lower layers than s8 it's cd8_only
                # in higher layers it's all the samples of cd8 present
    input_row = [maxdepth, s4, s8, all_seq] + cd8_values
    print(input_row)
    return input_row


###################################################################################################################
#### 1. Generating coefficients ###
# each probability is influenced by the s4 and s8 constants

## numerator ##
def coeff_generator(currentS4, currentS8, sample):
    if sample == "cd4":
        return currentS4 - 1, currentS8, currentS4
    else:
        return currentS4, currentS8 - 1, currentS8
    
# needs empty coefficient [] and visited () before starts
def dfs_coefficients(currentS4, currentS8, depth, maxDepth,coefficients,visited,coeff=1):
    if depth == maxDepth:
        if coeff not in visited:
            coefficients.append(coeff)
            visited.add(coeff)
        return 

    if currentS4 > 0:
        newS4, newS8, c = coeff_generator(currentS4, currentS8, "cd4")
        dfs_coefficients(newS4, newS8, depth + 1, maxDepth,coefficients,visited,coeff * c)
    
    if currentS8 > 0:
        newS4, newS8, c = coeff_generator(currentS4, currentS8, "cd8")
        dfs_coefficients(newS4, newS8, depth + 1, maxDepth,coefficients,visited, coeff * c)
    
    return 

# Generate all coefficients up to depth 5

def coeff_rac_produce(s4,s8,depth,coefficients):
    ## denominator ##
    samples = s4+s8
    denominator = 1
    for i in range(depth):
    #print(denominator)
    #print(samples-i)
        denominator = denominator * (samples-i)
    #print(denominator)

    ### list of racional numbers ###
    coeff_rac =[]
    for i in range(depth+1):
        coeff_rac.append(sp.Rational(coefficients[i],denominator))
    print(coeff_rac)
    return coeff_rac

###################################################################################################################
### 2. PASCAL'S TRIANGLE ###
## We use Pascal's triangle values in this code for the binomial theorem
# then, we use these values in equations later
# i.e. 
       # 1
      # 1 1
     # 1 2 1    # for 2nd layer
    # 1 3 3 1   # for 3rd layer

print("*"*30)
from math import factorial
def pascal(depth):
    coeff_pascal = []
    for j in range(depth+1):
        newCoeff = factorial(depth)//(factorial(j)*factorial(depth-j))
        print(newCoeff,end=" ")
        coeff_pascal.append(newCoeff)
    print("\n")
    return coeff_pascal

#################################################################################################################
#### 3. PROBABILITIES ####
#Pc, n4 = sp.symbols('Pc n4')


# For CD4  biased sequences
# print("******** CD4 BIASED SEQUENCES ********")
# coeff_pascal = pascal(maxdepth)
# coeff_rac = coeff_rac_produce(s4,s8,maxdepth)
# generating different probabilities for each group
def cd4_bias(c_pascal,c_rac,depth):
    Pc, n4 = sp.symbols('Pc n4')
    p4 = []
    p4_all = 0
    for i in range(depth,-1,-1):
    #print(i)
        newTerm = c_rac[depth-i] * Pc**i * (1-Pc)**(depth-i) * c_pascal[depth-i]
    #print(newTerm)
        p4.append(newTerm)
        p4_all = p4_all + newTerm
    print(p4_all)
    return p4, p4_all

# For CD8 biased sequences


def cd8_bias(c_pascal,c_rac,depth):
    Pc, n4 = sp.symbols('Pc n4')
    p8 = []
    p8_all = 0
    for i in range(depth,-1,-1):
    #print(i)
        newTerm = c_rac[depth-i] * Pc**(depth-i) * (1-Pc)**i * c_pascal[depth-i]
    #print(newTerm)
        p8.append(newTerm)
        p8_all = p8_all + newTerm
    print(p8_all)
    return p8,p8_all




####################################################################################################################
#### 4.  EQUATIONS  ####


def eqSolving(eq1,eq2,depth,label):
    Pc, n4 = sp.symbols('Pc n4')
    print("first solved")
    # Substitute from one equation to another
    a_expr = sp.solve(eq1, n4)[0] # solve eq1
    print("second solved")
    # Substitute into eq2
    eq2_sub = sp.simplify(eq2.subs(n4, a_expr))
    print("third solved")
    #eq3_sub = sp.simplify(eq3.subs(n4,a_expr))

    solutions1 = []
    df_sol = pd.DataFrame(solutions1)
    #labels = ["cd4_only, cd8_only","cd4_only, cd8_1_mix"]
    i = 0
    solutions = []
    try:
    # solve returns a list of solutions
        xv_solutions = sp.solve(eq2_sub, Pc)
        for xv in xv_solutions:
        # Filter for real, finite solutions
            if xv.is_real and xv.is_finite:
                av = a_expr.subs(Pc, xv)
                if av.is_real and av.is_finite:
                    solutions.append((float(xv), float(av)))
    except:
        pass

    print("Solutions 1 (Pc, n4):")
    for s in solutions:
        print(s)

    for guess in [0.01, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
        try:
            sol = sp.nsolve(eq2_sub, Pc,guess)
            av = a_expr.subs(Pc, sol)
            sol_Pc_float,sol_n4_float = float(sol),float(av)
            if 0 <= sol_Pc_float <= 1:
                print(f"  Pc = {sol_Pc_float:.10f}")
                if sol in solutions:
                    pass
                else:
                    sol3 = [{
                        'n4': sol_n4_float,
                        'Pc_val': sol_Pc_float,
                        'depth':depth,
                        'equation':label,
                        'expected_cd4':-1,
                        'expected_cd8':-1}]
                    df_new = pd.DataFrame(sol3)
                    df_sol = pd.concat([df_sol,df_new],ignore_index=True)
        except:
            pass
    print(df_sol)
    i = i +1

    print("\n" + "*"*50)
    print(solutions)
    print("*"*50)
    print(df_sol)
    print("*"*50)
    return solutions,df_sol

def main():
    df_final = pd.DataFrame([])
    print("Write the name of TCR chain: ")
    chain = input()
    for i in range(2,6):

        coefficients = []
        visited = set()
        # Generate all coefficients up to depth 5
        input_r = input_function(i,chain)
        #input_row = [maxdepth, s4, s8, all_seq] + cd8_values
        maxdepth = input_r[0]
        s4,s8 = input_r[1],input_r[2]
        all_seq = input_r[3]

        obs_val = input_r[4:]
        print("hiiiii")
        print(obs_val)
        for j in range(len(obs_val)):
            if obs_val[j] == "-":
                obs_val[j] = None
            else:
                obs_val[j] = int(obs_val[j])
        print(obs_val)
        #cd4_cd4,cd8_cd8,cd8_1_mix = int(input_r[4]),input_r[6],input_r[5]
        #all_seq = input_r[3]

        dfs_coefficients(s4, s8, 0, maxdepth,coefficients,visited)
        print(list(coefficients))
        Pc, n4 = sp.symbols('Pc n4')
        coeff_pascal = pascal(maxdepth)
        coeff_rac = coeff_rac_produce(s4,s8,maxdepth,coefficients)
        p4,p4_all = cd4_bias(coeff_pascal,coeff_rac,maxdepth)
        p8,p8_all = cd8_bias(coeff_pascal,coeff_rac,maxdepth)
        # Equations
        obs_val_index = 0
        print(len(obs_val))
        equations = []
        while obs_val_index < len(obs_val) and obs_val[obs_val_index] != None:
            print(obs_val_index)
            print(obs_val[obs_val_index])
            right_side_eq = obs_val[obs_val_index]
            eq = sp.Eq(n4 * (p4[obs_val_index] / p4_all) + (all_seq - n4) * (p8[obs_val_index] / p8_all), right_side_eq)
            equations.append(eq)
            obs_val_index += 1
        #eq1 = sp.Eq(n4 * (p4[0] / p4_all) + (all_seq - n4) * (p8[0] / p8_all), cd4_cd4)
        #eq2 = sp.Eq(n4 * (p4[maxdepth] / p4_all) + (all_seq - n4) * (p8[maxdepth] / p8_all), cd8_cd8)
        #eq3 = sp.Eq(n4 * (p4[1]/p4_all) + (all_seq - n4) * (p8[1] / p8_all),cd8_1_mix)
        #df_sol = []
        for cd8_val_eq in range(1,i+1):
            print(cd8_val_eq)
            method_name = "CD8_0, CD8_" + str(cd8_val_eq)
            solutions,df_sol = eqSolving(equations[0],equations[cd8_val_eq],maxdepth,method_name)
        #solu, df_sol1 = eqSolving(equations[0],equations[2],maxdepth,"CD8_0, CD8_1")
        # Apply testing function and get CD4/CD8 counts for each solution
            if df_sol.empty != True:
                cd4_counts, cd8_counts = testing(df_sol,coefficients,input_r)
        
        # Add the counts as new columns to df_sol
                df_sol['cd4_cd4_calculated'] = cd4_counts
                df_sol['cd8_cd8_calculated'] = cd8_counts
                df_final = pd.concat([df_final,df_sol],ignore_index=True)

        #df_final = pd.concat([df_final,df_sol],ignore_index=True)
        #testing(solutions,coefficients,input_r)

    print("\n" + "*"*60)
    print(df_final)
    filename = "2510200_"+chain[0]+"_depths_pc_n4_calculated.csv"
    df_final.to_csv(filename)

def testing(df_sol,coeffs,input_r):
# TESTING
    # input_row = [maxdepth, s4, s8, all_seq] + cd8_values
    maxdepth = input_r[0]
    s4,s8 = input_r[1],input_r[2]
    #cd4_cd4,cd8_cd8,cd8_1_mix = input_r[3],input_r[4],input_r[5]
    all_seq = input_r[3]

    obs_val = input_r[4:]
    for i in range(len(obs_val)):
        if obs_val[i] == "-" or obs_val[i] == "" or pd.isna(obs_val[i]):
            obs_val[i] = None
        else:
            obs_val[i] = float(obs_val[i])
    
    cd4_cd4 = obs_val[0] if obs_val[0] is not None else 0  # CD8_0 (CD4_only)
    cd8_cd8 = obs_val[maxdepth] if len(obs_val) > maxdepth and obs_val[maxdepth] is not None else 0  # CD8_maxdepth



    Pc, n4 = sp.symbols('Pc n4')
    print("******** CD4 BIASED SEQUENCES ********")
    coeff_pascal = pascal(maxdepth)
    coeff_rac = coeff_rac_produce(s4,s8,maxdepth,coeffs)
    p4,p4_all = cd4_bias(coeff_pascal,coeff_rac,maxdepth)
    print("******** CD8 BIASED SEQUENCES ********")
    p8,p8_all = cd8_bias(coeff_pascal,coeff_rac,maxdepth)
    print("\n" + "="*50)
    print("TESTING SOLUTIONS")
    print("="*50)
#solutions,df_sol = eqSolving(eq1,eq2,eq3)
    cd4_cd4_counts = []
    cd8_cd8_counts = []
    n4_val_all,Pc_val_all = df_sol['n4'], df_sol['Pc_val']
    df_sol['expected_cd4'], df_sol['expected_cd8'] = cd4_cd4, cd8_cd8 
    #i = 0
    for i in range(len(df_sol['n4'])):
        n4_val,Pc_val = n4_val_all[i],Pc_val_all[i]
        print(f"\nTesting Solution {i+1}: Pc = {Pc_val:.6f}, n4 = {n4_val:.6f}")
    
    # Substitute the solution values into the probability expressions
        p4_val = []
        p4_all_val = 0
        p8_val = []
        p8_all_val = 0
        for j in range(maxdepth+1):
            p4_val.append(p4[j].subs(Pc, Pc_val))
            p4_all_val = p4_all_val + p4_val[j]

        for k in range(maxdepth+1):
            p8_val.append(p8[k].subs(Pc, Pc_val))
            p8_all_val = p8_all_val + p8_val[k]
    
        cd4_cd4_count = n4_val * (p4_val[0] / p4_all_val) + (all_seq - n4_val) * (p8_val[0] / p8_all_val)
        cd8_cd8_count = n4_val * (p4_val[maxdepth] / p4_all_val) + (all_seq - n4_val) * (p8_val[maxdepth] / p8_all_val)
        cd4_cd4_counts.append(float(cd4_cd4_count))
        cd8_cd8_counts.append(float(cd8_cd8_count))

        print(f"  CD4-CD4 count: {float(cd4_cd4_count)} (expected: {cd4_cd4})")
        print(f"  CD8-CD8 count: {float(cd8_cd8_count)} (expected: {cd8_cd8})")
        print(f"  Error in CD4-CD4: {abs(float(cd4_cd4_count) - cd4_cd4)}")
        print(f"  Error in CD8-CD8: {abs(float(cd8_cd8_count) - cd8_cd8)}")
        #i = i+1

    return cd4_cd4_counts, cd8_cd8_counts
main()
input_function(2,"alpha")
input_function(2,"beta")
input_function(15,"alpha")