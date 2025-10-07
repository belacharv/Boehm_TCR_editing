import sympy as sp
import pandas as pd

Pc, n4 = sp.symbols('Pc n4')
###################################################################################################################
#### 0. Input

print("Insert depth:",end=" ")
maxdepth = int(input())

print("Number of CD4 samples:",end=" ")
s4 = int(input())
print("Number of CD8 samples:",end=" ")
s8 = int(input())
print("Number of CD4_only sequences:",end=" ")
cd4_cd4 = int(input())
print("Number of CD8_only sequences:",end=" ")
cd8_cd8 = int(input())
all = s4+s8
print("Number of sequences with exactly 1 CD8:",end=" ")
cd8_1_mix = int(input())
print("Number of all seq:",end=" ")
all_seq = int(input())

###################################################################################################################
#### 1. Generating coefficients ###

## numerator ##


def coeff_generator(currentS4, currentS8, sample):
    if sample == "cd4":
        return currentS4 - 1, currentS8, currentS4
    else:
        return currentS4, currentS8 - 1, currentS8

coefficients = []
visited = set()

def dfs_coefficients(currentS4, currentS8, depth, maxDepth, coeff=1):
    #print(currentS4,currentS8,depth)
    if depth == maxDepth:
        #print(coeff)
        if coeff not in visited:
            coefficients.append(coeff)
            visited.add(coeff)
        return 

    if currentS4 > 0:
        newS4, newS8, c = coeff_generator(currentS4, currentS8, "cd4")
        dfs_coefficients(newS4, newS8, depth + 1, maxDepth, coeff * c)
    
    if currentS8 > 0:
        newS4, newS8, c = coeff_generator(currentS4, currentS8, "cd8")
        dfs_coefficients(newS4, newS8, depth + 1, maxDepth, coeff * c)
    
    return 

# Generate all coefficients up to depth 5
dfs_coefficients(s4, s8, 0, maxdepth)
print(list(coefficients))


## denominator ##
samples = s4+s8
denominator = 1
for i in range(maxdepth):
    #print(denominator)
    #print(samples-i)
    denominator = denominator * (samples-i)
#print(denominator)

### list of racional numbers ###
coeff_rac =[]
for i in range(maxdepth+1):
    coeff_rac.append(sp.Rational(coefficients[i],denominator))
print(coeff_rac)

###################################################################################################################
### 2. PASCAL'S TRIANGLE
print("*"*30)
from math import factorial
coeff_pascal = []
for j in range(maxdepth+1):
    newCoeff = factorial(maxdepth)//(factorial(j)*factorial(maxdepth-j))
    print(newCoeff,end=" ")
    coeff_pascal.append(newCoeff)
print("\n")

#################################################################################################################
#### 3. PROBABILITIES ####
Pc, n4 = sp.symbols('Pc n4')

# For CD4  biased sequences
print("******** CD4 BIASED SEQUENCES ********")
p4 = []
p4_all = 0
for i in range(maxdepth,-1,-1):
    #print(i)
    newTerm = coeff_rac[maxdepth-i] * Pc**i * (1-Pc)**(maxdepth-i) * coeff_pascal[maxdepth-i]
    #print(newTerm)
    p4.append(newTerm)
    p4_all = p4_all + newTerm
print(p4_all)

# For CD8 biased sequences
print("******** CD8 BIASED SEQUENCES ********")
p8 = []
p8_all = 0
for i in range(maxdepth,-1,-1):
    #print(i)
    newTerm = coeff_rac[maxdepth-i] * Pc**(maxdepth-i) * (1-Pc)**i * coeff_pascal[maxdepth-i]
    #print(newTerm)
    p8.append(newTerm)
    p8_all = p8_all + newTerm
print(p8_all)




####################################################################################################################
#### 4.  EQUATIONS  ####

# Equations
eq1 = sp.Eq(n4 * (p4[0] / p4_all) + (all_seq - n4) * (p8[0] / p8_all), cd4_cd4)
eq2 = sp.Eq(n4 * (p4[maxdepth] / p4_all) + (all_seq - n4) * (p8[maxdepth] / p8_all), cd8_cd8)
eq3 = sp.Eq(n4 * (p4[1]/p4_all) + (all_seq - n4) * (p8[1] / p8_all),cd8_1_mix)
# print("=====================================")
# print(eq1)
# print(eq2)
# print("=====================================")

# Substitute from one equation to another
a_expr = sp.solve(eq1, n4)[0] # solve eq1
# Substitute into eq2
eq2_sub = sp.simplify(eq2.subs(n4, a_expr))
eq3_sub = sp.simplify(eq3.subs(n4,a_expr))

solutions1 = []
df_sol = pd.DataFrame(solutions1)
labels = ["cd4_only, cd8_only","cd4_only, cd8_1_mix"]
i = 0
for eq in (eq2_sub,eq3_sub):
    solutions = []
    try:
    # solve returns a list of solutions
        xv_solutions = sp.solve(eq, Pc)
    
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
            sol = sp.nsolve(eq, Pc,guess)
            av = a_expr.subs(Pc, sol)
            sol_Pc_float,sol_n4_float = float(sol),float(av)
            if 0 <= sol_Pc_float <= 1:
                print(f"  Pc = {sol_Pc_float:.10f}")
                if sol in solutions:
                    pass
                else:
                    sol3 = [{
                        'n4': sol_n4_float,
                        'value': sol_Pc_float,
                        'depth':maxdepth,
                        'equation':labels[i]}]
                    df_new = pd.DataFrame(sol3)
                    df_sol = pd.concat([df_sol,df_new],ignore_index=True)
        except:
            pass
    print(df_sol)
    i = i +1

print("\n" + "*"*50)
print(df_sol)
print("*"*50)

# TESTING
print("\n" + "="*50)
print("TESTING SOLUTIONS")
print("="*50)
len(solutions)
for i, (Pc_val, n4_val) in enumerate(solutions):
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
    # p4_44_val = p4[0].subs(Pc, Pc_val)
    # p4_48_val = p4[1].subs(Pc, Pc_val)
    # p4_88_val = p4[2].subs(Pc, Pc_val)
    
    # p8_44_val = p8[0].subs(Pc, Pc_val)
    # p8_48_val = p8[1].subs(Pc, Pc_val)
    # p8_88_val = p8[2].subs(Pc, Pc_val)
    
    # Calculate CD4-CD4 count (should equal cd4_cd4)
    cd4_cd4_count = n4_val * (p4_val[0] / p4_all_val) + (all_seq - n4_val) * (p8_val[0] / p8_all_val)
    
    # Calculate CD8-CD8 count (should equal cd8_cd8)
    cd8_cd8_count = n4_val * (p4_val[maxdepth] / p4_all_val) + (all_seq - n4_val) * (p8_val[maxdepth] / p8_all_val)
    
    print(f"  CD4-CD4 count: {float(cd4_cd4_count)} (expected: {cd4_cd4})")
    print(f"  CD8-CD8 count: {float(cd8_cd8_count)} (expected: {cd8_cd8})")
    print(f"  Error in CD4-CD4: {abs(float(cd4_cd4_count) - cd4_cd4)}")
    print(f"  Error in CD8-CD8: {abs(float(cd8_cd8_count) - cd8_cd8)}")


# currS4,currS8 = coeff_generator(currS4,currS8,"cd4")
# currS4,currS8 = coeff_generator(currS4,currS8,"cd8")
# currS4,currS8 = coeff_generator(currS4,currS8,"cd4")
# currS4,currS8 = coeff_generator(currS4,currS8,"cd4")
#print(generating_constants(3))