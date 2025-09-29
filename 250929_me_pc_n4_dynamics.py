# SUBSTITUTING N4 AND LOOKING AT PC
import sympy as sp

#n4 = 337
# Observed counts
cd4_cd4_obs = 337
cd4_cd8_obs = 39
cd8_cd8_obs = 10
all_seq = 386



def solvingEq(n4,chain):
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

    solutions = []
    # Method 1: Using CD4-CD4 eq
    #print("\n" + "-"*60)
    #print("Method 1: Using CD4-CD4 equation")
    expected_cd4_cd4 = n4 * (p4_44/p4_total) + (all_seq-n4) * (p8_44/p8_total)
    eq1 = sp.Eq(expected_cd4_cd4, cd4_cd4_obs)
    #print(f"Equation: {eq1}")
    solutions_1 = sp.solve(eq1, Pc)
    #print(f"\nSolutions from CD4-CD4 constraint:")
    for sol in solutions_1:
        try:
            sol_float = float(sol)
            if 0 <= sol_float <= 1:
                print(f"  Pc = {sol_float:.10f}")
                solutions.append(sol_float)
        except:
            print(float(sol))
            pass

    # Method 2: Using CD8-CD8 eq
    #print("\n" + "-"*60)
    #print("Method 2: Using CD8-CD8 equation")
    expected_cd8_cd8 = n4 * (p4_88/p4_total) + (all_seq-n4) * (p8_88/p8_total)
    eq2 = sp.Eq(expected_cd8_cd8, cd8_cd8_obs)
    #print(f"Equation: {eq2}")

    solutions_2 = sp.solve(eq2, Pc)
    #print(f"\nSolutions from CD8-CD8 constraint:")
    for sol in solutions_2:
        try:
            sol_float = float(sol)
            if 0 <= sol_float <= 1:
                print(f"  Pc = {sol_float:.10f}")
                solutions.append(sol_float)
        except:
            print("heree")
            pass

    # Method 3: Using CD4-CD8 count
    #print("\n" + "-"*60)
    #print("Method 3: Using CD4-CD8 equation")
    expected_cd4_cd8 = n4 * (p4_48/p4_total) + (all_seq-n4) * (p8_48/p8_total)
    eq3 = sp.Eq(expected_cd4_cd8, cd4_cd8_obs)
    #print(f"Equation: {eq3}")

    solutions_3 = sp.solve(eq3, Pc)
    #print(f"\nSolutions from CD4-CD8 constraint:")
    for sol in solutions_3:
        try:
            sol_float = float(sol)
            if 0 <= sol_float <= 1:
                print(f"  Pc = {sol_float:.10f}")
                solutions.append(sol_float)
        except:
            print("heree")
            pass
    print(n4, solutions)
    return solutions

def verify_solution(Pc_val, n4_val, n8_val, c44, c48, c88):
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
    print(f"  Expected CD4-CD4: {exp_44:.4f} (observed: {cd4_cd4_obs})")
    print(f"  Expected CD4-CD8: {exp_48:.4f} (observed: {cd4_cd8_obs})")
    print(f"  Expected CD8-CD8: {exp_88:.4f} (observed: {cd8_cd8_obs})")
    print(f"  Total: {exp_44+exp_48+exp_88:.4f}")
    
    return exp_44, exp_48, exp_88
"""
for n4 in range(337,386):
    n8 = all_seq - n4
    solutions_1 = solvingEq(n4,"a")
    # Verify valid solutions from each method
    if solutions_1:
        try:
            for sol in solutions_1:
                sol_float = float(sol)
                if 0 <= sol_float <= 1:
                    verify_solution(sol_float, n4, all_seq-n4)
        except:
            pass
"""
soloo = solvingEq(386,"a")
