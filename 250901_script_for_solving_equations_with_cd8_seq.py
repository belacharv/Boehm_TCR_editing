# EQUATIONS WHERE WE COUNT WITH TRUE CD8 SEQ
import sympy as sp

Pc, n4 = sp.symbols('Pc n4')
# constants for each combination
    # c44 = constant CD4-CD4
c44 = sp.Rational(10*9, 15*14)   # 3/7
c48 = sp.Rational(10*5, 15*14)   # 5/21
c88 = sp.Rational(5*4, 15*14)    # 2/21
# Probability to go by one of the 3 possible cases 
    # cd4-cd4 marked as 44
    # cd4-cd8 marked as 48
    # cd8-cd8 marked as 88
# For CD4 biased sequences
p4_44 = Pc**2 * c44 # both of the T cells chose the right fate (according to the bias of the sequence)
p4_48 = 2*Pc*(1-Pc) * c48 # one was correct, one was wrong
p4_88 = (1-Pc)**2 * c88 # both were wrong
# For CD8 biased sequences
p8_44 = (1-Pc)**2 * c44
p8_48 = 2*Pc*(1-Pc) * c48
p8_88 = Pc**2 * c88

# Equations
eq1 = sp.Eq(n4 * (p4_44 / (p4_44 + p4_48 + p4_88)) + (386 - n4) * (p8_44 / (p8_44 + p8_48 + p8_88)), 337)
eq2 = sp.Eq(n4 * (p4_88 / (p4_44 + p4_48 + p4_88)) + (386 - n4) * (p8_88 / (p8_44 + p8_48 + p8_88)), 10)

print("=====================================")
print(eq1)
print(eq2)
print("=====================================")


# Substitute from one equation to another
a_expr = sp.solve(eq1, n4)[0] # solve eq1
# Substitute into eq2
eq2_sub = sp.simplify(eq2.subs(n4, a_expr))

# Use nsolve with different starting guesses for x
# solutions = []
# for guess in [0, 1]: # nsolve needs an interval in which it will look for solutions
#     try:
#         xv = sp.nsolve(eq2_sub, Pc, guess)
#         av = a_expr.subs(Pc, xv)
#         solutions.append((float(xv), float(av)))
#     except:
#         pass

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

print("Solutions (Pc, n4):")
for s in solutions:
    print(s)


# TESTING
print("\n" + "="*50)
print("TESTING SOLUTIONS")
print("="*50)
len(solutions)
for i, (Pc_val, n4_val) in enumerate(solutions):
    print(f"\nTesting Solution {i+1}: Pc = {Pc_val:.6f}, n4 = {n4_val:.6f}")
    
    # Substitute the solution values into the probability expressions
    p4_44_val = p4_44.subs(Pc, Pc_val)
    p4_48_val = p4_48.subs(Pc, Pc_val)
    p4_88_val = p4_88.subs(Pc, Pc_val)
    
    p8_44_val = p8_44.subs(Pc, Pc_val)
    p8_48_val = p8_48.subs(Pc, Pc_val)
    p8_88_val = p8_88.subs(Pc, Pc_val)
    
    # Calculate CD4-CD4 count (should equal 337)
    cd4_cd4_count = n4_val * (p4_44_val / (p4_44_val + p4_48_val + p4_88_val)) + (386 - n4_val) * (p8_44_val / (p8_44_val + p8_48_val + p8_88_val))
    
    # Calculate CD8-CD8 count (should equal 10)
    cd8_cd8_count = n4_val * (p4_88_val / (p4_44_val + p4_48_val + p4_88_val)) + (386 - n4_val) * (p8_88_val / (p8_44_val + p8_48_val + p8_88_val))
    
    print(f"  CD4-CD4 count: {float(cd4_cd4_count)} (expected: 337)")
    print(f"  CD8-CD8 count: {float(cd8_cd8_count)} (expected: 10)")
    print(f"  Error in CD4-CD4: {abs(float(cd4_cd4_count) - 337)}")
    print(f"  Error in CD8-CD8: {abs(float(cd8_cd8_count) - 10)}")