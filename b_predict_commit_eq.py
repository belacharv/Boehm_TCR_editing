import numpy as np
#import matplotlib.pyplot as plt
import sympy as sp
import pandas as pd
from scipy import stats

class TcellLineageModel:
    def __init__(self, Pc, n4_total, total_sequences, chain):
        self.Pc = Pc  # Commitment precision
        self.n4_total = n4_total  # Number of CD4-biased sequences
        self.n8_total = total_sequences - n4_total  # Number of CD8-biased sequences
        self.total_sequences = total_sequences
        #print(Pc)
        #print(total_sequences)
        #print(self.n8_total)
        
        # Constants from your original model
        # for alpha
        if chain == "a":
            self.c44 = 10*9/(15*14)  # 3/7
            self.c48 = 10*5/(15*14)  # 5/21
            self.c88 = 5*4/(15*14)   # 2/21

        # for beta
        if chain == "b":
            self.c44 = 21*20/(31*30)
            self.c48 = 21*10/(31*30)
            self.c88 = 10*9/(31*30)
    
    def simulate_dataset(self, seed=None):
        """
        Each sequence appears in exactly 2 T-cells, so we simulate 386 sequence pairs.
        Each sequence is either CD4-biased or CD8-biased, determining the probabilities
        for the pair outcomes (CD4-CD4, CD4-CD8, CD8-CD8).
        """
        if seed is not None:
            np.random.seed(seed)
        
        cd4_cd4_count = 0
        cd4_cd8_count = 0
        cd8_cd8_count = 0
        
        # Simulate all 386 sequences (each appearing in exactly 2 T-cells)
        for seq_idx in range(self.total_sequences):
            
            # Determine if this sequence is CD4-biased or CD8-biased
            if seq_idx < self.n4_total:
                # This is a CD4-biased sequence
                p_44 = self.Pc**2 * self.c44
                p_48 = 2*self.Pc*(1-self.Pc) * self.c48
                p_88 = (1-self.Pc)**2 * self.c88
            else:
                # This is a CD8-biased sequence
                p_44 = (1-self.Pc)**2 * self.c44
                p_48 = 2*self.Pc*(1-self.Pc) * self.c48
                p_88 = self.Pc**2 * self.c88
            
            # Normalize probabilities
            total_p = p_44 + p_48 + p_88
            probs = [p_44/total_p, p_48/total_p, p_88/total_p]
            
            # Sample outcome for this sequence pair
            outcome = np.random.choice([0, 1, 2], p=probs)  # 0=CD4-CD4, 1=CD4-CD8, 2=CD8-CD8
            
            if outcome == 0:
                cd4_cd4_count += 1
            elif outcome == 1:
                cd4_cd8_count += 1
            else:
                cd8_cd8_count += 1
        
        # Verify total count
        total_count = cd4_cd4_count + cd4_cd8_count + cd8_cd8_count
        assert total_count == self.total_sequences, f"Total count {total_count} != {self.total_sequences}"
        
        return cd4_cd4_count, cd4_cd8_count, cd8_cd8_count
    
        
    # don't think should be used, it's basically substitution to the equation
    def expected_counts(self):
        """Calculate expected counts based on model"""
        # For CD4-biased sequences
        p4_44 = self.Pc**2 * self.c44
        p4_48 = 2*self.Pc*(1-self.Pc) * self.c48
        p4_88 = (1-self.Pc)**2 * self.c88
        p4_total = p4_44 + p4_48 + p4_88
        
        # For CD8-biased sequences
        p8_44 = (1-self.Pc)**2 * self.c44
        p8_48 = 2*self.Pc*(1-self.Pc) * self.c48
        p8_88 = self.Pc**2 * self.c88
        p8_total = p8_44 + p8_48 + p8_88
        
        expected_cd4_cd4 = self.n4_total * (p4_44/p4_total) + self.n8_total * (p8_44/p8_total)
        expected_cd4_cd8 = self.n4_total * (p4_48/p4_total) + self.n8_total * (p8_48/p8_total)
        expected_cd8_cd8 = self.n4_total * (p4_88/p4_total) + self.n8_total * (p8_88/p8_total)
        print(f"Expected results: {expected_cd4_cd4}, {expected_cd4_cd8}, {expected_cd8_cd8}")
        return expected_cd4_cd4, expected_cd4_cd8, expected_cd8_cd8

def fit_model_to_data(cd4_cd4_obs, cd8_cd8_obs, total_sequences,chain):
    """Fit the model to observed data using sympy.solve for exact solutions"""
    
    # Define symbolic variables
    Pc, n4 = sp.symbols('Pc n4', real=True)
    
    # Constants - using exact fractions for better symbolic computation
    # c44 = sp.Rational(21*20, 31*30)
    # c48 = sp.Rational(21*10, 31*30)
    # c88 = sp.Rational(10*9, 31*30)


    if chain == "a":
        c44 = sp.Rational(10*9, 15*14) 
        c48 = sp.Rational(10*5, 15*14) 
        c88 = sp.Rational(5*4, 15*14)

        # for beta
    if chain == "b":
        c44 = sp.Rational(21*20, 31*30)
        c48 = sp.Rational(21*10, 31*30)
        c88 = sp.Rational(10*9, 31*30)
    
    # Probabilities for CD4-biased sequences
    p4_44 = Pc**2 * c44
    p4_48 = 2*Pc*(1-Pc) * c48
    p4_88 = (1-Pc)**2 * c88
    p4_total = p4_44 + p4_48 + p4_88
    
    # Probabilities for CD8-biased sequences
    p8_44 = (1-Pc)**2 * c44
    p8_48 = 2*Pc*(1-Pc) * c48
    p8_88 = Pc**2 * c88
    p8_total = p8_44 + p8_48 + p8_88
    
    # Equations to solve
    eq1 = n4 * (p4_44/p4_total) + (total_sequences - n4) * (p8_44/p8_total) - cd4_cd4_obs
    eq2 = n4 * (p4_88/p4_total) + (total_sequences - n4) * (p8_88/p8_total) - cd8_cd8_obs
    
    # Solve the system of equations
    solutions = sp.solve([eq1, eq2], [Pc, n4])
    
    # Filter valid solutions
    valid_solutions = []
    
    # Handle different types of solutions
    if isinstance(solutions, list):
        # Multiple solutions
        for sol in solutions:
            try:
                if isinstance(sol, dict):
                    Pc_val = float(sol[Pc])
                    n4_val = float(sol[n4])
                elif isinstance(sol, tuple):
                    Pc_val = float(sol[0])
                    n4_val = float(sol[1])
                else:
                    continue
                
                # Check validity constraints
                if (0 <= Pc_val <= 1 and 0 <= n4_val <= total_sequences):
                    valid_solutions.append((Pc_val, n4_val))
            except (ValueError, TypeError):
                continue
                
    elif isinstance(solutions, dict):
        # Single solution as dictionary
        try:
            Pc_val = float(solutions[Pc])
            n4_val = float(solutions[n4])
            
            if (0 <= Pc_val <= 1 and 0 <= n4_val <= total_sequences):
                valid_solutions.append((Pc_val, n4_val))
        except (ValueError, TypeError):
            pass
    
    return valid_solutions
"""
# Validation Study 1: Parameter Recovery
print("VALIDATION STUDY 1: Parameter Recovery")
print("="*60)

# Test with known parameters
true_Pc = 0.9151860764850289 #0.8311920832252657 
true_n4 = 372 #107
true_n8 = 14 #7
total_seq = true_n4 + true_n8
model = TcellLineageModel(true_Pc, true_n4, total_seq,"a")

print(f"True parameters: Pc = {true_Pc}, n4 = {true_n4}")

# Generate multiple datasets and see if we can recover parameters
recovery_results = []
n_simulations = 100

for i in range(n_simulations):
    # Simulate data
    cd4_cd4, cd4_cd8, cd8_cd8 = model.simulate_dataset(seed=i)
    print(cd4_cd4, cd4_cd8, cd8_cd8)
    
    # Fit model to simulated data
    solutions = fit_model_to_data(cd4_cd4, cd8_cd8, total_seq)
    
    if len(solutions) >= 1:
        # Take the first valid solution
        fitted_Pc, fitted_n4 = solutions[1]
        fitted_n8 = model.total_sequences - fitted_n4
        
        print(f"n4: {fitted_n4}, n8: {fitted_n8}")
        print(fitted_Pc, fitted_n4)
        
        recovery_results.append({
            'sim_id': i,
            'true_Pc': true_Pc,
            'true_n4': true_n4,
            'fitted_Pc': fitted_Pc,
            'fitted_n4': fitted_n4,
            'cd4_cd4': cd4_cd4,
            'cd8_cd8': cd8_cd8,
            'error_Pc': abs(fitted_Pc - true_Pc),
            'error_n4': abs(fitted_n4 - true_n4),
            'error_cd4_cd4': abs(cd4_cd4 - 337)
        })
    else:
        print("*"*50)
        print(f"No valid solutions found for simulation {i}")
        print("*"*50)
        print(cd4_cd4, cd4_cd8, cd8_cd8)
        print("*"*50)

df_recovery = pd.DataFrame(recovery_results)
print(df_recovery)
print(f"Successful recoveries: {len(df_recovery)}/{n_simulations}")

if len(df_recovery) > 0:
    print(f"Mean Pc error: {df_recovery['error_Pc'].mean():.4f} ± {df_recovery['error_Pc'].std():.4f}")
    print(f"Mean n4 error: {df_recovery['error_n4'].mean():.4f} ± {df_recovery['error_n4'].std():.4f}")
    print(f"Mean cd4_cd4 error: {df_recovery['error_cd4_cd4'].mean():.4f} ± {df_recovery['error_cd4_cd4'].std():.4f}")
    df_recovery.to_csv("250912_b_parameter_recovery_pc_estimation.csv", index=False)

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
else:
    print("No successful parameter recoveries found!")
"""

# Validation Study 2: Your actual data
# print(f"\nVALIDATION STUDY 2: Your Observed Data")
# print("="*60)

# # Your observed data
# observed_cd4_cd4 = 337
# observed_cd8_cd8 = 10
# observed_total = 386

# # Fit to your data
# solutions = fit_model_to_data(observed_cd4_cd4, observed_cd8_cd8, observed_total)

# if len(solutions) > 0:
#     fitted_Pc, fitted_n4 = solutions[0]
#     print(f"Observed data: CD4-CD4 = {observed_cd4_cd4}, CD8-CD8 = {observed_cd8_cd8}")
#     print(f"Fitted parameters: Pc = {fitted_Pc:.6f}, n4 = {fitted_n4:.6f}")

#     # Create model with fitted parameters
#     fitted_model = TcellLineageModel(fitted_Pc, fitted_n4, observed_total)
#     expected_counts = fitted_model.expected_counts()

#     print(f"Expected counts from model:")
#     print(f"  CD4-CD4: {expected_counts[0]:.2f} (observed: {observed_cd4_cd4})")
#     print(f"  CD4-CD8: {expected_counts[1]:.2f} (observed: {observed_total - observed_cd4_cd4 - observed_cd8_cd8})")
#     print(f"  CD8-CD8: {expected_counts[2]:.2f} (observed: {observed_cd8_cd8})")
# else:
#     print("No valid solutions found for the observed data!")