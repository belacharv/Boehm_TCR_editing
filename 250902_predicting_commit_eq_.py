import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.optimize import fsolve
import pandas as pd
from scipy import stats

class TcellLineageModel:
    def __init__(self, Pc, n4_total, total_sequences=386):
        self.Pc = Pc  # Commitment precision
        self.n4_total = n4_total  # Number of CD4-biased sequences
        self.n8_total = total_sequences - n4_total  # Number of CD8-biased sequences
        self.total_sequences = total_sequences
        print(total_sequences)
        print(self.n8_total)
        
        # Constants from your original model
        self.c44 = 10*9/(15*14)  # 3/7
        self.c48 = 10*5/(15*14)  # 5/21
        self.c88 = 5*4/(15*14)   # 2/21
    
    def simulate_dataset(self, seed=None):
        """Simulate a dataset based on the model parameters"""
        if seed is not None:
            np.random.seed(seed)
        
        cd4_cd4_count = 0
        cd4_cd8_count = 0
        cd8_cd8_count = 0
        
        # Simulate CD4-biased sequences
        for _ in range(self.n4_total):
            # For CD4-biased sequences
            p_44 = self.Pc**2 * self.c44
            p_48 = 2*self.Pc*(1-self.Pc) * self.c48
            p_88 = (1-self.Pc)**2 * self.c88
            
            # Normalize probabilities
            total_p = p_44 + p_48 + p_88
            probs = [p_44/total_p, p_48/total_p, p_88/total_p]
            
            # Sample outcome
            outcome = np.random.choice([0, 1, 2], p=probs)  # 0=CD4-CD4, 1=CD4-CD8, 2=CD8-CD8
            
            if outcome == 0:
                cd4_cd4_count += 1
            elif outcome == 1:
                cd4_cd8_count += 1
            else:
                cd8_cd8_count += 1

        print(f"Only CD4 biased\nCD4-CD4: {cd4_cd4_count}\nCD4-CD8: {cd4_cd8_count}\nCD8-CD8: {cd8_cd8_count}")
        
        # Simulate CD8-biased sequences
        for _ in range(self.n8_total):
            # For CD8-biased sequences
            p_44 = (1-self.Pc)**2 * self.c44
            p_48 = 2*self.Pc*(1-self.Pc) * self.c48
            p_88 = self.Pc**2 * self.c88
            
            # Normalize probabilities
            total_p = p_44 + p_48 + p_88
            probs = [p_44/total_p, p_48/total_p, p_88/total_p]
            
            # Sample outcome
            outcome = np.random.choice([0, 1, 2], p=probs)
            
            if outcome == 0:
                cd4_cd4_count += 1
            elif outcome == 1:
                cd4_cd8_count += 1
            else:
                cd8_cd8_count += 1
        print(f"Only CD8 biased\nCD4-CD4: {cd4_cd4_count}\nCD4-CD8: {cd4_cd8_count}\nCD8-CD8: {cd8_cd8_count}")
        return cd4_cd4_count, cd4_cd8_count, cd8_cd8_count
        

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
        
        return expected_cd4_cd4, expected_cd4_cd8, expected_cd8_cd8

def fit_model_to_data(cd4_cd4_obs, cd8_cd8_obs, total_sequences=386):
    """Fit the model to observed data using the same approach as your original code"""
    
    def equations(vars):
        Pc, n4 = vars
        
        # Constants
        c44 = 10*9/(15*14)
        c48 = 10*5/(15*14)
        c88 = 5*4/(15*14)
        
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
        
        # Expected counts
        eq1 = n4 * (p4_44/p4_total) + (total_sequences - n4) * (p8_44/p8_total) - cd4_cd4_obs
        eq2 = n4 * (p4_88/p4_total) + (total_sequences - n4) * (p8_88/p8_total) - cd8_cd8_obs
        
        return [eq1, eq2]
    
    # Try different initial guesses
    for guess in [(0.5, 200), (0.8, 300), (0.2, 100)]:
        try:
            solution = fsolve(equations, guess)
            if 0 <= solution[0] <= 1 and 0 <= solution[1] <= total_sequences:
                return solution[0], solution[1]
        except:
            continue
    
    return None, None

# Validation Study 1: Parameter Recovery
print("VALIDATION STUDY 1: Parameter Recovery")
print("="*60)

# Test with known parameters
true_Pc = 0.915186
true_n4 = 372
model = TcellLineageModel(true_Pc, true_n4)

print(f"True parameters: Pc = {true_Pc}, n4 = {true_n4}")

# Generate multiple datasets and see if we can recover parameters
recovery_results = []
n_simulations = 100

for i in range(n_simulations):
    # Simulate data
    cd4_cd4, cd4_cd8, cd8_cd8 = model.simulate_dataset(seed=i)
    
    # Fit model to simulated data
    fitted_Pc, fitted_n4 = fit_model_to_data(cd4_cd4, cd8_cd8)
    
    if fitted_Pc is not None:
        recovery_results.append({
            'sim_id': i,
            'true_Pc': true_Pc,
            'true_n4': true_n4,
            'fitted_Pc': fitted_Pc,
            'fitted_n4': fitted_n4,
            'cd4_cd4': cd4_cd4,
            'cd8_cd8': cd8_cd8,
            'error_Pc': abs(fitted_Pc - true_Pc),
            'error_n4': abs(fitted_n4 - true_n4)
        })

df_recovery = pd.DataFrame(recovery_results)
print(f"Successful recoveries: {len(df_recovery)}/{n_simulations}")
print(f"Mean Pc error: {df_recovery['error_Pc'].mean():.4f} ± {df_recovery['error_Pc'].std():.4f}")
print(f"Mean n4 error: {df_recovery['error_n4'].mean():.4f} ± {df_recovery['error_n4'].std():.4f}")

# Validation Study 2: Your actual data
print(f"\nVALIDATION STUDY 2: Your Observed Data")
print("="*60)

# Your observed data
observed_cd4_cd4 = 337
observed_cd8_cd8 = 10
observed_total = 386

# Fit to your data
fitted_Pc, fitted_n4 = fit_model_to_data(observed_cd4_cd4, observed_cd8_cd8, observed_total)

print(f"Observed data: CD4-CD4 = {observed_cd4_cd4}, CD8-CD8 = {observed_cd8_cd8}")
print(f"Fitted parameters: Pc = {fitted_Pc:.6f}, n4 = {fitted_n4:.6f}")

# Create model with fitted parameters
fitted_model = TcellLineageModel(fitted_Pc, fitted_n4, observed_total)
expected_counts = fitted_model.expected_counts()

print(f"Expected counts from model:")
print(f"  CD4-CD4: {expected_counts[0]:.2f} (observed: {observed_cd4_cd4})")
print(f"  CD4-CD8: {expected_counts[1]:.2f} (observed: {observed_total - observed_cd4_cd4 - observed_cd8_cd8})")
print(f"  CD8-CD8: {expected_counts[2]:.2f} (observed: {observed_cd8_cd8})")