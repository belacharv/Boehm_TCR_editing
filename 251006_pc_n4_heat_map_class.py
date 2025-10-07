import sympy as sp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


class TcellSeqModel:
    def __init__(self,chain):
        self.chain = chain
        if chain == "a":
            self.c44 = sp.Rational(10*9, 15*14)   # 3/7
            self.c48 = sp.Rational(10*5, 15*14)   # 5/21
            self.c88 = sp.Rational(5*4, 15*14)
        elif chain == "b":
            self.c44 = sp.Rational(21*20, 31*30)
            self.c48 = sp.Rational(21*10, 31*30)
            self.c88 = sp.Rational(10*9, 31*30)
        else:
            raise Exception("Wrong input, uncorrect TCR chain selection")
    
    def returnCoeff(self):
        return self.c44,self.c48,self.c88

    def returnProbab(self):
        Pc = sp.symbols('Pc')
        c44 = self.c44
        c48 = self.c48
        c88 = self.c88

        # For CD4 biased sequences  
        p4_44 = Pc**2 * c44 # both of the T cells chose the right fate (according to the bias of the sequence)
        p4_48 = 2*Pc*(1-Pc) * c48 # one was correct, one was wrong
        p4_88 = (1-Pc)**2 * c88 # both were wrong
        p4_total = p4_44 + p4_48 + p4_88
        P4 = [p4_44,p4_48,p4_88]

        # For CD8 biased sequences
        p8_44 = (1-Pc)**2 * c44
        p8_48 = 2*Pc*(1-Pc) * c48
        p8_88 = Pc**2 * c88
        p8_total = p8_44 + p8_48 + p8_88
        P8 = [p8_44,p8_48,p8_88]
        return P4,P8
    
def calculate_exp(Pc_val, n4_val, all_seq, chain):
    obj = TcellSeqModel(chain)
    p4,p8 = obj.returnProbab()
    n8_val = all_seq - n4_val
  
    exp_44 = n4_val * (p4[0]/(p4[0]+p4[1]+p4[2])) + n8_val * (p8[0]/(p8[0]+p8[1]+p8[2]))
    exp_48 = n4_val * (p4[1]/(p4[0]+p4[1]+p4[2])) + n8_val * (p8[1]/(p8[0]+p8[1]+p8[2]))
    exp_88 = n4_val * (p4[2]/(p4[0]+p4[1]+p4[2])) + n8_val * (p8[2]/(p8[0]+p8[1]+p8[2]))
    
    return exp_44,exp_48,exp_88

def create_heatmap(all_seq, chain, observed):
    """Create heatmap of CD4-CD4 counts across n4 and Pc values"""
    # Create grid - use continuous values for both axes
    n4_values = np.linspace(0, all_seq, 1000)  # Continuous n4 values
    pc_values = np.linspace(0.0, 1.0, 1000)  # Continuous Pc values
    
    # Calculate CD4-CD4 difference for each combination
    cd4_cd4_diff_grid = np.zeros((len(n4_values), len(pc_values)))
    cd4_cd8_diff_grid = np.zeros((len(n4_values), len(pc_values)))
    cd8_cd8_diff_grid = np.zeros((len(n4_values), len(pc_values)))
    
    diff_grids = [cd4_cd4_diff_grid,cd4_cd8_diff_grid,cd8_cd8_diff_grid]
    label = ["CD4-CD4","CD4-CD8","CD8-CD8"]

    for i, n4 in enumerate(n4_values):
        for j, pc in enumerate(pc_values):
            val44,val48,val88 = calculate_exp(pc, n4, all_seq, chain)
            print(type(val44))
            cd4_cd4_diff_grid[i, j] = abs(observed[0] - val44)
            cd4_cd8_diff_grid[i, j] = abs(observed[1] - val48)
            cd8_cd8_diff_grid[i, j] = abs(observed[2] - val88)
    
    # Create heatmap
    plt.figure(figsize=(12, 8))
    
    # Use colormap where 0 (perfect match) is dark blue/green and large differences are red
    for i in range(0,len(diff_grids)):
        im = plt.imshow(diff_grids[i], aspect='auto', origin='lower',
                    extent=[pc_values.min(), pc_values.max(), n4_values.min(), n4_values.max()],
                    cmap='coolwarm_r', vmin=0, vmax=diff_grids[i].max())
    
        plt.colorbar(im, label=f'|Observed - Expected| {label[i]} Count')
        plt.xlabel('Pc', fontsize=12)
        plt.ylabel('True CD4 (n4)', fontsize=12)
        plt.title(f'CD4-CD4 Difference from Observed value {observed[i]}', fontsize=14)
    
    # Add contour line for zero difference (exact match)
    contour = plt.contour(pc_values, n4_values, cd4_cd4_diff_grid, 
                         levels=[0], colors='black', linewidths=2)
    plt.clabel(contour, inline=True, fontsize=10, fmt='Perfect match')
    
    plt.tight_layout()
    plt.show()

    
 
##### ALPHA ############################################################################################
obs_alpha = [337,39,10]
all_seq_a = 386

# Create heatmap
create_heatmap(all_seq_a, "a", obs_alpha)