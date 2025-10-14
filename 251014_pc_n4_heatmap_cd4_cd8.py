# SUBSTITUTING N4 AND LOOKING AT PC
import sympy as sp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm     

def calculate_exp(Pc_val, n4_val, all_seq, chain):
    """Calculate expected CD4-CD4 count for given Pc and n4"""
    if chain == "a":
        c44 = 10*9 / (15*14)
        c48 = 10*5 / (15*14)
        c88 = 5*4 / (15*14)
    elif chain == "b":
        c44 = 21*20 / (31*30)
        c48 = 21*10 / (31*30)
        c88 = 10*9 / (31*30)
    
    n8_val = all_seq - n4_val
    
    # For CD4 biased sequences
    p4_44_val = Pc_val**2 * c44
    p4_48_val = 2*Pc_val*(1-Pc_val) * c48
    p4_88_val = (1-Pc_val)**2 * c88
    p4_total_val = p4_44_val + p4_48_val + p4_88_val
    
    # For CD8 biased sequences
    p8_44_val = (1-Pc_val)**2 * c44
    p8_48_val = 2*Pc_val*(1-Pc_val) * c48
    p8_88_val = Pc_val**2 * c88
    p8_total_val = p8_44_val + p8_48_val + p8_88_val
    
    exp_44 = n4_val * (p4_44_val/p4_total_val) + n8_val * (p8_44_val/p8_total_val)
    exp_48 = n4_val * (p4_48_val/(p4_44_val+p4_48_val+p4_88_val)) + n8_val * (p8_48_val/(p8_44_val+p8_48_val+p8_88_val))
    exp_88 = n4_val * (p4_88_val/(p4_44_val+p4_48_val+p4_88_val)) + n8_val * (p8_88_val/(p8_44_val+p8_48_val+p8_88_val))
    
    return exp_44,exp_48,exp_88

def create_heatmap(all_seq, chain, observed, parameter):
    """Create heatmap of CD4-CD4 counts across n4 and Pc values"""
    # Create grid - use continuous values for both axes
    n4_values = np.linspace(0, all_seq, 1000)  # Continuous n4 values
    pc_values = np.linspace(0.5, 1.0, 1000)  # Continuous Pc values

    cd4_cd4_diff_grid = np.zeros((len(n4_values), len(pc_values)))
    cd4_cd8_diff_grid = np.zeros((len(n4_values), len(pc_values)))
    cd8_cd8_diff_grid = np.zeros((len(n4_values), len(pc_values)))
    
    diff_grids = [cd4_cd4_diff_grid,cd4_cd8_diff_grid,cd8_cd8_diff_grid]
    label = ["CD4-CD4","CD4-CD8","CD8-CD8"]

    for i, n4 in enumerate(n4_values):
        for j, pc in enumerate(pc_values):
            val44,val48,val88 = calculate_exp(pc, n4, all_seq, chain)
            #print(type(val44))
            cd4_cd4_diff_grid[i, j] = abs(observed[0] - val44)
            cd4_cd8_diff_grid[i, j] = abs(observed[1] - val48)
            cd8_cd8_diff_grid[i, j] = abs(observed[2] - val88)
    
    if parameter == "true_cd4":
        y_axis = n4_values
        org = "lower"
        ylabel = "True CD4"
    else: 
        if parameter == "true_cd8":
            y_axis = all_seq - n4_values
            org = "upper"
            ylabel = "True CD8"
        else:
            raise Exception("Parameter not valid")

    # Create heatmap
    #plt.figure(figsize=(12, 8))
    #print(cd4_cd4_diff_grid[:,1])
    # Use colormap where 0 (perfect match) is dark blue/green and large differences are red
    for i in range(0,len(diff_grids)):
        contour = plt.contour(pc_values, y_axis, diff_grids[i], 
                         levels=15, colors='black', linewidths=2)
        plt.clabel(contour, inline=True, fontsize=10)#, #fmt=f'{observed[i]}')
        contour = plt.contour(pc_values, y_axis, diff_grids[i], 
                         levels=[0.05], colors='black', linewidths=2)
        plt.clabel(contour, inline=True, fontsize=10)
        im = plt.imshow(diff_grids[i], 
                        aspect='auto',
                        origin=org,
                        extent=[pc_values.min(), pc_values.max(), y_axis.min(), y_axis.max()],
                    cmap='coolwarm_r', 
                    #vmin=0, vmax=diff_grids[i].max()
                    )
        print(diff_grids[i].min(), diff_grids[i].max())
        plt.colorbar(im, label=f'|Observed - Expected| {label[i]} Count')
        plt.xlabel('Pc', fontsize=12)
        plt.ylabel(ylabel, fontsize=12)
        plt.title(f'{label[i]} Difference from Observed value {observed[i]}', fontsize=14)
    
        # Add contour line for zero difference (exact match)
        
        #plt.clabel(contour, inline=True, fontsize=10, fmt='Perfect match')
    
        plt.tight_layout()
        #plt.show()
        name = str(chain)+"_0"+str(i+1)+"_"+parameter+"_heatmap.png"
        plt.savefig(name)
        plt.show()

##### ALPHA ############################################################################################
obs_alpha = [337,39,10]
all_seq_a = 386

# Generate original plots
#df_alpha = generate_table(all_seq_a,"a",obs_alpha)
#p#rint(df_alpha)
#draw_plots(df_alpha)

# Create heatmap
#n4_values = np.linspace(0, all_seq_a, 1000)

create_heatmap(all_seq_a, "a", obs_alpha,"true_cd4")
create_heatmap(all_seq_a, "a", obs_alpha,"true_cd8")
##### BETA #############################################################################################
obs_beta = [89,21,4]
all_seq_b = 114

#df_beta = generate_table(all_seq_b,"b",obs_beta)
#print(df_beta)
#draw_plots(df_beta)
create_heatmap(all_seq_b, "b", obs_beta,"true_cd4")
create_heatmap(all_seq_b, "b", obs_beta,"true_cd8")