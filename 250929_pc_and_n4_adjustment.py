import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def derive_n4_from_pc_constraint(chain="a"):
    """
    Derive the relationship between n4 and Pc such that 
    expected CD4-CD4 count equals 337 (for total_sequences = 386)
    
    Returns a symbolic expression and a numerical function
    """
    # Define symbolic variables
    Pc, n4 = sp.symbols('Pc n4', real=True, positive=True)
    
    # Constants based on chain type
    if chain == "a":
        c44 = sp.Rational(10*9, 15*14)  # 3/7
        c48 = sp.Rational(10*5, 15*14)  # 5/21
        c88 = sp.Rational(5*4, 15*14)   # 2/21
    elif chain == "b":
        c44 = sp.Rational(21*20, 31*30)
        c48 = sp.Rational(21*10, 31*30)
        c88 = sp.Rational(10*9, 31*30)
    
    total_sequences = 386
    target_cd4_cd4 = 337
    
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
    
    # Expected CD4-CD4 count equation
    expected_cd4_cd4 = n4 * (p4_44/p4_total) + (total_sequences - n4) * (p8_44/p8_total)
    
    # Solve for n4 as a function of Pc
    constraint_eq = expected_cd4_cd4 - target_cd4_cd4
    n4_solution = sp.solve(constraint_eq, n4)
    
    print("Symbolic solution for n4 as function of Pc:")
    print(f"n4 = {n4_solution[0]}")
    print("\nSimplified:")
    n4_simplified = sp.simplify(n4_solution[0])
    print(f"n4 = {n4_simplified}")
    
    # Create a numerical function
    n4_func = sp.lambdify(Pc, n4_solution[0], 'numpy')
    
    return n4_solution[0], n4_func


def calculate_n4_for_pc(Pc_value, chain="a"):
    """
    Calculate the required n4 value for a given Pc to achieve 
    expected CD4-CD4 count of 337
    """
    _, n4_func = derive_n4_from_pc_constraint(chain)
    n4_value = n4_func(Pc_value)
    return n4_value


def plot_pc_n4_relationship(chain="a"):
    """
    Plot the relationship between Pc and n4 for the constraint
    """
    _, n4_func = derive_n4_from_pc_constraint(chain)
    
    # Generate Pc values
    pc_values = np.linspace(0.5, 1.0, 200)
    n4_values = n4_func(pc_values)
    
    # Filter valid n4 values (between 0 and 386)
    valid_mask = (n4_values >= 0) & (n4_values <= 386)
    
    plt.figure(figsize=(10, 6))
    plt.plot(pc_values[valid_mask], n4_values[valid_mask], 'b-', linewidth=2)
    
    # Mark your original point
    original_pc = 0.9151860764850289
    original_n4 = 372
    plt.plot(original_pc, original_n4, 'ro', markersize=10, 
             label=f'Original: Pc={original_pc:.4f}, n4={original_n4}')
    
    plt.xlabel('Commitment Precision (Pc)', fontsize=12)
    plt.ylabel('Number of CD4-biased sequences (n4)', fontsize=12)
    plt.title('Relationship between Pc and n4 for CD4-CD4 count = 337', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    return pc_values[valid_mask], n4_values[valid_mask]


def verify_constraint(Pc_value, n4_value, chain="a"):
    """
    Verify that a given (Pc, n4) pair produces the expected CD4-CD4 count of 337
    """
    if chain == "a":
        c44 = 10*9/(15*14)
        c48 = 10*5/(15*14)
        c88 = 5*4/(15*14)
    elif chain == "b":
        c44 = 21*20/(31*30)
        c48 = 21*10/(31*30)
        c88 = 10*9/(31*30)
    
    total_sequences = 386
    n8_value = total_sequences - n4_value
    
    # For CD4-biased sequences
    p4_44 = Pc_value**2 * c44
    p4_48 = 2*Pc_value*(1-Pc_value) * c48
    p4_88 = (1-Pc_value)**2 * c88
    p4_total = p4_44 + p4_48 + p4_88
    
    # For CD8-biased sequences
    p8_44 = (1-Pc_value)**2 * c44
    p8_48 = 2*Pc_value*(1-Pc_value) * c48
    p8_88 = Pc_value**2 * c88
    p8_total = p8_44 + p8_48 + p8_88
    
    expected_cd4_cd4 = n4_value * (p4_44/p4_total) + n8_value * (p8_44/p8_total)
    expected_cd4_cd8 = n4_value * (p4_48/p4_total) + n8_value * (p8_48/p8_total)
    expected_cd8_cd8 = n4_value * (p4_88/p4_total) + n8_value * (p8_88/p8_total)
    
    print(f"\nVerification for Pc={Pc_value:.6f}, n4={n4_value:.2f}, n8={n8_value:.2f}")
    print(f"Expected CD4-CD4: {expected_cd4_cd4:.4f} (target: 337)")
    print(f"Expected CD4-CD8: {expected_cd4_cd8:.4f}")
    print(f"Expected CD8-CD8: {expected_cd8_cd8:.4f}")
    print(f"Total: {expected_cd4_cd4 + expected_cd4_cd8 + expected_cd8_cd8:.4f}")
    
    return expected_cd4_cd4, expected_cd4_cd8, expected_cd8_cd8


# Example usage
if __name__ == "__main__":
    print("="*60)
    print("Deriving n4 as a function of Pc (chain alpha)")
    print("="*60)
    
    n4_symbolic, n4_func = derive_n4_from_pc_constraint(chain="a")
    
    print("\n" + "="*60)
    print("Testing with different Pc values:")
    print("="*60)
    
    test_pc_values = [0.6,0.7,0.8,0.88,0.89,0.9,0.91, 0.9151860764850289, 0.92, 0.93]
    
    for pc in test_pc_values:
        n4 = calculate_n4_for_pc(pc, chain="a")
        print(f"\nPc = {pc:.6f} → n4 = {n4:.4f}")
        verify_constraint(pc, n4, chain="a")
    
    print("\n" + "="*60)
    print("Plotting relationship...")
    print("="*60)
    plot_pc_n4_relationship(chain="a")