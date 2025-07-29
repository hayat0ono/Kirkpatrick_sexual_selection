# 様々な初期値からGAシミュレーションを実行し、最終的な平衡状態をプロットする

import random
from pathlib import Path
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from model import Model

def run_simulation(t1_initial, p1_initial):
    model = Model(N, male_female_ratio, t1_initial, p1_initial, male_death_prob, s_ratio, female_death_prob, a0_coeff, a1_coeff, end_time, lifetime, num_child, mutation_rate)
    model.build_objects()
    while model.check_to_stop():
        model.step()
    return model.get_t1_ratio(), model.get_p1_ratio()

def plot_results(results, v1, v2, save_path=None):
        """
        Args:
            save_path (str): Save path (display only if None)
        """
        # make output directory if not exists
        save_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Extract data
        t1_finals = [r[0] for r in results]
        p1_finals = [r[1] for r in results]
        
        # Create figure
        plt.figure(figsize=(10, 8))
        
        # P1 vs T1 scatter plot
        plt.scatter(p1_finals, t1_finals, alpha=0.7, s=50, marker='s', 
                     facecolors='none', edgecolors='black', linewidth=1.5)
        
        # display theoretical equilibrium values
        p1_range = np.linspace(0, 1, 100)
        t1_equilibrium = []
        
        for p1 in p1_range:
            if p1 <= v1:
                t1_equilibrium.append(0.0)
            elif p1 < v2:
                t1_eq = (p1 - v1) / (v2 - v1)
                t1_equilibrium.append(max(0.0, min(1.0, t1_eq)))
            else:
                t1_equilibrium.append(1.0)
        
        # Plot theoretical equilibrium curve
        plt.plot(p1_range, t1_equilibrium, 'r--', linewidth=2, label='Theoretical Equilibrium')
        
        # Graph settings
        plt.xlabel('P1 Gene Frequency')
        plt.ylabel('T1 Gene Frequency')
        plt.title('P1 vs T1 Final Gene Frequencies')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        
        plt.show()

N = 131072
male_female_ratio = 0.5
male_death_prob = 0.0
s_ratio = 0.2
female_death_prob = 0.0
a0_coeff = 2.0
a1_coeff = 3.0
end_time = 500
lifetime = 1
num_child = 2
mutation_rate = 0.00001

results = []
for _ in range(51):
    t1_initial = random.uniform(0.1, 0.9)
    p1_initial = random.uniform(0.1, 0.9)
    t1_final, p1_final = run_simulation(t1_initial, p1_initial)
    results.append((t1_final, p1_final))

v1 = (a0_coeff + s_ratio - 1) / ((a0_coeff * a1_coeff - 1) * (1.0 - s_ratio))
v2 = a1_coeff * (a0_coeff + s_ratio - 1) / (a0_coeff * a1_coeff - 1)

plot_results(results, v1, v2, save_path=Path("output/scatter_plot_of_equilibrium_state.png"))