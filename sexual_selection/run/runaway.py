# 数個の初期値からGAシミュレーションを実行し、ランナウェイ効果を確認する

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
    result_t1_vs_p1 = set()
    result_male_survival_rate = set()
    while model.check_to_stop():
        if model.time % 50 == 0:
            t1_ratio = model.get_t1_ratio()
            p1_ratio = model.get_p1_ratio()
            result_t1_vs_p1.add((model.time, t1_ratio, p1_ratio))
        male_survival_rate = (model.count_t_male(0) + model.count_t_male(1) * (1 - model.s_ratio)) / (model.count_t_male(0) + model.count_t_male(1))
        result_male_survival_rate.add((model.time, male_survival_rate))
        model.step()
    return result_t1_vs_p1, result_male_survival_rate

def plot_results(results, v1, v2, save_dir=None):
    """
    Args:
        save_path (str): Save path (display only if None)
    """
    # make output directory if not exists
    save_dir.mkdir(parents=True, exist_ok=True)

    results_plt_t1_vs_p1 = []
    results_plt_male_survival_rate = []

    for result_t1_vs_p1, result_male_survival_rate, p1_initial in results:
        sorted_result_t1_vs_p1 = sorted(result_t1_vs_p1, key=lambda x: x[0])
        sorted_result_male_survival_rate = sorted(result_male_survival_rate, key=lambda x: x[0])
        t1_ratios = [r[1] for r in sorted_result_t1_vs_p1]
        p1_ratios = [r[2] for r in sorted_result_t1_vs_p1]
        generations = [r[0] for r in sorted_result_male_survival_rate]
        male_survival_rates = [r[1] for r in sorted_result_male_survival_rate]
        results_plt_t1_vs_p1.append((t1_ratios, p1_ratios, p1_initial))
        results_plt_male_survival_rate.append((generations, male_survival_rates, p1_initial))

    plt.figure(figsize=(10, 8))
    for t1_ratios, p1_ratios, p1_initial in results_plt_t1_vs_p1:
        plt.scatter(p1_ratios, t1_ratios, alpha=0.7, s=50, marker='s', 
                     facecolors='none', edgecolors='black', linewidth=1.5)
        plt.plot(p1_ratios, t1_ratios, label=f'Initial P1: {p1_initial}')
    
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
    
    if save_dir:
        save_path = save_dir / "runaway.png"
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    
    plt.show()
    plt.close()

    plt.figure(figsize=(10, 8))
    for generations, male_survival_rates, p1_initial in results_plt_male_survival_rate:
        plt.plot(generations, male_survival_rates, label=f'Initial P1: {p1_initial}')
    
    plt.xlabel('Generation')
    plt.ylabel('Male Survival Rate')
    plt.title('Male Survival Rate vs Generation')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.xlim(0, end_time)
    plt.ylim(0, 1.0)
    plt.tight_layout()

    if save_dir:
        save_path = save_dir / "male_survival_rate.png"
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    
    plt.show()
    plt.close()

N = 13107
male_female_ratio = 0.5
male_death_prob = 0.0
s_ratio = 0.2
female_death_prob = 0.0
a0_coeff = 2.5
a1_coeff = 2.0
end_time = 500
lifetime = 1
num_child = 2
mutation_rate = 0.0001

initial_values = [
    (0.0, 0.6),
    (0.0, 0.7),
    (0.0, 0.8),
]

results = []
for t1_initial, p1_initial in initial_values:
    result_t1_vs_p1, result_male_survival_rate = run_simulation(t1_initial, p1_initial)
    results.append((result_t1_vs_p1, result_male_survival_rate, p1_initial))

v1 = (a0_coeff + s_ratio - 1) / ((a0_coeff * a1_coeff - 1) * (1.0 - s_ratio))
v2 = a1_coeff * (a0_coeff + s_ratio - 1) / (a0_coeff * a1_coeff - 1)

plot_results(results, v1, v2, save_dir=Path("output"))