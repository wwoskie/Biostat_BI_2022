# -*- coding: utf-8 -*-
# Pandas понадобится нам для чтения денных
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest

def ci_maker(population, sample_size, samples_mean):
    std = np.std(population)
    se = std / np.sqrt(sample_size)

    mean_mean = np.mean(samples_mean)
    left_b = mean_mean - 1.96 * se
    right_b = mean_mean + 1.96 * se

    return left_b, right_b

def check_intervals_intersect(first_ci, second_ci):
    if min(second_ci) < min(first_ci):
        second_ci, first_ci = first_ci, second_ci
    are_intersect = min(second_ci) < max(first_ci)
    return are_intersect # True or False

def check_dge_with_ci(first_table, second_table):
    # dge - differential gene expression
    ci_test_results = []

    for gene in first_table.select_dtypes('number').columns:
        first_ci = st.t.interval(alpha=0.95, # 95% доверительный интервал
              df=len(first_table[gene]) - 1, # число степеней свободы - 1
              loc=np.mean(first_table[gene]), # Среднее
              scale=st.sem(first_table[gene])) # Стандартная ошибка среднего

        second_ci = st.t.interval(alpha=0.95, # 95% доверительный интервал
              df=len(second_table[gene]) - 1, # число степеней свободы - 1
              loc=np.mean(second_table[gene]), # Среднее
              scale=st.sem(second_table[gene])) # Стандартная ошибка среднего

        ci_test_results.append(check_intervals_intersect(first_ci, second_ci))

    return ci_test_results

def check_dge_with_ztest(first_table, second_table):
    z_test_p_values = [round(ztest(first_table[gene], second_table[gene])[1], 4) for gene in first_table.select_dtypes('number').columns]
    z_test_results = np.array(z_test_p_values) < 0.05
    return z_test_p_values, z_test_results

def mean_diff_counter(first_table, second_table):
    mean_diff = [round(np.mean(first_table[gene]) - np.mean(second_table[gene]), 1) for gene in first_table.select_dtypes('number').columns]
    return mean_diff

def dge_ci_z_test(first_cell_type_expressions_path = 'data/first_table.csv', 
                  second_cell_type_expressions_path = 'data/second_table.csv', 
                  save_results_table = 'processed_data/results_table.csv'):
    
    first_table = pd.read_csv(first_cell_type_expressions_path, index_col=0)
    second_table = pd.read_csv(second_cell_type_expressions_path, index_col=0)

    ci_test_results = check_dge_with_ci(first_table, second_table)
    z_test_p_values, z_test_results = check_dge_with_ztest(first_table, second_table)
    mean_diff = mean_diff_counter(first_table, second_table)
    
    results = {
    "ci_test_results": ci_test_results,
    "z_test_results": z_test_results,
    "z_test_p_values": z_test_p_values,
    "mean_diff": mean_diff
    }

    results = pd.DataFrame(results)
    results.to_csv(save_results_table)

# вызов с импутом для реально неленивых людей


dge_ci_z_test(input('Type path to your first table: '), 
                  input('Type path to your second table: '), 
                  input('Type path to your results table: '))