#!/usr/bin/python

import evo_q
import time
import math
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

N_TESTS = 5
N_DIMS = 2
TEST_MULTI_IND = 0
TEST_LAMBDA_IND = 1
TEST_CONVERGE_IND = 2
TEST_NOISE_IND = 3
TEST_COMPARE_LATIN = 4
run_tests = [True for i in range(N_TESTS)]

if len(sys.argv) > 1:
    run_tests = [False for i in range(N_TESTS)]

for arg in sys.argv:
    if "multi" in arg:
        run_tests[TEST_MULTI_IND] = True
    if "lambda" in arg:
        run_tests[TEST_LAMBDA_IND] = True
    if "converge" in arg:
        run_tests[TEST_CONVERGE_IND] = True
    if "noise" in arg:
        run_tests[TEST_NOISE_IND] = True
    if "latin" in arg:
        run_tests[TEST_COMPARE_LATIN] = True

N_GENS = 20
PRINT_PERIOD = 5

org_xs = []
org_ys = []

final_objective_coords = [[], []]
first_front_coords = [[], []]
start_org_xs = []
start_org_ys = []
start_org_zs = []
end_org_xs = []
end_org_ys = []
end_org_zs = []

A = 10
noise_std_dev = A/5

class TableFormatter:
    def __init__(self, header):
        self.header = header
        self.head_lengths = [len(val) for val in header]
        self.data_rows = []

    def add_row(self, data):
        if len(data) != len(self.head_lengths):
            raise ValueError('length of row does not match length of header')
        self.data_rows.append( data.copy() )
        for i, length in enumerate(self.head_lengths):
            if len(data[i]) > length:
                self.head_lengths[i] = len(data[i])
                if (self.head_lengths[i] - len(self.header[i])) % 2 == 1:
                    self.header[i] += " "
                    self.head_lengths[i] += 1

    def fill_spaces(self, col_i, val):
        ret = ""
        dl = self.head_lengths[col_i] - len(val)
        dl_half = int(dl/2)
        for j in range(max(0, dl_half)):
            ret += " "
        ret += val
        for j in range(max(0, dl_half)):
            ret += " "
        if dl % 2 == 1:
            ret += " "
        return ret

    def get_col_separator(self, spaces_between, print_bars):
        ret = ""
        for i in range(spaces_between):
            ret += " "
        if print_bars:
            ret += "|"
        for i in range(spaces_between):
            ret += " "
        return ret

    def get_row_separator(self, spaces_between, print_bars):
        ret = ""
        for i, length in enumerate(self.head_lengths):
            for i in range(length):
                ret += "-"
            if i != len(self.head_lengths) - 1:
                for j in range(spaces_between):
                    ret += "-"
                if print_bars:
                    ret += "+"
                for j in range(spaces_between):
                    ret += "-"
        return ret

    def show(self, spaces_between=1, print_col_bars=False, print_row_bars=False):
        pstr = ""
        for i, val in enumerate(self.header):
            pstr += self.fill_spaces(i, val)
            if i != len(self.header) - 1:
                pstr += self.get_col_separator(spaces_between, print_col_bars)
        print(pstr)
        if print_row_bars:
            print(self.get_row_separator(spaces_between, print_col_bars))
        for row_arr in self.data_rows:
            pstr = ""
            for i, val in enumerate(row_arr):
                pstr += self.fill_spaces(i, val)
                if i != len(self.header) - 1:
                    pstr += self.get_col_separator(spaces_between, print_col_bars)
            print(pstr)
            if print_row_bars:
                print(self.get_row_separator(spaces_between, print_col_bars))

def schaffer_fitness(org):
    parents = pop.get_parents()
    for org in parents:
        x1 = org.read_real(0)
        fits = [-x1*x1, -(x1 - 2)*(x1 - 2)]
        org.set_fitness(fits)

def rastrigin_2d(x, y):
    ret = 2*A
    ret += x*x - A*math.cos(2*math.pi*x)
    ret += y*y - A*math.cos(2*math.pi*y)
    return ret

def noisy_sphere(org):
    ret = np.random.normal(0, noise_std_dev)
    for i in range(N_DIMS):
        ret -= org.read_real(i)*org.read_real(i)
    return [ret]

def rastrigin_2d_grad(x, y):
    h = 0.01
    step_dist = 0.1
    for t in range(10):
        df_dx = (rastrigin_2d(x + h, y) - rastrigin_2d(x - h, y))/(2*h)
        df_dy = (rastrigin_2d(x, y + h) - rastrigin_2d(x, y - h))/(2*h)
        norm = math.sqrt(df_dx*df_dx + df_dy*df_dy)
        x -= step_dist*df_dx/norm
        y -= step_dist*df_dy/norm
    return rastrigin_2d(x,y)

def rastrigin_2d_np(x, y):
    return 2*A + x ** 2 - A*np.cos(2*math.pi*x) + y ** 2 - A*np.cos(2*math.pi*y)

schaffer = lambda org: [-org.read_real(0)*org.read_real(0), - (org.read_real(0) - 2)*(org.read_real(0) - 2)]
rastrigin = lambda org: [-rastrigin_2d(org.read_real(0), org.read_real(1))]
rastrigin_cost = lambda org: [rastrigin_2d(org.read_real(0), org.read_real(1))]

if (run_tests[TEST_MULTI_IND]):
    print("Now trying multi-objective optimization")
    prob = evo_q.Problem(32, 1, 2)
    print("successfully called constructor")
    prob.set_phenotype_parameters(["real_(-10.0, 10.0)"])
    print("successfully loaded phenotype")
    prob.set_template_parameter(0, 2.1)
    print("successfully set template parameter")
    pop = prob.initialize_population("ga.conf")
    prob.set_fitness_function(schaffer)
    print("successfully initialized the population")
    parents = pop.get_parents()
    print("successfully retrieved the parents in the population len = " + str(len(parents)))

    for i in range(N_GENS):
        #evaluate fitness for each organism
        pop.evaluate()
        if i % PRINT_PERIOD == 0:
            first_front = pop.get_best()
            first_front_coords = [[], []]
            for org in first_front:
                fits = org.get_fitness_list()
                first_front_coords[0].append(org.read_real(0)*org.read_real(0))
                first_front_coords[1].append((org.read_real(0) - 2)*(org.read_real(0) - 2))
            plt.scatter(first_front_coords[0], first_front_coords[1], s=20, c='blue', alpha=0.5)
            title_str = 'Objective space of schaffer optimization after ' + str(i) + ' generations'
            plt.title(title_str)
            plt.xlabel(r'$x^2$')
            plt.ylabel(r'$(x-2)^2$')
            plt.show()
        pop.iterate()
    pop.evaluate()
    parents = pop.get_best()
    for org in parents:
        fits = org.get_fitness_list()
        final_objective_coords[0].append(fits[0])
        final_objective_coords[1].append(fits[1])

    plt.scatter(final_objective_coords[0], final_objective_coords[1], s=20, c='red', alpha=0.5)
    plt.title('Objective space of schaffer optimization after optimization')
    plt.xlabel(r'$x^2$')
    plt.ylabel(r'$(x-2)^2$')
    plt.show()

    print("successfully used python bound evolution")

if run_tests[TEST_LAMBDA_IND]:
    print("Now running a test using lambda functions")
    #setup the rastrigin problem
    prob_rast = evo_q.Problem(32, 2, 1)
    prob_rast.set_phenotype_parameters(["real_(-2.56, 2.56)", "real_(-2.56, 2.56)"])
    pop_rast = prob_rast.initialize_population("ga.conf")

    prob_rast.set_fitness_function(rastrigin)
    parents = pop_rast.get_parents()
    for org in parents:
        start_org_xs.append(org.read_real(0))
        start_org_ys.append(org.read_real(1))
        start_org_zs.append( -rastrigin_2d(org.read_real(0), org.read_real(1)) + 0.75)

    prev_max_fitness = -1000000
    table = TableFormatter(["Generation", "Phenotype", "Max Fitness"])
    for i in range(pop_rast.get_num_gens()):
        pop_rast.evaluate()
        if i % PRINT_PERIOD == 0:
            org = pop_rast.get_best()[0]
            table.add_row([str(i), str(org.get_phenotype()), str(pop_rast.get_max_fitness())])
        pop_rast.iterate()
    table.show()

    pop_rast.evaluate()
    parents = pop_rast.get_parents()
    for org in parents:
        end_org_xs.append(org.read_real(0))
        end_org_ys.append(org.read_real(1))
        end_org_zs.append( -rastrigin_2d(org.read_real(0), org.read_real(1)) + 0.75)

    x = np.linspace(-2.56, 2.56)
    y = np.linspace(-2.56, 2.56)
    xx_land, yy_land = np.meshgrid(x, y)
    xx_start, yy_start = np.meshgrid(start_org_xs, start_org_ys)
    xx_end, yy_end = np.meshgrid(end_org_xs, end_org_ys)

    zz_land = -rastrigin_2d_np(xx_land, yy_land)
    zz_start = -rastrigin_2d_np(xx_start, yy_start) + 0.75
    zz_end = -rastrigin_2d_np(xx_end, yy_end) + 0.75

    fs = 25
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.grid(False)
    ax.set_xlabel(r'$\mathrm{x}$', fontsize=fs, labelpad=12)
    ax.set_ylabel(r'$\mathrm{y}$', fontsize=fs, labelpad=12)
    ax.set_zlabel(r'$C$', fontsize=fs, labelpad=12)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.plot_surface(xx_land, yy_land, zz_land, cmap='winter', alpha=0.24)
    ax.scatter3D(start_org_xs, start_org_ys, start_org_zs, c='red', s=[20 for n in range(len(xx_start))])
    plt.show()

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.grid(False)
    ax.set_xlabel(r'$\mathrm{x}$', fontsize=fs, labelpad=12)
    ax.set_ylabel(r'$\mathrm{y}$', fontsize=fs, labelpad=12)
    ax.set_zlabel(r'$C$', fontsize=fs, labelpad=12)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.plot_surface(xx_land, yy_land, zz_land, cmap='winter', alpha=0.5)
    ax.scatter3D(end_org_xs, end_org_ys, end_org_zs, c='red', s=[20 for n in range(len(xx_end))])
    plt.show()

    print("successfully used lambda function version")

if run_tests[TEST_CONVERGE_IND]:
    conv = evo_q.PlateauCutoff(0.0001, 5)
    prob_rast = evo_q.Problem(32, 2, 1)
    prob_rast.set_phenotype_parameters(["real_(-2.56, 2.56)", "real_(-2.56, 2.56)"])
    prob_rast.set_fitness_function(rastrigin)
    pop_rast = prob_rast.initialize_population("ga.conf")
    results = pop_rast.run(conv, store_intermediate=True)
    #print("generation " + str(results["Generation"]) + ", " + str(results["Solution"]) + ", max_fitness = " + str(results["Fitness"]))
    print("Generation\tFitness\tSolution")
    table = TableFormatter(["Generation", "Fitness", "Solution"])
    for i in range(results["Generations"]):
        table.add_row([str(i), "{:.2f}".format(results["Fitness"][i]), str(results["Solution"][i].get_phenotype())])
    table.show()
    print("successfully used plateau-based convergence checking")

HIST_BINS = 10
HIST_RANGE = [-100, 100]
N_SAMPLES = 100
if run_tests[TEST_NOISE_IND]:
    conv = evo_q.VarianceCutoff(0.01)
    prob = evo_q.Problem(32, 2, 1)
    prob.set_phenotype_parameters(["real_(-{}, {})".format(A, A) for i in range(N_DIMS)])
    prob.set_fitness_function(noisy_sphere)
    pop = prob.initialize_population("ga.conf")
    pop.set_noise_compensation(5)
    #results = pop.run(conv, store_intermediate=True)
    #results = pop.run(store_intermediate=True)
    print()
    results_noisy = []
    results_actual = []
    table = TableFormatter(["Generation", "Fitness", "Noiseless", "Solution"])
    for i in range(N_GENS):
        pop.evaluate()
        org = pop.get_best()[0]
        x = org.read_real(0)
        y = org.read_real(1)
        table.add_row([str(i), "{:.2f}".format(pop.get_max_fitness()), "{:.2f}".format(x*x + y*y), "{:.2f}, {:.2f}".format(x, y)])
        pop.iterate()
        results_noisy.append(-org.get_fitness())
        results_actual.append(x*x + y*y)
    table.show()
    plt.scatter(range(len(results_noisy)), results_noisy, label='organism fitness from noisy sampling')
    plt.scatter(range(len(results_actual)), results_actual, label='organism fitness without noise')
    plt.title('Noisy Sphere Convergence Plot')
    plt.legend()
    plt.xlabel(r'$x^2$')
    plt.ylabel(r'$y^2$')
    plt.show()
    print("successfully evaluated noisy fitness function")

if run_tests[TEST_COMPARE_LATIN]:
    conv_non_latin = evo_q.PlateauCutoff(0.0001, 5)
    conv_with_latin = evo_q.PlateauCutoff(0.0001, 5)
    prob_rast = evo_q.Problem(32, 2, 1)
    prob_rast.set_phenotype_parameters(["real_(-2.56, 2.56)", "real_(-2.56, 2.56)"])
    prob_rast.set_fitness_function(rastrigin)
    pop_rast_non_latin = prob_rast.initialize_population("ga.conf", False)
    pop_rast_with_latin = prob_rast.initialize_population("ga.conf", True)
    #pop_rast_non_latin.set_cost()
    #pop_rast_with_latin.set_cost()
    results_non_latin = pop_rast_non_latin.run(conv_non_latin, store_intermediate=True)
    results_with_latin = pop_rast_with_latin.run(conv_with_latin, store_intermediate=True)
    #print("generation " + str(results["Generation"]) + ", " + str(results["Solution"]) + ", max_fitness = " + str(results["Fitness"]))
    plt.scatter(range(results_non_latin["Generations"]), results_non_latin["Fitness"], label='without latin-hypercube initialization')
    plt.scatter(range(results_with_latin["Generations"]), results_with_latin["Fitness"], label='with latin-hypercube initialization')
    plt.title('Comparison between latin-hypercube and uniform random initialization')
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.legend()
    plt.show()
