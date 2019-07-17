#!/usr/bin/python

import evo_p
import time
import math
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

N_TESTS = 3
TEST_MULTI_IND = 0
TEST_LAMBDA_IND = 1
TEST_CONVERGE_IND = 2
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

if (run_tests[TEST_MULTI_IND]):
    print("Now trying multi-objective optimization")
    prob = evo_p.Problem(32, 1, 2)
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
    prob_rast = evo_p.Problem(32, 2, 1)
    prob_rast.set_phenotype_parameters(["real_(-2.56, 2.56)", "real_(-2.56, 2.56)"])
    pop_rast = prob_rast.initialize_population("ga.conf")

    prob_rast.set_fitness_function(rastrigin)
    parents = pop_rast.get_parents()
    for org in parents:
        start_org_xs.append(org.read_real(0))
        start_org_ys.append(org.read_real(1))
        start_org_zs.append( -rastrigin_2d(org.read_real(0), org.read_real(1)) + 0.75)

    prev_max_fitness = -1000000
    for i in range(pop_rast.get_num_gens()):
        pop_rast.evaluate()
        if i % PRINT_PERIOD == 0:
            org = pop_rast.get_best()[0]
            print("generation " + str(i) + ", " + str(org.get_phenotype()) + ", max_fitness = " + str(pop_rast.get_max_fitness()))
        pop_rast.iterate()

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
    conv = evo_p.PlateauCutoff(0.0001, 5)
    prob_rast = evo_p.Problem(32, 2, 1)
    prob_rast.set_phenotype_parameters(["real_(-2.56, 2.56)", "real_(-2.56, 2.56)"])
    prob_rast.set_fitness_function(rastrigin)
    pop_rast = prob_rast.initialize_population("ga.conf")
    results = pop_rast.run(conv, store_intermediate=True)
    #print("generation " + str(results["Generation"]) + ", " + str(results["Solution"]) + ", max_fitness = " + str(results["Fitness"]))
    print("Generation\tFitness\tSolution")
    for i in range(results["Generations"]):
        print( str(i) + "\t\t" + "{0:.2f}".format(results["Fitness"][i]) + "\t" + str(results["Solution"][i].get_phenotype()) )
    print("successfully used plateau-based convergence checking")
