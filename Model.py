#!/usr/bin/env python
#-*- coding: utf-8 -*-

from random import expovariate, gauss
import sys
import time
from optparse import OptionParser
from optparse import OptionGroup
import matplotlib.pyplot as plt

#there is a dictionary of rates, called rates_dict[name] = rate, eg.
# rates_dict['dna_k'] = 1

#{{{ rates_dict  &  pop_cont_dict
class rates_dict:
    def __init__(self):
    ##Ruty's parameters which appear magically from data
        self.dict={}
        self.dict['dna_k'] = 0.011
        self.dict['rna_k'] = 0.22
        self.dict['protein_k'] = 0.11
        self.dict['dna_l'] = 0.055
        self.dict['rna_l'] = 0.02
        self.dict['protein_l'] = 0.0001
    def change_param(self, param, val):
        self.dict[param] = float(val)

def init_pop_dict(result):
    dna_status = result[1]
    rna_pop    = result[2]
    protein_pop = result[3]

    d = {}
    d['dna_k'] = 1
    d['dna_l'] = 1
    d['rna_k'] = dna_status
    d['rna_l'] = rna_pop
    d['protein_k'] = rna_pop
    d['protein_l'] = protein_pop
    return d
#}}}

#{{{ first_event
def first_event(pop_const_dict, rates_dict):
    #returns the first event among the rates
    min = 0
    min_param = 'NULL'
    for param in rates_dict:    #param, eg. 'dna_k'
        val = pop_const_dict[param] * rates_dict[param]
        if val != 0:
            sample = expovariate(val)
        else:
            sample = 0
        if (min > sample or min == 0) and sample != 0:
            min = sample
            min_param = param
    return min, min_param
#}}}

#{{{ step_aux  &  step
def step_aux(last_result, rates_dict):
    #since each step only depends on what came before it
    pop_const_dict = init_pop_dict(last_result)

    dna_status = last_result[1]
    rna_pop = last_result[2]
    protein_pop = last_result[3]

    #call first_event
    winning_rate = first_event(pop_const_dict, rates_dict)
    next_time = winning_rate[0]
    winning_rate = winning_rate[1]

    if winning_rate == 'dna_k' and dna_status == 0:
        next_result = (next_time, 1, rna_pop, protein_pop)
    elif winning_rate == 'dna_l' and dna_status == 1:
        next_result = (next_time, 0, rna_pop, protein_pop)
    elif winning_rate == 'rna_k' and dna_status == 1:
        next_result = (next_time, dna_status, rna_pop + 1, protein_pop)
    elif winning_rate == 'rna_l' and rna_pop > 0:
        next_result = (next_time, dna_status, rna_pop - 1, protein_pop)
    elif winning_rate == 'protein_k'and rna_pop != 0:
        next_result = (next_time, dna_status, rna_pop, protein_pop + 1)
    elif winning_rate == 'protein_l' and protein_pop > 0:
        next_result = (next_time, dna_status, rna_pop, protein_pop - 1)
    else:
        next_result = 'NULL'
    return next_result

def step(last_result, rates_dict):
    next_result = 'NULL'
    while next_result == 'NULL':
        next_result = step_aux(last_result, rates_dict)
    return next_result
#}}}

#{{{ single_cell (max_time, rates_dict, asym=0.5)
def format(a,b,c,d):
    return "%f\t%d\t%d\t%d\t" %(a,b,c,d)

def single_cell(max_time, rates_dict, asym=0.5):
    last_result = (0,0,0,2000)  # using initial protein = 2000 based on results from running the simulation with initial protein = 0 and asym = 1.
    time = 0
    division_t = gauss(90,10)

    while time < max_time:
        yield time, last_result[1], last_result[2], last_result[3]
        #no division
        #change pop
        last_result = step(last_result, rates_dict)
        time += last_result[0]

        if time >= division_t:
            #there is a division
            #rescale by asym
            last_result = (division_t, last_result[1], round(last_result[2] * asym), round(last_result[3] * asym))
            time = division_t
            division_t += gauss(90,10)
#}}}

#{{{ multiple_cell (max_time, param_name, param_range, rates_dict)
def multiple_cell(max_time, param_name, param_range, rates_dict):
    #This is designed to vary one paramater at a time
    #param_range is a list usually of the form, [ x * step_size for x in range(n) ]

    for param_instance in param_range:
        #ugly exception but this code usually won't be run more than a hundred or two hundred time depending on length of param_range
        if param_name == 'asym':
            yield single_cell(max_time, rates_dict.dict, param_instance)
        else:
            rates_dict.change_param(param_name, param_instance)
            yield single_cell(max_time, rates_dict.dict, asym=0.5)

#}}}

#{{{protein_avg(simulations_list, param_range):
def protein_avg(simulations_list, param_range):
    # consolidate multiple_cell data in averages and return [(param_value, average protein level)]
    # simulations_list = list of generators for simulation data

    sum = 0
    dataPoint_count = 0
    protein_averages = []
    for sim in simulations_list:
        for dataPoint in sim:
            sum += dataPoint[3]
            dataPoint_count += 1
        average = sum / dataPoint_count
        protein_averages.append(average)

    return zip(param_range, protein_averages)       # assuming that the simulations run correspond to ascending list of parameters
#}}}

#{{{ protein_CofVar(sim_gen_gen, param_range):
def mean(list_or_gen):
    total = 0
    N = 0
    for data in list_or_gen:
        total += data
        N += 1
    return total / N

def var(sim_gen):
    l_sqrs = []
    l = []
    for dataPoint in sim_gen:
        protein_level = dataPoint[3]
        l_sqrs.append(protein_level ** 2)
        l.append(protein_level)
    return mean(l_sqrs) - mean(l)


def protein_CofVar(sim_gen_gen, param_range):
    CofVars = []
    means = []
    for gen in sim_gen_gen:
        sim_list = [dataPoint for dataPoint in gen]
        mu = mean([dataPoint[3] for dataPoint in sim_list])
        CofVars.append(var(sim_list) / mu)
        means.append(mu)
    return zip(param_range, means)
#}}}

#{{{ plot_data_tuple(data_list)
def plot_data_tuple(data_list, xLabel="time (min)", yLabel="Proteins (no. of molecules)"):
    # data_list = a list of tuples:    (next_time, dna_status, rna_pop, protein_pop)

    times = [float(dataPoint[0]) for dataPoint in data_list]
    proteins = [dataPoint[3] for dataPoint in data_list]

    plt.plot(times, proteins, 'r.')
    plt.axis([0, times[-1], 0, proteins[-1]*2])
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.show()
#}}}

#{{{ plot_data_stats(sim_list, xLabel, yLabel="Average Protein Count (no. of molecules)"):
def plot_data_stats(sim_list, xLabel, yLabel="Average Protein Count (no. of molecules)"):
    param_vals = [tuple[0] for tuple in sim_list]
    avg_vals = [tuple[1] for tuple in sim_list]

    param_vals, avg_vals = zip(*sim_list)
    param_vals = list(param_vals)
    avg_vals = list(avg_vals)

    plt.plot(param_vals, avg_vals, 'r.')
    plt.axis([0, param_vals[-1], 0, avg_vals[-1]*2])
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.show()
#}}}

#{{{ main
def main():
    d = rates_dict()

    usage = "%prog [-s -m] max_time(minutes) [options depend on simulation type see below]"
    parser = OptionParser(usage)

    parser.add_option("-s", "--single", action="store_true", dest="single",
            help="single simulation.      Usage: max_time(minutes) [parameters]")
    parser.add_option("-m", "--multiple", action="store_true", dest="multiple",
            help="multiple simulation.    Usage: max_time(minutes) --param_name --min --max --step_size.")

    parser.add_option("-g", "--graph", action="store_true", dest="graph",
            help="graph the results from the simulation.")

    single = OptionGroup(parser, "single simulation", "Run a single simulation for a set of parameters (given below --dna_k, etc.).  Choose parameters to differ from Ruty Rinott's original numbers: dna_k = 0.011 rna_k = 0.22 protein_k = 0.11 dna_l = 0.055 rna_l = 0.02 protein_l = 0.0001.  Asym is the amount of asymmetry in cell division (0 <= asym <= 1)")
    single.add_option("--dna_k")
    single.add_option("--rna_k")
    single.add_option("--protein_k")
    single.add_option("--dna_l")
    single.add_option("--rna_l")
    single.add_option("--protein_l")
    single.add_option("--asym")

    multiple = OptionGroup(parser, "multiple simulation", "Run multiple single cell simulations over a range of a given parameter.  The rest of the parameters are as given by Ruty Rinott (see above).  NB: The \"options\" given below are required")
    multiple.add_option("--param_name", help="Parameter to vary over")
    multiple.add_option("--min", help="Minimum value for vary parameter")
    multiple.add_option("--max", help="Maximum value for vary parameter")
    multiple.add_option("--step_size", help="Step size to vary parameter by between the minimum value and the maximum value given above")
    multiple.add_option("--cofvar", action="store_true", help="For graphing only: plot the coefficient of variation (var / mean) against mean protein expression values.  Default is to plot average protein levels against parameter value")

    parser.add_option_group(single)
    parser.add_option_group(multiple)

    (options, args) = parser.parse_args()
    if len(args) == 0:
        parser.print_help()
    if options.single and options.multiple:
        parser.error("Can't run --single and --multiple args mutually exclusive")

    if options.single:
        if options.dna_k:
            d.change_param("dna_k", options.dna_k)
        if options.rna_k:
            d.change_param("rna_k", options.rna_k)
        if options.protein_k:
            d.change_param("protein_k", options.protein_k)
        if options.dna_l:
            d.change_param("dna_l", options.dna_l)
        if options.rna_l:
            d.change_param("rna_l", options.rna_l)
        if options.protein_l:
            d.change_param("protein_l", options.protein_l)
        asym = 0.5
        if options.asym:
            asym = float(options.asym)
        if options.graph:
            plot_data_tuple([x for x in single_cell(int(args[0]), d.dict, asym)])
            return
        for x in single_cell(int(args[0]), d.dict, asym):
            print x

    if options.multiple:
        if not options.step_size:
            parser.error("See required options. Sorry for the nonsense \"required options\"")
        if not options.min:
            parser.error("See required options. Sorry for the nonsense \"required options\"")
        if not options.max:
            parser.error("See required options. Sorry for the nonsense \"required options\"")
        if not options.param_name:
            parser.error("See required options. Sorry for the nonsense \"required options\"")

        param_range = [x * float(options.step_size) for x in range(int(float(options.min)/float(options.step_size)), int(float(options.max)/float(options.step_size)))]

        if options.graph:
            if options.cofvar:
                plot_data_stats(protein_CofVar(multiple_cell(int(args[0]), options.param_name, param_range, d), param_range),"mean expression", "CV ( mean / var )")
                return
            else:
                plot_data_stats(protein_avg([x for x in multiple_cell(int(args[0]), options.param_name, param_range, d)], param_range),options.param_name)
                return

        # otherwise print out the data
        for x in multiple_cell(int(args[0]), options.param_name, param_range, d):
            for y in x:
                print y
#}}}


if __name__ == "__main__":
    main()
