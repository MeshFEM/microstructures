#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import itertools

import sys

import numpy as np


def print_experiment_info(params):
    print("\nExperimenting with parameters: \n")
    print(str(params))
    for i in range(0, len(params)):
        print("p" + str(i+1) + ": " + str(params[i]))


def create_name(pattern, params):
    name = partition_folder + "/" + pattern + "_"
    for i in range(0, len(params) - 1):
        name += "p" + str(i+1) + "-" + str(round(params[i], 3)) + "_"

    name += "p" + str(len(params)) + "-" + str(round(params[len(params) - 1], 3))

    return name


def parameters_string(params):
    params_string = ""
    for i in range(0, len(params) - 1):
        params_string += str(params[i]) + ", "

    params_string += str(params[len(params) - 1])

    return params_string


def parameter_table(params):
    params_string = str(len(params)) + " "
    for i in range(0, len(params) - 1):
        params_string += str(params[i]) + " "

    params_string += str(params[len(params) - 1])

    return params_string


def inflate_pattern(pattern, params, mesh_name, inflation_graph_radius=5):
    pathname = os.path.dirname(sys.argv[0])
    inflator_executable_path = pathname + '/../../cmake-build-release/isosurface_inflator/isosurface_cli'

    if os.path.isfile(mesh_name):
        print("Already computed")
    else:
        params_string = parameters_string(params)
        pattern_path = pathname + "/../../data/patterns/3D/reference_wires/pattern" + pattern + ".wire"

        cmd = [inflator_executable_path, args.symmetry, pattern_path, '--params', params_string, '-m',
               pathname + '/coarser_3d_meshing_opts.opt', '--cheapPostprocessing',
               '--inflation_graph_radius', str(inflation_graph_radius), mesh_name]
        print(cmd)
        try:
            result = subprocess.call(cmd, timeout=100)
            return bool(result == 0)
        except KeyboardInterrupt:
            raise
        except:
            print("Could not build mesh. Go to the following!")
            return False


def simulate_mesh(mesh_name, output_log, args):
    homogenization_executable_path = script_directory + '/../../../MeshFEM/cmake-build-release/MeshFEM/PeriodicHomogenization_cli'

    if args.material == 'B9Creator':
        material = script_directory + '/../../data/materials/B9Creator.material'
    else:
        material = script_directory + '/' + args.material

    cmd = [homogenization_executable_path, mesh_name, '-m', material]

    if args.symmetry == 'orthotropic' or args.symmetry == 'cubic':
        cmd += ['--ortho']

    print(cmd)
    with open(output_log, 'w') as out_log:
        try:
            result = subprocess.call(cmd, stdout=out_log)
            return bool(result == 0)
        except KeyboardInterrupt:
            raise
        except:
            print("Could not run simulation. Go to the following!")
            return False

    return True


def add_to_partition_table(mesh_name, log_file, params_table_string):
    out_log = open(log_file, 'r')

    float_pattern = re.compile(r'\-?\d+\.?\d*e?-?\d*')  # Compile a pattern to capture float values

    wire_content = out_log.readlines()
    for l, line in enumerate(wire_content):
        if line.startswith("Elasticity tensor:") or line.startswith("Homogenized elasticity tensor:"):
            tensor_values = []
            for r in range(1,7):
                row = [float(i) for i in float_pattern.findall(wire_content[l + r])]
                for v in row:
                    tensor_values.append(v)

    # Infer pattern given name
    path, filename = os.path.split(mesh_name)
    pattern_name = filename[:4]

    tensor_string = ' '.join(map(str, tensor_values))

    # Finally, write to table
    print('{} {} {}\n'.format(pattern_name, tensor_string, params_table_string))
    table_file.write('{} {} {}\n'.format(pattern_name, tensor_string, params_table_string))


def run_experiment(pattern, params, args):
    print_experiment_info(params)
    name = create_name(pattern, params)
    mesh_name = name + '.msh'
    simulation_log = name + '.sim'
    params_table_string = parameter_table(params)

    if params_table_string in partition_data:
        print("Already computed. Skip!")
        return

    result = inflate_pattern(pattern, params, mesh_name)

    if result and os.path.isfile(mesh_name):
        good_simulation = simulate_mesh(mesh_name, simulation_log, args)

        if good_simulation and os.path.isfile(simulation_log):
            add_to_partition_table(mesh_name, simulation_log, params_table_string)


    # Cleaning
    if os.path.isfile(mesh_name):
        os.remove(mesh_name)

        if os.path.isfile(simulation_log):
            os.remove(simulation_log)


if __name__== "__main__":
    parser = argparse.ArgumentParser(description='Sweep through different parameter sets of 3D topologies')
    parser.add_argument('topology', help='path to topology in which the sweep will be made')
    parser.add_argument('range_file', help='file with ranges to be used for each parameter')
    parser.add_argument('--symmetry', default='cubic', help='symmetry to be used in inflation of pattern')
    parser.add_argument('--partition', type=int, default=0, help='partition of experiments to be run')
    parser.add_argument('--total-partitions', type=int, default=1, help='number of experiment partitions to be used')
    parser.add_argument('--material', default='Default.material', help='choose material')
    args = parser.parse_args()

    pattern = args.topology
    partition_folder = "instances" + str(args.partition)
    if not os.path.exists(partition_folder):
        os.makedirs(partition_folder)

    pathname = os.path.dirname(sys.argv[0])
    script_directory = os.path.abspath(pathname)

    parameters_values = []

    # Read file, line by line, and add all values for each parameter
    range_file = open(args.range_file, 'r')
    float_pattern = re.compile(r'\-?\d+\.?\d*e?-?\d*')  # Compile a pattern to capture float values

    range_content = range_file.readlines()
    for l, line in enumerate(range_content):
        p_values = [float(i) for i in float_pattern.findall(line)]  # Convert strings to float

        parameters_values.append(np.array(p_values))

    # Collecting all experiments
    experiments = []
    for e in itertools.product(*parameters_values):
        experiments.append(e)

    num_experiments = len(experiments)
    each_partition = num_experiments / args.total_partitions

    print("Total number of experiments: " + str(num_experiments))
    print("Average number of experiments per partition: " + str(each_partition))

    # Find which range of experiment indices we should run
    start = args.partition * each_partition
    stop_before = start + each_partition

    start = int(start)
    stop_before = int(stop_before)

    # Create table and leave it open
    partition_table = "table" + str(args.partition) + ".txt"
    partition_data = []
    if os.path.isfile(partition_table):
        partition_data = open(partition_table).read()

    table_file = open(partition_table, 'a', buffering=1)

    print("Number of experiments for this partition: " + str(stop_before-start))
    print("Starting at index: " + str(start))
    print("Stopping at index: " + str(stop_before))

    # Inflate and simulate the corresponding set of parameters
    for e in range(start, stop_before):
        print("Index: " + str(e))
        experiment = list(experiments[e])
        print("Params: " + str(experiment))
        run_experiment(pattern, experiment, args)

    table_file.close()