#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import textwrap
import random
import statistics
from random import randint
import matplotlib
import networkx as nx
import matplotlib.pyplot as plt
random.seed(9001)
matplotlib.use("Agg")

__author__ = "Meryam Boulayat"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Meryam Boulayat"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Meryam Boulayat"
__email__ = "meryam.boulayat@etu.u-paris.fr"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file doesn't exist

    :return: (str) Path
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that yields the read sequences.
    """
    with open(fastq_file, "r") as file:
        for line in file:
            seq = next(file).strip()
            yield seq
            next(file)
            next(file)


def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.

    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    seq_length = len(read)
    for i in range(seq_length-(kmer_size-1)):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    list_sequences = list(read_fastq(fastq_file))
    dico = {}
    for seq in list_sequences:
        list_kmer = list(cut_kmer(seq, kmer_size))
        for kmer in list_kmer:
            if kmer not in dico :
                dico[kmer] = list_kmer.count(kmer)
    return dico


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graph = nx.DiGraph()

    for kmer, weight in kmer_dict.items():
        prefix = kmer[:-1]
        suffix = kmer[1:]

        if graph.has_edge(prefix, suffix):
            graph[prefix][suffix]['weight'] += weight
        else:
            graph.add_edge(prefix, suffix, weight=weight)

    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:

        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)

        elif not delete_entry_node and not delete_sink_node:
            graph.remove_nodes_from(path[1:(len(path)-1)])

        elif delete_entry_node:
            graph.remove_nodes_from(path[:(len(path)-1)])

        elif delete_sink_node:
            graph.remove_nodes_from(path[1:])

    return graph


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    if len(set(weight_avg_list)) > 1:
        weight_stdev = statistics.stdev(weight_avg_list)
    else:
        weight_stdev = 0

    if len(set(path_length)) > 1:
        length_stdev = statistics.stdev(path_length)
    else:
        length_stdev = 0

    if weight_stdev > 0:
        max_weight = max(weight_avg_list)
        best_paths = [path_list[i] for i,
                      weight in enumerate(weight_avg_list) if weight == max_weight]
    elif length_stdev > 0:
        max_length = max(path_length)
        best_paths = [path_list[i] for i, length in enumerate(path_length) if length == max_length]
    else:
        random_index = randint(0, len(path_list) - 1)
        best_paths = [path_list[random_index]]


    paths_to_remove = [path for path in path_list if path not in best_paths]

    graph = remove_paths(graph, paths_to_remove, delete_entry_node, delete_sink_node)

    return graph


def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    all_paths = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))
    path_lengths = [len(path) for path in all_paths]
    weight_avg_list = [path_average_weight(graph, path) for path in all_paths]
    cleaned_graph = select_best_path(graph, all_paths, path_lengths, weight_avg_list)

    return cleaned_graph


def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False

    for node_n in graph:
        predecessors_list = list(graph.predecessors(node_n))

        if len(predecessors_list) > 1:
            for i, pred_i in enumerate(predecessors_list[:-1]):
                for pred_j in predecessors_list[i + 1:]:
                    ancestor_node = nx.lowest_common_ancestor(graph, pred_i, pred_j)

                    if ancestor_node is not None:
                        bubble = True
                        break

                if bubble:
                    break

            if bubble:
                break

    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, ancestor_node, node_n))

    return graph


def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    list_nodes = list(graph.nodes())
    path_list = []

    for node in list_nodes:
        predecessors = list(graph.predecessors(node))

        if len(predecessors) > 1:
            for start_node in starting_nodes:
                paths = list(nx.all_simple_paths(graph, start_node, node))
                if paths:
                    path_list.append(list(paths)[0])

            if len(path_list) > 1:
                len_list = [len(path) for path in path_list]
                path_weights = []
                for i, path in enumerate(path_list):
                    if len_list[i] > 1:
                        weight = path_average_weight(graph, path)
                    else:
                        weight = graph[path[0]][path[1]]["weight"]
                    path_weights.append(weight)

                new_graph = select_best_path(
                    graph, path_list, len_list, path_weights,
                    delete_entry_node=True, delete_sink_node=False
                )

                if new_graph != graph:
                    return solve_entry_tips(new_graph, starting_nodes)
                else:
                    return graph
    return graph


def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    list_nodes = list(graph.nodes())
    path_list = []

    for node in list_nodes:
        successors = list(graph.successors(node))

        if len(successors) > 1:
            for ending_node in ending_nodes:
                paths = list(nx.all_simple_paths(graph, node, ending_node))
                if paths:
                    path_list.append(list(paths)[0])

            if len(path_list) > 1:
                len_list = [len(path) for path in path_list]
                path_weights = []
                for i, path in enumerate(path_list):
                    if len_list[i] > 1:
                        weight = path_average_weight(graph, path)
                    else:
                        weight = graph[path[0]][path[1]]["weight"]
                    path_weights.append(weight)

                new_graph = select_best_path(
                    graph, path_list, len_list, path_weights,
                    delete_entry_node=False, delete_sink_node=True
                )

                if new_graph != graph:
                    return solve_out_tips(new_graph, ending_nodes)
                else:
                    return graph

    return graph


def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    list_graph = list(graph.nodes)
    list_node = []
    for node in list_graph:
        list_pred = list(graph.predecessors(node))
        if len(list_pred) == 0:
            list_node.append(node)
    return list_node


def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    list_graph = list(graph.nodes)
    list_node = []
    for node in list_graph:
        list_pred = list(graph.successors(node))
        if len(list_pred) == 0:
            list_node.append(node)

    return list_node


def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    list_contigs = []
    for node_start in starting_nodes:
        for node_end in ending_nodes:
            if nx.has_path(graph, node_start, node_end):
                paths = nx.all_simple_paths(graph, node_start, node_end)
                for path in paths:
                    sequence = path[0]
                    for node in path[1:]:
                        sequence += node[-1]
                    list_contigs.append([sequence, len(sequence)])

    return list_contigs


def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, 'w') as file:
        for index, (contig, length) in enumerate(contigs_list):
            header = f">contig_{index} len={length}"
            sequence = textwrap.fill(contig, width=80)
            file.write(header + "\n")
            file.write(sequence + "\n")


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """
    #fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    graph = simplify_bubbles(graph)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    new_starting_nodes = get_starting_nodes(graph)
    new_ending_nodes = get_sink_nodes(graph)
    contigs_list = get_contigs(graph, new_starting_nodes, new_ending_nodes)
    save_contigs(contigs_list, args.output_file)


    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()
    