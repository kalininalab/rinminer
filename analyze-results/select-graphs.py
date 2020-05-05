#!/usr/bin/python3

import argparse
import sys, os
import igraph

class IGraph:
    def __init__(self, graph):
        self.igraph = igraph.Graph(n=len(graph.vertices), directed=False)
        self.igraph.add_edges([(e[0], e[1]) for e in graph.edges])
        self.v_labels = [v[1] for v in graph.vertices]
        self.e_labels = [e[2] for e in graph.edges]

    def subg_iso(self, other):
        return self.igraph.subisomorphic_vf2(other.igraph,
                                             color1=self.v_labels,
                                             color2=other.v_labels,
                                             edge_color1=self.e_labels,
                                             edge_color2=other.e_labels)

class Graph:
    def __init__(self, id, offset):
        self.id = id
        self.offset = offset
        self.vertices = []
        self.edges = []
        self.parents = []
        self.size = 0
        self.score = 0
        self.igraph = None

class GraphsCollection:
    def __init__(self, graphs_file_path):
        self.id_to_graph = {}
        self.database = []

        self._parse_graphs(graphs_file_path)

    def _add_graph(self, g):
        self.id_to_graph[g.id] = g

    def _parse_graphs (self, graphs_file_path):
        graph = None
        header = True
        last_id = -1

        with open(graphs_file_path, 'r', encoding='ascii') as graphs_file:
            offset = graphs_file.tell()
            line = graphs_file.readline()
            while line:
                if line[0] == 'v':
                    splt = [int(x) for x in line[2:].split()]
                    id_, aa_label = splt[0], splt[1]
                    graph.vertices.append((id_, aa_label))

                elif line[0] == 'e':
                    splt = [int(x) for x in line[2:].split()]
                    id_A, id_B, interaction_label = splt[0], splt[1], splt[2]
                    graph.edges.append((id_A, id_B, interaction_label))
                    graph.size += 1

                elif line[0] == 't':
                    header = False

                    if graph is not None:
                        self._add_graph(graph)

                    try:
                        graph_id = int(line.split()[1])
                    except:
                        graph_id = last_id
                        last_id += 1

                    graph = Graph(graph_id, offset)
                    graph.parent_ids = []

                elif line[0] == '#':
                    if not header:
                        splt = line.strip().split(' ')
                        if splt[1] == 'p':
                            graph.parent_ids = [int(x) for x in splt[2:]]
                        elif splt[1] == 's':
                            graph.score = int(splt[2])
                        elif splt[1] == 'i':
                            graph.stats = [int(x) for x in splt[2:]]

                elif line == '\n':
                    pass

                else:
                    sys.stderr.write('Unhandled input format: {0}\n'.format(line))

                offset = graphs_file.tell()
                line = graphs_file.readline()

            if graph is not None:
                self._add_graph(graph)

        # Establish parent relationships between graphs
        for g_id, g in self.id_to_graph.items():
            for p_id in g.parent_ids:
                try:
                    p = self.id_to_graph[p_id]
                    g.parents.append(p)
                except:
                    pass
            del g.parent_ids

    def _parents(self, graph, degree):
        parents = set([graph])
        for depth in range(degree):
            new_parents = set()
            for parent in parents:
                new_parents.update(parent.parents)
            # We reached the minimum included graph size
            if len(new_parents) == 0:
                return parents
            parents = new_parents

        return parents
        
    def get_offsets(self):
        results = {}
        for g_id, graph in self.id_to_graph.items():
            results[g_id] = graph.offset
        return results

    def select_graphs(self, num_graphs, min_size, max_size):
        results = []

        for graph in sorted(self.id_to_graph.values(), key=lambda x: x.score, reverse=True):
            include = True

            if graph.size < min_size or graph.size > max_size:
                continue
            
            graph.igraph = IGraph(graph)

            for included in results:
                if len(included.vertices) >= len(graph.vertices) \
                and included.igraph.subg_iso(graph.igraph):
                    include = False
                    break

            if include:
                results.append(graph)
                if len(results) == num_graphs:
                    break
            else:
                del graph. igraph

        return results

parser = argparse.ArgumentParser(description='Select a set of the highest scoring motifs while excluding too similar motifs')
parser.add_argument('-b', '--batch', action="store_true", help='no interactive input required, only requested output')
parser.add_argument('-f', '--force', action="store_true", help='force clearing the directory first')
parser.add_argument('-m', '--min-size', type=int, default=1, help='only include subgraphs at least this size')
parser.add_argument('-x', '--max-size', type=int, default=sys.maxsize, help='only include subgraphs at most this size')
parser.add_argument('-l', '--limit', type=int, default=15, help='limit output to LIMIT best graphs')

parser.add_argument('in_file', metavar='result.lg', help='input file')
parser.add_argument('out_index_file', metavar='result_index.txt', help='output graph index file')
parser.add_argument('out_file', metavar='output.lst', help='output file')
args = parser.parse_args()

try:
    if os.path.exists(args.out_index_file):
        if not os.path.isfile(args.out_index_file):
            print("Graph index file ({}) not a file, exiting.".format(args.out_index_file))
            sys.exit(0)
        if args.force:
            if not args.batch:
                print("Graph index file ({}) exists, overwriting.".format(args.out_index_file))
        elif not args.batch:
            sys.stdout.write("Graph index file ({}) exists. Overwrite [Y/n]? ".format(args.out_index_file))
            sys.stdout.flush()
            inp = sys.stdin.readline().strip()
            if not inp.upper() == 'Y' and not inp == '':
                sys.exit(0)
        else:
            sys.exit(0)

    if os.path.exists(args.out_file):
        if not os.path.isfile(args.out_file):
            print("Output file ({}) not a file, exiting.".format(args.out_file))
            sys.exit(0)
        if args.force:
            if not args.batch:
                print("Output file ({}) exists, overwriting.".format(args.out_file))
        elif not args.batch:
            sys.stdout.write("Output file ({}) exists. Overwrite [Y/n]? ".format(args.out_file))
            sys.stdout.flush()
            inp = sys.stdin.readline().strip()
            if not inp.upper() == 'Y' and not inp == '':
                sys.exit(0)
        else:
            sys.exit(0)

    collection = GraphsCollection(args.in_file)
    graphs = collection.select_graphs(args.limit, args.min_size, args.max_size)
    offsets = collection.get_offsets()
    with open(args.out_index_file, 'w', encoding='ascii') as out:
        for g_id, offset in sorted(offsets.items()):
            out.write('{} {}\n'.format(g_id, offset))
    with open(args.out_file, 'w', encoding='ascii') as out:
        for g in [ str(g.id) for g in graphs ]:
            out.write('{}\n'.format(g))

except FileNotFoundError:
    if args.batch:
        pass
    else:
        raise