#!/usr/bin/python3

import argparse

import gzip
from Bio.Data.SCOPData import protein_letters_3to1
import Bio.PDB

import numpy as np
import math

import sys, os
import re

parser = argparse.ArgumentParser(description='Create a graph database file from a group of .sif files')
parser.add_argument('-g', '--grouping', type=int, choices=[0, 1, 2, 3], default=0, help='residue type grouping (1 = least specific grouping, ..., 3 = most specific grouping, 0 = no grouping)')
parser.add_argument('-t', '--disable-type', action="store_true", help='disable use of interaction type labels')
parser.add_argument('-d', '--disable-distance', action="store_true", help='disable use of distance labels')
parser.add_argument('-s', '--enable-secondary', action="store_true", help='enable use of secondary structure labels')
parser.add_argument('-a', '--enable-accessibility', action="store_true", help='enable use of surface accessibility labels')
parser.add_argument('-b', '--enable-angles', action="store_true", help='disable use of ca-cb angle labels')
parser.add_argument('-x', '--suffix', type=str, default=None, help="read segment ranges from PDB_suffix.seg instead of PDB.seg")
parser.add_argument('-c', '--classification', action="store_true", help='generate database for subgraph based classification significance calculation')
parser.add_argument('out_file', metavar='database', help='name of graph database output file')
parser.add_argument('in_files', metavar='input.sif', nargs='+', help='list of .sif files')
args = parser.parse_args()

# Single letter codes sorted according to their three letter codes
amino_acid_order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

class FixedStructureBuilder(Bio.PDB.StructureBuilder.StructureBuilder):
    # Allcaps for elements. This has been properly fixed in Biopython many years ago, but that version may not be available on all machines...
    def init_atom (self, name, coord, b_factor, occupancy, altloc, fullname, serial_number=None, element=None):
        if element != None:
            element = element.upper()
        Bio.PDB.StructureBuilder.StructureBuilder.init_atom(self, name, coord, b_factor, occupancy, altloc, fullname, serial_number, element)

def create_group (tmp_grps):
    result = {}
    for i, grp in enumerate(tmp_grps):
        for res_type in grp:
            result[res_type] = i
    return result

# mappings from amino acid to numerical label; line graphs only support these

# group similar amino acids (broad)
type_groups1 = create_group([['G', 'A', 'S', 'T'], ['C', 'V', 'I', 'L', 'P', 'F', 'Y', 'M', 'W'], ['N', 'Q', 'H'], ['D', 'E'], ['K', 'R']])

# group similar amino acids (fine)
type_groups2 = create_group([['G', 'A'], ['S', 'T'], ['C'], ['V', 'I', 'L'], ['P'], ['F', 'Y', 'W'], ['M'], ['N', 'Q'], ['H'], ['D', 'E'], ['K', 'R']])

# group similar amino acids (finer)
type_groups3 = create_group([['G'], ['A'], ['S', 'T'], ['C'], ['V', 'I', 'L'], ['P'], ['F', 'Y'], ['W'], ['M'], ['N'], ['Q'], ['H'], ['D', 'E'], ['K', 'R']])

# each type gets its own label, no grouping (finest)
individual_types = {y: x for (x, y) in enumerate(amino_acid_order)}

if args.grouping == 0:
    grouping_scheme = individual_types
elif args.grouping == 1:
    grouping_scheme = type_groups1
elif args.grouping == 2:
    grouping_scheme = type_groups2
elif args.grouping == 3:
    grouping_scheme = type_groups3


include_distances = not args.disable_distance
include_types = not args.disable_type
include_secondary_structure = args.enable_secondary
include_solvent_accessibility = args.enable_accessibility
include_angles = args.enable_angles

output_file_path = args.out_file

if output_file_path.endswith('.sif'):
    sys.stderr.write('Please specify a output file\n')
    exit(1)

output_dir = os.path.dirname(output_file_path)

output_message_path = output_dir + '/graph_gen_output.txt'
output_message_file = open(output_message_path, 'w')
output_error_path = output_dir + '/graph_gen_errors.txt'
output_error_file = open(output_error_path, 'w')

def print_message (msg):
    sys.stdout.write(msg)
    sys.stdout.flush()
    output_message_file.write(msg)

def print_error (msg):
    sys.stderr.write(msg)
    sys.stderr.flush()
    output_error_file.write(msg)

if not (include_distances or include_types):
    print_error('No edge label format specified\n')
    exit(1)

print_message('Suffix: {}\n'.format(args.suffix))
print_message('Grouping: {}\n'.format(args.grouping))
print_message('Include interaction type labels: {}\n'.format(include_types))
print_message('Include secondary structure labels: {}\n'.format(include_secondary_structure))
print_message('Include surface accessibility labels: {}\n'.format(include_solvent_accessibility))
print_message('Include distance labels: {}\n'.format(include_distances))
print_message('Include angle labels: {}\n'.format(include_angles))

print_message('\n')

####################################

def three_to_one (aa):
    error = None

    try:
        result = protein_letters_3to1[aa]
        if result == 'X':
            error == '\'any\' type residue'
    except:
        error = 'Invalid three-letter code \'{}\''.format(aa)
        result = 'X'

    # Residues of type X are ignored
    if len(result) != 1:
        error = 'Invalid three-letter code mapping \'{}\'->\'{}\''.format(aa, result)
        result = 'X'

    if result == 'Z':
        # Wikipedia:
        # Glutamic acid (E) or glutamine (Q)
        # A placeholder when either amino acid may occupy a position
        result = 'E'

    return result, error

# Used to convert 'A:123B' to ('A', 123, 'B')
def string_to_res(res_text):
    try:
        convert_regex = re.compile('([a-zA-Z0-9]+):(-?[0-9]+)([a-zA-Z]*)')
        groups = convert_regex.match(res_text).groups()
    except:
        print('Error when converting residue string: {}'.format(res_text))
        raise
    return (groups[0], int(groups[1]), groups[2])

def res_to_string(res):
    return '{}:{}{}'.format(res[0], res[1], res[2])

def parse_sif_res_string(res_string):
    orig_res_chain, orig_res_id, orig_icode, res_type = res_string.split(':')
    orig_res_id = int(orig_res_id)
    if orig_icode == '_':
        orig_icode = ''
    res_id = (orig_res_chain, orig_res_id, orig_icode)
    res = (res_id, res_type)
    return res

def load_segments (path):
    if not os.path.exists(path):
        print_error('File not found: {0}\n'.format(path))
        return None

    f = open(path, 'r')

    segments = {}

    for line in f:
        if line == "\n":
            continue

        _, chain, frm, to = map(str.strip, line.split(','))
        if chain not in segments:
            segments[chain] = []

        frm_icode = None
        if frm != "_":
            frm, frm_icode = re.match("(-?\d+)(\D*)", frm).groups()
            frm = int(frm)

        to_icode = None
        if to != "_":
            to, to_icode = re.match("(-?\d+)(\D*)", to).groups()
            to = int(to)

        segments[chain].append((frm, frm_icode, to, to_icode))

    f.close()

    return segments

def is_in_segments(segments, res):
    chain, res_id, icode = res
    if chain not in segments:
        return False
    for range_ in segments[chain]:
        frm, frm_icode, to, to_icode = range_
        if (frm == "_" or res_id > frm or (res_id == frm and icode >= frm_icode)) \
           and (to == "_" or res_id < to or (res_id == to and icode <= to_icode)):
            return True
    return False

def run_dssp (structure, path):
    ss_type_map = {}
    acc_map = {}

    model = structure[0]

    try:
        dssp = Bio.PDB.DSSP(model, path)
    except:
        print_error('Error while running dssp on {}:\n'.format(path))
        raise

    for biopy_chain, biopy_res_id in dssp.keys():
        values = dssp[(biopy_chain, biopy_res_id)]

        icode = biopy_res_id[2]
        if icode == ' ':
            icode = ''
        res = (biopy_chain, biopy_res_id[1], icode)

        # 0: Loop
        # 1: Helix
        # 2: Sheet
        ss_type = 0

        ss_type_string = values[2]

        # Helices
        if ss_type_string in ['H', 'G', 'I']:
            ss_type = 1
        # Beta-Strand
        # consider 'B' a loop, they are too isolated
        elif (ss_type_string == 'E'):
            ss_type = 2

        # 0: accessible
        # 1: inaccessible
        acc_type = 0

        rel_sol_surf_acc_area = values[3]

        if rel_sol_surf_acc_area < 0.16:
            acc_type = 1

        ss_type_map[res] = ss_type
        acc_map[res] = acc_type

    return ss_type_map, acc_map

def gen_biopy_res_map(structure):
    res_map = {}
    try:
        for chain in structure[0]:
            chain_id = chain.id
            for res in chain.get_residues():
                het, res_id, icode = res.id
                if icode == ' ':
                    icode = ''
                res_map[(chain_id, res_id, icode)] = res.id
    except:
        print_error('Could not generate biopython residue map for: {}\n'.format(structure))
        raise
    return res_map

def get_biopy_residue(structure, res, res_map):
    chain, res_id, icode = res
    biopy_res_id = res_map[res]
    biopy_res = structure[0][chain][biopy_res_id]
    return biopy_res

def ca_is_in_structure(biopy_res):
    result = 'CA' in biopy_res
    return result

def get_ca_coords(biopy_res):
    coords = biopy_res['CA'].get_coord()
    return coords

def cb_is_in_structure(biopy_res):
    result = 'CB' in biopy_res or \
             ('N' in biopy_res and 'C' in biopy_res and 'CA' in biopy_res)
    return result

def get_cb_coords(biopy_res):
    try:
        cb = biopy_res['CB'].get_coord()
    except KeyError:
        # From: http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
        # get atom coordinates as vectors
        n = biopy_res['N'].get_vector()
        c = biopy_res['C'].get_vector()
        ca = biopy_res['CA'].get_vector()
        # center at origin
        n = n - ca
        c = c - ca
        # find rotation matrix that rotates n -120 degrees along the ca-c vector
        rot = Bio.PDB.rotaxis(-math.pi*120.0/180.0, c)
        # apply rotation to ca-n vector
        cb_at_origin = n.left_multiply(rot)
        # put on top of ca atom
        cb = cb_at_origin + ca
        cb = cb.get_array()
    return cb

def convert_file_to_graph (sif_path, suffix):
    if not os.path.exists(path):
        print_error('File not found: {0}\n'.format(sif_path))
        exit(1)

    if not path.endswith('.sif'):
        print_error('File not in sif format : {0}\n'.format(sif_path))
        exit(1)

    if suffix:
        seg_path = sif_path.replace(".sif", "_" + suffix + ".seg")
    else:
        seg_path = sif_path.replace(".sif", ".seg")
    segments = load_segments(seg_path)

    fam = os.path.basename(os.path.dirname(os.path.dirname(sif_path)))
    fams = set([fam])

    pdb_id = os.path.basename(sif_path)[0:-4]
    pdb_parser = Bio.PDB.PDBParser(QUIET = True, structure_builder=FixedStructureBuilder())
    pdb_path = None
    pdb_gzipped = False
    for suffix in [".pdb", ".pdb.gz", ".ent", ".ent.gz"]:
        test_path = sif_path.replace(".sif", suffix)
        if os.path.isfile(test_path):
            pdb_path = test_path
            if pdb_path.endswith(".gz"):
                pdb_gzipped = True
            break
    if pdb_path == None:
        print_error('No associated PDB file found: {0}\n'.format(sif_path))
        exit(1)

    if pdb_gzipped:
        pdb_file = gzip.open(pdb_path, 'rt')
    else:
        pdb_file = open(pdb_path, 'rt')

    structure = pdb_parser.get_structure(pdb_id, pdb_file)
    pdb_file.close()

    res_map = gen_biopy_res_map(structure)

    residue_types = {}
    interactions = set([])

    # for recoverable errors
    problems = []

    f = open(sif_path, 'r')

    for line in f:
        temp = line.split()
        #edges
        if len(temp) == 3:
            (res_string_A, edge, res_string_B) = temp

            # ignore edge type
            if not edge == "combi:all_all":
                continue

            res_A, res_type_A = parse_sif_res_string(res_string_A)
            res_B, res_type_B = parse_sif_res_string(res_string_B)

            # allow using just a subset of residues based on segments files
            if not is_in_segments(segments, res_A):
                continue
            if not is_in_segments(segments, res_B):
                continue

            residue_types[res_A] = res_type_A
            residue_types[res_B] = res_type_B

            # additional safety measure against accidental double edges
            if res_A < res_B:
                interactions.add((res_A, res_B, 1))
            else:
                interactions.add((res_B, res_A, 1))

        #nodes
        elif len(temp) == 1:
            res_string = temp[0]
            res, res_type = parse_sif_res_string(res_string)

            if not is_in_segments(segments, res):
                continue

            residue_types[res] = res_type

        else:
            print_error('Unhandled input format\n')
            exit(1)

    f.close()

    # remove entries with missing C-alpha/beta atoms or 'any' (X) type residues
    removed_residues = []
    for res, res_type in sorted(residue_types.items()):
        biopy_res = get_biopy_residue(structure, res, res_map)

        if include_distances and not ca_is_in_structure (biopy_res):
            problems.append('Residue with missing C-alpha atom removed ({})\n'.format(res_to_string(res)))
            removed_residues.append(res)
            continue

        if include_angles and not cb_is_in_structure (biopy_res):
            problems.append('Residue with missing C-beta atom (or atoms required to approximate it) removed ({})\n'.format(res_to_string(res)))
            removed_residues.append(res)
            continue

        aa_code, error = three_to_one(res_type)
        if aa_code == 'X':
            problems.append('Residue removed because of type issue: {} ({})\n'.format(error, res_to_string(res)))
            removed_residues.append(res)
            continue

        if res[2] != '':
            problems.append('Residue with insertion code found ({})\n'.format(res_to_string(res)))
            # Not removed

    for res in removed_residues:
        del residue_types[res]

    # remove interactions containing removed residues
    interactions = [i for i in interactions if not i[0] in removed_residues and not i[1] in removed_residues]

    # the residues list need to be sorted by residue id for the next step
    residues = sorted(residue_types.keys())

    # Add covalent bonds between consecutive residues
    for i in range(len(residues) - 1):
        res_A = residues[i]
        res_B = residues[i + 1]
        chain_A, res_id_A, icode_A = res_A
        chain_B, res_id_B, icode_B = res_B

        # Determine if residues are consecutive
        # Some insert codes start at A, some start at B
        if chain_A == chain_B \
           and (res_id_A + 1 == res_id_B
                or (res_id_A == -1 and res_id_B == 1)
                or (res_id_A == res_id_B and icode_A == '' and icode_B == 'A')
                or (res_id_A == res_id_B and icode_A == '' and icode_B == 'B')
                or (res_id_A == res_id_B and icode_A != '' and chr(ord(icode_A) + 1) == icode_B)):

            if (res_A, res_B, 1) in interactions:
                interactions.remove((res_A, res_B, 1))

            interactions.append((res_A, res_B, 0))
        else:
            problems.append('Non-consecutive residue IDs ({}, {})\n'.format(res_to_string(res_A), res_to_string(res_B)))

    # map residues to a numerical id starting at 0 since gaston expects that
    res_to_vid = {}
    vid_to_res = {}
    for i, res in enumerate(residues):
        res_to_vid[res] = i
        vid_to_res[i] = res


    if include_secondary_structure or include_solvent_accessibility:
        ss_map, acc_map = run_dssp(structure, pdb_path)

    # prepare the final vertex list
    vertices = []
    for res in residues:
        res_type = residue_types[res]
        aa_code, error = three_to_one(res_type)
        aa_label = grouping_scheme[aa_code]

        ss_label = None
        if include_secondary_structure:
            try:
                ss_label = ss_map[res]
            except KeyError:
                ss_label = 0
                problems.append('Missing secondary structure information, falling back to turn ({})\n'.format(res_to_string(res)))

        acc_label = None
        if include_solvent_accessibility:
            try:
                acc_label = acc_map[res]
            except KeyError:
                acc_label = 0
                problems.append('Missing surface accessibility information, falling back to accessible ({})\n'.format(res_to_string(res)))

        vertices.append((res_to_vid[res], aa_label, ss_label, acc_label, res, residue_types[res]))

    # calculate distances for each edge
    ca_coords = {}
    cb_coords = {}
    for res in residues:
        biopy_res = get_biopy_residue(structure, res, res_map)

        if include_distances or include_angles:
            ca_coords[res] = get_ca_coords(biopy_res)
        else:
            ca_coords[res] = None

        if include_angles:
            cb_coords[res] = get_cb_coords(biopy_res)
        else:
            cb_coords[res] = None

    edges = []
    for res_A, res_B, interaction_type in sorted(interactions):
        if include_distances:
            distance = np.linalg.norm(ca_coords[res_A] - ca_coords[res_B])
            distance = int(distance * 1000)
        else:
            distance = 0

        if include_angles:
            res_A_vec_ab = cb_coords[res_A] - ca_coords[res_A]
            res_A_vec_ab = res_A_vec_ab / np.linalg.norm(res_A_vec_ab)
            res_B_vec_ab = cb_coords[res_B] - ca_coords[res_B]
            res_B_vec_ab = res_B_vec_ab / np.linalg.norm(res_B_vec_ab)
            angle = np.arccos (np.dot(res_A_vec_ab, res_B_vec_ab))
            angle = int(angle * 1000)
        else:
            angle = 0


        edges.append((res_to_vid[res_A], res_to_vid[res_B], interaction_type, distance, angle))


    result = (vertices, edges, pdb_id, fams, sif_path, pdb_path, seg_path)

    if len(vertices) == 0:
        result = None
        problems.append('Skipped: Empty graph\n')

    if problems:
        print_error('Potential problems when processing {}:\n'.format(sif_path))
        for problem in problems:
            print_error(problem)
        print_error('\n')

    return result

graphs_for_pdb_id = {}

for path in args.in_files:
    if not path.endswith('.sif'):
        print_error('Input file doesn\'t end in .sif: {0}\n'.format(path))
        exit(1)

    g = convert_file_to_graph(os.path.abspath(path), args.suffix)

    if g:
        pdb_id = g[2]
        try:
            graphs_for_pdb_id[pdb_id].append(g)
        except:
            graphs_for_pdb_id[pdb_id] = []
            graphs_for_pdb_id[pdb_id].append(g)


# combine families of duplicate entries with matching residue sets
all_graphs = []
print_new_line = False
for pdb_id, graphs in sorted(graphs_for_pdb_id.items()):
    added_graphs = []
    merged_families = []
    for graph in graphs:
        vertices = graph[0]
        residue_set = set()
        for vid, aa_label, ss_label, acc_label, res, residue_type in vertices:
            residue_set.add(res)
        new_graph = True
        for i in range(len(added_graphs)):
            existing_residue_set, existing_graph = added_graphs[i]
            if residue_set.issubset(existing_residue_set):
                # add new family to existing graph and don't add new graph
                # this covers the case for identical residue sets as well
                existing_graph[3].update(graph[3])
                merged_families[i].update(graph[3])
                new_graph = False
            elif existing_residue_set.issubset(residue_set):
                # replace existing graph with new one and copy old families
                graph[3].update(existing_graph[3])
                added_graphs[i] = (residue_set, graph)
                merged_families[i].update(graph[3])
                new_graph = False
        if new_graph:
            added_graphs.append((residue_set, graph))
            merged_families.append(graph[3])
    for m in merged_families:
        if len(m) > 1:
            print_error('Overlapping entries for {} in families {} merged\n'.format(pdb_id, ', '.join(sorted(m))))
            print_new_line = True
    all_graphs.extend([g[1] for g in added_graphs])
if print_new_line:
    print_error('\n')

#gaston
f = open(output_file_path, 'w')
output_dir = os.path.abspath(os.path.dirname(output_file_path))
for (i, g) in enumerate(all_graphs):
    f.write('t # {}\n'.format(i))
    (vertices, edges, pdb_id, fams, sif_path, pdb_path, seg_path) = g

    # sif file path
    sif_path = os.path.relpath(sif_path, output_dir)
    f.write('# sif {}\n'.format(sif_path))

    # pdb file path
    pdb_path = os.path.relpath(pdb_path, output_dir)
    f.write('# pdb {}\n'.format(pdb_path))

    # segments file path
    seg_path = os.path.relpath(seg_path, output_dir)
    f.write('# seg {}\n'.format(seg_path))

    # pdb id
    f.write('# id {}\n'.format(pdb_id))

    # family
    if args.classification:
        f.write('# fam {}\n'.format(':'.join(fams)))

    # residue mapping
    f.write('# r')
    for vid, aa_label, ss_label, acc_label, res, residue_type in vertices:
        f.write(' {}'.format(res_to_string(res)))
    f.write('\n')

    # nodes
    for vid, aa_label, ss_label, acc_label, res, residue_type in vertices:
        labels = [aa_label]
        if include_secondary_structure:
            labels.append(ss_label)
        if include_solvent_accessibility:
            labels.append(acc_label)
        f.write('v {} {}\n'.format(vid, ' '.join(map(str, labels))))

    # edges
    for node_id_A, node_id_B, edge_label, distance_label, angle_label in edges:
        labels = []
        if include_types:
            labels.append(edge_label)
        if include_distances:
            labels.append(distance_label)
        if include_angles:
            labels.append(angle_label)

        f.write('e {} {} {}\n'.format(node_id_A, node_id_B, ' '.join(map(str, labels))))
            
f.close()


#statistics
amino_acids = sorted(amino_acid_order)

aa_count = {}
for aa in amino_acids:
    aa_count[aa] = 0

total_residues = 0

for g in all_graphs:
    vertices = g[0]
    for vid, aa_label, ss_label, acc_label, res, residue_type in vertices:
        short_res_type, error = three_to_one(residue_type)
        aa_count[short_res_type] += 1
        total_residues += 1

print_message("Absolute:\n")
for aa in amino_acids:
    print_message(aa + ": " + str(aa_count[aa]) + "\n")
print_message("Total: " + str(total_residues) + "\n")
print_message("\n")

print_message("Relative:\n")
for aa in amino_acids:
    print_message(aa + ": " + str(aa_count[aa]/float(total_residues)) + "\n")

output_message_file.close()
output_error_file.close()
