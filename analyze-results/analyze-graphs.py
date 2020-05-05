#!/usr/bin/python3

import argparse
import subprocess
import sys, os
import regex as re
import math as m
import random
from scipy.stats import fisher_exact
from Bio.Data.SCOPData import protein_letters_3to1
import Bio.PDB
import gzip
from itertools import combinations

random.seed('THIS NEEDS TO BE STATIC ACROSS DIFFERENT ITERATIONS')

class FixedStructureBuilder(Bio.PDB.StructureBuilder.StructureBuilder):
    # Ignore hetero flag to allow addressing HETATM residues via their ID+icode
    def init_residue(self, resname, field, resseq, icode):
        Bio.PDB.StructureBuilder.StructureBuilder.init_residue(self, resname, ' ', resseq, icode)

    # Allcaps for elements. This works around a bug that has been fixed many
    # years ago. Thanks Debian.
    def init_atom (self, name, coord, b_factor, occupancy, altloc, fullname, serial_number=None, element=None):
        if element != None:
            element = element.upper()
        Bio.PDB.StructureBuilder.StructureBuilder.init_atom(self, name, coord, b_factor, occupancy, altloc, fullname, serial_number, element)

pdb_parser = Bio.PDB.PDBParser(QUIET = True, structure_builder=FixedStructureBuilder())

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
individual_types = {y: x for (x, y) in enumerate(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])}

ss_types = ['L', 'A', 'B']

def gen_reverse_mapping(grouping, ss):
    grouping_scheme_num_types = len(set(grouping_scheme.values()))

    ss_num_types = 1
    if ss:
        ss_num_types = 3

    result = {}

    for x in range(ss_num_types * grouping_scheme_num_types):
        grp_id = x % grouping_scheme_num_types
        ss_type = (x - grp_id) / grouping_scheme_num_types
        res_type = []
        for res, grp in grouping.items():
            if grp == grp_id:
                res_type.append(res)
        res_type = sorted(res_type)
        temp = ""
        if ss:
            temp += "["
        temp += ", ".join(res_type)
        if ss:
            temp += "] (" + ss_types[ss_type] + ")"
        result[x] = temp

    return result

def load_segments (path):
    if not os.path.exists(path):
        sys.stderr.write('File not found: {0}\n'.format(path))
        sys.exit(1)

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

def is_in_segments(segments, chain, res_id, icode):
    if chain not in segments:
        return False

    for range_ in segments[chain]:
        frm, frm_icode, to, to_icode = range_

        if (frm == "_" or res_id > frm or (res_id == frm and icode >= frm_icode)) \
           and (to == "_" or res_id < to or (res_id == to and icode <= to_icode)):
            return True

    return False

def parse_pattern(pattern):
    items = pattern.split('-')
    item_numbers = []
    scores_binary = []
    scores_percent = []
    regex = r''
    for item_number, item in enumerate(items):
        item_score = None

        # These are not really enforceable, because the sequence in the
        # structure might contain additional residues
        item = item.replace('<', '').replace('>', '')

        # Replace cases like [<A] with A after removing < or >
        item = re.sub(r'\[(.)\]', r'\1', item)

        min_size = 1
        max_size = 1
        size = re.match(r'(.*)\((.*)\)', item)
        if size:
            item = size.group(1)
            size = [int(x) for x in size.group(2).split(',')]
            try:
                min_size, max_size = size
            except:
                min_size = size[0]
                max_size = size[0]

        group = re.match(r'\[(.*)\]', item)
        excl_group = re.match(r'\{(.*)\}', item)

        if group:
            # Not very specific, but still somehat important
            regex_item = '[{}]'.format(group.group(1))
            item_score_binary = 0.
            # Residue has to have at least one amino acid type, so divide by 19 instead of 20
            item_score_percent = float(20 - len(group.group(1))) / 19.
        elif excl_group:
            # Probably not too important, can make only a small contribution
            regex_item = '[^{}]'.format(excl_group.group(1))
            item_score_binary = 0.
            # Residue has to have at least one amino acid type, so divide by 19 instead of 20
            item_score_percent = float(len(excl_group.group(1))) / 19.
        else:
            if item == 'x':
                # Anything could match, should not be cosidered for scoring
                regex_item = '.'
                item_score_binary = 0.
                item_score_percent = 0.
            else:
                # Exact match, should always be considered for scoring with the highest weight
                regex_item = item
                item_score_binary = 1.
                item_score_percent = 1.

        required_repeats = min_size
        min_size = 0
        max_size = max_size - required_repeats

        for _ in range(required_repeats):
            regex += r'({})'.format(regex_item)
            scores_binary.append(item_score_binary)
            scores_percent.append(item_score_percent)
            item_numbers.append(item_number)

        for _ in range(max_size):
            regex += r'({}?)'.format(regex_item)
            scores_binary.append(item_score_binary)
            scores_percent.append(item_score_percent)
            item_numbers.append(item_number)

    return (regex, item_numbers, scores_binary, scores_percent)

# Used to convert 'A:123B' to ('A', 123, 'B')
def convert_text_to_res(res_text):
    try:
        convert_regex = re.compile('([a-zA-Z0-9]+):(-?[0-9]+)([a-zA-Z]*)')
        groups = convert_regex.match(res_text).groups()
    except:
        print('Error when converting residue string: {}'.format(res_text))
        raise
    return (groups[0], int(groups[1]), groups[2])

def convert_res_to_text(res):
    return '{}:{}{}'.format(res[0], res[1], res[2])

class Graph:
    def __init__(self, id):
        self.id = id
        self.parents = []
        self.vertices = []
        self.edges = []
        self.size = 0
        self.score = 0
        self.mappings = None
        self.class_counts = None
        self.is_interesting = False
        self.is_significant = False
        self.is_top = False # Across the entire parameter set for the database
        self.is_best = False # Within the parameter set for the database

class DatabaseEntry:
    def __init__(self, pretty_id, pretty_vertices, seg_file_path, pdb_file_path, pattern_data):
        self.pretty_id = pretty_id
        self.pretty_vertices = pretty_vertices
        self.seg_file_path = seg_file_path
        self.pdb_file_path = pdb_file_path
        self.pattern_matches_binary = None
        self.pattern_matches_percent = None
        self._structure = None
        self.warned_res_names = set()

        # This relies on a static seed used across the different iterations of this script
        self.shuffled_vertices = pretty_vertices[:]
        random.shuffle(self.shuffled_vertices)

        if pattern_data:
            self.pattern_matches_binary, self.pattern_matches_percent = self.build_pattern_matches(pattern_data)

    def get_structure(self):
        if self._structure != None:
            return self._structure

        if self.pdb_file_path.endswith('.gz'):
            with gzip.open(self.pdb_file_path, 'rt') as pdb_file:
                self._structure = pdb_parser.get_structure(self.pretty_id, pdb_file)
        else:
            with open(self.pdb_file_path, 'rt') as pdb_file:
                self._structure = pdb_parser.get_structure(self.pretty_id, pdb_file)

        return self._structure

    def build_pattern_matches(self, pattern):
        pattern_regex, item_numbers, scores_binary, scores_percent = pattern
        structure = self.get_structure()
        segments = load_segments(self.seg_file_path)
        sequence = ""
        sequence_to_id = []
        last_valid_res_id = 0

        for res in structure.get_residues():
            chain, (hetatm, res_id, icode) = res.get_full_id()[2:]
            if hetatm == " " \
               and is_in_segments(segments, chain, res_id, icode):
                try:
                    sequence += protein_letters_3to1[res.get_resname()]
                    sequence_to_id.append((chain, res_id, icode.strip()))
                    last_valid_res_id = res_id
                except:
                    if res_id == last_valid_res_id + 1 and res.get_resname() not in self.warned_res_names:
                        print('Unknown residue type "{}", skipping'.format(res.get_resname()))
                        self.warned_res_names.add(res.get_resname())

        pattern_matches_binary = []
        pattern_matches_percent = []
        for match in re.finditer(pattern_regex, sequence, overlapped=True):
            match_scores_binary = {}
            match_scores_percent = {}
            for g in range(0, len(item_numbers)):
                item = item_numbers[g]
                score_binary = scores_binary[g]
                score_percent = scores_percent[g]
                span = match.span(g + 1)
                for x in range(span[0], span[1]):
                    res = sequence_to_id[x]
                    match_scores_binary[res] = (score_binary, item)
                    match_scores_percent[res] = (score_percent, item)
            pattern_matches_binary.append(match_scores_binary)
            pattern_matches_percent.append(match_scores_percent)

        return pattern_matches_binary, pattern_matches_percent


class GraphsCollection:
    def __init__(self, graphs_file_path, offsets, graph_ids, reverse_mapping, pattern, pattern_data, structman_classes, base_dir, score_difference, significance_threshold, top):
        self.database = []
        self.graphs = []
        self.id_to_graph = {}
        self.reverse_mapping = reverse_mapping
        self.pattern = pattern
        self.pattern_data = pattern_data
        self.score_difference = score_difference
        self.structman_classes = structman_classes

        if significance_threshold == None:
            self.use_significance = False
            self.significance_threshold = 0.0
        else:
            self.use_significance = True
            self.significance_threshold = significance_threshold

        # candidates
        self._parse_graphs(graphs_file_path, offsets, graph_ids, base_dir)

        if len(self.graphs) > 0:
            self.graphs[0].is_top = top
            self.graphs[0].is_best = True

    def _alignment_selectors(self, pdb_id, residues):
        residues_selectors = []
        for chain, res_id, icode in residues:
            # Use "first" instead of "CA" so pair_fit won't complain about
            # mismatching numbers of atoms if CA is missing or occupancy != 1.0
            residues_selector = "first {0}//{1}/{2}{3}/".format(pdb_id, chain, res_id, icode)
            residues_selectors.append(residues_selector)
        return residues_selectors

    def _residue_selector(self, pdb_id, residues):
        residues_selectors = []
        current_chain = None
        for chain, res_id, icode in residues:
            if chain != current_chain:
                residues_selector = "{0} and chain {1} and resi {2}{3}".format(pdb_id, chain, res_id, icode)
                residues_selectors.append(residues_selector)
                current_chain = chain
                continue
            residues_selector = residues_selector + "+{}{}".format(res_id, icode)
            residues_selectors[-1] = residues_selector
        residues_selector_string = " or ".join(residues_selectors)
        return residues_selector_string

    def generate_pymol_script (self, f, graph):
        object_ids = {}
        sorted_db_ids = []

        first_structure = None
        target_selectors = None

        # sorted by score contribution; target = exact map = original structure
        for db_id, _ in sorted(graph.mappings.items(), key=lambda x: x[1][0][1], reverse=True):
            pretty_id = self.database[db_id].pretty_id
            object_id = pretty_id
            counter = 1
            while object_id in object_ids.values():
                counter += 1
                object_id = '{}_{}'.format(pretty_id, counter)
            object_ids[db_id] = object_id
            sorted_db_ids.append(db_id)

            f.write('fetch {}, name={}, async=0\n'.format(pretty_id, object_id))

            selectors = self._alignment_selectors(object_id, graph.mappings[db_id][0][3])

            if not first_structure:
                first_structure = object_id
                target_selectors = selectors
            else:
                atoms = ', '.join(['{}, {}'.format(mobile, target) for mobile, target in zip(selectors, target_selectors)])
                f.write('pair_fit {}\n'.format(atoms))

        all_mappings = []
        all_first_mappings = []

        for db_id in sorted_db_ids:
            object_id = object_ids[db_id]
            number = 1

            for (exact_, score_contribution_, mapping_, pretty_mapping_, random_mapping_) in graph.mappings[db_id]:
                pretty_mapping = pretty_mapping_
                selection = self._residue_selector(object_id, pretty_mapping)
                name = "{}_mapping_{}".format(object_id, number)
                all_mappings.append(name)
                if number == 1:
                    all_first_mappings.append(name)


                # q  occupancy
                f.write('create {}, ({}) and q > 0.499, 0, 1\n'.format(name, selection))
                # copy of the entire structure as first state to ensure proper alignment in the sequence view
                f.write('create {}, ({}) and q > 0.499, 0, 2\n'.format(name, object_id))

                number += 1

        # Don't include hydrogens in motif mappings
        f.write('remove ({}) and symbol H\n'.format(' or '.join(all_mappings)))

        f.write('hide everything, all\n')
        f.write('show cartoon, {}\n'.format(' or '.join(object_ids.values())))
        f.write('show sticks, {}\n'.format(' or '.join(all_mappings)))
        f.write('disable all\n')
        f.write('enable {}\n'.format(' or '.join(all_first_mappings)))
        f.write('enable {}\n'.format(first_structure))

        f.write('set cartoon_transparency, 0.5\n')
        f.write('color grey80, {}\n'.format(first_structure))
        f.write('color chartreuse, {}\n'.format(' or '.join(all_first_mappings)))

        # Center on the originating structure
        f.write('center {}, -1\n'.format(all_mappings[0]))

    def calc_graph_statistics (self, graph):
        if graph.class_counts == None:
            return 0, 0, 0, 0, 0., 0., 0.

        tp = graph.class_counts[0]
        fn = graph.class_counts[1] - graph.class_counts[0]
        fp = graph.class_counts[2]
        tn = graph.class_counts[3] - graph.class_counts[2]

        # One-sided fisher's exact test, we only care about the probability of better matching contingency tables
        # Odds ratio: *unconditional* maximum likelihood estimate: (tp * tn) / (fp * fn)
        # R uses conditional maximum likelihood estimate instead:
        # http://scipy.github.io/devdocs/generated/scipy.stats.fisher_exact.html
        # https://github.com/wch/r-source/blob/trunk/src/library/stats/R/fisher.test.R
        odds_ratio, p_value = fisher_exact([[tp, fn], [fp, tn]], 'greater')

        # Matthews correlation coefficient
        mcc = (tp * tn - fp * fn) / m.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

        return tp, fn, fp, tn, mcc, odds_ratio, p_value

    def get_origin_mapping (self, graph):
        for db_id, db_mappings in graph.mappings.items():
            for (exact, score_contribution, mapping, pretty_mapping, random_mapping) in db_mappings:
                if exact:
                    return db_id, pretty_mapping

    def get_origin_mapping_and_pattern (self, graph, score_variant):
        origin_pretty_mapping = None
        origin_db_id = None
        best_pattern_score = 0.
        best_pattern_unmapped = 0
        best_pattern_match = None

        # Use the first exact mapping we find for pretty vertex names.
        # No point in using more than that
        for db_id, db_mappings in graph.mappings.items():
            for (exact, score_contribution, mapping, pretty_mapping, random_mapping) in db_mappings:
                if exact:
                    # Always at least have a pretty mapping if there is an exact mapping
                    if origin_pretty_mapping == None:
                        origin_pretty_mapping = pretty_mapping
                        origin_db_id = db_id
                        best_pattern_unmapped = len(pretty_mapping)

                    # Find best pattern match for the origin structure
                    if score_variant == 'percent':
                        pattern_matches = self.database[db_id].pattern_matches_percent
                    elif score_variant == 'binary':
                        pattern_matches = self.database[db_id].pattern_matches_binary
                    else:
                        raise NameError

                    # Sometimes the structure can't be matched to the sequence pattern
                    # because of missing residues. Ignore these structures when scoring.
                    if not pattern_matches or len(pattern_matches) < 1:
                        continue

                    for pattern_match in pattern_matches:
                        max_score = sum([x[0] for x in pattern_match.values()])
                        if max_score == 0.:
                            continue
                        score = 0.
                        unmapped = 0
                        for res in pretty_mapping:
                            try:
                                score += pattern_match[res][0]
                            except:
                                unmapped += 1
                        score /= max_score
                        if score > best_pattern_score or (score == best_pattern_score and unmapped <= best_pattern_unmapped):
                            origin_pretty_mapping = pretty_mapping
                            origin_db_id = db_id
                            best_pattern_score = score
                            best_pattern_unmapped = unmapped
                            best_pattern_match = pattern_match


        return (origin_db_id, origin_pretty_mapping, best_pattern_match)
        
    def calc_pattern_score (self, graph, score_variant):
        scores = {}
        unmappeds = {}

        for db_id, db_mappings in graph.mappings.items():
            if score_variant == 'percent':
                pattern_matches = self.database[db_id].pattern_matches_percent
            elif score_variant == 'binary':
                pattern_matches = self.database[db_id].pattern_matches_binary
            else:
                raise NameError

            # Sometimes the structure can't be matched to the sequence pattern
            # because of missing residues. Ignore these structures when scoring.
            if not pattern_matches or len(pattern_matches) < 1:
                continue

            scores[db_id] = 0.
            unmappeds[db_id] = len(graph.vertices)
            for pattern_match in pattern_matches:
                max_score = sum([x[0] for x in pattern_match.values()])
                if max_score == 0.:
                    continue
                for (_, _, _, pretty_mapping, _) in db_mappings:
                    score = 0.
                    unmapped = 0
                    for res in pretty_mapping:
                        try:
                            score += pattern_match[res][0]
                        except:
                            unmapped += 1
                    score /= max_score
                    if score > scores[db_id] or (score == scores[db_id] and unmapped < unmappeds[db_id]):
                        scores[db_id] = score
                        unmappeds[db_id] = unmapped

        score = sum(scores.values()) / len(scores.values())
        unmapped = sum(unmappeds.values()) / len(unmappeds.values())

        return (score, unmapped)

    def count_strictly_conserved (self, graph):
        # Consider all pattern entries here not just those present in the structure
        pattern_scores_binary = self.pattern_data[2]
        num_strictly_conserved = int(sum(pattern_scores_binary))

        scores = {}
        unmappeds = {}

        for db_id, db_mappings in graph.mappings.items():
            pattern_matches = self.database[db_id].pattern_matches_binary
            scores[db_id] = 0.
            unmappeds[db_id] = len(graph.vertices)
            for pattern_match in pattern_matches:
                for (_, _, _, pretty_mapping, _) in db_mappings:
                    score = 0.
                    unmapped = 0
                    for res in pretty_mapping:
                        try:
                            score += pattern_match[res][0]
                        except:
                            unmapped += 1
                    # no division by max score
                    if score > scores[db_id] or (score == scores[db_id] and unmapped < unmappeds[db_id]):
                        scores[db_id] = score
                        unmappeds[db_id] = unmapped

        # Consider all structures here, not just the ones with mappings
        sum_strictly_conserved = len(self.database) * num_strictly_conserved
        sum_strictly_conserved_matched = int(sum(scores.values()))

        return (sum_strictly_conserved, sum_strictly_conserved_matched)

    def calc_rmsd (self, graph):
        origin_db_id, origin_pretty_mapping = self.get_origin_mapping(graph)
        origin_structure = self.database[origin_db_id].get_structure()
        rmsds = []
        for db_id, db_mappings in graph.mappings.items():
            if db_id == origin_db_id:
                continue

            structure = self.database[db_id].get_structure()

            # Only consider the first/highest scoring mapping
            pretty_mapping = db_mappings[0][3]

            fixed = []
            moving = []

            for v_id in range(len(graph.vertices)):
                origin_chain, origin_res_id, origin_icode = origin_pretty_mapping[v_id]
                try:
                    origin_res = origin_structure[0][origin_chain][(' ', origin_res_id, ' ' if origin_icode == '' else origin_icode)]
                except:
                    sys.stderr.write('Missing residue in {}({}): ({},{})\n'.format(origin_structure, origin_chain, origin_res_id, origin_icode))
                    raise

                chain, res_id, icode = pretty_mapping[v_id]
                try:
                    res = structure[0][chain][(' ', res_id, ' ' if icode == '' else icode)]
                except:
                    sys.stderr.write('Missing residue in {}({}): ({},{})\n'.format(structure, chain, res_id, icode))
                    raise

                atom_names = [ a.get_name() for a in res.get_atoms() ]

                for origin_atom in origin_res.get_atoms():
                    origin_atom_name = origin_atom.get_name()
                    if origin_atom_name in atom_names:
                        fixed.append(origin_atom)
                        moving.append(res[origin_atom_name])

            superimposer = Bio.PDB.Superimposer()
            superimposer.set_atoms(fixed, moving)
            rmsds.append(superimposer.rms)

        avg_rmsd = sum(rmsds) / len(rmsds)
        return avg_rmsd


    def generate_pattern_html (self, graph, score_variant):
        origin_db_id, origin_pretty_mapping, best_pattern_match = self.get_origin_mapping_and_pattern(graph, score_variant)
        
        try:
            pattern_match_range = (min(best_pattern_match.keys()), max(best_pattern_match.keys()))
        except:
            pattern_match_range = None

        # Sort subgraph residues by their residue IDs
        res_list = []
        for v_id, labels in graph.vertices:
            res_id = origin_pretty_mapping[v_id]
            aa_type = self.reverse_mapping[labels[0]]
            res_list.append((res_id, aa_type))

        res_list = sorted(res_list)

        patterns = []
        patterns.append([])

        i = 0
        match_start_prefix = '<u>'
        match_end_suffix = '</u>'
        match_open = False
        chain_split = False
        while i < len(res_list):

            counter = 1
            res = res_list[i][0]
            res_chain = res[0]
            res_num = res[1]
            res_aa_type = res_list[i][1]

            item_label = res_aa_type

            if i < len(res_list) - 1:
                res2 = res_list[i + 1][0]
            else:
                res2 = res

            res2_chain = res2[0]
            res2_num = res2[1]

            # This is -1 for the last residue
            distance = res2_num - res_num - 1

            # Split patterns if they span multiple chains
            if res_chain != res2_chain:
                distance = 0
                chain_split = True

            if distance == 1:
                distance_label = 'x'
            elif distance > 1:
                distance_label = 'x({})'.format(distance)

            if pattern_match_range:
                # Determine matched groups
                try:
                    item_pattern_match_group = str(best_pattern_match[res][1])
                    item_label = '<span data-groups="{1}">{0}</span>'.format(item_label, item_pattern_match_group)
                except:
                    pass

                if distance >= 1:
                    distance_pattern_match_groups = set()
                    for x in best_pattern_match.keys():
                        if x > res and x < res2:
                            distance_pattern_match_groups.add(str(best_pattern_match[x][1]))
                    if distance_pattern_match_groups:
                        distance_label = '<span data-groups="{1}">{0}</span>'.format(distance_label, ','.join(sorted(distance_pattern_match_groups)))

                # Open match
                if not match_open:
                    # Match starts at residue
                    if res >= pattern_match_range[0] and res <= pattern_match_range[1]:
                        match_open = True
                        item_label = match_start_prefix + item_label
                    # Match somewhere in the following distance
                    elif distance >= 1 and res2 > pattern_match_range[0] and res2 <= pattern_match_range[1]:
                        match_open = True
                        distance_label = match_start_prefix + distance_label
                # Close match
                if match_open:
                    # Current residue is the end of the pattern string or match
                    if distance < 0 or res == pattern_match_range[1]:
                        match_open = False
                        item_label = item_label + match_end_suffix
                    # Match ends somewhere in the following distance
                    elif distance >= 1 and res2 > pattern_match_range[1]:
                        match_open = False
                        distance_label = distance_label + match_end_suffix

            patterns[-1].append(item_label)

            if distance >= 1:
                patterns[-1].append(distance_label)

            # Start new pattern for next chain
            if chain_split:
                patterns.append([])
                chain_split = False

            i += 1

        return ' + '.join(['-'.join(pattern) for pattern in patterns])

    def calc_graph_properties (self, graph):
        graph.support = len(graph.mappings.keys())
        graph.is_interesting = self.graph_is_interesting(graph)
        graph.pattern_html = self.generate_pattern_html(graph, 'percent')
        graph.rmsd = self.calc_rmsd(graph)

        for db_id, db_mappings in graph.mappings.items():
            for (exact, score_contribution, mapping, pretty_mapping, random_mapping) in db_mappings:
                if exact:
                    graph.origin_id = db_id
                    break

        graph.stats = self.calc_graph_statistics(graph)
        graph.is_significant = graph.stats[6] < self.significance_threshold

        graph.pattern_score_percent = 0.0
        graph.pattern_score_binary = 0.0
        graph.pattern_unmapped = 0.0
        graph.pattern_strictly_conserved = 0
        graph.pattern_strictly_conserved_matched = 0

        if self.pattern_data != None:
            graph.pattern_score_percent, graph.pattern_unmapped = self.calc_pattern_score(graph, 'percent')
            graph.pattern_score_binary = self.calc_pattern_score(graph, 'binary')[0]
            graph.pattern_strictly_conserved, graph.pattern_strictly_conserved_matched = self.count_strictly_conserved(graph)



    def generate_csv_string (self, f, graph, name):
        # 'id' 'name' 'size' 'vertices' 'support' 'score' 'interesting' 'significant'
        out = [graph.id, name, graph.size, len(graph.vertices), graph.support, graph.score, graph.is_interesting, graph.is_significant]
        # 'tp' 'fn' 'fp' 'tn' 'mcc' 'odds_ratio' 'p_value'
        out.extend(graph.stats)
        # 'pattern_score_percent' 'pattern_score_binary' 'pattern_exceed' 'pattern_conserved' 'pattern_conserved_matched' 'rmsd'
        out.extend([graph.pattern_score_percent, graph.pattern_score_binary, graph.pattern_unmapped, graph.pattern_strictly_conserved, graph.pattern_strictly_conserved_matched, graph.rmsd])

        if self.structman_classes and graph.id in self.structman_classes:
            core, lig, prot, surf, dna, rna, metal, ion = self.structman_classes[graph.id]
        else:
            core, lig, prot, surf, dna, rna, metal, ion = [ 0., 0., 0., 0., 0., 0., 0., 0. ]

        out.extend([ core, lig, prot, surf, dna, rna, metal, ion ])

        csv = '|'.join(map(str, out))
        f.write(csv + '\n')

    def generate_dot_string (self, graph, f, score_variant):
        origin_db_id, origin_pretty_mapping, best_pattern_match = self.get_origin_mapping_and_pattern(graph, score_variant)
        f.write(b'graph 0 {\n')

        f.write(b'mode = ipsep\n')
        f.write(b'overlap = ipsep\n')
        f.write(b'newrank = true\n')
        #f.write(b'nodesep = 2.0\n')
        f.write(b'mclimit = 2.0\n')
        f.write(b'bgcolor = "transparent"\n')
        f.write(b'remincross = true\n')
        f.write(b'levelsgap = -5\n')
        #f.write(b'ratio = 2.0\n')
        f.write(b'splines = true\n')
        f.write(b'node [fontname = "Arial"]\n')

        for v_id, labels in graph.vertices:
            label = self.reverse_mapping[labels[0]]
            color = 'white'
            if origin_pretty_mapping:
                mapped_residue = origin_pretty_mapping[v_id]
                #label = '<{}<BR/><FONT POINT-SIZE="8">{}</FONT>>'.format(label, convert_res_to_text(mapped_residue))
                label = '{}'.format(label)
                if best_pattern_match and mapped_residue in best_pattern_match:
                    color = 'lightblue'
            f.write('  {} [label={}, shape="circle" style="filled" fillcolor="{}"];\n'.format(v_id, label, color).encode('utf-8'))

        for v_id_1, v_id_2, labels in graph.edges:
            edge_style = ' [len={}]'.format(float(labels[1]) / 2000)
            # mark intractions as dashed lines
            if labels[0] == 1:
                edge_style += ' [style=dashed]'
            if labels[1] != 0:
                edge_style += ' [label=<<table cellpadding="1" border="0" cellborder="0"><tr><td>{:.3f}</td></tr></table>>]'.format(float(labels[1]) / 1000)
            f.write('  {} -- {}{};\n'.format(v_id_1, v_id_2, edge_style).encode('utf-8'))

        f.write(b'}')

    def generate_html_overview_header (self, f, title):
        f.write(('<!DOCTYPE html>\n'
                 '<html>\n'
                 '<head>\n'
                 '<title>{0}</title>\n'
                 '<style>\n'
                 '    body {{\n'
                 '        font-family: sans-serif;\n'
                 '    }}\n'
                 '    .graph-container {{\n'
                 '        display: flex;\n'
                 '        flex-wrap: wrap;\n'
                 '    }}\n'
                 '    .graph-container > div {{\n'
                 '        display: flex;\n'
                 '        flex-direction: column;\n'
                 '        justify-content: space-between;\n'
                 '        margin: 3px;\n'
                 '        padding: 5px;\n'
                 '        background-color: rgb(245, 245, 245);\n'
                 '        text-align: center;\n'
                 '    }}\n'
                 '    .graph-container > div.interesting {{\n'
                 '        background-color: rgb(136, 235, 140);\n'
                 '    }}\n'
                 '    .graph-container > div.significant {{\n'
                 '        box-shadow: 0px 0px 0 3px #fe3 inset;\n'
                 '    }}\n'
                 '    table.motif-properties {{\n'
                 '        margin-left: auto;\n'
                 '        margin-right: auto;\n'
                 '    }}\n'
                 '    table.motif-properties > tbody > tr > td:nth-child(1) {{\n'
                 '        text-align: left;\n'
                 '    }}\n'
                 '    table.motif-properties > tbody > tr > td:nth-child(2) {{\n'
                 '        text-align: right;\n'
                 '        padding-left: 20px;\n'
                 '    }}\n'
                 '    table.motif-properties > tbody > tr > td.pattern {{\n'
                 '        font-family: monospace;\n'
                 '        font-size: larger;\n'
                 '        font-weight: bold;\n'
                 '        text-align: center;\n'
                 '        padding-right: initial;\n'
                 '    }}\n').format(title))

        if self.pattern_data != None:
            f.write(('    #floating-pattern {\n'
                     '        background-color: rgb(250, 250, 250);\n'
                     '        border: 1px solid rgb(155, 155, 155);\n'
                     '        border-top: 0;\n'
                     '        position: fixed;\n'
                     '        padding: 5px;\n'
                     '        top: 0;\n'
                     '        left: 50%;\n'
                     '        transform: translateX(-50%);\n'
                     '        font-size: x-large;\n'
                     '    }\n'
                     '    #floating-pattern > .highlight {\n'
                     '        text-decoration: underline;\n'
                     '    }\n'
                     '    .visible {\n'
                     '        visibility: visible;\n'
                     '        opacity: 1;\n'
                     '        transition: opacity 0.3s ease;\n'
                     '    }\n'
                     '    .hidden {\n'
                     '        visibility: hidden;\n'
                     '        opacity: 0;\n'
                     '        transition: visibility 0s 0.3s, opacity 0.3s ease;\n'
                     '    }\n'))

        if self.structman_classes:
            f.write(('    .sm-classes {\n'
                     '        display: inline-block;\n'
                     '        white-space: nowrap;\n'
                     '        width: 100%;\n'
                     '        border: 0;\n'
                     '    }\n'
                     '    .sm-classes div {\n'
                     '        box-sizing: border-box;\n'
                     '        display: inline-block;\n'
                     '        border-width: 1px 0 1px 1px;\n'
                     '        border-color: black;\n'
                     '        border-style: solid;\n'
                     '        padding: 3px;\n'
                     '        overflow: hidden;\n'
                     '        white-space: nowrap;\n'
                     '        text-align: center;\n'
                     '    }\n'
                     '    .sm-classes div:last-child {\n'
                     '        border-width: 1px 1px 1px 1px;\n'
                     '    }\n'
                     '    .sm-classes .surf {\n'
                     '        background-color: rgb(250, 200, 120);\n'
                     '    }\n'
                     '    .sm-classes .core {\n'
                     '        background-color: rgb(230, 130, 100);;\n'
                     '    }\n'
                     '    .sm-classes .dna {\n'
                     '        background-color: rgb(255, 240, 50);\n'
                     '    }\n'
                     '    .sm-classes .rna {\n'
                     '        background-color: rgb(250, 135, 180);\n'
                     '    }\n'
                     '    .sm-classes .lig {\n'
                     '        background-color: rgb(110, 190, 230);\n'
                     '    }\n'
                     '    .sm-classes .prot {\n'
                     '        background-color: rgb(136, 235, 140);\n'
                     '    }\n'
                     '    .sm-classes .metal {\n'
                     '        background-color: rgb(255, 255, 255);\n'
                     '    }\n'
                     '    .sm-classes .ion {\n'
                     '        background-color: rgb(140, 240, 230);\n'
                     '    }\n'))


        f.write(('</style>\n'
                 '</head>\n'
                 '<body>\n'
                 '<h1>{0}</h1>\n').format(title))

        if self.pattern != None:
            pattern = []
            for group, entry in enumerate(self.pattern.strip().split('-')):
                pattern.append('<span id="group-{}">{}</span>'.format(group, entry))
            f.write('<div id="floating-pattern">{}</div>\n'.format('-'.join(pattern)))

        f.write('<div class="graph-container">\n')

    def generate_html_overview_footer (self, f, output):
        f.write('</div>\n')

        if os.path.isfile(output):
            with open(output, 'r') as output_file:
                f.write('<p>\n<span style="font-size: larger;"><b>Output:</b></span>\n<pre>')
                for line in output_file:
                    f.write(line)
                f.write('</pre>')

        if self.pattern_data != None:
            f.write(('<script>\n'
                     'function highlight_groups() {\n'
                     'var groups = this.dataset.groups.split(\',\')\n'
                     '    for (group_id in groups) {\n'
                     '        group_element = document.getElementById("group-" + groups[group_id]);\n'
                     '        if (group_element)\n'
                     '            group_element.classList.add(\'highlight\');\n'
                     '    }\n'
                     '}\n'
                     '\n'
                     'function unhighlight_groups() {\n'
                     '    var groups = this.dataset.groups.split(\',\')\n'
                     '    for (group_id in groups) {\n'
                     '        group_element = document.getElementById("group-" + groups[group_id]);\n'
                     '        if (group_element)\n'
                     '            group_element.classList.remove(\'highlight\');\n'
                     '    }\n'
                     '}\n'
                     '\n'
                     'for (span_id in document.getElementsByTagName("span")) {\n'
                     '    var span = document.getElementsByTagName("span")[span_id];\n'
                     '    if (span.dataset && span.dataset.groups) {\n'
                     '        span.addEventListener(\'mouseover\', highlight_groups);\n'
                     '        span.addEventListener(\'mouseout\', unhighlight_groups);\n'
                     '    }\n'
                     '}\n'
                     '\n'
                     'function hide_pattern() {\n'
                     '    var pattern = document.getElementById("floating-pattern")\n'
                     '    pattern.classList.remove(\'visible\');\n'
                     '    pattern.classList.add(\'hidden\');\n'
                     '}\n'
                     '\n'
                     'function unhide_pattern() {\n'
                     '    var pattern = document.getElementById("floating-pattern")\n'
                     '    pattern.classList.remove(\'hidden\');\n'
                     '    pattern.classList.add(\'visible\');\n'
                     '}\n'
                     '\n'
                     'document.getElementsByTagName("h1")[0].addEventListener(\'mouseover\', hide_pattern);\n'
                     'document.getElementsByTagName("h1")[0].addEventListener(\'mouseout\', unhide_pattern);\n'
                     '</script>\n'))

        f.write('</body>\n')
        f.write('</html>\n')

    def generate_html_overview_entry (self, f, graph, svg_file, pml_file):
        style_classes = []
        if graph.is_interesting:
            style_classes.append('interesting')
        if graph.is_significant:
            style_classes.append('significant')
        style_classes_property = ''
        if len(style_classes) > 0:
            style_classes_property = ' class="{}"'.format(' '.join(style_classes))

        f.write('<div{}>\n'.format(style_classes_property))

        image = '<img src="{}">'.format(svg_file)
        if pml_file:
            image = '<a href="{}">{}</a>'.format(pml_file, image)
        f.write(image + '\n')

        f.write('<table class="motif-properties">\n')
        f.write('<tr><td colspan="2" class="pattern">{}</td></tr>\n'.format(graph.pattern_html))
        f.write('<tr><td>Score:</td><td>{}</td></tr>\n'.format(graph.score))
        f.write('<tr><td>Support:</td><td>{}/{}</td></tr>\n'.format(graph.support, len(self.database)))

        if self.use_significance:
            tp, fn, fp, tn, mcc, odds_ratio, p_value = graph.stats

            tpr = tp / (tp + fn)
            f.write('<tr><td>TPR:</td><td>{:.1f} ({:d}/{:d})</td></tr>\n'.format(tpr, tp, tp + fn))

            fpr = fp / (fp + tn)
            f.write('<tr><td>FPR:</td><td>{:.1f} ({:d}/{:d})</td></tr>\n'.format(fpr, fp, fp + tn))

            f.write('<tr><td>p-value:</td><td>{:.3e}</td></tr>\n'.format(p_value))
            f.write('<tr><td>Odds ratio:</td><td>{:.1f}</td></tr>\n'.format(odds_ratio))
            f.write('<tr><td>MCC:</td><td>{:.5f}</td></tr>\n'.format(mcc))

        f.write('<tr><td>Origin:</td><td>{}</td></tr>\n'.format(self.database[graph.origin_id].pretty_id))
        f.write('<tr><td>Size:</td><td>{:d}</td></tr>\n'.format(graph.size))

        f.write('<tr><td>RMSD:</td><td>{:.3f}</td></tr>\n'.format(graph.rmsd))

        if self.pattern_data != None:
            f.write('<tr><td>Pattern score (percent):</td><td>{:.5f}</td></tr>\n'.format(graph.pattern_score_percent))
            f.write('<tr><td>Pattern score (binary):</td><td>{:.5f}</td></tr>\n'.format(graph.pattern_score_binary))
            f.write('<tr><td>Non-pattern matches:</td><td>{:.5f}</td></tr>\n'.format(graph.pattern_unmapped))

        if self.structman_classes and graph.id in self.structman_classes:
            core, lig, prot, surf, dna, rna, metal, ion = self.structman_classes[graph.id]
            f.write('<tr><td colspan="2"><div class="sm-classes">')
            for sm_class_id, percentage in zip(['surf', 'core', 'dna', 'rna', 'lig', 'prot', 'metal', 'ion'], [surf, core, dna, rna, lig, prot, metal, ion]):
                if percentage > 0.0:
                    f.write('<div class="{0}" style="width: {1:.2%};" title="{0}: {1:.2%}">{1:.2%}</div>'.format(sm_class_id, percentage))
            f.write('</div></td></tr>\n')

        f.write('</table>\n')
        f.write('</div>\n')

    def generate_smlf_string (self, f, f_random, graph, family, prefix, additional_tags):
        tags = [ 'size-' + str(graph.size) ]

        if family:
            tags.append('fam-' + family)

        if graph.is_interesting:
            tags.append('interesting')
        else:
            tags.append('not-interesting')

        if graph.is_significant:
            tags.append('significant')
        else:
            tags.append('not-significant')

        if graph.is_best:
            tags.append('best')

        if graph.is_top:
            tags.append('top')

        if graph.rmsd <= 3.0:
            tags.append('rmsd-3')
        elif graph.rmsd <= 6.0:
            tags.append('rmsd-3-6')
        elif graph.rmsd <= 9.0:
            tags.append('rmsd-6-9')
        else:
            tags.append('rmsd-9')

        if additional_tags:
            tags.extend(additional_tags)

        tag_combos = []

        for size in range(1, len(tags) + 1):
            for combined_tag in combinations(tags, size):
                tag_combos.append('_'.join(combined_tag))

        for db_id, db_mappings in graph.mappings.items():
            # Use the highest scoring mapping for each structure
            exact, score_contribution, mapping, pretty_mapping, random_mapping = sorted(db_mappings, key=lambda x: (x[0], x[1]), reverse=True)[0]
            pretty_db_id = self.database[db_id].pretty_id
            for v_id, labels in graph.vertices:
                res_type = self.reverse_mapping[labels[0]]
                chain, res_id, icode = pretty_mapping[v_id]
                f.write('{0}:{1}\t{2}{3}{4}\t\t{5}\n'.format(pretty_db_id, chain, res_type, res_id, icode, ','.join(tag_combos + [ prefix + str(graph.id) ])))

                chain, res_id, icode = random_mapping[v_id]
                f_random.write('{0}:{1}\t?{2}{3}\t\t{4}\n'.format(pretty_db_id, chain, res_id, icode, ','.join(tag_combos)))


    # Simple heuristic to indicate which subgraphs might be more worth manual
    # investigation
    def graph_is_interesting(self, graph):
        if graph.size < 3:
            return False

        if graph.score < (graph.size * self.score_difference / 2):
            return False

        boring_counter = 0
        for vertex in graph.vertices:
            v_label = self.reverse_mapping[vertex[1][0]]
            if 'I' in v_label or 'L' in v_label or 'V' in v_label or 'A' in v_label:
                boring_counter += 1

        if boring_counter >= graph.size / 2:
            return False

        return True

    def _parse_graphs (self, graphs_file_path, offsets, graph_ids, base_dir):
        with open(graphs_file_path, 'r', encoding='ascii') as graphs_file:
            # Parse the database header
            db_pretty_vertices = []
            db_pretty_ids = []
            db_seg_file_paths = []
            db_pdb_file_paths = []

            for line in graphs_file:
                if line[0] == '#':
                    splt = line.strip().split(' ')
                    if splt[1] == 'r':
                        pretty_vertices = [convert_text_to_res(x) for x in splt[2:]]
                        db_pretty_vertices.append(pretty_vertices)
                    elif splt[1] == 'id':
                        pdb_id = splt[2]
                        db_pretty_ids.append(pdb_id)
                    elif splt[1] == 'pdb':
                        pdb_path = splt[2]
                        if pdb_path[0] != '/':
                            pdb_path = base_dir + '/' + pdb_path
                        db_pdb_file_paths.append(pdb_path)
                    elif splt[1] == 'seg':
                        seg_path = splt[2]
                        if seg_path[0] != '/':
                            seg_path = base_dir + '/' + seg_path
                        db_seg_file_paths.append(seg_path)

                elif line == '\n':
                    pass

                elif line[0] == 't':
                    # End of header
                    break

                else:
                    sys.stderr.write('Unhandled input format: {0}\n'.format(line))

            for pretty_id, pretty_vertices, seg_file_path, pdb_file_path in zip(db_pretty_ids, db_pretty_vertices, db_seg_file_paths, db_pdb_file_paths):
                self.database.append(DatabaseEntry(pretty_id, pretty_vertices, seg_file_path, pdb_file_path, self.pattern_data))

            # Parse graphs using their offsets
            for graph_id in graph_ids:
                graph = Graph(graph_id)
                mappings = {}

                num_t = 0

                graphs_file.seek(offsets[graph_id])
                line = graphs_file.readline()

                while line:
                    if line[0] == 't':
                        num_t += 1
                        if num_t > 1:
                            break

                    elif line[0] == 'v':
                        splt = [int(x) for x in line[2:].split()]
                        id_, labels = splt[0], splt[1:]
                        graph.vertices.append((id_, labels))

                    elif line[0] == 'e':
                        splt = [int(x) for x in line[2:].split()]
                        id_A, id_B, labels = splt[0], splt[1], splt[2:]
                        graph.edges.append((id_A, id_B, labels))
                        graph.size += 1

                    elif line[0] == '#':
                        splt = line.strip().split()
                        if splt[1] == 'm':
                            try:
                                db_id, exact, score_contribution, mapping = line[4:].split(' : ')
                                exact = bool(int(exact))
                                score_contribution = int(score_contribution)
                            except:
                                db_id, mapping = line[4:].split(' : ')
                                exact = False
                                score_contribution = 0

                            db_id = int(db_id)

                            if db_id not in mappings:
                                mappings[db_id] = []

                            mapping = [int(x) for x in mapping.split()]

                            pretty_mapping = None
                            if db_id < len(self.database):
                                pretty_mapping = [self.database[db_id].pretty_vertices[x] for x in mapping]
                                random_mapping = [self.database[db_id].shuffled_vertices[x] for x in mapping]

                            mappings[db_id].append((exact, score_contribution, mapping, pretty_mapping, random_mapping))
                        elif splt[1] == 'p':
                            graph.parents = [int(x) for x in splt[2:]]
                        elif splt[1] == 's':
                            graph.score = int(splt[2])
                        elif splt[1] == 'i':
                            graph.class_counts = [int(x) for x in splt[2:]]

                    elif line == '\n':
                        pass

                    else:
                        sys.stderr.write('Unhandled input format: {0}\n'.format(line))

                    line = graphs_file.readline()

                # sort mappings with highest score first
                for db_id in mappings.keys():
                    mappings[db_id].sort(key=lambda x: x[1], reverse=True)

                graph.mappings = mappings

                self.calc_graph_properties(graph)

                self.graphs.append(graph)



parser = argparse.ArgumentParser(description='Create dot files from a rinminer result file')
parser.add_argument('-g', '--grouping', type=int, choices=[0, 1, 2, 3], default=0, help='residue type grouping (1 = least specific grouping, ..., 3 = most specific grouping, 0 = no grouping)')
parser.add_argument('-s', '--structure', action="store_true", help='use secondary structure information in the labels')
parser.add_argument('-f', '--force', action="store_true", help='force clearing the directory first')
parser.add_argument('-b', '--batch', action="store_true", help='no interactive input required, only requested output')
parser.add_argument('-i', '--interesting', action="store_true", help='give measure of interestingness of results to \'interestingness\'')
parser.add_argument('-p', '--pymol', action="store_true", help='generate pymol scripts')
parser.add_argument('-d', '--score-diff', type=int, default=900, help='minimum score difference required in rinminer; used to calculate interestingness')
parser.add_argument('-c', '--csv', action="store_true", help='csv output of all graph statistics')
parser.add_argument('-a', '--pattern', default=None, help='pattern used for evaluation')
parser.add_argument('-y', '--base-dir', default="/", help='base directory used for relative paths')
parser.add_argument('-t', '--title', default="Subgraphs", help='title used in the overview html')
parser.add_argument('-n', '--significance-threshold', default=None, type=float, help='threshold used to determine significance of classification based on subgraph')
parser.add_argument('-o', '--output', default=None, help='path to the rinminer terminal output')
parser.add_argument('-l', '--family', default='', help='family identifier for SMLF tags')
parser.add_argument('-x', '--prefix', default='', help='prefix used for the unique SMLF ID')
parser.add_argument('-m', '--smlf', action="store_true", help='generate SMLF files')
parser.add_argument('-w', '--structman', default=None, help='path to structman simple classification (requires correct prefix specified)')
parser.add_argument('-z', '--top', action="store_true", help='the highest scoring subgraph from the input should be marked as top result')
parser.add_argument('-u', '--tags', default=None, help='additional tags to be used in te SMLF output')

parser.add_argument('in_file', metavar='result.lg', help='name of .lg input file')
parser.add_argument('in_index_file', metavar='result_index.txt', help='graph index file')
parser.add_argument('out_dir', metavar='output_dir', help='path to the directory for the dot files')
parser.add_argument('ids', metavar='ids', nargs='*', type=int, help='candidate IDs')
args = parser.parse_args()

svg_files = []
pml_files = []

try:
    try:
        os.makedirs(args.out_dir)
    except:
        pass

    offsets = {}
    if len(args.ids) > 0:
        with open(args.in_index_file, 'r') as index_file:
            for line in index_file:
                g_id, offset = [ int(x) for x in line.strip().split() ]
                offsets[g_id] = offset

    if not os.path.isdir(args.out_dir):
        sys.stderr.write('Invalid output directory\n')
        sys.exit(1)

    if os.listdir(args.out_dir):
        clear = False
        if args.force:
            if not args.batch:
                print("Target directory ({}) not empty.".format(args.out_dir))
            clear = True
        elif not args.batch:
            sys.stdout.write("Target directory ({}) not empty. Clear directory [Y/n]? ".format(args.out_dir))
            sys.stdout.flush()
            inp = sys.stdin.readline().strip()
            if inp.upper() == 'Y' or inp == '':
                clear = True

        if clear:
            if not args.batch:
                print("Clearing {}".format(args.out_dir))
            for f in os.listdir(args.out_dir):
                file_path = os.path.join(args.out_dir, f)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                except Exception as e:
                    sys.stderr.write(e)

    if args.grouping == 0:
        grouping_scheme = individual_types
    elif args.grouping == 1:
        grouping_scheme = type_groups1
    elif args.grouping == 2:
        grouping_scheme = type_groups2
    elif args.grouping == 3:
        grouping_scheme = type_groups3

    rm = gen_reverse_mapping(grouping_scheme, args.structure)

    if args.pattern:
        pattern_data = parse_pattern(args.pattern)
    else:
        pattern_data = None

    additional_tags = None
    if args.tags:
        additional_tags = args.tags.split(',')

    class_headers = dict((x,i) for i,x in enumerate(['Core', 'Ligand Interaction', 'Protein Interaction', 'Surface', 'DNA Interaction', 'RNA Interaction', 'Metal Interaction', 'Ion Interaction'], start=0))
    structman_classes = None
    if args.structman:
        structman_classes = {}
        with open(args.structman, 'r') as sm_file:
            prefix = args.prefix
            class_mapping = {}
            remaining_classes = []
            input_classes = next(sm_file).rstrip().split('\t')[1:]
            for i, input_class in enumerate(input_classes, start=0):
                if input_class in class_headers:
                    class_mapping[i] = class_headers[input_class]
                else:
                    remaining_classes.append(i)
            for line in sm_file:
                if line.startswith(prefix):
                    splt = line.strip().split('\t')
                    tag = splt[0]
                    g_id = int(splt[0][len(prefix):])
                    parsed_classes = [ float(x) for x in splt[1:] ]
                    orderd_classes = [ 0. ] * len(class_headers)
                    for frm, to in class_mapping.items():
                        orderd_classes[to] = parsed_classes[frm]
                    for x in remaining_classes:
                        if parsed_classes[x] != 0.:
                            sys.stderr.write('\'{}\' non-zero for tag \'{}\': {}\n'.format(input_classes[x], tag, parsed_classes[x]))
                            raise ValueError
                    structman_classes[g_id] = tuple(orderd_classes)

    collection = GraphsCollection(args.in_file, offsets, args.ids, rm, args.pattern, pattern_data, structman_classes, args.base_dir, args.score_diff, args.significance_threshold, args.top)

    # too lazy to solve this properly :)
    max_len = len(str(len(collection.graphs)))

    overview_path = '{}/overview.html'.format(args.out_dir)
    overview_file = open(overview_path, 'w')

    collection.generate_html_overview_header(overview_file, args.title)

    if args.csv:
        csv_path = '{}/overview.csv'.format(args.out_dir)
        csv_file = open(csv_path, 'w')

    if args.smlf:
        smlf_path = '{}/overview.smlf'.format(args.out_dir)
        smlf_file = open(smlf_path, 'w')
        smlf_random_path = '{}/overview_random.smlf'.format(args.out_dir)
        smlf_random_file = open(smlf_random_path, 'w')

    for number, graph in enumerate(collection.graphs, start=1):
        name = '{}_{}_{}'.format(str(number).zfill(max_len), graph.size, graph.id)

        svg_path = '{}/{}.svg'.format(args.out_dir, name)
        proc = subprocess.Popen(['dot', '-Tsvg', '-o{}'.format(svg_path)], stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        collection.generate_dot_string(graph, proc.stdin, 'percent')
        proc.stdin.close()

        dot_path = '{}/{}.dot'.format(args.out_dir, name)
        with open(dot_path, 'wb') as dot_file:
            collection.generate_dot_string(graph, dot_file, 'percent')

        pml_path = None
        if args.pymol:
            pml_path = '{}/{}.pml'.format(args.out_dir, name)
            pml_file = open(pml_path, 'w')
            collection.generate_pymol_script(pml_file, graph)
            pml_file.close()

        if args.smlf:
            collection.generate_smlf_string(smlf_file, smlf_random_file, graph, args.family, args.prefix, additional_tags)

        if args.csv:
            collection.generate_csv_string(csv_file, graph, name)

        if pml_path:
            pml_path = os.path.relpath(pml_path, start=os.path.dirname(overview_path))

        svg_path = os.path.relpath(svg_path, start=os.path.dirname(overview_path))

        collection.generate_html_overview_entry(overview_file, graph, svg_path, pml_path)

    if args.smlf:
        smlf_file.close()

    if args.csv:
        csv_file.close()

    collection.generate_html_overview_footer(overview_file, args.output)

except FileNotFoundError:
    if args.batch:
        pass
    else:
        raise
