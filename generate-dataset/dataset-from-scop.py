#!/usr/bin/python3

import argparse
import os, sys
import gzip
import errno
import re
from urllib.request import urlretrieve
import pisces

parser = argparse.ArgumentParser(description='Generate a dataset based on SCOP families')
parser.add_argument('blastdb_path', metavar='blastdb_path', help='Path to PISCES BLASTDB directory')
parser.add_argument('scop_des', metavar='scop_des_file', help='Path to dir.des.scop.*.txt')
parser.add_argument('sfamilies', metavar='sfamilies_file', help='Path to file containing the selected SCOP (super) families')
parser.add_argument('target_dir', metavar='dataset_dir', help='Directory in which the data set is created')
parser.add_argument('-m', '--max-res', type=float, default=3.0, help='Maximum resolution')
parser.add_argument('-r', '--max-r-value', type=float, default=1.0, help='Maximum R-value')
parser.add_argument('-l', '--min-length', type=int, default=0, help='Minimum length')
parser.add_argument('-a', '--aux', type=int, default=50, help='Auxiliary chain length threshold. Auxiliary chains are not subject to sequence identity thresholds.')
parser.add_argument('-i', '--max-identity', type=int, default=30, help='Maximum identity')
parser.add_argument('-e', '--min-entries', type=int, default=3, help='Minimum number of entries in a family')
parser.add_argument('-p', '--pdb-mirror', help='Local PDB mirror path instead of downloading files from the RCSB PDB webserver.')
args = parser.parse_args()

def create_pdb_file(pdb_id, path):
    if args.pdb_mirror:
        mirror_path = args.pdb_mirror + '/data/structures/divided/pdb/' + pdb_id.lower()[1:3] + '/pdb' + pdb_id.lower() + '.ent.gz'
        if os.path.isfile(mirror_path):
            try:
                os.remove(path)
            except OSError as exc:
                if exc.errno != errno.ENOENT:
                    raise
                pass
            os.symlink(mirror_path, path)
            return
        else:
            print(mirror_path)
            sys.stderr.write('PDB ID \'{}\' not found in local database, attempting download.\n'.format(pdb_id))
    urlretrieve ('http://files.rcsb.org/download/' + pdb_id + '.pdb.gz', path)

# Read names and IDs of all super families used to generate the data set
families = {}
family_names = {}
with open(args.sfamilies, 'r') as sfamilies_file:
    for line in sfamilies_file:
        fam_id, fam_name = line.strip().split(';', 1)
        families[fam_id] = []
        family_names[fam_id] = fam_name

# Add entries from the scop database into their respective families
with open(args.scop_des, 'r') as scop_des_file:
    for line in scop_des_file:
        if line[0] == '#':
            continue
        line = line.strip().split('\t')
        if line[1] == 'px':
            family = None
            for f in families:
                if line[2].startswith(f):
                    family = f
                    break
            if family == None:
                continue
            entry = line[4]
            pdb_id, temp_locations = entry.split(' ')
            pdb_id = pdb_id.upper()
            if not pdb_id[0].isdigit():
                sys.stderr.write('Skipping protein \'{}\': Invalid PDB ID\n'.format(pdb_id))
                continue
            locations = []
            for x in temp_locations.split(','):
                chain, rnge = x.split(':')
                if rnge == '':
                    rnge = None
                locations.append((chain, rnge))
            families[family].append((pdb_id, locations))

# Filter out similar structures
pisces = pisces.Pisces(args.blastdb_path, args.max_res, args.max_r_value, args.min_length, args.aux, args.max_identity)
filtered_families = {}
removed_families = []
for family, entries in families.items():
    if len(entries) < args.min_entries:
        removed_families.append(family)
        continue

    fam_pdb_chain_ids = set([ (x[0], tuple([ l[0] for l in x[1] ])) for x in entries ])

    keep, stats = pisces.filter(fam_pdb_chain_ids)
    keep = set(keep)

    # Now filter all entries in the family using the non-redundant set of kept structures
    filtered_entries = []
    for entry in entries:
        chain_ids = (entry[0], tuple(sorted([ l[0] for l in entry[1] ])))
        if chain_ids in keep:
            filtered_entries.append(entry)

    # Remove families with too few protein structures remaining
    if len(filtered_entries) < args.min_entries:
        removed_families.append(family)
        continue

    filtered_families[family] = filtered_entries


for family in removed_families:
    sys.stderr.write('Skipping family \'{}\': Too many structures culled\n'.format(family))


# Create dataset
for family, entries in filtered_families.items():
    family_dir = args.target_dir + '/' + family
    try:
        os.makedirs(family_dir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    name_file = open(family_dir + '/family_name.txt', 'w')
    name_file.write(family_names[family])
    name_file.close()   

    used_names = []

    for entry in entries:
        pdb_id = entry[0]
        locations = entry[1]

        name = pdb_id + '_' + locations[0][0]
        full_name = name
        counter = 2

        while full_name in used_names:
            full_name = name + str(counter)
            counter += 1

        used_names.append(full_name)

        dir_path = family_dir + '/' + full_name

        try:
            os.makedirs(dir_path)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        pdb_file_path = dir_path + '/' + pdb_id + '.pdb.gz'
        try:
            create_pdb_file (pdb_id, pdb_file_path)
        except Exception as e:
            sys.stderr.write('\'{}\' removed - could not obtain PDB structure: {}\n'.format(pdb_id, e))
            try:
                os.rmdir(dir_path)
            except:
                pass
            continue

        segments = []
        for location in locations:
            chain = location[0]
            try:
                m = re.match('(-?[0-9]+[a-zA-Z]?)-(-?[0-9]+[a-zA-Z]?)', location[1])
                rnge = [m.group(1), m.group(2)]
            except:
                rnge = ['_', '_']
            segments.append('{},{},{},{}'.format(0, chain, rnge[0], rnge[1]))

        segments_file = open(dir_path + '/' + pdb_id + '.seg', 'w')
        segments_file.write('\n'.join(segments))
        segments_file.close()
