#!/usr/bin/python3

import argparse
import os, sys
import urllib
import gzip
import errno
from urllib.request import urlretrieve

parser = argparse.ArgumentParser(description='Generate a dataset based on a list of PDB IDs (with chains and ranges specified in SCOP format)')
parser.add_argument('pdb_ids_file', nargs='+', metavar='PDB_IDs_FILE', help='Path to files containing PDB IDs')
parser.add_argument('target_dir', metavar='dataset_dir', help='Target directory')
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


try:
    os.makedirs(args.target_dir)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass

for pdb_ids_path in args.pdb_ids_file:
    pdb_ids_file = open(pdb_ids_path, 'rt')
    family = os.path.basename(pdb_ids_path).split('.')[0]
    family_dir = args.target_dir + '/' + family

    try:
        os.makedirs(family_dir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    name_file = open(family_dir + '/family_name.txt', 'w')
    name_file.write(family)
    name_file.close()   

    entries = []

    for line in pdb_ids_file:
        if line[0] == '#':
            continue
        pdb_id, temp_locations = line.strip().split(' ')
        pdb_id = pdb_id.upper()
        locations = []
        for x in temp_locations.split(','):
            chain, rnge = x.split(':')
            if rnge == '':
                rnge = None
            locations.append((chain, rnge))
        entries.append((pdb_id, locations))

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
                m = re.match('(-?[0-9]*)-(-?[0-9]*)', location[1])
                rnge = [m.group(1), m.group(2)]
            except:
                rnge = ['_', '_']
            segments.append('{},{},{},{}'.format(0, chain, rnge[0], rnge[1]))

        segments_file = open(dir_path + '/' + pdb_id + '.seg', 'w')
        segments_file.write('\n'.join(segments))
        segments_file.close()
