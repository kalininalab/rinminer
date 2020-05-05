#!/usr/bin/python3

import argparse
import os, sys
import re
import gzip
import errno
from Bio import SeqIO
import pisces
import tempfile
import shutil
from urllib.request import urlretrieve

def pattern_to_regex(pattern):
    pattern = pattern.replace('-','')
    pattern = re.sub('\[([^[]*)>([^]]*)\]', '([\g<1>\g<2>]|$)', pattern)
    pattern = re.sub('\[([^[]*)<([^]]*)\]', '([\g<1>\g<2>]|^)', pattern)
    pattern = pattern.replace('<','^').replace('>','$')
    pattern = pattern.replace('{','[^').replace('}',']')
    pattern = re.sub('\(([0-9,]*)\)', '{\g<1>}', pattern)
    pattern = pattern.replace('x','.')
    pattern = '.*' + pattern + '.*'
    return pattern

def get_chains_for_pattern(pattern, pdb_path):
    regex = pattern_to_regex(pattern)
    r = re.compile(regex)

    include_chains = []

    with gzip.open(pdb_path, 'rt') as pdb_file:
        for record in SeqIO.parse(pdb_file, "pdb-seqres"):
            chain = record.annotations["chain"]
            chain_seq = str(record.seq).strip()
            if r.match(chain_seq):
                include_chains.append(chain)

    return include_chains

def get_chains_for_refs(refs, pdb_path):
    include_chains = []

    with gzip.open(pdb_path, 'rt') as pdb_file:
        for line in pdb_file:
            line = line.strip()
            if line.startswith('DBREF'):
                chain = line[12]
                ref = line[33:42].strip()
                if ref in refs:
                    include_chains.append(chain)
    return include_chains

parser = argparse.ArgumentParser(description='Generate a dataset based on Prosite patterns')
parser.add_argument('blastdb_path', metavar='blastdb_path', help='Path to PISCES BLASTDB directory')
parser.add_argument('prosite_path', metavar='prosite.dat', help='Path to prosite.dat')
parser.add_argument('target_dir', metavar='dataset_dir', help='Directory in which the data set is created')
parser.add_argument('-m', '--max-res', type=float, default=3.0, help='Maximum resolution')
parser.add_argument('-r', '--max-r-value', type=float, default=1.0, help='Maximum R-value')
parser.add_argument('-l', '--min-length', type=int, default=0, help='Maximum length')
parser.add_argument('-i', '--max-identity', type=int, default=50, help='Maximum identity')
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

# used to filter out similar structures
pisces = pisces.Pisces(args.blastdb_path, args.max_res, args.max_r_value, args.min_length, args.max_identity)

# parse prosite entries
entries = []
with open(args.prosite_path, 'rt') as prosite_file:
    entry_id = None
    entry_type = None
    for line in prosite_file:
        line = line.strip()
        line_type = line[0:2]
        if line_type == 'ID':
            entry = line[2:].split(';')
            entry_id = entry[0].strip()
            entry_type = entry[1].rstrip('.').strip()
            pattern = ''
            structures = []
            refs = []
        elif line_type == 'AC':
            acc = line[3:].rstrip(';').strip()
        elif line_type == 'PA':
            pattern += line[3:].rstrip('.').strip()
        elif line_type == 'DR':
            entry = filter(None, line[3:].split(';'))
            for x in entry:
                ref = [r.strip() for r in x.split(',')]
                # only use true positives
                if ref[2] == 'T':
                    refs.append(ref[0])
        elif line_type == '3D':
            structures = list(filter(None, map(str.strip, line[2:].split(';'))))
        elif line_type == '//':
            if entry_type == 'PATTERN' and len(structures) > 0 and len(refs) > 0:
                entries.append((entry_id, acc, pattern, set(refs), set(structures)))
            entry_id = None
            entry_type = None
    if entry_type == 'PATTERN' and len(structures) > 0 and len(refs) > 0:
        entries.append((entry_id, acc, pattern, set(refs), set(structures)))

tmp_dir = tempfile.TemporaryDirectory()
pdb_files = {}

for num, entry in enumerate(entries, 1):
    name, accession_number, pattern, refs, structures =  entry
    chains = []

    print('{}/{}: {} ({}); {}'.format(num, len(entries), name, accession_number, pattern))

    # determine matching chains, because prosite only tells us the pdb id
    for pdb_id in structures:
        if pdb_id not in pdb_files:
            try:
                pdb_file_path = tmp_dir.name + '/' + pdb_id + '.pdb.gz'
                create_pdb_file (pdb_id, pdb_file_path)
                pdb_files[pdb_id] = pdb_file_path
            except Exception as e:
                sys.stderr.write('Structure \'{}\' could not be found: {}\n'.format(pdb_id, e))
                continue
        else:
            pdb_file_path = pdb_files[pdb_id]
        try:
            #struct_chains_pat = get_chains_for_pattern(pattern, pdb_id)
            #if len(struct_chains_pat) == 0:
                #sys.stderr.write('Pattern \'{}\' not found in structure \'{}\'\n'.format(pattern, pdb_id))
                #sys.exit(1)
            struct_chains = get_chains_for_refs(refs, pdb_file_path)
            #if len(struct_chains) == 0:
                #sys.stderr.write('No matching chains found in structure \'{}\'\n'.format(pdb_id))
                #sys.exit(1)
            #if len(struct_chains) != len(struct_chains_pat):
                #sys.stderr.write('Pattern/references mismatch ({}/{}) for pattern \'{}\' in structure \'{}\'\n'.format(len(struct_chains_pat), len(struct_chains), pattern, pdb_id))
                #sys.exit(1)
            struct_chains = [pdb_id + struct_chain for struct_chain in struct_chains]
            chains.extend(struct_chains)
        except IOError:
            sys.stderr.write('%s not found, skipping\n' % pdb_id)

    keep, stats = pisces.filter(chains)

    print('{}/{} chains with at most {}% identity. {} culled for experimental reasons, {} culled for identity'.format(len(keep), len(chains), args.max_identity, stats[0], stats[1]))

    if len(keep) < args.min_entries:
        sys.stderr.write('Skipping \'{}\': Too many structures culled\n'.format(accession_number))
        print('----------')
        continue


    # create data set
    family_dir = args.target_dir + '/' + accession_number
    try:
        os.makedirs(family_dir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    name_file = open(family_dir + '/family_name.txt', 'w')
    name_file.write(name)
    name_file.close()

    pattern_file = open(family_dir + '/pattern.txt', 'w')
    pattern_file.write(pattern)
    pattern_file.close()

    for struct in keep:
        pdb_id = struct[:4]
        chain = struct[4:]

        full_name = pdb_id + '_' + chain
        dir_path = family_dir + '/' + full_name

        try:
            os.makedirs(dir_path)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        segments_file = open(dir_path + '/' + pdb_id + '.seg', 'w')
        segments_file.write('0,{},_,_'.format(chain))
        segments_file.close()

        pdb_file_path = dir_path + '/' + pdb_id + '.pdb.gz'
        shutil.copyfile(pdb_files[pdb_id], pdb_file_path, follow_symlinks=False)

    print('----------')

tmp_dir.cleanup()
