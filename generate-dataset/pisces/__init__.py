import os, sys

# Uses BLASTDB.tar.gz from http://dunbrack.fccc.edu/Guoli/pisces_download.php

class Pisces:
    def __init__(self, blastdb_path, max_res=3.0, max_r_value=1.0, min_length=50, aux_length=25, max_identity=30):
        self.max_identity = max_identity
        self.aux_length = aux_length

        # prepare list of all chain IDs in the DB sorted by "quality"
        with open(blastdb_path + '/resolution.dat', 'r', encoding='ascii') as res_file:
            to_sort = []
            for line in res_file:
                pdb_chain_id, method, resolution_, r_value_, free_r_value_, ca_only, length_ = line.strip().split()

                if method != "XRAY":
                    continue

                if ca_only != "no":
                    continue

                resolution = 100.0
                try:
                    resolution = float(resolution_)
                except:
                    pass
                if resolution > max_res:
                    continue

                r_value = 1.0
                try:
                    r_value = float(r_value_)
                except:
                    pass
                if r_value > max_r_value:
                    continue

                length = int(length_)
                if length < min_length:
                    continue

                to_sort.append((resolution, r_value, -length, pdb_chain_id))

            to_sort.sort()
            self.lengths = { entry[3]: -entry[2] for entry in to_sort }
            self.pdb_sort_keys = { entry[3]: i for (i, entry) in enumerate(to_sort) }

        # Prepare redundancy mapping
        replaced = {}
        with open(blastdb_path + '/MakeNonRedundant.log.pdb', 'r', encoding='ascii') as redundancy_file:
            for line in redundancy_file:
                if line[0] == '#':
                    continue
                orig, replacement = line.strip().split(' by ')
                # The replacement might not be in our set of chains, but is still needed. So don't filter it here
                if orig in self.pdb_sort_keys:
                    replaced[orig] = replacement

        chain_ids_non_redundant = (self.pdb_sort_keys.keys() - replaced.keys()) | set(replaced.values())
        self.struct_chain_to_idx = { struct_chain: i for (i, struct_chain) in enumerate(chain_ids_non_redundant) }
        for (replaced, replacement) in replaced.items():
            self.struct_chain_to_idx[replaced] = self.struct_chain_to_idx[replacement]

        # Prepare identity matrix
        self.identical = [ set() for i in range(len(chain_ids_non_redundant)) ]

        with open(blastdb_path + '/pdbaa.align', 'r', encoding='ascii') as align_file:
            for line in align_file:
                hit, query, identity, _ = line.split()

                if hit not in chain_ids_non_redundant \
                or query not in chain_ids_non_redundant \
                or int(identity, 10) <= self.max_identity:
                    continue

                h_idx = self.struct_chain_to_idx[hit]
                q_idx = self.struct_chain_to_idx[query]
                self.identical[q_idx].add(h_idx)
                self.identical[h_idx].add(q_idx)

    # Sort instances by their structure properties and structure id
    def _sort_instances(self, instance):
        struct_id, chains, aux_chains = instance
        sort_keys_chains = [ self.pdb_sort_keys[struct_id + chain] for chain in chains ]
        sort_keys_chains.sort()
        return sort_keys_chains

    def filter(self, instances):
        single_chain = False
        if isinstance(next(iter(instances)), str):
            single_chain = True
            instances = [ (x[:4], (x[4])) for x in instances ]
        initial_size = len(instances)
        filtered_instances = self._filter_included(instances)
        after_experimental_filtering = len(filtered_instances)
        filtered_instances = self._filter_identity(filtered_instances)
        after_identity_filtering = len(filtered_instances)
        if single_chain:
            filtered_instances = [ x[0] + x[1][0] for x in filtered_instances ]
        stats = (initial_size - after_experimental_filtering, after_experimental_filtering - after_identity_filtering)
        return (filtered_instances, stats)

    # This filters out chains from entries that are not included in the PISCES
    # database.
    def _filter_included(self, instances):
        filtered_instances = []
        for struct_id, chains in instances:
            filtered_chains = [ chain for chain in chains if struct_id + chain in self.pdb_sort_keys and self.lengths[struct_id + chain] >= self.aux_length ]
            if len(filtered_chains) > 0:
                # Small chains that are not considered for identity
                aux_chains = [ chain for chain in chains if struct_id + chain not in self.lengths or self.lengths[struct_id + chain] < self.aux_length ]
                filtered_instances.append((struct_id, tuple(filtered_chains), tuple(aux_chains)))
        return filtered_instances

    # This filters instances consisting of one or more chains such that the resulting
    # instances contain no chains with a sequence identity above the threshold in
    # any other instance than itself. Chains within an instance however are
    # allowed to have a higher identity to account for homo k-mers. Multiple instances
    # interacting with exactly the same chains are also *not* filtered.
    def _filter_identity(self, instances):
        chain_ids_flat = set()
        for struct_id, chains, aux_chains in instances:
            chain_ids_flat.update([ self.struct_chain_to_idx[struct_id + chain] for chain in chains ])
        instances.sort(key=self._sort_instances)

        # Determine the set of structures to keep while removing redundant structures
        kept_instances = set()
        skip = set()

        for instance in instances:
            struct_id, chains, aux_chains = instance

            all_ids = { self.struct_chain_to_idx[struct_id + chain] for chain in chains }

            if not all_ids.isdisjoint(skip):
                continue

            kept_instances.add(instance)

            # This relies on the identity matrix containing all combinations with a significant identity
            skip.update(all_ids)
            skip.update(set().union(*[self.identical[i] for i in all_ids]) & chain_ids_flat)

        result = []
        for struct_id, chains, aux_chains in kept_instances:
            chains = tuple(sorted(chains + aux_chains))
            result.append((struct_id, chains))

        return result
