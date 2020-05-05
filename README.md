# RINminer

RINminer is a frequent subgraph mining tool based on the gSpan algorithm with
additional features and optimizations for residue interaction networks (RINs).
Its code has been published under the MIT license.

### Input format

The RINminer input format is compatible with the gSpan reference implementation
which uses a transaction based graph database representation. The input file
contains all graphs represented as individual transactions

A transaction corresponding to a graph starts with a `t`-line. In RINminer this
can be followed by comment lines giving further information about this graph.
These comments with a `#` at the beginning of the line.

This is followed by a vertex list where each vertex is represented
by `v VERTEX_ID LABEL` where the vertex IDs start at 0 and must not contain any
gaps, however their order can be arbitrary. The vertex label is a single
integer and in the case of RINs is used to indicate the residue type.

The vertex list is followed by an edge list of the format
`e VERTEX_ID_1 VERTEX_ID_2 LABEL DISTANCE_LABEL` where the two vertex IDs are
referring to the vertices between which the edge is found. In the case of RINs
the first label is an integer used to indicate the interaction type (contact or
sequence) and the second (optional) label is used to indicate the distance as
integer in milliangstrom. This distance label will be matched approximately and
only be used if the corresponding commandline options are used.

Additionally RINminer supports comment lines before the first transaction. Other
than comments within transactions, these comments at the beginning are copied to
the output file and can be used to store mapping information from vertices to
residues or other meta data.

##### Example input

```
t # 0
# this is the first graph
v 0 19
v 1 10
e 0 1 0 4510
t # 1
# this is the second graph
v 0 18
v 1 19
v 2 10
e 0 1 1 3210
e 1 2 0 5098
e 2 0 0 4490
```

### Output format

The output format follows the same graph representation as the input format with
additional metadata encoded as comments. These metadata lines are:

`# p` followed by a space seprated list of the graph ids of the parents of a
subgraph.

`# t` indicates the type of subgraph: 0 = path, 1 = tree, 2 = cyclic.

`# s` corresponds to the score of the subgraph.

`# m <input graph id> : <exact match?> : <mapping score> : <list of vertex IDs in the input graph>`
describes a single mapping of a subgraph onto an input graph. The first vertex
in the list corresponds to the first vertex in the subgraph.

### Commandline options:
```
rinminer [OPTIONS] <support threshold> <input file> [<output file>]
```

The support threshold indicates the minimum number of database graphs a subgraph
has to be found in.

The input file is in the above described format.

If an output path is specified an output file will be generated. This file will
contain the subgraphs in a format very similar to the input format with
additional metadata in the form of comment lines. The size of this file is
limited and it will only contain the largest (number of edges) subgraphs. It
will contain only complete sets of subgraphs of a certain size. The truncation
ensures that no incomplete sets will be written to the file.

| Option      | Description |
| ----------- | ----------- |
| -t N        | Number of threads. By default RINminer only uses 1 thread. |
| -o N        | Minimum size of subgraphs in output. This helps avoid disk writes for smaller subgraphs. |
| -f N[m,g,t] | Maximum output file size in (mega/giga/tera) bytes. |
| -l database | File path to graph database used to determine the classification ability of each subgraph. This file is in the same format as the input format with the addition of a `# fam FAMILYNAME` comment line in each transaction, where family name is an arbitrary string. |
| -c classes  | Which classes in the significance database are considered positive. The argument is a colon separated string and can contain multiple classes. (Only used in combination with -l) |
| -r          | If specified, the labels will be considered directed. The last edge label will be used to determine the direction of an edge: 0 - undirected, 1 - from V1 to V2, 2 - from V2 to V1. |
| -a N        | Minimum approximate support threshold. Approximate support can be used to still include subgraphs, even if their support value falls below the support threshold, but at a lower rate. This defines the limit up to which support threshold this should be considered. |
| -d N        | Maximum approximate support distance. This sets the maximum support difference between a subgraph and its parent subgraphs for which a subgraph should still be included even if below the support threshold. (Only used in combination with -a) |
| -i N        | Maximum distance difference when matching edges with distance labels. Per mille of the average distance of the compared edges. |
| -s N        | Minimum score increase when determining whether to include a subgraph. |
| -y 1 or 0   | Input file contains interaction type labels (1) or not (0). This can be used to run RINminer with input files only containing distance edge labels. |


### Compiling RINminer

RINminer itself is written in C has no external dependencies. It can be compiled using:

```
cd rinminer
make
```

# Datasets, batch processing and analysis

In addition to the RINminer source code this repository also contains the code
used to generate the datasets used in the paper, as well as for batch processing
and result analysis.

### Dependencies

* [Python 3](https://www.python.org/)
* [Biopython](https://biopython.org/)
* [SciPy](https://www.scipy.org/)
* [python-igraph](https://igraph.org/)
* [Python regex](https://pypi.org/project/regex/)
* [graphviz](https://www.graphviz.org/)
* [RINerator](https://bitbucket.org/doncheva/rinerator/src/master/) (needs to be installed to `rinerator/` within this directory)

RINerator depdends on [probe](https://github.com/rlabduke/probe/) and [reduce](https://github.com/rlabduke/reduce). For the paper we used modified versions of both to fix crashes encountered with some structures. These fixes have been submitted for review.

`reduce` requires `reduce_wwPDB_het_dict.txt` from the source code repository to be copied to `/usr/local/reduce_wwPDB_het_dict.txt` or the `REDUCE_HET_DICT` environment variable to point to that file.

To install RINerator run `git clone https://bitbucket.org/doncheva/rinerator.git` from the root directory of this repository and then modify the probe/reduce paths in `rinerator/Source/get_segments.py`.

Generating datasets from SCOP or PROSITE additionally requires the following files:
* [BLASTDB.tar.gz from PISCES](http://dunbrack.fccc.edu/Guoli/pisces_download.php) (extracted)
* [prosite.dat from PROSITE](ftp://ftp.expasy.org/databases/prosite/prosite.dat)
* [dir.des.scop.txt from SCOP](https://scop.berkeley.edu/downloads/parse/dir.des.scop.1.75.txt)

### Dataset generation

There are three scripts in `generate-dataset/` directory for generating datasets
from different sources (SCOP superfamilies, PROSITE patterns and PDB ID lists).
These scripts create the required files for the batch RIN generation, processing
and analysis. Running any of these scripts with the `-h` option describes their
usage.

### Graph database preparation

The `prepare-database/` directory contains scripts used to generate RINs for the
datasets. The main script is `make-databases.sh` which accepts a dataset
directory as its only parameter.

### Batch processing

Once the graph databases for a dataset have been generated, RINminer can be run
on that dataset. The batch processing for an entire dataset can be started
with the `run-dataset.sh` command from the `run` directory. Running this
command without arguments will provide usage information.

### Analysis

The analysis of the found subgraphs happens in two steps. First the 15 highest
scoring subgraphs without topological duplicates are selected and then in the
second step they are analyzed and visualized. These responsible scripts can be
found in `analyze-results/` and are run automatically from `run-dataset.sh`.

### Example

```
# Generate dataset based on selected SCOP families
cd generate-dataset
./dataset-from-pdb-ids.py example-inputs/pdb-ids/rdrps.txt /tmp/rdrp-dataset
cd ..

# Convert structures into graph databases
cd prepare-database
./make-databases.sh /tmp/rdrp-dataset
cd ..

# Run RINminer on the dataset
cd run
./run-dataset.sh /tmp/rdrp-dataset /tmp/rdrp-results /tmp/rdrp-report
```

After the run command has finished a HTML report showing the highest scoring
subgraphs with additional information under `/tmp/rdrp-report/report.html`.
