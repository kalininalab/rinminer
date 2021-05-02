#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <pthread.h>
#include <unistd.h>
#include <locale.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

typedef int vertex_label_t;

#define VERTEX_LABEL_MAX (INT_MAX)

typedef struct edge_label_t {
    int interaction_type;
    int distance;
} edge_label_t;

#define EDGE_LABEL_MAX ((edge_label_t) { INT_MAX, INT_MAX })

typedef uint16_t vertex_t;

#define VERTEX_INVALID ((vertex_t) UINT16_MAX)
#define VERTEX_MAX ((vertex_t) VERTEX_INVALID - 1)

typedef enum direction_t direction_t;
enum direction_t {
    UNDIRECTED = 0,
    FORWARD    = 1,
    BACKWARD   = 2,

    DIRECTION_COUNT
};

typedef struct dfs_code_t {
    vertex_t from;
    vertex_t to;
    vertex_label_t from_label;
    vertex_label_t to_label;
    edge_label_t edge_label;
    direction_t direction;
} dfs_code_t;

#define DFS_CODE_MAX ((dfs_code_t) { VERTEX_MAX, VERTEX_MAX - 1, VERTEX_LABEL_MAX, VERTEX_LABEL_MAX, EDGE_LABEL_MAX, BACKWARD })

typedef struct neighbor_t {
    vertex_t vertex_id;
    edge_label_t edge_label;
    direction_t direction;
} neighbor_t;

typedef enum graph_type_t graph_type_t;
enum graph_type_t {
    PATH,
    TREE,
    CYCLE,

    GRAPH_TYPE_COUNT
};

typedef struct mapping_info_t mapping_info_t;
struct mapping_info_t {
    int mapping_score;
    float mismatch_score;
};

typedef struct mapped_vertex_t mapped_vertex_t;
struct mapped_vertex_t {
    mapped_vertex_t *parent;
    volatile unsigned int reference_count;
    vertex_t database_vertex_id;
};

typedef struct mapping_t mapping_t;
struct mapping_t {
    mapping_t *next;
    mapped_vertex_t *mapped_vertex;
    int exact;
    mapping_info_t info;
};

typedef struct db_mapping_t db_mapping_t;
struct db_mapping_t {
    mapping_t *mappings;
    unsigned int db_idx;
    vertex_t last_mapped_vertex_id;
    int exact;
};

typedef struct embedding_store_t embedding_store_t;
struct embedding_store_t {
    mapped_vertex_t *vertex_free_list;
    mapped_vertex_t **vertex_blocks;
    mapping_t *mapping_free_list;
    mapping_t **mapping_blocks;
    unsigned int vertex_num_blocks;
    unsigned int vertex_remaining_in_block;
    unsigned int vertex_current_block_index;
    unsigned int mapping_num_blocks;
    unsigned int mapping_remaining_in_block;
    unsigned int mapping_current_block_index;
};

#define EMBEDDING_STORE_BLOCK_SIZE (1024 * 1024)

typedef struct subgraph_t subgraph_t;
struct subgraph_t {
    dfs_code_t *dfs_code;
    vertex_label_t *vertex_labels;
    vertex_t *rightmost_path;
    db_mapping_t *db_mappings;
    db_mapping_t *sig_db_mappings;
    int *supported_by;

    subgraph_t **parents;
    unsigned int *parent_ids;

    unsigned int id;
    uint32_t hash;
    unsigned int support;
    int score;

    unsigned int num_mappings;
    unsigned int rightmost_path_len;
    unsigned int dfs_code_len;
    unsigned int num_parents;

    unsigned int sig_support;
    unsigned int sig_good_support;

    vertex_t rightmost_vertex;

    graph_type_t graph_type;
    subgraph_t *next;
};

typedef struct extension_t extension_t;
struct extension_t {
    extension_t *next;
    int *supported_by;
    unsigned int support;
    dfs_code_t code;
};

typedef struct extension_store_t extension_store_t;
struct extension_store_t {
    extension_t *extensions;
    extension_t **hash_map;
    unsigned int *max_extensions_for_length;
    unsigned int num_extensions;
    unsigned int hash_map_size;
    unsigned int extensions_size;
};

#define EXTENSION_STORE_HASH_MAP_SIZE (2048)

typedef struct rem_edge_store_t rem_edge_store_t;
struct rem_edge_store_t {
    unsigned int *cycles_closed_by;
    unsigned int *cycles_open_at;
    unsigned int *edge_count;
};

typedef struct database_entry_t database_entry_t;
struct database_entry_t {
    vertex_label_t *vertex_labels;
    unsigned int *neighbors_len;
    neighbor_t **neighbors;
    unsigned int id;
    int classification;
    unsigned int num_vertices;
    unsigned int num_edges;
};

typedef struct database_t database_t;
struct database_t {
    database_entry_t *entries;
    unsigned int num_entries;
    unsigned int num_classified_good;
    unsigned int num_label_types;
};

typedef struct mismatch_matrix_t mismatch_matrix_t;
struct mismatch_matrix_t {
    float *matrix;
    unsigned int num_entries;
};

typedef struct config_t config_t;
struct config_t {
    unsigned int direct_min_support;
    unsigned int approx_max_support_difference;
    unsigned int approx_min_support;
    unsigned int max_distance_difference;
    unsigned int min_score_difference;
    float max_mismatch_score;
    mismatch_matrix_t mismatch_matrix;
};

typedef struct output_t output_t;
struct output_t {
    char *buffer;
    int fd;
    off_t buffer_size;
    off_t cursor;
    off_t bytes_written;
    off_t *start_barriers;
    unsigned int starts_len;
    pthread_mutex_t mutex;
    unsigned int min_size;
    unsigned int header_bytes_used;
    int use_distances;
    int use_directions;
};

typedef struct subgraph_store_t subgraph_store_t;
struct subgraph_store_t {
    subgraph_t *subgraphs;
    subgraph_t **hash_map;
    unsigned int num_subgraphs;
    unsigned int hash_map_size;
};

typedef struct thread_data_t thread_data_t;
struct thread_data_t {
    subgraph_store_t *in_subgraphs;
    volatile unsigned int *current_item;
    config_t *config;
    database_t *database;
    database_t *sig_database;
    embedding_store_t embedding_store;
    embedding_store_t sig_embedding_store;
    embedding_store_t canon_embedding_store;
    extension_store_t extension_store;
    extension_store_t *shared_extension_store;
    rem_edge_store_t rem_edge_store;
    output_t *output;
    char *output_buffer;
    unsigned int output_buffer_size;
    subgraph_store_t out_subgraphs;
    unsigned int item_count;
};

static inline void free_extension (extension_t *extension)
{
    free (extension->supported_by);
}

static inline void free_extension_store (extension_store_t *extension_store)
{
    free (extension_store->hash_map);
    free (extension_store->extensions);
}

static inline void free_rem_edge_store (rem_edge_store_t *rem_edge_store)
{
    free (rem_edge_store->cycles_open_at);
    free (rem_edge_store->cycles_closed_by);
    free (rem_edge_store->edge_count);
}

static inline void free_subgraph_labels (subgraph_t *subgraph)
{
    free (subgraph->vertex_labels);
    subgraph->vertex_labels = NULL;
}

static inline void free_subgraph_rightmost_path (subgraph_t *subgraph)
{
    free (subgraph->rightmost_path);
    subgraph->rightmost_path = NULL;
}

static inline void free_subgraph_parents (subgraph_t *subgraph)
{
    free (subgraph->parents);
    subgraph->parents = NULL;
}

static inline void free_subgraph_parent_ids (subgraph_t *subgraph)
{
    free (subgraph->parent_ids);
    subgraph->parent_ids = NULL;
}

static inline void unref_mappings (embedding_store_t *embedding_store,
                                   db_mapping_t      *db_mapping);

static inline void free_subgraph_mappings (embedding_store_t *embedding_store,
                                           subgraph_t        *subgraph)
{
    for (unsigned int db_iterator = 0; db_iterator < subgraph->support; db_iterator++) {
        db_mapping_t *db_mapping = subgraph->db_mappings + db_iterator;
        unref_mappings (embedding_store, db_mapping);
    }
}

static inline void free_subgraph_sig_mappings (embedding_store_t *embedding_store,
                                               subgraph_t        *subgraph)
{
    for (unsigned int sig_db_iterator = 0; sig_db_iterator < subgraph->sig_support; sig_db_iterator++) {
        db_mapping_t *sig_db_mapping = subgraph->sig_db_mappings + sig_db_iterator;
        unref_mappings (embedding_store, sig_db_mapping);
    }
}

static inline void free_subgraph (subgraph_t *subgraph)
{
    free (subgraph->vertex_labels);
    free (subgraph->rightmost_path);
    free (subgraph->dfs_code);
    free (subgraph->parents);
    free (subgraph->parent_ids);
    free (subgraph->db_mappings);
    free (subgraph->supported_by);
    free (subgraph->sig_db_mappings);
    // free_subgraph_mappings needs to be called manually where appropriate
}

static inline void free_embedding_store (embedding_store_t *embedding_store)
{
#ifndef NDEBUG
    // Check against leaks - does not work with the current approach, because
    // free lists can contain entries allocated by other threads.
    /*unsigned int f = 0;
    mapped_vertex_t *v = embedding_store->vertex_free_list;
    while (v) {
        f++;
        v = v->parent;
    }
    unsigned int u = embedding_store->vertex_remaining_in_block;
    unsigned int r = (embedding_store->vertex_num_blocks - embedding_store->vertex_current_block_index - 1) * EMBEDDING_STORE_BLOCK_SIZE;
    unsigned int t = f + u + r;
    fprintf (stderr, "v: %d/%d (%d + %d + %d)\n", t, embedding_store->vertex_num_blocks * EMBEDDING_STORE_BLOCK_SIZE, f, u, r);
    assert (t == embedding_store->vertex_num_blocks * EMBEDDING_STORE_BLOCK_SIZE);

    f = 0;
    mapping_t *m = embedding_store->mapping_free_list;
    while (m) {
        f++;
        m = m->next;
    }
    u = embedding_store->mapping_remaining_in_block;
    r = (embedding_store->mapping_num_blocks - embedding_store->mapping_current_block_index - 1) * EMBEDDING_STORE_BLOCK_SIZE;
    t = f + u + r;
    fprintf (stderr, "m: %d/%d (%d + %d + %d)\n", t, embedding_store->mapping_num_blocks * EMBEDDING_STORE_BLOCK_SIZE, f, u, r);
    assert (t == embedding_store->mapping_num_blocks * EMBEDDING_STORE_BLOCK_SIZE);
    */
#endif

    for (unsigned int i = 0; i < embedding_store->vertex_num_blocks; i++) {
        free (embedding_store->vertex_blocks[i]);
    }
    free (embedding_store->vertex_blocks);

    for (unsigned int i = 0; i < embedding_store->mapping_num_blocks; i++) {
        free (embedding_store->mapping_blocks[i]);
    }
    free (embedding_store->mapping_blocks);
}

static inline void free_canonical_database_entry (database_entry_t *database_entry)
{
    // The first entry in neighbors points to the address of the memory allocated
    // for all vertices.
    free (database_entry->neighbors[0]);
    free (database_entry->neighbors);
    free (database_entry->neighbors_len);
}

static inline void free_database_entry (database_entry_t *database_entry)
{
    free (database_entry->vertex_labels);
    free_canonical_database_entry (database_entry);
}

static inline void free_database (database_t *database)
{
    for (unsigned int i = 0; i < database->num_entries; i++) {
        free_database_entry (&database->entries[i]);
    }
    free (database->entries);
}

static inline void free_subgraph_store (subgraph_store_t *subgraph_store)
{
    for (unsigned int i = 0; i < subgraph_store->num_subgraphs; i++) {
        free_subgraph (&subgraph_store->subgraphs[i]);
    }
    free (subgraph_store->subgraphs);
    free (subgraph_store->hash_map);
}

static inline direction_t invert_direction (direction_t direction)
{
    if (direction == FORWARD) {
        return BACKWARD;
    }
    if (direction == BACKWARD) {
        return FORWARD;
    }
    return UNDIRECTED;
}

static inline float mismatch_score_between (mismatch_matrix_t  matrix,
                                            vertex_label_t     first_label,
                                            vertex_label_t     second_label);

// Returns a negative number if a is smaller, a positive number if b is smaller
// and 0 otherwise
static inline int vertex_label_compare (vertex_label_t *a,
                                        vertex_label_t *b)
{
    return *a - *b;
}

// Returns 1 if two labels are compatible and can be matched in a non-exact
// isomorphism and 0 otherwise
static inline unsigned int vertex_label_match (config_t       *config,
                                               mapping_info_t *mapping_info,
                                               vertex_label_t *a,
                                               vertex_label_t *b,
                                               mapping_info_t *new_mapping_info)
{
    float mismatch_score = mapping_info->mismatch_score
                           + mismatch_score_between (config->mismatch_matrix, *a, *b);
    if (new_mapping_info){
        new_mapping_info->mismatch_score = mismatch_score;
    }
    return mismatch_score <= config->max_mismatch_score;
}

// Returns a negative number if a is smaller, a positive number if b is smaller
// and 0 otherwise
static inline int edge_label_compare (edge_label_t *a,
                                      edge_label_t *b)
{
    if (a->interaction_type < b->interaction_type) {
        return -1;
    }

    if (a->interaction_type > b->interaction_type) {
        return 1;
    }

    if (a->distance < b->distance) {
        return -1;
    }

    if (a->distance > b->distance) {
        return 1;
    }

    return 0;
}

// Returns 1 if two labels are compatible and can be matched in a non-exact
// isomorphism and 0 otherwise
static inline unsigned int edge_label_match (config_t       *config,
                                             mapping_info_t *mapping_info,
                                             edge_label_t   *a,
                                             edge_label_t   *b,
                                             mapping_info_t *new_mapping_info)
{
    if (a->interaction_type != b->interaction_type) {
        return 0;
    }

    int normalize = MAX ((a->distance + b->distance) / 2, 1);
    int distance_difference = (abs (a->distance - b->distance) * 1000) / normalize;

    if (distance_difference <= config->max_distance_difference) {
        if (new_mapping_info) {
            int edge_score = 1000 - distance_difference;
            new_mapping_info->mapping_score = mapping_info->mapping_score + edge_score;
        }
        return 1;
    }

    return 0;
}

// Jenkins hash function
// http://www.burtleburtle.net/bob/hash/doobs.html
static inline uint32_t dfs_code_hash (dfs_code_t   *code,
                                      unsigned int  length)
{
    uint32_t hash = 0;

    for (unsigned int i = 0; i < length; i++) {
        hash += code[i].from;
        hash += (hash << 10);
        hash ^= (hash >> 6);

        hash += code[i].from_label;
        hash += (hash << 10);
        hash ^= (hash >> 6);

        hash += code[i].to;
        hash += (hash << 10);
        hash ^= (hash >> 6);

        hash += code[i].to_label;
        hash += (hash << 10);
        hash ^= (hash >> 6);

        hash += code[i].edge_label.interaction_type;
        hash += (hash << 10);
        hash ^= (hash >> 6);

        hash += code[i].edge_label.distance;
        hash += (hash << 10);
        hash ^= (hash >> 6);

        hash += code[i].direction;
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }

    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);

    return hash;
}

static inline void calc_subgraph_hash (subgraph_t *subgraph)
{
    uint32_t hash = dfs_code_hash (subgraph->dfs_code, subgraph->dfs_code_len);
    subgraph->hash = hash;
}

static inline unsigned int labels_to_index (mismatch_matrix_t matrix,
                                            vertex_label_t    vertex_label_a,
                                            vertex_label_t    vertex_label_b)
{
    unsigned int index = vertex_label_a * matrix.num_entries + vertex_label_b;
    return index;
}

mismatch_matrix_t generate_default_mismatch_matrix (unsigned int num_vertex_label_types)
{
    mismatch_matrix_t mismatch_matrix = {};
    unsigned int size = num_vertex_label_types * num_vertex_label_types;
    mismatch_matrix.num_entries = num_vertex_label_types;
    mismatch_matrix.matrix = malloc (size * sizeof(mismatch_matrix.matrix[0]));
    for (unsigned int i = 0; i < num_vertex_label_types; i++) {
        for (unsigned int j = 0; j < num_vertex_label_types; j++) {
            unsigned int index = labels_to_index (mismatch_matrix, i, j);
            if (i == j) {
                mismatch_matrix.matrix[index] = 0.0;
            }
            else {
                mismatch_matrix.matrix[index] = 0.1;
            }
        }
    }
    return mismatch_matrix;
}

mismatch_matrix_t parse_mismatch_matrix (unsigned int  num_vertex_label_types,
                                         char         *file_path)
{
    FILE* f = fopen (file_path, "r");
    if (!f) {
        fprintf (stderr, "Could not open file: %s\n", file_path);
        exit (1);
    }

    mismatch_matrix_t mismatch_matrix = {};
    unsigned int size = num_vertex_label_types * num_vertex_label_types;
    mismatch_matrix.num_entries = num_vertex_label_types;
    mismatch_matrix.matrix = malloc (size * sizeof(mismatch_matrix.matrix[0]));
    unsigned int entries = 0;

    while (!feof (f)) {
        char c = fgetc (f);
        if ((c >= '0' && c <= '9') || c == '.') {
            ungetc (c, f);
            fscanf (f,"%f", &mismatch_matrix.matrix[entries++]);
        }
        if (entries > size) {
            size *= 2;
            mismatch_matrix.matrix = realloc (mismatch_matrix.matrix, size * sizeof(mismatch_matrix.matrix[0]));
        }
    }

    if (entries < size) {
        fprintf (stderr, "Database contains labels larger than in mismatch scoring matrix: %d\n", num_vertex_label_types - 1);
        exit (1);
    }

    return mismatch_matrix;
}

static inline void free_mismatch_matrix (mismatch_matrix_t mismatch_matrix) {
    free (mismatch_matrix.matrix);
}

static inline float mismatch_score_between (mismatch_matrix_t matrix,
                                            vertex_label_t    first_label,
                                            vertex_label_t    second_label)
{
    unsigned int index = labels_to_index (matrix, first_label, second_label);
    return matrix.matrix[index];
}

static inline int rightmost_path_search (vertex_t     *list,
                                         unsigned int  length,
                                         vertex_t      target)
{
    for (unsigned int i = 0; i < length; i++) {
        if (list[i] == target) {
            return i;
        }
    }
    return -1;
}

static inline int code_matches (dfs_code_t *a,
                                dfs_code_t *b)
{
    return a->from == b->from
           && a->to == b->to
           && a->direction == b->direction
           && vertex_label_compare (&a->from_label, &b->from_label) == 0
           && vertex_label_compare (&a->to_label, &b->to_label) == 0
           && edge_label_compare (&a->edge_label, &b->edge_label) == 0;
}

static inline int compare_code (dfs_code_t *a,
                                dfs_code_t *b)
{
    // condition 1:
    if (a->from < a->to && b->from < b->to) {
        if (a->to < b->to || (a->to == b->to && a->from > b->from)) {
            return -1;
        }
    }

    // condition 2:
    if (a->to < a->from && b->to < b->from) {
        if (a->from < b->from || (a->from == b->from && a->to < b->to)) {
            return -1;
        }
    }

    // condition 3:
    if (a->from < a->to && b->to < b->from) {
        if (a->to <= b->from) {
            return -1;
        }
    }

    // condition 4:
    if (a->to < a->from && b->from < b->to) {
        if (a->from < b->to) {
            return -1;
        }
    }

    if (a->from == b->from && a->to == b->to) {
        int from_comparison = vertex_label_compare(&a->from_label, &b->from_label);
        int to_comparison = vertex_label_compare(&a->to_label, &b->to_label);
        int edge_comparison = edge_label_compare(&a->edge_label, &b->edge_label);

        if (from_comparison < 0) {
            return -1;
        }
        else if (from_comparison == 0) {
            if (to_comparison < 0) {
                return -1;
            }
            else if (to_comparison == 0) {
                if (edge_comparison < 0) {
                    return -1;
                }
                else if (edge_comparison == 0) {
                    if (a->direction < b->direction) {
                        return -1;
                    }
                    else if (a->direction == b->direction) {
                        return 0;
                    }
                    else {
                        return 1;
                    }
                }
                else {
                    return 1;
                }
            }
            else {
                return 1;
            }
        }
        else {
            return 1;
        }
    }

    return 1;
}

static inline int contains_edge (subgraph_t   *subgraph,
                                 vertex_t      vertex_a,
                                 vertex_t      vertex_b,
                                 edge_label_t *edge_label,
                                 direction_t   direction)
{

    for (unsigned int i = subgraph->dfs_code_len; i-- > 0;) {
        dfs_code_t *cur = &subgraph->dfs_code[i];
        if (((cur->from == vertex_a && cur->to == vertex_b && cur->direction == direction)
             || (cur->from == vertex_b && cur->to == vertex_a && cur->direction == invert_direction(direction)))
            && edge_label_compare (&cur->edge_label, edge_label) == 0) {
            return 1;
        }
    }
    return 0;
}

static inline void init_extension_store (extension_store_t *extension_store,
                                         unsigned int       num_initial_extensions,
                                         unsigned int      *max_extensions_for_length)
{
    extension_store->hash_map_size = EXTENSION_STORE_HASH_MAP_SIZE;
    extension_store->hash_map = malloc (extension_store->hash_map_size * sizeof(extension_store->hash_map[0]));
    extension_store->extensions_size = num_initial_extensions;
    extension_store->extensions = malloc (num_initial_extensions * sizeof(extension_store->extensions[0]));
    extension_store->max_extensions_for_length = max_extensions_for_length;
}

static inline void prepare_extension_store (extension_store_t *extension_store,
                                            unsigned int       num_mappings,
                                            unsigned int       rightmost_path_length)
{
    memset (extension_store->hash_map, 0, extension_store->hash_map_size * sizeof(extension_store->hash_map[0]));
    extension_store->num_extensions = 0;
    unsigned int max_extensions = num_mappings * extension_store->max_extensions_for_length[rightmost_path_length];
    if (max_extensions > extension_store->extensions_size) {
        extension_store->extensions_size = max_extensions;
        free (extension_store->extensions);
        extension_store->extensions = malloc (max_extensions * sizeof(extension_store->extensions[0]));
    }
}

static inline void init_embedding_store (embedding_store_t *embedding_store)
{
    embedding_store->vertex_num_blocks = 1;
    embedding_store->vertex_blocks = malloc (sizeof(embedding_store->vertex_blocks[0]));
    embedding_store->vertex_blocks[0] = malloc (EMBEDDING_STORE_BLOCK_SIZE * sizeof(embedding_store->vertex_blocks[0][0]));
    embedding_store->vertex_remaining_in_block = EMBEDDING_STORE_BLOCK_SIZE;
    embedding_store->vertex_current_block_index = 0;
    embedding_store->vertex_free_list = NULL;

    embedding_store->mapping_num_blocks = 1;
    embedding_store->mapping_blocks = malloc (sizeof(embedding_store->mapping_blocks[0]));
    embedding_store->mapping_blocks[0] = malloc (EMBEDDING_STORE_BLOCK_SIZE * sizeof(embedding_store->mapping_blocks[0][0]));
    embedding_store->mapping_remaining_in_block = EMBEDDING_STORE_BLOCK_SIZE;
    embedding_store->mapping_current_block_index = 0;
    embedding_store->mapping_free_list = NULL;
}


static inline void clear_embedding_store (embedding_store_t *embedding_store)
{
    embedding_store->vertex_remaining_in_block = EMBEDDING_STORE_BLOCK_SIZE;
    embedding_store->vertex_current_block_index = 0;
    embedding_store->vertex_free_list = NULL;

    embedding_store->mapping_remaining_in_block = EMBEDDING_STORE_BLOCK_SIZE;
    embedding_store->mapping_current_block_index = 0;
    embedding_store->mapping_free_list = NULL;
}

static inline void init_rem_edge_store (rem_edge_store_t *rem_edge_store,
                                        unsigned int      max_vertices)
{
    rem_edge_store->cycles_open_at = malloc (sizeof(rem_edge_store->cycles_open_at[0]) * max_vertices);
    rem_edge_store->cycles_closed_by = malloc (sizeof(rem_edge_store->cycles_closed_by[0]) * max_vertices);
    rem_edge_store->edge_count = malloc (sizeof(rem_edge_store->edge_count[0]) * max_vertices);
}

mapped_vertex_t *request_mapped_vertex (embedding_store_t *embedding_store,
                                        int                use_refcounts)
{
    mapped_vertex_t *result;
    if (use_refcounts && embedding_store->vertex_free_list != NULL) {
        result = embedding_store->vertex_free_list;
        embedding_store->vertex_free_list = result->parent;
    }
    else {
        if (embedding_store->vertex_remaining_in_block == 0) {
            embedding_store->vertex_current_block_index++;
            embedding_store->vertex_remaining_in_block = EMBEDDING_STORE_BLOCK_SIZE;

            // If all blocks are full, allocate a new block
            if (embedding_store->vertex_current_block_index == embedding_store->vertex_num_blocks){
                embedding_store->vertex_num_blocks++;
                embedding_store->vertex_blocks = realloc (embedding_store->vertex_blocks,
                                                          embedding_store->vertex_num_blocks * sizeof(embedding_store->vertex_blocks[0]));
                embedding_store->vertex_blocks[embedding_store->vertex_current_block_index] = malloc (EMBEDDING_STORE_BLOCK_SIZE * sizeof(embedding_store->vertex_blocks[0][0]));
            }
        }
        mapped_vertex_t *current_block = embedding_store->vertex_blocks[embedding_store->vertex_current_block_index];
        result = &current_block[EMBEDDING_STORE_BLOCK_SIZE - embedding_store->vertex_remaining_in_block];
        embedding_store->vertex_remaining_in_block--;
    }

    return result;
}

mapping_t *request_mapping (embedding_store_t *embedding_store,
                            int                use_refcounts)
{
    mapping_t *result;
    if (use_refcounts && embedding_store->mapping_free_list != NULL) {
        result = embedding_store->mapping_free_list;
        embedding_store->mapping_free_list = result->next;
    }
    else {
        if (embedding_store->mapping_remaining_in_block == 0) {
            embedding_store->mapping_current_block_index++;
            embedding_store->mapping_remaining_in_block = EMBEDDING_STORE_BLOCK_SIZE;

            // If all blocks are full, allocate a new block
            if (embedding_store->mapping_current_block_index == embedding_store->mapping_num_blocks){
                embedding_store->mapping_num_blocks++;
                embedding_store->mapping_blocks = realloc (embedding_store->mapping_blocks,
                                                          embedding_store->mapping_num_blocks * sizeof(embedding_store->mapping_blocks[0]));
                embedding_store->mapping_blocks[embedding_store->mapping_current_block_index] = malloc (EMBEDDING_STORE_BLOCK_SIZE * sizeof(embedding_store->mapping_blocks[0][0]));
            }
        }
        mapping_t *current_block = embedding_store->mapping_blocks[embedding_store->mapping_current_block_index];
        result = &current_block[EMBEDDING_STORE_BLOCK_SIZE - embedding_store->mapping_remaining_in_block];
        embedding_store->mapping_remaining_in_block--;
    }

    return result;
}

static inline void unref_mapping_helper (mapping_t        *mapping,
                                         mapped_vertex_t **out_free_list_start,
                                         mapped_vertex_t **out_free_list_end)

{
    mapped_vertex_t *mapped_vertex = mapping->mapped_vertex;
    unsigned int refs = __sync_sub_and_fetch(&mapped_vertex->reference_count, 1);

    mapped_vertex_t *free_list_start = NULL;
    mapped_vertex_t *free_list_end = NULL;

    if (refs == 0) {
        free_list_start = mapped_vertex;
        free_list_end = mapped_vertex;


        while (free_list_end->parent) {
            unsigned int parent_refs = __sync_sub_and_fetch(&free_list_end->parent->reference_count, 1);
            if (parent_refs > 0) {
                break;
            }
            free_list_end = free_list_end->parent;
        }
    }

    *out_free_list_start = free_list_start;
    *out_free_list_end = free_list_end;
}

static inline void unref_mappings (embedding_store_t *embedding_store,
                                   db_mapping_t      *db_mappings)
{
    mapping_t *mapping_start = db_mappings->mappings;
    mapping_t *mapping_end = NULL;
    mapped_vertex_t *vertex_free_list_start = NULL;
    mapped_vertex_t *vertex_free_list_end = NULL;

    mapping_t *mapping = mapping_start;
    while (mapping) {
        mapping_end = mapping;
        mapped_vertex_t *new_free_list_start;
        mapped_vertex_t *new_free_list_end;
        unref_mapping_helper (mapping, &new_free_list_start, &new_free_list_end);
        if (new_free_list_start) {
            new_free_list_end->parent = vertex_free_list_start;
            vertex_free_list_start = new_free_list_start;
            if (!vertex_free_list_end) {
                vertex_free_list_end = new_free_list_end;
            }
        }
        mapping = mapping->next;
    }

    if (vertex_free_list_start) {
        vertex_free_list_end->parent = embedding_store->vertex_free_list;
        embedding_store->vertex_free_list = vertex_free_list_start;
    }

    mapping_end->next = embedding_store->mapping_free_list;
    embedding_store->mapping_free_list = mapping_start;
}

static inline mapped_vertex_t *add_mapped_vertex (embedding_store_t *embedding_store,
                                                  mapped_vertex_t   *parent,
                                                  vertex_t           db_vertex_id,
                                                  int                use_refcounts)
{
    mapped_vertex_t *vertex = request_mapped_vertex (embedding_store, use_refcounts);

    vertex->parent = parent;
    vertex->database_vertex_id = db_vertex_id;

    if (use_refcounts) {
        if (parent) {
            vertex->reference_count = 1;
            __sync_add_and_fetch(&parent->reference_count, 1);
        }
        else {
            // The first mapped vertex is coupled to its children
            vertex->reference_count = 0;
        }
    }

    return vertex;
}

static inline mapping_t *add_mapping_helper (embedding_store_t *embedding_store,
                                             db_mapping_t      *db_mappings,
                                             int                use_refcounts)
{
    mapping_t *mapping = request_mapping (embedding_store, use_refcounts);
    mapping->next = db_mappings->mappings;
    db_mappings->mappings = mapping;

    return mapping;
}

static inline mapping_t *add_mapping (embedding_store_t *embedding_store,
                                      db_mapping_t      *db_mappings,
                                      mapped_vertex_t   *parent,
                                      vertex_t           db_vertex_id,
                                      int                exact,
                                      mapping_info_t    *mapping_info,
                                      int                use_refcounts)
{
    mapped_vertex_t *vertex = add_mapped_vertex (embedding_store, parent, db_vertex_id, use_refcounts);
    mapping_t *mapping = add_mapping_helper (embedding_store, db_mappings, use_refcounts);

    mapping->mapped_vertex = vertex;
    mapping->exact = exact;
    mapping->info = *mapping_info;

    return mapping;
}

static inline mapping_t * copy_map_with_new_mapping_info (embedding_store_t *embedding_store,
                                                          mapping_t         *source_mapping,
                                                          db_mapping_t      *target_db_mappings,
                                                          int                exact,
                                                          mapping_info_t    *mapping_info,
                                                          int                use_refcounts)
{
    mapped_vertex_t *vertex = source_mapping->mapped_vertex;
    mapping_t *mapping = add_mapping_helper (embedding_store, target_db_mappings, use_refcounts);

    mapping->mapped_vertex = vertex;
    mapping->exact = exact;
    mapping->info = *mapping_info;

    if (use_refcounts) {
        __sync_add_and_fetch(&vertex->reference_count, 1);
    }

    return mapping;
}

void calc_subgraph_score (config_t   *config,
                          subgraph_t *subgraph)
{
    int score = 0;

    unsigned int num_considered_mappings = config->approx_min_support;

    int *highest_mapping_scores = malloc (num_considered_mappings * sizeof(highest_mapping_scores[0]));

    for (unsigned int i = 0; i < num_considered_mappings; i++) {
        highest_mapping_scores[i] = 0;
    }

    int lowest_high_score_index = 0;
    int lowest_high_score = highest_mapping_scores[lowest_high_score_index];

    for (unsigned int db_iterator = 0; db_iterator < subgraph->support; db_iterator++) {
        db_mapping_t *db_mapping = subgraph->db_mappings + db_iterator;
        int highest_mapping_score = 0;

        // Determine highest score of all mappings into this db graph
        for (mapping_t *mapping = db_mapping->mappings;
             mapping;
             mapping = mapping->next)
        {
            mapping_info_t *mapping_info = &mapping->info;
            if (mapping_info->mapping_score > highest_mapping_score) {
                highest_mapping_score = mapping_info->mapping_score;
            }
        }

        // Add to highest scoring list if higher than the current lowest score in that list
        if (highest_mapping_score > lowest_high_score) {
            highest_mapping_scores[lowest_high_score_index] = highest_mapping_score;
            lowest_high_score_index = 0;
            lowest_high_score = highest_mapping_scores[lowest_high_score_index];

            for (unsigned int i = 0; i < num_considered_mappings; i++) {
                if (highest_mapping_scores[i] < lowest_high_score) {
                    lowest_high_score_index = i;
                    lowest_high_score = highest_mapping_scores[i];
                }
            }
        }
    }

    for (unsigned int i = 0; i < num_considered_mappings; i++) {
        score += highest_mapping_scores[i];
    }

    score /= num_considered_mappings;
    subgraph->score = score;

    free (highest_mapping_scores);
}

static inline unsigned int map (db_mapping_t    *db_mapping,
                                mapped_vertex_t *mapped_vertex,
                                vertex_t         subgraph_vertex_id_to_map)
{
    assert (subgraph_vertex_id_to_map <= db_mapping->last_mapped_vertex_id);

    for (vertex_t distance = db_mapping->last_mapped_vertex_id - subgraph_vertex_id_to_map;
         distance > 0;
         distance--) {
        mapped_vertex = mapped_vertex->parent;
    }

    return mapped_vertex->database_vertex_id;
}

static inline vertex_t reverse_map (db_mapping_t    *db_mapping,
                                    mapped_vertex_t *mapped_vertex,
                                    vertex_t         database_vertex_id_to_map)
{
    for (vertex_t vertex_id = db_mapping->last_mapped_vertex_id + 1; vertex_id-- > 0;) {
        if (mapped_vertex->database_vertex_id == database_vertex_id_to_map) {
            return vertex_id;
        }
        mapped_vertex = mapped_vertex->parent;
    }

    return VERTEX_INVALID;
}

static inline int _init_subgraph_isomorphisms_helper (config_t          *config,
                                                      database_t        *database,
                                                      embedding_store_t *embedding_store,
                                                      database_entry_t  *db_entry,
                                                      subgraph_t        *subgraph,
                                                      db_mapping_t      *db_mapping,
                                                      int                exact_mode,
                                                      int                use_refcounts)
{
    vertex_label_t *vertex_v_label = &subgraph->dfs_code[0].to_label;
    edge_label_t *edge_label = &subgraph->dfs_code[0].edge_label;
    direction_t direction = subgraph->dfs_code[0].direction;

    mapping_info_t empty_mapping_info = (mapping_info_t) {};

    unsigned int num_db_vertices = db_entry->num_vertices;
    vertex_label_t *db_labels = db_entry->vertex_labels;

    unsigned int num_new_embedding_entries = 0;
    int exact_support = 0;

    // This is the end of the list, new mappings are prepended to it
    db_mapping->mappings = NULL;

    // Initialize mappings for the first vertex
    for (vertex_t isomorphic_vertex_u = 0;
         isomorphic_vertex_u < num_db_vertices;
         isomorphic_vertex_u++)
    {
        int u_exact_match = vertex_label_compare (&db_labels[isomorphic_vertex_u], &subgraph->vertex_labels[0]) == 0;

        if (exact_mode && !u_exact_match) {
            continue;
        }

        mapping_info_t u_mapping_info = (mapping_info_t) {};

        if (!exact_mode &&
            !vertex_label_match (config,
                                 &empty_mapping_info,
                                 &db_labels[isomorphic_vertex_u],
                                 &subgraph->vertex_labels[0],
                                 &u_mapping_info))
        {
            continue;
        }

        int u_exact = exact_mode || u_exact_match;

        mapped_vertex_t *parent = NULL;

        neighbor_t *neighbors = db_entry->neighbors[isomorphic_vertex_u];
        unsigned int neighbors_len = db_entry->neighbors_len[isomorphic_vertex_u];

        // Mappings for the second vertex. Always a forward edge.
        for (unsigned int k = 0; k < neighbors_len; k++) {
            neighbor_t *neighbor = &neighbors[k];

            int v_exact_match = vertex_label_compare (&db_labels[neighbor->vertex_id], vertex_v_label) == 0;
            int e_exact_match = edge_label_compare (&neighbor->edge_label, edge_label) == 0;

            if (exact_mode && (!v_exact_match || !e_exact_match)) {
                continue;
            }

            if (direction != neighbor->direction) {
                continue;
            }

            mapping_info_t v_mapping_info = (mapping_info_t) {};

            if (!exact_mode &&
                !vertex_label_match (config,
                                     &u_mapping_info,
                                     &db_labels[neighbor->vertex_id],
                                     vertex_v_label,
                                     &v_mapping_info))
            {
                continue;
            }

            if (!exact_mode &&
                !edge_label_match (config,
                                   &u_mapping_info,
                                   &neighbor->edge_label,
                                   edge_label,
                                   &v_mapping_info))
            {
                continue;
            }

            int v_exact = exact_mode || (u_exact && v_exact_match && e_exact_match);
            exact_support |= v_exact;

            if (!parent) {
                parent = add_mapped_vertex (embedding_store, NULL, isomorphic_vertex_u, use_refcounts);
            }
            add_mapping (embedding_store, db_mapping, parent, neighbor->vertex_id, v_exact, &v_mapping_info, use_refcounts);
            num_new_embedding_entries++;
        }
    }

    db_mapping->exact = exact_support;
    db_mapping->last_mapped_vertex_id = 1;

    return num_new_embedding_entries;
}

static inline void _init_subgraph_isomorphisms (config_t          *config,
                                                database_t        *database,
                                                embedding_store_t *embedding_store,
                                                subgraph_t        *subgraph,
                                                int                exact_mode,
                                                int                use_refcounts)
{
    db_mapping_t *db_mappings = malloc (subgraph->support * sizeof(db_mapping_t));
    unsigned int new_mappings_counter = 0;

    for (unsigned int db_idx = 0; db_idx < database->num_entries; db_idx++)
    {
        if (!subgraph->supported_by[db_idx]) {
            continue;
        }

        db_mapping_t *db_mapping = &db_mappings[new_mappings_counter];
        database_entry_t *db_entry = &database->entries[db_idx];

        db_mapping->db_idx = db_idx;

        int num_new_mappings = _init_subgraph_isomorphisms_helper (config,
                                                                   database,
                                                                   embedding_store,
                                                                   db_entry,
                                                                   subgraph,
                                                                   db_mapping,
                                                                   exact_mode,
                                                                   use_refcounts);

        assert (num_new_mappings > 0);
        subgraph->num_mappings += num_new_mappings;

        new_mappings_counter++;
    }

    subgraph->db_mappings = db_mappings;
}

__attribute__((noinline)) void init_canonical_subgraph_isomorphisms (database_t        *database,
                                                                     embedding_store_t *alloc,
                                                                     subgraph_t        *subgraph)
{
    _init_subgraph_isomorphisms (NULL, database, alloc, subgraph, 1, 0);
}

__attribute__((noinline)) void init_subgraph_isomorphisms (config_t          *config,
                                                           database_t        *database,
                                                           embedding_store_t *alloc,
                                                           subgraph_t        *subgraph)
{
    _init_subgraph_isomorphisms (config, database, alloc, subgraph, 0, 1);
}

static inline void init_sig_subgraph_isomorphisms (config_t          *config,
                                                   database_t        *sig_database,
                                                   embedding_store_t *embedding_store,
                                                   subgraph_t        *subgraph)
{
    db_mapping_t *sig_db_mappings = malloc (sig_database->num_entries * sizeof(db_mapping_t));
    unsigned int new_mappings_counter = 0;
    unsigned int new_good_mappings_counter = 0;

    for (unsigned int sig_db_idx = 0;
         sig_db_idx < sig_database->num_entries;
         sig_db_idx++)
    {
        db_mapping_t *sig_db_mapping = &sig_db_mappings[new_mappings_counter];
        database_entry_t *sig_db_entry = &sig_database->entries[sig_db_idx];

        sig_db_mapping->db_idx = sig_db_idx;

        int num_new_mappings = _init_subgraph_isomorphisms_helper (config,
                                                                   sig_database,
                                                                   embedding_store,
                                                                   sig_db_entry,
                                                                   subgraph,
                                                                   sig_db_mapping,
                                                                   0,
                                                                   1);

        if (num_new_mappings) {
            new_mappings_counter++;
            if (sig_db_entry->classification) {
                new_good_mappings_counter++;
            }
        }
    }

    subgraph->sig_db_mappings = sig_db_mappings;
    subgraph->sig_support = new_mappings_counter;
    subgraph->sig_good_support = new_good_mappings_counter;
}

static inline int _update_subgraph_isomorphisms_helper (config_t          *config,
                                                        database_t        *database,
                                                        embedding_store_t *embedding_store,
                                                        database_entry_t  *db_entry,
                                                        db_mapping_t      *db_mapping,
                                                        db_mapping_t      *parent_db_mapping,
                                                        dfs_code_t        *code_to_map,
                                                        int                exact_mode,
                                                        int                use_refcounts)
{
    vertex_t vertex_u = code_to_map->from;
    vertex_t vertex_v = code_to_map->to;
    vertex_label_t *vertex_v_label = &code_to_map->to_label;
    edge_label_t *edge_label = &code_to_map->edge_label;
    direction_t direction = code_to_map->direction;

    vertex_label_t *db_labels = db_entry->vertex_labels;

    unsigned int num_new_embedding_entries = 0;
    int exact_support = 0;

    // This is the end of the list, new mappings are prepended to it
    db_mapping->mappings = NULL;

    // Forward edge
    if (vertex_v > vertex_u) {
        // A new vertex will be mapped after this
        db_mapping->last_mapped_vertex_id = parent_db_mapping->last_mapped_vertex_id + 1;

        for (mapping_t *mapping = parent_db_mapping->mappings;
             mapping;
             mapping = mapping->next)
        {
            mapped_vertex_t *mapped_vertex = mapping->mapped_vertex;

            vertex_t isomorphic_vertex_u = map (parent_db_mapping, mapped_vertex, vertex_u);

            neighbor_t *neighbors = db_entry->neighbors[isomorphic_vertex_u];
            unsigned int neighbors_len = db_entry->neighbors_len[isomorphic_vertex_u];

            for (unsigned int k = 0; k < neighbors_len; k++) {
                neighbor_t *neighbor = &neighbors[k];

                int v_exact_match = vertex_label_compare (&db_labels[neighbor->vertex_id], vertex_v_label) == 0;
                int e_exact_match = edge_label_compare (&neighbor->edge_label, edge_label) == 0;

                if (exact_mode && (!v_exact_match || !e_exact_match)) {
                    continue;
                }

                if (direction != neighbor->direction) {
                    continue;
                }

                mapping_info_t new_mapping_info = (mapping_info_t) {};

                if (!exact_mode &&
                    !vertex_label_match (config,
                                         &mapping->info,
                                         &db_labels[neighbor->vertex_id],
                                         vertex_v_label,
                                         &new_mapping_info))
                {
                    continue;
                }

                if (!exact_mode &&
                    !edge_label_match (config,
                                       &mapping->info,
                                       &neighbor->edge_label,
                                       edge_label,
                                       &new_mapping_info))
                {
                    continue;
                }

                // Only consider this vertex if it hasn't been mapped yet
                if (reverse_map(parent_db_mapping, mapped_vertex, neighbor->vertex_id) != VERTEX_INVALID)
                {
                    continue;
                }

                int exact = exact_mode || (mapping->exact && v_exact_match && e_exact_match);
                exact_support |= exact;
                add_mapping (embedding_store, db_mapping, mapped_vertex, neighbor->vertex_id, exact, &new_mapping_info, use_refcounts);
                num_new_embedding_entries++;
            }
        }
    }
    // Backward edge
    else {
        // No new vertices for backward edges
        db_mapping->last_mapped_vertex_id = parent_db_mapping->last_mapped_vertex_id;

        for (mapping_t *mapping = parent_db_mapping->mappings;
             mapping;
             mapping = mapping->next)
        {
            mapped_vertex_t *mapped_vertex = mapping->mapped_vertex;
            mapping_info_t *mapping_info = &mapping->info;

            vertex_t isomorphic_vertex_v = map (parent_db_mapping, mapped_vertex, vertex_v);
            vertex_t isomorphic_vertex_u = map (parent_db_mapping, mapped_vertex, vertex_u);

            neighbor_t *neighbors = db_entry->neighbors[isomorphic_vertex_u];
            unsigned int neighbors_len = db_entry->neighbors_len[isomorphic_vertex_u];

            for (unsigned int k = 0; k < neighbors_len; k++) {
                neighbor_t *neighbor = &neighbors[k];

                // The vertex connected by the backward edge has to be isomorphic
                // to the one in the subgraph in the current mapping
                if (isomorphic_vertex_v != neighbor->vertex_id)
                {
                    continue;
                }

                if (direction != neighbor->direction) {
                    continue;
                }

                int e_exact_match = edge_label_compare (&neighbor->edge_label, edge_label) == 0;

                if (exact_mode && !e_exact_match) {
                    continue;
                }

                mapping_info_t new_mapping_info = *mapping_info;

                if (!exact_mode &&
                    !edge_label_match (config,
                                       mapping_info,
                                       &neighbor->edge_label,
                                       edge_label,
                                       &new_mapping_info))
                {
                    continue;
                }

                int exact = exact_mode || (mapping->exact && e_exact_match);
                exact_support |= exact;
                // Copy previous generations map for this generation since there are no
                // new vertices and also update exact state as the new backwards edge
                // might have caused an inexact match
                copy_map_with_new_mapping_info (embedding_store, mapping, db_mapping, exact, &new_mapping_info, use_refcounts);
                num_new_embedding_entries++;

                // No need to check the other neighboring vertices, we already found the isomorphic one
                break;
            }
        }
    }

    db_mapping->exact = exact_support;

    return num_new_embedding_entries;
}

static inline void _update_subgraph_isomorphisms (config_t          *config,
                                                  database_t        *database,
                                                  embedding_store_t *embedding_store,
                                                  subgraph_t        *subgraph,
                                                  db_mapping_t      *parent_db_mappings_list,
                                                  unsigned int       num_parent_db_mappings,
                                                  int                exact_mode,
                                                  int                use_refcounts)
{
    unsigned int code_pos = subgraph->dfs_code_len - 1;
    dfs_code_t *code_to_map = &subgraph->dfs_code[code_pos];
    db_mapping_t *db_mappings = malloc (subgraph->support * sizeof(db_mapping_t));
    unsigned int new_mappings_counter = 0;

    for (unsigned int i = 0; i < num_parent_db_mappings; i++)
    {
        db_mapping_t *db_mapping = &db_mappings[new_mappings_counter];
        db_mapping_t *parent_db_mapping = &parent_db_mappings_list[i];
        unsigned int db_idx = parent_db_mapping->db_idx;
        database_entry_t *db_entry = &database->entries[db_idx];

        if (!subgraph->supported_by[db_idx]) {
            continue;
        }

        db_mapping->db_idx = db_idx;

        int num_new_mappings = _update_subgraph_isomorphisms_helper(config,
                                                                    database,
                                                                    embedding_store,
                                                                    db_entry,
                                                                    db_mapping,
                                                                    parent_db_mapping,
                                                                    code_to_map,
                                                                    exact_mode,
                                                                    use_refcounts);

        assert (num_new_mappings > 0);
        subgraph->num_mappings += num_new_mappings;

        new_mappings_counter++;
    }

    subgraph->db_mappings = db_mappings;
}

__attribute__((noinline)) void update_canonical_subgraph_isomorphisms (database_t        *database,
                                                                       embedding_store_t *embedding_store,
                                                                       subgraph_t        *subgraph,
                                                                       db_mapping_t      *all_parent_sig_db_mappings,
                                                                       unsigned int       num_parent_db_mappings)
{
    _update_subgraph_isomorphisms (NULL, database, embedding_store, subgraph, all_parent_sig_db_mappings, num_parent_db_mappings, 1, 0);
}

__attribute__((noinline)) void update_subgraph_isomorphisms (config_t          *config,
                                                             database_t        *database,
                                                             embedding_store_t *embedding_store,
                                                             subgraph_t        *subgraph,
                                                             db_mapping_t      *all_parent_sig_db_mappings,
                                                             unsigned int       num_parent_db_mappings)
{
    _update_subgraph_isomorphisms (config, database, embedding_store, subgraph, all_parent_sig_db_mappings, num_parent_db_mappings, 0, 1);
}

static inline void update_sig_subgraph_isomorphisms (config_t          *config,
                                                     database_t        *sig_database,
                                                     embedding_store_t *embedding_store,
                                                     subgraph_t        *subgraph,
                                                     db_mapping_t      *all_parent_sig_db_mappings,
                                                     unsigned int       num_parent_sig_db_mappings)
{
    db_mapping_t *sig_db_mappings = malloc (num_parent_sig_db_mappings * sizeof(db_mapping_t));
    unsigned int new_mappings_counter = 0;
    unsigned int new_good_mappings_counter = 0;

    unsigned int code_pos = subgraph->dfs_code_len - 1;
    dfs_code_t *code_to_map = &subgraph->dfs_code[code_pos];

    for (unsigned int i = 0; i < num_parent_sig_db_mappings; i++)
    {
        db_mapping_t *sig_db_mapping = &sig_db_mappings[new_mappings_counter];
        db_mapping_t *parent_sig_db_mapping = &all_parent_sig_db_mappings[i];
        unsigned int sig_db_idx = parent_sig_db_mapping->db_idx;
        database_entry_t *sig_db_entry = &sig_database->entries[sig_db_idx];

        sig_db_mapping->db_idx = sig_db_idx;

        int num_new_mappings = _update_subgraph_isomorphisms_helper(config,
                                                                    sig_database,
                                                                    embedding_store,
                                                                    sig_db_entry,
                                                                    sig_db_mapping,
                                                                    parent_sig_db_mapping,
                                                                    code_to_map,
                                                                    0,
                                                                    1);

        if (num_new_mappings > 0) {
            new_mappings_counter++;
            if (sig_db_entry->classification) {
                new_good_mappings_counter++;
            }
        }
    }

    subgraph->sig_db_mappings = sig_db_mappings;
    subgraph->sig_support = new_mappings_counter;
    subgraph->sig_good_support = new_good_mappings_counter;
}

static inline void update_extensions (database_t         *database,
                                      unsigned int        current_db_entry,
                                      extension_store_t  *extension_store,
                                      dfs_code_t         *extension_code)
{
    unsigned int index = dfs_code_hash (extension_code, 1) % extension_store->hash_map_size;
    extension_t *hash_map_entry = extension_store->hash_map[index];

    // only add extensions we haven't found yet but also update the support for those we have already
    while (hash_map_entry) {
        if (code_matches (&hash_map_entry->code, extension_code)) {
            if (!hash_map_entry->supported_by[current_db_entry]) {
                hash_map_entry->supported_by[current_db_entry] = 1;
                hash_map_entry->support++;
            }
            return;
        }
        hash_map_entry = hash_map_entry->next;
    }

    extension_t *extension = extension_store->extensions + extension_store->num_extensions;
    extension_store->num_extensions++;

    extension->code = *extension_code;

    extension->supported_by = calloc (database->num_entries, sizeof(extension->supported_by[0]));
    extension->supported_by[current_db_entry] = 1;

    extension->support = 1;

    extension->next = extension_store->hash_map[index];
    extension_store->hash_map[index] = extension;
}

void count_fuzzy_support (config_t          *config,
                          database_t        *database,
                          subgraph_t        *subgraph,
                          extension_store_t *extension_store)
{
    for (unsigned int db_iterator = 0; db_iterator < subgraph->support; db_iterator++) {
        db_mapping_t *db_mapping = subgraph->db_mappings + db_iterator;
        unsigned int db_idx = db_mapping->db_idx;
        vertex_label_t *db_labels = database->entries[db_idx].vertex_labels;

        unsigned int remaining_db_entries = subgraph->support - db_iterator;

        for (unsigned int i = 0; i < extension_store->num_extensions; i++) {
            extension_t *extension = &extension_store->extensions[i];

            // Has already been counted via collisions in update_extensions
            if (extension->supported_by[db_idx]) {
                continue;
            }

            // Extension can't possibly be supported if this is the case, so skip any further checks
            if (extension->support + remaining_db_entries < config->approx_min_support) {
                continue;
            }

            vertex_t vertex_u = extension->code.from;
            vertex_t vertex_v = extension->code.to;
            vertex_label_t *vertex_v_label = &extension->code.to_label;
            edge_label_t *edge_label = &extension->code.edge_label;
            direction_t direction = extension->code.direction;

            int found = 0;

            // Forward extension
            if (vertex_v > vertex_u) {
                for (mapping_t *mapping = db_mapping->mappings;
                     mapping;
                     mapping = mapping->next)
                {
                    mapped_vertex_t *mapped_vertex = mapping->mapped_vertex;
                    mapping_info_t *mapping_info = &mapping->info;

                    vertex_t isomorphic_vertex_u = map (db_mapping, mapped_vertex, vertex_u);

                    neighbor_t *neighbors = database->entries[db_idx].neighbors[isomorphic_vertex_u];
                    unsigned int neighbors_len = database->entries[db_idx].neighbors_len[isomorphic_vertex_u];

                    for (unsigned int k = 0; k < neighbors_len; k++) {
                        neighbor_t *neighbor = &neighbors[k];

                        if (direction != neighbor->direction) {
                            continue;
                        }

                        if (!vertex_label_match (config,
                                                 mapping_info,
                                                 &db_labels[neighbor->vertex_id],
                                                 vertex_v_label,
                                                 NULL))
                        {
                            continue;
                        }

                        if (!edge_label_match (config,
                                               mapping_info,
                                               &neighbor->edge_label,
                                               edge_label,NULL))
                        {
                            continue;
                        }

                        if (reverse_map(db_mapping,
                                        mapped_vertex,
                                        neighbor->vertex_id) != VERTEX_INVALID)
                        {
                            continue;
                        }

                        found = 1;
                        break;
                    }

                    if (found) {
                        break;
                    }
                }
            }
            // Backward extension
            else {
                for (mapping_t *mapping = db_mapping->mappings;
                     mapping;
                     mapping = mapping->next)
                {
                    mapped_vertex_t *mapped_vertex = mapping->mapped_vertex;
                    mapping_info_t *mapping_info = &mapping->info;

                    vertex_t isomorphic_vertex_v = map (db_mapping, mapped_vertex, vertex_v);
                    vertex_t isomorphic_vertex_u = map (db_mapping, mapped_vertex, vertex_u);

                    neighbor_t *neighbors = database->entries[db_idx].neighbors[isomorphic_vertex_u];
                    unsigned int neighbors_len = database->entries[db_idx].neighbors_len[isomorphic_vertex_u];

                    for (unsigned int k = 0; k < neighbors_len; k++) {
                        neighbor_t *neighbor = &neighbors[k];

                        if (isomorphic_vertex_v != neighbor->vertex_id) {
                            continue;
                        }

                        if (direction != neighbor->direction) {
                            continue;
                        }

                        if (!edge_label_match (config,
                                               mapping_info,
                                               &neighbor->edge_label,
                                               edge_label, NULL))
                        {
                            continue;
                        }

                        found = 1;
                        break;
                    }

                    if (found) {
                        break;
                    }
                }
            }

            if (found) {
                extension->supported_by[db_idx] = 1;
                extension->support++;
            }
        }
    }
}

// Sets the labels of code a to code b if b is smaller
static inline void minimize_labels (dfs_code_t *a,
                                    dfs_code_t *b)
{
        vertex_label_t *smaller_vertex_label, *bigger_vertex_label;
        direction_t direction;

        if (vertex_label_compare (&b->from_label, &b->to_label) < 0) {
            smaller_vertex_label = &b->from_label;
            bigger_vertex_label = &b->to_label;
            direction = b->direction;
        }
        else {
            smaller_vertex_label = &b->to_label;
            bigger_vertex_label = &b->from_label;
            // We switched from and to, so we need to change direction as well
            direction = invert_direction(b->direction);
        }

        int first_label_comparison = vertex_label_compare (smaller_vertex_label, &a->from_label);

        if (first_label_comparison < 0) {
            a->from_label = *smaller_vertex_label;
            a->to_label = *bigger_vertex_label;
            a->edge_label = b->edge_label;
            a->direction = direction;

            return;
        }

        int second_label_comparison = vertex_label_compare (bigger_vertex_label, &a->to_label);

        if (first_label_comparison == 0
            && second_label_comparison < 0) {
            a->to_label = *bigger_vertex_label;
            a->edge_label = b->edge_label;
            a->direction = direction;

            return;
        }

        int third_label_comparison = edge_label_compare (&b->edge_label, &a->edge_label);

        if (first_label_comparison == 0
            && second_label_comparison == 0
            && third_label_comparison < 0) {
            a->edge_label = b->edge_label;
            a->direction = direction;
        }

        if (first_label_comparison == 0
            && second_label_comparison == 0
            && third_label_comparison == 0
            && direction < a->direction) {
            a->direction = direction;
        }
}

static inline int canonicality_pre_check (subgraph_t *subgraph,
                                          dfs_code_t *extension)
{
    // Early partial canonicality check to avoid expensive support calculations
    dfs_code_t potential_new_start = DFS_CODE_MAX;
    potential_new_start.from = 0;
    potential_new_start.to = 1;
    minimize_labels (&potential_new_start, extension);

    if (compare_code (&potential_new_start, &subgraph->dfs_code[0]) < 0) {
        return 0;
    }

    return 1;
}

void rightmost_path_extensions (config_t          *config,
                                database_t        *database,
                                extension_store_t *extension_store,
                                subgraph_t        *subgraph)
{
    prepare_extension_store (extension_store, subgraph->num_mappings, subgraph->rightmost_path_len);

    for (unsigned int db_iterator = 0; db_iterator < subgraph->support; db_iterator++) {
        db_mapping_t *db_mapping = subgraph->db_mappings + db_iterator;

        // Only extend based on exact isomorphisms, otherwise the number of subgraphs grows too fast
        if (!db_mapping->exact) {
            continue;
        }

        unsigned int db_idx = db_mapping->db_idx;
        database_entry_t *db_entry = &database->entries[db_idx];

        // This part of the dfs code doesn't change in the loop
        dfs_code_t code;
        code.from = subgraph->rightmost_vertex;
        code.from_label = subgraph->vertex_labels[subgraph->rightmost_vertex];

        // backward extensions from rightmost child
        for (mapping_t *mapping = db_mapping->mappings;
             mapping;
             mapping = mapping->next)
        {
            // We are only interested in exact extensions
            if (!mapping->exact) {
                continue;
            }

            mapped_vertex_t *mapped_vertex = mapping->mapped_vertex;

            vertex_t isomorphic_vertex_u_r = map (db_mapping, mapped_vertex, subgraph->rightmost_vertex);
            neighbor_t *neighbors = db_entry->neighbors[isomorphic_vertex_u_r];
            unsigned int neighbors_len = db_entry->neighbors_len[isomorphic_vertex_u_r];

            for (unsigned int k = 0; k < neighbors_len; k++) {
                neighbor_t *neighbor = &neighbors[k];
                vertex_t vertex_v = reverse_map (db_mapping, mapped_vertex, neighbor->vertex_id);

                // neighbor vertex has to exist in the subgraph already if we want to add a backwards edge
                if (vertex_v == VERTEX_INVALID) {
                    continue;
                }

                // only consider vertices that are on the rightmost path
                if (rightmost_path_search(subgraph->rightmost_path, subgraph->rightmost_path_len, vertex_v) == -1) {
                    continue;
                }

                // don't add extensions for edges already in the subgraph
                if (contains_edge (subgraph, subgraph->rightmost_vertex, vertex_v, &neighbor->edge_label, neighbor->direction)) {
                    continue;
                }

                code.to = vertex_v;
                code.to_label = subgraph->vertex_labels[vertex_v];
                code.edge_label = neighbor->edge_label;
                code.direction = neighbor->direction;

                if (!canonicality_pre_check (subgraph, &code)) {
                    continue;
                }

                update_extensions (database, db_idx, extension_store, &code);
            }
        }

        // forward extensions from vertices on rightmost path
        for (mapping_t *mapping = db_mapping->mappings;
             mapping;
             mapping = mapping->next)
        {
            // We are only interested in exact extensions
            if (!mapping->exact) {
                continue;
            }

            mapped_vertex_t *mapped_vertex = mapping->mapped_vertex;

            for (unsigned int k = subgraph->rightmost_path_len; k-- > 0;) {
                vertex_t vertex_u = subgraph->rightmost_path[k];
                vertex_t isomorphic_vertex_u = map (db_mapping, mapped_vertex, vertex_u);
                neighbor_t *neighbors = db_entry->neighbors[isomorphic_vertex_u];
                unsigned int neighbors_len = db_entry->neighbors_len[isomorphic_vertex_u];

                // This part of the dfs code doesn't change in the loop
                code.from = vertex_u;
                code.from_label = subgraph->vertex_labels[vertex_u];

                for (unsigned int l = 0; l < neighbors_len; l++) {
                    neighbor_t *neighbor = &neighbors[l];
                    vertex_t vertex_x = neighbor->vertex_id;

                    // New vertex from forward extension can't already exist in subgraph
                    if (reverse_map (db_mapping, mapped_vertex, vertex_x) != VERTEX_INVALID) {
                        continue;
                    }

                    code.to = subgraph->rightmost_vertex + 1;
                    code.to_label = db_entry->vertex_labels[vertex_x];
                    code.edge_label = neighbor->edge_label;
                    code.direction = neighbor->direction;

                    if (!canonicality_pre_check (subgraph, &code)) {
                        continue;
                    }

                    update_extensions (database, db_idx, extension_store, &code);
                }
            }
        }
    }
    count_fuzzy_support (config, database, subgraph, extension_store);
}

subgraph_t new_subgraph_from_extension (subgraph_t  *subgraph,
                                        extension_t *extension)
{
    assert (extension->code.from <= subgraph->rightmost_vertex);
    assert (extension->code.to <= subgraph->rightmost_vertex + 1);

    subgraph_t result = {};

    result.dfs_code_len = subgraph->dfs_code_len + 1;
    result.dfs_code = malloc (result.dfs_code_len * sizeof(result.dfs_code[0]));
    memcpy (result.dfs_code,
            subgraph->dfs_code,
            subgraph->dfs_code_len * sizeof(result.dfs_code[0]));
    result.dfs_code[result.dfs_code_len - 1] = extension->code;

    result.support = extension->support;
    result.supported_by = extension->supported_by;

    // Forward edge
    if (extension->code.to > extension->code.from) {
        assert (extension->code.to == subgraph->rightmost_vertex + 1);

        // reuse the old right most path up to the last vertex before we added our new vertex
        // vertex IDs on the rightmost path SHOULD be always increasing, so we can use binary search here
        int pos = rightmost_path_search (subgraph->rightmost_path, subgraph->rightmost_path_len, extension->code.from);
        result.rightmost_path = malloc ((pos + 2) * sizeof(result.rightmost_path[0]));
        memcpy (result.rightmost_path,
                subgraph->rightmost_path,
                (pos + 1) * sizeof(subgraph->rightmost_path[0]));
        result.rightmost_path[pos + 1] = extension->code.to;

        result.rightmost_path_len = pos + 2;

        result.rightmost_vertex = extension->code.to;

        unsigned int num_vertices = result.rightmost_vertex + 1;
        result.vertex_labels = malloc (num_vertices * sizeof(result.vertex_labels[0]));
        memcpy (result.vertex_labels,
                subgraph->vertex_labels,
                (num_vertices - 1) * sizeof(result.vertex_labels[0]));
        result.vertex_labels[num_vertices - 1] = extension->code.to_label;

        // Forward extensions can only turn paths into trees, everything else is
        // guaranteed to stay the same
        if (subgraph->graph_type == PATH) {
            result.graph_type = PATH;
            // 0 is the only vertex which can occur twice as "from" in a path
            if (extension->code.from == 0) {
                // the first edge always start at 0
                for (unsigned int i = 1; i < subgraph->dfs_code_len; i++) {
                    if (subgraph->dfs_code[i].from == 0) {
                        result.graph_type = TREE;
                        break;
                    }
                }
            }
            // Extensions of vertices that are neither the rightmost vertex nor 0 are
            // guaranteed to create branches
            else if (extension->code.from != subgraph->rightmost_vertex) {
                result.graph_type = TREE;
            }
        }
        else {
            result.graph_type = subgraph->graph_type;
        }
    }
    // Backward edge
    else {
        result.rightmost_path = malloc (subgraph->rightmost_path_len * sizeof(result.rightmost_path[0]));
        memcpy (result.rightmost_path,
                subgraph->rightmost_path,
                subgraph->rightmost_path_len * sizeof(subgraph->rightmost_path[0]));

        result.rightmost_path_len = subgraph->rightmost_path_len;

        result.rightmost_vertex = subgraph->rightmost_vertex;

        unsigned int num_vertices = subgraph->rightmost_vertex + 1;
        result.vertex_labels = malloc (num_vertices * sizeof(result.vertex_labels[0]));
        memcpy (result.vertex_labels,
                subgraph->vertex_labels,
                num_vertices * sizeof(result.vertex_labels[0]));

        result.graph_type = CYCLE;
    }

    return result;
}

database_entry_t canonical_database_entry (subgraph_t *subgraph) {
    database_entry_t db_entry;
    unsigned int num_vertices = subgraph->rightmost_vertex + 1;
    db_entry.num_vertices = num_vertices;
    db_entry.num_edges = subgraph->dfs_code_len;
    db_entry.vertex_labels = subgraph->vertex_labels;
    db_entry.neighbors = malloc (num_vertices * sizeof(db_entry.neighbors[0]));
    db_entry.neighbors_len = calloc (num_vertices, sizeof(db_entry.neighbors_len[0]));
    neighbor_t *neighbors_alloc = malloc (2 * subgraph->dfs_code_len * sizeof(neighbors_alloc[0]));

    for (dfs_code_t *cur = subgraph->dfs_code; cur < subgraph->dfs_code + subgraph->dfs_code_len; cur++) {
        // Count number of neighbors
        db_entry.neighbors_len[cur->from]++;
        db_entry.neighbors_len[cur->to]++;
    }

    unsigned int index = 0;
    for (unsigned int i = 0; i < num_vertices; i++) {
        db_entry.neighbors[i] = &neighbors_alloc[index];
        index += db_entry.neighbors_len[i];
    }

    // reuse this memory as a counter for the next step
    memset (db_entry.neighbors_len, 0, num_vertices * sizeof(db_entry.neighbors_len[0]));

    for (dfs_code_t *cur = subgraph->dfs_code; cur < subgraph->dfs_code + subgraph->dfs_code_len; cur++) {
        db_entry.neighbors[cur->from][db_entry.neighbors_len[cur->from]].vertex_id = cur->to;
        db_entry.neighbors[cur->from][db_entry.neighbors_len[cur->from]].edge_label = cur->edge_label;
        db_entry.neighbors[cur->from][db_entry.neighbors_len[cur->from]].direction = cur->direction;
        db_entry.neighbors[cur->to][db_entry.neighbors_len[cur->to]].vertex_id = cur->from;
        db_entry.neighbors[cur->to][db_entry.neighbors_len[cur->to]].edge_label = cur->edge_label;
        db_entry.neighbors[cur->to][db_entry.neighbors_len[cur->to]].direction = invert_direction(cur->direction);
        db_entry.neighbors_len[cur->from]++;
        db_entry.neighbors_len[cur->to]++;
    }

    return db_entry;
}

subgraph_t canonical_subgraph_template (subgraph_t *subgraph)
{
    subgraph_t result = {};
    unsigned int num_vertices = subgraph->rightmost_vertex + 1;

    result.dfs_code_len = subgraph->dfs_code_len;
    result.dfs_code = malloc (result.dfs_code_len * sizeof(result.dfs_code[0]));
    result.rightmost_path = malloc (num_vertices * sizeof(result.rightmost_path[0]));
    result.vertex_labels = malloc (num_vertices * sizeof(result.vertex_labels[0]));

    result.support = 1;
    result.supported_by = malloc (1 * sizeof(result.supported_by[0]));
    result.supported_by[0] = 1;

    return result;
}

void free_canonical_subgraph (subgraph_t *subgraph)
{
    free (subgraph->dfs_code);
    free (subgraph->rightmost_path);
    free (subgraph->vertex_labels);
    free (subgraph->supported_by);
    free (subgraph->db_mappings);
}

void reset_canonical_subgraph_template (subgraph_t *subgraph,
                                        dfs_code_t *initial_entry)
{
    subgraph->dfs_code_len = 1;
    subgraph->dfs_code[0] = *initial_entry;
    subgraph->rightmost_path[0] = 0;
    subgraph->rightmost_path[1] = 1;
    subgraph->rightmost_vertex = 1;
    subgraph->rightmost_path_len = 2;
    subgraph->vertex_labels[0] = initial_entry->from_label;
    subgraph->vertex_labels[1] = initial_entry->to_label;
    free (subgraph->db_mappings);
    subgraph->db_mappings = NULL;
    subgraph->num_mappings = 0;
}

void update_canonical_subgraph_template (subgraph_t  *subgraph,
                                         extension_t *extension)
{
    subgraph->dfs_code[subgraph->dfs_code_len] = extension->code;
    subgraph->dfs_code_len++;

    // Forward edge
    if (extension->code.to > extension->code.from) {
        int pos = rightmost_path_search (subgraph->rightmost_path, subgraph->rightmost_path_len, extension->code.from);
        subgraph->rightmost_path[pos + 1] = extension->code.to;
        subgraph->rightmost_path_len = pos + 2;
        subgraph->rightmost_vertex = extension->code.to;
        subgraph->vertex_labels[extension->code.to] = extension->code.to_label;
    }
    // No need to do anything special for backward edges
}

extension_t find_canonical_extension (database_t *database,
                                      subgraph_t *subgraph)
{
    int extension_found = 0;
    extension_t result;

    result.code = DFS_CODE_MAX;
    result.code.from = 0;

    db_mapping_t *db_mapping = &subgraph->db_mappings[0];

    for (mapping_t *mapping = db_mapping->mappings;
         mapping;
         mapping = mapping->next)
    {
        mapped_vertex_t *mapped_vertex = mapping->mapped_vertex;

        // backward extensions from rightmost child
        vertex_t isomorphic_vertex_u_r = map (db_mapping, mapped_vertex, subgraph->rightmost_vertex);
        neighbor_t *neighbors = database->entries[0].neighbors[isomorphic_vertex_u_r];
        unsigned int neighbors_len = database->entries[0].neighbors_len[isomorphic_vertex_u_r];

        // This part of the dfs code doesn't change in the loop
        dfs_code_t code;
        code.from = subgraph->rightmost_vertex;
        code.from_label = subgraph->vertex_labels[subgraph->rightmost_vertex];

        for (unsigned int k = 0; k < neighbors_len; k++) {
            neighbor_t *neighbor = &neighbors[k];
            vertex_t vertex_v = reverse_map (db_mapping, mapped_vertex, neighbor->vertex_id);

            // neighbor vertex has to exist in the subgraph already if we want to add a backwards edge
            if (vertex_v == VERTEX_INVALID) {
                continue;
            }

            // only consider vertices that are on the rightmost path
            if (rightmost_path_search(subgraph->rightmost_path, subgraph->rightmost_path_len, vertex_v) == -1) {
                continue;
            }

            // don't add extensions for edges already in the subgraph
            if (contains_edge (subgraph, subgraph->rightmost_vertex, vertex_v, &neighbor->edge_label, neighbor->direction)) {
                continue;
            }

            code.to = vertex_v;
            code.to_label = subgraph->vertex_labels[vertex_v];
            code.edge_label = neighbor->edge_label;
            code.direction = neighbor->direction;

            if (compare_code (&code, &result.code) == -1) {
                result.code = code;
                extension_found = 1;
            }
        }

        // No forward edge can have a smaller DFS code than a backward edge in
        // between two smaller vertices, regardless of the mapping
        if (extension_found) {
            continue;
        }

        // forward extensions from vertices on rightmost path
        for (unsigned int k = subgraph->rightmost_path_len; k-- > 0;) {
            vertex_t vertex_u = subgraph->rightmost_path[k];
            vertex_t isomorphic_vertex_u = map (db_mapping, mapped_vertex, vertex_u);
            neighbors = database->entries[0].neighbors[isomorphic_vertex_u];
            neighbors_len = database->entries[0].neighbors_len[isomorphic_vertex_u];

            // This part of the dfs code doesn't change in the loop
            code.from = vertex_u;
            code.from_label = database->entries[0].vertex_labels[isomorphic_vertex_u];

            for (unsigned int l = 0; l < neighbors_len; l++) {
                neighbor_t *neighbor = &neighbors[l];
                vertex_t vertex_x = neighbor->vertex_id;
                if (reverse_map (db_mapping, mapped_vertex, vertex_x) == VERTEX_INVALID) {
                    code.to = subgraph->rightmost_vertex + 1;
                    code.to_label = database->entries[0].vertex_labels[vertex_x];
                    code.edge_label = neighbor->edge_label;
                    code.direction = neighbor->direction;

                    if (compare_code (&code, &result.code) == -1) {
                        result.code = code;
                    }
                }
            }
        }
    }

    assert (result.code.to != DFS_CODE_MAX.to);
    assert (vertex_label_compare (&result.code.to_label, &DFS_CODE_MAX.to_label) < 0);
    assert (vertex_label_compare (&result.code.from_label, &DFS_CODE_MAX.from_label) < 0);
    assert (edge_label_compare (&result.code.edge_label, &DFS_CODE_MAX.edge_label) < 0);

    return result;
}

// This function assumes that the first n-1 DFS code entries are already canonical
int is_canonical (subgraph_t        *subgraph,
                  database_t        *canonical_db,
                  embedding_store_t *canon_test_embedding_store)
{
    dfs_code_t initial_code = subgraph->dfs_code[0];

    // Find smallest label, this can be either the current first DFS code or the new extension
    minimize_labels (&initial_code, &subgraph->dfs_code[subgraph->dfs_code_len - 1]);

    subgraph_t canonical_subgraph = canonical_subgraph_template (subgraph);
    reset_canonical_subgraph_template (&canonical_subgraph, &initial_code);

    // Generate mapping for the first 2 vertices of the initial DFS code
    init_canonical_subgraph_isomorphisms (canonical_db, canon_test_embedding_store, &canonical_subgraph);

    for (unsigned int i = 1; i < subgraph->dfs_code_len; i++) {
        extension_t canonical_extension = find_canonical_extension (canonical_db, &canonical_subgraph);
        int canonical = (compare_code (&subgraph->dfs_code[i], &canonical_extension.code) < 1);

        if (canonical) {
            db_mapping_t *prev_mappings = canonical_subgraph.db_mappings;
            update_canonical_subgraph_template (&canonical_subgraph, &canonical_extension);
            update_canonical_subgraph_isomorphisms (canonical_db, canon_test_embedding_store, &canonical_subgraph, prev_mappings, 1);
            free (prev_mappings);
        }
        else {
            free_canonical_subgraph (&canonical_subgraph);
            clear_embedding_store (canon_test_embedding_store);
            return 0;
        }
    }
    free_canonical_subgraph (&canonical_subgraph);
    clear_embedding_store (canon_test_embedding_store);
    return 1;
}

// Note that parent can be NULL when using approximate matching
subgraph_t *find_canonical_parent (subgraph_t         *subgraph,
                                   database_t         *parent_db,
                                   subgraph_t         *canonical_subgraph,
                                   unsigned int        edge_to_remove,
                                   subgraph_t        **subgraphs_hash_map,
                                   unsigned int        hash_map_size,
                                   embedding_store_t  *canon_test_embedding_store)
{
    database_entry_t *parent_db_entry = &parent_db->entries[0];
    parent_db_entry->num_edges = subgraph->dfs_code_len - 1;
    parent_db_entry->num_vertices = subgraph->rightmost_vertex + 1;

    // Find smallest possible starting entry
    dfs_code_t initial_code = DFS_CODE_MAX;
    initial_code.from = 0;
    initial_code.to = 1;
    for (unsigned int i = 0; i < subgraph->dfs_code_len; i++) {
        // Skip the removed edge
        if (i == edge_to_remove) {
            continue;
        }

        dfs_code_t *cur = subgraph->dfs_code + i;
        minimize_labels (&initial_code, cur);
    }

    // Move neighbor entries from the removed entry to the last position and reduce the neighbor length to hide them
    dfs_code_t *removed_entry = &subgraph->dfs_code[edge_to_remove];
    unsigned int rem_from = removed_entry->from;
    unsigned int rem_to = removed_entry->to;
    neighbor_t *last_neighbor;
    neighbor_t temp;
    for (unsigned int i = 0; i < parent_db_entry->neighbors_len[rem_from] - 1; i++) {
        if (parent_db_entry->neighbors[rem_from][i].vertex_id == rem_to) {
            last_neighbor = &parent_db_entry->neighbors[rem_from][parent_db_entry->neighbors_len[rem_from] - 1];
            temp = *last_neighbor;
            *last_neighbor = parent_db_entry->neighbors[rem_from][i];
            parent_db_entry->neighbors[rem_from][i] = temp;
            break;
        }
    }
    parent_db_entry->neighbors_len[rem_from]--;
    for (unsigned int i = 0; i < parent_db_entry->neighbors_len[rem_to] - 1; i++) {
        if (parent_db_entry->neighbors[rem_to][i].vertex_id == rem_from) {
            last_neighbor = &parent_db_entry->neighbors[rem_to][parent_db_entry->neighbors_len[rem_to] - 1];
            temp = *last_neighbor;
            *last_neighbor = parent_db_entry->neighbors[rem_to][i];
            parent_db_entry->neighbors[rem_to][i] = temp;
            break;
        }
    }
    parent_db_entry->neighbors_len[rem_to]--;

    reset_canonical_subgraph_template (canonical_subgraph, &initial_code);

    // Generate mapping for the first 2 vertices of the initial DFS code
    init_canonical_subgraph_isomorphisms (parent_db, canon_test_embedding_store, canonical_subgraph);

    for (unsigned int i = 1; i < parent_db_entry->num_edges; i++) {
        extension_t canonical_extension = find_canonical_extension (parent_db, canonical_subgraph);
        db_mapping_t *prev_mappings = canonical_subgraph->db_mappings;
        update_canonical_subgraph_template (canonical_subgraph, &canonical_extension);
        update_canonical_subgraph_isomorphisms (parent_db, canon_test_embedding_store, canonical_subgraph, prev_mappings, 1);
        free (prev_mappings);
    }

    calc_subgraph_hash (canonical_subgraph);

    unsigned int slot = canonical_subgraph->hash % hash_map_size;
    subgraph_t *hash_map_entry = subgraphs_hash_map[slot];
    subgraph_t *parent = NULL;

    while (hash_map_entry) {
        if (hash_map_entry->hash == canonical_subgraph->hash) {
            parent = hash_map_entry;
            for (unsigned int i = 0; i < parent_db_entry->num_edges; i++) {
                if (!code_matches (canonical_subgraph->dfs_code + i,
                                   parent->dfs_code + i)){
                    parent = NULL;
                    break;
                }
            }

            if (parent != NULL) {
                break;
            }
        }
        hash_map_entry = hash_map_entry->next;
    }

    // Restore original neighbor length to add the removed edge back
    parent_db_entry->neighbors_len[rem_from]++;
    parent_db_entry->neighbors_len[rem_to]++;

    clear_embedding_store (canon_test_embedding_store);

    return parent;
}

// Find all edges that can be removed without resulting in a disconnected graph
void find_removable_edges (subgraph_t        *subgraph,
                           rem_edge_store_t  *rem_edge_store,
                           char             **out_removable_edges)
{
    char *removable_edges = calloc (subgraph->dfs_code_len, sizeof(removable_edges));

    unsigned int num_vertices = subgraph->rightmost_vertex + 1;
    unsigned int *cycles_closed_by = rem_edge_store->cycles_closed_by;
    unsigned int *cycles_open_at = rem_edge_store->cycles_open_at;
    unsigned int *edge_count = rem_edge_store->edge_count;
    memset (cycles_closed_by, 0, sizeof(cycles_closed_by[0]) * num_vertices);
    memset (cycles_open_at, 0, sizeof(cycles_open_at[0]) * num_vertices);
    memset (edge_count, 0, sizeof(edge_count[0]) * num_vertices);

    for (unsigned int i = subgraph->dfs_code_len; i-- > 0;) {
        dfs_code_t *current_code = subgraph->dfs_code + i;

        // Forward edges
        if (current_code->from < current_code->to) {
            // Update closed cycles
            cycles_open_at[current_code->to] -= cycles_closed_by[current_code->to];
            // Inherit opened cycles from child vertex
            cycles_open_at[current_code->from] += cycles_open_at[current_code->to];
            // Count the number of edges going out from the "from" vertex
            // There is no need to count the incoming edges of the "to" vertex here
            edge_count[current_code->from]++;
            edge_count[current_code->to]++;

            // No more edges to "to" after this DFS code entry, so check if
            // there are any open cycles at this point.
            // Also check the number of edges to determine if this is a terminal
            // vertex.
            if (cycles_open_at[current_code->to] > 0
                || edge_count[current_code->to] == 1)
            {
                removable_edges[i] = 1;
            }
        }

        // Backward edges
        else {
            // Backward edges are always part of a cycle
            removable_edges[i] = 1;

            // Remember that a new cycle started and when it will close again
            cycles_open_at[current_code->from]++;
            cycles_closed_by[current_code->to]++;

            // Count the edge for both vertices
            edge_count[current_code->from]++;
            edge_count[current_code->to]++;
        }
    }

    // The first DFS code entry needs a special case for the starting vertex
    if (edge_count[0] == 1) {
        removable_edges[0] = 1;
    }

    *out_removable_edges = removable_edges;
}

int parent_check (config_t           *config,
                  database_t         *database,
                  subgraph_t         *subgraph,
                  database_t         *canonicality_db,
                  subgraph_t        **subgraphs_hash_map,
                  unsigned int        hashmap_size,
                  embedding_store_t  *canon_test_embedding_store,
                  rem_edge_store_t   *rem_edge_store,
                  int                 directly_supported)
{
    unsigned int max_allowed_parent_support = directly_supported ? database->num_entries : subgraph->support + config->approx_max_support_difference;
    unsigned int max_allowed_parent_score = subgraph->score - config->min_score_difference;
    subgraph_t canonical_subgraph = canonical_subgraph_template (subgraph);
    int non_supported_parent_found = 0;
    char *removable_edges;
    unsigned int num_parents = 0;
    subgraph_t **parents = malloc (subgraph->dfs_code_len * sizeof(parents[0]));

    find_removable_edges (subgraph, rem_edge_store, &removable_edges);

    for (unsigned int i = 0; i < subgraph->dfs_code_len; i++) {
        if (removable_edges[i]) {
            subgraph_t *parent;

            parent = find_canonical_parent (subgraph,
                                            canonicality_db,
                                            &canonical_subgraph,
                                            i,
                                            subgraphs_hash_map,
                                            hashmap_size,
                                            canon_test_embedding_store);
            if (parent) {
                // Check if parent has already been found (in case of symmetries)
                int found = 0;
                for (unsigned int j = 0; j < num_parents; j++) {
                    if (parents[j] == parent) {
                        found = 1;
                        break;
                    }
                }
                if (!found) {
                    if (parent->support <= max_allowed_parent_support && parent->score <= max_allowed_parent_score) {
                        parents[num_parents++] = parent;
                    }
                    else {
                        non_supported_parent_found = 1;
                        break;
                    }
                }
            }
            else {
                non_supported_parent_found = 1;
                break;
            }
        }
        if (non_supported_parent_found) {
            break;
        }
    }

    if (!non_supported_parent_found) {
        parents = realloc (parents, num_parents * sizeof(parents[0]));
    }

    free_canonical_subgraph (&canonical_subgraph);

    subgraph->parents = parents;
    subgraph->num_parents = num_parents;

    subgraph->parent_ids = malloc (num_parents * sizeof(subgraph->parent_ids[0]));

    for (unsigned int i = 0; i < num_parents; i++) {
        subgraph->parent_ids[i] = parents[i]->id;
    }

    free (removable_edges);

    return !non_supported_parent_found;
}

void init_output (output_t     *output,
                  char         *path,
                  off_t         max_file_size,
                  unsigned int  min_output_size,
                  int           use_distances,
                  int           use_directions)
{
    *output = (output_t) {};

    if (!path) {
        output->fd = -1;
        return;
    }

    output->fd = open (path, O_RDWR | O_CREAT | O_TRUNC, 00644);
    if (output->fd == -1) {
        fprintf (stderr, "Could not open output file: %s\n", strerror(errno));
        exit (1);
    }

    ftruncate (output->fd, max_file_size);

    output->buffer = mmap(NULL, max_file_size, PROT_READ | PROT_WRITE, MAP_SHARED, output->fd, 0);
    if (output->buffer == MAP_FAILED) {
        fprintf (stderr, "Could not mmap output file: %s\n", strerror(errno));
        exit (1);
    }

    output->buffer_size = max_file_size;

    output->min_size = min_output_size;

    output->use_distances = use_distances;
    output->use_directions = use_directions;

    pthread_mutex_init (&output->mutex, NULL);
}

void write_header_entry (output_t *output,
                         char     *entry,
                         off_t     size)
{
    if (output->fd == -1) {
        return;
    }

    if (size > output->buffer_size) {
        fprintf (stderr, "Output file too small for header\n");
        exit (1);
    }

    memcpy (output->buffer, entry, size);
    output->buffer += size;
    output->buffer_size -= size;
    output->header_bytes_used += size;
}

static inline void swap_range (char  *buffer,
                               off_t  first_start,
                               off_t  second_start,
                               off_t  size)
{
    // No overlaps allowed
    assert (labs(first_start - second_start) >= size);

    int temp_size = 16 * 1024 * 1024; // 16MiB
    char *temp = malloc(temp_size);
    off_t to_swap = size;
    off_t swapped = 0;
    while (swapped < to_swap) {
        off_t swap_size = MIN(temp_size, to_swap - swapped);
        memcpy (temp,
                buffer + first_start + swapped,
                swap_size);
        memcpy (buffer + first_start + swapped,
                buffer + second_start + swapped,
                swap_size);
        memcpy (buffer + second_start + swapped,
                temp,
                swap_size);
        swapped += swap_size;
    }
    free (temp);
}

void finalize_output (output_t *output)
{
    if (output->fd == -1) {
        return;
    }

    off_t start = output->cursor;
    off_t size = 0;

    unsigned int start_barrier_index = output->min_size - 1;
    while (start_barrier_index < output->starts_len) {
        off_t size_required = output->bytes_written - output->start_barriers[start_barrier_index];
        if (size_required <= output->buffer_size){
            start = output->start_barriers[start_barrier_index] % output->buffer_size;
            size = size_required;
            break;
        }

        fprintf (stderr,
                 "Skipping output of size %u: output file too small, would need at least %ld Bytes\n",
                 start_barrier_index + 1,
                 size_required + output->header_bytes_used);
        start_barrier_index++;
    }

    free (output->start_barriers);

    // buffer already linear, just needs moving to the beginning of the buffer
    if (start <= output->cursor) {
        memmove (output->buffer, output->buffer + start, output->cursor - start);
    }

    // in place linearization
    else {
        off_t first = 0;
        off_t middle = start;
        off_t last = output->buffer_size;
        off_t next = middle;
        off_t rotation_size;
        while (first < size) {
            rotation_size = MIN(middle - first, last - middle);
            swap_range (output->buffer, first, next, rotation_size);
            first += rotation_size;
            next += rotation_size;
            if (next == last) {
                next = middle;
            }
            else if (first == middle) {
                middle = next;
            }

        }
    }

    output->buffer -= output->header_bytes_used;
    output->buffer_size += output->header_bytes_used;
    size += output->header_bytes_used;

    munmap (output->buffer, output->buffer_size);
    ftruncate (output->fd, size);
    close (output->fd);

}

static inline void mark_output_start (output_t *output)
{
    if (output->fd == -1) {
        return;
    }

    if (output->starts_len % 10 == 0) {
        output->start_barriers = realloc (output->start_barriers, (output->starts_len + 10) * sizeof(output->start_barriers[0]));
    }

    output->start_barriers[output->starts_len] = output->bytes_written;
    output->starts_len++;
}

static inline void write_output (output_t *output,
                                 char     *out_string,
                                 off_t     out_size)
{
    // split output string into two parts if close to the end of the ring buffer
    off_t first_part_size = MIN(out_size, output->buffer_size - output->cursor);
    off_t second_part_size = out_size - first_part_size;

    memcpy (output->buffer + output->cursor, out_string, first_part_size);
    output->cursor += first_part_size;

    if (second_part_size > output->buffer_size) {
        fprintf (stderr, "Output file too small\n");
        exit (1);
    }

    if (second_part_size) {
        output->cursor = 0;
        memcpy (output->buffer, out_string + first_part_size, second_part_size);
        output->cursor += second_part_size;
    }

    output->bytes_written += out_size;
}

static inline void write_thread_buffer (char         **out_buffer,
                                        unsigned int  *out_buffer_size,
                                        unsigned int  *out_buffer_alloc,
                                        char          *in_buffer,
                                        unsigned int   in_buffer_size)
{
    int changed_alloc = 0;
    while (*out_buffer_size + in_buffer_size > *out_buffer_alloc) {
        *out_buffer_alloc *= 2;
        changed_alloc = 1;
    }
    if (changed_alloc) {
        *out_buffer = realloc (*out_buffer, *out_buffer_alloc);
    }
    memcpy (*out_buffer + *out_buffer_size, in_buffer, in_buffer_size);
    *out_buffer_size += in_buffer_size;
}

void write_subgraph (database_t    *database,
                     database_t    *sig_database,
                     subgraph_t    *subgraph,
                     thread_data_t *thread_data)
{
    output_t *output = thread_data->output;

    if (output->fd == -1) {
        return;
    }

    if (subgraph->dfs_code_len < output->min_size) {
        return;
    }

    unsigned int out_buffer_size = 0;
    unsigned int out_buffer_alloc = 16 * 1024;
    char *out_buffer = malloc (out_buffer_alloc);
    char line_buffer[4096];
    char *line_buf_pos = line_buffer;
    unsigned int line_size = 0;

    line_size = sprintf (line_buffer, "t %u\n", subgraph->id);
    write_thread_buffer (&out_buffer, &out_buffer_size, &out_buffer_alloc, line_buffer, line_size);

    if (subgraph->num_parents > 0) {
        line_size = sprintf (line_buffer, "# p");
        line_buf_pos = line_buffer + line_size;

        for (unsigned int i = 0; i < subgraph->num_parents; i++) {
            line_size += sprintf (line_buf_pos, " %u", subgraph->parent_ids[i]);
            line_buf_pos = line_buffer + line_size;
        }

        line_size += sprintf (line_buf_pos, "\n");
        line_buf_pos = line_buffer + line_size;

        write_thread_buffer (&out_buffer, &out_buffer_size, &out_buffer_alloc, line_buffer, line_size);
    }

    line_size = sprintf (line_buffer, "# t %u\n", subgraph->graph_type);
    write_thread_buffer (&out_buffer, &out_buffer_size, &out_buffer_alloc, line_buffer, line_size);

    line_size = sprintf (line_buffer, "# s %d\n", subgraph->score);
    write_thread_buffer (&out_buffer, &out_buffer_size, &out_buffer_alloc, line_buffer, line_size);

    if (sig_database) {
        unsigned int num_classified_bad = sig_database->num_entries - sig_database->num_classified_good;
        unsigned int support_bad = subgraph->sig_support - subgraph->sig_good_support;
        line_size = sprintf (line_buffer, "# i %u %u %u %u\n",
                             subgraph->sig_good_support,
                             sig_database->num_classified_good,
                             support_bad,
                             num_classified_bad);
        write_thread_buffer (&out_buffer, &out_buffer_size, &out_buffer_alloc, line_buffer, line_size);
    }

    unsigned int *temp = malloc ((subgraph->rightmost_vertex + 1) * sizeof(temp[0]));

    for (unsigned int db_iterator = 0; db_iterator < subgraph->support; db_iterator++) {
        db_mapping_t *db_mapping = subgraph->db_mappings + db_iterator;
        unsigned int db_id = database->entries[db_mapping->db_idx].id;

        for (mapping_t *mapping = db_mapping->mappings;
             mapping;
             mapping = mapping->next)
        {
            mapped_vertex_t *mapped_vertex = mapping->mapped_vertex;
            mapping_info_t *mapping_info = &mapping->info;

            int score_contribution = mapping_info->mapping_score;
            int exact = mapping->exact;

            for (vertex_t k = db_mapping->last_mapped_vertex_id + 1; k-- > 0;) {
                temp[k] = map (db_mapping, mapped_vertex, k);
            }

            line_size = sprintf (line_buffer, "# m %u : %d : %d :", db_id, exact, score_contribution);
            line_buf_pos = line_buffer + line_size;

            for (vertex_t k = 0; k <= db_mapping->last_mapped_vertex_id; k++) {
                line_size += sprintf (line_buf_pos, " %u", temp[k]);
                line_buf_pos = line_buffer + line_size;
            }

            line_size += sprintf (line_buf_pos, "\n");
            line_buf_pos = line_buffer + line_size;

            write_thread_buffer (&out_buffer, &out_buffer_size, &out_buffer_alloc, line_buffer, line_size);
        }
    }

    free (temp);

    unsigned int num_labels = subgraph->rightmost_vertex + 1;
    for (unsigned int i = 0; i < num_labels; i++) {
        line_size = sprintf (line_buffer, "v %u %u\n", i, subgraph->vertex_labels[i]);
        write_thread_buffer (&out_buffer, &out_buffer_size, &out_buffer_alloc, line_buffer, line_size);
    }

    for (unsigned int i = 0; i < subgraph->dfs_code_len; i++) {
        line_size = sprintf (line_buffer,
                             "e %u %u %d",
                             subgraph->dfs_code[i].from,
                             subgraph->dfs_code[i].to,
                             subgraph->dfs_code[i].edge_label.interaction_type);
        line_buf_pos = line_buffer + line_size;

        if (output->use_distances) {
            line_size += sprintf (line_buf_pos,
                                  " %d",
                                  subgraph->dfs_code[i].edge_label.distance);
            line_buf_pos = line_buffer + line_size;
        }

        if (output->use_directions) {
            line_size += sprintf (line_buf_pos,
                                  " %d",
                                  subgraph->dfs_code[i].direction);
            line_buf_pos = line_buffer + line_size;
        }

        line_size += sprintf (line_buf_pos, "\n");
        line_buf_pos = line_buffer + line_size;

        write_thread_buffer (&out_buffer, &out_buffer_size, &out_buffer_alloc, line_buffer, line_size);
    }

    pthread_mutex_lock (&output->mutex);
    write_output (output, out_buffer, out_buffer_size);
    pthread_mutex_unlock (&output->mutex);
    free (out_buffer);
}

void *gspan_step_helper (void *data)
{
    thread_data_t *thread_data = data;
    config_t *config = thread_data->config;
    database_t *database = thread_data->database;
    database_t *sig_database = thread_data->sig_database;
    extension_store_t *extension_storage = &thread_data->extension_store;
    embedding_store_t *embedding_store = &thread_data->embedding_store;
    embedding_store_t *sig_embedding_store = &thread_data->sig_embedding_store;
    embedding_store_t *canon_embedding_store = &thread_data->canon_embedding_store;
    subgraph_store_t *in_subgraphs = thread_data->in_subgraphs;
    subgraph_store_t *out_subgraphs = &thread_data->out_subgraphs;
    rem_edge_store_t *rem_edge_store = &thread_data->rem_edge_store;
    unsigned int num_new_subgraphs = 0;
    subgraph_t *new_subgraphs = NULL;

    unsigned int current_item;

    while (current_item = __sync_fetch_and_add (thread_data->current_item, 1),
           current_item < thread_data->in_subgraphs->num_subgraphs)
    {
        subgraph_t *parent_subgraph = &in_subgraphs->subgraphs[current_item];

        rightmost_path_extensions (config, database, extension_storage, parent_subgraph);

        write_subgraph (database, sig_database, parent_subgraph, thread_data);

        for (unsigned int j = 0; j < extension_storage->num_extensions; j++) {
            extension_t *extension = &extension_storage->extensions[j];

            if (extension->support < config->approx_min_support) {
                free_extension (extension);
                continue;
            }

            subgraph_t test_subgraph = new_subgraph_from_extension (parent_subgraph, extension);
            database_entry_t canonicality_db_entry = canonical_database_entry (&test_subgraph);
            database_t canonicality_db = {};
            canonicality_db.num_entries = 1;
            canonicality_db.entries = &canonicality_db_entry;

            if (!is_canonical (&test_subgraph,
                               &canonicality_db,
                               canon_embedding_store))
            {
                free_canonical_database_entry (&canonicality_db_entry);
                free_subgraph (&test_subgraph);
                continue;
            }

            int directly_supported = test_subgraph.support >= config->direct_min_support;

            update_subgraph_isomorphisms (config,
                                          database,
                                          embedding_store,
                                          &test_subgraph,
                                          parent_subgraph->db_mappings,
                                          parent_subgraph->support);

            calc_subgraph_score (config, &test_subgraph);

            // Subgraph is discarded if at least one parent is unsupported
            if (!parent_check (config,
                               database,
                               &test_subgraph,
                               &canonicality_db,
                               in_subgraphs->hash_map,
                               in_subgraphs->hash_map_size,
                               canon_embedding_store,
                               rem_edge_store,
                               directly_supported))
            {
                free_subgraph_mappings (embedding_store, &test_subgraph);
                free_canonical_database_entry (&canonicality_db_entry);
                free_subgraph (&test_subgraph);
                continue;
            }

            free_canonical_database_entry (&canonicality_db_entry);

            calc_subgraph_hash (&test_subgraph);

            if (sig_database) {
                update_sig_subgraph_isomorphisms (config,
                                                  sig_database,
                                                  sig_embedding_store,
                                                  &test_subgraph,
                                                  parent_subgraph->sig_db_mappings,
                                                  parent_subgraph->sig_support);
            }

            num_new_subgraphs++;
            if (num_new_subgraphs % 100 == 1) {
                new_subgraphs = realloc (new_subgraphs, (num_new_subgraphs + 99) * sizeof(new_subgraphs[0]));
            }
            new_subgraphs[num_new_subgraphs - 1] = test_subgraph;
        }
        free_subgraph_labels (parent_subgraph);
        free_subgraph_rightmost_path (parent_subgraph);
        free_subgraph_parents (parent_subgraph);
        free_subgraph_parent_ids (parent_subgraph);
        free_subgraph_mappings (embedding_store, parent_subgraph);
        if (sig_database) {
            free_subgraph_sig_mappings (sig_embedding_store, parent_subgraph);
        }
        // Can't free the DFS code yet, as it is still needed for parent detection
    }

    out_subgraphs->subgraphs = new_subgraphs;
    out_subgraphs->num_subgraphs = num_new_subgraphs;

    return NULL;
}

void gspan_step (database_t        *database,
                 database_t        *sig_database,
                 unsigned int       num_threads,
                 thread_data_t     *thread_data,
                 subgraph_store_t  *in_subgraphs,
                 subgraph_store_t  *out_subgraphs,
                 unsigned int      *output_stats,
                 subgraph_t       **out_best_score_subgraph,
                 output_t          *output)
{
    unsigned int num_in_subgraphs = in_subgraphs->num_subgraphs;

    // Subgraphs can not be split across multiple threads
    num_threads = MIN(num_threads, num_in_subgraphs);

    pthread_t *threads = malloc (num_threads * sizeof(threads[0]));
    volatile unsigned int current_item = 0;

    for (unsigned int i = 0; i < num_threads; i++) {
        thread_data[i].current_item = &current_item;
        thread_data[i].item_count = num_in_subgraphs;
        thread_data[i].in_subgraphs = in_subgraphs;
    }

    mark_output_start (output);

    for (unsigned int i = 0; i < num_threads; i++) {
        pthread_create (&threads[i], NULL, gspan_step_helper, &thread_data[i]);
    }

    unsigned int num_new_subgraphs = 0;

    for (unsigned int i = 0; i < num_threads; i++) {
        pthread_join (threads[i], NULL);
        num_new_subgraphs += thread_data[i].out_subgraphs.num_subgraphs;
    }

    subgraph_t *new_subgraphs = malloc (num_new_subgraphs * sizeof(new_subgraphs[0]));

    unsigned int copied = 0;

    for (unsigned int i = 0; i < num_threads; i++) {
        unsigned int size = thread_data[i].out_subgraphs.num_subgraphs;
        memcpy (new_subgraphs + copied,
                thread_data[i].out_subgraphs.subgraphs,
                size * sizeof(new_subgraphs[0]));
        free (thread_data[i].out_subgraphs.subgraphs);
        copied += size;
    }

    memset(output_stats, 0, GRAPH_TYPE_COUNT * sizeof(output_stats[0]));

    int best_score = INT_MIN;
    subgraph_t *best_score_subgraph;

    unsigned int subgraph_id = in_subgraphs->subgraphs[num_in_subgraphs - 1].id + 1;
    for (unsigned int i = 0; i < num_new_subgraphs; i++) {
        new_subgraphs[i].id = subgraph_id++;
        output_stats[new_subgraphs[i].graph_type]++;
        if (new_subgraphs[i].score > best_score) {
            best_score = new_subgraphs[i].score;
            best_score_subgraph = &new_subgraphs[i];
        }
    }

    unsigned int new_subgraphs_hash_map_size = MAX(1024, num_new_subgraphs / 3);
    subgraph_t **new_subgraphs_hash_map = calloc (new_subgraphs_hash_map_size, sizeof(new_subgraphs_hash_map[0]));
    for (unsigned int i = 0; i < num_new_subgraphs; i++) {
        unsigned int hash = new_subgraphs[i].hash;
        unsigned int slot = hash % new_subgraphs_hash_map_size;
        subgraph_t *prev_subgraph = new_subgraphs_hash_map[slot];
        new_subgraphs_hash_map[slot] = &new_subgraphs[i];
        new_subgraphs[i].next = prev_subgraph;
    }

    free (threads);

    out_subgraphs->subgraphs = new_subgraphs;
    out_subgraphs->num_subgraphs = num_new_subgraphs;
    out_subgraphs->hash_map = new_subgraphs_hash_map;
    out_subgraphs->hash_map_size = new_subgraphs_hash_map_size;

    *out_best_score_subgraph = best_score_subgraph;
}

void *gspan_init_helper (void *data)
{
    thread_data_t *thread_data = data;
    config_t *config = thread_data->config;
    database_t *database = thread_data->database;
    database_t *sig_database = thread_data->sig_database;
    extension_store_t *extension_store = thread_data->shared_extension_store;
    embedding_store_t *embedding_store = &thread_data->embedding_store;
    embedding_store_t *sig_embedding_store = &thread_data->sig_embedding_store;
    subgraph_store_t *out_subgraphs = &thread_data->out_subgraphs;

    unsigned int subgraphs_alloc = 0;
    unsigned int num_subgraphs = 0;
    subgraph_t *subgraphs = NULL;

    // Assume empty graph as parent
    unsigned int min_support = MIN(database->num_entries - config->approx_max_support_difference, config->direct_min_support);
    unsigned int min_score = config->min_score_difference;

    unsigned int current_item;

    while (current_item = __sync_fetch_and_add (thread_data->current_item, 1),
           current_item < thread_data->item_count)
    {
        extension_t *extension = extension_store->extensions + current_item;

        // Determine support of the extension
        for (unsigned int db_idx = 0; db_idx < database->num_entries; db_idx++) {
            unsigned int remaining_db_entries = database->num_entries - db_idx;

            // Already counted via hash collision
            if (extension->supported_by[db_idx]) {
                continue;
            }

            // Extension can't possibly be supported if this is the case, so skip any further checks
            if (extension->support + remaining_db_entries < config->approx_min_support) {
                continue;
            }

            database_entry_t *current_db_entry = database->entries + db_idx;
            unsigned int num_db_vertices = current_db_entry->num_vertices;
            vertex_label_t *db_labels = current_db_entry->vertex_labels;

            vertex_label_t *vertex_u_label = &extension->code.from_label;
            vertex_label_t *vertex_v_label = &extension->code.to_label;
            edge_label_t *edge_label = &extension->code.edge_label;
            direction_t direction = extension->code.direction;

            int found = 0;

            for (vertex_t isomorphic_vertex_u = 0;
                 isomorphic_vertex_u < num_db_vertices;
                 isomorphic_vertex_u++)
            {
                neighbor_t *neighbors = database->entries[db_idx].neighbors[isomorphic_vertex_u];
                unsigned int neighbors_len = database->entries[db_idx].neighbors_len[isomorphic_vertex_u];

                mapping_info_t mapping_info = (mapping_info_t) {};

                if (!vertex_label_match (config,
                                         &mapping_info,
                                         &db_labels[isomorphic_vertex_u],
                                         vertex_u_label,
                                         &mapping_info))
                {
                    continue;
                }

                for (unsigned int k = 0; k < neighbors_len; k++) {
                    neighbor_t *neighbor = &neighbors[k];

                    if (direction != neighbor->direction){
                        continue;
                    }

                    if (!vertex_label_match (config,
                                             &mapping_info,
                                             &db_labels[neighbor->vertex_id],
                                             vertex_v_label,
                                             NULL))
                    {
                        continue;
                    }

                    if (!edge_label_match (config,
                                           &mapping_info,
                                           &neighbor->edge_label,
                                           edge_label,
                                           NULL))
                    {
                        continue;
                    }

                    found = 1;
                    break;
                }

                if (found) {
                    break;
                }
            }

            if (found) {
                extension->supported_by[db_idx] = 1;
                extension->support++;
            }
        }

        if (extension->support < min_support) {
            free_extension (extension);
            continue;
        }

        // Generate a new initial subgraph from the extension
        num_subgraphs++;
        if (num_subgraphs > subgraphs_alloc) {
            subgraphs_alloc += 2048;
            subgraphs = realloc (subgraphs, subgraphs_alloc * sizeof(subgraphs[0]));
        }

        subgraph_t *current_subgraph = subgraphs + num_subgraphs - 1;
        *current_subgraph = (subgraph_t) {};

        current_subgraph->vertex_labels = malloc (2 * sizeof(current_subgraph->vertex_labels[0]));
        current_subgraph->vertex_labels[0] = extension->code.from_label;
        current_subgraph->vertex_labels[1] = extension->code.to_label;

        current_subgraph->rightmost_path = malloc (2 * sizeof(current_subgraph->rightmost_path[0]));
        current_subgraph->rightmost_path[0] = 0;
        current_subgraph->rightmost_path[1] = 1;

        current_subgraph->rightmost_path_len = 2;

        current_subgraph->rightmost_vertex = 1;

        current_subgraph->dfs_code_len = 1;
        current_subgraph->dfs_code = malloc (1 * sizeof(current_subgraph->dfs_code[0]));
        current_subgraph->dfs_code[0] = extension->code;

        current_subgraph->support = extension->support;
        current_subgraph->supported_by = extension->supported_by;

        current_subgraph->graph_type = PATH;

        // Map the first two vertices
        init_subgraph_isomorphisms (config, database, embedding_store, current_subgraph);

        calc_subgraph_score (config, current_subgraph);

        if (current_subgraph->score < min_score) {
            num_subgraphs--;
            free_subgraph_mappings (embedding_store, current_subgraph);
            free_subgraph (current_subgraph);
            continue;
        }

        calc_subgraph_hash (current_subgraph);

        if (sig_database) {
            init_sig_subgraph_isomorphisms (config, sig_database, sig_embedding_store, current_subgraph);

        }
    }

    out_subgraphs->num_subgraphs = num_subgraphs;
    out_subgraphs->subgraphs = subgraphs;

    return NULL;
}

void gspan_init (database_t        *database,
                 database_t        *sig_database,
                 unsigned int       num_threads,
                 thread_data_t     *thread_data,
                 subgraph_store_t  *out_subgraphs,
                 unsigned int      *output_stats,
                 subgraph_t       **out_best_score_subgraph)
{
    dfs_code_t code;
    code.from = 0;
    code.to = 1;

    // For initialization the number of extensions does not depend on the number
    // of mappings/databases, so use 1 here.
    extension_store_t *extension_store = thread_data[0].shared_extension_store;
    prepare_extension_store (extension_store, 1, 2);

    for (unsigned int i = 0; i < database->num_entries; i++) {
        database_entry_t * current_db_entry = database->entries + i;
        for (unsigned int j = 0; j < current_db_entry->num_vertices; j++) {
            code.from_label = current_db_entry->vertex_labels[j];
            for (unsigned int k = 0; k < current_db_entry->neighbors_len[j]; k++) {
                neighbor_t *neighbor = &current_db_entry->neighbors[j][k];
                int comparison;
                code.to_label = current_db_entry->vertex_labels[neighbor->vertex_id];
                comparison = vertex_label_compare (&code.from_label, &code.to_label);

                // only consider canonical DFS codes
                if (comparison > 0) {
                    continue;
                }
                if (comparison == 0) {
                    if (neighbor->direction > invert_direction (neighbor->direction)) {
                        continue;
                    }
                }

                code.edge_label = neighbor->edge_label;
                code.direction = neighbor->direction;

                update_extensions (database, i, extension_store, &code);
            }
        }
    }

    num_threads = MIN(num_threads, extension_store->num_extensions);
    pthread_t *threads = malloc (num_threads * sizeof(threads[0]));
    volatile unsigned int current_item = 0;

    for (unsigned int i = 0; i < num_threads; i++) {
        thread_data[i].current_item = &current_item;
        thread_data[i].item_count = extension_store->num_extensions;
    }

    for (unsigned int i = 0; i < num_threads; i++) {
        pthread_create (&threads[i], NULL, gspan_init_helper, &thread_data[i]);
    }

    unsigned int num_new_subgraphs = 0;

    for (unsigned int i = 0; i < num_threads; i++) {
        pthread_join (threads[i], NULL);
        num_new_subgraphs += thread_data[i].out_subgraphs.num_subgraphs;
    }

    subgraph_t *new_subgraphs = malloc (num_new_subgraphs * sizeof(new_subgraphs[0]));
    unsigned int copied = 0;

    for (unsigned int i = 0; i < num_threads; i++) {
        unsigned int size = thread_data[i].out_subgraphs.num_subgraphs;
        memcpy (new_subgraphs + copied,
                thread_data[i].out_subgraphs.subgraphs,
                size * sizeof(new_subgraphs[0]));
        free (thread_data[i].out_subgraphs.subgraphs);
        copied += size;
    }

    memset(output_stats, 0, GRAPH_TYPE_COUNT * sizeof(output_stats[0]));
    output_stats[PATH] = num_new_subgraphs;

    int best_score = INT_MIN;
    subgraph_t *best_score_subgraph;
    unsigned int new_subgraphs_hash_map_size = MAX(1024, num_new_subgraphs / 3);
    subgraph_t **new_subgraphs_hash_map = calloc (new_subgraphs_hash_map_size, sizeof(new_subgraphs_hash_map[0]));

    unsigned int subgraph_id = 0;
    for (unsigned int i = 0; i < num_new_subgraphs; i++) {
        unsigned int hash = new_subgraphs[i].hash;
        unsigned int slot = hash % new_subgraphs_hash_map_size;
        subgraph_t *prev_subgraph = new_subgraphs_hash_map[slot];

        new_subgraphs[i].id = subgraph_id++;

        if (new_subgraphs[i].score > best_score) {
            best_score = new_subgraphs[i].score;
            best_score_subgraph = &new_subgraphs[i];
        }

        new_subgraphs_hash_map[slot] = &new_subgraphs[i];
        new_subgraphs[i].next = prev_subgraph;
    }

    free (threads);

    out_subgraphs->num_subgraphs = num_new_subgraphs;
    out_subgraphs->subgraphs = new_subgraphs;
    out_subgraphs->hash_map = new_subgraphs_hash_map;
    out_subgraphs->hash_map_size = new_subgraphs_hash_map_size;

    *out_best_score_subgraph = best_score_subgraph;
}

int database_entry_vertex_size_compare (const void *a,
                                      const void *b)
{
    const database_entry_t *entry_a = a;
    const database_entry_t *entry_b = b;
    return entry_a->num_vertices - entry_b->num_vertices;
}

int vertex_degree_compare_desc (const void *a,
                              const void *b)
{
    return *(int *)b - *(int *)a;
}

void read_database_file (FILE       *db_file,
                         database_t *database,
                         int         read_type_labels,
                         int         read_distance_labels,
                         int         read_direction,
                         int         include_comments,
                         char       *good_classes_str,
                         output_t   *output)
{
    char *line;
    char *line_buffer = NULL;
    size_t line_length = 0;
    size_t line_buffer_size = 0;
    size_t line_bytes_read = 0;
    const size_t line_buffer_increment = 4096;

    unsigned int line_count = 0;

    char line_type = 0;
    int mode = -1;

    int parsed = 0;
    int parsed_length;

    unsigned int vid = 0;
    vertex_label_t vlabel = 0;
    unsigned int efrom = 0;
    unsigned int eto = 0;
    edge_label_t elabel = {};
    direction_t direction;

    database_entry_t *current_entry = NULL;

    int classification = 0;
    vertex_label_t max_vertex_label = 0;
    unsigned int highest_vertex_id;
    unsigned int num_vertices = 0;
    unsigned int num_edges = 0;
    unsigned int vertex_alloc_size = 0;
    unsigned int old_vertex_alloc_size = 0;
    unsigned int neighbors_alloc_size = 0;
    neighbor_t **neighbors = NULL;
    neighbor_t *neighbors_alloc = NULL;
    unsigned int *neighbors_len = NULL;
    vertex_label_t *vertex_labels = NULL;

    database->num_classified_good = 0;
    database->num_entries = 0;
    database->entries = NULL;

    const char *class_delimiters = ":;, ";
    char *good_classes_cpy = NULL;
    char **good_classes = NULL;
    unsigned int num_good_classes = 0;
    if (good_classes_str) {
        good_classes_cpy = strdup(good_classes_str);
        char *good_class = strtok(good_classes_cpy, class_delimiters);
        while(good_class) {
            num_good_classes++;
            good_classes = realloc(good_classes, num_good_classes * sizeof(good_classes[0]));
            good_classes[num_good_classes - 1] = good_class;
            good_class = strtok(NULL, class_delimiters);
        }
    }

    while (!feof (db_file)) {
        line_length = 0;
        do {
            if (line_length + line_buffer_increment > line_buffer_size) {
                line_buffer_size += line_buffer_increment;
                line_buffer = realloc (line_buffer, line_buffer_size);
            }

            fgets (line_buffer + line_length, line_buffer_increment, db_file);
            line_bytes_read = strlen(line_buffer + line_length);
            line_length += line_bytes_read;
        } while (line_bytes_read > 0 && line_buffer[line_length - 1] != '\n');

        // Fix Windows encoded line endings on non-Windows machines
        if (line_length >= 2 && line_buffer[line_length - 2] == '\r' && line_buffer[line_length - 1] == '\n') {
            line_buffer[line_length - 2] = '\n';
            line_buffer[line_length - 1] = '\0';
            line_length--;
        }

        line = line_buffer;
        line_type = line[0];
        line_count++;

        switch (line_type) {
        case 't':
            mode = 0;
            database->num_entries++;
            classification = 0;
            highest_vertex_id = 0;
            num_vertices = 0;
            num_edges = 0;
            break;
        case 'v':
            if (!(mode == 0 || mode == 1)) {
                fprintf (stderr, "Missing transaction before vertex section. Line: %u\n", line_count);
                exit (1);
            }
                               
            mode = 1;

            parsed = sscanf (line, "v %u %d", &vid, &vlabel);

            if (parsed < 2) {
                fprintf (stderr, "Error parsing vertex. Line: %u\n", line_count);
                exit (1);
            }

            if (vid > VERTEX_MAX) {
                fprintf (stderr, "Vertex ID too large. Line: %u\n", line_count);
                exit (1);
            }

            num_vertices++;

            if (vid >= highest_vertex_id) {
                highest_vertex_id = vid;

                if (vertex_alloc_size <= highest_vertex_id) {
                    vertex_alloc_size = highest_vertex_id + 1000;
                    vertex_labels = realloc (vertex_labels, vertex_alloc_size * sizeof(vertex_labels[0]));
                }
            }

            vertex_labels[vid] = vlabel;

            if (vlabel > max_vertex_label) {
                max_vertex_label = vlabel;
            }
            break;
        case 'e':
            if (!(mode == 1 || mode == 2)) {
                fprintf (stderr, "Missing vertex section before edge section. Line: %u\n", line_count);
                exit (1);
            }
            
            mode = 2;

            parsed = sscanf (line,
                             "e %u %u %d%n",
                             &efrom,
                             &eto,
                             &elabel.interaction_type,
                             &parsed_length);
            line += parsed_length;

            if (parsed < 3) {
                fprintf (stderr, "Error parsing edge. Line: %u\n", line_count);
                exit (1);
            }

            if (efrom >= num_vertices) {
                fprintf (stderr, "Unspecified vertex %u used in graph %u. Line: %u\n", efrom, database->num_entries - 1, line_count);
                exit (1);
            }

            if (eto >= num_vertices) {
                fprintf (stderr, "Unspecified vertex %u used in graph %u. Line: %u\n", eto, database->num_entries - 1, line_count);
                exit (1);
            }

            if (efrom == eto) {
                fprintf (stderr, "Loop on vertex %u in graph %u. Line: %u\n", efrom, database->num_entries - 1, line_count);
                exit (1);
            }

            if (!read_type_labels) {
                elabel.interaction_type = 0;
            }

            elabel.distance = 0;

            if (read_distance_labels) {
                parsed = sscanf (line,
                                 "%d %n",
                                 &elabel.distance,
                                 &parsed_length);
                line += parsed_length;

                if (parsed < 1) {
                    fprintf (stderr, "Distance labels requested but none given in edge (%u, %u) in graph %u. "
                                     "Falling back to 0. Line: %u\n", efrom, eto, database->num_entries - 1, line_count);
                }
            }

            direction = UNDIRECTED;

            if (read_direction) {
                parsed = sscanf (line,
                                 " %u%n",
                                 &direction,
                                 &parsed_length);
                line += parsed_length;

                if (parsed < 1) {
                    fprintf (stderr, "Direction labels requested but none given in edge (%u, %u) in graph %u. "
                                     "Falling back to undirected. Line: %u\n", efrom, eto, database->num_entries - 1, line_count);
                }
            }


            for (unsigned int i = 0; i < neighbors_len[efrom]; i++) {
                neighbor_t *neighbor = &neighbors[efrom][i];
                if (neighbor->vertex_id == eto
                    && edge_label_compare (&neighbor->edge_label, &elabel) == 0
                    && direction == neighbor->direction)
                {
                    fprintf (stderr, "Duplicate edge between vertices %u and %u in graph %u. Line: %u\n", efrom, eto, database->num_entries - 1, line_count);
                    exit (1);
                }
            }
            neighbors[efrom][neighbors_len[efrom]].vertex_id = eto;
            neighbors[efrom][neighbors_len[efrom]].edge_label = elabel;
            neighbors[efrom][neighbors_len[efrom]].direction = direction;

            neighbors[eto][neighbors_len[eto]].vertex_id = efrom;
            neighbors[eto][neighbors_len[eto]].edge_label = elabel;
            neighbors[eto][neighbors_len[eto]].direction = invert_direction(direction);

            neighbors_len[efrom]++;
            neighbors_len[eto]++;

            num_edges++;
            break;
        case '#':
            {
                const char *fam_tag = "# fam ";
                char fams[4096];
                if (good_classes && strncmp (line, fam_tag, strlen (fam_tag)) == 0) {
                    strncpy (fams, line + strlen (fam_tag), 4096);
                    // Remove \n
                    fams[strlen (fams) - 1] = 0;
                    char *fam = strtok(fams, class_delimiters);
                    while(fam) {
                        for (unsigned int good_classes_iter = 0;
                             good_classes_iter < num_good_classes;
                             good_classes_iter++)
                        {
                            if (strcmp(good_classes[good_classes_iter], fam) == 0) {
                                classification = 1;
                                database->num_classified_good++;
                                break;
                            }
                        }
                        if (classification) {
                            break;
                        }
                        fam = strtok(NULL, class_delimiters);
                    }
                }

                // Copy all comments to the beginning of the output file
                if (include_comments) {
                    write_header_entry (output, line, line_length);
                }
            }
            break;
        case '\n':
            break;
        case '\r':
            break;
        default:
            fprintf (stderr, "Unrecognized input sequence '%s'. Line: %u\n", line, line_count);
            exit (1);
        }

        // Mode switches depending on the next line type
        line_type = fgetc (db_file);
        ungetc(line_type, db_file);

        // We are done reading vertices, prepare edges
        if (mode == 1 && line_type != 'v' && line_type != '#') {
            if (highest_vertex_id + 1 != num_vertices) {
                fprintf (stderr, "Gap between vertex IDs in graph %u: %u vertices found and the highest ID is %u\n",
                         database->num_entries - 1, num_vertices, highest_vertex_id);
                exit (1);
            }

            // Pre-allocate all data structures related to edges
            size_t new_neighbors_alloc_size = num_vertices * num_vertices;
            if (new_neighbors_alloc_size > neighbors_alloc_size) {
                free (neighbors_alloc);
                neighbors_alloc_size = new_neighbors_alloc_size;
                neighbors_alloc = malloc (neighbors_alloc_size * sizeof(neighbors_alloc[0]));
            }

            if (vertex_alloc_size > old_vertex_alloc_size) {
                free (neighbors);
                free (neighbors_len);
                old_vertex_alloc_size = vertex_alloc_size;
                neighbors = malloc (vertex_alloc_size * sizeof(neighbors[0]));
                neighbors_len = malloc (vertex_alloc_size * sizeof(neighbors_len[0]));
            }
            memset (neighbors_len, 0, num_vertices * sizeof(neighbors_len[0]));

            size_t index = 0;
            for (unsigned int i = 0; i < num_vertices; i++) {
                neighbors[i] = &neighbors_alloc[index];
                index += num_vertices;
            }

            mode = 2;
        }

        // We are done reading edges, put everything into the database entry now
        else if (mode == 2 && line_type != 'e' && line_type != '#') {
            database->entries = realloc (database->entries, database->num_entries * sizeof(database->entries[0]));
            current_entry = database->entries + database->num_entries - 1;
            current_entry->id = database->num_entries - 1;
            current_entry->classification = classification;
            current_entry->num_vertices = num_vertices;
            current_entry->num_edges = num_edges;
            current_entry->vertex_labels = malloc (num_vertices * sizeof(current_entry->vertex_labels[0]));
            memcpy (current_entry->vertex_labels, vertex_labels, num_vertices * sizeof(current_entry->vertex_labels[0]));
            current_entry->neighbors_len = malloc (num_vertices * sizeof(current_entry->neighbors_len[0]));
            memcpy (current_entry->neighbors_len, neighbors_len, num_vertices * sizeof(current_entry->neighbors_len[0]));

            current_entry->neighbors = malloc (num_vertices * sizeof(current_entry->neighbors[0]));
            neighbor_t *entry_neighbors_alloc = malloc (num_edges * 2 * sizeof(entry_neighbors_alloc[0]));
            unsigned int index = 0;
            for (unsigned int i = 0; i < num_vertices; i++) {
                current_entry->neighbors[i] = &entry_neighbors_alloc[index];
                memcpy(&entry_neighbors_alloc[index], neighbors[i], neighbors_len[i] * sizeof(entry_neighbors_alloc[0]));
                index += neighbors_len[i];
            }

            mode = 0;
        }
    }

    free (vertex_labels);
    free (neighbors);
    free (neighbors_len);
    free (neighbors_alloc);
    free (good_classes);
    free (good_classes_cpy);
    free (line_buffer);

    qsort (database->entries, database->num_entries, sizeof(database->entries[0]), database_entry_vertex_size_compare);

    database->num_label_types = max_vertex_label + 1;
}

thread_data_t *setup_thread_data (unsigned int       num_threads,
                                  config_t          *config,
                                  database_t        *database,
                                  database_t        *sig_database,
                                  extension_store_t *shared_extension_store,
                                  output_t          *output,
                                  unsigned int       max_subgraph_vertices,
                                  unsigned int       num_initial_extensions,
                                  unsigned int      *max_extensions_for_length)
{
    thread_data_t *thread_data = malloc (num_threads * sizeof(thread_data[0]));

    for (unsigned int i = 0; i < num_threads; i++) {
        thread_data[i].config = config;
        thread_data[i].database = database;
        thread_data[i].sig_database = sig_database;
        init_embedding_store (&thread_data[i].embedding_store);
        if (sig_database) {
            init_embedding_store (&thread_data[i].sig_embedding_store);
        }
        init_embedding_store (&thread_data[i].canon_embedding_store);
        init_extension_store (&thread_data[i].extension_store, num_initial_extensions, max_extensions_for_length);
        thread_data[i].shared_extension_store = shared_extension_store;
        init_rem_edge_store (&thread_data[i].rem_edge_store, max_subgraph_vertices);
        thread_data[i].output = output;
        thread_data[i].output_buffer_size = 16 * 1024;
        thread_data[i].output_buffer = malloc (thread_data[i].output_buffer_size);
    }

    return thread_data;
}

void free_thread_data (thread_data_t *thread_data,
                       unsigned int   num_threads,
                       int            use_sig_database)
{
    for (unsigned int i = 0; i < num_threads; i++) {
        free_embedding_store (&thread_data[i].embedding_store);
        if (use_sig_database) {
            free_embedding_store (&thread_data[i].sig_embedding_store);
        }
        free_embedding_store (&thread_data[i].canon_embedding_store);
        free_extension_store (&thread_data[i].extension_store);
        free_rem_edge_store (&thread_data[i].rem_edge_store);
        free (thread_data[i].output_buffer);
    }

    free (thread_data);
}


void print_result (unsigned int  subgraph_size,
                   unsigned int  num_subgraphs,
                   unsigned int *output_stats,
                   unsigned int  total_num_subgraphs,
                   subgraph_t   *best_score_subgraph,
                   database_t   *sig_database)
{
    if (num_subgraphs <= 0) {
        return;
    }

    if (subgraph_size == 1) {
        printf ("\n");
        printf (" Size |    Found |    Paths    Trees   Cyclic |    Total | Best Score \n");
        printf ("----------------------------------------------------------------------\n");
    }

    printf (" %4u | %8u | %8u %8u %8u | %8u | %d (%u; %u",
            subgraph_size,
            num_subgraphs,
            output_stats[PATH], output_stats[TREE], output_stats[CYCLE],
            total_num_subgraphs,
            best_score_subgraph->score,
            best_score_subgraph->id,
            best_score_subgraph->support);

    if (sig_database) {
        unsigned int num_classified_bad = sig_database->num_entries - sig_database->num_classified_good;
        unsigned int support_bad = best_score_subgraph->sig_support - best_score_subgraph->sig_good_support;

        printf ("; %u/%u - %u/%u",
                best_score_subgraph->sig_good_support,
                sig_database->num_classified_good,
                support_bad,
                num_classified_bad);
    }

    printf (")\n");

    fflush (stdout);
}


void print_usage (char *executable)
{
    fprintf (stderr, "%s [OPTION] <direct support threshold> <input file> [<output file>]\n"
                     "Options:\n"
                     "  -t <number of threads>\n"
                     "  -o <minimum size of subgraphs in output>\n"
                     "  -f <maximum output file size>\n"
                     "  -l <significance database>\n"
                     "  -c <good classes>\n"
                     "  -r (use directions)\n"
                     "  -a <minimum approximate support threshold>\n"
                     "  -d <maximum approximate support distance>\n"
                     "  -i <maximum distance difference>\n"
                     "  -s <minimum score difference>\n"
                     "  -m <maximum mismatch score>\n"
                     "  -x <mismatch scoring matrix>\n"
                     "  -y <type labels 1 (on; default) or 0 (off)>\n",
                     executable);
}

int main (int   argc,
          char *argv[])
{
    char *db_path = NULL;
    char *sig_db_path = NULL;
    char *out_path = NULL;
    char *mismatch_matrix_file_path = NULL;
    FILE *sig_db_file;
    FILE *db_file;

    char *good_classes = NULL;

    int num_threads = 1;
    int read_type_labels = 1;
    int directions = 0;
    int direct_min_support = 0;
    int approx_min_support = 0;
    int approx_max_support_difference = 0;
    int max_distance_difference = 0;
    int min_score_difference = 0;
    float max_mismatch_score = -1.0f;

    int output_min_size = 1;
    long int output_max_file_size = 500l * 1000l * 1000l;

    // For parsing floats properly under all locales
    setlocale(LC_NUMERIC, "C");

    int opt;
    opterr = 0;
    while ((opt = getopt (argc, argv, "t:o:f:c:l:ra:d:i:s:m:x:y")) != -1) {
        switch (opt) {
            case 't':
                num_threads = atoi (optarg);
                break;
            case 'o':
                output_min_size = atoi (optarg);
                break;
            case 'f':
                output_max_file_size = atol (optarg);
                int len = strlen (optarg);
                for (int x = 0; x < len; x++) {
                    switch (optarg[x]) {
                    case 'G':
                    case 'g':
                        output_max_file_size *= 1000l;
                    /* Fallthrough */
                    case 'M':
                    case 'm':
                        output_max_file_size *= 1000l;
                    /* Fallthrough */
                    case 'K':
                    case 'k':
                        output_max_file_size *= 1000l;
                    }
                }
                break;
            case 'c':
                good_classes = optarg;
                break;
            case 'l':
                sig_db_path = optarg;
                break;
            case 'r':
                directions = 1;
                break;
            case 'a':
                approx_min_support = atoi (optarg);
                break;
            case 'd':
                approx_max_support_difference = atoi (optarg);
                break;
            case 'i':
                max_distance_difference = atoi (optarg);
                break;
            case 's':
                min_score_difference = atoi (optarg);
                break;
            case 'm':
                {
                    float temp = atof (optarg);
                    if (temp < 0.0f) {
                        max_mismatch_score = -2.0f;
                    }
                    else {
                        max_mismatch_score = temp;
                    }
                }
                break;
            case 'x':
                mismatch_matrix_file_path = optarg;
                break;
            case 'y':
                read_type_labels = atoi (optarg);
                break;
            case '?':
                print_usage (argv[0]);
                return 1;
            default:
                abort ();
          }
    }

    if (num_threads < 1) {
        fprintf (stderr, "Invalid number of threads: %d\n", num_threads);
        fprintf (stderr, "There has to be at least 1 thread\n");
        exit (1);
    }

    if (output_max_file_size < 0) {
        fprintf (stderr, "Invalid maximum output file size: %ld\n", output_max_file_size);
        fprintf (stderr, "The maximum output file size needs to be positive\n");
        exit (1);
    }

    if (output_min_size < 1) {
        fprintf (stderr, "Invalid minimum output size: %d\n", output_min_size);
        fprintf (stderr, "The output size has to be at least 1\n");
        exit (1);
    }

    if (max_distance_difference < 0) {
        fprintf (stderr, "Invalid maximum distance difference: %d\n", max_distance_difference);
        fprintf (stderr, "This value can not be negative\n");
        exit (1);
    }

    if (max_mismatch_score == -2.0f) {
        fprintf (stderr, "Invalid max mismatch score: %lf\n", max_mismatch_score);
        fprintf (stderr, "This score can not be negative\n");
        exit (1);
    }

    if (max_mismatch_score != -1.0f && mismatch_matrix_file_path == NULL) {
        fprintf (stderr, "Max mismatch score specified but no mismatch scoring matrix given\n");
        exit (1);
    }

    if (max_mismatch_score == -1.0f && mismatch_matrix_file_path != NULL) {
        fprintf (stderr, "Mismatch scoring matrix given but no max mismatch score specified\n");
        exit (1);
    }

    if (min_score_difference < 0) {
        fprintf (stderr, "Invalid minimum score difference: %d\n", min_score_difference);
        fprintf (stderr, "This value can not be negative\n");
        exit (1);
    }

    if (sig_db_path && !good_classes) {
        fprintf (stderr, "Significance database given but no good classes specified\n");
        exit (1);
    }

    if (!sig_db_path && good_classes) {
        fprintf (stderr, "Good classes specified but no significance database given\n");
        exit (1);
    }

    if (max_mismatch_score == -1.0f) {
        max_mismatch_score = 0.0f;
    }

    if (argc - optind < 2
        || argc - optind > 3) {
        print_usage (argv[0]);
        return 1;
    }

    direct_min_support = atoi (argv[optind]);
    db_path = argv[optind + 1];

    if (argc - optind == 3) {
        out_path = argv[optind + 2];
    }

    db_file = fopen (db_path, "r");
    if (!db_file) {
        fprintf (stderr, "Could not open file: %s\n", db_path);
        exit (1);
    }

    output_t output;
    init_output (&output,
                 out_path,
                 output_max_file_size,
                 output_min_size,
                 !!max_distance_difference,
                 directions);

    database_t database = {};
    read_database_file (db_file,
                        &database,
                        !!read_type_labels,
                        !!max_distance_difference,
                        directions,
                        1,
                        NULL,
                        &output);

    fclose (db_file);

    // negative or zero threshold values are shorthands relative to the corresponding maximum values
    if (direct_min_support < 1) {
        direct_min_support = database.num_entries + direct_min_support;
    }

    if (approx_min_support < 1) {
        approx_min_support = direct_min_support + approx_min_support;
    }

    unsigned int max_num_label_types = database.num_label_types;

    database_t sig_database = {};
    database_t *sig_database_ptr = NULL;

    if (sig_db_path && good_classes) {
        sig_db_file = fopen (sig_db_path, "r");
        if (!sig_db_file) {
            fprintf (stderr, "Could not open file: %s\n", sig_db_path);
            exit (1);
        }

        read_database_file (sig_db_file,
                            &sig_database,
                            !!read_type_labels,
                            !!max_distance_difference,
                            directions,
                            0,
                            good_classes,
                            NULL);

        fclose (sig_db_file);
        sig_database_ptr = &sig_database;

        max_num_label_types = MAX(max_num_label_types, sig_database.num_label_types);
    }

    if (direct_min_support < 1 || direct_min_support > database.num_entries) {
        fprintf (stderr, "Invalid minimum direct support threshold: %d\n", direct_min_support);
        exit (1);
    }

    if (approx_min_support < 1 || approx_min_support > database.num_entries) {
        fprintf (stderr, "Invalid minimum approximate support threshold: %d\n", approx_min_support);
        exit (1);
    }

    if (direct_min_support < approx_min_support) {
        fprintf (stderr, "Minimum direct support threshold smaller than minimum approximate support threshold\n");
        exit (1);
    }

    if (approx_min_support > direct_min_support - approx_max_support_difference) {
        fprintf (stderr, "Approximation distance (%u) is larger than the distance "
                         "between the direct support threshold (%u) and the "
                         "approximate support threshold (%u).\n",
                         approx_max_support_difference,
                         direct_min_support,
                         approx_min_support);
        exit (1);
    }

    printf ("Input file: %s\n", db_path);

    printf ("Number of database graphs: %d\n", database.num_entries);
    printf ("Minimum direct support threshold: %d\n", direct_min_support);

    if (approx_min_support < direct_min_support && approx_max_support_difference > 0) {
        printf ("Minimum approximate support threshold: %d\n", approx_min_support);
        printf ("Approximation distance: %d\n", approx_max_support_difference);
    }

    if (mismatch_matrix_file_path) {
        printf ("Mismatch matrix: %s\n", mismatch_matrix_file_path);
        printf ("Maximum mismatch score: %f\n", max_mismatch_score);
    }

    if (sig_db_path) {
        printf ("Significance database file: %s\n", sig_db_path);
        printf ("Good significance classes: %s\n", good_classes);
        printf ("Number of graphs classified as good: %d/%d\n",
                sig_database.num_classified_good, sig_database.num_entries);
    }

    if (out_path) {
        printf ("Output file: %s\n", out_path);
        // SI units are the best units.
        printf ("Maximum output file size: %ld MB (%ld Bytes)\n",
                output_max_file_size / (1000l * 1000l),
                output_max_file_size);
        if (output_min_size > 1) {
            printf ("Mininmum output size: %d\n", output_min_size);
        }
    }

    fflush(stdout);

    config_t config = {};
    config.direct_min_support = direct_min_support;
    config.approx_min_support = approx_min_support;
    config.approx_max_support_difference = approx_max_support_difference;
    config.max_mismatch_score = max_mismatch_score;
    config.max_distance_difference = max_distance_difference;
    config.min_score_difference = min_score_difference;

    if (mismatch_matrix_file_path) {
        config.mismatch_matrix = parse_mismatch_matrix (max_num_label_types, mismatch_matrix_file_path);
    }
    else {
        config.mismatch_matrix = generate_default_mismatch_matrix (max_num_label_types);
    }

    // Determine upper bound for the maximum number of extensions per mapping

    unsigned int max_num_vertices_db_entry_for_lowest_support = database.num_entries - direct_min_support;
    // This requires the database to be sorted by number of vertices (ascending)
    unsigned int max_subgraph_vertices = database.entries[max_num_vertices_db_entry_for_lowest_support].num_vertices;
    unsigned int *max_extensions_for_length = calloc (max_subgraph_vertices + 1, sizeof(max_extensions_for_length[0]));
    unsigned int num_initial_extensions = 0;

    for (unsigned int db_idx = 0; db_idx < database.num_entries; db_idx++) {
        database_entry_t *db_entry = &database.entries[db_idx];
        // The extensions of the empty graph can potentially include all edges
        // from all databases.
        num_initial_extensions += db_entry->num_edges;

        // In all other cases extensions are based on an existing graph. This
        // existing graph can allow at most as many extensions as the highest
        // sum of degrees of the vertices of a potential rightmost path in a
        // database graph minus the edges already included in the rightmost
        // path. Only consider the largest sum here, it will later be multiplied
        // by the number of mappings. This is just an upper bound!
        unsigned int *degrees = malloc (db_entry->num_vertices * sizeof(degrees[0]));
        assert(sizeof(degrees[0]) == sizeof(db_entry->neighbors_len[0]));
        memcpy (degrees, db_entry->neighbors_len, db_entry->num_vertices * sizeof(degrees[0]));
        qsort (degrees, db_entry->num_vertices, sizeof(degrees[0]), vertex_degree_compare_desc);
        for (unsigned int length = 2; length <= max_subgraph_vertices; length++) {
            if (length > db_entry->num_vertices) {
                break;
            }
            unsigned int degree_sum = 0;
            for (unsigned int i = 0; i < length; i++) {
                degree_sum += degrees[i];
            }

            // All internal vertices on the rightmost path are using two edges already
            // the two terminal vertices only one.
            unsigned int rightmost_path_edges = 2 * (length - 1);

            // On disconnected graphs with multiple connected components, the maximum
            // lenght for the rightmost path does not cover all vertices. In these
            // cases the maximum extensions for this graph at a length for the rightmost
            // path that is longer than its largest connected component is 0.
            if (degree_sum < rightmost_path_edges) {
                degree_sum = 0;
            }
            else {
                degree_sum -= rightmost_path_edges;
            }

            if (degree_sum > max_extensions_for_length[length]) {
                max_extensions_for_length[length] = degree_sum;
            }
        }
        free (degrees);
    }

    extension_store_t shared_extension_store;
    init_extension_store (&shared_extension_store,
                          num_initial_extensions,
                          max_extensions_for_length);

    thread_data_t *thread_data = setup_thread_data (num_threads,
                                                    &config,
                                                    &database,
                                                    sig_database_ptr,
                                                    &shared_extension_store,
                                                    &output,
                                                    max_subgraph_vertices,
                                                    num_initial_extensions,
                                                    max_extensions_for_length);

    subgraph_store_t prev_subgraphs;
    subgraph_store_t subgraphs;
    subgraph_t *best_score_subgraph;
    unsigned int output_stats[GRAPH_TYPE_COUNT];

    unsigned int total = 0;

    gspan_init (&database,
                sig_database_ptr,
                num_threads,
                thread_data,
                &prev_subgraphs,
                output_stats,
                &best_score_subgraph);

    free_extension_store (&shared_extension_store);

    total = prev_subgraphs.num_subgraphs;

    unsigned int subgraph_size = 1;

    print_result (subgraph_size,
                  prev_subgraphs.num_subgraphs,
                  output_stats,
                  total,
                  best_score_subgraph,
                  sig_database_ptr);

    while (prev_subgraphs.num_subgraphs) {
        subgraph_size++;

        gspan_step (&database,
                    sig_database_ptr,
                    num_threads,
                    thread_data,
                    &prev_subgraphs,
                    &subgraphs,
                    output_stats,
                    &best_score_subgraph,
                    &output);

        total += subgraphs.num_subgraphs;

        print_result (subgraph_size,
                      subgraphs.num_subgraphs,
                      output_stats,
                      total,
                      best_score_subgraph,
                      sig_database_ptr);

        free_subgraph_store (&prev_subgraphs);
        prev_subgraphs = subgraphs;
    }

    free_subgraph_store (&prev_subgraphs);

    free_thread_data (thread_data, num_threads, !!sig_database_ptr);

    free (max_extensions_for_length);

    printf ("\n");
    printf ("Found %d supported subgraphs\n", total);
    fflush(stdout);

    finalize_output (&output);

    free_database (&database);
    if (sig_database_ptr)
    {
        free_database (sig_database_ptr);
    }

    free_mismatch_matrix (config.mismatch_matrix);

    return 0;
}
