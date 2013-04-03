/**
 * @file uFind.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef UFIND_H_
#define UFIND_H_

#include "loadseq.h"

typedef struct {
    int parent;
    int rank;
} ufind_t;

typedef struct link {
    int id;
    struct link *next;
} link_t;

typedef struct {
    int count;
    char singleton;      /* 1 if singleton 0 if no */
    struct link *head;
} cluster_t;

extern char **n2gidHash; /* i -> gid */

ufind_t *init_union(int size);
void free_union(ufind_t *uf);
int find(ufind_t *ufSet, int inputIndex);
void union_elems(ufind_t *ufSet, int elem1, int elem2);

void disp(ufind_t *ufSet, int size);
int disp_all_clusters(ufind_t *uf, int size, int *singletons, char *dir, sequence_t *seqs);

#endif /* end of uFind.h */
