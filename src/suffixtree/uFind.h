#ifndef UFIND_H_
#define UFIND_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "type.h"
#include "elib.h"

typedef struct ufind {
  int parent;
  int rank;
}UF;

struct link {
  int id;
  struct link *next;
};

struct clusters {
  int count;
  char singleton;      /* 1 if singleton 0 if no */
  struct link *head;
};

extern char **n2gidHash; /* i -> gid */

struct ufind *init_union(int size);
void free_union(struct ufind *uf);
int find(struct ufind *ufSet,int inputIndex);
void union_elems(struct ufind *ufSet,int elem1,int elem2);
void merge_elems(struct ufind *ufSet,int elem1,int elem2);

void disp(struct ufind *ufSet,int size);
int disp_all_clusters(struct ufind *uf, int size, int *singletons, char *dir, SEQ *seqs);

#endif /* end of uFind.h */
