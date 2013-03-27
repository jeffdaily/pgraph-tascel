#include "uFind.h"

#pragma mta parallel on
struct ufind *init_union(int size) {
	int i;
	struct ufind *ufSet;

	ufSet = emalloc(size*sizeof(struct ufind));

	/* init all set elems */
	for (i=0; i<size; i++) {
		ufSet[i].parent = i;
		ufSet[i].rank = 0;
	}

	return ufSet;
}

#pragma mta parallel off
void free_union(struct ufind *uf){
    free(uf);
}

#pragma mta inline
int find(struct ufind *ufSet, int elemIndex) {
	int myParent;
	myParent = ufSet[elemIndex].parent;

	if (ufSet[myParent].parent != myParent) {
        /* path compression part */
		myParent = ufSet[elemIndex].parent = find(ufSet, myParent);

        /* elemIndex's parent is updated, so the rank should decrease */
		ufSet[myParent].rank--;
		return myParent;
	} else {
		return myParent;
	}
}


#pragma mta inline
void union_elems(struct ufind *ufSet, int elem1, int elem2) {
    merge_elems(ufSet, find(ufSet, elem1), find(ufSet, elem2));
}

/* can never been called outside of this file */
#pragma mta inline
void merge_elems(struct ufind *ufSet, int elem1, int elem2) {
	if (elem1!=elem2) {
		if (ufSet[elem1].rank > ufSet[elem2].rank) {
			ufSet[elem2].parent = elem1;
		} else {
			ufSet[elem1].parent = elem2;
			if (ufSet[elem1].rank==ufSet[elem2].rank) {
				ufSet[elem2].rank++;
			}
		}
	}
}

int disp_all_clusters(struct ufind *uf, int size, int *singletons, char *dir, SEQ *seqs) {
	int i=0, iClusters=0;
	struct clusters *clust;
	struct link *temp, *temp_fix;
	int root;
	char s[200];
	FILE *fp, *fp1, *tmp;

	clust = emalloc(sizeof(struct clusters)*size);

	for (i=0; i<size; i++) {
		clust[i].count=1;
		clust[i].singleton='1';
		clust[i].head=NULL;
	}

	for (i=0; i<size; i++) {
		root = find(uf, i);
		if (i!=root) {
			temp = emalloc(sizeof(struct link));

			temp->id=i;
			temp->next=clust[root].head;
			clust[root].head=temp;
			clust[root].count++;
			clust[root].singleton='0';
			clust[i].singleton = '0';
			clust[i].head = NULL;
		}
	} /* end for */

	sprintf(s, "%s/vertexClustSize.%d.SGL", dir, size);
    
    /* to detect existing file */
    while(1){
        if((tmp=fopen(s, "r"))) {
            strcat(s, "d");
            fclose(tmp);
        } else break;
    }
	fp1 = efopen(s, "w");

	sprintf(s, "%s/vertexClust.%d.SGL", dir, size);
    
    /* to detect existing file */
    while(1){
        if((tmp=fopen(s, "r"))) {
            strcat(s, "d");
            fclose(tmp);
        } else break;
    }
	fp = efopen(s, "w");

	iClusters=0;
	*singletons=0;
	for (i=0; i<size; i++) {
		if (clust[i].singleton=='0' && clust[i].head!=NULL) {

			fprintf(fp, "{Cluster#}   %d\n", iClusters);
			fprintf(fp, "{Member#}  %s\n", seqs[i].gid);
			temp = clust[i].head;
			while (temp!=NULL) {
				fprintf(fp, "{Member#}  %s\n", seqs[temp->id].gid);
				temp = temp->next;
			}
			fprintf(fp1, "%d %d\n", iClusters, clust[i].count);
			iClusters++;
		} else {
			if (clust[i].singleton=='1') {
				fprintf(fp, "{Cluster#}   %d\n", iClusters);
				fprintf(fp, "{Member#}  %s\n", seqs[i].gid);
				fprintf(fp1, "%d %d\n", iClusters, clust[i].count);

				iClusters++;
				(*singletons)++;
			}
		}
	} /* end for 2 */

	fclose(fp);
	fclose(fp1);

	for (i=0; i<size; i++) {
		if (clust[i].head) {
			temp = clust[i].head;
			while (temp->next != 0) {
				temp_fix = temp->next;
				free(temp);
				temp = temp_fix;
			}
			free(temp);
		}
	}
	if (clust) free(clust);

	return iClusters;

}

void disp(struct ufind *ufSet, int size) {
	int i;
	for (i=0; i<size; i++) {
		if (i!=ufSet[i].parent)
			printf("(%d,%d) ", i, ufSet[i].parent);
		fflush(stdout);
	}
	printf("\n");
}
