/**
 * @file SuffixArray.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <utility>
#include <vector>

#include "sais.h"

using ::std::cout;
using ::std::endl;
using ::std::make_pair;
using ::std::numeric_limits;
using ::std::pair;
using ::std::set;
using ::std::size_t;
using ::std::vector;

#include "Parameters.hpp"
#include "Sequence.hpp"
#include "SequenceDatabase.hpp"
#include "SuffixBuckets.hpp"
#include "SuffixArray.hpp"

#ifndef SIZE_MAX
#define SIZE_MAX (size_t(-1))
#endif

namespace pgraph {

const size_t SuffixArray::npos(-1);

static int sufcheck(const unsigned char *T, const int *SA, int n, int verbose) {
    int C[256];
    int i, p, q, t;
    int c;

    if(verbose) { fprintf(stderr, "sufcheck: "); }
    if(n == 0) {
        if(verbose) { fprintf(stderr, "Done.\n"); }
        return 0;
    }

    /* Check arguments. */
    if((T == NULL) || (SA == NULL) || (n < 0)) {
        if(verbose) { fprintf(stderr, "Invalid arguments.\n"); }
        return -1;
    }

    /* check range: [0..n-1] */
    for(i = 0; i < n; ++i) {
        if((SA[i] < 0) || (n <= SA[i])) {
            if(verbose) {
                fprintf(stderr, "Out of the range [0,%d].\n"
                        "  SA[%d]=%d\n",
                        n - 1, i, SA[i]);
            }
            return -2;
        }
    }

    /* check first characters. */
    for(i = 1; i < n; ++i) {
        if(T[SA[i - 1]] > T[SA[i]]) {
            if(verbose) {
                fprintf(stderr, "Suffixes in wrong order.\n"
                        "  T[SA[%d]=%d]=%d > T[SA[%d]=%d]=%d\n",
                        i - 1, SA[i - 1], T[SA[i - 1]], i, SA[i], T[SA[i]]);
            }
            return -3;
        }
    }

    /* check suffixes. */
    for(i = 0; i < 256; ++i) { C[i] = 0; }
    for(i = 0; i < n; ++i) { ++C[T[i]]; }
    for(i = 0, p = 0; i < 256; ++i) {
        t = C[i];
        C[i] = p;
        p += t;
    }

    q = C[T[n - 1]];
    C[T[n - 1]] += 1;
    for(i = 0; i < n; ++i) {
        p = SA[i];
        if(0 < p) {
            c = T[--p];
            t = C[c];
        } else {
            c = T[p = n - 1];
            t = q;
        }
        if((t < 0) || (p != SA[t])) {
            if(verbose) {
                fprintf(stderr, "Suffix in wrong position.\n"
                        "  SA[%d]=%d or\n"
                        "  SA[%d]=%d\n",
                        t, (0 <= t) ? SA[t] : -1, i, SA[i]);
            }
            return -4;
        }
        if(t != q) {
            ++C[c];
            if((n <= C[c]) || (T[SA[C[c]]] != c)) { C[c] = -1; }
        }
    }

    if(1 <= verbose) { fprintf(stderr, "Done.\n"); }
    return 0;
}

SuffixArray::SuffixArray(
        SequenceDatabase *sequences, Bucket *bucket, const Parameters &param, int k)
    :   sequences(sequences)
    ,   bucket(bucket)
    ,   param(param)
    ,   SIGMA(param.alphabet.size())
    ,   DOLLAR(param.alphabet_dollar)
    ,   BEGIN(param.alphabet_begin)
    ,   alphabet(param.alphabet)
    ,   alphabet_table(numeric_limits<unsigned char>::max(), npos)
    ,   window_size(k >= 0 ? k : param.window_size)
    ,   T(NULL)
    ,   SA(NULL)
    ,   LCP(NULL)
    ,   BWT(NULL)
    ,   SID(NULL)
    ,   POS(NULL)
    ,   END()
    ,   SIDmap()
    ,   n(0)
    ,   sid(0)
    ,   sentinal(0)
    ,   sequences_cache()
{
    /* convert the given bucket into a form we can process, keeping only the
     * longest suffix for each sequence represented */
    map<size_t,size_t> input_subset;
    for (Suffix *p=bucket->suffixes; p!=NULL; p=p->next) {
#if 1
        if (input_subset.count(p->sid)) {
            if (input_subset[p->sid] > p->pid) {
                input_subset[p->sid] = p->pid;
            }
        }
        else {
            input_subset[p->sid] = p->pid;
        }
#else
        input_subset[p->sid] = 0;
#endif
    }
    size_t uncompressed_size = 0;
    size_t compressed_size = 0;
    for (map<size_t,size_t>::iterator it=input_subset.begin();
            it!=input_subset.end(); ++it) {
        Sequence &sequence = get_sequence(it->first);
        uncompressed_size += sequence.size();
        compressed_size += sequence.size() - it->second;
    }
    n = compressed_size;
#if DEBUG
    cout << "uncompressed_size = " << uncompressed_size << endl;
    cout << "  compressed_size = " <<   compressed_size << endl;
#endif

    /* Allocate 9n bytes of memory. */
    T = (unsigned char *)malloc((size_t)(n+1) * sizeof(unsigned char));
    SA = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for computing LCP */
    LCP = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for lcp tree */
    BWT = (unsigned char *)malloc((size_t)(n+1) * sizeof(unsigned char));
    SID = (int *)malloc((size_t)n * sizeof(int));
    POS = (int *)malloc((size_t)n * sizeof(int));
    if((T == NULL)
            || (SA == NULL)
            || (LCP == NULL)
            || (BWT == NULL)
            || (SID == NULL)
            || (POS == NULL))
    {
        fprintf(stderr, "sais: Cannot allocate memory.\n");
        exit(EXIT_FAILURE);
    }

    bool uses_delim = (sequences->get_delimiter() != '\0');
    /* copy sequences to compressed input T */
    size_t offset = 0;
    for (map<size_t,size_t>::iterator it=input_subset.begin();
            it!=input_subset.end(); ++it) {
        Sequence &sequence = get_sequence(it->first);
        size_t length = sequence.size() - it->second;
        if (uses_delim) {
            length -= 1;
        }
        (void)memcpy(&T[offset], &sequence[it->second], length);
        offset += length;
        T[offset] = '$';
        offset += 1;
        SIDmap.push_back(it->first);
    }

    T[n]='\0'; /* so we can print it */
    //if (n < 256) printf("%s\n", T);
    //printf("%s\n", T);

    /* determine sentinal */
    sentinal = '$';
#if DEBUG
    printf("sentinal='%c'\n", sentinal);
#endif

#if DEBUG
    printf("original n=%d\n", n);
#endif
    /* determine actual end of file (last char) */
    {
        int off = 0;
        while (!isgraph(T[n-off])) {
            ++off;
        }
        n = n - off + 1;
    }
#if DEBUG
    printf("adjusted n=%d\n", n);
#endif

    /* scan T from left to right to build auxiliary info stuff... yep */
    sid = 0;
    int pos = 0;
    for (int i=0; i<n; ++i) {
        SID[i] = sid;
        POS[i] = pos++;
        if (T[i] == sentinal) {
            END.push_back(i);
            ++sid;
            pos = 0;
        }
    }
    if (0 == sid) { /* no sentinal found, that's okay */
#if DEBUG
        printf("no sentinal(%c) found in input\n", sentinal);
#endif
        END.push_back(n-1);
    }

    /* Construct the suffix array. */
#if DEBUG
    fprintf(stderr, "sais: %d bytes ... \n", n);
#endif
    double start = MPI_Wtime();
    if(sais(T, SA, LCP, (int)n) != 0) {
        fprintf(stderr, "sais: Cannot allocate memory.\n");
        exit(EXIT_FAILURE);
    }
    double finish = MPI_Wtime();
#if DEBUG
    fprintf(stderr, "induced SA: %.4f sec\n", finish - start);
#endif

    /* naive BWT: */
    /* also "fix" the LCP array to clamp LCP's that are too long */
    /* while we're at it, determine bucket dilenations */
    start = MPI_Wtime();
    for (int i = 0; i < n; ++i) {
        int len = END[SID[SA[i]]] - SA[i]; // don't include sentinal
        if (LCP[i] > len) LCP[i] = len;
        BWT[i] = (SA[i] > 0) ? T[SA[i]-1] : sentinal;
    }
    finish = MPI_Wtime();
#if DEBUG
    fprintf(stderr, " naive BWT: %.4f sec\n", finish - start);
#endif

    bool validate = false;
    if (validate) {
        /* check SA: */
        if (sufcheck(T, SA, (int)n, 1) != 0) { exit(EXIT_FAILURE); }

        /* check LCP: */
        fprintf(stderr, "LCP check: ");
        for (int i = 1; i < n; ++i) {
            int l = 0;
            while (T[SA[i]+l]==T[SA[i-1]+l]
                    && (END[SID[SA[i-1]]]-SA[i-1])>l) ++l;
            if (l != LCP[i]) {
                printf("Error at position %i\n", i);
                printf("%i vs. %i\n", l, LCP[i]);
                for (int j = 0; j < 10; j++) printf("%c", T[SA[i]+j]); printf("\n");
                for (int j = 0; j < 10; j++) printf("%c", T[SA[i-1]+j]); printf("\n");
                exit(-1);
            }
        }
        fprintf(stderr, "Done.\n");
    }
#if DEBUG
    else {
        printf("skipping validation\n");
    }
#endif
    LCP[n] = 0; /* doesn't really exist, but for the root */
}


SuffixArray::~SuffixArray()
{
    free(T);
    free(SA);
    free(LCP);
    free(BWT);
    free(SID);
    free(POS);
    for (map<size_t,Sequence*>::iterator it=sequences_cache.begin();
            it!=sequences_cache.end(); ++it) {
        delete it->second;
    }
}


void SuffixArray::print()
{
#define CUTOFF (20)
    printf("Index\tSA\tLCP\tBWT\tT\tSeqID\tSeqPos\tSuffix\n");
    for(int i=0; i<n; ++i) {
        int len = END[SID[SA[i]]] - SA[i] + 1;
        if (len > CUTOFF) len = CUTOFF;
        //int len = n - SA[i];
        printf("%d\t%d\t%d\t%c\t%c\t%d\t%d\t%.*s",
                i, SA[i], LCP[i], BWT[i], T[i], SID[SA[i]], POS[SA[i]], len, T+SA[i]);
        if (len == CUTOFF) printf("...");
        printf("\n");
    }
}

}; /* namespace pgraph */

