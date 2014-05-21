/*
 * suftest.cpp for sais-lite
 * Copyright (c) 2008-2010 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
#include "config.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <iostream>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <stack>
#include <utility>
#include <vector>

#include "sais.h"

using ::std::cout;
using ::std::endl;
using ::std::list;
using ::std::make_pair;
using ::std::map;
using ::std::ostream;
using ::std::pair;
using ::std::set;
using ::std::size_t;
using ::std::stack;
using ::std::vector;

typedef map<char,list<int> > pset_t;

struct quad {
    int lcp;
    int lb;
    int rb;
    vector<quad> children;

    quad()
        : lcp(0), lb(0), rb(INT_MAX), children() {}
    quad(int lcp, int lb, int rb)
        : lcp(lcp), lb(lb), rb(rb), children() {}
    quad(int lcp, int lb, int rb, vector<quad> children)
        : lcp(lcp), lb(lb), rb(rb), children(children) {}

    bool empty() { return rb == INT_MAX; }
};

ostream& operator << (ostream &os, const quad &q) {
    os << "quad("
        << q.lcp
        << "," << q.lb
        << "," << q.rb
#if 0
        << ",{";
    if (q.children.size() > 0) {
        os << q.children[0];
        for (size_t i=1; i<q.children.size(); ++i) {
            os << "," << q.children[i];
        }
    }
    os << "})";
#else
       << "," << q.children.empty()
       << ")";
#endif
    return os;
}

static int count;
static int count_generated;
static set<pair<int,int> > pairs;

inline static void pair_check(
        const int &i,
        const int &j,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const char &sentinal);

inline static void process(const quad &q);

inline static void process(
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const char &sentinal,
        const int &cutoff);


/* Checks the suffix array SA of the string T. */
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

static void print_help(const char *progname, int status) {
    fprintf(stderr, "usage: %s [-p] [-c cutoff] [-s sentinal] FILE\n\n", progname);
    exit(status);
}

int main(int argc, const char *argv[]) {
    FILE *fp = NULL;
    const char *fname = NULL;
    unsigned char *T = NULL;
    int *SA = NULL;
    int *LCP = NULL;
    unsigned char *BWT = NULL;
    int *SID = NULL;
    int *END = NULL;
    int n = 0;
    clock_t start = 0;
    clock_t finish = 0;
    int i = 0;
    int j = 0;
    int l = 0;
    int sid = 0;
    char sentinal = 0;
    int print = 0;
    int validate = 0;
    int lb = 0;
    int cutoff = 1;

    /* Check arguments. */
    i = 1;
    while (i < argc) {
        if ((strncmp(argv[i], "-h", 2) == 0)
                || (strncmp(argv[i], "--help", 6) == 0)) {
            fprintf(stderr, "help requested\n");
            print_help(argv[0], EXIT_SUCCESS);
        }
        else if (strncmp(argv[i], "-c", 2) == 0) {
            if (i+1 < argc && strlen(argv[i+1]) == 1) {
                cutoff = atoi(argv[i+1]);
                ++i;
            }
            else {
                if (i+1 >= argc) {
                    fprintf(stderr, "-s takes a parameter\n");
                }
                else {
                    fprintf(stderr, "bad argument to -s: %s\n", argv[i+1]);
                }
                print_help(argv[0], EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "-s", 2) == 0) {
            if (i+1 < argc && strlen(argv[i+1]) == 1) {
                sentinal = argv[i+1][0];
                ++i;
            }
            else {
                if (i+1 >= argc) {
                    fprintf(stderr, "-s takes a parameter\n");
                }
                else {
                    fprintf(stderr, "bad argument to -s: %s\n", argv[i+1]);
                }
                print_help(argv[0], EXIT_FAILURE);
            }
        }
        else if (strncmp(argv[i], "-p", 2) == 0) {
            print = 1;
        }
        else if (strncmp(argv[i], "-x", 2) == 0) {
            validate = 1;
        }
        else if (strncmp(argv[i], "-", 1) != 0) {
            /* filename */
            if (i+1 != argc) { /* last argument */
                print_help(argv[0], EXIT_FAILURE);
            }
            fname = argv[i];
        }
        else {
            print_help(argv[0], EXIT_FAILURE);
        }
        ++i;
    }
    if (fname == NULL) {
        fprintf(stderr, "missing input file\n");
        print_help(argv[0], EXIT_FAILURE);
    }

    /* Open a file for reading. */
    if((fp = fopen(fname, "rb")) == NULL) {
        fprintf(stderr, "%s: Cannot open file `%s': ", argv[0], fname);
        perror(NULL);
        exit(EXIT_FAILURE);
    }

    /* Get the file size. */
    if(fseek(fp, 0, SEEK_END) == 0) {
        n = ftell(fp);
        rewind(fp);
        if(n < 0) {
            fprintf(stderr, "%s: Cannot ftell `%s': ", argv[0], fname);
            perror(NULL);
            exit(EXIT_FAILURE);
        }
    } else {
        fprintf(stderr, "%s: Cannot fseek `%s': ", argv[0], fname);
        perror(NULL);
        exit(EXIT_FAILURE);
    }

    /* Allocate 9n bytes of memory. */
    T = (unsigned char *)malloc((size_t)(n+1) * sizeof(unsigned char));
    SA = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for computing LCP */
    LCP = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for lcp tree */
    BWT = (unsigned char *)malloc((size_t)(n+1) * sizeof(unsigned char));
    SID = (int *)malloc((size_t)n * sizeof(int));
    if (print) {
        END = (int *)malloc((size_t)n * sizeof(int));
    }
    if((T == NULL)
            || (SA == NULL)
            || (LCP == NULL)
            || (BWT == NULL)
            || (SID == NULL)
            || (END == NULL && print))
    {
        fprintf(stderr, "%s: Cannot allocate memory.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* Read n bytes of data. */
    if(fread(T, sizeof(unsigned char), (size_t)n, fp) != (size_t)n) {
        fprintf(stderr, "%s: %s `%s': ",
                argv[0],
                (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
                fname);
        perror(NULL);
        exit(EXIT_FAILURE);
    }
    fclose(fp);

    T[n]='\0'; /* so we can print it */
    //if (n < 256) printf("%s\n", T);
    if (n < 256) printf("%s\n", T);

    /* determine sentinal */
    if (sentinal == 0) {
        int off = 0;
        while (!isgraph(T[n-off])) {
            ++off;
        }
        sentinal = T[n-off];
    }
    printf("sentinal='%c'\n", sentinal);

    printf("original n=%d\n", n);
    /* determine actual end of file (last char) */
    {
        int off = 0;
        while (!isgraph(T[n-off])) {
            ++off;
        }
        n = n - off + 1;
    }
    printf("adjusted n=%d\n", n);

    /* scan T from left to right to build auxiliary info stuff... yep */
    if (print) {
        sid = 0;
        for (i=0; i<n; ++i) {
            SID[i] = sid;
            if (T[i] == sentinal) {
                END[sid] = i;
                ++sid;
            }
        }
    }
    else {
        sid = 0;
        for (i=0; i<n; ++i) {
            SID[i] = sid;
            if (T[i] == sentinal) {
                ++sid;
            }
        }
    }
    if (0 == sid) { /* no sentinal found, that's okay */
        printf("no sentinal(%c) found in input\n", sentinal);
        if (print) END[0] = n-1;
    }

    /* Construct the suffix array. */
    fprintf(stderr, "%s: %d bytes ... \n", fname, n);
    start = clock();
    if(sais(T, SA, LCP, (int)n) != 0) {
        fprintf(stderr, "%s: Cannot allocate memory.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    finish = clock();
    fprintf(stderr, "induced SA: %.4f sec\n", (double)(finish - start) / (double)CLOCKS_PER_SEC);

    /* naive BWT: */
    start = clock();
    for (i = 0; i < n; ++i) {
        BWT[i] = (SA[i] > 0) ? T[SA[i]-1] : sentinal;
    }
    finish = clock();
    fprintf(stderr, " naive BWT: %.4f sec\n", (double)(finish - start) / (double)CLOCKS_PER_SEC);

    if (validate) {
        /* check SA: */
        if (sufcheck(T, SA, (int)n, 1) != 0) { exit(EXIT_FAILURE); }

        /* check LCP: */
        for (i = 1; i < n; ++i) {
            l = 0;
            while (T[SA[i]+l]==T[SA[i-1]+l]) ++l;
            if (l != LCP[i]) {
                printf("Error at position %i\n", i);
                printf("%i vs. %i\n", l, LCP[i]);
                for (j = 0; j < 10; j++) printf("%c", T[SA[i]+j]); printf("\n");
                for (j = 0; j < 10; j++) printf("%c", T[SA[i-1]+j]); printf("\n");
                exit(-1);
            }
        }
    }
    else {
        printf("skipping validation\n");
    }

    if (print) {
        printf("Index\tSA\tLCP\tBWT\tSeqID\tSuffix\n");
        for(i=0; i<n; ++i) {
            int len = END[SID[SA[i]]] - SA[i] + 1;
            printf("%d\t%d\t%d\t%c\t%d\t%.*s\n", i, SA[i], LCP[i], BWT[i], SID[SA[i]], len, T+SA[i]);
        }
        /* this is the last time END is used */
        free(END);
    }

    /* we don't even need the input file any longer */
    free(T);

    stack<quad> the_stack;
    start = clock();
    count = 0;
    count_generated = 0;
    quad last_interval;
    the_stack.push(quad());
    LCP[n] = 0;
    for (i = 1; i <= n; ++i) {
        lb = i - 1;
        while (LCP[i] < the_stack.top().lcp) {
            the_stack.top().rb = i - 1;
            last_interval = the_stack.top();
            the_stack.pop();
            process(last_interval, SA, BWT, SID, sentinal, cutoff);
            lb = last_interval.lb;
            if (LCP[i] <= the_stack.top().lcp) {
                last_interval.children.clear();
                the_stack.top().children.push_back(last_interval);
                last_interval = quad();
            }
        }
        if (LCP[i] > the_stack.top().lcp) {
            if (!last_interval.empty()) {
                last_interval.children.clear();
                the_stack.push(quad(LCP[i],lb,INT_MAX,vector<quad>(1, last_interval)));
                last_interval = quad();
            }
            else {
                the_stack.push(quad(LCP[i],lb,INT_MAX));
            }
        }
    }
    the_stack.top().rb = n - 1;
    process(the_stack.top(), SA, BWT, SID, sentinal, cutoff);
    finish = clock();
    fprintf(stderr, "SA has %d internal nodes\n", count);
    fprintf(stderr, "processing: %.4f sec\n", (double)(finish - start) / (double)CLOCKS_PER_SEC);

#if 0
    for (set<pair<int,int> >::iterator pit=pairs.begin();
            pit!=pairs.end(); ++pit) {
        cout << "pair " << pit->first << "," << pit->second << endl;
    }
#endif
    cout << "   unique pairs = " << pairs.size() << endl;
    cout << "generated pairs = " << count_generated << endl;

    /* Deallocate memory. */
    free(SA);
    free(LCP);
    free(BWT);
    free(SID);

    return 0;
}


inline static void process(const quad &q)
{
    ++count;
    cout << q << endl;
}

inline static void pair_check(
        const int &i,
        const int &j,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const char &sentinal)
{
    const int sidi = SID[SA[i]];
    const int sidj = SID[SA[j]];
    if ((BWT[i] != BWT[j] || BWT[i] == sentinal) && sidi != sidj) {
        ++count_generated;
        if (sidi < sidj) {
            pairs.insert(make_pair(sidi,sidj));
        }
        else {
            pairs.insert(make_pair(sidj,sidi));
        }
    }
}

/* this naive pair generation actually works but generates an extreme number of
 * duplicates */
#if 0
inline static void process(
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const char &sentinal,
        const int &cutoff)
{
    ++count;
    //cout << q << endl;

    if (q.lcp < cutoff) return;

    /* naive impl ??? */
    for (int i=q.lb; i<=q.rb; ++i) {
        for (int j=i+1; j<=q.rb; ++j) {
            pair_check(i, j, SA, BWT, SID, sentinal);
        }
    }
}

#else

/* try to reduce number of duplicate pairs generated */
/* we observe that l-intervals (i.e. internal nodes) always have at least two
 * children, but these children could be singleton l-intervals,
 * e.g., [i..j]=[1..1], in addition to l-intervals with non-singleton
 * ranges/quads. For each l-interval, we take the cross product of its child
 * l-intervals. Naively, as above, we could take the cross product of the
 * entire lb/rb range of the l-interval, but this generates too many duplicate
 * pairs. Instead, the complexity should be bounded by the number of exact
 * matches...
 */
inline static void process(
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const char &sentinal,
        const int &cutoff)
{
    const int n_children = q.children.size();
    int child_index = 0;

    ++count;

    if (q.lcp < cutoff) return;

    for (int i=q.lb; i<=q.rb; ++i) {
        int j = i+1;
        if (child_index < n_children) {
            if (i >= q.children[child_index].lb) {
                j = q.children[child_index].rb+1;
                if (i >= q.children[child_index].rb) {
                    ++child_index;
                }
            }
        }
        for (/*nope*/; j<=q.rb; ++j) {
            pair_check(i, j, SA, BWT, SID, sentinal);
        }
    }
}
#endif

