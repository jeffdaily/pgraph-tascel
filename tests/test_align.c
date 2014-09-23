#include "config.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align/align.h"
#include "blosum/blosum62.h"
#include "timer.h"


static float pct(unsigned long long orig_, unsigned long long new_)
{
    float orig = (float)orig_;
    float new = (float)new_;
    return 100.0*(orig - new)/orig;
}


int main(int argc, char **argv)
{
    const char *seqA =
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS"
        "SLPSMRADSFTKELMEKISS";
    const char *seqB =
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV"
        "MTNKICIYAISKNEEKFV";
    const int lena = strlen(seqA);
    const int lenb = strlen(seqB);
    const int longest = MAX(lena,lenb) + 16 /* +16 for woz padding */;
    int score;
    int matches;
    int length;
    unsigned long long timer;
    unsigned long long timer_ref;
    size_t limit = 10000;
    size_t i;
    int * tbl_pr = malloc(sizeof(int) * longest);
    int * del_pr = malloc(sizeof(int) * longest);
    int * mch_pr = malloc(sizeof(int) * longest);
    int * len_pr = malloc(sizeof(int) * longest);

    timer_init();
    printf("%s timer\n", timer_name());
    printf("alg\t\t\ttime\t%%imp\tscore\tmatches\tlength\n");

#if 1
    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("nw reference\t\t%llu\t\t%d\n", timer_ref/limit, score);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_woz(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("nw wozniak\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_striped(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("nw striped\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);
#endif

#if 1
    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("nw stats reference\t%llu\t\t%d\t%d\t%d\n", timer_ref/limit, score, matches, length);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_woz(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("nw stats wozniak\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_striped(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length);
    }
    timer = timer_end(timer);
    printf("nw stats striped\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if 1
    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("sg reference\t\t%llu\t\t%d\n", timer_ref/limit, score);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_woz(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("sg wozniak\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_striped(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("sg striped\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);
#endif

#if 1
    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("sg stats reference\t%llu\t\t%d\t%d\t%d\n", timer_ref/limit, score, matches, length);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_woz(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("sg stats wozniak\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_striped(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sg stats striped\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if 1
    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("sw reference\t\t%llu\t\t%d\n", timer_ref/limit, score);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_woz(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("sw wozniak\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_striped(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("sw striped\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);
#endif

#if 1
    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("sw stats reference\t%llu\t\t%d\t%d\t%d\n", timer_ref/limit, score, matches, length);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats_woz(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("sw stats wozniak\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if 1
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats_striped(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sw stats striped\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

    free(tbl_pr);
    free(del_pr);
    free(mch_pr);
    free(len_pr);

    return 0;
}
