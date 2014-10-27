#include "config.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align/align.h"
#include "align/align_wozniak_128_16.h"
#include "align/align_wozniak_128_8.h"
#include "align/align_striped_128_16.h"
#include "align/align_striped_128_8.h"
#include "align/align_scan_128_16.h"
#include "align/align_scan_128_8.h"
#include "blosum/blosum62.h"
#include "timer.h"

#define USE_PERCENT_IMPROVED 0

static float pct(unsigned long long orig_, unsigned long long new_)
{
    float orig = (float)orig_;
    float new = (float)new_;
#if USE_PERCENT_IMPROVED
    return 100.0*(orig - new)/orig;
#else
    return orig / new;
#endif
}


int main(int argc, char **argv)
{
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
    //const char *seqA = "MEFYDVAVTV";
    //const char *seqB = "AALGVAARAGFLAAGFASSS";
    const int lena = strlen(seqA);
    const int lenb = strlen(seqB);
    const int longest = MAX(lena,lenb) + 32 /* +32 for woz padding */;
    int score;
    int matches;
    int length;
    unsigned long long timer;
    unsigned long long timer_ref;
    size_t limit = 1000;
    size_t i;
    int * tbl_pr = malloc(sizeof(int) * longest);
    int * del_pr = malloc(sizeof(int) * longest);
    int * mch_pr = malloc(sizeof(int) * longest);
    int * len_pr = malloc(sizeof(int) * longest);

    timer_init();
    printf("%s timer\n", timer_name());
#if USE_PERCENT_IMPROVED
    printf("alg\t\t\t\ttime\t%%imp\tscore\tmatches\tlength\n");
#else
    printf("alg\t\t\t\ttime\tx_imp\tscore\tmatches\tlength\n");
#endif

    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("nw reference\t\t\t%llu\t\t%d\n", timer_ref/limit, score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_scan_row(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("nw scan row\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_scan(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("nw scan\t\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_scan_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("nw scan 128 16\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_scan_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("nw scan 128 8\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_wozniak_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("nw wozniak 128 16\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_wozniak_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("nw wozniak 128 8\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_striped_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("nw striped 128 16\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_striped_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("nw striped 128 8\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("nw stats reference\t\t%llu\t\t%d\t%d\t%d\n", timer_ref/limit, score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_scan(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("nw stats scan reference\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_scan_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length);
    }
    timer = timer_end(timer);
    printf("nw stats scan 128 16\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_wozniak_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("nw stats wozniak 128 16\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_wozniak_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("nw stats wozniak 128 8\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_striped_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length);
    }
    timer = timer_end(timer);
    printf("nw stats striped 128 16\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_striped_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length);
    }
    timer = timer_end(timer);
    printf("nw stats striped 128 8\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("sg reference\t\t\t%llu\t\t%d\n", timer_ref/limit, score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_scan(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("sg scan\t\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_scan_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("sg scan 128 16\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_scan_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("sg scan 128 8\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_wozniak_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("sg wozniak 128 16\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_wozniak_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("sg wozniak 128 8\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_striped_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("sg striped 128 16\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_striped_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("sg striped 128 8\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("sg stats reference\t\t%llu\t\t%d\t%d\t%d\n", timer_ref/limit, score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_wozniak_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("sg stats wozniak 128 16\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_wozniak_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("sg stats wozniak 128 8\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_striped_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sg stats striped 128 16\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("sw reference\t\t\t%llu\t\t%d\n", timer_ref/limit, score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_scan(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("sw scan\t\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_scan_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("sw scan 128 16\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_scan_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("sw scan 128 8\t\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_wozniak_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr);
    }
    timer = timer_end(timer);
    printf("sw wozniak 128 16\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_striped_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__);
    }
    timer = timer_end(timer);
    printf("sw striped 128 16\t\t%llu\t%4.1f\t%d\n", timer/limit, pct(timer_ref,timer), score);

    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer_ref = timer_end(timer_ref);
    printf("sw stats reference\t\t%llu\t\t%d\t%d\t%d\n", timer_ref/limit, score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats_wozniak_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("sw stats wozniak 128 16\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats_striped_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sw stats striped 128 16\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

    free(tbl_pr);
    free(del_pr);
    free(mch_pr);
    free(len_pr);

    return 0;
}
