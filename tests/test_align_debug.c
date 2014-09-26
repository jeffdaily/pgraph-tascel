#include "config.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align/align_debug.h"
#include "blosum/blosum62.h"

static float pct(unsigned long long orig_, unsigned long long new_)
{
    float orig = (float)orig_;
    float new = (float)new_;
    return 100.0*(orig - new)/orig;
}


static void print_array(
        const char * filename,
        const int * const restrict array,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len)
{
    int i;
    int j;
    FILE *f = fopen(filename, "w");
    fprintf(f, " ");
    for (j=0; j<s2Len; ++j) {
        fprintf(f, "%3c", s2[j]);
    }
    fprintf(f, "\n");
    for (i=0; i<s1Len; ++i) {
        fprintf(f, "%c", s1[i]);
        for (j=0; j<s2Len; ++j) {
            fprintf(f, "%3d", array[i*s2Len + j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}



int main(int argc, char **argv)
{
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
    const int lena = strlen(seqA);
    const int lenb = strlen(seqB);
    const int longest = MAX(lena,lenb) + 16 /* +16 for woz padding */;
    const int tbl_size = (16+lena)*(16+lenb);
    int score;
    int matches;
    int length;
    int * tbl_pr = malloc(sizeof(int) * longest);
    int * del_pr = malloc(sizeof(int) * longest);
    int * mch_pr = malloc(sizeof(int) * longest);
    int * len_pr = malloc(sizeof(int) * longest);
    int * scr_tbl = malloc(sizeof(int) * tbl_size);
    int * mch_tbl = malloc(sizeof(int) * tbl_size);
    int * len_tbl = malloc(sizeof(int) * tbl_size);

    score = nw_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr, scr_tbl);
    print_array("nw_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = nw_woz_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr, scr_tbl);
    print_array("nw_woz_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = nw_striped_debug(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("nw_stp_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = nw_stats_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, &matches, &length, tbl_pr, del_pr, mch_pr, len_pr, scr_tbl, mch_tbl, len_tbl);
    print_array("nw_stats_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_ref_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_ref_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = nw_stats_woz_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, &matches, &length, tbl_pr, del_pr, mch_pr, len_pr, scr_tbl, mch_tbl, len_tbl);
    print_array("nw_stats_woz_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_woz_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_woz_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = nw_stats_striped_debug(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("nw_stats_stp_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_stp_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_stp_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = sg_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr, scr_tbl);
    print_array("sg_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = sg_woz_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr, scr_tbl);
    print_array("sg_woz_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = sg_striped_debug(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("sg_stp_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = sg_stats_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, &matches, &length, tbl_pr, del_pr, mch_pr, len_pr, scr_tbl, mch_tbl, len_tbl);
    print_array("sg_stats_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_ref_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_ref_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = sg_stats_woz_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, &matches, &length, tbl_pr, del_pr, mch_pr, len_pr, scr_tbl, mch_tbl, len_tbl);
    print_array("sg_stats_woz_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_woz_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_woz_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = sg_stats_striped_debug(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("sg_stats_stp_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_stp_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_stp_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = sw_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr, scr_tbl);
    print_array("sw_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = sw_woz_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, tbl_pr, del_pr, scr_tbl);
    print_array("sw_woz_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = sw_striped_debug(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("sw_stp_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = sw_stats_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, &matches, &length, tbl_pr, del_pr, mch_pr, len_pr, scr_tbl, mch_tbl, len_tbl);
    print_array("sw_stats_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_ref_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_ref_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = sw_stats_woz_debug(seqA, lena, seqB, lenb, 10, 1, blosum62, &matches, &length, tbl_pr, del_pr, mch_pr, len_pr, scr_tbl, mch_tbl, len_tbl);
    print_array("sw_stats_woz_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_woz_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_woz_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = sw_stats_striped_debug(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("sw_stats_stp_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_stp_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_stp_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    free(tbl_pr);
    free(del_pr);
    free(mch_pr);
    free(len_pr);
    free(scr_tbl);
    free(mch_tbl);
    free(len_tbl);

    return 0;
}
