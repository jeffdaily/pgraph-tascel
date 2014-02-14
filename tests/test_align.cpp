/**
 * @file test_align.cpp
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "alignment.hpp"
#include "alignment_defs.hpp"
#include "blosum/blosum62.h"
#include "timer.h"

using namespace pgraph;

static const char seq0[] = 
"MMVLEIRNVAVIGAGSMGHAIAEVVAIHGFNVKLMDVSEDQLKRAMEKIEEGLRKSYERGYISEDPEKVLKRIEATADLIEVAKDADLVIEAIPEIFDLKKKVFSEIEQYCPDHTIFATNTSSLSITKLAEATKRPEKFIGMHFFNPPKILKLLEIVWGEKTSEETIRIVEDFARKIDRIIIHVRKDVPGFIVNRIFVTMSNEASWAVEMGEGTIEEIDSAVKYRLGLPMGLFELHDVLGGGSVDVSYHVLEYYRQTLGESYRPSPLFERLFKAGHYGKKTGKGFYDWSEGKTNEVPLRAGANFDLLRLVAPAVNEAAWLIEKGVASAEEIDLAVLHGLNYPRGLLRMADDFGIDSIVKKLNELYEKYNGEERYKVNPVLQKMVEEGKLGRTTGEGFYKYGD";

static const char seq19[] =
"MKVKKVCVLGAGAMGSGIAQVCATAGYEVWVRDIKQEFLDRGKAAIEKNLQKAVSKGKMTEEKAKEILSRIHFTLDMEEAVKDADLVIEAVPEIMDLKKQVFAEVQKYAKPECIFASNTSGLSITELGNATDRPEKFLGLHFFNPPPVMALVEVIKGEKTSDETIKFGVEFVKSLGKVPVVVKKDVAGFIVNRILVPYLVLAIDDVEKGVATKEEIDATMMYKYGFPMGPIELSDFVGLDILYHASQQWDIVPQSKLLEEKFKANELGMKTGKGFYDWSAGRPKIPQELAGKYDAIRLIAPMVNIAADLIAMGVADAKDIDTAMKLGTNMPKGPCELGDEIGLDVILAKVEELYKEKGFEILKPSEHLKKMVSEGKLGEKSGEGFYSYGK";

static const char seq47[] = 
"MLLLPSRGKLYVVGIGPGKEELMTLKAKRAIEEADYIVGYQTYVDRISHLIEGKKVVTTPMRKELDRVKIALELAKEHVVALISGGDPSIYGILPLVIEYAVEKKVDVEIEAIPGVTAASAASSLLGSAISGDFAVVSLSDLLVPWSVVEKRLLYALSGDFVVAIYNPSSRRRKENFRKAMEIVRRFRGDAWVGVVRNAGREGQQVEIRRVSEVDEVDMNTILIVGNSETKVVDGKMFTPRGYSNKYNIG";

static const char seq1576[] =
"DQADMVRRVKNYENGFINNPIVISPTTTVGEAKSMKEKYGFAGFPVTADGKRNAKLVGAITSRDIQFVEDNSLLVQDVMTKNPVTGAQGITLSEGNEILKKIKKGRLLVVDEKGNLVSMLSRTDLMKNQKYPLASKSANTKQLLWGASIGTMDADKERLRLLVKAGLDVVILDS";

static const char seq1800[] =
"NPIVISPTTTVGEAKSMKERFGFSGFPVTEDGKRNGKLMGIVTSRDIQFVEDNSLLVQDVMTKNPVTGAQGITLSEGNEILKKIKKGKLLIVDDNGNLVSMLSRTDLMKN";

static const char seq_original1[] =
"MGDSSKKVKDSFDTISEPDSFDEPKGVPISMEPVFSTAAGIRIDVKQESIDKSKKMLNSDLKSKSSSKGG"
"FSSPLVRKNNGSSAFVSPFRREGTSSTTTKRPASGGFEDFEAPPAKKSTSSSSKKSKKHSKKEKKKEFKE"
"IHADVLRVSRIYEKDKFRIILQESSSTPLILATCSYNRGSDIKFGDRIHVDAEVCKKSSSGDVTEIYIDR"
"VLKNKENGAKSGIRRHSIAKKPFCIKPRFIHELSDTKIKKTVVQVNLLDLNLDFYAGCSKCKHSLPEAAN"
"QCEFCKDSQGKSELSMYSRVRVMDFSGQMFINVTTKNMKKLLDLLGYEGFDNWFRFKDPQERQNYVFRPV"
"MVEIEKSNDEWECTDVAEVDWKDFGSYLKHKEDKKKRRSKKKHP";

static const char seq_original2[] =
"MADVALRITETVARLQKELKCGICCSTYKDPILSTCFHIFCRSCINACFERKRKVQCPICRSVLDKRSCR"
"DTYQITMAVQNYLKLSEAFKKDIENMNTFKSLPPEKMFMESQMPLDITIIPENDGKRCAPDFAIPFLPVR"
"RKRPSRPQPPSAFAEEPAEPVEPPEPATKQPVELQSRVFPLEKLKKDVETSTETYKISREELKNVDIEEY"
"INTLRENSTEIDEIDALFQLMPTMRQFLRNNINQLMEKFHVAPPKKSEKPANRRVSFASSQDLENIKIMT"
"ASESLETPPEPIQKLAQKPEVFKSTQNLIDLNLNTAVKKPVVVASDDDEVVEDSEGELQIDEDDLANVTC"
"ATSVRNIGKSLCAEYIREGRSISQKSTAYLYAIARKCVIVGRQWLVDCITTGLLLSEADYTITSCSSTIP"
"VKIPPSIGSEMGWLRSRNDEHGKLFAGRRFMILRKFTMNPYFDYKQLIELVQQCGGEILSCYENLSPEKL"
"YIIFSKHSKAIEESKNIENLYKCDVVTMEWVLDSISEYLILPTQPYKAVDSIGCLQD";

static void print_result(cell_t result, const char *msg) {
    printf("---------------------------------------------------\n");
    printf("%s\n", msg);
    printf("result->score=%d\n", result.score);
    printf("result->ndig=%d\n", result.matches);
    printf("result->alen=%d\n", result.length);
}


static int callback(char ch1, char ch2) {
    return blosum62_[MAP_BLOSUM[ch1-'A']][MAP_BLOSUM[ch2-'A']];
}


int test(const char *seq1, const char *seq2)
{
    size_t i = 0;
    cell_t **tbl = NULL;
    int **del = NULL;
    int **ins = NULL;
    size_t seq1_len = strlen(seq1);
    size_t seq2_len = strlen(seq2);
    int max_seq_len = seq1_len > seq2_len ? seq1_len : seq2_len;
    cell_t result = {0,0,0};
    bool is_edge_answer = false;
    int sscore = 0;
    size_t maxLen = 0;
    const int AOL = 8;
    const int SIM = 4;
    const int OS = 3;
    const int OPEN = -10;
    const int GAP = -1;
    unsigned long long t = timer_start();
#define NROW 2

    tbl = allocate_cell_table(NROW, max_seq_len);
    del = allocate_int_table(NROW, max_seq_len);
    ins = allocate_int_table(NROW, max_seq_len);
    select_blosum(62);

    result = affine_gap_align(
            seq1, seq1_len,
            seq2, seq2_len,
            tbl, del, ins,
            OPEN, GAP,
            blosum62_, MAP_BLOSUM, 'A');
    print_result(result, "custom map");

    result = affine_gap_align(
            seq1, seq1_len,
            seq2, seq2_len,
            tbl, del, ins,
            OPEN, GAP,
            callback);
    print_result(result, "callback");

    result = affine_gap_align_blosum(
            seq1, seq1_len,
            seq2, seq2_len,
            tbl, del, ins,
            OPEN, GAP);
    print_result(result, "blosum");

    result = align_semi_affine(
            seq1, seq1_len,
            seq2, seq2_len,
            blosum62_, MAP_BLOSUM, 'A',
            OPEN, GAP,
            tbl, del, ins);
    print_result(result, "new semi custom map");

    result = align_semi_affine(
            seq1, seq1_len,
            seq2, seq2_len,
            callback,
            OPEN, GAP,
            tbl, del, ins);
    print_result(result, "new semi callback");

    result = align_semi_affine(
            seq1, seq1_len,
            seq2, seq2_len,
            OPEN, GAP,
            tbl, del, ins);
    print_result(result, "new semi blosum");

    result = align_global_affine(
            seq1, seq1_len,
            seq2, seq2_len,
            blosum62_, MAP_BLOSUM, 'A',
            OPEN, GAP,
            tbl, del, ins);
    print_result(result, "new global custom map");

    result = align_global_affine(
            seq1, seq1_len,
            seq2, seq2_len,
            callback,
            OPEN, GAP,
            tbl, del, ins);
    print_result(result, "new global callback");

    result = align_global_affine(
            seq1, seq1_len,
            seq2, seq2_len,
            OPEN, GAP,
            tbl, del, ins);
    print_result(result, "new global blosum");

    result = align_local_affine(
            seq1, seq1_len,
            seq2, seq2_len,
            blosum62_, MAP_BLOSUM, 'A',
            OPEN, GAP,
            tbl, del, ins);
    print_result(result, "new local custom map");

    result = align_local_affine(
            seq1, seq1_len,
            seq2, seq2_len,
            callback,
            OPEN, GAP,
            tbl, del, ins);
    print_result(result, "new local callback");

    result = align_local_affine(
            seq1, seq1_len,
            seq2, seq2_len,
            OPEN, GAP,
            tbl, del, ins);
    print_result(result, "new local blosum");

#if HAVE_EMMINTRIN_H
    result = align_local_affine_ssw(
            seq1, seq1_len,
            seq2, seq2_len,
            OPEN, GAP);
    print_result(result, "new local blosum ssw");
#endif

    printf("---------------------------------------------------\n");

    t = timer_start();
    for (i = 0; i < 10000U; ++i) {
        result = affine_gap_align(
                seq1, seq1_len,
                seq2, seq2_len,
                tbl, del, ins,
                OPEN, GAP,
                blosum62_, MAP_BLOSUM, 'A');
    }
    t = timer_end(t);
    printf("    custom map %s timer took %llu units\n", timer_name(), t);

    t = timer_start();
    for (i = 0; i < 10000U; ++i) {
        result = affine_gap_align(
                seq1, seq1_len,
                seq2, seq2_len,
                tbl, del, ins,
                OPEN, GAP,
                callback);
    }
    t = timer_end(t);
    printf("      callback %s timer took %llu units\n", timer_name(), t);

    t = timer_start();
    for (i = 0; i < 10000U; ++i) {
        result = affine_gap_align_blosum(
                seq1, seq1_len,
                seq2, seq2_len,
                tbl, del, ins,
                OPEN, GAP);
    }
    t = timer_end(t);
    printf("        blosum %s timer took %llu units\n", timer_name(), t);

    t = timer_start();
    for (i = 0; i < 10000U; ++i) {
        result = align_semi_affine(
                seq1, seq1_len,
                seq2, seq2_len,
                blosum62_, MAP_BLOSUM, 'A',
                OPEN, GAP,
                tbl, del, ins);
    }
    t = timer_end(t);
    printf("new custom map %s timer took %llu units\n", timer_name(), t);

    t = timer_start();
    for (i = 0; i < 10000U; ++i) {
        result = align_semi_affine(
                seq1, seq1_len,
                seq2, seq2_len,
                callback,
                OPEN, GAP,
                tbl, del, ins);
    }
    t = timer_end(t);
    printf(" new  callback %s timer took %llu units\n", timer_name(), t);

    t = timer_start();
    for (i = 0; i < 10000U; ++i) {
        result = align_semi_affine(
                seq1, seq1_len,
                seq2, seq2_len,
                OPEN, GAP,
                tbl, del, ins);
    }
    t = timer_end(t);
    printf("   new  blosum %s timer took %llu units\n", timer_name(), t);

    t = timer_start();
    for (i = 0; i < 10000U; ++i) {
        result = align_local_affine_ssw(
                seq1, seq1_len,
                seq2, seq2_len,
                OPEN, GAP);
    }
    t = timer_end(t);
    printf("new ssw blosum %s timer took %llu units\n", timer_name(), t);

    free_cell_table(tbl, NROW);
    free_int_table(del, NROW);
    free_int_table(ins, NROW);

    return 0;
}


int main(int argc, char **argv)
{
    //test(seq_original1, seq_original2);
    //test(seq0, seq19);
    //test(seq0, seq47);
    test(seq1576, seq1800);
    return 0;
}
