/**
 * @file test_align.c
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#include "config.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "align.h"
#include "blosum/blosum62.h"
#include "timer.h"

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

int test(const char *seq1, const char *seq2)
{
    size_t i = 0;
    pgraph_cell_t **tbl = NULL;
    int **del = NULL;
    int **ins = NULL;
    size_t seq1_len = strlen(seq1);
    size_t seq2_len = strlen(seq2);
    int max_seq_len = seq1_len > seq2_len ? seq1_len : seq2_len;
    pgraph_cell_t result = {0,0,0};
    int is_edge_answer = 0;
    int sscore = 0;
    size_t maxLen = 0;
    const int AOL = 8;
    const int SIM = 4;
    const int OS = 3;
    const int OPEN = -10;
    const int GAP = -1;
#define NROW 2

    tbl = pgraph_malloc_cell_table(NROW, max_seq_len);
    del = pgraph_malloc_int_table(NROW, max_seq_len);
    ins = pgraph_malloc_int_table(NROW, max_seq_len);

    pgraph_affine_gap_align(
            seq1, seq1_len,
            seq2, seq2_len,
            &result,
            tbl, del, ins,
            OPEN, GAP,
            _blosum62, PGRAPH_MAP_BLOSUM, 'A');

    is_edge_answer = pgraph_is_edge(
            result, seq1, seq1_len, seq2, seq2_len,
            AOL, SIM, OS, &sscore, &maxLen,
            _blosum62, PGRAPH_MAP_BLOSUM, 'A');

    printf("---------------------------------------------------\n");
    printf("result->score=%d\n", result.score);
    printf("result->ndig=%d\n", result.matches);
    printf("result->alen=%d\n", result.length);
    printf("is_edge_answer=%d\n", is_edge_answer);
    printf("alen/maxLen=%f\n", 1.0*result.length/maxLen);
    printf("nmatch/alen=%f\n", 1.0*result.matches/result.length);
    printf("score/sscore=%f\n", 1.0*result.score/sscore);
    printf("timing 10000 calls to pgraph_affine_gap_align\n");
    unsigned long long t = timer_start();
    for (i = 0; i < 10000U; ++i) {
        pgraph_affine_gap_align(
                seq1, seq1_len,
                seq2, seq2_len,
                &result,
                tbl, del, ins,
                OPEN, GAP,
                _blosum62, PGRAPH_MAP_BLOSUM, 'A');
    }
    t = timer_end(t);
    printf("%s timer took %llu units\n", timer_name(), t);

    pgraph_free_cell_table(tbl, NROW);
    pgraph_free_int_table(del, NROW);
    pgraph_free_int_table(ins, NROW);

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
