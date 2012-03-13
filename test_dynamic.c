#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "dynamic.h"

static const char seq1[] =
"MGDSSKKVKDSFDTISEPDSFDEPKGVPISMEPVFSTAAGIRIDVKQESIDKSKKMLNSDLKSKSSSKGG"
"FSSPLVRKNNGSSAFVSPFRREGTSSTTTKRPASGGFEDFEAPPAKKSTSSSSKKSKKHSKKEKKKEFKE"
"IHADVLRVSRIYEKDKFRIILQESSSTPLILATCSYNRGSDIKFGDRIHVDAEVCKKSSSGDVTEIYIDR"
"VLKNKENGAKSGIRRHSIAKKPFCIKPRFIHELSDTKIKKTVVQVNLLDLNLDFYAGCSKCKHSLPEAAN"
"QCEFCKDSQGKSELSMYSRVRVMDFSGQMFINVTTKNMKKLLDLLGYEGFDNWFRFKDPQERQNYVFRPV"
"MVEIEKSNDEWECTDVAEVDWKDFGSYLKHKEDKKKRRSKKKHP";

static const char seq2[] =
"MADVALRITETVARLQKELKCGICCSTYKDPILSTCFHIFCRSCINACFERKRKVQCPICRSVLDKRSCR"
"DTYQITMAVQNYLKLSEAFKKDIENMNTFKSLPPEKMFMESQMPLDITIIPENDGKRCAPDFAIPFLPVR"
"RKRPSRPQPPSAFAEEPAEPVEPPEPATKQPVELQSRVFPLEKLKKDVETSTETYKISREELKNVDIEEY"
"INTLRENSTEIDEIDALFQLMPTMRQFLRNNINQLMEKFHVAPPKKSEKPANRRVSFASSQDLENIKIMT"
"ASESLETPPEPIQKLAQKPEVFKSTQNLIDLNLNTAVKKPVVVASDDDEVVEDSEGELQIDEDDLANVTC"
"ATSVRNIGKSLCAEYIREGRSISQKSTAYLYAIARKCVIVGRQWLVDCITTGLLLSEADYTITSCSSTIP"
"VKIPPSIGSEMGWLRSRNDEHGKLFAGRRFMILRKFTMNPYFDYKQLIELVQQCGGEILSCYENLSPEKL"
"YIIFSKHSKAIEESKNIENLYKCDVVTMEWVLDSISEYLILPTQPYKAVDSIGCLQD";

int main(int argc, char **argv)
{
    cell_t **tbl = NULL;
    int **del = NULL;
    int **ins = NULL;
    size_t seq1_len = strlen(seq1);
    size_t seq2_len = strlen(seq2);
    int max_seq_len = seq1_len > seq2_len ? seq1_len : seq2_len;
    cell_t result;
    is_edge_param_t param;
    int is_edge_answer = 0;

    assert(NROW == 2);
    tbl = alloc_tbl(NROW, max_seq_len);
    del = alloc_int(NROW, max_seq_len);
    ins = alloc_int(NROW, max_seq_len);

    affine_gap_align(seq1, seq1_len, seq2, seq2_len, &result, tbl, del, ins);
    printf("result->score=%d\n", result.score);
    printf("result->ndig=%d\n", result.ndig);
    printf("result->alen=%d\n", result.alen);

    param.AOL = 8;
    param.SIM = 4;
    param.OS = 3;
    is_edge_answer = is_edge(result, seq1, seq1_len, seq2, seq2_len, param);
    printf("is_edge_answer=%d\n", is_edge_answer);

    free_tbl(tbl, NROW);
    free_int(del, NROW);
    free_int(ins, NROW);

    return 0;
}
