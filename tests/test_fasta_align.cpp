#include "config.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "alignment.hpp"
#include "blosum/blosum62.h"

/* from fasta36 */
extern "C" {
#include "defs.h"
#include "param.h"
#include "dropgsw2.h"
#include "smith_waterman_sse2.h"
}

/* from Complete-Striped-Smith-Waterman-Library */
extern "C" {
#include "ssw.h"
}
struct _profile{
    __m128i* profile_byte;  // 0: none
    __m128i* profile_word;  // 0: none
    const int8_t* read;
    const int8_t* mat;
    int32_t readLen;
    int32_t n;
    uint8_t bias;
};

int main(int argc, char **argv)
{
    const char * s1 = "SLPSMRADSFTKELMEKISSVRTSTLTFAPEAGTPRLRDIINKNITEEEILRACRVAYEAGKNQIKLYFMDGLPGETYEDIAGIAALASHVVDEYYRTPGRNKARQPQVTLSVACFIPKPHTPFQWERQNAPEELADKQAFLSGKITDRKVRHNYHDAKVSRIEAVFARGDRRLGRALEEAARRHVRFDAWEDCFDYDGWMDIFETVGIDPAFYANRTIPDDEILPWDMISCGVTKSFLLSERHKAQQAIATPACRDQCSGCGVNRLVDKRYCRWCPGHPESSDSAGRITSDREIRKKPEETSAQKGNVKPARQIRIRFRKYGAMLYISHLDLAKTVMRSIVRSGLPVYYSEGFNPKPKLVFGTPLSVGCGGEAEVLDIRLMKAVSNAEITEKLKAVMPNGVEVTQVYEQKGKLTDVKWAENVIEWRNTDVSPELAEKTEALFQSPVVMMKKSKSGEKEVDITSYIRSLRAEALDGGLRITAVTAAEQENYLNPEYIVQAAERAFGISGENGWHVITRTRLLLADGETDFA";
    size_t s1_len = strlen(s1);
    const char * s2 = "MTNKICIYAISKNEEKFVEKWYDSMKEADAVVVLDTGSTDNTVEKLRKLGATVEVKKIDPWRFDVARNESLKLVPDDCNILMSTDLDEWLEPGWSKPLREKWIEGVHERGVYKYSWSHLKDGSSGRIFRYDKIHSRKWKWMAPVHELLCDEAGSNEYYYDQILDLFDDIHLHHYPDPNKSRGSYLPLLELRAKENPEDWYGLIYLAHEYFYRGKNEKAIALLKRILSEYKDHYSILEKASCYLFMGDGYKAIGDMCEDEEERNKNYGLAKLAYLNAIRTEPSYIEPYLDLSKVYFEEKDFDVAETYIKRGLQNSYRHFTWLERDTS";
    size_t s2_len = strlen(s2);

    /* prep needed for ssw_init */
    s_profile *the_s_profile = NULL;
    int8_t *s1_num = new int8_t[s1_len];
    int8_t *s2_num = new int8_t[s2_len];
    s_align *result = NULL;

    /* This table is used to transform amino acid letters into numbers. */
    static const int8_t table[128] = {
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
    };

    /* initialize score matrix */
    for (size_t m = 0; m < s1_len; ++m) s1_num[m] = table[(int)s1[m]];
    for (size_t m = 0; m < s2_len; ++m) s2_num[m] = table[(int)s2[m]];

    /* call into ssw routines */
    the_s_profile = ssw_init(s1_num, s1_len, blosum62__, 24, 2);
    result = ssw_align(the_s_profile, s2_num, s2_len, 10, 1, 2, 0, 0, s1_len/2);
    ::std::cout << result->score1 << ::std::endl;

    ::pgraph::cell_t result_cell = ::pgraph::align_local_affine_ssw(s1, s1_len, s2, s2_len, -10, -1);
    ::std::cout << result_cell.score
        << " " << result_cell.matches
        << " " << result_cell.length
        << ::std::endl;
#if 1
    int n0 = s1_len;
    int e;
    int f;
    int i;
    int j;
    int nsq = 24;
    int l;
    int data;
    int bias;
    unsigned char *  pc;
    unsigned short * ps;
    int overflow;
    int n_count;
    int col_len;

    /* prep needed for fasta */
    struct f_struct *f_str = new struct f_struct;
    f_str->ss = NULL;
    f_str->waa_s = NULL;
    f_str->waa_a = NULL;
    f_str->pam2p[0] = NULL;
    f_str->pam2p[1] = NULL;
    f_str->max_res = 0;
    ::std::fill(f_str->aa0_f, f_str->aa0_f+MAXSQ, 0.0);
    f_str->kar_p = NULL;
    f_str->e_cut = 0.0;
    f_str->show_ident = 0;
    f_str->max_repeat = 0;
    f_str->bias = 0;
    f_str->ceiling = 0;
    f_str->word_score = NULL;
    f_str->byte_score = NULL;
    f_str->workspace = NULL;
    f_str->alphabet_size = 0;
    f_str->word_score_memory = NULL;
    f_str->byte_score_memory = NULL;
    f_str->workspace_memory = NULL;
    f_str->try_8bit = 0;
    f_str->done_8bit = 0;
    f_str->done_16bit = 0;

    f_str->workspace_memory  = (void *)malloc(3*16*(MAXTST+MAXLIB+32)+256);
    f_str->workspace  = (void *) ((((size_t) f_str->workspace_memory) + 255) & (~0xff));
    f_str->word_score_memory = (void *)malloc((n0 + 32) * sizeof (short) * (nsq + 1) + 256);
    f_str->byte_score_memory = (void *)malloc((n0 + 32) * sizeof (char) * (nsq + 1) + 256);
    f_str->word_score = (unsigned short *) ((((size_t) f_str->word_score_memory) + 255) & (~0xff));
    f_str->byte_score = (unsigned char *) ((((size_t) f_str->byte_score_memory) + 255) & (~0xff));

    overflow = 0;

    {
#if 1
        /* Classical simple substitution matrix */
        /* Find the bias to use in the substitution matrix */
        bias = 127;
        for (i = 0; i < nsq ; i++) {
            for (j = 0; j < nsq ; j++) {
                data = blosum62[i][j];
                if (data < -128) {
                    fprintf(stderr,"*** ERROR *** data out of range: %d[%d,%d]\n",
                            data, i, j);
                }
                if (data < bias) {
                    bias = data;
                }
            }
        }
#endif

#if 0
        /* Fill our specially organized byte- and word-size scoring arrays. */
        ps = f_str->word_score;
        col_len = (n0 + 7) / 8;
        n_count = (n0 + 7) & 0xfffffff8;
        for (f = 0; f < n_count; ++f) {
            *ps++ = 0;
        }
        for (f = 1; f < nsq ; f++) {
            for (e = 0; e < col_len; e++) {
                for (i = e; i < n_count; i += col_len) {
                    if (i >= n0) {
                        data = 0;
                    } else {
                        data = ppst->pam2[ip][aa0[i]][f];
                    }
                    *ps++ = (unsigned short)data;
                }
            }
        }
#endif

#if 1
        pc = f_str->byte_score;
        col_len = (n0 + 15) / 16;
        n_count = (n0 + 15) & 0xfffffff0;
        for (f = 0; f < n_count; ++f) {
            *pc++ = 0;
        }
        for (f = 0; f < nsq ; f++) {
            for (e = 0; e < col_len; e++) {
                for (i = e; i < n_count; i += col_len) {
                    if (i >= n0) {
                        data = -bias;
                    } else {
                        data = blosum62[table[s1[i]]][f] - bias;
                    }
                    if (data > 255) {
                        printf("Fatal error. data: %d bias: %d, position: %d/%d, "
                                "Score out of range for 8-bit SSE2 datatype.\n",
                                data, bias, f, e);
                        exit(1);
                    }
                    *pc++ = (unsigned char)data;
                }
            }
        }
#endif
    }

    f_str->bias = (unsigned char) (-bias);
    //f_str->alphabet_size = nsq+1;

    /* Some variable to keep track of how many 8-bit runs we need to rerun
     * in 16-bit accuracy. If there are too many reruns it can be faster
     * to use 16-bit alignments directly.
     */

    /* We can only do 8-bit alignments if the scores were small enough. */
    //f_str->try_8bit = (overflow == 0) ? 1 : 0;
    f_str->try_8bit = 0;

#if 1
    int score = smith_waterman_sse2_byte(
            (const unsigned char *)s1,
            //(unsigned char *)(the_s_profile->profile_byte),
            f_str->byte_score,
            s1_len,
            (const unsigned char *)s2,
            s2_len,
            0,
            10,
            1,
            f_str);
#else
    int score = smith_waterman_sse2_word(
            (const unsigned char *)s1,
            (unsigned short *)(the_s_profile->profile_byte),
            s1_len,
            (const unsigned char *)s2,
            s2_len,
            10,
            1,
            f_str);
#endif

    ::std::cout << score << ::std::endl;
#endif

    delete [] s1_num;
    delete [] s2_num;
    free(f_str->workspace_memory);
    delete f_str;

    return 0;
}
