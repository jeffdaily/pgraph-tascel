#include "config.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "alignment.hpp"
#include "blosum/blosum62.h"
#include "timer.h"


int main(int argc, char **argv)
{
    //const char * s1 = "SLPSMRADSFTKELMEKISSVRTSTLTFAPEAGTPRLRDIINKNITEEEILRACRVAYEAGKNQIKLYFMDGLPGETYEDIAGIAALASHVVDEYYRTPGRNKARQPQVTLSVACFIPKPHTPFQWERQNAPEELADKQAFLSGKITDRKVRHNYHDAKVSRIEAVFARGDRRLGRALEEAARRHVRFDAWEDCFDYDGWMDIFETVGIDPAFYANRTIPDDEILPWDMISCGVTKSFLLSERHKAQQAIATPACRDQCSGCGVNRLVDKRYCRWCPGHPESSDSAGRITSDREIRKKPEETSAQKGNVKPARQIRIRFRKYGAMLYISHLDLAKTVMRSIVRSGLPVYYSEGFNPKPKLVFGTPLSVGCGGEAEVLDIRLMKAVSNAEITEKLKAVMPNGVEVTQVYEQKGKLTDVKWAENVIEWRNTDVSPELAEKTEALFQSPVVMMKKSKSGEKEVDITSYIRSLRAEALDGGLRITAVTAAEQENYLNPEYIVQAAERAFGISGENGWHVITRTRLLLADGETDFA";
    //const char * s1 = "SLPSMRADSFTKELMEKISS";
    const char * s1 = "SLPSMRADSFTKELMEKISSSLPSMRADSFTKELMEKISSSLPSMRADSFTKELMEKISSSLPSMRADSFTKELMEKISSSLPSMRADSFTKELMEKISSSLPSMRADSFTKELMEKISSSLPSMRADSFTKELMEKISSSLPSMRADSFTKELMEKISS";
    size_t s1_len = strlen(s1);
    //const char * s2 = "MTNKICIYAISKNEEKFVEKWYDSMKEADAVVVLDTGSTDNTVEKLRKLGATVEVKKIDPWRFDVARNESLKLVPDDCNILMSTDLDEWLEPGWSKPLREKWIEGVHERGVYKYSWSHLKDGSSGRIFRYDKIHSRKWKWMAPVHELLCDEAGSNEYYYDQILDLFDDIHLHHYPDPNKSRGSYLPLLELRAKENPEDWYGLIYLAHEYFYRGKNEKAIALLKRILSEYKDHYSILEKASCYLFMGDGYKAIGDMCEDEEERNKNYGLAKLAYLNAIRTEPSYIEPYLDLSKVYFEEKDFDVAETYIKRGLQNSYRHFTWLERDTS";
    //const char * s2 = "MTNKICIYAISKNEEKFV";
    const char * s2 = "MTNKICIYAISKNEEKFVMTNKICIYAISKNEEKFVMTNKICIYAISKNEEKFVMTNKICIYAISKNEEKFVMTNKICIYAISKNEEKFVMTNKICIYAISKNEEKFVMTNKICIYAISKNEEKFVMTNKICIYAISKNEEKFV";
    size_t s2_len = strlen(s2);
    size_t longest = ::std::max(s1_len,s2_len) + 16;

    int open = -10;
    int gap = -1;
#if USE_SIMILARITIES
    ::pgraph::cell_t result_cell = {0,0,0,0};
#else
    ::pgraph::cell_t result_cell = {0,0,0};
#endif
    ::pgraph::cell_t **tbl = ::pgraph::allocate_cell_table(2, longest);
    ::pgraph::tbl_t **all = ::pgraph::allocate_tbl_table(2, longest);
    int **scr = ::pgraph::allocate_int_table(2, longest);
    int **mat = ::pgraph::allocate_int_table(2, longest);
#if USE_SIMILARITIES
    int **sim = ::pgraph::allocate_int_table(2, longest);
#endif
    int **len = ::pgraph::allocate_int_table(2, longest);
    int **del = ::pgraph::allocate_int_table(2, longest);
    int **ins = ::pgraph::allocate_int_table(2, longest);
    int *woz_scr = new int[longest];
    int *woz_del = new int[longest];
    int *woz_mat = new int[longest];
    int *woz_len = new int[longest];
    unsigned long long timer;
    size_t i = 0;
    size_t limit = 1000;
    //size_t limit = 1;

    timer_init();
    ::std::cout << timer_name() << " timer" << ::std::endl;

    ::std::cout << ::std::endl;
    ::std::cout << "AOS tests" << ::std::endl;
    ::std::cout << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_global_affine(
                s1, s1_len, s2, s2_len, open, gap, tbl, del, ins);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_global_affine"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_semi_affine(
                s1, s1_len, s2, s2_len, open, gap, tbl, del, ins);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_semi_affine"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_local_affine(
                s1, s1_len, s2, s2_len, open, gap, tbl, del, ins);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_local_affine"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    ::std::cout << ::std::endl;
    ::std::cout << "SOA tests" << ::std::endl;
    ::std::cout << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_global_affine(
#if USE_SIMILARITIES
                s1, s1_len, s2, s2_len, open, gap, scr, mat, sim, len, del, ins);
#else
                s1, s1_len, s2, s2_len, open, gap, scr, mat, len, del, ins);
#endif
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_global_affine"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_semi_affine(
#if USE_SIMILARITIES
                s1, s1_len, s2, s2_len, open, gap, scr, mat, sim, len, del, ins);
#else
                s1, s1_len, s2, s2_len, open, gap, scr, mat, len, del, ins);
#endif
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_semi_affine"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_local_affine(
#if USE_SIMILARITIES
                s1, s1_len, s2, s2_len, open, gap, scr, mat, sim, len, del, ins);
#else
                s1, s1_len, s2, s2_len, open, gap, scr, mat, len, del, ins);
#endif
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_local_affine"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    ::std::cout << ::std::endl;
    ::std::cout << "AOS full tests" << ::std::endl;
    ::std::cout << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_global_affine(
                s1, s1_len, s2, s2_len, open, gap, tbl);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_global_affine"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_semi_affine(
                s1, s1_len, s2, s2_len, open, gap, tbl);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_semi_affine"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_local_affine(
                s1, s1_len, s2, s2_len, open, gap, tbl);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_local_affine"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    ::std::cout << ::std::endl;
    ::std::cout << "Wozniak SSE tests" << ::std::endl;
    ::std::cout << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_global_affine_sse(
                s1, s1_len, s2, s2_len, open, gap,
                woz_scr, woz_del, woz_mat, woz_len);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_global_affine_sse"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_semi_affine_sse(
                s1, s1_len, s2, s2_len, open, gap,
                woz_scr, woz_del, woz_mat, woz_len);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_semi_affine_sse"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_local_affine_sse(
                s1, s1_len, s2, s2_len, open, gap,
                woz_scr, woz_del, woz_mat, woz_len);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_global_local_sse"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;

#if 1
    ::std::cout << ::std::endl;
    ::std::cout << "SSW tests" << ::std::endl;
    ::std::cout << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_local_affine_ssw(s1, s1_len, s2, s2_len, open, gap);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_local_affine_ssw"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;
#endif

#if 0
    ::std::cout << ::std::endl;
    ::std::cout << "FASTA tests" << ::std::endl;
    ::std::cout << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_global_affine_fasta(s1, s1_len, s2, s2_len, open, gap);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_local_affine_ssw"
        << "\t" << timer/limit
        << "\t" << result_cell.score
        << "\t" << result_cell.matches
#if USE_SIMILARITIES
        << "\t" << result_cell.similarities
#endif
        << "\t" << result_cell.length
        << ::std::endl;
#endif

    ::pgraph::free_cell_table(tbl, 2);
    ::pgraph::free_int_table(scr, 2);
    ::pgraph::free_int_table(mat, 2);
#if USE_SIMILARITIES
    ::pgraph::free_int_table(sim, 2);
#endif
    ::pgraph::free_int_table(len, 2);
    ::pgraph::free_int_table(del, 2);
    ::pgraph::free_int_table(ins, 2);
    delete [] woz_scr;
    delete [] woz_del;
    delete [] woz_mat;
    delete [] woz_len;

    return 0;
}
