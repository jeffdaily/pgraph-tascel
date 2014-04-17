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
    const char * s1 = "SLPSMRADSFTKELMEKISSVRTSTLTFAPEAGTPRLRDIINKNITEEEILRACRVAYEAGKNQIKLYFMDGLPGETYEDIAGIAALASHVVDEYYRTPGRNKARQPQVTLSVACFIPKPHTPFQWERQNAPEELADKQAFLSGKITDRKVRHNYHDAKVSRIEAVFARGDRRLGRALEEAARRHVRFDAWEDCFDYDGWMDIFETVGIDPAFYANRTIPDDEILPWDMISCGVTKSFLLSERHKAQQAIATPACRDQCSGCGVNRLVDKRYCRWCPGHPESSDSAGRITSDREIRKKPEETSAQKGNVKPARQIRIRFRKYGAMLYISHLDLAKTVMRSIVRSGLPVYYSEGFNPKPKLVFGTPLSVGCGGEAEVLDIRLMKAVSNAEITEKLKAVMPNGVEVTQVYEQKGKLTDVKWAENVIEWRNTDVSPELAEKTEALFQSPVVMMKKSKSGEKEVDITSYIRSLRAEALDGGLRITAVTAAEQENYLNPEYIVQAAERAFGISGENGWHVITRTRLLLADGETDFA";
    size_t s1_len = strlen(s1);
    const char * s2 = "MTNKICIYAISKNEEKFVEKWYDSMKEADAVVVLDTGSTDNTVEKLRKLGATVEVKKIDPWRFDVARNESLKLVPDDCNILMSTDLDEWLEPGWSKPLREKWIEGVHERGVYKYSWSHLKDGSSGRIFRYDKIHSRKWKWMAPVHELLCDEAGSNEYYYDQILDLFDDIHLHHYPDPNKSRGSYLPLLELRAKENPEDWYGLIYLAHEYFYRGKNEKAIALLKRILSEYKDHYSILEKASCYLFMGDGYKAIGDMCEDEEERNKNYGLAKLAYLNAIRTEPSYIEPYLDLSKVYFEEKDFDVAETYIKRGLQNSYRHFTWLERDTS";
    size_t s2_len = strlen(s2);

    int open = -10;
    int gap = -1;
    ::pgraph::cell_t result_cell = {0,0,0,0};
    ::pgraph::cell_t **tbl = ::pgraph::allocate_cell_table(2, s1_len);
    int **del = ::pgraph::allocate_int_table(2, s1_len);
    int **ins = ::pgraph::allocate_int_table(2, s1_len);
    unsigned long long timer;
    size_t i = 0;
    size_t limit = 100;

    timer_init();
    ::std::cout << timer_name() << " timer" << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_global_affine(
                s1, s1_len, s2, s2_len, open, gap, tbl, del, ins);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_global_affine" << ::std::endl
        << result_cell.score
        << " " << result_cell.matches
        << " " << result_cell.similarities
        << " " << result_cell.length
        << ::std::endl;
    ::std::cout << timer/limit << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_semi_affine(
                s1, s1_len, s2, s2_len, open, gap, tbl, del, ins);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_semi_affine" << ::std::endl
        << result_cell.score
        << " " << result_cell.matches
        << " " << result_cell.similarities
        << " " << result_cell.length
        << ::std::endl;
    ::std::cout << timer/limit << ::std::endl;

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_local_affine(
                s1, s1_len, s2, s2_len, open, gap, tbl, del, ins);
    }
    timer = timer_end(timer);
    ::std::cout
        << "align_local_affine" << ::std::endl
        << result_cell.score
        << " " << result_cell.matches
        << " " << result_cell.similarities
        << " " << result_cell.length
        << ::std::endl;
    ::std::cout << timer/limit << ::std::endl;

#if 0
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result_cell = ::pgraph::align_global_affine_sse(s1, s1_len, s2, s2_len, open, gap);
    }
    timer = timer_end(timer);
    ::std::cout
        << result_cell.score
        << " " << result_cell.matches
        << " " << result_cell.similarities
        << " " << result_cell.length
        << ::std::endl;
    ::std::cout << timer/limit << ::std::endl;
#endif

    ::pgraph::free_cell_table(tbl, 2);
    ::pgraph::free_int_table(del, 2);
    ::pgraph::free_int_table(ins, 2);

    return 0;
}
