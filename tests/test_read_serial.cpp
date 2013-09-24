/**
 * @author jeff.daily@pnnl.gov
 */
#include "config.h"

#include <stdio.h>
#include <sys/stat.h>

#include <cassert>
#include <climits>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#define ARG_LEN_MAX 1024

int rank = 0;
int nprocs = 0;


int main(int argc, char **argv)
{
    long file_size = -1;
    char *file_buffer = NULL;

    /* sanity check that we got the correct number of arguments */
    if (argc <= 1 || argc >= 3) {
        if (0 == rank) {
            if (argc <= 1) {
                printf("missing input file\n");
            }
            else if (argc >= 3) {
                printf("too many arguments\n");
            }
            printf("usage: test_read sequence_file\n");
        }
        return 1;
    }

    /* we don't want all procs to stat the file! */
    if (0 == rank) {
        struct stat st;
        assert(0 == stat(argv[1], &st));
        file_size = st.st_size;
    }

    if (0 == rank) {
        printf("file_size=%ld\n", file_size);
    }

    /* allocate a buffer for the file, of the entire size */
    file_buffer = new char[file_size];

    if (0 == rank) {
        FILE *file = NULL;
        size_t read_count = 0;
        unsigned long my_seg_count = 0;

        file = fopen(argv[1], "r");
        if (NULL == file) {
            perror("fopen");
            printf("unable to open file on process 0\n");
            assert(0);
        }
        printf("process 0 opened file\n");

        (void)memset(file_buffer, 0, file_size);
        printf("process 0 memset buffer\n");
        read_count = fread(file_buffer, file_size, 1, file);
        if (0 == read_count) {
            printf("unable to read file on process 0\n");
            assert(0);
        }
        printf("process 0 read file\n");
        
        assert(file_buffer[0] == '>');
        ++my_seg_count;
        for (int i=1; i<file_size; ++i) {
            if (file_buffer[i] == '>' && file_buffer[i-1] == '\n') {
                ++my_seg_count;
            }
        }

        printf("[0] my_seg_count=%lu\n", my_seg_count);
    }


    /* clean up */
    delete [] file_buffer;

    return 0;
}
