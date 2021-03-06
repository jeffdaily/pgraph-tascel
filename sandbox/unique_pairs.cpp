#include "config.h"

#include <dirent.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cassert>
#include <cerrno>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;


/* check if a string ends with another string */
bool endswith(string const &fullString, string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}


/* get *.bin files from given dir */
int listdir(string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        if (endswith(dirp->d_name, ".bin")) {
            files.push_back(string(dirp->d_name));
        }
    }
    closedir(dp);
    return 0;
}


size_t get_filesize(const string &file)
{
    struct stat buf;
    if (0 == stat(file.c_str(), &buf)) {
        return buf.st_size;
    }
    else {
        perror("stat");
        return 0;
    }
}


char* read_file(const string &filename)
{
    char *buffer = NULL;
    ifstream file(filename.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()) {
        size_t size = file.tellg();
        buffer = new char[size];
        file.seekg(0, ios::beg);
        file.read(buffer, size);
        file.close();
    }
    return buffer;
}


int main(int argc, char **argv)
{
    if (argc != 3) {
        cout << "usage: <progname> input_dir output_dir" << endl;
        return 1;
    }

    string dir = string(argv[1]);
    vector<string> files = vector<string>();
    listdir(dir,files);

    size_t file_count = files.size();
    size_t width = static_cast<size_t>(log10(file_count)+0.5);
    size_t pc_total = 0;
    size_t new_total = 0;

    cout << fixed << setprecision(3);

    for (size_t i=0; i<file_count; ++i) {
        set<pair<int,int> > pairs;
        size_t pc = 0;
        string infile = dir+'/'+files[i];
        size_t size = get_filesize(infile);
        double count_f = 1.0*size/9.0;
        assert(floor(count_f) == count_f);
        int count = static_cast<int>(count_f);
        cout << "processing " << setw(width) << i+1 << "/" << file_count
            << " " << infile << endl;
        cout << "\t0%" << flush;
        char *buffer = read_file(infile);
        assert(buffer);
        for (size_t j=0; j<size; j+=9) {
            int id1 = *((int*)(&buffer[j]));
            int id2 = *((int*)(&buffer[j+4]));
            char newline = buffer[j+8];
            assert(id1 >= 0);
            assert(id2 >= 0);
            assert(newline == '\n');
            ++pc;
            if (pc%100000 == 0) {
                cout << "\r" << flush << "\t" << (100.0*pc/count) << "%" << flush;
            }
            if (id1 > id2) {
                int tmp = id1;
                id1 = id2;
                id2 = tmp;
            }
            pair<int,int> pair(id1,id2);
            pairs.insert(pair);
        }
        delete [] buffer;
        cout << "\r" << flush;
        pc_total += pc;

        size_t pairs_size = pairs.size();
        string outfile = string(argv[2])+'/'+files[i];
        if (pc) {
            cout << "writing " << outfile
                << " (removed " << (pc-pairs_size)
                << " or " << (1.0*(pc-pairs_size)/pc) << "%)" << endl;
        }
        else {
            cout << "writing " << outfile << " (empty)" << endl;
        }
        cout << "\t0%" << flush;
        new_total += pairs_size;
        pc = 0;
        ofstream out(outfile.c_str(), ios::binary|ios::trunc);
        if (!out) {
            cout << "error opening output file" << endl;
            return 1;
        }
        for (set<pair<int,int> >::iterator it=pairs.begin();
                it!=pairs.end(); ++it) {
            pc += 1;
            if (pc%10000 == 0) {
                cout << "\r" << flush << "\t" << (100.0*pc/pairs_size) << "%" << flush;
            }
            const char *newline = "\n";
            out.write(reinterpret_cast<const char*>(&(it->first)), sizeof(int));
            out.write(reinterpret_cast<const char*>(&(it->second)), sizeof(int));
            out.write(newline, 1);
        }
        cout << "\r" << flush;
        out.close();
        pairs.clear();
    }
    cout << "NEW TOTAL PAIRS = " << new_total << endl;
    cout << "OLD TOTAL PAIRS = " << pc_total << endl;

    return 0;
}
