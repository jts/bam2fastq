/*
Copyright 2010, HudsonAlpha Institute for Biotechnology

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

bam2fastq.cpp

Written by Phillip Dexheimer
*/

#include "sam.h"
#include <getopt.h>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cctype>

using namespace std;

const char version[] = "1.1.0";
const char shortopts[] = "o:vhfqs";

int save_aligned = 1;
int save_unaligned = 1;
int save_filtered = 1;
int overwrite_files = 0;
int stdout_pairs = 0;
int stdout_all = 0;
int print_msgs = 1;
int strict = 0;

static struct option longopts[] = {
    { "help",            no_argument,       NULL,           'h' },
    { "version",         no_argument,       NULL,           'v' },
    { "output",          required_argument, NULL,           'o' },
    { "force",           no_argument,       NULL,           'f' },
    { "quiet",           no_argument,       NULL,           'q' },
    { "strict",          no_argument,       NULL,           's' },
    { "overwrite",       no_argument,       &overwrite_files,0  },
    { "aligned",         no_argument,       &save_aligned,   1  },
    { "no-aligned",      no_argument,       &save_aligned,   0  },
    { "unaligned",       no_argument,       &save_unaligned, 1  },
    { "no-unaligned",    no_argument,       &save_unaligned, 0  },
    { "filtered",        no_argument,       &save_filtered,  1  },
    { "no-filtered",     no_argument,       &save_filtered,  0  },
    { "pairs-to-stdout", no_argument,       &stdout_pairs,   1  },
    { "all-to-stdout",   no_argument,       &stdout_all,     1  },
    { NULL,           0,                 NULL,            0  }
};

void usage(int error=1) {
    cerr << "bam2fastq v" << version << " - extract sequences from a BAM file" << endl
         << endl
         << "Usage: bam2fastq [options] <bam file>" << endl
         << endl
         << "Options:" << endl
         << "  -o FILENAME, --output FILENAME" << endl
         << "       Specifies the name of the FASTQ file(s) that will be generated.  May" << endl
         << "       contain the special characters % (replaced with the lane number) and" << endl
         << "       # (replaced with _1 or _2 to distinguish PE reads, _M for unpaired reads)." << endl
         << "       [Default: s_%#_sequence.txt]" << endl << endl
         << "  --pairs-to-stdout" << endl
         << "       Write the paired reads to stdout " << endl << endl
         << "  --all-to-stdout" << endl
         << "       Write all reads to stdout, ignoring pairing" << endl << endl
         << "  -f, --force, --overwrite" << endl
         << "       Create output files specified with --output, overwriting existing" << endl
         << "       files if necessary [Default: exit program rather than overwrite files]" << endl << endl
         << "  --aligned" << endl
         << "  --no-aligned" << endl
         << "       Reads in the BAM that are aligned will (will not) be extracted." << endl
         << "       [Default: extract aligned reads]" << endl << endl
         << "  --unaligned" << endl
         << "  --no-unaligned" << endl
         << "       Reads in the BAM that are not aligned will (will not) be extracted." << endl
         << "       [Default: extract unaligned reads]" << endl << endl
         << "  --filtered" << endl
         << "  --no-filtered" << endl
         << "       Reads that are marked as failing QC checks will (will not) be extracted." << endl
         << "       [Default: extract filtered reads]" << endl << endl
         << "  -q, --quiet" << endl
         << "       Suppress informational messages [Default: print messages]" << endl << endl
         << "  -s, --strict" << endl
         << "       Keep bam2fastq's processing to a minimum, assuming that the BAM strictly"
         << "       meets specifications. [Default: allow some errors in the BAM]" << endl << endl
         << endl;
    exit(error);
}

map<int, char> bases;

const string get_pair_name(const bam1_t *b) {
    return string(bam1_qname(b));
}

//Returns 0 for read 1 and 1 for read 2.
//Counterintuitive, but good for lookups
const int get_read_idx(const bam1_t *b) {
    return !(b->core.flag & BAM_FREAD1);
}

const string get_read_name(const bam1_t *b) {
    static const char suffix[2][3] = { "/1", "/2" };
    string name = get_pair_name(b);
    if (b->core.flag & BAM_FPAIRED)
        name += suffix[get_read_idx(b)];
    return name;
}

int get_lane_id(const bam1_t *b) {
    string name(get_pair_name(b));
    size_t start = name.find(":");
    if (start == string::npos)
        return 0;
    start++;
    size_t stop = name.find(":", start);
    if (stop == string::npos || stop == start)
        return 0;
    istringstream ist(string(name.begin()+start, name.begin()+stop));
    int lane;
    ist >> lane;
    return lane;
}

const string get_sequence(const bam1_t *b) {
    const uint8_t *seq = bam1_seq(b);
    size_t len = b->core.l_qseq;
    string sequence("");
    sequence.reserve(len);
    uint8_t offset = (b->core.flag & BAM_FREVERSE) ? 16 : 0;
    for (size_t i=0; i<len; i++) {
        sequence += bases[bam1_seqi(seq, i) + offset];
    }
    if (offset)
        sequence = string(sequence.rbegin(), sequence.rend());
    return sequence;
}

const string get_qualities(const bam1_t *b) {
    const uint8_t *qual = bam1_qual(b);
    size_t len = b->core.l_qseq;
    string quality("");
    quality.reserve(len);
    for (size_t i=0; i<len; i++) {
        quality += char_traits<char>::to_char_type(char_traits<char>::to_int_type(qual[i])+33);
    }
    if (b->core.flag & BAM_FREVERSE)
        quality = string(quality.rbegin(), quality.rend());
    return quality;
}

void mangle(string &name) {
    if (name.length() < 3)
        return;
    if (isdigit(name[name.length()-1]) && !isdigit(name[name.length()-2]))
        name.erase(name.length()-2);
}

vector<ostream *> initialize_all_stdout() {
    vector<ostream *> files;
    files.push_back(&std::cout);
    return files;
}

vector<ostream *> initialize_paired_stdout() {
    vector<ostream *> files;
    files.push_back(&std::cout);
    files.push_back(&std::cout);
    files.push_back(new ofstream("unpaired_reads.fastq", ios_base::out));
    return files;
}

vector<ostream *> initialize_output(const string &out_template, int lane) {
    
    vector<ostream *> files;
    string output(out_template);
    
    //First, replace % in the filename with the lane number
    size_t laneMarker = output.find('%');
    if (laneMarker != string::npos) {
        if (lane == 0) {
            cerr << "The lane could not be determined from the reads.  Specify output files" << endl
                 << "(using --output) that do not include the lane number (%)" << endl;
            return files;
        }
        ostringstream laneStr;
        laneStr << lane;
        output.replace(laneMarker, 1, laneStr.str());
    }
    
    //Replace # with read number and open ofstreams
    size_t readMarker = output.find('#');

    //Replace
    if (readMarker == string::npos) {
        std::cerr << "Plesae ensure that the output filename (--output) includes" << endl
                  << "a # symbol to be replaced with the read number" << endl;
        return files;
    }

    string file1(output);
    file1.replace(readMarker, 1, "_1");
    string file2(output);
    file2.replace(readMarker, 1, "_2");
    string file3(output);
    file3.replace(readMarker, 1, "_M");

    if (print_msgs) {
        cerr << "This looks like paired data from lane " << lane << "." << endl
             << "Output will be in " << file1 << " and " << file2 << endl
             << "Single-end reads will be in " << file3 << endl;
    }

    //If we're not going to overwrite, check to see if the files exist
    if (!overwrite_files) {
        ifstream test;
        test.open(file1.c_str());
        if (test.is_open()) {
            cerr << "ERROR: " << file1 << " already exists.  Specify --force to overwrite" << endl;
            return files;
        }
        //test is not open, so don't try closing it
        test.open(file2.c_str());
        if (test.is_open()) {
            cerr << "ERROR: " << file2 << " already exists.  Specify --force to overwrite" << endl;
            return files;
        }
    }

    files.push_back(new ofstream(file1.c_str(), ios_base::out));
    files.push_back(new ofstream(file2.c_str(), ios_base::out));
    files.push_back(new ofstream(file3.c_str(), ios_base::out));

    return files;
}

//Effective STL, Item 7
//Except, of course, that I'm not using smart pointers
struct DeleteObject {
    template<typename T>
    void operator()(const T *ptr) const {
        delete ptr;
    }
};

void parse_bamfile(const char *bam_filename, const string &output_template) {
    samfile_t *sam = samopen(bam_filename, "rb", NULL);
    if (sam == NULL) {
        cerr << "Could not open " << bam_filename << endl;
        return;
    }
    bam1_t *read = bam_init1();
    size_t exported = 0;
    size_t all_seen = 0;

    bam_read1(sam->x.bam, read);
    int lane = get_lane_id(read);

    vector<ostream *> output;
    if(stdout_pairs) {
        output = initialize_paired_stdout();
    } else if(stdout_all) {
        output = initialize_all_stdout();
    } else {
        output = initialize_output(output_template, lane);
    }

    if (output.empty())
        return;

    map<string, string> unPaired;
    map<string, string>::iterator position;

    do {
        all_seen++;
        if (!save_aligned && !(read->core.flag & BAM_FUNMAP))
            continue;
        if (!save_unaligned && (read->core.flag & BAM_FUNMAP))
            continue;
        if (!save_filtered && (read->core.flag & BAM_FQCFAIL))
            continue;

        exported++;
        ostringstream ostr;
        ostr << "@" << get_read_name(read) << endl
             << get_sequence(read) << endl
             << "+" << endl
             << get_qualities(read) << endl;

        //Paired-end is complicated, because both members of the pair
        //have to be output at the same position of the two files

        // If there is only one output filehandle we don't care about
        // pairing and can write immediately
        if(output.size() == 1) {
            *output[0] << ostr.str();
        } else if( !(read->core.flag & BAM_FPAIRED) ) {
            // Is this an unpaired read in a BAM with pairs? write to the _M file
            *output[2] << ostr.str();
        } else {

            // Search for the pair in the map
            string pairName(get_pair_name(read));
            if (!strict)
                mangle(pairName);
            position = unPaired.find(pairName);
            if (position == unPaired.end()) {
                //I haven't seen the other member of this pair, so just save it
                unPaired[pairName] = ostr.str();
            } else {
                //Aha!  This will be the second of the two.  Dump them both,
                //then clean up
                int r_idx = get_read_idx(read);

                // Since we want to output interleaved pairs in stdout mode
                // we need to take care to write them in the correct order
                if(r_idx == 0) {
                    *output[0] << ostr.str();
                    *output[1] << position->second;
                } else {
                    *output[0] << position->second;
                    *output[1] << ostr.str();
                }
                unPaired.erase(position);
            }
        }

    } while (bam_read1(sam->x.bam, read) > 0);

    //The documentation for bam_read1 says that it returns the number of
    //bytes read - which is true, unless it doesn't read any.  It returns
    //-1 for normal EOF and -2 for unexpected EOF.  So don't just wait for
    //it to return 0...
    bam_destroy1(read);
    samclose(sam);

    // Write the remaining unpaired file to the single-end file
    for(map<string, string>::iterator iter = unPaired.begin(); 
            iter != unPaired.end(); ++iter) {
        *output[2] << iter->second;
    }

    // Clean up filehandles
    for(size_t i = 0; i < output.size(); ++i) {
        if(output[i] != &std::cout)
            delete output[i];
    }

    if (print_msgs) {
        cerr << all_seen << " sequences in the BAM file" << endl;
        cerr << exported << " sequences exported" << endl;
    }
}

int main (int argc, char *argv[]) {
    bases[1] = 'A';
    bases[2] = 'C';
    bases[4] = 'G';
    bases[8] = 'T';
    bases[15] = 'N';
    //Complements (original index + 16)
    bases[17] = 'T';
    bases[18] = 'G';
    bases[20] = 'C';
    bases[24] = 'A';
    bases[31] = 'N';
    string output_template("s_%#_sequence.txt");
    int ch;
    while ((ch = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1)
        switch (ch) {
            case 'v':
                cerr << "bam2fastq v" << version << endl;
                exit(0);
            case 'h' :
                usage(0);
            case 'o' :
                output_template = string(optarg);
                break;
            case 'f' :
                overwrite_files = 1;
                break;
            case 'q' :
                print_msgs = 0;
                break;
            case 's' :
                strict = 1;
                break;
            case '?' : //Unrecognized option
                usage(2);
            //The remaining options will set the appropriate variable themselves
        };
    argc -= optind;
    argv += optind;
    if (argc == 0)
        usage(1);
    parse_bamfile(argv[0], output_template);
    return 0;
}
