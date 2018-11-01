#include <stdio.h>
#include <stdlib.h>


#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <string>


typedef std::unordered_map<std::string, std::string> readname_seq;
//typedef std::map<std::string, std::string> readname_seq;

int main(int argc, char *argv[])
{
    if (argc != 3) {printf("Usage: %s <readname file by sort> <fasta>\n\n", argv[0]); return -1;}

    readname_seq reads;
    const char *readfile = argv[1];
    const char *fastfile = argv[2];

    std::fstream infa(fastfile);
    std::string line, name, seqs;

    while(getline(infa, line))
    {
        const char *cline = line.c_str();
        if ('>' == cline[0]) {
            if (seqs.size() > 0) reads.insert(std::make_pair(name, seqs)); 
            name = line;
            seqs.clear();
        } else {
            seqs.append(line).append("\n");
        }
    }
    //last read
    reads.insert(std::make_pair(name, seqs));
    infa.close();
    std::cout<<"Read fasta done!"<<std::endl;

    std::fstream inrd(readfile);
    std::ofstream outf(std::string(fastfile).append(".sort"));
    while(getline(inrd, line))
    {
        auto it = reads.find(line);
        if (reads.end() != it) {
            outf.write(line.c_str(), line.size());
            outf.write("\n", 1);
            outf.write(it->second.c_str(), it->second.size());
            reads.erase(it);
        }
    }
    inrd.close();
    outf.close();


    return 0;
}
