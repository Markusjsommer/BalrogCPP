#include <iostream>
#include <cctype>
#include "cxxopts.hpp"
#include "FastaReader.h"
#include "GeneFinder.h"

int main(int argc, char* argv[]) {
    // parse command line options
    cxxopts::Options options("Balrog", "Balrog is a prokaryotic gene finder based on a temporal convolutional network");
    options.add_options()
            ("i, in", "Path to input fasta or gzipped fasta", cxxopts::value<std::string>())
            ("o, out", "Path to output annotation", cxxopts::value<std::string>())
            ("g, gene-model", "Path to pretrained LibTorch gene model", cxxopts::value<std::string>())
            ("t, TIS-model", "Path to pretrained LibTorch TIS model", cxxopts::value<std::string>())
            ("max-overlap", "Maximum allowable overlap between genes in nucleotides", cxxopts::value<int>()->default_value("60"))
            ("min-length", "Minimum allowable gene length in nucleotides", cxxopts::value<int>()->default_value("60"))
            ("table", "Nucleotide to amino acid translation table. 11 for most bacteria/archaea, 4 for Mycoplasma/Spiroplasma.",
                    cxxopts::value<int>()->default_value("11"))
            ("max-connections", "Maximum number of forward connections in the directed acyclic graph used to find a set of coherent genes in each genome.",
                    cxxopts::value<int>()->default_value("50"))
            ("gene-batch-size", "Batch size for the temporal convolutional network used to score genes.",
                    cxxopts::value<int>()->default_value("200"))
            ("TIS-batch-size", "Batch size for the temporal convolutional network used to score TIS.",
                    cxxopts::value<int>()->default_value("5000"))
            ("verbose", "Verbose output", cxxopts::value<bool>()->default_value("true"))
            ("h,help", "Print Balrog usage")
            ;
    auto result = options.parse(argc, argv);

    // check validity and display help
    if (result.count("help") or not result.count("in")){
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // check translation table
    int table = result["table"].as<int>();
    if (table != 11 and table != 4){
        std::cout << "Only translation tables 11 and 4 are currently implemented. Please open a GitHub issue if you need another." << std::endl;
        exit(0);
    }

    // load LibTorch jit traced gene model
    if (result["verbose"].as<bool>()){
        std::cout << "Importing gene model..." << std::endl;
    }
    std::string model_path = result["gene-model"].as<std::string>();
    torch::jit::script::Module gene_model;
    try {
        gene_model = torch::jit::load(model_path);
    }
    catch (const c10::Error& e) {
        std::cerr << "error loading the LibTorch gene model\n";
        return -1;
    }

    // load LibTorch jit traced TIS model
    if (result["verbose"].as<bool>()){
        std::cout << "Importing TIS model..." << std::endl;
    }
    std::string model_path_TIS = result["TIS-model"].as<std::string>();
    torch::jit::script::Module TIS_model;
    try {
        TIS_model = torch::jit::load(model_path_TIS);
    }
    catch (const c10::Error& e) {
        std::cerr << "error loading the LibTorch TIS model\n";
        return -1;
    }

    // read fasta
    if (result["verbose"].as<bool>()){
        std::cout << "Reading fasta..." << std::endl;
    }
    std::vector<std::string> seq_vec;
    std::vector<std::string> contigname_vec;

    std::string in_path = result["in"].as<std::string>();
    FastaReader io;
    io.read_fasta(in_path, seq_vec, contigname_vec);

    // capitalize all nucelotides
    for (auto &seq : seq_vec){
        for(auto &c: seq){
            c = toupper(c);
        }
    }

    // find genes on all contigs
    std::vector<std::vector<std::pair<int, int>>> gene_coord_vec;
    std::vector<std::vector<bool>> gene_strand_vec;
    std::vector<std::vector<std::string>> gene_nucseq_vec;
    std::vector<std::vector<std::string>> gene_protseq_vec;

    int i = 0;
    for (std::string seq : seq_vec){
        ++i;

        GeneFinder gf(gene_model, TIS_model);
        if (result["verbose"].as<bool>()) {
            std::cout << std::endl << "contig " << i << " of " << seq_vec.size() << " : length " << seq.length() << std::endl;
        }

        std::vector<std::pair<int, int>> gene_coord;
        std::vector<bool> gene_strand;
        std::vector<std::string> gene_nucseq;
        std::vector<std::string> gene_protseq;


        gf.find_genes(seq,
                      gene_coord,
                      gene_strand,
                      gene_nucseq,
                      gene_protseq,
                      table,
                      result["min-length"].as<int>(),
                      result["max-overlap"].as<int>(),
                      result["verbose"].as<bool>(),
                      result["gene-batch-size"].as<int>(),
                      result["TIS-batch-size"].as<int>());


        gene_coord_vec.emplace_back(gene_coord);
        gene_strand_vec.emplace_back(gene_strand);
        gene_nucseq_vec.emplace_back(gene_nucseq);
        gene_protseq_vec.emplace_back(gene_protseq);
    }

    // export genes to gff annotation file
    if (result["verbose"].as<bool>()) {
        std::cout << "Writing predicted genes to " << result["out"].as<std::string>() << std::endl;
    }
    std::ofstream out;
    out.open(result["out"].as<std::string>());
    out << "##gff-version 3" << std::endl;
    std::string contigname;
    for (int j=0; j < contigname_vec.size(); ++j) {
        // extract name up to first space
        contigname = contigname_vec[j].substr(0, contigname_vec[j].find(' '));
        contigname.erase(std::remove(contigname.begin(), contigname.end(), '>'), contigname.end());
        // write all 1-indexed sequence region names and coords
        out << "##sequence-region " << contigname << " " << 1 << " " << seq_vec[j].length() << std::endl;
    }
    for (int j=0; j < contigname_vec.size(); ++j) {
        contigname = contigname_vec[j].substr(0, contigname_vec[j].find(' '));
        contigname.erase(std::remove(contigname.begin(), contigname.end(), '>'), contigname.end());
        // write CDS features
        for (int k=0; k < gene_coord_vec[j].size(); ++k){
            int coord0;
            int coord1;
            if (gene_strand_vec[j][k]){
                coord0 = gene_coord_vec[j][k].first + 1;
                coord1 = gene_coord_vec[j][k].second + 3;
            } else{
                coord0 = gene_coord_vec[j][k].second + 1;
                coord1 = gene_coord_vec[j][k].first + 3;
            }
            out << contigname << "\tBalrog\tCDS\t" << coord0 << "\t" << coord1 << "\t" << "." << "\t";
            if (gene_strand_vec[j][k]){
                out << "+\t";
            } else{
                out << "-\t";
            }
            out << "0\tinference=ab initio prediction:Balrog;product=hypothetical protein" << std::endl;
        }
    }
    out.close();

    if (result["verbose"].as<bool>()) {
        std::cout << "Done...\n" << std::endl;
    }

    return 0;
}
