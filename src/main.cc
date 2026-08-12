#include "Result.hh"
#include "cmdline.hh"
#include "pseudo_loop.hh"
#include "HFold.hh"
#include "Hotspot.hh"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>


int main(int argc, char* argv[]) {
    args_info a;
    if (cmdline_parser(argc, argv, &a) != 0) return 1;

	/* ────── Command‑line → variables ────── */
    std::string sequence   = getSequence(a);
    std::string restricted = a.input_structure_given ? a.input_structure_arg   : ""; // Initial structure
    std::string fileI      = a.input_file_given      ? a.input_file_arg        : "";
    std::string fileO      = a.output_file_given     ? a.output_file_arg       : "";
    int suboptCount        = a.opt_given             ? a.opt_arg               : 1;
    int dangles            = a.dangles_given         ? a.dangles_arg           : 1;
    std::string paramFile  = a.paramFile_given       ? a.paramFile_arg         : std::string(PARAMS_DIR) + "/rna_DirksPierce09.par";
    std::string shapeFile  = a.shape_given           ? a.shape_arg             : "";
    bool verbose           = a.verbose_flag; 

    std::vector<RNAEntry> inputs = get_all_inputs(fileI, sequence, restricted);

    for (RNAEntry& current : inputs){
        preprocess_sequence(current.sequence, current.structure, a.noConv_given);
        load_energy_parameters(paramFile, current.sequence,a.paramFile_given);
        SHAPEData ShapeData(shapeFile,current.sequence.length());

        std::vector<Hotspot> hotspots = build_hotspots(current.sequence, current.structure, ShapeData, suboptCount);

        std::vector<Result>  results  = fold_hotspots(current.sequence, hotspots,ShapeData, dangles, a.input_structure_given);

        output_results(current.sequence, results, fileO, suboptCount, current.name, inputs.size(),false,verbose);
    }

    cmdline_parser_free(&a);
    return 0;
}
