// Iterative HFold files
#include "HFold.hh"
#include "cmdline.hh"
#include "pseudo_loop.hh"
#include "h_globals.hh"
#include "cmdline.hh"
#include "Result.hh"
#include "Hotspot.hh"
// a simple driver for the HFold
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>

// filesystem::exists not supported in older macOS which is needed for the Conda Build
bool file_exists(const std::string& name) { 
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

std::string getSequence(args_info a) {
	if (a.inputs_num > 0) {
		return a.inputs[0];
	} else if (!a.input_file_given) {
		std::cout << "Sequence: ";
		std::string seq;
		std::getline(std::cin, seq);
		return seq;
	}
    return "";
}

void trim(std::string& s) {
    /**
     * @brief Trims leading and trailing whitespace from a string.
     * 
     * This function modifies the input string in-place to remove any leading
     * and trailing whitespace characters, including spaces, tabs, newlines,
     * carriage returns, form feeds, and vertical tabs.
     * 
     * @param s The string to be trimmed.
     * @return void
     */
    s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
    s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
}

std::vector<RNAEntry> get_all_file_entries(const std::string& file){
    if(!file_exists(file)){
        std::cerr << "Error: Input file not found: " << file << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // state machine to parse the file
    #define UNINITIALIZED -1
    #define NAME 0
    #define SEQUENCE 1
    #define STRUCTURE 2

    std::ifstream in(file.c_str());
    std::string line;
    RNAEntry current;
    std::vector<RNAEntry> entries;
    int state = UNINITIALIZED;
    int line_number = 0;

    while(getline(in, line)){
        ++line_number;
        trim(line);
        if (line.empty()) continue;

        // Check if the line is the name of the entry
        if ((state == NAME || state == UNINITIALIZED) && (line[0] != '>')) {
            std::cerr << "Error: Expected '>' at the beginning of the line: " << line << ". Line number: " << line_number <<std::endl;
            exit(EXIT_FAILURE);
        }

        if (line[0] == '>'){

            if (!current.name.empty() && !current.sequence.empty() && current.structure.empty()) {
                current.structure = "";
            }

            if (state != UNINITIALIZED) { // valid entry, save it
                if (current.sequence.empty() && current.structure.empty()) {
                    std::cerr << "Warning: Sequence and structure are empty for entry: " << current.name << ". Line number: " << line_number << ". Skipping..."<<  std::endl;
                }
                entries.push_back(current);
                current = RNAEntry();
            }

            current.name = line.substr(1);
            state = SEQUENCE;

        } else if (state == SEQUENCE){
            if (!validateSequence(line, false)) {
                std::cerr << "Error: Sequence is invalid for entry: " << current.name  << ". Line number: " << line_number << std::endl;
                exit(EXIT_FAILURE);
            }
            current.sequence = line;
            state = STRUCTURE;
            
        } else if (state == STRUCTURE) {
            if (!validateStructure(current.sequence, line, false)) {
                std::cerr << "Error: Structure is invalid for entry: " << current.name << ". Line number: " << line_number  << std::endl;
                exit(EXIT_FAILURE);
            }
            current.structure = line;
            state = NAME;
        } else {
            // Should never reach here
            std::cerr << "Error: Unexpected state at line " << line_number << ": " << line << ". Line number: " << line_number  << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // Handle the last entry
    if (!current.name.empty() && !current.sequence.empty()) {
        if (current.structure.empty()) {
            current.structure = "";
        }
        entries.push_back(current);
    }

    return entries;
}

std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq, const std::string& restricted) {
    std::vector<RNAEntry> entries;
    if (!seq.empty()) {
        entries.emplace_back("Console Sequence", seq, restricted);
    }
    if (!fileI.empty()){
        std::vector<RNAEntry> file_entries = get_all_file_entries(fileI);
        entries.insert(entries.end(), file_entries.begin(), file_entries.end());
    }
    if (entries.empty()) throw std::runtime_error("Sequence is missing");
    return entries;
}

// check length and if any characters other than ._() and if the base pairs are valid
bool validateStructure(std::string &seq, std::string &structure, bool exit_on_invalid) {
    cand_pos_t n = structure.length();
    std::vector<int> pairs;
    for (int j = 0; j < n; ++j) {
        if (structure[j] == '(') pairs.push_back(j);
        if (structure[j] == ')') {
            if (pairs.empty()) {
                if (exit_on_invalid) {
                    std::cerr << "Error: Incorrect input: More right parentheses than left" << std::endl;
                    exit(EXIT_FAILURE);
                }
                return false;
            } else {
                cand_pos_t i = pairs.back();
                pairs.pop_back();
                if (seq[i] == 'A' && seq[j] == 'U') {
                } else if (seq[i] == 'C' && seq[j] == 'G') {
                } else if ((seq[i] == 'G' && seq[j] == 'C') || (seq[i] == 'G' && seq[j] == 'U')) {
                } else if ((seq[i] == 'U' && seq[j] == 'G') || (seq[i] == 'U' && seq[j] == 'A')) {
                } else {
                    if (exit_on_invalid) {
                        std::cerr << "Error: Incorrect input: " << seq[i] << " does not pair with " << seq[j] << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    return false;
                }
            }
        }
    }
    if (!pairs.empty()) {
        if (exit_on_invalid) {
            std::cerr << "Error: Incorrect input: More left parentheses than right" << std::endl;
            exit(EXIT_FAILURE);
        }
        return false;
    }
    return true;
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
bool validateSequence(std::string& sequence, bool exit_on_invalid){ 
	if(sequence.length() == 0){
        if (exit_on_invalid) {
            std::cerr << "Error: Sequence is missing." << std::endl;
            exit(EXIT_FAILURE);
        }
        return false;
	}

  // return false if any characters other than GCAUT -- future implement check based on type
  for(char c : sequence) {
    if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T' || c == 'N')) {
        if (exit_on_invalid) {
            std::cerr  << "Error: Sequence contains invalid character " << c << ". Allowed: G, C, A, U, T, N." << std::endl;
            exit(EXIT_FAILURE);
        }
        return false;
    }
  }
    return true;
}

std::string hfold(std::string seq,std::string res, double &energy, SHAPEData &ShapeData, bool pk_free, bool pk_only, int dangle){
	sparse_tree tree(res, seq.size());
    pseudo_loop min_fold(seq, res,tree,ShapeData, pk_free, pk_only, dangle);
    energy = min_fold.hfold();
    std::string structure = min_fold.structure;
    return structure;
}

void seqtoRNA(std::string &sequence){
	/**
	 * @brief Converts a DNA sequence to RNA by replacing T with U.
	 * 
	 * This function modifies the input string in-place and logs the 
	 * conversion if any T's are found.
	 * 
	 * @param sequence The DNA sequence to be converted (e.g. "ATCGT").
	 *                 The result (RNA) will be stored in the same string.
	 */
	bool isRNA = true;
	std::string original_sequence = sequence;
    for (char &c : sequence) {
      	if (c == 'T') {
			c = 'U';
			isRNA = false;
		}
    }
	if(!isRNA){
		std::cout << "Input sequence contains T's, converting to U's" << std::endl;
		std::cout << "Original sequence: " << original_sequence << std::endl;
		std::cout << "Converted sequence: " << sequence << std::endl;
	}
}

void preprocess_sequence(std::string& seq, std::string& restricted, bool noConv) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper); // convert sequence to uppercase
    if (!noConv) seqtoRNA(seq);
    validateSequence(seq);
    if (!restricted.empty()) validateStructure(seq, restricted);
}

void load_energy_parameters(const std::string& paramFile, const std::string& seq, bool param_given) {
    if (paramFile.empty()) return; // No parameter file provided
    
    if(param_given){
        if (file_exists(paramFile)) {
            vrna_params_load(paramFile.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
        } else {
            std::cerr << "Error: Input file not found: " << paramFile << std::endl;
        }
    } else{
        if (seq.find('T') != std::string::npos) {
            vrna_params_load_DNA_Mathews2004();
        } else if (file_exists(paramFile)){
            vrna_params_load(paramFile.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
        } else{
            std::cerr << "Error: Input file not found: " << paramFile << std::endl;
        }
    }
}

std::vector<Hotspot> build_hotspots(const std::string& seq, const std::string& restricted, SHAPEData &ShapeData, int suboptCount) {
    std::vector<Hotspot> hotspots;
    vrna_param_s* params = vrna_params(NULL);

    if (!restricted.empty()) {
        hotspots.emplace_back(1, restricted.length(), restricted.length() + 1);
        hotspots.back().set_structure(restricted);
    }
    if (static_cast<int>(hotspots.size()) < suboptCount) {
        get_hotspots(seq, hotspots,ShapeData, suboptCount, params);
    }

    free(params);
    return hotspots;
}

std::vector<Result> fold_hotspots(
    const std::string& seq, 
    const std::vector<Hotspot>& hotspots, SHAPEData &ShapeData,
    int dangles, bool input_structure_given
) {
    std::vector<Result> results;
    results.reserve(hotspots.size()); // Pre-allocate memory for results

    for (const Hotspot& hs : hotspots) {
        std::string final_structure = ""; double final_energy = 0.0; 
        double method1_energy = 0.0, method2_energy = 0.0, method3_energy = 0.0, method4_energy = 0.0;
        int method_chosen = -1;
        std::string method1_structure = hfold(seq, hs.get_structure(), method1_energy = 0.0,ShapeData,false,false, dangles);
        if(method1_energy < final_energy){
            final_energy = method1_energy;
            final_structure=method1_structure;
            method_chosen = 1;
		}
        std::string method2_structure = method2(seq,hs.get_structure(),method2_energy,ShapeData,dangles);
		if(method2_energy < final_energy){
            final_energy = method2_energy;
            final_structure=method2_structure;
            method_chosen = 2;
		}

        //Method3
		std::string pk_free = hfold(seq,hs.get_structure(),method3_energy,ShapeData,true,false,dangles);
		std::string relaxed = obtainRelaxedStems(hs.get_structure(),pk_free);
		cand_pos_t n = hs.get_structure().length();
		for(cand_pos_t i =0; i< n;++i) if(hs.get_structure()[i] == 'x') relaxed[i] = 'x';
		std::string method3_structure = method2(seq,relaxed,method3_energy,ShapeData,dangles);
		if(method3_energy < final_energy){
			final_energy = method3_energy;
			final_structure=method3_structure;
			method_chosen = 3;
		}

        //Method4
        std::vector< std::pair<int,int> > disjoint_substructure_index;
		find_disjoint_substructure(hs.get_structure(),disjoint_substructure_index);
		std::string disjoint_structure = hs.get_structure();
		for(auto current_substructure_index : disjoint_substructure_index){
			cand_pos_t i = current_substructure_index.first;
			cand_pos_t j = current_substructure_index.second;
			double energy = INF;

			std::string subsequence = seq.substr(i,j-i+1);
			std::string substructure = hs.get_structure().substr(i,j-i+1);

			std::string pk_free = hfold(subsequence,substructure,energy,ShapeData,true,false,dangles);
			std::string relaxed = obtainRelaxedStems(substructure,pk_free);
			cand_pos_t sub_n = substructure.length();
			for(cand_pos_t i =0; i< sub_n;++i) if(substructure[i] == 'x') relaxed[i] = 'x';
			disjoint_structure.replace(i,j-i+1,relaxed);
		}
		std::string method4_structure = method2(seq,disjoint_structure,method4_energy,ShapeData,dangles);
		if(method4_energy < final_energy){
			final_energy = method4_energy;
			final_structure=method4_structure;
			method_chosen = 4;
		}
        std::cout << method1_structure << " " << method1_energy << std::endl << method2_structure << " " << method2_energy << std::endl << method3_structure << " " << method3_energy << std::endl << method4_structure << " " << method4_energy << std::endl;


		if (!input_structure_given && final_energy > 0.0) {
			final_energy = 0.0;
			final_structure = std::string(seq.length(), '.');
		}
        results.emplace_back(seq, hs.get_structure(), hs.get_energy(), final_structure, final_energy, method_chosen);
    }

    std::sort(results.begin(), results.end(), Result::Result_comp{});
    return results;
}

void output_results(
    const std::string& seq,
    const std::vector<Result>& results,
    const std::string& fileO,
    int suboptCount,
    const std::string& name,
    const std::size_t input_count,
    const bool skip_duplicates,
    const bool verbose
) {
    int number_of_output = std::min(results.size(), static_cast<std::size_t>(suboptCount));

    if (!fileO.empty()) { //output to file
        std::ofstream out(fileO, std::ios::app);
        if (!name.empty()){
            out << ">" << name << std::endl; 
        }
        out << seq << std::endl;
        if (number_of_output == 1) {
            out << "Restricted" << ": " << results[0].get_restricted() << std::endl;
            out << "Result" << ":     " << results[0].get_final_structure()
                << " (" << results[0].get_final_energy() << ")" << std::endl;
        } else {
            for (int i = 0; i < number_of_output; i++) {
                if (skip_duplicates && i && results[i].get_final_structure() == results[i - 1].get_final_structure()){ // skip duplicates
                    continue;
                }
                out << "Restricted_" << i << ": " << results[i].get_restricted() << std::endl;
                out << "Result_" << i << ":     " << results[i].get_final_structure()
                    << " (" << results[i].get_final_energy() << ")" << std::endl;
                if(verbose){
				    out << "Method: " << results[i].get_method_chosen() << std::endl;
			    }
            }
        }
    } else { // output to console
        if (!name.empty() && input_count > 1){
            std::cout << ">" << name << std::endl; 
        }
        
        if (results.size() == 1) { 
            std::cout << "Sequence:        " << seq << std::endl;
			std::cout << "Restricted:      " << results[0].get_restricted() << std::endl;
            std::cout << "Final Structure: " << results[0].get_final_structure()
                      << " (" << results[0].get_final_energy() << ")" << std::endl;
            if(verbose){
				std::cout << "Method: " << results[0].get_method_chosen() << std::endl;
			}
        } else {
            int alignment = std::floor(std::log10(number_of_output));
            std::cout << "Sequence:     " << std::string(alignment, ' ') << seq << std::endl;
            for (int i = 0; i < number_of_output; i++) {
                alignment = std::floor(std::log10(number_of_output)) - (std::floor(std::log10(i !=0 ? i : 1)));
                if (skip_duplicates && i && results[i].get_final_structure() == results[i - 1].get_final_structure()){ // skip duplicates
                    continue;
                }
                std::cout << "Restricted_" << i << ": " << std::string(alignment, ' ') << results[i].get_restricted() << std::endl;
                std::cout << "Result_" << i << ":     " << std::string(alignment, ' ') << results[i].get_final_structure()
                          << " (" << results[i].get_final_energy() << ")" << std::endl;
                if(verbose){
				    std::cout << "Method: " << results[0].get_method_chosen() << std::endl;
			    }
            }
        }
    }
}
// Iterative HFold Util
std::string remove_structure_intersection(const std::string &restricted, std::string structure){
	assert(restricted.length() == structure.length());
    cand_pos_t length = structure.length();
    std::replace(structure.begin(), structure.end(), '(', '.');
    std::replace(structure.begin(), structure.end(), ')', '.');
	for(cand_pos_t i=0; i< length; ++i){
		if(restricted[i] == 'x') structure[i] = 'x';
	}
    std::replace(structure.begin(), structure.end(), '[', '(');
    std::replace(structure.begin(), structure.end(), ']', ')');
	return structure;
}
std::string remove_x(std::string structure){
	std::replace(structure.begin(), structure.end(), 'x', '.');
	return structure; 
}
std::string method2(const std::string &seq, const std::string &restricted, double &method2_energy, SHAPEData &ShapeData, int dangles){

	std::string pk_only_output = hfold(seq,restricted,method2_energy,ShapeData,false,true,dangles);
	std::string pk_free_removed = remove_structure_intersection(restricted,pk_only_output);
	std::string no_x_restricted = remove_x(restricted);

	if(pk_only_output != no_x_restricted) return hfold(seq,pk_free_removed,method2_energy,ShapeData,false,false,dangles);
	else return pk_only_output;
}
/**
 * @brief Fills the pair array
 * p_table will contain the index of each base pair
 * X or x tells the program the base cannot pair and . sets it as unpaired but can pair
 * @param structure Input structure
 * @param p_table Restricted array
 */
void detect_pairs(const std::string &structure, std::vector<cand_pos_t> &p_table){
	cand_pos_t i, j, length = structure.length();
	std::vector<cand_pos_t>  pairs;
	pairs.push_back(length);

	for (i=length-1; i >=0; --i){
		if ((structure[i] == 'x') || (structure[i] == 'X'))
			p_table[i] = -1;
		else if (structure[i] == '.')
			p_table[i] = -2;
		if (structure[i] == ')'){
			pairs.push_back(i);
		}
		if (structure[i] == '('){
			j = pairs[pairs.size()-1];
			pairs.erase(pairs.end()-1);
			p_table[i] = j;
			p_table[j] = i;
		}
	}
	pairs.pop_back();
	if (pairs.size() != 0) std::cout << pairs[0] << std::endl << structure << std::endl;
	if (pairs.size() != 0)
	{
		fprintf (stderr, "The given structure is not valid: more left parentheses than right parentheses: \n");
		exit (1);
	}
}
// /**
//  * @brief returns a vector of pairs which represent the start and end indices for each disjoint substructure in the structure
//  * 
//  * @param CL_ Candidate list
//  * @return total number of candidates
//  */
void find_disjoint_substructure(std::string structure, std::vector< std::pair<int,int> > &pair_vector){
	cand_pos_t n = structure.length();
	cand_pos_t count = 0;
	cand_pos_t i = 0;
	for(cand_pos_t k=0; k<n;++k){
		if(structure[k] == '('){
			if(count == 0) i = k;
			count++;

		}else if(structure[k] == ')'){
			count--;
			if(count == 0){
				std::pair <int,int> ij_pair (i,k);
				pair_vector.push_back(ij_pair);
			}
		}
	}
}
// Loop-size offsets (a, b)
const std::vector<std::pair<cand_pos_t, cand_pos_t>> kLoopOffsets = {
    {1, 1},         // stacked pair
    {2, 1}, {1, 2}, // bulge of size 1
    {2, 2},         // 1x1 interior loop
    {3, 2}, {2, 3}, // 1x2 / 2x1 interior loop
    {3, 3},         // 2x2 interior loop
    {4, 3}, {3, 4}  // 3x2 / 2x3 interior loop
};
// sign = +1 looks outward from (i, j) (i-a, j+b); sign = -1 looks inward (i+a, j-b).
bool hasNeighborInG1(cand_pos_t i, cand_pos_t j, int sign, const std::vector<cand_pos_t> &G1_pair){
    for (const auto &[a, b] : kLoopOffsets) {
        if (paired_structure(i - sign * a, j + sign * b, G1_pair)) return true;
    }
    return false;
}
void tryRelaxPair(cand_pos_t i, cand_pos_t j, int sign, const std::string &restricted, const std::string &pkfree_structure, std::string &relaxed, std::vector<cand_pos_t> &G1_pair) {
    if (restricted[i] == pkfree_structure[i]) return; // already in G1
    if (!hasNeighborInG1(i, j, sign, G1_pair)) return;

    relaxed[i] = pkfree_structure[i];
    relaxed[j] = pkfree_structure[j];
    G1_pair[i] = j;
    G1_pair[j] = i;
}
/**
 * @brief Takes the input constraint structure as the base output structure. Finds the stacking bases, bulges of size 1, internal loops of size 1x1, 2x1, and 1x2 on the pk_free structure
 * which have these substructures forming on the base pairs and adds them to the output structure.
 * 
 * @param restricted input constraint structure used as the base output structure
 * @param pkfree_structure structure post pseudoknot-free base pair filling.
*/
std::string obtainRelaxedStems(const std::string &restricted, const std::string &pkfree_structure){
	cand_pos_t n = restricted.length();

	//Gresult <- G1
	std::string relaxed = restricted;
	
	std::vector<cand_pos_t> G1_pair;
	std::vector<cand_pos_t> G2_pair;
	G1_pair.resize(n,-2);
	G2_pair.resize(n,-2);
	detect_pairs(restricted,G1_pair);
	detect_pairs(pkfree_structure,G2_pair);

	// Outward sweep: extend stems by checking just outside each G2 pair.
    for (cand_pos_t k = 0; k < n; ++k) {
        if (G2_pair[k] <= -1 || k >= G2_pair[k]) continue; // handle each pair once, as (i < j)
        tryRelaxPair(k, G2_pair[k], /*sign=*/+1, restricted, pkfree_structure, relaxed, G1_pair);
    }
    // Inward sweep: extend stems by checking just inside each G2 pair.
    for (cand_pos_t k = n - 1; k >= 0; --k) {
        if (G2_pair[k] <= -1 || k >= G2_pair[k]) continue;
        tryRelaxPair(k, G2_pair[k], /*sign=*/-1, restricted, pkfree_structure, relaxed, G1_pair);
    }
	return relaxed;
}
cand_pos_t paired_structure(cand_pos_t i, cand_pos_t j, const std::vector<cand_pos_t> &pair_index){
	cand_pos_t n = pair_index.size();
	return (i >= 0 && j < n && (pair_index[i] == j));
}