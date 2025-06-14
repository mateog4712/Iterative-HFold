// Iterative HFold files
#include "hotspot.hh"
#include "Result.hh"
#include "cmdline.hh"
#include "W_final.hh"
#include "h_globals.hh"
// a simple driver for the HFold
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>
#include <string>
#include <getopt.h>

bool exists (const std::string path) {
  struct stat buffer;   
  return (stat (path.c_str(), &buffer) == 0); 
}

void get_input(std::string file, std::string &sequence, std::string &structure){
	if(!exists(file)){
		std::cout << "Input file does not exist" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::ifstream in(file.c_str());
	std::string str;
	int i = 0;
	while(getline(in,str)){
		if(str[0] == '>') continue;
		if(i==0) sequence = str;
		if(i==1) structure = str;
		++i;
	}
	in.close();
}

//check length and if any characters other than ._()
void validateStructure(std::string &seq, std::string &structure){
	int n = structure.length();
	std::vector<int> pairs;
	for(int j = 0; j<n;++j){
		if(structure[j] == '(') pairs.push_back(j);
		if(structure[j] == ')'){
			if(pairs.empty()){
				std::cout << "Incorrect input: More left parentheses than right" << std::endl;
				exit(0);
			}
			else {
				int i = pairs.back();
				pairs.pop_back();
				if(seq[i] == 'A' && seq[j] == 'U'){}
				else if (seq[i] == 'C' && seq[j] == 'G'){}
				else if ((seq[i] == 'G' && seq[j] == 'C') || (seq[i] == 'G' && seq[j] == 'U')){}
				else if ((seq[i] == 'U' && seq[j] == 'G') || (seq[i] == 'U' && seq[j] == 'A')){}
				else{
					std::cout << "Incorrect input: " << seq[i] << " does not pair with " << seq[j] << std::endl;
					exit(0);
				}
			}
		}
	}
	if(!pairs.empty()){
		std::cout << "Incorrect input: More left parentheses than right" << std::endl;
		exit(0);
	}
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
void validateSequence(std::string sequence){

	if(sequence.length() == 0){
		std::cout << "sequence is missing" << std::endl;
		exit(EXIT_FAILURE);
	}
  // return false if any characters other than GCAUT -- future implement check based on type
  for(char c : sequence) {
    if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
		std::cout << "Sequence contains character " << c << " that is not G,C,A,U, or T." << std::endl;
		exit(EXIT_FAILURE);
    }
  }
}


std::string remove_structure_intersection(std::string restricted, std::string structure){
	cand_pos_t length = structure.length();
	for(cand_pos_t i=0; i< length; ++i){
		if(restricted[i] == '(' || restricted[i] == ')') structure[i] = '.';
		
		if (structure[i] == '[') structure[i] = '(';
		
		if (structure[i] == ']') structure[i] = ')';

		if(restricted[i] == 'x') structure[i] = 'x';
	}
	return structure;
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

cand_pos_t paired_structure(cand_pos_t i, cand_pos_t j, std::vector<cand_pos_t> &pair_index){
	cand_pos_t n = pair_index.size();
	return (i >= 0 && j < n && (pair_index[i] == j));
}

/**
 * @brief Takes the input constraint structure as the base output structure. Finds the stacking bases, bulges of size 1, internal loops of size 1x1, 2x1, and 1x2 on the pk_free structure
 * which have these substructures forming on the base pairs and adds them to the output structure.
 * 
 * @param restricted input constraint structure used as the base output structure
 * @param pkfree_structure structure post pseudoknot-free base pair filling.
*/
std::string obtainRelaxedStems(std::string restricted, std::string pkfree_structure){
	cand_pos_t n = restricted.length();

	//Gresult <- G1
	std::string relaxed = restricted;

	cand_pos_t i = 0;
	cand_pos_t j = 0;
	
	std::vector<cand_pos_t> G1_pair;
	std::vector<cand_pos_t> G2_pair;
	G1_pair.resize(n,-2);
	G2_pair.resize(n,-2);
	detect_pairs(restricted,G1_pair);
	detect_pairs(pkfree_structure,G2_pair);

	
	for(int k=0;k<n;++k){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (restricted[i] != pkfree_structure[i])){//if ij not in G1
					//include stacking base pairs
					if(paired_structure(i-1,j+1,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include bulges of size 1
					}else if(paired_structure(i-2,j+1,G1_pair) || paired_structure(i-1,j+2,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x1
					}else if(paired_structure(i-2,j+2,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x2 or 2x1
					}else if( paired_structure(i-3,j+2,G1_pair) || paired_structure(i-2,j+3,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					}
				}
			}
		}
	}

	for(int k=n-1;k>=0;--k){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (restricted[i] != pkfree_structure[i])){//if ij not in G1
					//include stacking base pairs
					if(paired_structure(i+1,j-1,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
						
					//include bulges of size 1
					}else if(paired_structure(i+1,j-2,G1_pair) || paired_structure(i+2,j-1,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x1
					}else if(paired_structure(i+2,j-2,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x2 or 2x1
					}else if(paired_structure(i+2,j-3,G1_pair) || paired_structure(i+3,j-2,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					}
				}
			}
		}
	}
	return relaxed;
}
void seqtoRNA(std::string &sequence){
	bool DNA = false;
    for (char &c : sequence) {
      	if (c == 'T' || c == 't') {
			c = 'U';
			DNA = true;
		}
    }
	noGU = DNA;
}
std::string remove_x(std::string structure){
	for(char &c: structure) {if(c == 'x') c = '.';}
	return structure; 
}

std::string hfold(std::string seq,std::string res, energy_t &energy, bool pk_free, bool pk_only, int dangles){
	sparse_tree tree(res,res.length());
	W_final min_fold(seq,res, pk_free, pk_only,dangles);
	energy = min_fold.hfold(tree);
    std::string structure = min_fold.structure;
    return structure;
}

std::string method2(std::string &seq, std::string &restricted, energy_t &method2_energy, int dangles){

	std::string pk_only_output = hfold(seq,restricted,method2_energy,false,true,dangles);
	std::string pk_free_removed = remove_structure_intersection(restricted,pk_only_output);
	std::string no_x_restricted = remove_x(restricted);

	if(pk_only_output != no_x_restricted) return hfold(seq,pk_free_removed,method2_energy,false,false,dangles);
	else return pk_only_output;
}


int main (int argc, char *argv[])
{
    args_info args_info;

	// get options (call getopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string seq;
	if (args_info.inputs_num>0) {
	seq=args_info.inputs[0];
	} else {
		if(!args_info.input_file_given) std::getline(std::cin,seq);
	}	

	std::string restricted;
    args_info.input_structure_given ? restricted = input_struct : restricted = "";

	std::string fileI;
    args_info.input_file_given ? fileI = input_file : fileI = "";

	std::string fileO;
    args_info.output_file_given ? fileO = output_file : fileO = "";	

	int number_of_suboptimal_structure = args_info.subopt_given ? subopt : 1;

	int dangles = args_info.dangles_given ? dangle_model : 2;

	if(fileI != ""){
		
		if(exists(fileI)){
			get_input(fileI,seq,restricted);
		}
		if(seq == ""){
			std::cout << "sequence is missing from file" << std::endl; 
		}
		
	}
	std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	
	validateSequence(seq);
	if(restricted != "") validateStructure(seq,restricted);

	std::string file= args_info.paramFile_given ? parameter_file : "params/rna_DirksPierce09.par";
	if(exists(file)){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}
	else if (seq.find('T') != std::string::npos){
		vrna_params_load_DNA_Mathews2004();
	}

	std::vector<Hotspot> hotspot_list;

	// Hotspots

	vrna_param_s *params;
	params = scale_parameters();
	if(restricted != ""){
		Hotspot hotspot(1,restricted.length(),restricted.length()+1);
		hotspot.set_structure(restricted);
		hotspot_list.push_back(hotspot);
	}
	if((number_of_suboptimal_structure-hotspot_list.size())>0) {
		get_hotspots(seq, hotspot_list,number_of_suboptimal_structure,params);
	}
	free(params);

	// Data structure for holding the output
	std::vector<Result> result_list;

	// Iterate through all hotspots or the single given input structure
	cand_pos_t size = hotspot_list.size();
	for(int i = 0;i<size;++i){
		energy_t method1_energy = INF;
		energy_t method2_energy = INF;
		energy_t method3_energy = INF;
		energy_t method4_energy = INF;
		energy_t final_en = INF;
		int method_chosen = -1;
		std::string final_structure;
		std::string res = hotspot_list[i].get_structure();

		//Method1
		std::string method1_structure = hfold(seq,res,method1_energy,false,false,dangles);
		if(method1_energy < final_en){
		final_en = method1_energy;
		final_structure=method1_structure;
		method_chosen = 1;
		}

		//Method2
		std::string method2_structure = method2(seq,res,method2_energy,dangles);
		if(method2_energy < final_en){
		final_en = method2_energy;
		final_structure=method2_structure;
		method_chosen = 2;
		}

		//Method3
		std::string pk_free = hfold(seq,res,method3_energy,true,false,dangles);
		std::string relaxed = obtainRelaxedStems(res,pk_free);
		cand_pos_t n = res.length();
		for(int i =0; i< n;++i) if(res[i] == 'x') relaxed[i] = 'x';
		std::string method3_structure = method2(seq,relaxed,method3_energy,dangles);
		if(method3_energy < final_en){
			final_en = method3_energy;
			final_structure=method3_structure;
			method_chosen = 3;
		}




		//Method4
		std::vector< std::pair<int,int> > disjoint_substructure_index;
		find_disjoint_substructure(res,disjoint_substructure_index);
		std::string disjoint_structure = res;
		for(auto current_substructure_index : disjoint_substructure_index){
			cand_pos_t i = current_substructure_index.first;
			cand_pos_t j = current_substructure_index.second;
			energy_t energy = INF;

			std::string subsequence = seq.substr(i,j-i+1);
			std::string substructure = res.substr(i,j-i+1);

			std::string pk_free = hfold(subsequence,substructure,energy,true,false,dangles);
			std::string relaxed = obtainRelaxedStems(substructure,pk_free);
			cand_pos_t sub_n = substructure.length();
			for(int i =0; i< sub_n;++i) if(substructure[i] == 'x') relaxed[i] = 'x';
			disjoint_structure.replace(i,j-i+1,relaxed);
		}
		std::string method4_structure = method2(seq,disjoint_structure,method4_energy,dangles);
		if(method4_energy < final_en){
			final_en = method4_energy;
			final_structure=method4_structure;
			method_chosen = 4;
		}
		double final_energy = final_en/100.0;

		if(!args_info.input_structure_given && final_energy>0.0){
			final_energy = 0.0;
			final_structure = std::string(n,'.');
		}
		
		Result result(seq,hotspot_list[i].get_structure(),hotspot_list[i].get_energy(),final_structure,final_energy,method_chosen);
		result_list.push_back(result);
	}

    

	Result::Result_comp result_comp;
	std::sort(result_list.begin(), result_list.end(),result_comp );

	int number_of_output = 1;

	if(number_of_suboptimal_structure != 1){
			number_of_output = std::min( (int) result_list.size(),number_of_suboptimal_structure);
	}

	//output to file
	if(fileO != ""){
		std::ofstream out(fileO);
		out << seq << std::endl;
		std::cout << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
		std::cout << "Result_" << 0 << ":     " << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;
		for (int i=1; i < number_of_output; i++) {
			if(result_list[i].get_final_structure() == result_list[i-1].get_final_structure()) continue;
			std::cout << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
			std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
			if(args_info.verbose_given){
				std::cout << "Method: " << result_list[i].get_method_chosen() << std::endl;
			}
			
		}

	}else{
		//kevin: june 22 2017
		//Mateo: Sept 13 2023
		//changed format for ouptut to stdout
		std::cout << seq << std::endl;
		if(result_list.size() == 1){
			// std::cout << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
			std::cout << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;
			if(args_info.verbose_given){
				std::cout << "Method: " << result_list[0].get_method_chosen() << std::endl;
			}
		}
		else{
			std::cout << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
			std::cout << "Result_" << 0 << ":     " << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;
			for (int i=1; i < number_of_output; i++) {
				if(result_list[i].get_final_structure() == result_list[i-1].get_final_structure()) continue;
				std::cout << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
				std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
				if(args_info.verbose_given){
					std::cout << "Method: " << result_list[i].get_method_chosen() << std::endl;
				}
				
			}
		}
	}
	cmdline_parser_free(&args_info);

    return 0;
}
