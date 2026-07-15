#include <Rcpp.h>
#include <vector>
#include <string>

using namespace Rcpp;

//' Assign a variant type from ALLELE_LENGTHS field
//' @title Assign a variant type from ALLELE_LENGTHS field
//' @param al_lens a vector of allele lengths formatted as in ALLELE_LENGTHS, e.g. 1,0,10/56
//' @return a vector of characters with the variant types
//' @author Jean Monlong
//' @keywords internal
// [[Rcpp::export]]
CharacterVector assign_type_from_lengths(CharacterVector al_lens){
    int n_snarls = al_lens.size();
    CharacterVector var_types(n_snarls);
    
    // process every input string
    for (int idx = 0; idx < n_snarls; idx++) {
        // quickly handle SNPs (vast majority)
        if(al_lens[idx] == "1,1") {
            var_types[idx] = "SNP";
            continue;
        }
        String al_len_s = al_lens[idx];
        std::string al_len = al_len_s.get_cstring();
        // to remember if any allele was a MNP or SV
        bool any_sv = false;
        bool any_non_snv = false;
        // start looking for the delimiters (very ugly)
        auto f_start = 0U;
        auto f_end = std::min(al_len.find(','), al_len.find('/'));
        int cur_al_len;
        while (f_end != std::string::npos) {
            cur_al_len = std::stoi(al_len.substr(f_start, f_end - f_start).c_str());
            if (cur_al_len >= 50) {
                any_sv = true;
                break;
            } else if (cur_al_len != 1) {
                any_non_snv = true;
            }
            f_start = f_end + 1;
            f_end = std::min(al_len.find(',', f_start), al_len.find('/', f_start));
        }
        // process the last element
        cur_al_len = std::stoi(al_len.substr(f_start, al_len.size()).c_str());
        if (cur_al_len >= 50) {
            any_sv = true;
        } else if (cur_al_len != 1) {
            any_non_snv = true;
        }
        // assign the variant type for this snarl
        if (any_sv) {
            var_types[idx] = "SV";
        } else if(any_non_snv) {
            var_types[idx] = "MNP";
        } else {
            var_types[idx] = "SNP";
        }
    }
    
    return var_types;
}
