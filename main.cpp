#include "swr.h"


//usage // <executable name> <input file containing both s1 and s2> <0: global, 1: local> <optional: path to parameters config file>
int main(int argc, char *argv[]){
    //parse input parameters 
    string a_sequence; 
    string b_sequence;
    ifstream seq_file (argv[1]);
    ifstream params (argv[3]);
    int match = 1;
    int mismatch = -2;
    int h = -5;
    int g = -1;
    //gets sequences from a sequence file or defaults 
    if (seq_file.is_open()){
        getline (seq_file,a_sequence);
        getline (seq_file,b_sequence);
    }
    else  {
        a_sequence = "acatgctacacgtactccgataccccgtaaccgataacgatacacagacctcgtacgcttgctacaacgtactctataaccgagaacgattgacatgcctcgtacacatgctacacgtactccgatgaccccgt";
        b_sequence = "acattctacgaacctctcgataaccccataaccgataacgattgacacctcgtacgctttctacaacttactctctcgataaccccataaccgataacgattgacacctcgtacacatggtacatacgtactctcgataccccgt";
    }
    cout << "****************Params*Start***************" << endl;
    //print sequences 
    cout << a_sequence << endl;
    cout << b_sequence << endl;
    //reads parameter.config file
    if (params.is_open()){
        string line;

        //reads line to value while converting string to int 
        //I should clean this up eventually 
        while(getline(params,line)){
            string field;
            string value;
            bool flip = 0;
            int converted_value = 0;
            for(char c : line){
                if(c == ' '){
                    flip = 1;
                }
                if(flip == 0){
                    field += c;
                }
                else{
                    value += c;
                }
            }

            cout << field << value << endl;
            std::istringstream(value) >> converted_value;
            if(field.compare("match") == 0){
                match = converted_value;
            }
            else if(field.compare("mismatch") == 0){
                mismatch = converted_value;
            }
            else if(field.compare("h") == 0){
                h = converted_value;
            }
            else if(field.compare("g") == 0){
                g = converted_value;
            }

        }
    }
    cout << "****************End*Params*****************" << endl << endl;

    cout << "*computing swr grid*" << endl;
    swr *comp = new swr(a_sequence, b_sequence, match, mismatch);
    
    comp->compute_grid();
    
    //comp->print_grid();
    point result = comp->largest_alignment_index(comp->determine_largest());

    cout << endl << "****************Results********************" << endl;
    comp->print_local_alignment(result);
    
    return 0;
}