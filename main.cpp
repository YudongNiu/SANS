#include "ReadAnyBURLRules.cpp"

using namespace std;

int main(int argc, char* argv[]){
    /** cod */
    std::string dataset(argv[1]);
    std::string method(argv[2]);
    
    if(method == "peele") Effective_hg_DSD(dataset, "edge");
    else if(method == "peelt") Effective_hg_DSD(dataset, "triangle");
    else if(method == "scanse") Effective_prop_DSD(dataset, "edge");
    else if(method == "scanst") Effective_prop_DSD(dataset, "triangle");
}