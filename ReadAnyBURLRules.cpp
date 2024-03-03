#include "effectiveness.cpp"
#include <chrono>
#include <thread>

std::vector<std::string> split(const std::string& str, const std::string &pattern, bool app){
    std::vector<std::string> res;
    if(str.empty()) return res;

    std::string::size_type temp_pos;
    std::string appstr;
    if(app) appstr = str + pattern;
    else appstr = str;

    int size = appstr.size();
    for(unsigned int i=0;i<size;i++){
        temp_pos = appstr.find(pattern, i);
        if(temp_pos < size){
            std::string substr = appstr.substr(i, temp_pos - i);
            res.push_back(substr);
            i = temp_pos;
        }
    }
    return res;
}

void Effective_prop_DSD(const std::string &choice, const std::string &metric){
    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    unsigned int rule_count = 0;
    unsigned int dom_count = 0;

    getline(rules_in, rules_line);
    rules_in.close();

    double avg_time = 0.0;

    int state = -1;
    // state == 0: next int represents a variable rule;
    // state == 1: next int represents a instance rule;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                qp->EDirect.pop_back();
                qp->ETypes.pop_back();
                qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }

                    peeling::COD_prop(qp, metric, &g, visited, back_visited, dom_count);

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    peeling::COD_prop(qp, metric, &g, visited, back_visited, dom_count);
                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    qp->clear();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;
}

void Effective_hg_DSD(const std::string &choice, const std::string &metric){
    HeterGraph g(choice);

    std::string rules_path = (choice + "/"+choice+"-cod-global-rules.dat");
    std::ifstream rules_in;
    rules_in.open(rules_path);
    std::string rules_line;

    auto visited = new std::vector<std::vector<bool>*>();
    auto back_visited = new std::vector<std::vector<bool>*>();
    auto qp = new Pattern();

    getline(rules_in, rules_line);
    rules_in.close();

    unsigned int rule_count = 0;
    double avg_construct_time = 0.0;
    double avg_peel_time = 0.0;
    //double avg_goodness_density = 0.0;

    int state = -1;
    // state == 0: next int represents a variable rule;
    // state == 1: next int represents a instance rule;

    std::string::size_type temp_pos;
    int size = rules_line.size();

    for(unsigned int i = 0; i<size;i++) {
        temp_pos = rules_line.find(' ', i);
        if (temp_pos < size) {
            int sub = stoi(rules_line.substr(i, temp_pos - i));

            if (sub == -1) state = 0; // next int indicates a variable rule
            else if (sub == -2) qp->EDirect.push_back(1); // ->
            else if (sub == -3) qp->EDirect.push_back(-1); // <-
            else if (sub == -4) { // pop
                qp->EDirect.pop_back();
                qp->ETypes.pop_back();
                qp->NTypes.pop_back();
            } else if (sub == -5) state = 1; // next int indicates a instance rule
            else {
                if (state == 1) { // instance rule
                    qp->instance = sub;
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    auto valid = peeling::COD_hg(qp, metric, &g, visited, back_visited, avg_construct_time, avg_peel_time);
                    if(valid) rule_count++;

                    qp->instance = -1;
                    state = -1;
                } else if (state == 0) { // variable rule
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                    while(visited->size()< qp->ETypes.size()+1){
                        visited->push_back(new std::vector<bool>(g.NT.size(), false));
                        back_visited->push_back(new std::vector<bool>(g.NT.size(), false));
                    }
                    auto valid = peeling::COD_hg(qp, metric, &g, visited, back_visited, avg_construct_time, avg_peel_time);
                    if(valid) rule_count++;

                    state = -1;
                } else {
                    qp->ETypes.push_back(sub);
                    qp->NTypes.push_back(-1);
                }
            }
            i = temp_pos;
        }
    }
    qp->clear();

    for (auto &it : *visited) delete it;
    for (auto &it : *back_visited) delete it;
    delete visited;
    delete back_visited;

    avg_construct_time /= rule_count;
    avg_peel_time /= rule_count;
    //avg_goodness_density /= rule_count;

    std::cout<<"~time_construct_per_rule:"<<avg_construct_time<<" s"<<std::endl;
    std::cout<<"~time_peel_per_rule:"<<avg_peel_time<<" s"<<std::endl;
    //std::cout<<"~(peel/peelpp)_density_ratio:"<<avg_goodness_density<<std::endl;
    std::cout<<"rule_count:"<<rule_count<<std::endl;
}