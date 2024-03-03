#include <iostream>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <map>
#include <chrono>
#include "hin.cpp"
#include "pattern.cpp"
#include "peer.cpp"
#include "synopse.cpp"
#include "join_synopse.cpp"
#include "param.h"
#include "hg.cpp"
#include "hidden_edge.cpp"
#include "compare.cpp"
#include "updatable_priority_queue.h"
#include "flowgraph.cpp"
#include "maxflow.cpp"


namespace peeling{
    double density(std::vector<std::vector<unsigned int> *> *g, std::vector<unsigned int> *s){
        std::vector<bool> indication(g->size(), false);
        for(unsigned int i: *s) indication[i] = true;
        double senum=0.0;
        for(unsigned int n: *s) for (unsigned int i : *g->at(n)) if (indication[i]) senum += 1.0;
        senum /= s->size();
        return senum / 2.0;
    }
    
    std::pair<std::vector<std::unordered_map<unsigned int, unsigned int>*>*, std::vector<unsigned int>*>
            triangle_hg_construction(std::vector<std::vector<unsigned int>*>* hg){
        auto triangle_deg = new std::vector<unsigned int>(hg->size(), 0);
        auto triangle_hg = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        for (auto &i : *hg) {
            if (!i->empty()) triangle_hg->push_back(new std::unordered_map<unsigned int, unsigned int>());
            else triangle_hg->push_back(nullptr);
        }

        for(unsigned int i=0;i<hg->size();i++){
            if(!hg->at(i)->empty()){
                for(unsigned int j=0;j<hg->at(i)->size();j++){
                    unsigned int nbr = hg->at(i)->at(j);
                    if(nbr >= i) break;
                    for(unsigned int jj=0;jj<hg->at(nbr)->size();jj++){
                        unsigned int nnbr = hg->at(nbr)->at(jj);
                        if(nnbr >= nbr) break;
                        bool triangle_find = false;
                        int left = 0, right = hg->at(i)->size() - 1;
                        while(left <= right){
                            int mid = left + (right - left) / 2;
                            if(hg->at(i)->at(mid) == nnbr) {
                                triangle_find = true; break;
                            }
                            if(hg->at(i)->at(mid) < nnbr) left = mid + 1;
                            else right = mid - 1;
                        }

                        if(triangle_find){
                            if(triangle_hg->at(i)->find(nbr) == triangle_hg->at(i)->end())
                                triangle_hg->at(i)->insert({nbr, 1});
                            else triangle_hg->at(i)->at(nbr)++;

                            if(triangle_hg->at(i)->find(nnbr) == triangle_hg->at(i)->end())
                                triangle_hg->at(i)->insert({nnbr, 1});
                            else triangle_hg->at(i)->at(nnbr) = (triangle_hg->at(i)->at(nnbr) + 1);

                            if(triangle_hg->at(nbr)->find(i) == triangle_hg->at(nbr)->end())
                                triangle_hg->at(nbr)->insert({i, 1});
                            else triangle_hg->at(nbr)->at(i)++;

                            if(triangle_hg->at(nbr)->find(nnbr) == triangle_hg->at(nbr)->end())
                                triangle_hg->at(nbr)->insert({nnbr, 1});
                            else triangle_hg->at(nbr)->at(nnbr) = (triangle_hg->at(nbr)->at(nnbr) + 1);

                            if(triangle_hg->at(nnbr)->find(i) == triangle_hg->at(nnbr)->end())
                                triangle_hg->at(nnbr)->insert({i, 1});
                            else triangle_hg->at(nnbr)->at(i)++;

                            if(triangle_hg->at(nnbr)->find(nbr) == triangle_hg->at(nnbr)->end())
                                triangle_hg->at(nnbr)->insert({nbr, 1});
                            else triangle_hg->at(nnbr)->at(nbr) = (triangle_hg->at(nnbr)->at(nbr) + 1);

                            triangle_deg->at(i)++;
                            triangle_deg->at(nbr)++;
                            triangle_deg->at(nnbr)++;
                        }}}}}
        return {triangle_hg, triangle_deg};
    }
    
    double triangle_density(std::vector<std::vector<unsigned int>*>* g, std::vector<unsigned int>*s){
        double triangle_count = 0;
        std::vector<bool> indication(g->size(), false);
        for(unsigned int i: *s) indication[i] = true;

        auto induced_g = new std::vector<std::vector<unsigned int>*>();
        for(unsigned int i=0;i<g->size();i++) induced_g->push_back(new std::vector<unsigned int>());
        for(unsigned int i=0;i<g->size();i++){
            if(indication[i])
            for(unsigned int j=0;j<g->at(i)->size();j++){
                unsigned int n = g->at(i)->at(j);
                if(indication[n] && n != i) induced_g->at(i)->push_back(n);
            }
        }
        
        auto triangle_hg_pair = triangle_hg_construction(induced_g);
        auto temp = triangle_hg_pair.first;
        auto triangle_deg = triangle_hg_pair.second;
        for(auto n: *s) {
            triangle_count += triangle_deg->at(n);
        }

        for(auto &i: *induced_g) delete i;
        delete induced_g;
        delete triangle_deg;
        return triangle_count / (3.0 * s->size());
    }
    
    std::pair<double, double> hg_peel(std::vector<std::vector<unsigned int>*>* hg){
        better_priority_queue::updatable_priority_queue<int, int> pQ;
        std::vector<unsigned int> remained_deg;
        double current_enum = 0.0;

        unsigned int hg_size = 0;

        for(unsigned int i=0;i<hg->size();i++){
            if(!hg->at(i)->empty()){
                pQ.push(i, -hg->at(i)->size());
                remained_deg.push_back(hg->at(i)->size());
                current_enum += hg->at(i)->size();
                hg_size++;
            }
            else remained_deg.push_back(0);
        }
        
        std::vector<unsigned int> rm_nodes;
        double peel_den = current_enum / hg_size;
        double peel_size = hg_size;

        while(pQ.size()>1){
            int next = pQ.pop_value().key;
            current_enum -= (2 * remained_deg[next] - 1);

            rm_nodes.push_back((unsigned int)next);
            double current_den = current_enum / (hg_size-rm_nodes.size());
            if(current_den > peel_den) {
                peel_den = current_den;
                peel_size = hg_size-rm_nodes.size();
            }

            remained_deg[next] = 0;
            for (int nbr : *hg->at((unsigned int) next)) {
                if(remained_deg[nbr] > 0){
                    remained_deg[nbr] -= 1;
                    pQ.update(nbr, -remained_deg[nbr]);
                }
            }
        }
        return {peel_den / 2.0 , peel_size};
    }
    
    double parresketch_peel(std::vector<std::vector<unsigned int>*>* g){
        std::random_device dev;
        std::mt19937 generator(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);

        unsigned int qualified_size = g->size();
        std::vector<bool> qualified_nodes(g->size(), true);

        auto ractive = new std::vector<std::vector<unsigned int>*>();
        ractive->push_back(new std::vector<unsigned int>());
        for(unsigned int i=0;i<g->size();i++) ractive->at(0)->push_back(i + 1);

        auto com = new std::vector<unsigned int>();

        better_priority_queue::updatable_priority_queue<int, double> pQ;
        bool resketch = true;

        auto remained_deg = new std::unordered_map<unsigned int, double>();
        double current_enum = 0.0;
        int removed_index = -1;
        std::vector<unsigned int> rm_nodes;
        double peel_den = 0.0;

        auto synopses = new std::vector<jsy::Synopse>();
        for(unsigned int i=0;i<g->size();i++) {
            jsy::Synopse s; jsy::init_synopse(&s, i);
            synopses->push_back(s);
        }

        auto reverse_synopses = new std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>();
        for(unsigned int i=0;i<L;i++) reverse_synopses->push_back(new std::unordered_map<unsigned int, std::vector<unsigned int>*>());
        auto rand2ps = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        auto p2rands = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        for(unsigned int l=0;l<L;l++){
            rand2ps->push_back(new std::unordered_map<unsigned int, unsigned int>());
            p2rands->push_back(new std::unordered_map<unsigned int, unsigned int>());
        }

        for(unsigned int l=0;l<L;l++){
            for(unsigned int i=0;i<g->size();i++){
                int rand = distribute(generator);
                while(rand2ps->at(l)->count(rand) > 0) rand = distribute(generator);
                rand2ps->at(l)->insert({rand, i});
                p2rands->at(l)->insert({i, rand});
            }
        }
        
        while(qualified_size > 1){
            if(resketch){
                current_enum = 0.0;
                for(unsigned int i=0;i<L;i++){
                    for (auto &j : *reverse_synopses->at(i)) delete j.second;
                    reverse_synopses->at(i)->clear();
                }
                for (auto &j : *synopses) clear(&j, true);
                remained_deg->clear();
                
                for(unsigned int l = 0; l < L; l++)
                    for(unsigned int n=0;n<g->size();n++)
                        if(qualified_nodes[n])
                            for(unsigned int nn=0;nn<g->at(n)->size();nn++)
                                if(qualified_nodes[g->at(n)->at(nn)]){
                                    jsy::add(&synopses->at(n), l, p2rands->at(l)->at(g->at(n)->at(nn)));
                                }
                for(unsigned int l=0;l < L;l++) for(unsigned int n=0;n<g->size();n++) synopses->at(n).rs[l] = p2rands->at(l)->at(n);
                jsy::set_overflow(synopses);

                for(unsigned int i=0;i<g->size();i++){
                    auto synopse = synopses->at(i);
                    for(unsigned int l=0;l<L;l++)
                        for(unsigned int j=0;j<synopse.sizes[l];j++) {
                            if(reverse_synopses->at(l)->find(rand2ps->at(l)->at(synopse.minks[l*K+j])) == reverse_synopses->at(l)->end())
                                reverse_synopses->at(l)->insert({rand2ps->at(l)->at(synopse.minks[l*K+j]), new std::vector<unsigned int>()});
                            reverse_synopses->at(l)->at(rand2ps->at(l)->at(synopse.minks[l*K+j]))->push_back(i);}
                }

                for(unsigned int i=0;i<g->size();i++)
                    if(qualified_nodes[i]) remained_deg->insert({i, jsy::dyn_estimate(&synopses->at(i))});

                for(auto &iter: *remained_deg){
                    bool success = pQ.push(iter.first, -iter.second);
                    if(!success) pQ.update(iter.first, -iter.second);
                    current_enum += iter.second;
                }
                double current_den = current_enum / qualified_size;
                if(current_den > peel_den){
                    peel_den = current_den;
                    removed_index = g->size() - qualified_size - 1;
                }
            }
            auto next = pQ.pop_value();
            qualified_size--;
            rm_nodes.push_back((unsigned int)next.key);
            current_enum += next.priority;

            qualified_nodes[(unsigned int)next.key] = false;

            auto updated_nodes = new std::unordered_map<unsigned int, bool>();
            resketch = jsy::remove((unsigned int)(next.key), synopses, ractive, &qualified_nodes, reverse_synopses, updated_nodes);

            if(!resketch){
                for (auto &updated_node : *updated_nodes) {
                    auto updated = updated_node.first;
                    auto priority_pair = pQ.get_priority(updated);
                    if(priority_pair.first){
                        double pre_est = -priority_pair.second;
                        double temp_den = jsy::dyn_estimate(&synopses->at(updated));
                        current_enum -= (pre_est - temp_den);
                        pQ.update(updated, -temp_den);
                    }
                }
                remained_deg->at((unsigned int)next.key) = 0;
                double current_den = current_enum / qualified_size;
                if(current_den > peel_den){
                    peel_den = current_den;
                    removed_index = g->size() - qualified_size - 1;
                }
            }
        }

        std::vector<bool> preserved(ractive->at(0)->size(), false);
        for(unsigned int i=0;i<preserved.size();i++) if(ractive->at(0)->at(i) > 0) preserved[i] = true;
        for(int i=0;i<=removed_index;i++) preserved[rm_nodes[i]] = false;
        for(unsigned int i=0;i<preserved.size();i++) if(preserved[i]) com->push_back(i);

        auto den = density(g, com);

        auto hg_pair = hg_peel(g);

        return den / hg_pair.first;
    }
    
    double sketch_peel(std::vector<std::vector<unsigned int>*>* g){
        std::random_device dev;
        std::mt19937 generator(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);

        unsigned int qualified_size = g->size();
        std::vector<bool> qualified_nodes(g->size(), true);

        auto ractive = new std::vector<std::vector<unsigned int>*>();
        ractive->push_back(new std::vector<unsigned int>());
        for(unsigned int i=0;i<g->size();i++) ractive->at(0)->push_back(i + 1);

        auto com = new std::vector<unsigned int>();

        better_priority_queue::updatable_priority_queue<int, double> pQ;
        bool resketch = true;

        auto remained_deg = new std::unordered_map<unsigned int, double>();
        double current_enum = 0.0;
        int removed_index = -1;
        std::vector<unsigned int> rm_nodes;
        double peel_den = 0.0;

        auto synopses = new std::vector<jsy::Synopse>();
        for(unsigned int i=0;i<g->size();i++) {
            jsy::Synopse s; jsy::init_synopse(&s, i);
            synopses->push_back(s);
        }

        auto reverse_synopses = new std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>();
        for(unsigned int i=0;i<L;i++) reverse_synopses->push_back(new std::unordered_map<unsigned int, std::vector<unsigned int>*>());
        auto rand2ps = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        auto p2rands = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        for(unsigned int l=0;l<L;l++){
            rand2ps->push_back(new std::unordered_map<unsigned int, unsigned int>());
            p2rands->push_back(new std::unordered_map<unsigned int, unsigned int>());
        }

        while(qualified_size > 1){
            if(resketch){
                current_enum = 0.0;
                for(unsigned int i=0;i<L;i++){
                    rand2ps->at(i)->clear();
                    p2rands->at(i)->clear();
                    for (auto &j : *reverse_synopses->at(i)) delete j.second;
                    reverse_synopses->at(i)->clear();
                }
                for (auto &j : *synopses) clear(&j, true);
                remained_deg->clear();

                for(unsigned int l=0;l<L;l++){
                    for(unsigned int i=0;i<g->size();i++){
                        int rand = distribute(generator);
                        while(rand2ps->at(l)->count(rand) > 0) rand = distribute(generator);
                        rand2ps->at(l)->insert({rand, i});
                        p2rands->at(l)->insert({i, rand});
                    }
                }

                for(unsigned int l = 0; l < L; l++)
                    for(unsigned int n=0;n<g->size();n++) 
                        if(qualified_nodes[n])
                            for(unsigned int nn=0;nn<g->at(n)->size();nn++) 
                                if(qualified_nodes[g->at(n)->at(nn)]){
                                    jsy::add(&synopses->at(n), l, p2rands->at(l)->at(g->at(n)->at(nn)));
                                }
                for(unsigned int l=0;l < L;l++) for(unsigned int n=0;n<g->size();n++) synopses->at(n).rs[l] = p2rands->at(l)->at(n);
                jsy::set_overflow(synopses);

                for(unsigned int i=0;i<g->size();i++){
                    auto synopse = synopses->at(i);
                    for(unsigned int l=0;l<L;l++)
                        for(unsigned int j=0;j<synopse.sizes[l];j++) {
                            if(reverse_synopses->at(l)->find(rand2ps->at(l)->at(synopse.minks[l*K+j])) == reverse_synopses->at(l)->end())
                                reverse_synopses->at(l)->insert({rand2ps->at(l)->at(synopse.minks[l*K+j]), new std::vector<unsigned int>()});
                            reverse_synopses->at(l)->at(rand2ps->at(l)->at(synopse.minks[l*K+j]))->push_back(i);}
                }

                for(unsigned int i=0;i<g->size();i++)
                    if(qualified_nodes[i]) remained_deg->insert({i, jsy::dyn_estimate(&synopses->at(i))});
                
                for(auto &iter: *remained_deg){
                    bool success = pQ.push(iter.first, -iter.second);
                    if(!success) pQ.update(iter.first, -iter.second);
                    current_enum += iter.second;
                }
                double current_den = current_enum / qualified_size;
                if(current_den > peel_den){
                    peel_den = current_den;
                    removed_index = g->size() - qualified_size - 1;
                }
            }
            auto next = pQ.pop_value();
            qualified_size--;
            rm_nodes.push_back((unsigned int)next.key);
            current_enum += next.priority;

            qualified_nodes[(unsigned int)next.key] = false;

            auto updated_nodes = new std::unordered_map<unsigned int, bool>();
            resketch = jsy::remove((unsigned int)(next.key), synopses, ractive, &qualified_nodes, reverse_synopses, updated_nodes);

            if(!resketch){
                for (auto &updated_node : *updated_nodes) {
                    auto updated = updated_node.first;
                    auto priority_pair = pQ.get_priority(updated);
                    if(priority_pair.first){
                        double pre_est = -priority_pair.second;
                        double temp_den = jsy::dyn_estimate(&synopses->at(updated));
                        current_enum -= (pre_est - temp_den);
                        pQ.update(updated, -temp_den);
                    }
                }
                remained_deg->at((unsigned int)next.key) = 0;
                double current_den = current_enum / qualified_size;
                if(current_den > peel_den){
                    peel_den = current_den;
                    removed_index = g->size() - qualified_size - 1;
                }
            }
        }

        std::vector<bool> preserved(ractive->at(0)->size(), false);
        for(unsigned int i=0;i<preserved.size();i++) if(ractive->at(0)->at(i) > 0) preserved[i] = true;
        for(int i=0;i<=removed_index;i++) preserved[rm_nodes[i]] = false;
        for(unsigned int i=0;i<preserved.size();i++) if(preserved[i]) com->push_back(i);

        auto den = density(g, com);

        auto hg_pair = hg_peel(g);
        
        return den / (hg_pair.first);
    }

    double nonresketch_peel(std::vector<std::vector<unsigned int>*>* g){
        std::random_device dev;
        std::mt19937 generator(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);

        auto reverse_synopses = new std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>();
        for(unsigned int l=0;l<L;l++)
            reverse_synopses->push_back(new std::unordered_map<unsigned int, std::vector<unsigned int>*>());

        auto rand2ps = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        auto p2rands = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        for(unsigned int l=0;l<L;l++){
            rand2ps->push_back(new std::unordered_map<unsigned int, unsigned int>());
            p2rands->push_back(new std::unordered_map<unsigned int, unsigned int>());
        }
        for(unsigned int l=0;l<L;l++) {
            for (unsigned int i = 0; i < g->size(); i++) {
                int rand = distribute(generator);
                while (rand2ps->at(l)->count(rand) > 0) rand = distribute(generator);
                rand2ps->at(l)->insert({rand, i});
                p2rands->at(l)->insert({i, rand});
            }
        }
        auto synopses = new std::vector<jsy::Synopse>();
        for(unsigned int i=0;i<g->size();i++) {
            jsy::Synopse s; jsy::init_synopse(&s, i);
            synopses->push_back(s);
        }
        for(unsigned int l = 0; l < L; l++)
            for(unsigned int n=0;n<g->size();n++)
                for(unsigned int nn=0;nn<g->at(n)->size();nn++)
                    jsy::add(&synopses->at(n), l, p2rands->at(l)->at(g->at(n)->at(nn)));
        for(unsigned int l=0;l < L;l++) for(unsigned int n=0;n<g->size();n++) synopses->at(n).rs[l] = p2rands->at(l)->at(n);
        jsy::set_overflow(synopses);

        for(unsigned int i=0;i<g->size();i++){
            auto synopse = synopses->at(i);
            for(unsigned int l=0;l<L;l++)
                for(unsigned int j=0;j<synopse.sizes[l];j++) {
                    if(reverse_synopses->at(l)->find(rand2ps->at(l)->at(synopse.minks[l*K+j])) == reverse_synopses->at(l)->end())
                        reverse_synopses->at(l)->insert({rand2ps->at(l)->at(synopse.minks[l*K+j]), new std::vector<unsigned int>()});
                    reverse_synopses->at(l)->at(rand2ps->at(l)->at(synopse.minks[l*K+j]))->push_back(i);
                }
        }

        auto com = new std::vector<unsigned int>();
        better_priority_queue::updatable_priority_queue<int, double> pQ;
        std::vector<bool> qualified_nodes(g->size(), true);
        unsigned int qualified_size = g->size();

        auto remained_deg = new std::unordered_map<unsigned int, double>();
        double current_enum = 0.0;

        for(unsigned int i=0;i<synopses->size();i++)
            remained_deg->insert({i, jsy::dyn_estimate(&synopses->at(i))});
        for(auto &iter: *remained_deg){
            bool success = pQ.push(iter.first, -iter.second);
            if(!success) pQ.update(iter.first, -iter.second);
            current_enum += iter.second;
        }

        double peel_den = current_enum / g->size();
        int removed_index = -1;
        std::vector<double> dens;
        std::vector<unsigned int> rm_nodes;

        auto ractive = new std::vector<std::vector<unsigned int>*>();
        ractive->push_back(new std::vector<unsigned int>());
        for(unsigned int i=0;i<g->size();i++) ractive->at(0)->push_back(i + 1);

        while(qualified_size > 1){
            auto next = pQ.pop_value();
            qualified_nodes[(unsigned int)next.key] = false;

            current_enum += next.priority;
            auto updated_nodes = new std::unordered_map<unsigned int, bool>();
            jsy::remove_nonresketch((unsigned int)(next.key), synopses, ractive, &qualified_nodes, reverse_synopses, updated_nodes);

            for (auto &updated_node : *updated_nodes) {
                auto updated = updated_node.first;
                auto priority_pair = pQ.get_priority(updated);
                if(priority_pair.first){
                    double pre_est = -priority_pair.second;
                    double temp_den = jsy::dyn_estimate(&synopses->at(updated));
                    current_enum -= (pre_est - temp_den);
                    pQ.update(updated, -temp_den);
                }
            }

            remained_deg->at((unsigned int)next.key) = 0;

            qualified_size--;
            rm_nodes.push_back((unsigned int)next.key);

            double current_den = current_enum / qualified_size;
            if(current_den > peel_den){
                peel_den = current_den;
                removed_index = g->size() - qualified_size - 1;
            }
        }

        std::vector<bool> preserved(g->size(), true);
        for(int i=0;i<=removed_index;i++) preserved[rm_nodes[i]] = false;
        for(unsigned int i=0;i<preserved.size();i++) if(preserved[i]) com->push_back(i);
        
        auto den = density(g, com);

        auto hg_pair = hg_peel(g);
        return den / (hg_pair.first);
    }

    unsigned int resketch_analysis(std::vector<std::vector<unsigned int>*>* g, const std::string & metric){
        std::random_device dev;
        std::mt19937 generator(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);

        auto reverse_synopses = new std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>();
        for(unsigned int l=0;l<L;l++)
            reverse_synopses->push_back(new std::unordered_map<unsigned int, std::vector<unsigned int>*>());

        auto rand2ps = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        auto p2rands = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        for(unsigned int l=0;l<L;l++){
            rand2ps->push_back(new std::unordered_map<unsigned int, unsigned int>());
            p2rands->push_back(new std::unordered_map<unsigned int, unsigned int>());
        }
        for(unsigned int l=0;l<L;l++) {
            for (unsigned int i = 0; i < g->size(); i++) {
                int rand = distribute(generator);
                while (rand2ps->at(l)->count(rand) > 0) rand = distribute(generator);
                rand2ps->at(l)->insert({rand, i});
                p2rands->at(l)->insert({i, rand});
            }
        }
        auto synopses = new std::vector<jsy::Synopse>();
        for(unsigned int i=0;i<g->size();i++) {
            jsy::Synopse s; jsy::init_synopse(&s, i);
            synopses->push_back(s);
        }
        for(unsigned int l = 0; l < L; l++)
            for(unsigned int n=0;n<g->size();n++)
                for(unsigned int nn=0;nn<g->at(n)->size();nn++)
                    jsy::add(&synopses->at(n), l, p2rands->at(l)->at(g->at(n)->at(nn)));
        jsy::set_overflow(synopses);

        for(unsigned int i=0;i<g->size();i++){
            auto synopse = synopses->at(i);
            for(unsigned int l=0;l<L;l++)
                for(unsigned int j=0;j<synopse.sizes[l];j++) {
                    if(reverse_synopses->at(l)->find(synopse.minks[l*K+j]) == reverse_synopses->at(l)->end())
                        reverse_synopses->at(l)->insert({synopse.minks[l*K+j], new std::vector<unsigned int>()});
                    reverse_synopses->at(l)->at(synopse.minks[l * K + j])->push_back(i);
                }
        }

        unsigned int resketch_num = 0;

        if(metric == "edge") {
            better_priority_queue::updatable_priority_queue<int, int> pQ;
            std::vector<unsigned int> remained_deg;

            unsigned int hg_size = 0;
            for (unsigned int i = 0; i < g->size(); i++) {
                if (!g->at(i)->empty()) {
                    pQ.push(i, -g->at(i)->size());
                    remained_deg.push_back(g->at(i)->size());
                    hg_size++;
                } else remained_deg.push_back(0);
            }

            std::vector<bool> removed(g->size(), false);

            while (pQ.size() > 1) {
                auto next = (unsigned int) pQ.pop_value().key;
                removed[next] = true;

                //remove next from sketches
                for (unsigned int l = 0; l < L; l++) {
                    auto to_be_resketch = jsy::remove(p2rands->at(l)->at(next), l, synopses, &removed, reverse_synopses);
                    if (to_be_resketch) {
                        resketch_num++;

                        for (auto &j: *reverse_synopses->at(l)) delete j.second;
                        reverse_synopses->at(l)->clear();
                        for (auto &i: *synopses) clear(&i, l);
                        rand2ps->at(l)->clear();
                        p2rands->at(l)->clear();

                        for (unsigned int i = 0; i < g->size(); i++) {
                            int rand = distribute(generator);
                            while (rand2ps->at(l)->count(rand) > 0) rand = distribute(generator);
                            rand2ps->at(l)->insert({rand, i});
                            p2rands->at(l)->insert({i, rand});
                        }
                        for (unsigned int n = 0; n < g->size(); n++)
                            for (unsigned int nn = 0; nn < g->at(n)->size(); nn++)
                                if (!removed[n] && !removed[g->at(n)->at(nn)])
                                    jsy::add(&synopses->at(n), l, p2rands->at(l)->at(g->at(n)->at(nn)));
                        jsy::set_overflow(synopses, l);

                        for (unsigned int i = 0; i < g->size(); i++) {
                            auto synopse = synopses->at(i);
                            if (!removed[i])
                                for (unsigned int j = 0; j < synopse.sizes[l]; j++) {
                                    if (reverse_synopses->at(l)->find(synopse.minks[l * K + j]) ==
                                        reverse_synopses->at(l)->end())
                                        reverse_synopses->at(l)->insert(
                                                {synopse.minks[l * K + j], new std::vector<unsigned int>()});
                                    reverse_synopses->at(l)->at(synopse.minks[l * K + j])->push_back(i);
                                }
                        }
                    }
                }

                remained_deg[next] = 0;
                for (int nbr: *g->at(next)) {
                    if (remained_deg[nbr] > 0) {
                        remained_deg[nbr] -= 1;
                        pQ.update(nbr, -remained_deg[nbr]);
                    }
                }
            }
        }
        else if (metric == "triangle"){
            better_priority_queue::updatable_priority_queue<int, int> pQ;
            unsigned int hg_size = 0;

            auto triangle_hg_pair = triangle_hg_construction(g);
            auto triangle_hg = triangle_hg_pair.first;
            auto remained_deg = triangle_hg_pair.second;

            for(unsigned int i=0;i<g->size();i++){
                if(!g->at(i)->empty()){
                    pQ.push(i, -remained_deg->at(i));
                    hg_size++;
                }
            }

            std::vector<bool> removed(g->size(), false);
            while(pQ.size() > 2){
                auto next = (unsigned int)pQ.pop_value().key;
                removed[next] = true;

                //remove next from sketches
                for (unsigned int l = 0; l < L; l++) {
                    auto to_be_resketch = jsy::remove(p2rands->at(l)->at(next), l, synopses, &removed, reverse_synopses);
                    if (to_be_resketch) {
                        resketch_num++;
                        for (auto &j: *reverse_synopses->at(l)) delete j.second;
                        reverse_synopses->at(l)->clear();
                        for (auto &i: *synopses) clear(&i, l);
                        rand2ps->at(l)->clear();
                        p2rands->at(l)->clear();

                        for (unsigned int i = 0; i < g->size(); i++) {
                            int rand = distribute(generator);
                            while (rand2ps->at(l)->count(rand) > 0) rand = distribute(generator);
                            rand2ps->at(l)->insert({rand, i});
                            p2rands->at(l)->insert({i, rand});
                        }
                        for (unsigned int n = 0; n < g->size(); n++)
                            for (unsigned int nn = 0; nn < g->at(n)->size(); nn++)
                                if (!removed[n] && !removed[g->at(n)->at(nn)])
                                    jsy::add(&synopses->at(n), l, p2rands->at(l)->at(g->at(n)->at(nn)));
                        jsy::set_overflow(synopses, l);

                        for (unsigned int i = 0; i < g->size(); i++) {
                            auto synopse = synopses->at(i);
                            if (!removed[i])
                                for (unsigned int j = 0; j < synopse.sizes[l]; j++) {
                                    if (reverse_synopses->at(l)->find(synopse.minks[l * K + j]) ==
                                        reverse_synopses->at(l)->end())
                                        reverse_synopses->at(l)->insert(
                                                {synopse.minks[l * K + j], new std::vector<unsigned int>()});
                                    reverse_synopses->at(l)->at(synopse.minks[l * K + j])->push_back(i);
                                }
                        }
                    }
                }

                remained_deg->at(next) = 0;
                for (auto &it : *triangle_hg->at(next)) {
                    auto tnbr = it.first;
                    if(remained_deg->at(tnbr) > 0){
                        remained_deg->at(tnbr) -= it.second;
                        pQ.update(tnbr, -remained_deg->at(tnbr));
                        for(auto &itt: *triangle_hg->at(tnbr)){
                            auto ttnbr = itt.first;
                            if(triangle_hg->at(next)->find(ttnbr) != triangle_hg->at(next)->end()){
                                triangle_hg->at(tnbr)->at(ttnbr) -= 1;
                            }
                        }
                    }
                }
                triangle_hg->at(next)->clear();
                delete triangle_hg->at(next);
                triangle_hg->at(next) = nullptr;
            }
        }
        return resketch_num;
    }

    std::vector<unsigned int>* FlowComputation(std::vector<std::vector<unsigned int>*>* g,
            double rho, int edge_num){
        auto source_seg = new std::vector<unsigned int>();
        auto fg = new FlowGraph<double, double, double>(g->size(), (int)edge_num);
        for(unsigned int i=0;i<g->size();i++) fg->add_node();
        for(unsigned int i=0;i<g->size();i++)
            fg->add_tweights(i, edge_num, edge_num + 2*rho - g->at(i)->size());
        for(unsigned int i=0;i<g->size();i++){
            for(unsigned int j=0;j<g->at(i)->size();j++){
                unsigned int n = g->at(i)->at(j);
                if(n > i){
                    fg->add_edge(i, n, 1, 1);
                }
            }
        }
        double flow = fg->maxflow();
        for(unsigned int i=0;i<g->size();i++){
            if(fg->what_segment(i) == FlowGraph<double, double, double>::SOURCE)
                source_seg->push_back(i);
        }
        delete fg;
        return source_seg;
    }

    std::pair<double, std::vector<unsigned int>*> Goldberg(std::vector<std::vector<unsigned int>*>* g){
        unsigned int edge_num = 0;
        for (auto &i : *g) edge_num += i->size();

        auto subgraph = new std::vector<unsigned int>();

        double min_deg = 0;
        double max_deg = 0;
        for (auto &i : *g) {
            if(i->size() > max_deg)
                max_deg = i->size() - 1;
        }

        while(max_deg - min_deg >= 0.00001){
            double least_density = (max_deg +  min_deg) / 2;
            std::vector<unsigned int>* source_segment = FlowComputation(g, least_density, edge_num);
            if(source_segment->empty()){
                max_deg = least_density;
            }
            else{
                subgraph->clear();
                min_deg = least_density;
                for(unsigned int i: *source_segment) subgraph->push_back(i);
            }
            delete source_segment;
        }
        return {max_deg, subgraph};
    }

    double triangle_estimate(std::vector<jsy::Synopse>* synopses,
                            std::vector<std::vector<unsigned int>*>* ractive,
                             std::vector<std::unordered_map<unsigned int, unsigned int>*>* rand2ps,
                            unsigned int v){
        double nbr_size = dyn_estimate(&synopses->at(ractive->at(0)->at(v)-1));
        double I_sum = 0.0;
        double sample_count = 0.0;
        for(unsigned int l=0;l<L;l++){
            for(unsigned int i=0;i < synopses->at(ractive->at(0)->at(v)-1).sizes[l];i++){
                unsigned int r = synopses->at(ractive->at(0)->at(v)-1).minks[i + K * l];
                unsigned int nbr_sample = rand2ps->at(l)->at(r);
                double I = jsy::intersect_estimate(&synopses->at(ractive->at(0)->at(v)-1), &synopses->at(ractive->at(0)->at(nbr_sample)-1), l);
                I_sum += (I - 2);
                sample_count += 1.0;
            }
        }
        double est = ((nbr_size - 1) * I_sum) / (sample_count * 2.0);
        return est;
    }

    std::pair<std::vector<unsigned int>*, double> sketch_triangle_peel(unsigned int meta_layer,
            std::vector<HiddenEdge>* hidden_edges, std::vector<std::vector<unsigned int> *>* ractive,
            std::vector<std::vector<jsy::Synopse> *>* synopses){
        
        double resketch_count = 0.0;

        unsigned int qualified_size = 0;
        std::vector<bool> qualified_nodes(ractive->at(0)->size(), false);
        for(unsigned int i=0;i<qualified_nodes.size();i++)
            if(ractive->at(0)->at(i) > 0) {
                qualified_nodes[i] = true;
                qualified_size++;
            }
        unsigned int hg_size = qualified_size;

        auto com = new std::vector<unsigned int>();
        better_priority_queue::updatable_priority_queue<int, double> pQ;
        bool resketch = true;

        auto remained_deg = new std::unordered_map<unsigned int, double>();
        double current_enum = 0.0;
        int removed_index = -1;
        std::vector<double> dens;
        std::vector<unsigned int> rm_nodes;
        double peel_den = 0.0;

        auto reverse_synopses = new std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>();
        for(unsigned int i=0;i<L;i++) reverse_synopses->push_back(new std::unordered_map<unsigned int, std::vector<unsigned int>*>());

        std::vector<std::unordered_map<unsigned int, unsigned int>*>* rand2ps = nullptr;
        while(qualified_size > 2){
            if(resketch){
                resketch_count += 1.0;
                current_enum = 0.0;
                remained_deg->clear();

                for(unsigned int i=0;i<L;i++){
                    for (auto &j : *reverse_synopses->at(i)) delete j.second;
                    reverse_synopses->at(i)->clear();
                }

                for (unsigned int l = 0; l <= meta_layer; l++) for (auto &i : *synopses->at(l)) clear(&i, true);
                rand2ps = jsy::gnn_synopses(meta_layer, hidden_edges, synopses, &qualified_nodes);
                jsy::set_overflow(synopses->at(0), &qualified_nodes, false);

                for (auto &i : *synopses->at(0)){
                    unsigned int _p = i.p;
                    if(qualified_nodes[_p]) for(unsigned int l=0;l<L;l++)
                            reverse_synopses->at(l)->insert({_p, new std::vector<unsigned int>()});
                }

                for(unsigned int i=0;i<synopses->at(0)->size();i++){
                    auto synopse = synopses->at(0)->at(i);
                    if(qualified_nodes[synopse.p])
                        for(unsigned int l=0;l < L; l++) for(unsigned int j=0; j< synopse.sizes[l];j++)
                                reverse_synopses->at(l)->at(rand2ps->at(l)->at(synopse.minks[l * K + j]))->push_back(i);
                }

                for(auto &i: *synopses->at(0)){
                    unsigned int _p = i.p;
                    if(qualified_nodes[_p]){
                        double est = triangle_estimate(synopses->at(0), ractive, rand2ps, _p);
                        remained_deg->insert({_p, est});
                    }
                }
                for(auto &iter: *remained_deg){
                    bool success = pQ.push(iter.first, -iter.second);
                    if(!success) pQ.update(iter.first, -iter.second);
                    current_enum += iter.second;
                }
                double current_den = current_enum / qualified_size;
                if(current_den > peel_den){
                    peel_den = current_den;
                    removed_index = hg_size - qualified_size - 1;
                }
            }
            auto next = pQ.pop_value();
            current_enum += next.priority;

            auto updated_nodes = new std::unordered_map<unsigned int, bool>();
            resketch = jsy::remove((unsigned int)(next.key), synopses->at(0), ractive, &qualified_nodes, reverse_synopses, updated_nodes);

            if(!resketch){
                for (auto &updated_node : *updated_nodes) {
                    auto updated = updated_node.first;
                    auto priority_pair = pQ.get_priority(updated);
                    if(priority_pair.first){
                        double pre_est = -priority_pair.second;
                        double temp_den = triangle_estimate(synopses->at(0), ractive, rand2ps, updated);
                        current_enum -= (pre_est - temp_den);
                        pQ.update(updated, -temp_den);
                    }
                }
                remained_deg->at((unsigned int)next.key) = 0;
                qualified_nodes[(unsigned int)next.key] = false;
                qualified_size--;
                rm_nodes.push_back((unsigned int)next.key);
                
                double current_den = current_enum / qualified_size;
                if(current_den > peel_den){
                    peel_den = current_den;
                    removed_index = hg_size - qualified_size - 1;
                }
            }
        }

        std::vector<bool> preserved(ractive->at(0)->size(), false);
        for(unsigned int i=0;i<preserved.size();i++) if(ractive->at(0)->at(i) > 0) preserved[i] = true;
        for(int i=0;i<=removed_index;i++) preserved[rm_nodes[i]] = false;
        for(unsigned int i=0;i<preserved.size();i++) if(preserved[i]) com->push_back(i);

        for (unsigned int l = 0; l <= meta_layer; l++) for (auto &i : *synopses->at(l)) clear(&i, true);
        return {com, resketch_count};
    }

    std::vector<unsigned int>* nonresketch_peel(unsigned int meta_layer,
                                                std::vector<HiddenEdge>* hidden_edges,
                                                std::vector<std::vector<unsigned int>*>* ractive,
                                                std::vector<std::vector<jsy::Synopse>*>* synopses){
        unsigned int qualified_size = 0;
        std::vector<bool> qualified_nodes(ractive->at(0)->size(), false);
        for(unsigned int i=0;i<qualified_nodes.size();i++)
            if(ractive->at(0)->at(i) > 0) {
                qualified_nodes[i] = true;
                qualified_size++;
            }
        unsigned int hg_size = qualified_size;

        auto reverse_synopses = new std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>();
        for(unsigned int i=0;i<L;i++) reverse_synopses->push_back(new std::unordered_map<unsigned int, std::vector<unsigned int>*>());

        auto rand2ps = jsy::gnn_synopses(meta_layer, hidden_edges, synopses, &qualified_nodes);
        for(auto &i: *synopses->at(0)){
            unsigned int _p = i.p;
            if(qualified_nodes[_p]) for(unsigned int l=0;l<L;l++)
                    reverse_synopses->at(l)->insert({_p, new std::vector<unsigned int>()});
        }
        for(unsigned int i=0;i<synopses->at(0)->size();i++){
            auto synopse = synopses->at(0)->at(i);
            if(qualified_nodes[synopse.p])
                for(unsigned int l=0;l < L; l++) for(unsigned int j=0; j< synopse.sizes[l];j++)
                        reverse_synopses->at(l)->at(rand2ps->at(l)->at(synopse.minks[l * K + j]))->push_back(i);
        }
        jsy::set_overflow(synopses->at(0), &qualified_nodes);

        auto com = new std::vector<unsigned int>();
        better_priority_queue::updatable_priority_queue<int, double> pQ;

        auto remained_deg = new std::unordered_map<unsigned int, double>();
        double current_enum = 0.0;

        for(auto &i: *synopses->at(0)){
            unsigned int _p = i.p;
            if(qualified_nodes[_p]){
                remained_deg->insert({_p, jsy::dyn_estimate(&i)});
            }
        }
        for(auto &iter: *remained_deg){
            bool success = pQ.push(iter.first, -iter.second);
            if(!success) pQ.update(iter.first, -iter.second);
            current_enum += iter.second;
        }

        double peel_den = current_enum / qualified_size;
        int removed_index = -1;
        std::vector<double> dens;
        std::vector<unsigned int> rm_nodes;

        while(qualified_size > 1){
            auto next = pQ.pop_value();
            qualified_nodes[(unsigned int)next.key] = false;

            current_enum += next.priority;
            auto updated_nodes = new std::unordered_map<unsigned int, bool>();
            jsy::remove_nonresketch((unsigned int)(next.key), synopses->at(0), ractive, &qualified_nodes, reverse_synopses, updated_nodes);

            for (auto &updated_node : *updated_nodes) {
                auto updated = updated_node.first;
                auto priority_pair = pQ.get_priority(updated);
                if(priority_pair.first){
                    double pre_est = -priority_pair.second;
                    double temp_den = jsy::dyn_estimate(&synopses->at(0)->at(ractive->at(0)->at(updated)-1));
                    current_enum -= (pre_est - temp_den);
                    pQ.update(updated, -temp_den);
                }
            }

            remained_deg->at((unsigned int)next.key) = 0;

            qualified_size--;
            rm_nodes.push_back((unsigned int)next.key);

            double current_den = current_enum / qualified_size;
            if(current_den > peel_den){
                peel_den = current_den;
                removed_index = hg_size - qualified_size - 1;
            }
        }

        std::vector<bool> preserved(ractive->at(0)->size(), false);
        for(unsigned int i=0;i<preserved.size();i++) if(ractive->at(0)->at(i) > 0) preserved[i] = true;
        for(int i=0;i<=removed_index;i++) preserved[rm_nodes[i]] = false;
        for(unsigned int i=0;i<preserved.size();i++) if(preserved[i]) com->push_back(i);

        for (unsigned int l = 0; l <= meta_layer; l++) for (auto &i : *synopses->at(l)) clear(&i, true);
        return com;
    }

    std::pair<std::vector<unsigned int>*, double> sketch_peel(unsigned int meta_layer,
            std::vector<HiddenEdge>* hidden_edges,
            std::vector<std::vector<unsigned int> *>* ractive,
            std::vector<std::vector<jsy::Synopse> *>* synopses){
        double resketch_count = 0.0;
        
        unsigned int qualified_size = 0;
        std::vector<bool> qualified_nodes(ractive->at(0)->size(), false);
        for(unsigned int i=0;i<qualified_nodes.size();i++)
            if(ractive->at(0)->at(i) > 0) {
                qualified_nodes[i] = true;
                qualified_size++;
            }
        unsigned int hg_size = qualified_size;

        auto com = new std::vector<unsigned int>();
        better_priority_queue::updatable_priority_queue<int, double> pQ;
        bool resketch = true;

        auto remained_deg = new std::unordered_map<unsigned int, double>();
        double current_enum = 0.0;
        int removed_index = -1;
        std::vector<double> dens;
        std::vector<unsigned int> rm_nodes;
        double peel_den = 0.0;

        auto reverse_synopses = new std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>();
        for(unsigned int i=0;i<L;i++) reverse_synopses->push_back(new std::unordered_map<unsigned int, std::vector<unsigned int>*>());

        while(qualified_size > 1){
            if(resketch){
                resketch_count += 1.0;
                current_enum = 0.0;
                remained_deg->clear();

                for(unsigned int i=0;i<L;i++){
                    for (auto &j : *reverse_synopses->at(i)) delete j.second;
                    reverse_synopses->at(i)->clear();
                }

                for (unsigned int l = 0; l <= meta_layer; l++) for (auto &i : *synopses->at(l)) clear(&i, true);
                auto rand2ps = jsy::gnn_synopses(meta_layer, hidden_edges, synopses, &qualified_nodes);
                jsy::set_overflow(synopses->at(0), &qualified_nodes, false);

                for (auto &i : *synopses->at(0)){
                    unsigned int _p = i.p;
                    if(qualified_nodes[_p]) for(unsigned int l=0;l<L;l++)
                            reverse_synopses->at(l)->insert({_p, new std::vector<unsigned int>()});
                }

                for(unsigned int i=0;i<synopses->at(0)->size();i++){
                    auto synopse = synopses->at(0)->at(i);
                    if(qualified_nodes[synopse.p])
                        for(unsigned int l=0;l < L; l++) for(unsigned int j=0; j< synopse.sizes[l];j++)
                            reverse_synopses->at(l)->at(rand2ps->at(l)->at(synopse.minks[l * K + j]))->push_back(i);
                }

                for(auto &i: *synopses->at(0)){
                    unsigned int _p = i.p;
                    if(qualified_nodes[_p]){
                        remained_deg->insert({_p, jsy::dyn_estimate(&i)});
                    }
                }
                for(auto &iter: *remained_deg){
                    bool success = pQ.push(iter.first, -iter.second);
                    if(!success) pQ.update(iter.first, -iter.second);
                    current_enum += iter.second;
                }
                double current_den = current_enum / qualified_size;
                if(current_den > peel_den){
                    peel_den = current_den;
                    removed_index = hg_size - qualified_size - 1;
                }
            }

            auto next = pQ.pop_value();
            current_enum += next.priority;

            auto updated_nodes = new std::unordered_map<unsigned int, bool>();
            resketch = jsy::remove((unsigned int)(next.key), synopses->at(0), ractive, &qualified_nodes, reverse_synopses, updated_nodes);
            if(!resketch){
                for (auto &updated_node : *updated_nodes) {
                    auto updated = updated_node.first;
                    auto priority_pair = pQ.get_priority(updated);
                    if(priority_pair.first){
                        double pre_est = -priority_pair.second;
                        double temp_den = jsy::dyn_estimate(&synopses->at(0)->at(ractive->at(0)->at(updated)-1));
                        current_enum -= (pre_est - temp_den);
                        pQ.update(updated, -temp_den);
                    }
                }
                remained_deg->at((unsigned int)next.key) = 0;
                qualified_nodes[(unsigned int)next.key] = false;
                qualified_size--;
                rm_nodes.push_back((unsigned int)next.key);

                double current_den = current_enum / qualified_size;
                if(current_den > peel_den){
                    peel_den = current_den;
                    removed_index = hg_size - qualified_size - 1;
                }
            }
        }

        std::vector<bool> preserved(ractive->at(0)->size(), false);
        for(unsigned int i=0;i<preserved.size();i++) if(ractive->at(0)->at(i) > 0) preserved[i] = true;
        for(int i=0;i<=removed_index;i++) preserved[rm_nodes[i]] = false;
        for(unsigned int i=0;i<preserved.size();i++) if(preserved[i]) com->push_back(i);
        for (unsigned int l = 0; l <= meta_layer; l++) for (auto &i : *synopses->at(l)) clear(&i, true);
        return {com, resketch_count};
    }

    std::pair<double, double> hg_triangle_peel(std::vector<std::vector<unsigned int>*>* hg){
        better_priority_queue::updatable_priority_queue<int, int> pQ;
        double current_tnum=0.0;
        unsigned int hg_size = 0;

        auto triangle_hg_pair = triangle_hg_construction(hg);
        auto triangle_hg = triangle_hg_pair.first;
        auto remained_deg = triangle_hg_pair.second;

        for(unsigned int i=0;i<hg->size();i++){
            if(!hg->at(i)->empty()){
                pQ.push(i, -remained_deg->at(i));
                current_tnum += remained_deg->at(i);
                hg_size++;
            }
        }

        std::vector<unsigned int> rm_nodes;
        std::vector<double> dens;
        double peel_den = current_tnum / hg_size;
        double peel_size = hg_size;

        while(pQ.size() > 2){
            auto next = (unsigned int)pQ.pop_value().key;
            
            current_tnum -= 3 * remained_deg->at(next);
            rm_nodes.push_back(next);
            double current_den = current_tnum / (hg_size-rm_nodes.size());
            if(current_den > peel_den) {
                peel_den = current_den;
                peel_size = hg_size - rm_nodes.size();
            }

            remained_deg->at(next) = 0;
            for (auto &it : *triangle_hg->at(next)) {
                auto tnbr = it.first;
                if(remained_deg->at(tnbr) > 0){
                    remained_deg->at(tnbr) -= it.second;
                    pQ.update(tnbr, -remained_deg->at(tnbr));
                    for(auto &itt: *triangle_hg->at(tnbr)){
                        auto ttnbr = itt.first;
                        if(triangle_hg->at(next)->find(ttnbr) != triangle_hg->at(next)->end()){
                            triangle_hg->at(tnbr)->at(ttnbr) -= 1;
                        }
                    }
                }
            }
            triangle_hg->at(next)->clear();
            delete triangle_hg->at(next);
            triangle_hg->at(next) = nullptr;
        }
        return {peel_den / 3.0, peel_size};
    }

    double hg_peelpp(std::vector<std::vector<unsigned int>*>* hg, unsigned int T){
        std::vector<unsigned int> load(hg->size(), 0);
        double init_current_enum = 0.0;
        std::vector<unsigned int> remained_deg(hg->size(), 0);
        unsigned int hg_size = 0;
        for (auto &i : *hg) {
            if(!i->empty()){
                hg_size++;
                init_current_enum += i->size();
            }
        }

        double peel_den = init_current_enum / hg_size;

        for(unsigned int t = 0; t < T; t++){
            better_priority_queue::updatable_priority_queue<int, int> pQ;
            double current_enum = init_current_enum;
            for(unsigned int i=0;i<hg->size();i++){
                if(!hg->at(i)->empty()){
                    pQ.push(i, -hg->at(i)->size()-load[i]);
                    remained_deg[i] = hg->at(i)->size();
                }
                else remained_deg[i] = 0;
            }

            std::vector<bool> removed(hg->size(), false);
            std::vector<unsigned int> rm_nodes;
            while(pQ.size() > 1){
                int next = pQ.pop_value().key;
                current_enum -= (2 * remained_deg[next] - 1);
                load[next] += remained_deg[next];

                rm_nodes.push_back((unsigned int)next);
                double current_den = current_enum / (hg_size - rm_nodes.size());
                if(current_den > peel_den) peel_den = current_den;

                remained_deg[next] = 0;
                for(int nbr: *hg->at((unsigned int) next)){
                    if(remained_deg[nbr] > 0){
                        remained_deg[nbr] -= 1;
                        pQ.update(nbr, -remained_deg[nbr]-load[nbr]);
                    }
                }
            }
        }
        return peel_den;
    }

    void COD_casestudy(Pattern *qp, const std::string &method, HeterGraph *g, std::vector<std::vector<bool>*>* visited,
                       std::vector<std::vector<bool>*>* back_visited){
        double running_time = 0.0;
        double resketch_count = 0.0;

        auto start = std::chrono::steady_clock::now();

        auto frontiers = new std::vector<unsigned int>();
        auto peers = new std::vector<unsigned int>();
        Peers(qp, g, frontiers, peers, visited);

        std::vector<std::set<unsigned int> *> *active = ActiveMidNodes(peers, frontiers, qp, g, visited, back_visited);
        auto ractive = new std::vector<std::vector<unsigned int> *>();
        for (unsigned int l = 0; l <= qp->ETypes.size(); l++) ractive->push_back(new std::vector<unsigned int>(g->NT.size(), 0));
        for (unsigned int n = 0; n < g->NT.size(); n++) for (unsigned int l : *active->at(n)) ractive->at(l)->at(n) = 1;

        auto synopses_combined = new std::vector<std::vector<jsy::Synopse> *>();
        for (unsigned int l = 0; l <= qp->ETypes.size(); l++) {
            synopses_combined->push_back(new std::vector<jsy::Synopse>());
            for (unsigned int p = 0; p < g->NT.size(); p++)
                if (ractive->at(l)->at(p) > 0) {
                    jsy::Synopse s; jsy::init_synopse(&s, p);
                    synopses_combined->at(l)->push_back(s);
                    ractive->at(l)->at(p) = synopses_combined->at(l)->size();
                }
        }

        auto hidden_edges = new std::vector<HiddenEdge>();
        for (unsigned int l = 0; l < qp->ETypes.size(); l++) {
            for (unsigned int p = 0; p < g->NT.size(); p++) {
                if (ractive->at(l)->at(p) > 0) {
                    if (qp->EDirect[l] == 1) {
                        for (unsigned int nbr: *(g->EL[p])) {
                            if (ractive->at(l + 1)->at(nbr) > 0) {
                                HiddenEdge he{.s=ractive->at(l)->at(p) - 1, .t=ractive->at(l + 1)->at(nbr) - 1, .l=l};
                                hidden_edges->push_back(he);
                            }
                        }
                    } else {
                        for (unsigned int nbr: *(g->rEL[p])) {
                            if (ractive->at(l + 1)->at(nbr) > 0) {
                                HiddenEdge he{.s=ractive->at(l)->at(p) - 1, .t=ractive->at(l + 1)->at(nbr) - 1, .l=l};
                                hidden_edges->push_back(he);
                            }}}}}}
        unsigned int max_meta_layer;
        if(qp->instance == -1) max_meta_layer = qp->ETypes.size();
        else max_meta_layer = qp->ETypes.size() - 1;

        auto end = std::chrono::steady_clock::now();
        running_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        for (unsigned int meta_layer = 1; meta_layer <= max_meta_layer; meta_layer++) {
            if(method == "scanse") {
                auto sketch_peel_pair = sketch_peel(meta_layer, hidden_edges, ractive, synopses_combined);
                auto dense_subgraph = sketch_peel_pair.first;
                resketch_count += sketch_peel_pair.second;

                for(unsigned int detected: *dense_subgraph) std::cout<<detected<<" "; std::cout<<std::endl;

                std::vector<std::vector<unsigned int> *> *hg = hidden_graph_construction(meta_layer, qp, g, visited, ractive);
                std::vector<bool> reserve(hg->size(), false);
                for(unsigned int detected: *dense_subgraph) reserve[detected] = true;

                for(unsigned int i=0;i<hg->size();i++){
                    if(reserve[i]){
                        for(unsigned int j=0;j<hg->at(i)->size();j++){
                            if(reserve[hg->at(i)->at(j)] && hg->at(i)->at(j) > i){
                                std::cout<<i<<" "<<hg->at(i)->at(j)<<std::endl;
                            }
                        }
                    }
                }

                delete dense_subgraph;
                for (auto &i : *hg) delete i;
                delete hg;

            }
            else if(method == "scanst"){
                auto sketch_peel_pair = sketch_triangle_peel(meta_layer, hidden_edges, ractive, synopses_combined);
                resketch_count += sketch_peel_pair.second;
                auto dense_subgraph = sketch_peel_pair.first;

                for(unsigned int detected: *dense_subgraph) std::cout<<detected<<" "; std::cout<<std::endl;
                
                std::vector<std::vector<unsigned int> *> *hg = hidden_graph_construction(meta_layer, qp, g, visited, ractive);
                std::vector<bool> reserve(hg->size(), false);
                for(unsigned int detected: *dense_subgraph) reserve[detected] = true;

                for(unsigned int i=0;i<hg->size();i++){
                    if(reserve[i]){
                        for(unsigned int j=0;j<hg->at(i)->size();j++){
                            if(reserve[hg->at(i)->at(j)] && hg->at(i)->at(j) > i){
                                std::cout<<i<<" "<<hg->at(i)->at(j)<<std::endl;
                            }
                        }
                    }
                }
                
                delete dense_subgraph;
                for (auto &i : *hg) delete i;
                delete hg;
            }
        }
        resketch_count /= max_meta_layer;

        delete frontiers;
        delete peers;
        for (auto &i : *active) delete i;
        delete active;
        for (auto &i : *ractive) delete i;
        delete ractive;

        running_time /= 1000000000;
        std::cout<<"# time:"<<running_time<<std::endl;
        std::cout<<"# resketch_count:"<<resketch_count<<std::endl;
    }
    
    double Mem_prop(Pattern *qp, HeterGraph *g, std::vector<std::vector<bool>*>* visited,
                    std::vector<std::vector<bool>*>* back_visited){
        if(qp->instance >= 0 && qp->NTypes.size() <= 2) return 0;
        auto frontiers = new std::vector<unsigned int>();
        auto peers = new std::vector<unsigned int>();
        Peers(qp, g, frontiers, peers, visited);
        unsigned int peer_size = peers->size();

        if (peer_size <= 2) {
            delete frontiers;
            delete peers;
            return 0;
        }

        std::vector<std::set<unsigned int> *> *active = ActiveMidNodes(peers, frontiers, qp, g, visited, back_visited);
        auto ractive = new std::vector<std::vector<unsigned int> *>();
        for (unsigned int l = 0; l <= qp->ETypes.size(); l++) ractive->push_back(new std::vector<unsigned int>(g->NT.size(), 0));
        for (unsigned int n = 0; n < g->NT.size(); n++) for (unsigned int l : *active->at(n)) ractive->at(l)->at(n) = 1;

        auto synopses_combined = new std::vector<std::vector<jsy::Synopse> *>();
        for (unsigned int l = 0; l <= qp->ETypes.size(); l++) {
            synopses_combined->push_back(new std::vector<jsy::Synopse>());
            for (unsigned int p = 0; p < g->NT.size(); p++)
                if (ractive->at(l)->at(p) > 0) {
                    jsy::Synopse s; jsy::init_synopse(&s, p);
                    synopses_combined->at(l)->push_back(s);
                    ractive->at(l)->at(p) = synopses_combined->at(l)->size();
                }
        }
        double mem_count = 0.0;
        for(unsigned int l=0;l<=qp->ETypes.size();l++){
            mem_count += synopses_combined->at(l)->size();
        }
        mem_count *= (K * L * 4.0);
        mem_count /= 1024;
        return mem_count;
    }

    double Mem_hg(Pattern *qp, HeterGraph *g, std::vector<std::vector<bool>*>* visited,
            std::vector<std::vector<bool>*>* back_visited){
        if(qp->instance >= 0 && qp->NTypes.size() <= 2) return 0;
        auto frontiers = new std::vector<unsigned int>();
        auto peers = new std::vector<unsigned int>();
        Peers(qp, g, frontiers, peers, visited);
        unsigned int peer_size = peers->size();

        if (peer_size <= 2) {
            delete frontiers;
            delete peers;
            return 0;
        }

        std::vector<std::set<unsigned int> *> *active = ActiveMidNodes(peers, frontiers, qp, g, visited, back_visited);
        auto ractive = new std::vector<std::vector<unsigned int> *>();
        for (unsigned int l = 0; l <= qp->ETypes.size(); l++) ractive->push_back(new std::vector<unsigned int>(g->NT.size(), 0));
        for (unsigned int n = 0; n < g->NT.size(); n++) for (unsigned int l : *active->at(n)) ractive->at(l)->at(n) = 1;

        unsigned int max_meta_layer;
        if(qp->instance == -1) max_meta_layer = qp->ETypes.size();
        else max_meta_layer = qp->ETypes.size() - 1;

        double mem_count = 0.0;

        for(unsigned int meta_layer = 1; meta_layer <= max_meta_layer; meta_layer++){
            double temp_mem_count = 0.0;
            auto hg = hidden_graph_construction(meta_layer, qp, g, visited, ractive);
            for (auto &n : *hg) temp_mem_count += (n->size() + 1);
            if(temp_mem_count > mem_count){
                mem_count = temp_mem_count;
            }
        }
        mem_count *= 4.0;
        mem_count /= 1024;
        return mem_count;
    }

    void COD_prop(Pattern *qp, const std::string &metric, HeterGraph *g, std::vector<std::vector<bool>*>* visited,
                    std::vector<std::vector<bool>*>* back_visited, unsigned int &count){
        if(qp->instance >= 0 && qp->NTypes.size() <= 2) return;

        double running_time = 0.0;
        double resketch_count = 0.0;

        auto start = std::chrono::steady_clock::now();

        auto frontiers = new std::vector<unsigned int>();
        auto peers = new std::vector<unsigned int>();
        Peers(qp, g, frontiers, peers, visited);
        unsigned int peer_size = peers->size();

        auto end = std::chrono::steady_clock::now();

        if (peer_size <= 2) {
            delete frontiers;
            delete peers;
            return;
        }
        qp->print();

        running_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        start = std::chrono::steady_clock::now();

        std::vector<std::set<unsigned int> *> *active = ActiveMidNodes(peers, frontiers, qp, g, visited, back_visited);
        auto ractive = new std::vector<std::vector<unsigned int> *>();
        for (unsigned int l = 0; l <= qp->ETypes.size(); l++) ractive->push_back(new std::vector<unsigned int>(g->NT.size(), 0));
        for (unsigned int n = 0; n < g->NT.size(); n++) for (unsigned int l : *active->at(n)) ractive->at(l)->at(n) = 1;

        auto synopses_combined = new std::vector<std::vector<jsy::Synopse> *>();
        for (unsigned int l = 0; l <= qp->ETypes.size(); l++) {
            synopses_combined->push_back(new std::vector<jsy::Synopse>());
            for (unsigned int p = 0; p < g->NT.size(); p++)
                if (ractive->at(l)->at(p) > 0) {
                    jsy::Synopse s; jsy::init_synopse(&s, p);
                    synopses_combined->at(l)->push_back(s);
                    ractive->at(l)->at(p) = synopses_combined->at(l)->size();
                }
        }

        auto hidden_edges = new std::vector<HiddenEdge>();
        for (unsigned int l = 0; l < qp->ETypes.size(); l++) {
            for (unsigned int p = 0; p < g->NT.size(); p++) {
                if (ractive->at(l)->at(p) > 0) {
                    if (qp->EDirect[l] == 1) {
                        for (unsigned int nbr: *(g->EL[p])) {
                            if (ractive->at(l + 1)->at(nbr) > 0) {
                                HiddenEdge he{.s=ractive->at(l)->at(p) - 1, .t=ractive->at(l + 1)->at(nbr) - 1, .l=l};
                                hidden_edges->push_back(he);
                            }
                        }
                    } else {
                        for (unsigned int nbr: *(g->rEL[p])) {
                            if (ractive->at(l + 1)->at(nbr) > 0) {
                                HiddenEdge he{.s=ractive->at(l)->at(p) - 1, .t=ractive->at(l + 1)->at(nbr) - 1, .l=l};
                                hidden_edges->push_back(he);
                            }}}}}}
        unsigned int max_meta_layer;
        if(qp->instance == -1) max_meta_layer = qp->ETypes.size();
        else max_meta_layer = qp->ETypes.size() - 1;

        end = std::chrono::steady_clock::now();
        running_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        for (unsigned int meta_layer = 1; meta_layer <= max_meta_layer; meta_layer++) {
            if(metric == "edge") {
                start = std::chrono::steady_clock::now();
                auto sketch_peel_pair = sketch_peel(meta_layer, hidden_edges, ractive, synopses_combined);
                auto dense_subgraph = sketch_peel_pair.first;
                end = std::chrono::steady_clock::now();
                running_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
                
                std::cout<<"subgraph:";
                for(unsigned int u: *dense_subgraph) std::cout<<u<<" "; std::cout<<std::endl;

                delete dense_subgraph;
            }
            else if(metric == "triangle"){
                start = std::chrono::steady_clock::now();
                auto sketch_peel_pair = sketch_triangle_peel(meta_layer, hidden_edges, ractive, synopses_combined);
                auto dense_subgraph = sketch_peel_pair.first;
                end = std::chrono::steady_clock::now();
                running_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
                
                std::cout<<"subgraph:";
                for(unsigned int u: *dense_subgraph) std::cout<<u<<" "; std::cout<<std::endl;
                
                delete dense_subgraph;
            }
        }
        count += max_meta_layer;
        
        delete frontiers;
        delete peers;
        for (auto &i : *active) delete i;
        delete active;
        for (auto &i : *ractive) delete i;
        delete ractive;

        running_time /= 1000000000;
    }

    bool COD_hg(Pattern *qp, const std::string &metric, HeterGraph *g, std::vector<std::vector<bool>*>* visited,
                std::vector<std::vector<bool>*>* back_visited, double & construct_time, double& peel_time){
        double temp_construct_time = 0.0, temp_peel_time = 0.0;
        if(qp->instance >= 0 && qp->NTypes.size() <= 2) return false;

        auto start = std::chrono::steady_clock::now();

        auto frontiers = new std::vector<unsigned int>();
        auto peers = new std::vector<unsigned int>();
        Peers(qp, g, frontiers, peers, visited);
        unsigned int peer_size = peers->size();

        auto end = std::chrono::steady_clock::now();

        if (peer_size <= 2) {
            delete frontiers;
            delete peers;
            return false;
        }
        qp->print();

        temp_construct_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        start = std::chrono::steady_clock::now();

        std::vector<std::set<unsigned int> *> *active = ActiveMidNodes(peers, frontiers, qp, g, visited, back_visited);
        auto ractive = new std::vector<std::vector<unsigned int> *>();
        for (unsigned int l = 0; l <= qp->ETypes.size(); l++) ractive->push_back(new std::vector<unsigned int>(g->NT.size(), 0));
        for (unsigned int n = 0; n < g->NT.size(); n++) for (unsigned int l : *active->at(n)) ractive->at(l)->at(n) = 1;

        unsigned int max_meta_layer;
        if(qp->instance == -1) max_meta_layer = qp->ETypes.size();
        else max_meta_layer = qp->ETypes.size() - 1;

        end = std::chrono::steady_clock::now();
        temp_construct_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        
        for (unsigned int meta_layer = 1; meta_layer <= max_meta_layer; meta_layer++) {
            std::cout << "# l=" << meta_layer << std::endl;

            start = std::chrono::steady_clock::now();
            auto hg = hidden_graph_construction(meta_layer, qp, g, visited, ractive);
            end = std::chrono::steady_clock::now();
            temp_construct_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

            if(metric == "edge") {
                start = std::chrono::steady_clock::now();
                auto peel_pair = hg_peel(hg);
                end = std::chrono::steady_clock::now();
                temp_peel_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

                std::cout << "hg_peel_den:" << peel_pair.first << std::endl;
                std::cout << "hg_peel_size:" << peel_pair.second << std::endl;
            }
            else if(metric == "triangle"){
                start = std::chrono::steady_clock::now();
                for(unsigned int i=0;i<hg->size();i++)
                    if(!hg->at(i)->empty()) std::sort(hg->at(i)->begin(), hg->at(i)->end());
                end = std::chrono::steady_clock::now();
                temp_construct_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
                
                start = std::chrono::steady_clock::now();
                auto peel_pair = hg_triangle_peel(hg);
                end = std::chrono::steady_clock::now();
                temp_peel_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

                std::cout << "hg_peel_den:" << peel_pair.first << std::endl;
                std::cout << "hg_peel_size:" << peel_pair.second << std::endl;
            }
        }

        delete frontiers;
        delete peers;
        for (auto &i : *active) delete i;
        delete active;
        for (auto &i : *ractive) delete i;
        delete ractive;

        temp_construct_time /= 1000000000;
        temp_peel_time /= 1000000000;
        construct_time += temp_construct_time;
        peel_time += temp_peel_time;
        std::cout<<"construct_time:"<<temp_construct_time<<"s"<<std::endl;
        std::cout<<"peel_time:"<<temp_peel_time<<"s"<<std::endl;
        return true;
    }
}