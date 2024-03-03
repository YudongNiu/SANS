#ifndef JOINTSYNOPSE
#define JOINTSYNOPSE

#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <algorithm>
#include <queue>
#include <fstream>
#include <set>
#include <math.h>
#include <sstream>
#include <unordered_map>
#include <map>
#include <chrono>
#include "pattern.cpp"
#include "hin.cpp"
#include "param.h"
#include "hidden_edge.cpp"

namespace jsy {  //join-synopses
    struct Synopse {
        unsigned int minks[K * L];
        unsigned int sizes[L];
        bool overflows[L]; // if the node has more than k neighbors /* DSD */
        unsigned int p;

        unsigned int rs[L]; // the random values corresponding to "p"

        double clique_est;
    };
    
    void synopse_print(Synopse *s){
        std::cout<<"p:"<<s->p<<std::endl;
        std::cout<<"minks:";
        for(unsigned int l=0;l < L;l++){
            for(unsigned int k=0;k<K;k++) {
                std::cout << s->minks[l * K + k]<<" ";
            }
            std::cout<<" | ";
        }
        std::cout<<std::endl<<"sizes:";
        for (unsigned int size : s->sizes)
            std::cout<< size <<" ";
        std::cout<<std::endl<<"rs:";
        for (unsigned int r : s->rs)
            std::cout<< r <<" ";
        std::cout<<std::endl<<"overflow:";
        for (bool overflow : s->overflows)
            std::cout<< overflow <<" ";
        std::cout<<std::endl<<std::endl;
    }

    void init_synopse(Synopse *s, unsigned int _p) {
        for (unsigned int &mink : s->minks) mink = 0;
        for (unsigned int &size : s->sizes) size = 0;
        for (bool &overflow: s->overflows) overflow = false;
        s->p = _p;
        for(unsigned int &r: s->rs) r = 0;
        s->clique_est = 0;
    }

    void clear(Synopse *s, bool clear_rs) {
        for (unsigned int i = 0; i < L; i++) {
            for (unsigned int j = 0; j < s->sizes[i]; j++) s->minks[i * K + j] = 0;
            s->sizes[i] = 0;
            if(clear_rs) s->rs[i] = 0;
            s->overflows[i] = false;
        }
    }
    
    void clear(Synopse *s, unsigned int l){
        for(unsigned int j=0;j<s->sizes[l];j++) s->minks[l * K + j] = 0;
        s->sizes[l] = 0;
        s->overflows[l] = false;
    }

    void layer_clear(Synopse *s, unsigned int pl, bool clear_rs){
        for(unsigned int j=0;j<s->sizes[pl];j++) s->minks[pl * K + j] = 0;
        s->sizes[pl] = 0;
        if(clear_rs) s->rs[pl] = 0;
    }

    // choice a random nbr of "p" from its last synopse
    unsigned int choice(Synopse *s, unsigned int rand_num,
            std::vector<std::unordered_map<unsigned int, unsigned int>*>* rand2ps){
        unsigned int pos = rand_num % s->sizes[L-1];
        return rand2ps->at(L-1)->at(s->minks[K * (L-1) + pos]);
    }

//    unsigned int bundle_choice(Synopse *s, unsigned int rand_bundle, unsigned int rand_num,
//            std::vector<std::unordered_map<unsigned int, unsigned int>*>* rand2ps){
//        unsigned int bundle = rand_bundle % L;
//        unsigned int pos = rand_num % s->sizes[bundle];
//        return rand2ps->at(bundle)->at(s->minks[K * bundle + pos]);
//    }

    std::vector<unsigned int>* traverse(Synopse *s, unsigned int layer,
                                        std::vector<std::unordered_map<unsigned int, unsigned int>*>* rand2ps){
        auto sampl = new std::vector<unsigned int>();
        for(unsigned int i=0;i<s->sizes[layer];i++) {
            unsigned int value = s->minks[i + K * layer];
            sampl->push_back(rand2ps->at(layer)->at(value));
        }
        return sampl;
    }

    std::vector<unsigned int>* sample(Synopse *s, unsigned int num_sampl, unsigned int layer,
                                      std::vector<std::unordered_map<unsigned int, unsigned int>*>* rand2ps){
        std::random_device dev;
        std::mt19937 rg(dev());

        std::vector<unsigned int> values;
        for(unsigned int i=0;i<s->sizes[layer];i++) values.push_back(s->minks[i + K * layer]);
        std::shuffle(values.begin(), values.end(), rg);

        auto sampl = new std::vector<unsigned int>();
        if(num_sampl >= s->sizes[layer]){
            for(unsigned int i=0;i<s->sizes[layer];i++){
                unsigned int value = values[i];
                sampl->push_back(rand2ps->at(layer)->at(value));
            }
        }
        else{
            for(unsigned int i=0;i<num_sampl;i++){
                unsigned int value = values[i];
                sampl->push_back(rand2ps->at(layer)->at(value));
            }
        }
        return sampl;
    }

    void layer_add(Synopse *s, Synopse *t, unsigned int pl){
        for(unsigned int i=0;i<t->sizes[pl];i++){
            unsigned int x = t->minks[pl * K + i];
            bool exists = false;
            for(unsigned int j=0;j<s->sizes[pl];j++)
                if(s->minks[pl * K + j] == x){
                    exists = true;
                    break;
                }
            if(exists) continue;
            if(s->sizes[pl] < K){
                s->minks[pl * K + s->sizes[pl]] = x;
                s->sizes[pl]++;
                exists = true;
            }
            if(exists) continue;
            int top_pos = 0, top = 0;
            for(unsigned int j = 0;j<s->sizes[pl];j++){
                if(s->minks[pl * K + j] > top){
                    top = s->minks[pl * K + j];
                    top_pos = j;
                }
            }
            if(top > x) s->minks[pl * K + top_pos] = x;
        }
    }
    
    void add(Synopse *s, unsigned int l, unsigned int x){
        bool exists = false;
        for (unsigned int j = 0; j < s->sizes[l]; j++)
            if (s->minks[l * K + j] == x) {
                exists = true;
                break;
            }
        if (exists) return;
        if (s->sizes[l] < K) {
            s->minks[l * K + s->sizes[l]] = x;
            s->sizes[l]++;
            exists = true;
        }
        if (exists) return;
        int top_pos = 0, top = 0;
        for (unsigned int j = 0; j < s->sizes[l]; j++) {
            if (s->minks[l * K + j] > top) {
                top = s->minks[l * K + j];
                top_pos = j;
            }
        }
        if (top > x) s->minks[l * K + top_pos] = x;
    }

    void add(Synopse *s, Synopse *t) {
        // add Synopse t into Synopse s
        for (unsigned int l = 0; l < L; l++) {
            for (unsigned int i = 0; i < t->sizes[l]; i++) {
                unsigned int x = t->minks[l * K + i];
                bool exists = false;
                for (unsigned int j = 0; j < s->sizes[l]; j++)
                    if (s->minks[l * K + j] == x) {
                        exists = true;
                        break;
                    }
                if (exists) continue;
                if (s->sizes[l] < K) {
                    s->minks[l * K + s->sizes[l]] = x;
                    s->sizes[l]++;
                    exists = true;
                }
                if (exists) continue;
                int top_pos = 0, top = 0;
                for (unsigned int j = 0; j < s->sizes[l]; j++) {
                    if (s->minks[l * K + j] > top) {
                        top = s->minks[l * K + j];
                        top_pos = j;
                    }
                }
                if (top > x) s->minks[l * K + top_pos] = x;
            }
        }
    }

    double estimate(Synopse *s) {
        double est = 0;
        for (unsigned int l = 0; l < L; l++) {
            if (s->sizes[l] < K) est += s->sizes[l];
            else {
                int top = 0;
                for (unsigned int i = 0; i < s->sizes[l]; i++) if (s->minks[l * K + i] > top) top = s->minks[l * K + i];
                double r = static_cast<double>(top) / RAND_MAX;
                est += ((K - 1) / r);
                std::cout<<"temp_est:"<<((K - 1) / r)<<std::endl;//////////////////////////////////////
            }
        }
        est /= L;
        return est;
    }

    double intersect_estimate(Synopse* s, Synopse* t, unsigned int l){
        double est = 0;
        if(!s->overflows[l] && !t->overflows[l]){
            for(unsigned int i=0;i<t->sizes[l];i++){
                unsigned int x = t->minks[l * K + i];
                bool exists = false;
                for(unsigned int j=0;j<s->sizes[l];j++)
                    if(s->minks[l * K + j] == x){
                        exists = true;
                        break;
                    }
                if(exists) est++;
            }
        }
        else{
            unsigned int uset[2 * K];
            unsigned int uset_size=0;
            if(s->overflows[l] && t->overflows[l]){
                if(s->sizes[l]<= t->sizes[l]){
                    for(unsigned int i=0;i<s->sizes[l];i++) uset[i] = s->minks[l * K + i];
                    uset_size = s->sizes[l];
                    for(unsigned int i=0;i<t->sizes[l];i++){
                        unsigned int x = t->minks[l * K + i];
                        bool exists = false;
                        for(unsigned int j=0;j<uset_size;j++) if(uset[j] == x){
                                exists = true;
                                break;
                            }
                        if(exists) continue;
                        int top_pos = 0, top = 0;
                        for(unsigned int j=0;j<uset_size;j++) if(uset[j] > top){
                                top = uset[j];
                                top_pos = j;
                            }
                        if(top > x) uset[top_pos] = x;
                    }
                }
                else{
                    for(unsigned int i=0;i<t->sizes[l];i++) uset[i] = t->minks[l * K + i];
                    uset_size = t->sizes[l];
                    for(unsigned int i=0;i<s->sizes[l];i++){
                        unsigned int x = s->minks[l * K + i];
                        bool exists = false;
                        for(unsigned int j=0;j<uset_size;j++) if(uset[j] == x){
                                exists = true;
                                break;
                            }
                        if(exists) continue;
                        int top_pos = 0, top = 0;
                        for(unsigned int j=0;j<uset_size;j++) if(uset[j] > top){
                                top = uset[j];
                                top_pos = j;
                            }
                        if(top > x) uset[top_pos] = x;
                    }
                }
            }
            else if(s->overflows[l] && !t->overflows[l]){
                for(unsigned int i=0;i<s->sizes[l];i++) uset[i] = s->minks[l * K + i];
                uset_size = s->sizes[l];
                for(unsigned int i=0;i<t->sizes[l];i++){
                    unsigned int x = t->minks[l * K + i];
                    bool exists = false;
                    for(unsigned int j=0;j<uset_size;j++) if(uset[j] == x){
                            exists = true;
                            break;
                        }
                    if(exists) continue;
                    int top = 0;
                    for(unsigned int j=0;j<uset_size;j++) if(uset[j] > top) top = uset[j];
                    if(top > x){
                        uset[uset_size] = x;
                        uset_size++;
                    }
                }
            }
            else if(t->overflows[l] && !s->overflows[l]){
                for(unsigned int i=0;i<t->sizes[l];i++) uset[i] = t->minks[l * K + i];
                uset_size = t->sizes[l];
                for(unsigned int i=0;i<s->sizes[l];i++){
                    unsigned int x = s->minks[l * K + i];
                    bool exists = false;
                    for(unsigned int j=0;j<uset_size;j++) if(uset[j] == x){
                            exists = true;
                            break;
                        }
                    if(exists) continue;
                    int top = 0;
                    for(unsigned int j=0;j<uset_size;j++) if(uset[j] > top) top = uset[j];
                    if(top > x){
                        uset[uset_size] = x;
                        uset_size++;
                    }
                }
            }

            int top = 0;
            for(unsigned int i=0;i<uset_size;i++) if(uset[i] > top) top = uset[i];
            double r = static_cast<double>(top) / RAND_MAX;

            std::vector<unsigned int> insec;
            for(unsigned int i=0;i<t->sizes[l];i++){
                unsigned int x = t->minks[l * K + i];
                for(unsigned int j=0;j<s->sizes[l];j++)
                    if(s->minks[l * K + j] == x){
                        insec.push_back(x);
                        break;
                    }
            }

            double intersect = 0;
            for(auto x: insec){
                for (unsigned int j : uset) {
                    if(j == x){
                        intersect++;
                        break;
                    }
                }
            }
            est += (intersect / uset_size) * ( ( uset_size - 1) / r);
        }
        return est;
    }

    double intersect_estimate(Synopse* s, Synopse* t){
        double est = 0;
        for(unsigned int l=0;l<L;l++) est += intersect_estimate(s, t, l);
        est /= L;
        return est;
    }
    
    void set_overflow(std::vector<Synopse>* synopses, unsigned int l){
        for(auto &i: *synopses) if(i.sizes[l] >= K) i.overflows[l] = true;
    }

    //void set_overflow(std::vector<Synopse>* synopses, std::vector<bool>* qualified_nodes = nullptr){
    //    for(auto &i: *synopses){
    //        unsigned int _p = i.p;
    //        if(qualified_nodes == nullptr || qualified_nodes->at(_p))
    //            for(unsigned int l = 0;l<L;l++) if(i.sizes[l] >= K) i.overflows[l] = true;
    //    }
    //}
    
    void set_overflow(std::vector<Synopse>* synopses, std::vector<bool>* qualified_nodes=nullptr, bool rmself = true){
        for(auto &i: *synopses){
            unsigned int _p = i.p;
            if(qualified_nodes == nullptr || qualified_nodes->at(_p)){
                for(unsigned int l=0;l<L;l++){
                    if(i.sizes[l] >= K) i.overflows[l] = true;
                    if(rmself) {
                        unsigned int x = i.rs[l];
                        bool exists = false;
                        unsigned int j = 0;
                        for (; j < i.sizes[l]; j++)
                            if (i.minks[l * K + j] == x) {
                                exists = true;
                                break;
                            }
                        if (exists) {
                            i.minks[l * K + j] = i.minks[l * K + i.sizes[l] - 1];
                            i.minks[l * K + i.sizes[l] - 1] = 0;
                            i.sizes[l]--;
                        }
                    }
                }
            }
        }
    }

    double dyn_estimate(Synopse *s){
        double est = 0;
        for (unsigned int l=0;l<L;l++){
            if(!s->overflows[l]) est += s->sizes[l];
            else if(s->sizes[l] > 1) {
                int top = 0;
                for(unsigned int i=0; i<s->sizes[l]; i++) if(s->minks[l*K+i]>top) top=s->minks[l*K+i];
                double r = static_cast<double>(top) / RAND_MAX;
                est += ((s->sizes[l] - 1) / r);
            }
        }
        est /= L;
        return est;
    }

    
    bool remove(unsigned int x, unsigned int l, std::vector<Synopse>* synopses, std::vector<bool>* removed,
                std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>* reverse_synopses){
            //std::cout<<"x:"<<x<<std::endl;/////////////////////////////////////
            //if(reverse_synopses->at(l)->count(x) == 0) std::cout<<"x does not exists in reverse_synopses"<<std::endl;///////////////////////
            if(reverse_synopses->at(l)->count(x) == 0) return false;
            for(unsigned int pos: *reverse_synopses->at(l)->at(x)){
                if(!removed->at(pos)){
                    int j=0;
                    for(;j<synopses->at(pos).sizes[l];j++) if(synopses->at(pos).minks[l * K + j] == x) break;
                    synopses->at(pos).minks[l * K + j] = synopses->at(pos).minks[l * K + synopses->at(pos).sizes[l] - 1];
                    synopses->at(pos).minks[l * K + synopses->at(pos).sizes[l] - 1] = 0;
                    synopses->at(pos).sizes[l] -= 1;
                    if(synopses->at(pos).sizes[l] < Kmin && synopses->at(pos).overflows[l]) return true;
                }
            }
        return false;
    }
    
    void remove_nonresketch(unsigned int node, std::vector<Synopse>* synopses,
                            std::vector<std::vector<unsigned int> *>* ractive,
                            std::vector<bool>* qualified_nodes,
                            std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>* reverse_synopses,
                            std::unordered_map<unsigned int, bool>* updated_nodes){
        for(unsigned int l = 0;l < L; l++){
            unsigned int r = synopses->at(ractive->at(0)->at(node)-1).rs[l];
            if(reverse_synopses->at(l)->count(node) == 0) continue;
            for(unsigned int pos: *reverse_synopses->at(l)->at(node)){
                if(qualified_nodes->at(synopses->at(pos).p) && synopses->at(pos).sizes[l] > 0){
                    int j=0;
                    for(;j<synopses->at(pos).sizes[l];j++) if(synopses->at(pos).minks[l * K + j] == r) break;
                    synopses->at(pos).minks[l * K + j] = synopses->at(pos).minks[l * K + synopses->at(pos).sizes[l] - 1];
                    synopses->at(pos).minks[l * K + synopses->at(pos).sizes[l] - 1] = 0;
                    synopses->at(pos).sizes[l] -= 1;
                    if(updated_nodes->count(synopses->at(pos).p)==0) updated_nodes->insert({synopses->at(pos).p, true});
                }
            }
        }
    }

    // return whether re-sketching is required
    bool remove(unsigned int node, std::vector<Synopse>* synopses,
                std::vector<std::vector<unsigned int> *>* ractive,
                std::vector<bool>* qualified_nodes,
                std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>* reverse_synopses,
                std::unordered_map<unsigned int, bool>* updated_nodes){
        for(unsigned int l = 0;l < L; l++){
            unsigned int r = synopses->at(ractive->at(0)->at(node)-1).rs[l];
            if(reverse_synopses->at(l)->count(node) == 0) continue;
            for(unsigned int pos: *reverse_synopses->at(l)->at(node)){
                if(qualified_nodes->at(synopses->at(pos).p) && synopses->at(pos).sizes[l] > 0){
                    int j=0;
                    for(;j<synopses->at(pos).sizes[l];j++) if(synopses->at(pos).minks[l * K + j] == r) break;
                    synopses->at(pos).minks[l * K + j] = synopses->at(pos).minks[l * K + synopses->at(pos).sizes[l] - 1];
                    synopses->at(pos).minks[l * K + synopses->at(pos).sizes[l] - 1] = 0;
                    synopses->at(pos).sizes[l] -= 1;
                    if(synopses->at(pos).sizes[l] < Kmin && synopses->at(pos).overflows[l]) return true;
                    if(updated_nodes->count(synopses->at(pos).p) == 0) updated_nodes->insert({synopses->at(pos).p, true});
                }
            }
        }
        return false;
    }

    // NIU YUDONG: not finished
    //bool triangle_remove(unsigned int node, std::vector<Synopse>* synopses,
    //                     std::vector<std::vector<unsigned int> *>* ractive,
    //                     std::vector<bool>* qualified_nodes,
    //                     std::vector<std::unordered_map<unsigned int, std::vector<unsigned int>*>*>* reverse_synopses,
    //                     std::unordered_map<unsigned int, bool>* updated_nodes){
    //    for(unsigned int l=0;l < L;l++){
    //        unsigned int r = synopses->at(ractive->at(0)->at(node)-1).rs[l];
    //        for(unsigned int pos: *reverse_synopses->at(l)->at(node)){
    //            if(qualified_nodes->at(synopses->at(pos).p)){
    //                int j = 0;
    //                for(;j<synopses->at(pos).sizes[l];j++) if(synopses->at(pos).minks[l * K + j] == r) break;
    //                synopses->at(pos).minks[l * K + j] = synopses->at(pos).minks[l * K + synopses->at(pos).sizes[l] - 1];
    //                synopses->at(pos).minks[l * K + synopses->at(pos).sizes[l] - 1] = 0;
    //                synopses->at(pos).sizes[l] -= 1;
    //                if(synopses->at(pos).sizes[l] < Kmin && synopses->at(pos).overflows[l]) return true;
    //
    //            }
    //        }
    //    }
    //    return false;
    //}

    std::map<unsigned int, double>* synopses_fp_deg(unsigned int meta_layer, std::vector<HiddenEdge> *hidden_edges,
                                                std::vector<std::vector<Synopse> *>* synopses, double q_deg,
                                                double topr, unsigned int peer_size,
                                                std::vector<bool>* qualified_nodes = nullptr){
        auto degs = new std::map<unsigned int, double>();
        std::random_device dev;
        std::mt19937 generator(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);

        int meta_layer_count = -1;
        for (unsigned int ec = 0; ec < hidden_edges->size(); ec++) {
            unsigned int p = hidden_edges->at(ec).s;
            unsigned int nbr = hidden_edges->at(ec).t;
            unsigned int l = hidden_edges->at(ec).l;
            if(l == 0 && synopses->at(0)->at(p).sizes[0] == 0){
                unsigned int n=synopses->at(0)->at(p).p;
                if(qualified_nodes == nullptr || qualified_nodes->at(n))
                    for (unsigned int ll = 0; ll < L; ll++) {
                        int rand = distribute(generator);
                        synopses->at(0)->at(p).minks[ll * K] = (unsigned int) rand;
                        synopses->at(0)->at(p).sizes[ll] = 1;
                    }
            }
            if (l >= meta_layer) {
                meta_layer_count = ec - 1;
                break;
            }
            add(&synopses->at(l + 1)->at(nbr), &synopses->at(l)->at(p));
        }

        for (auto &i: *synopses->at(meta_layer)) {
            double clique_size = estimate(&i);
            if(clique_size >= q_deg && clique_size > peer_size * topr) return nullptr;
        }

        for(unsigned int l=1;l<meta_layer;l++) for(auto &i: *synopses->at(l)) i.clique_est = estimate(&i);

        if (meta_layer_count < 0) meta_layer_count = hidden_edges->size() - 1;

        for (unsigned int l = 0; l < meta_layer; l++) for (auto &i: *synopses->at(l)) clear(&i, true);

        unsigned int pre_nbr = hidden_edges->at((unsigned int)meta_layer_count).s;
        unsigned int pre_l = 1 + hidden_edges->at((unsigned int)meta_layer_count).l;

        for (int ec = meta_layer_count; ec >= 0; ec--) {
            unsigned int p = hidden_edges->at((unsigned int) ec).t;
            unsigned int nbr = hidden_edges->at((unsigned int) ec).s;
            unsigned int l = 1 + hidden_edges->at((unsigned int) ec).l;

            if((pre_nbr != nbr || pre_l != l) && l > 1){
                if(synopses->at(pre_l-1)->at(pre_nbr).clique_est > peer_size * topr &&
                    estimate(&synopses->at(pre_l-1)->at(pre_nbr))>=q_deg) return nullptr;
                pre_nbr = nbr;
                pre_l = l;
            }

            add(&synopses->at(l - 1)->at(nbr), &synopses->at(l)->at(p));
        }

        for (auto &i : *synopses->at(0)) {
            unsigned int _p = i.p;
            if (qualified_nodes == nullptr || qualified_nodes->at(_p))
                degs->insert(std::pair<unsigned int, double>(_p, estimate(&i)));
        }
        for (unsigned int l = 0; l <= meta_layer; l++) for (auto &i : *synopses->at(l)) clear(&i, true);
        return degs;
    }

    std::map<unsigned int, double>* synopses_fp(unsigned int meta_layer, std::vector<HiddenEdge> *hidden_edges,
                                                std::vector<std::vector<Synopse> *>* synopses, double q_deg,
                                                double topr, unsigned int peer_size, double & mcs,
                                                std::vector<bool> *qualified_nodes = nullptr){
        auto degs = new std::map<unsigned int, double>();
        std::random_device dev;
        std::mt19937 generator(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);

        int meta_layer_count = -1;
        for (unsigned int ec = 0; ec < hidden_edges->size(); ec++) {
            unsigned int p = hidden_edges->at(ec).s;
            unsigned int nbr = hidden_edges->at(ec).t;
            unsigned int l = hidden_edges->at(ec).l;
            if(l == 0 && synopses->at(0)->at(p).sizes[0] == 0){
                unsigned int n=synopses->at(0)->at(p).p;
                if(qualified_nodes == nullptr || qualified_nodes->at(n))
                    for (unsigned int ll = 0; ll < L; ll++) {
                        int rand = distribute(generator);
                        synopses->at(0)->at(p).minks[ll * K] = (unsigned int) rand;
                        synopses->at(0)->at(p).sizes[ll] = 1;
                    }
            }
            if (l >= meta_layer) {
                meta_layer_count = ec - 1;
                break;
            }
            add(&synopses->at(l + 1)->at(nbr), &synopses->at(l)->at(p));
        }

        mcs = 0;

        for (auto &i: *synopses->at(meta_layer)) {
            unsigned int _p = i.p;
            double clique_size = estimate(&i);
            if(clique_size > mcs) mcs = clique_size;
            if(clique_size >= q_deg && clique_size > peer_size * topr) return nullptr;
        }

        if (meta_layer_count < 0) meta_layer_count = hidden_edges->size() - 1;

        for (unsigned int l = 0; l < meta_layer; l++) for (auto &i: *synopses->at(l)) clear(&i, true);

        for (int ec = meta_layer_count; ec >= 0; ec--) {
            unsigned int p = hidden_edges->at((unsigned int) ec).t;
            unsigned int nbr = hidden_edges->at((unsigned int) ec).s;
            unsigned int l = 1 + hidden_edges->at((unsigned int) ec).l;
            add(&synopses->at(l - 1)->at(nbr), &synopses->at(l)->at(p));
        }

        for (auto &i : *synopses->at(0)) {
            unsigned int _p = i.p;
            if (qualified_nodes == nullptr || qualified_nodes->at(_p))
                degs->insert(std::pair<unsigned int, double>(_p, estimate(&i)));
        }
        for (unsigned int l = 0; l <= meta_layer; l++) for (auto &i : *synopses->at(l)) clear(&i, true);
        return degs;
    }

//    std::vector<std::unordered_map<unsigned int, unsigned int>*>* gnn_synopses_bipartite(unsigned int meta_layer,
//                                                                    std::vector<HiddenEdge> *hidden_edges,
//                                                                    std::vector<std::vector<Synopse> *>* synopses){
//        std::random_device dev;
//        std::mt19937 generator(dev());
//        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);
//
//        auto rand2ps = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
//        for(unsigned int l=0;l<L;l++) rand2ps->push_back(new std::unordered_map<unsigned int, unsigned int>());
//
//        int meta_layer_count = -1;
//        for(unsigned int ec = 0;ec<hidden_edges->size();ec++){
//            unsigned int p = hidden_edges->at(ec).s;
//            unsigned int nbr = hidden_edges->at(ec).t;
//            unsigned int l = hidden_edges->at(ec).l;
//            if(l==0 && synopses->at(0)->at(p).sizes[0] == 0){
//                unsigned int n = synopses->at(0)->at(p).p;
//                for(unsigned int ll = 0; ll <L ; ll++){
//                    int rand = distribute(generator);
//                    while(rand2ps->at(ll)->count(rand) > 0) rand = distribute(generator);
//                    synopses->at(0)->at(p).minks[ll * K] = (unsigned int) rand;
//                    synopses->at(0)->at(p).sizes[ll] = 1;
//                    rand2ps->at(ll)->at(rand) = n;
//                    synopses->at(0)->at(p).rs[ll] = (unsigned int)rand;
//                }
//            }
//            if(l >= meta_layer){
//                meta_layer_count = ec - 1;
//                break;
//            }
//            add(&synopses->at(l+1)->at(nbr), &synopses->at(l)->at(p));
//        }
//
//        if(meta_layer_count < 0) meta_layer_count = hidden_edges->size()-1;
//
//        for(unsigned int l=0;l<meta_layer;l++) for(auto &i: *synopses->at(l)) clear(&i, false);
//
//        for(int ec = meta_layer_count;ec >= 0;ec--){
//            unsigned int p = hidden_edges->at((unsigned int) ec).t;
//            unsigned int nbr = hidden_edges->at((unsigned int) ec).s;
//            unsigned int l = 1+hidden_edges->at((unsigned int) ec).l;
//            if(l == meta_layer && synopses->at(meta_layer)->at(p).sizes[0] == 0){
//                unsigned int n = synopses->at(meta_layer)->at(p).p;
//                for(unsigned int ll = 0; ll < L;ll++){
//                    int rand = distribute(generator);
//                    while(rand2ps->at(ll)->count(rand) > 0) rand = distribute(generator);
//                    synopses->at(meta_layer)->at(p).minks[ll * K] = (unsigned int) rand;
//                    synopses->at(meta_layer)->at(p).sizes[ll] = 1;
//                    rand2ps->at(ll)->at(rand) = n;
//                    synopses->at(meta_layer)->at(p).rs[ll] = (unsigned int) rand;
//                }
//            }
//            add(&synopses->at(l-1)->at(nbr), &synopses->at(l)->at(p));
//        }
//        return rand2ps;
//    }

    std::pair<unsigned int, unsigned int> sub_matching_graph(unsigned int meta_layer,
            std::vector<HiddenEdge>* hidden_edges,
            std::vector<std::vector<bool>*>* qualified_nodes,
            std::vector<std::vector<Synopse>*>* synopses){
        unsigned int sub_matching_graph_size = 0, matching_graph_size = 0;
        for (auto &hidden_edge : *hidden_edges) {
            unsigned int p = hidden_edge.s;
            unsigned int nbr = hidden_edge.t;
            unsigned int l = hidden_edge.l;
            unsigned int _p = synopses->at(l)->at(p).p;
            unsigned int _nbr = synopses->at(l+1)->at(nbr).p;
            if(l >= meta_layer) break;
            std::cout<<"p:"<<p<<"  nbr:"<<nbr<<"  l:"<<l<<"  _p:"<<_p<<std::endl;//////////////////////
            if(qualified_nodes->at(l)->at(_p)){
                sub_matching_graph_size++;
                qualified_nodes->at(l+1)->at(_nbr) = true;
            }
            matching_graph_size++;
            
            
        }
        return {sub_matching_graph_size, matching_graph_size};
    }

    std::unordered_map<unsigned int, unsigned int>* layer_synopses(unsigned int meta_layer,
                                                                    std::vector<HiddenEdge> *hidden_edges,
                                                                    std::vector<std::vector<Synopse>*>* synopses,
                                                                    unsigned int pl,
                                                                    std::vector<bool> *qualified_nodes = nullptr){
        std::random_device dev;
        std::mt19937 generator(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);

        auto rand2p = new std::unordered_map<unsigned int, unsigned int>();

        int meta_layer_count = -1;
        for(unsigned int ec = 0;ec < hidden_edges->size(); ec++){
            unsigned int p = hidden_edges->at(ec).s;
            unsigned int nbr = hidden_edges->at(ec).t;
            unsigned int l = hidden_edges->at(ec).l;
            if(l == 0 && synopses->at(0)->at(p).sizes[pl] == 0){
                unsigned int n = synopses->at(0)->at(p).p;
                if(qualified_nodes == nullptr || qualified_nodes->at(n)){
                    int rand = distribute(generator);
                    while(rand2p->count(rand) > 0) rand = distribute(generator);
                    synopses->at(0)->at(p).minks[pl * K] = (unsigned int) rand;
                    synopses->at(0)->at(p).sizes[pl] = 1;
                    rand2p->insert(std::make_pair((unsigned int)rand, n));
                    synopses->at(0)->at(p).rs[pl] = (unsigned int) rand;
                }
            }
            if(l >= meta_layer){
                meta_layer_count = ec - 1;
                break;
            }
            layer_add(&synopses->at(l + 1)->at(nbr), &synopses->at(l)->at(p), pl);
        }

        if(meta_layer_count < 0) meta_layer_count = hidden_edges->size() - 1;

        for(unsigned int l = 0; l < meta_layer; l++) for(auto &i: *synopses->at(l)) layer_clear(&i, pl, false);

        for(int ec = meta_layer_count;ec >= 0;ec--){
            unsigned int p = hidden_edges->at((unsigned int) ec).t;
            unsigned int nbr = hidden_edges->at((unsigned int) ec).s;
            unsigned int l = 1 + hidden_edges->at((unsigned int) ec).l;
            layer_add(&synopses->at(l-1)->at(nbr), &synopses->at(l)->at(p), pl);
        }
        return rand2p;
    }

    std::vector<std::unordered_map<unsigned int, unsigned int>*> *gnn_synopses(unsigned int meta_layer,
                                                std::vector<HiddenEdge> *hidden_edges,
                                                std::vector<std::vector<Synopse> *> *synopses,
                                                std::vector<bool> *qualified_nodes = nullptr){
        std::random_device dev;
        std::mt19937 generator(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);

        auto rand2ps = new std::vector<std::unordered_map<unsigned int, unsigned int>*>();
        for(unsigned int l=0;l<L;l++) rand2ps->push_back(new std::unordered_map<unsigned int, unsigned int>());

        int meta_layer_count = -1;
        for (unsigned int ec = 0; ec < hidden_edges->size(); ec++) {
            unsigned int p = hidden_edges->at(ec).s;
            unsigned int nbr = hidden_edges->at(ec).t;
            unsigned int l = hidden_edges->at(ec).l;
            if(l == 0 && synopses->at(0)->at(p).sizes[0] == 0){
                unsigned int n=synopses->at(0)->at(p).p;
                if(qualified_nodes == nullptr || qualified_nodes->at(n))
                    for (unsigned int ll = 0; ll < L; ll++){
                        int rand = distribute(generator);
                        while(rand2ps->at(ll)->count(rand) > 0) rand = distribute(generator);
                        synopses->at(0)->at(p).minks[ll * K] = (unsigned int) rand;
                        synopses->at(0)->at(p).sizes[ll] = 1;
                        rand2ps->at(ll)->insert(std::make_pair((unsigned int)rand, n));
                        synopses->at(0)->at(p).rs[ll] = (unsigned int)rand;
                    }
            }
            if (l >= meta_layer) {
                meta_layer_count = ec - 1;
                break;
            }
            add(&synopses->at(l + 1)->at(nbr), &synopses->at(l)->at(p));
        }

        if (meta_layer_count < 0) meta_layer_count = hidden_edges->size() - 1;

        for (unsigned int l = 0; l < meta_layer; l++) for (auto &i: *synopses->at(l)) clear(&i, false);

        for (int ec = meta_layer_count; ec >= 0; ec--){
            unsigned int p = hidden_edges->at((unsigned int) ec).t;
            unsigned int nbr = hidden_edges->at((unsigned int) ec).s;
            unsigned int l = 1 + hidden_edges->at((unsigned int) ec).l;
            add(&synopses->at(l - 1)->at(nbr), &synopses->at(l)->at(p));
        }
        return rand2ps;
    }

    std::map<unsigned int, double> *synopses(unsigned int meta_layer, std::vector<HiddenEdge> *hidden_edges,
                                             std::vector<std::vector<Synopse> *> *synopses,
                                             std::vector<bool> *qualified_nodes = nullptr) {
        auto degs = new std::map<unsigned int, double>();
        std::random_device dev;
        std::mt19937 generator(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distribute(1, RAND_MAX);

        int meta_layer_count = -1;
        for (unsigned int ec = 0; ec < hidden_edges->size(); ec++) {
            unsigned int p = hidden_edges->at(ec).s;
            unsigned int nbr = hidden_edges->at(ec).t;
            unsigned int l = hidden_edges->at(ec).l;
            if(l == 0 && synopses->at(0)->at(p).sizes[0] == 0){
                unsigned int n=synopses->at(0)->at(p).p;
                if(qualified_nodes == nullptr || qualified_nodes->at(n))
                    for (unsigned int ll = 0; ll < L; ll++) {
                        int rand = distribute(generator);
                        synopses->at(0)->at(p).minks[ll * K] = (unsigned int) rand;
                        synopses->at(0)->at(p).sizes[ll] = 1;
                    }
            }
            if (l >= meta_layer) {
                meta_layer_count = ec - 1;
                break;
            }
            add(&synopses->at(l + 1)->at(nbr), &synopses->at(l)->at(p));
        }

        if (meta_layer_count < 0) meta_layer_count = hidden_edges->size() - 1;

        for (unsigned int l = 0; l < meta_layer; l++) for (auto &i: *synopses->at(l)) clear(&i, true);

        for (int ec = meta_layer_count; ec >= 0; ec--) {
            unsigned int p = hidden_edges->at((unsigned int) ec).t;
            unsigned int nbr = hidden_edges->at((unsigned int) ec).s;
            unsigned int l = 1 + hidden_edges->at((unsigned int) ec).l;
            add(&synopses->at(l - 1)->at(nbr), &synopses->at(l)->at(p));
        }

        for (auto &i : *synopses->at(0)) {
            unsigned int _p = i.p;
            if (qualified_nodes == nullptr || qualified_nodes->at(_p))
                degs->insert(std::pair<unsigned int, double>(_p, estimate(&i)));
        }
        for (unsigned int l = 0; l <= meta_layer; l++) for (auto &i : *synopses->at(l)) clear(&i, true);
        return degs;
    }
}
#endif