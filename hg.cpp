#ifndef SEMANTICGRAPH
#define SEMANTICGRAPH

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
#include <map>
#include <chrono>
#include "pattern.cpp"
#include "peer.cpp"
#include "hin.cpp"
#include "param.h"

std::vector<std::vector<unsigned int>*>* hidden_graph_construction_noselfloop(
        unsigned int meta_layer, Pattern *qp, HeterGraph *g, std::vector<std::vector<bool>*>* visited,
        std::vector<std::vector<unsigned int>*>* ractive){
    auto visited_mod = new std::vector<std::pair<unsigned int, unsigned int>>();
    auto hidden_graph = new std::vector<std::vector<unsigned int>*>();

    for(unsigned int u=0;u<g->NT.size();u++){
        hidden_graph->push_back(new std::vector<unsigned int>());
        if(ractive->at(0)->at(u) > 0){
            auto midfrontiers = MidHop(u, qp, g, meta_layer, visited, visited_mod, ractive);
            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();

            for(unsigned int f: *midfrontiers){
                auto peers_new = rMidHop(f, qp, g, meta_layer, visited, visited_mod, ractive);
                for(unsigned int p: *peers_new)
                    if(p != u) hidden_graph->at(u)->push_back(p);
                delete peers_new;
            }

            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();
            delete midfrontiers;
        }
    }
    delete visited_mod;
    return hidden_graph;
}

std::vector<std::vector<unsigned int>*>* hidden_graph_construction(unsigned int meta_layer, Pattern *qp, HeterGraph *g,
                                                                   std::vector<std::vector<bool>*>* visited,
                                                                   std::vector<std::vector<unsigned int>*>* ractive){
    auto visited_mod = new std::vector<std::pair<unsigned int, unsigned int>>();
    auto hidden_graph = new std::vector<std::vector<unsigned int>*>();

    for(unsigned int u=0;u<g->NT.size();u++){
        hidden_graph->push_back(new std::vector<unsigned int>());
        if(ractive->at(0)->at(u) > 0){
            auto midfrontiers = MidHop(u, qp, g, meta_layer, visited, visited_mod, ractive);
            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();

            for(unsigned int f: *midfrontiers){
                auto peers_new = rMidHop(f, qp, g, meta_layer, visited, visited_mod, ractive);
                for(unsigned int p: *peers_new) hidden_graph->at(u)->push_back(p);
                delete peers_new;
            }

            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();
            delete midfrontiers;
        }
    }
    delete visited_mod;
    return hidden_graph;
}

std::vector<unsigned int>* hidden_graph_deg_union(unsigned int meta_layer, Pattern *qp, HeterGraph *g,
                                                  std::vector<std::vector<bool>*>* visited,
                                                  std::vector<std::vector<unsigned int>*>* ractive){
    auto deg_vec = new std::vector<unsigned int>(g->NT.size(), 0);

    std::vector<std::vector<unsigned int>*> x0_2_xL(g->NT.size(), new std::vector<unsigned int>());
    std::vector<std::vector<unsigned int>*> xL_2_x0(g->NT.size(), new std::vector<unsigned int>());

    auto visited_mod = new std::vector<std::pair<unsigned int, unsigned int>>();
    for(unsigned int u=0;u<g->NT.size();u++){
        if(ractive->at(0)->at(u) > 0){
            auto midfrontiers = MidHop(u, qp, g, meta_layer, visited, visited_mod, ractive);
            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();

            for(unsigned int f: *midfrontiers){
                x0_2_xL[u]->push_back(f);
                xL_2_x0[f]->push_back(u);
            }
        }
    }
    for(unsigned int f=0;f<g->NT.size();f++)
        if(ractive->at(meta_layer)->at(f)>0) std::sort(xL_2_x0[f]->begin(), xL_2_x0[f]->end(), std::greater<>());

    std::priority_queue<std::pair<unsigned int, unsigned int>> pq;
    for(unsigned int u=0;u<g->NT.size();u++){
        if(ractive->at(0)->at(u) > 0){
            std::vector<unsigned int> image_union;

            std::vector<unsigned int> positions;
            for(unsigned int i=0;i<x0_2_xL[u]->size();i++) positions.push_back(0);
            for(unsigned int i=0;i<x0_2_xL[u]->size();i++) {
                unsigned int f = x0_2_xL[u]->at(i);
                pq.push(std::pair<unsigned int, unsigned int>(xL_2_x0[f]->at(positions[i]), i));
                positions[i]++;
            }
            while(!pq.empty()){
                auto next = pq.top();
                pq.pop();
                unsigned int n = next.first;
                unsigned int pos = next.second;
                unsigned int f = x0_2_xL[u]->at(pos);

                positions[pos]++;
                if(positions[pos] < xL_2_x0[f]->size()){
                    pq.push(std::pair<unsigned int, unsigned int>(xL_2_x0[f]->at(positions[pos]), pos));
                }
                if(image_union[image_union.size()-1] > n){
                    image_union.push_back(n);
                }
            }
            deg_vec->at(u) = image_union.size();
        }
    }
    return deg_vec;
}

std::vector<std::vector<unsigned int>*>* subgraph(unsigned int meta_layer, Pattern *qp, HeterGraph *g,
                                    std::vector<std::vector<bool>*>* visited,
                                    std::vector<std::vector<unsigned int>*>* ractive,
                                    std::vector<bool>* qualified_nodes){
    auto visited_mod = new std::vector<std::pair<unsigned int, unsigned int>>();
    auto subgraph = new std::vector<std::vector<unsigned int>*>();

    for(unsigned int u=0;u<g->NT.size();u++){
        subgraph->push_back(new std::vector<unsigned int>());
        if(qualified_nodes->at(u)){
            auto midfrontiers = MidHop(u, qp, g, meta_layer, visited, visited_mod, ractive);
            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();

            for(unsigned int f: *midfrontiers){
                auto peers_new = rMidHop(f, qp, g, meta_layer, visited, visited_mod, ractive);
                for(unsigned int p: *peers_new)
                    if(qualified_nodes->at(p)) subgraph->at(u)->push_back(p);
                delete peers_new;
            }
            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();
            delete midfrontiers;
        }
    }
    delete visited_mod;
    return subgraph;
}

std::vector<unsigned int>* hidden_graph_deg(unsigned int meta_layer, Pattern *qp, HeterGraph *g,
                                            std::vector<std::vector<bool>*>* visited,
                                            std::vector<std::vector<unsigned int>*>* ractive){
    auto deg_vec = new std::vector<unsigned int>(g->NT.size(), 0);

    auto visited_mod = new std::vector<std::pair<unsigned int, unsigned int>>();
    for(unsigned int u=0;u<g->NT.size();u++){
        if(ractive->at(0)->at(u) > 0){
            unsigned int degree_u = 0;
            auto midfrontiers = MidHop(u, qp, g, meta_layer, visited, visited_mod, ractive);
            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();

            for(unsigned int f: *midfrontiers){
                auto peers_new = rMidHop(f, qp, g, meta_layer, visited, visited_mod, ractive);
                degree_u += peers_new->size();
                delete peers_new;
            }

            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();
            delete midfrontiers;

            deg_vec->at(u) = degree_u;
        }
    }
    delete visited_mod;
    return deg_vec;
}

std::vector<unsigned int>* hidden_graph_rw(unsigned int qn, double alpha, std::vector<std::vector<unsigned int>*>* hg){
    auto rw = new std::vector<unsigned int>();
    rw->push_back(qn);

    std::random_device dev;
    std::mt19937 generator(dev());
    std::uniform_int_distribution<std::mt19937::result_type> distrib(0, RAND_MAX);

    double threshold = ((distrib(generator) * 1.0) / RAND_MAX);
    while(threshold > alpha){
        unsigned int u = rw->at(rw->size()-1);
        unsigned int rand_pos = distrib(generator) % hg->at(u)->size();
        rw->push_back(hg->at(u)->at(rand_pos));
        threshold = ((distrib(generator) * 1.0) / RAND_MAX);
    }
    return rw;
}

std::vector<unsigned int>* hidden_graph_rw(unsigned int meta_layer, Pattern *qp, unsigned int qn, double alpha,
        HeterGraph *g, std::vector<std::vector<bool>*>* visited, std::vector<std::vector<unsigned int>*>* ractive){
    /* alpha: teleportation */

    auto rw = new std::vector<unsigned int>();
    rw->push_back(qn);
    auto visited_mod = new std::vector<std::pair<unsigned int, unsigned int>>();
    std::random_device dev;
    std::mt19937 generator(dev());
    std::uniform_real_distribution<> dist1(0.0, 1.0);
    std::uniform_int_distribution<std::mt19937::result_type> intdist(0, RAND_MAX);

    double threshold = dist1(generator);
    while(threshold > alpha){
        std::cout<<"threshold:"<<threshold<<std::endl;

        unsigned int u = rw->at(rw->size()-1);
        std::vector<unsigned int> nbrs;
        auto midfrontiers = MidHop(u, qp, g, meta_layer, visited, visited_mod, ractive);
        for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
        visited_mod->clear();

        for(unsigned int f: *midfrontiers){
            auto peers_new = rMidHop(f, qp, g, meta_layer, visited, visited_mod, ractive);
            for(unsigned int p: *peers_new) nbrs.push_back(p);
            delete peers_new;
        }
        for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
        visited_mod->clear();
        delete midfrontiers;

        unsigned int nextstep = nbrs[intdist(generator) % nbrs.size()];
        rw->push_back(nextstep);
        threshold = dist1(generator);
    }
    std::cout<<"threshold:"<<threshold<<std::endl;
    return rw;
}

std::vector<unsigned int>* hidden_graph_hindex(unsigned int meta_layer, Pattern *qp, HeterGraph *g,
                                               std::vector<std::vector<bool>*>* visited,
                                               std::vector<std::vector<unsigned int>*>* ractive){
    auto deg_vec = hidden_graph_deg(meta_layer, qp, g, visited, ractive);

    auto hindex_vec = new std::vector<unsigned int>(g->NT.size(), 0);

    auto visited_mod = new std::vector<std::pair<unsigned int, unsigned int>>();
    for(unsigned int u=0;u<g->NT.size();u++){
        if(ractive->at(0)->at(u) > 0){
            std::vector<unsigned int> nbrs;
            auto midfrontiers = MidHop(u, qp, g, meta_layer, visited, visited_mod, ractive);
            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();

            for(unsigned int f: *midfrontiers){
                auto peers_new = rMidHop(f, qp, g, meta_layer, visited, visited_mod, ractive);
                for(unsigned int p: *peers_new) nbrs.push_back(p);
                delete peers_new;
            }

            for(auto &it: *visited_mod) visited->at(it.first)->at(it.second) = false;
            visited_mod->clear();
            delete midfrontiers;

            std::vector<unsigned int> degrees;
            for(unsigned int n: nbrs) degrees.push_back(deg_vec->at(n));
            std::sort(degrees.begin(), degrees.end(), std::greater<>());

            unsigned int hindex = 0;
            for(;hindex<degrees.size();hindex++){
                if(degrees[hindex]<(hindex + 1)) break;
            }
            hindex_vec->at(u) = hindex;
        }
    }
    delete deg_vec;
    delete visited_mod;
    return hindex_vec;
}

#endif