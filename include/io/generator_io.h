/*******************************************************************************
 * include/io/generator_io.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _GENERATOR_IO_H_
#define _GENERATOR_IO_H_

#include <mpi.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <vector>

#include "generator_config.h"

namespace kagen {

template <typename T>
struct identity {
  typedef T type;
};

template <typename Edge = std::tuple<SInt, SInt>>
class GeneratorIO {
 public:
  GeneratorIO(PGeneratorConfig& config) : config_(config), local_num_edges_(0) {
    dist_.resize(config_.dist_size);
  }

  inline void UpdateDist(SInt node_id) {
    // if ((CRCHash::hash(node_id) % config_.n) < dist_.size()) dist_[node_id]++;
    if (node_id < dist_.size()) dist_[node_id]++;
    local_num_edges_++;
  }

  void OutputDist() const {
    // Exchange local dist
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<SInt> global_dist(dist_.size(), 0);
    MPI_Reduce(&dist_[0], &global_dist[0], dist_.size(), MPI_LONG, MPI_SUM,
               ROOT, MPI_COMM_WORLD);
    if (rank == ROOT) {
      FILE* fout = fopen(config_.output_file.c_str(), "w+");
      for (SInt i = 0; i < global_dist.size(); ++i) {
        fprintf(fout, "%llu\n", global_dist[i]);
      }
      fclose(fout);
    }
  }

  void ReserveEdges(SInt num_edges) { edges_.reserve(num_edges); }

  template <typename... Args>
  inline void PushEdge(Args... args) {
    edges_.emplace_back(std::make_tuple(args...));
    local_num_edges_++;
  }

  //TODO: is this needed? if not, remove
  inline void PushEdge(SInt v1, SInt v2){
	if( v1 < v2){
		edges_.emplace_back( std::make_tuple( v1, v2) );
	}else{
		edges_.emplace_back( std::make_tuple( v2, v1) );
	}
	//v1 > v2 ? edges_.emplace_back( std::make_tuple( v1, v2) ) : edges_.emplace_back( std::make_tuple( v2, v1) );
    local_num_edges_++;
  }

  void OutputEdges() const { 
#ifdef SINGLE_LIST
    GatherPrint(identity<Edge>());
#else
    Print(identity<Edge>()); 
#endif
  }

  SInt NumEdges() const { 
    return edges_.size() > 0 ? edges_.size() : local_num_edges_/2; 
  }

 private:
  PGeneratorConfig &config_;

  std::vector<SInt> dist_;
  std::vector<Edge> edges_;

  SInt local_num_edges_;
  

  void GatherPrint(identity<std::tuple<LPFloat, LPFloat, LPFloat, SInt>>) const {
    throw std::logic_error("Not implemented, should not be called.");
  }

  typedef std::tuple<LPFloat, LPFloat> coord_2d;
  typedef std::tuple<LPFloat, LPFloat, LPFloat> coord_3d;


  //template<typename coord_dd>
  int GatherCoords( identity<std::tuple<LPFloat, LPFloat, SInt>>, std::vector<std::tuple<LPFloat, LPFloat>> &globalCoords, std::vector<SInt> &globalIds) const {
    // Exchange local dist
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //cannot remove duplicates locally as function is const  and cannot change edges_

    //break vector of tuples into two vectors
    unsigned int lSize = NumEdges();
    std::vector<std::tuple<LPFloat, LPFloat>> localCoords(lSize);
    std::vector<SInt> localIds(lSize);
    assert( lSize==edges_.size() );

    for ( unsigned int i=0; i<lSize; i++ ){
        auto edge = edges_[i];
        localCoords[i] = std::tuple<LPFloat, LPFloat>( std::get<0>(edge) , std::get<1>(edge) );
        localIds[i] = std::get<2>(edge);
    }

    int total_num_edges = 0;
    
    // Gather number of edges for each PE
    std::vector<int> displ(size);
    std::vector<int> num_edges(size);

    MPI_Gather(&lSize, 1, MPI_INT,
               num_edges.data(), 1, MPI_INT, 
               ROOT, MPI_COMM_WORLD);
    int current_displ = 0;
    if (rank == ROOT) {
      for (SInt i = 0; i < num_edges.size(); ++i) {
        displ[i] = current_displ;
        total_num_edges += num_edges[i];
        current_displ = total_num_edges;
      }
    }

    std::cout<< __FILE__ << ", " << __LINE__ << ": " << rank << ", local size= " << lSize <<  ", total_num_edges= " << total_num_edges << std::endl;
    globalCoords.resize(total_num_edges);
    globalIds.resize(total_num_edges);

    // Gather actual coordinates
    MPI_Datatype MPI_COORD;
    MPI_Type_vector(1, 2, 0, MPI_DOUBLE, &MPI_COORD);
    MPI_Type_commit(&MPI_COORD);
    
    MPI_Gatherv(localCoords.data(), lSize, MPI_COORD,
                globalCoords.data(), num_edges.data(), displ.data(), MPI_COORD, 
                ROOT, MPI_COMM_WORLD);  

    // Gather actual vertex indices
    MPI_Gatherv(localIds.data(), lSize, MPI_LONG,
                globalIds.data(), num_edges.data(), displ.data(), MPI_LONG, 
                ROOT, MPI_COMM_WORLD);

    assert( globalCoords.size()==globalIds.size() );

    return total_num_edges;
  }


  //2D point output
  void GatherPrint(identity<std::tuple<LPFloat, LPFloat, SInt>>) const {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<std::tuple<LPFloat, LPFloat>> globalCoords; 
    std::vector<SInt> globalIds; 

    GatherCoords( identity<Edge>(), globalCoords, globalIds);

    if (rank == ROOT) {
        //checks
        if( not std::is_sorted(globalIds.begin(), globalIds.end()) ){
            //throw std::logic_error("ERROR: For some reason, the coordinate indices are not sorted.");
            std::cout<< "WARNING: indices are not sorted, will sort now" << std::endl;

            //zip again to sort and remove duplicates
            std::vector<std::tuple<LPFloat, LPFloat, SInt>> forSort;
            for( unsigned long i=0; i<globalCoords.size(); i++){
                forSort.push_back( std::tuple<LPFloat, LPFloat, SInt>( std::get<0>(globalCoords[i]), std::get<1>(globalCoords[i]), globalIds[i] ) );
            }

            std::stable_sort( forSort.begin(), forSort.end(), []( auto x1, auto x2){
                return std::get<2>(x1) < std::get<2>(x2);
                });
            
            std::cout<< __FILE__ << ", " << __LINE__ << ": before removing duplicates: " << forSort.size() << std::endl;
            forSort.erase( unique(forSort.begin(), forSort.end()), forSort.end() );
            std::cout<< __FILE__ << ", " << __LINE__ << ": after: " << forSort.size() << std::endl;

            //unzip
            SInt newSize = forSort.size();
            globalIds.resize( newSize, 0);
            globalCoords.resize( newSize,std::tuple<LPFloat, LPFloat>(0.0, 0.0) );
            for ( unsigned int i=0; i<newSize; i++ ){
              auto edge = forSort[i];
              globalCoords[i] = std::tuple<LPFloat, LPFloat>( std::get<0>(edge) , std::get<1>(edge) );
              globalIds[i] = std::get<2>(edge);
            }
        }

      for( unsigned int i=0; i<globalIds.size()-1; i++){
        if( globalIds[i]==globalIds[i+1]){
          throw std::logic_error("ERROR: vertex id " + std::to_string(globalIds[i]) + " is duplicated.");
        }
      }

      // Output edges
      FILE* fout = fopen(config_.coord_file.c_str(), "w+");

      for (auto coord : globalCoords){
        fprintf(fout, "%f %f\n", std::get<0>(coord) , std::get<1>(coord) );
      }
      fclose(fout);
    }
      
  }//GatherPrint, 2D

  //for the edges
  void GatherPrint(identity<std::tuple<SInt, SInt>>) const {
    // Exchange local dist
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Gather number of edges for each PE
    std::vector<int> displ(size);
    std::vector<int> num_edges(size);
    int lSize = NumEdges();
    MPI_Gather(&lSize, 1, MPI_INT,
               num_edges.data(), 1, MPI_INT, 
               ROOT, MPI_COMM_WORLD);
    int current_displ = 0;
    int total_num_edges = 0;
    if (rank == ROOT) {
      for (SInt i = 0; i < num_edges.size(); ++i) {
        displ[i] = current_displ;
        total_num_edges += num_edges[i];
        current_displ = total_num_edges;
      }
    }

    // Gather actual edges
    MPI_Datatype MPI_EDGE;
    MPI_Type_vector(1, 2, 0, MPI_LONG, &MPI_EDGE);
    MPI_Type_commit(&MPI_EDGE);
    std::vector<Edge> edges(total_num_edges);
    MPI_Gatherv(edges_.data(), lSize, MPI_EDGE,
                edges.data(), num_edges.data(), displ.data(), MPI_EDGE, 
                ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
      // Sort edges and remove duplicates
      std::sort(std::begin(edges), std::end(edges));
      //SInt total_edges = edges.size();
      edges.erase(unique(edges.begin(), edges.end()), edges.end());
      
      // Output edges
      FILE* fout = fopen(config_.output_file.c_str(), "w+");
#ifndef OMIT_HEADER
      fprintf(fout, "%llu %lu\n", config_.n, edges.size());
#endif
      for (auto edge : edges) fprintf(fout, "%llu %llu\n", std::get<0>(edge) , std::get<1>(edge) );
      fclose(fout);
    }
  }

  //
  // distributed versions
  //

  // node id output
  void Print(identity<std::tuple<SInt, SInt>>) const {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SInt num_edges = edges_.size();
    SInt total_num_edges = 0;
    MPI_Allreduce(&num_edges, &total_num_edges, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    FILE* fout =
        fopen((config_.output_file + "_" + std::to_string(rank)).c_str(), "w+");
#ifndef OMIT_HEADER
    fprintf(fout, "%% %llu %llu\n", config_.n, total_num_edges);
#endif
    for (auto edge : edges_) {
      fprintf(fout, "%llu %llu\n", std::get<0>(edge) , std::get<1>(edge) );
    }
    fclose(fout);
    if( rank==0)
      std::cout << "Finished writing in file "<< config_.output_file << std::endl;
  };

  // 3D geometric point output
  void Print(identity<std::tuple<LPFloat, LPFloat, LPFloat, SInt>>) const {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);
    for (PEID i = 0; i < size; i++) {
      if (i == rank) {
        std::string mode = i == 0 ? "w+" : "a";
        FILE* fout = fopen(config_.coord_file.c_str(), mode.c_str());
        for (auto edge : edges_) {
          fprintf(fout, "%f %f %f\n", std::get<0>(edge), std::get<1>(edge),
                  std::get<2>(edge));
        }
        fclose(fout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if( rank==0)
      std::cout << "Finished writing in file "<< config_.coord_file << std::endl;
  };

  // 2D geometric point output
  void Print(identity<std::tuple<LPFloat, LPFloat, SInt>>) const {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);
    for (PEID i = 0; i < size; i++) {
      if (i == rank) {
        std::string mode = i == 0 ? "w+" : "a";
        FILE* fout = fopen(config_.coord_file.c_str(), mode.c_str());
        for (auto edge : edges_) {
          fprintf(fout, "%f %f\n", std::get<0>(edge), std::get<1>(edge));
        }
        fclose(fout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if( rank==0)
      std::cout << "Finished writing in file "<< config_.coord_file << std::endl;
  };

  // ABUSE: adjacency list output
  void Print(identity<std::tuple<SInt, std::vector<SInt>>>) const {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // fugly but saves memory as this is a const method
    auto& nodes = const_cast<std::vector<Edge>&>(edges_);

    // sort edges by node id
    std::sort(std::begin(nodes), std::end(nodes),
              [](const Edge& t1, const Edge& t2) {
                return std::get<0>(t1) <
                       std::get<0>(t2);  // or use a custom compare function
              });

    // compute edge count
    SInt edgeCount = std::accumulate(
        std::begin(nodes), std::end(nodes), SInt(0),
        [](SInt a, const Edge& b) { return a + std::get<1>(b).size(); });

    FILE* fout =
        fopen((config_.output_file + std::to_string(rank)).c_str(), "w+");
    fprintf(fout, "%lu %llu\n", edges_.size(), edgeCount);

    for (auto& node : nodes) {
      auto& edges = std::get<1>(node);
      edgeCount += edges.size();
      std::sort(std::begin(edges), std::end(edges));

      std::string sep = "";
      for (const auto& edge : edges) {
        fprintf(fout, "%s%llu", sep.c_str(), edge);
        sep = " ";
      }
      fprintf(fout, "\n");
    }

    fclose(fout);
  };
};

}
#endif
