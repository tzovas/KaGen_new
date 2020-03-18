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
#include "timer.h"
#include "mem_monitor.h"

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

  void OutputEdges( bool binary=false) const { 
#ifdef SINGLE_LIST
    GatherPrint(identity<Edge>(), binary);
#else
    Print(identity<Edge>()); 
#endif
  }

  SInt NumEdges() const { 
    return edges_.size() > 0 ? edges_.size() : local_num_edges_/2; 
  }
  
  void removeDuplicates(){
    PEID rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     // Sort edges and remove duplicates locally
    std::cout<< __FILE__ << ", " << __LINE__ << ", " << rank << ": local num edges before removing duplicates: " << edges_.size() << std::endl;
    std::sort(std::begin(edges_), std::end(edges_));
    edges_.erase(unique(edges_.begin(), edges_.end()), edges_.end());
    std::vector<Edge>(edges_).swap(edges_); //shrink to fit
    std::cout<< __FILE__ << ", " << __LINE__ << ", " << rank << ": after: " << edges_.size() << std::endl;

  }

 private:
  PGeneratorConfig &config_;

  std::vector<SInt> dist_;
  std::vector<Edge> edges_;

  SInt local_num_edges_;

//---------------------------
// 2D
//
#ifdef DIM_2D

  void GatherPrint(identity<std::tuple<LPFloat, LPFloat, LPFloat, SInt>>, [[maybe_unused]] bool binary=false) const {
    throw std::logic_error("Not implemented, should not be called. For the 3D generators recompile with KAGEN_DIMENSION_2D=OFF");
  }

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
  void GatherPrint(identity<std::tuple<LPFloat, LPFloat, SInt>>, [[maybe_unused]] bool binary=false) const {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<std::tuple<LPFloat, LPFloat>> globalCoords; 
    std::vector<SInt> globalIds; 

    Timer t;
    double local_time = 0.0;
    double total_time = 0.0;
    t.Restart();

    GatherCoords( identity<Edge>(), globalCoords, globalIds);

    local_time = t.Elapsed();
    MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, ROOT,
             MPI_COMM_WORLD);

    if (rank == ROOT) {

      std::cout << "time to gather the coordinates: " << total_time << std::endl;

      local_time = 0.0;
      t.Restart();
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

      //sanity check
      for( unsigned int i=0; i<globalIds.size()-1; i++){
        if( globalIds[i]==globalIds[i+1]){
          throw std::logic_error("ERROR: vertex id " + std::to_string(globalIds[i]) + " is duplicated.");
        }
        //do not allow isolated vertices
        //if( not (globalIds[i+1]==globalIds[i]+1) ){
        //  throw std::logic_error("ERROR: vertex id " + std::to_string(globalIds[i]+1) + " is missing.");
        //}
      }
      local_time = t.Elapsed();
      std::cout << "time to sort and check: " << local_time << std::endl;

      if (not binary){
        // Output edges     
        FILE* fout = fopen(config_.coord_file.c_str(), "w+");
        for (auto coord : globalCoords){
          fprintf(fout, "%f %f\n", std::get<0>(coord) , std::get<1>(coord) );
        }
        fclose(fout);
      }else{
        //for binary files we enforce the header
        auto fout = std::fstream(config_.coord_file, std::ios::out | std::ios::binary);

        const auto sizeOfCoord = sizeof std::get<0>(globalCoords[0]);
        const LPFloat zero = 0.0;
        for (auto coord : globalCoords){
          fout.write( reinterpret_cast<const char*>(&std::get<0>(coord)), sizeOfCoord );
          fout.write( reinterpret_cast<const char*>(&std::get<1>(coord)), sizeOfCoord );
          //TODO: geographer expects 3 dimensional coordinates for the binary reader
          fout.write( reinterpret_cast<const char*>(&zero), sizeOfCoord );
        }
        fout.close();
      }
    }//if (rank == ROOT) 

  }//GatherPrint, 2D


#elif DIM_3D

  void GatherPrint(identity<std::tuple<LPFloat, LPFloat, SInt>>, [[maybe_unused]] bool binary=false) const {
    throw std::logic_error("Not implemented, should not be called. For the 2D generators recompile with KAGEN_DIMENSION_2D=ON");
  }
//---------------------------
// 3D
//

  int GatherCoords( identity<std::tuple<LPFloat, LPFloat, LPFloat, SInt>>, std::vector<std::tuple<LPFloat, LPFloat, LPFloat>> &globalCoords, std::vector<SInt> &globalIds) const {
    // Exchange local dist
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //cannot remove duplicates locally as function is const  and cannot change edges_

    //break vector of tuples into two vectors
    unsigned int lSize = NumEdges();
    std::vector<std::tuple<LPFloat, LPFloat, LPFloat>> localCoords(lSize);
    std::vector<SInt> localIds(lSize);
    assert( lSize==edges_.size() );

    for ( unsigned int i=0; i<lSize; i++ ){
        auto edge = edges_[i];
        localCoords[i] = std::tuple<LPFloat, LPFloat, LPFloat>( std::get<0>(edge) , std::get<1>(edge), std::get<2>(edge) );
        localIds[i] = std::get<3>(edge);
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

    globalCoords.resize(total_num_edges);
    globalIds.resize(total_num_edges);

    // Gather actual coordinates
    MPI_Datatype MPI_COORD;
    MPI_Type_vector(1, 3, 0, MPI_DOUBLE, &MPI_COORD);
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


//3D point output
  void GatherPrint(identity<std::tuple<LPFloat, LPFloat, LPFloat, SInt>>, [[maybe_unused]] bool binary=false) const {
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<std::tuple<LPFloat, LPFloat, LPFloat>> globalCoords; 
    std::vector<SInt> globalIds; 

    Timer t;
    double local_time = 0.0;
    double total_time = 0.0;
    t.Restart();

    GatherCoords( identity<Edge>(), globalCoords, globalIds);

    local_time = t.Elapsed();
    MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, ROOT,
             MPI_COMM_WORLD);

    if (rank == ROOT) {

      std::cout << "time to gather the coordinates: " << total_time << std::endl;

      local_time = 0.0;
      t.Restart();
      //checks
      if( not std::is_sorted(globalIds.begin(), globalIds.end()) ){
        //throw std::logic_error("ERROR: For some reason, the coordinate indices are not sorted.");
        std::cout<< "WARNING: indices are not sorted, will sort now" << std::endl;

        //zip again to sort and remove duplicates
        std::vector<std::tuple<LPFloat, LPFloat, LPFloat, SInt>> forSort;
        for( unsigned long i=0; i<globalCoords.size(); i++){
          forSort.push_back( std::tuple<LPFloat, LPFloat, LPFloat, SInt>( 
            std::get<0>(globalCoords[i]), std::get<1>(globalCoords[i]), std::get<2>(globalCoords[i]), globalIds[i] ) );
        }

        std::stable_sort( forSort.begin(), forSort.end(), []( auto x1, auto x2){
          return std::get<3>(x1) < std::get<3>(x2);
        });
    
        std::cout<< __FILE__ << ", " << __LINE__ << ": before removing duplicates: " << forSort.size() << std::endl;
        forSort.erase( unique(forSort.begin(), forSort.end()), forSort.end() );
        std::cout<< __FILE__ << ", " << __LINE__ << ": after: " << forSort.size() << std::endl;

        //unzip
        SInt newSize = forSort.size();
        globalIds.resize( newSize, 0);
        globalCoords.resize( newSize,std::tuple<LPFloat, LPFloat, LPFloat>(0.0, 0.0, 0.0) );
        for ( unsigned int i=0; i<newSize; i++ ){
          auto edge = forSort[i];
          globalCoords[i] = std::tuple<LPFloat, LPFloat, LPFloat>( std::get<0>(edge), std::get<1>(edge), std::get<2>(edge) );
          globalIds[i] = std::get<3>(edge);
        }
      }

      //sanity check
      for( unsigned int i=0; i<globalIds.size()-1; i++){
        if( globalIds[i]==globalIds[i+1]){
          throw std::logic_error("ERROR: vertex id " + std::to_string(globalIds[i]) + " is duplicated.");
        }
      }
      local_time = t.Elapsed();
      std::cout << "time to sort and check: " << local_time << std::endl;

      if (not binary){
        // Output edges     
        FILE* fout = fopen(config_.coord_file.c_str(), "w+");
        for (auto coord : globalCoords){
          fprintf(fout, "%f %f %f\n", std::get<0>(coord), std::get<1>(coord), std::get<2>(coord) );
        }
        fclose(fout);
      }else{
        //for binary files we enforce the header
        auto fout = std::fstream(config_.coord_file, std::ios::out | std::ios::binary);

        const auto sizeOfCoord = sizeof std::get<0>(globalCoords[0]);
        for (auto coord : globalCoords){
          fout.write( reinterpret_cast<const char*>(&std::get<0>(coord)), sizeOfCoord );
          fout.write( reinterpret_cast<const char*>(&std::get<1>(coord)), sizeOfCoord );
          //TODO: geographer expects 3 dimensional coordinates for the binary reader
          fout.write( reinterpret_cast<const char*>(&std::get<2>(coord)), sizeOfCoord );
        }
        fout.close();
      }
    }//if (rank == ROOT) 

  }//GatherPrint, 2D

#else
  throw std::runtime_error(" must specify one dimension using the KAGEN_DIMENSION_2D at compile time");
#endif

  //TODO: consider making it non-const. It would save some time
  //for the edges
  void GatherPrint(identity<std::tuple<SInt, SInt>>, bool binary=false) const {
    // Exchange local dist
    PEID rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Gather number of edges for each PE
    std::vector<int> displ(size);
    std::vector<int> num_edges(size);
    SInt lSize = NumEdges();
    MPI_Gather(&lSize, 1, MPI_INT,
               num_edges.data(), 1, MPI_INT,
               ROOT, MPI_COMM_WORLD);
    SInt current_displ = 0;
    SInt total_num_edges = 0;
    if (rank == ROOT) {
      for (SInt i = 0; i < num_edges.size(); ++i) {
        displ[i] = current_displ;
        total_num_edges += num_edges[i];
        current_displ = total_num_edges;
      }
    }

    if (rank == ROOT){
        std::vector<Edge> edges_tmp;
        std::cout<< "rank "<< rank <<": about to request a vector of size " <<total_num_edges << std::endl;
        std::cout<< "possible ERROR: requested size is larger than the maximum allowed vector size: " << edges_tmp.max_size() << std::endl;
        std::cout<< "the vector will need at least " << total_num_edges*((double) sizeof(Edge))/(1024*1024) << " MBs of memory" << std::endl;
    }

    // Gather actual edges
    MPI_Datatype MPI_EDGE;
    MPI_Type_vector(1, 2, 0, MPI_LONG, &MPI_EDGE);
    MPI_Type_commit(&MPI_EDGE);
    std::vector<Edge> edges(total_num_edges); //this is needed because method is const

    if (rank == ROOT){
        std::cout<< "### memory before gather " << std::endl;
    }

    //memory consumption
    printMemUsage();

    MPI_Gatherv(edges_.data(), lSize, MPI_EDGE,
                edges.data(), num_edges.data(), displ.data(), MPI_EDGE, 
                ROOT, MPI_COMM_WORLD);

    //memory consumption
    printMemUsage();

    if (rank == ROOT) {

      //add all edges (i, i+1) to force graph to be connected
      for( unsigned int i=0; i<config_.n-1; i++){
        edges.push_back( std::tuple<SInt, SInt>(i, i+1) );
      }

      // Sort edges and remove duplicates
      std::cout<< __FILE__ << ", " << __LINE__ << ": num edges before removing duplicates: " << edges.size() << std::endl;        
      std::sort(std::begin(edges), std::end(edges));
      edges.erase(unique(edges.begin(), edges.end()), edges.end());
      std::vector<Edge>(edges).swap(edges); //shrink to fit
      std::cout<< __FILE__ << ", " << __LINE__ << ": after: " << edges.size() << std::endl;
      
      // Output edges
      if ( not binary){
        FILE* fout = fopen(config_.output_file.c_str(), "w+");
#ifndef OMIT_HEADER
        fprintf(fout, "%llu %lu\n", config_.n, edges.size());
#endif
        for (auto edge : edges) fprintf(fout, "%llu %llu\n", std::get<0>(edge) , std::get<1>(edge) );
        fclose(fout);
      }else{
        //for binary files we enforce the header
        const SInt edgesSize = edges.size();
        std::cout<< __FILE__ << ", " << __LINE__ << ": n= " << config_.n << ", m= " << edgesSize << std::endl;
        {
          //write header also in a separate text file
          std::string headerFile = config_.output_file+".info";
          auto fheader = std::fstream(headerFile, std::fstream::app);
          fheader <<  "n= " << config_.n << ", m= " << edgesSize << ", r= "<< config_.r<< std::endl;
        }
        auto fout = std::fstream(config_.output_file, std::ios::out | std::ios::binary);
        
        fout.write( reinterpret_cast<const char*>(&config_.n), sizeof config_.n );
        fout.write( reinterpret_cast<const char*>(&edgesSize), sizeof edgesSize );
        for (auto edge : edges){
          fout.write( reinterpret_cast<const char*>(&std::get<0>(edge)), sizeof std::get<0>(edge) );
          fout.write( reinterpret_cast<const char*>(&std::get<1>(edge)), sizeof std::get<1>(edge) );
        }
        fout.close();
      }
    }// if (rank == ROOT) 

  }// GatherPrint

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
//std::cout << rank << ": " << edges_.size() << std::endl;
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

};//class GeneratorIO 

}
#endif
