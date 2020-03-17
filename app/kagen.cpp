/*******************************************************************************
 * app/generate_kagen.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
//#define DEL_STATS 1

#include <mpi.h>

#include "benchmark.h"
#include "generator_config.h"
#include "io/generator_io.h"
#include "parse_parameters.h"
#include "timer.h"

#include "geometric/delaunay/delaunay_2d.h"
#include "geometric/delaunay/delaunay_3d.h"
#include "geometric/rgg/rgg_2d.h"
#include "geometric/rgg/rgg_3d.h"
#include "gnm/gnm_directed.h"
#include "gnm/gnm_undirected.h"
#include "gnp/gnp_directed.h"
#include "gnp/gnp_undirected.h"
#include "hyperbolic/hyperbolic.h"
#include "barabassi/barabassi.h"
#include "kronecker/kronecker.h"
#include "grid/grid_2d.h"
#include "grid/grid_3d.h"

using namespace kagen;

void OutputParameters(PGeneratorConfig &config, const PEID /* rank */,
                      const PEID size, std::ostream& out=std::cout ) {
  if (config.generator == "gnm_directed" ||
      config.generator == "gnm_undirected" ||
      config.generator == "gnp_directed" ||
      config.generator == "gnp_undirected")
    out << "generate graph (n=" << config.n << ", m=" << config.m
              << " (p=" << config.p << "), k=" << config.k
              << ", s=" << config.seed << ", P=" << size << ")" << std::endl;

  else if (config.generator == "rgg_2d" || config.generator == "rgg_3d")
    out << "generate graph (n=" << config.n << ", r=" << config.r
              << ", k=" << config.k << ", s=" << config.seed << ", P=" << size
              << ")" << std::endl;

  else if (config.generator == "rdg_2d" || config.generator == "rdg_3d")
    out << "generate graph (n=" << config.n << ", k=" << config.k
              << ", s=" << config.seed << ", P=" << size << ")" << std::endl;

  else if (config.generator == "rhg")
    out << "generate graph (n=" << config.n << ", d=" << config.avg_degree
              << ", gamma=" << config.plexp << ", k=" << config.k
              << ", s=" << config.seed << ", P=" << size << ")" << std::endl;

  else if (config.generator == "ba")
    out << "generate graph (n=" << config.n << ", d=" << config.min_degree
              << ", k=" << config.k << ", s=" << config.seed << ", P=" << size
              << ")" << std::endl;

  else if (config.generator == "rmat")
    out << "generate graph (n=" << config.n << ", m=" << config.m
              << ", k=" << config.k << ", s=" << config.seed << ", P=" << size
              << ")" << std::endl;

  else if (config.generator == "grid_2d")
    out << "generate graph (row=" << config.grid_x << ", col=" << config.grid_y
              << ", p=" << config.p << ", k=" << config.k << ", s=" << config.seed 
              << ", P=" << size << ")" << std::endl;

  else if (config.generator == "grid_3d")
    out << "generate graph (x=" << config.grid_x << ", y=" << config.grid_y << ", z=" << config.grid_z
              << ", p=" << config.p << ", k=" << config.k << ", s=" << config.seed 
              << ", P=" << size << ")" << std::endl;
}

template <typename Generator, typename EdgeCallback>
void RunGenerator(PGeneratorConfig &config, const PEID rank,
                  const PEID size , Statistics &stats, Statistics &edge_stats,
                  Statistics &edges, const EdgeCallback &cb) {
  // Start timers
  Timer t;
  double local_time = 0.0;
  double total_time = 0.0;
  t.Restart();

  // Chunk distribution
  Generator gen(config, rank, cb);
  gen.Generate();

  // Output
  local_time = t.Elapsed();
  MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, ROOT,
             MPI_COMM_WORLD);

  if (rank == ROOT) {
    stats.Push(total_time);
    edge_stats.Push(total_time / gen.NumberOfEdges());
    edges.Push(gen.NumberOfEdges());
  }

  unsigned long int localNumEdges = gen.NumberOfEdges();
  unsigned long int globalNumEdges  = 0;
  MPI_Reduce( &localNumEdges, &globalNumEdges, 1, MPI_UNSIGNED_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);

  if (rank == ROOT){
    std::cout   << "done generating graph, time: " << total_time
                << ", total number of edges: " << globalNumEdges/2
                << ", write output..." << std::endl;
  }
  
  local_time = 0.0;
  total_time = 0.0;
  t.Restart();

  gen.Output();
  local_time = t.Elapsed();
  MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, ROOT,
             MPI_COMM_WORLD);

#ifdef OUTPUT_EDGES
  if (rank == ROOT and config.output_format=="binary"){
    //file to store some data about the output files
    std::string headerFile = config.output_file+".info";
    auto fheader = std::fstream(headerFile, std::fstream::app);
    OutputParameters(config, rank, size, fheader);
    fheader << "time to write output: "<< total_time << std::endl;
    std::cout << "time to write output: "<< total_time << std::endl;
  }
#endif
}

int main(int argn, char **argv) {
  // Init MPI
  MPI_Init(&argn, &argv);
  PEID rank=0, size=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Read command-line args
  PGeneratorConfig generator_config;
  ParseParameters( argn, argv, rank, size, generator_config, std::cout);

  std::string callingCommand = "";
  for (int i = 0; i < argn; i++) {
    callingCommand += std::string(argv[i]) + " ";
  }  

  if (rank == ROOT){
    OutputParameters(generator_config, rank, size);
    std::cout<<"%%Calling command: " << callingCommand <<std::endl;
    std::cout<<"%%Called with " << size << " MPI processes" << std::endl;
  }

  // Statistics
  Statistics stats;
  Statistics edge_stats;
  Statistics edges;

  auto edge_cb = [](SInt, SInt){};
  ULONG user_seed = generator_config.seed;
  for (ULONG i = 0; i < generator_config.iterations; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);
    generator_config.seed = user_seed + i;
    if (generator_config.generator == "gnm_directed")
      RunGenerator<GNMDirected<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "gnm_undirected")
      RunGenerator<GNMUndirected<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "gnp_directed")
      RunGenerator<GNPDirected<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "gnp_undirected")
      RunGenerator<GNPUndirected<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rgg_2d")
      RunGenerator<RGG2D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rgg_3d")
      RunGenerator<RGG3D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rdg_2d")
      RunGenerator<Delaunay2D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rdg_3d")
      RunGenerator<Delaunay3D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rhg")
      RunGenerator<Hyperbolic<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "ba")
      RunGenerator<Barabassi<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "rmat")
      RunGenerator<Kronecker<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "grid_2d")
      RunGenerator<Grid2D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else if (generator_config.generator == "grid_3d")
      RunGenerator<Grid3D<decltype(edge_cb)>, decltype(edge_cb)>
        (generator_config, rank, size, stats, edge_stats, edges, edge_cb);
    else 
      if (rank == ROOT) std::cout << "generator not supported" << std::endl;
  }

  if (rank == ROOT) {
    std::cout << "RESULT runner=" << generator_config.generator
              << " time=" << stats.Avg() << " stddev=" << stats.Stddev()
              << " iterations=" << generator_config.iterations
              << " edges=" << edges.Avg()
              << " time_per_edge=" << edge_stats.Avg()
              << " r= " << generator_config.r
              << " output stored in file " << generator_config.output_file << std::endl;
  }
  //std::cout << rank << ": edges=" << edges.Avg() << std::endl;

  MPI_Finalize();
  return 0;
}
