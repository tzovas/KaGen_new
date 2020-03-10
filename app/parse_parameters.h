/*******************************************************************************
 * app/parse_parameters.cpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _PARSE_PARAMETERS_H_
#define _PARSE_PARAMETERS_H_

#include <string.h>
#include <cmath>

#include "generator_config.h"
#include "tools/arg_parser.h"

#include "definitions.h"

namespace kagen {

void ParseParameters(int argn, char **argv,
                     PEID rank , PEID size,
                     PGeneratorConfig &generator_config,
                     std::ostream& out=std::cout ) {
  ArgParser args(argn, argv);

  // Generator
  generator_config.generator = args.Get<std::string>("gen", "");

  if ( (args.IsSet("help") || argn < 2) && rank==ROOT ){
    if (generator_config.generator == "") {
      out << "%%================================================" << std::endl;
      out << "%%==================== KaGen =====================" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%Usage:\t\t\tmpirun -n <num_proc> ./kagen -gen <generator> [additional parameters]" << std::endl;
      out << "%%Generators:\t\tgnm_directed|gnm_undirected|gnp_directed|gnp_undirected|rgg_2d|rgg_3d|rdg_2d|rdg_3d|ba|rhg" << std::endl;
      out << "%%Additional help:\t./kagen -gen <generator> -help" << std::endl;
    }
    
    if (generator_config.generator == "gnm_undirected" || generator_config.generator == "gnm_directed") {
      out << "%%================================================" << std::endl;
      out << "%%========== Erdos-Renyi Graphs G(n,m) ===========" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%Parameters:" << std::endl;
      out << "%%-n\t\t<number of vertices as a power of two>" << std::endl;
      out << "%%-m\t\t<number of edges as a power of two>" << std::endl;
      out << "%%-k\t\t<number of chunks>" << std::endl;
      out << "%%-seed\t\t<seed for PRNGs>" << std::endl;
      out << "%%-output\t\t<output file>" << std::endl;
      out << "%%-self_loops" << std::endl;
      out << "%%\nExample:" << std::endl;
      out << "%%mpirun -n 16 ./build/app/kagen -gen gnm_directed -n 20 -m 22 -self_loops -output tmp" << std::endl;
    } else if (generator_config.generator == "gnp_undirected" || generator_config.generator == "gnp_directed") {
      out << "%%================================================" << std::endl;
      out << "%%========== Erdos-Renyi Graphs G(n,p) ===========" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%Parameters:" << std::endl;
      out << "%%-n\t\t<number of vertices as a power of two>" << std::endl;
      out << "%%-p\t\t<edge probability>" << std::endl;
      out << "%%-k\t\t<number of chunks>" << std::endl;
      out << "%%-seed\t\t<seed for PRNGs>" << std::endl;
      out << "%%-output\t\t<output file>" << std::endl;
      out << "%%-self_loops" << std::endl;
      out << "%%\nExample:" << std::endl;
      out << "%%mpirun -n 16 ./build/app/kagen -gen gnp_directed -n 20 -p 0.001 -self_loops -output tmp" << std::endl;
    } else if (generator_config.generator == "rgg_2d" || generator_config.generator == "rgg_3d") {
      out << "%%================================================" << std::endl;
      out << "%%======= Random Geometric Graphs RGG(n,d) ========" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%Parameters:" << std::endl;
      out << "%%-n\t\t<number of vertices as a power of two>" << std::endl;
      out << "%%-r\t\t<radius for vertices to be connected> (r <= 1.0)>" << std::endl;
      out << "%%-k\t\t<number of chunks>" << std::endl;
      out << "%%-seed\t\t<seed for PRNGs>" << std::endl;
      out << "%%-output\t\t<output file>" << std::endl;
      out << "%%\nExample:" << std::endl;
      out << "%%mpirun -n 16 ./build/app/kagen -gen rgg_3d -n 20 -r 0.00275 -output tmp" << std::endl;
    } else if (generator_config.generator == "rdg_2d" || generator_config.generator == "rdg_3d") {
      out << "%%================================================" << std::endl;
      out << "%%======== Random Delaunay Graphs RDG(n) =========" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%Parameters:" << std::endl;
      out << "%%-n\t\t<number of vertices as a power of two>" << std::endl;
      out << "%%-k\t\t<number of chunks>" << std::endl;
      out << "%%-seed\t\t<seed for PRNGs>" << std::endl;
      out << "%%-output\t\t<output file>" << std::endl;
      out << "%%\nExample:" << std::endl;
      out << "%%mpirun -n 16 ./build/app/kagen -gen rdg_3d -n 20 -output tmp" << std::endl;
    } else if (generator_config.generator == "ba") {
      out << "%%================================================" << std::endl;
      out << "%%======= Barabassi-Albert Graphs BA(n,d) ========" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%Parameters:" << std::endl;
      out << "%%-n\t\t<number of vertices as a power of two>" << std::endl;
      out << "%%-md\t\t<min degree for each vertex>" << std::endl;
      out << "%%-k\t\t<number of chunks>" << std::endl;
      out << "%%-seed\t\t<seed for PRNGs>" << std::endl;
      out << "%%-output\t\t<output file>" << std::endl;
      out << "%%\nExample:" << std::endl;
      out << "%%mpirun -n 16 ./build/app/kagen -gen ba -n 20 -md 4 -output tmp" << std::endl;
    } else if (generator_config.generator == "rhg") {
      out << "%%Parameters for Random Hyperbolic Graphs RHG(n,gamma,d)" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%=== Random Hyperbolic Graphs RHG(n,gamma,d) ====" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%Parameters:" << std::endl;
      out << "%%-n\t\t<number of vertices as a power of two>" << std::endl;
      out << "%%-gamma\t\t<power-law exponent>" << std::endl;
      out << "%%-d\t\t<average degree>" << std::endl;
      out << "%%-k\t\t<number of chunks>" << std::endl;
      out << "%%-seed\t\t<seed for PRNGs>" << std::endl;
      out << "%%-output\t\t<output file>" << std::endl;
      out << "%%\nExample:" << std::endl;
      out << "%%mpirun -n 16 ./build/app/kagen -gen rhg -n 20 -d 8 -gamma 2.2 -output tmp" << std::endl;
    } else if (generator_config.generator == "rmat") {
      out << "%%Parameters for Kronecker Graphs RMAT(n,m)" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%=========== Kronecker Graphs RMAT(n,m) =========" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%Parameters:" << std::endl;
      out << "%%-n\t\t<number of vertices as a power of two>" << std::endl;
      out << "%%-m\t\t<number of edges as a power of two>" << std::endl;
      out << "%%-k\t\t<number of chunks>" << std::endl;
      out << "%%-seed\t\t<seed for PRNGs>" << std::endl;
      out << "%%-output\t\t<output file>" << std::endl;
      out << "%%\nExample:" << std::endl;
      out << "%%mpirun -n 16 ./build/app/kagen -gen rmat -n 20 -m 22 -output tmp" << std::endl;
    } else if (generator_config.generator == "grid") {
      out << "%%Parameters for 2D/3D Grid Graphs G(x,y(,z),periodic)" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%=========== Grid Graphs G(x,y(,z)) ================" << std::endl;
      out << "%%================================================" << std::endl;
      out << "%%Parameters:" << std::endl;
      out << "%%-x\t\t<size of first dimension>" << std::endl;
      out << "%%-y\t\t<size of second dimension>" << std::endl;
      out << "%%-z\t\t<size of third dimension>" << std::endl;
      out << "%%-p\t\t<probability of edge insertion>" << std::endl;
      out << "%%-periodic\t\t<use periodic boundary condition>" << std::endl;
      out << "%%-k\t\t<number of chunks>" << std::endl;
      out << "%%-seed\t\t<seed for PRNGs>" << std::endl;
      out << "%%-output\t\t<output file>" << std::endl;
      out << "%%\nExample:" << std::endl;
      out << "%%mpirun -n 16 ./build/app/kagen -gen grid -x 16 -y 16 -output tmp" << std::endl;
    }
    exit(0);
  }

  // Nodes
  bool exact_n = args.IsSet("exact_n");
  if (exact_n)
    generator_config.n = args.Get<ULONG>("n", 100);
  else
    generator_config.n = (ULONG)1 << args.Get<ULONG>("n", 3);

  // Blocks
  generator_config.k = args.Get<ULONG>("k", size);

  // RNG
  generator_config.seed = args.Get<ULONG>("seed", 1);
  generator_config.hash_sample = args.Get<bool>("hash_sample", false);
  generator_config.use_binom = args.IsSet("binom");

  // I/O
  generator_config.output_file = args.Get<std::string>("output", "out");
  generator_config.coord_file = args.Get<std::string>("coord_output", "coords");
  generator_config.output_dir = args.Get<std::string>("output_dir", "./");
  generator_config.output_format = args.Get<std::string>("output_format", "text");

  // if no output parameter is given, create file name automatically
  bool isOutfileGiven = args.IsSet("output");
  if( not isOutfileGiven){
    std::string fileExt=".edgl";

    if( generator_config.output_format == "binary"){
        fileExt = ".bgf";
    }
    generator_config.output_file = generator_config.generator + "_" + std::to_string(args.Get<ULONG>("n", 3))+fileExt;
    generator_config.coord_file = generator_config.output_file+".xyz";
  }
  
  bool givenDir = args.IsSet("output_dir");
  if( givenDir){
    generator_config.output_file = generator_config.output_dir + "/" + generator_config.output_file;
    generator_config.coord_file = generator_config.output_dir+ "/" + generator_config.coord_file;
  }

  generator_config.debug_output = args.Get<std::string>("debug", "dbg");
  generator_config.dist_size = args.Get<ULONG>("dist", 10);

  // Edges
  bool exact_m = args.IsSet("exact_m");
  if (exact_m)
    generator_config.m = args.Get<ULONG>("m", 0);
  else
    generator_config.m = (ULONG)1 << args.Get<ULONG>("m", 0);
  generator_config.p = args.Get<double>("p", 0.0);
  generator_config.self_loops = args.IsSet("self_loops");

  // Radius/Edges
  generator_config.r = args.Get<double>("r", 0.125);

  // Average degree
  generator_config.avg_degree = args.Get<double>("d", 5.0);
  generator_config.plexp = args.Get<double>("gamma", 2.6);

  if(  args.IsSet("d") ){
    double newR = generator_config.r;
    if(generator_config.generator == "rgg_2d"){
       newR = std::sqrt( (double) 2*generator_config.avg_degree /(3.1415*generator_config.n) );
    }else if (generator_config.generator == "rgg_3d"){
        newR = std::pow( (double) 3*generator_config.avg_degree /(2*3.1415*generator_config.n), 3);
    }

    if(rank==ROOT and generator_config.r!=newR){
      out << "%% Setting/overwriting r to " << generator_config.r << std::endl;
    }
    generator_config.r = newR;
  }

  // RHG 
  generator_config.thres = args.Get<ULONG>("t", 0);
  generator_config.query_both = args.Get<bool>("qb", false);

  // BA
  generator_config.min_degree = args.Get<ULONG>("md", 4);

  // GRID
  generator_config.grid_x = args.Get<ULONG>("x", 1);
  generator_config.grid_y = args.Get<ULONG>("y", 1);
  generator_config.grid_z = args.Get<ULONG>("z", 1);
  generator_config.periodic = args.IsSet("periodic");

  // Floating-point precision
  generator_config.precision = args.Get<ULONG>("prec", 32);

  // Sampling algorithm
  generator_config.base_size = (ULONG)1 << args.Get<ULONG>("sk", 8);
  generator_config.hyp_base = (ULONG)1 << args.Get<ULONG>("hk", 8);

  // Benchmarks
  generator_config.iterations = args.Get<ULONG>("i", 1);
}//void ParseParameters

}
#endif
