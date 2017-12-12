/******************************************************************************
 * t
 * barabassi.h
 *
 * Source of the graph generator
 ******************************************************************************
 * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef _BARABASSI_H_
#define _BARABASSI_H_

#include <iostream>
#include <vector>

#include "definitions.h"
#include "generator_config.h"
#include "generator_io.h"
#include "tools/spooky_hash.h"

class Barabassi {
 public:
  Barabassi(const PGeneratorConfig &config, const PEID rank)
      : config_(config),
        rank_(rank),
        io_(config),
        min_degree_(config_.min_degree),
        total_degree_(2 * config_.min_degree) {
    PEID size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Init variables
    from_ = rank * ceil(config_.n / (LPFloat)size);
    to_ = std::min((SInt)((rank + 1) * ceil(config_.n / (LPFloat)size) - 1),
                   config_.n - 1);
  }

  void Generate() {
    GenerateEdges();

    // Additional stats
    if (rank_ == ROOT) {
      HPFloat terra_bytes = 128 * config_.n * min_degree_;
      terra_bytes /= (8);
      terra_bytes /= (1024);
      terra_bytes /= (1024);
      terra_bytes /= (1024);
      terra_bytes /= (1024);

      std::cout << "memory of graph in tera bytes  " << std::setprecision(40)
                << terra_bytes << std::endl;
    }
  }

  void Output() const {
#ifdef OUTPUT_EDGES
    io_.OutputEdges();
#else
    io_.OutputDist();
#endif
  }

  SInt NumberOfEdges() const { return io_.NumEdges(); }

 private:
  // Config
  PGeneratorConfig config_;
  PEID rank_;

  // I/O
  GeneratorIO<> io_;

  // Constants and variables
  SInt min_degree_;
  SInt total_degree_;
  SInt from_, to_;

  void GenerateEdges() {
    for (SInt v = from_; v <= to_; v++) {
      for (SInt i = 0; i < min_degree_; i++) {
        SInt r = 2 * (v * min_degree_ + i) + 1;
        do {
          // compute hash h(r)
          SInt hash = Spooky::Hash(r);
          r = hash % r;
        } while (r % 2 == 1);
        SInt w = r / total_degree_;
#ifdef OUTPUT_EDGES
        io_.PushEdge(v, w);
        io_.PushEdge(w, v);
#else
        io_.UpdateDist(v);
        io_.UpdateDist(w);
#endif
      }
    }
  };
};
#endif
