/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "./density_matrix.hpp"
#include <mpi/vector.hpp>

#include <iomanip>

namespace triqs_cthyb {

  measure_density_matrix::measure_density_matrix(qmc_data const &data, std::vector<matrix_t> &density_matrix) : data(data), block_dm(density_matrix) {
    block_dm.resize(data.imp_trace.get_density_matrix().size());
    for (int i = 0; i < block_dm.size(); ++i) {
      block_dm[i]   = data.imp_trace.get_density_matrix()[i].mat;
      block_dm[i]() = 0;
    }
  }
  // --------------------

  void measure_density_matrix::accumulate(mc_weight_t s) {
    if (data.updated) {
      if (!flag) {
        z += data.n_acc * old_z;
        mc_weight_t s_temp = data.n_acc * old_s;
        int size = block_dm.size();
        for (int i = 0; i < size; ++i) 
          if (data.imp_trace.get_density_matrix()[i].is_valid) { block_dm[i] += s_temp * data.imp_trace.get_density_matrix()[i].mat; }
      }	
      data.imp_trace.compute(-1,0,true);
      old_z = s * data.atomic_reweighting;
      old_s = s / data.atomic_weight;
    }
    flag = false;
  }

  // ---------------------------------------------

  void measure_density_matrix::collect_results(mpi::communicator const &c) {
    
    z += data.n_acc * old_z;
    mc_weight_t s_temp = data.n_acc * old_s;
    int size = block_dm.size();
    for (int i = 0; i < size; ++i)
      if (data.imp_trace.get_density_matrix()[i].is_valid) { block_dm[i] += s_temp * data.imp_trace.get_density_matrix()[i].mat; }
  
    z                          = mpi::all_reduce(z, c);
    block_dm                   = mpi::all_reduce(block_dm, c);
    for (auto &b : block_dm){
        // Normalize
        b /= real(z);
        
        // Enforce hermiticity
        b = make_regular(0.5*(b + dagger(b)));
    }

    if (c.rank() != 0) return;

    // Check: the trace of the density matrix must be 1 by construction
    h_scalar_t tr = 0;
    for (auto &b : block_dm) tr += trace(b);
    if (std::abs(tr - 1) > 1.e-13)
      std::cerr << "Warning :: Trace of the density matrix is " << std::setprecision(13) << tr << std::setprecision(6) << " instead of 1"
                << std::endl;

  }
}
