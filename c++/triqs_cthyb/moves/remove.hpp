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
#pragma once
#include <algorithm>
#include <triqs/mc_tools.hpp>
#include "../qmc_data.hpp"

namespace triqs_cthyb {

  // Removal of C, C^dagger operator
  class move_remove_c_cdag {

    qmc_data &data;
    configuration &config;
    mc_tools::random_generator &rng;
    int block_index, block_size;
    std::vector<int> *count_remove;
    histogram *histo_proposed, *histo_accepted; // Analysis histograms
    double dtau;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;
    time_pt tau1, tau2;
    const std::vector<double> *hist_insert, *hist_remove;
    bool meas_wr;
    int Nmax;
    const uint64_t step_i;
    const time_pt t1;
    bool use_improved_sampling;
    std::vector<int> bins;
    std::vector<double> *wr_remove;

    histogram *add_histo(std::string const &name, histo_map_t *histos, int nbins);

    public:
    move_remove_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data, mc_tools::random_generator &rng,
                       histo_map_t *histos, int nbins, std::vector<double> const *hist_insert, std::vector<double> const *hist_remove,
		       std::vector<double> *wr_remove, std::vector<int> *count_remove);

    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };
}
