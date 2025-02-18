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

#include "./mirror.hpp"
#include <triqs/mc_tools.hpp>

namespace triqs_cthyb {

  move_mirror::move_mirror(qmc_data &data, mc_tools::random_generator &rng)
     : data(data),
       config(data.config),
       rng(rng) {}

  mc_weight_t move_mirror::attempt() {

    if (config.size() == 0) return 0;  // Nothing to exchange !

    // Build updated ops
    updated_ops.clear();
    auto it = config.end();
    for (auto &[tau, op] : config) {
      --it;
      auto new_tau  = -it->first;
      auto new_op   = it->second;
      new_op.dagger = !new_op.dagger;
      updated_ops[tau] = std::make_pair(new_tau,new_op);
    }

    // --- Modify the tree
    data.imp_trace.try_mirror(updated_ops);

    mc_weight_t det_ratio = 1.;
    for (int block_index : range(data.delta.size())) det_ratio *= data.dets[block_index].try_mirror();

    // for quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(det_ratio / data.atomic_weight);

    // --- Compute the atomic_weight ratio
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);

    if (new_atomic_weight == 0.0) return 0;
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

    // --- Compute the weight
    mc_weight_t p = atomic_weight_ratio * det_ratio;

    return p;
  }

  mc_weight_t move_mirror::accept() {

    data.imp_trace.min_tau = time_pt(0,config.beta());
    data.imp_trace.max_tau = time_pt(time_pt::Nmax,config.beta());
    data.updated = true;

    // Update the tree
    data.imp_trace.confirm_mirror();

    for (int block_index : range(data.delta.size())) data.dets[block_index].complete_operation();

    // Update the configuration
    config.mirror();
    config.finalize();

    // Update the determinant
    data.update_sign();

    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;

    auto result = data.current_sign / data.old_sign;

    return result;
  }

  void move_mirror::reject() {

    config.finalize();
    for (int block_index : range(data.delta.size())) data.dets[block_index].reject_last_try();
    data.imp_trace.cancel_mirror();

  }
}
