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

#include "./insert.hpp"

namespace triqs_cthyb {

  histogram *move_insert_c_cdag::add_histo(std::string const &name, histo_map_t *histos, int nbins) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, config.beta(), nbins}});
    return &(new_histo.first->second);
  }

  move_insert_c_cdag::move_insert_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data,
                                         mc_tools::random_generator &rng, histo_map_t *histos, int nbins,
                                         std::vector<double> const *hist_insert, std::vector<double> const *hist_remove, 
					 std::vector<double> *wr_insert, std::vector<int> *count_insert)
     : data(data),
       config(data.config),
       rng(rng),
       block_index(block_index),
       block_size(block_size),
       count_insert(count_insert),
       histo_proposed(add_histo("insert_length_proposed_" + block_name, histos, nbins)),
       histo_accepted(add_histo("insert_length_accepted_" + block_name, histos, nbins)),
       hist_insert(hist_insert),
       hist_remove(hist_remove),
       meas_wr(wr_insert),
       step_d(config.beta() / double(nbins)),
       step_i(time_pt::Nmax / nbins),
       t1(time_pt(1, config.beta())),
       use_improved_sampling(hist_insert && hist_remove),
       wr_insert(wr_insert) {} 

  mc_weight_t move_insert_c_cdag::attempt() {

#ifdef EXT_DEBUG
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "In config " << config.get_id() << std::endl;
    std::cerr << "* Attempt for move_insert_c_cdag (block " << block_index << ")" << std::endl;
#endif

    // Pick up the value of alpha and choose the operators
    auto rs1 = rng(block_size), rs2 = rng(block_size);
    op1 = op_desc{block_index, rs1, true, data.linindex[std::make_pair(block_index, rs1)]};
    op2 = op_desc{block_index, rs2, false, data.linindex[std::make_pair(block_index, rs2)]};

    auto &det    = data.dets[block_index];
    int det_size = det.size();

    // Choice of times for insertion. Find the time as double and them put them on the grid.
    if (!use_improved_sampling) tau1 = data.tau_seg.get_random_pt(rng);
    tau2 = data.tau_seg.get_random_pt(rng);
    double fac = 1.;
    if (use_improved_sampling) {
      // first choose the bin, each bin being weighted by the probability hist_insert[bin]*length(bin)
      double ran  = double(rng(time_pt::Nmax)) / double(time_pt::Nmax - 1); // random number in [0,1]
      double csum = 0;
      int nbins = (*hist_insert).size();
      int ibin = 0;
      for (int i = 0; i < nbins; ++i) {
        csum += (*hist_insert)[i] * step_d;
        if (csum >= ran || i == (nbins - 1) ) {
          ibin = i;
          break;
        }
      }

      time_pt bin_start  = time_pt(step_i * ibin, config.beta());
      time_pt bin_finish = time_pt(step_i * (ibin+1), config.beta());
      if (ibin == nbins - 1) bin_finish = time_pt(time_pt::Nmax, config.beta());
      
      // now draw a time point uniformly within this bin
      tau1 = tau2 + data.tau_seg.get_random_pt(rng, bin_start, bin_finish);

      // compute the probability of proposing the current config from the trial one
      // this is simply hist_remove[bin(tau1-tau2)] / sum_i(hist_remove(bin(tau_i - tau2)))
      // where the sum is performed over all creation operators of the current block, including the trial one
      int ind;
      double s = (*hist_remove)[ibin]; // normalization constant = sum_i(hist_remove(bin(tau_i - tau2)))
      time_pt dtau_r; 

      for (int i = 0; i < det_size; ++i) {
        dtau_r = det.get_x(i).first - tau2;
        ind = floor_div(dtau_r, t1) / step_i;
        s += (*hist_remove)[ind];
      }

      // corrective factor for t_ratio
      if (std::abs(s) < 1.e-15 || (*hist_remove)[ibin] == 0.0) return 0; // quick return
      fac = double(data.dets[block_index].size() + 1) * (*hist_remove)[ibin] / (s * config.beta() * (*hist_insert)[ibin]);
    }

#ifdef EXT_DEBUG
    std::cerr << "* Proposing to insert:" << std::endl;
    std::cerr << op1 << " at " << tau1 << std::endl;
    std::cerr << op2 << " at " << tau2 << std::endl;
#endif

    // record the length of the proposed insertion
    dtau = double(tau2 - tau1);
    if (histo_proposed) *histo_proposed << dtau;

    // Insert the operators op1 and op2 at time tau1, tau2
    // 1- In the very exceptional case where the insert has failed because an operator is already sitting here
    // (cf std::map doc for insert return), we reject the move.
    // 2- If ok, we store the iterator to the inserted operators for later removal in reject if necessary
    try {
      data.imp_trace.try_insert(tau1, op1);
      data.imp_trace.try_insert(tau2, op2);
    } catch (rbt_insert_error const &) {
      std::cerr << "Insert error : recovering ... " << std::endl;
      data.imp_trace.cancel_insert();
      return 0;
    }

    // Computation of det ratio

    // Find the position for insertion in the determinant
    // NB : the determinant stores the C in decreasing time order.
    int num_c_dag, num_c;
    for (num_c_dag = 0; num_c_dag < det_size; ++num_c_dag) {
      if (det.get_x(num_c_dag).first < tau1) break;
    }
    for (num_c = 0; num_c < det_size; ++num_c) {
      if (det.get_y(num_c).first < tau2) break;
    }

    // Insert in the det. Returns the ratio of dets (Cf det_manip doc).
    auto det_ratio = det.try_insert(num_c_dag, num_c, {tau1, op1.inner_index}, {tau2, op2.inner_index});

    // proposition probability
    mc_weight_t t_ratio = std::pow(block_size * config.beta() / double(det.size() + 1), 2);
    if (use_improved_sampling) t_ratio *= fac;

    // For quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(t_ratio * det_ratio / data.atomic_weight);

    // computation of the new trace after insertion
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
      std::cerr << "atomic_weight == 0" << std::endl;
#endif
      return 0;
    }
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(insert) trace_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;

    if (meas_wr) { 
      int ibin = floor_div(tau1 - tau2, t1) / step_i; 
      (*wr_insert)[ibin] += std::abs(p);
      (*count_insert)[ibin] ++;
    }

#ifdef EXT_DEBUG
    std::cerr << "Atomic ratio: " << atomic_weight_ratio << '\t';
    std::cerr << "Det ratio: " << det_ratio << '\t';
    std::cerr << "Prefactor: " << t_ratio << '\t';
    std::cerr << "Weight: " << p * t_ratio << std::endl;
    std::cerr << "p_yee * newtrace: " << p_yee * new_atomic_weight << std::endl;
#endif

    if (!isfinite(p * t_ratio)) {
      std::cerr << "Insert move info:\n";
      std::cerr << "Atomic ratio: " << atomic_weight_ratio << '\t';
      std::cerr << "Det ratio: " << det_ratio << '\t';
      std::cerr << "Prefactor: " << t_ratio << '\t';
      std::cerr << "Weight: " << p * t_ratio << std::endl;
      std::cerr << "p_yee * newtrace: " << p_yee * new_atomic_weight << std::endl;
      
      TRIQS_RUNTIME_ERROR << "(insert) p * t_ratio not finite p : " << p << " t_ratio : " << t_ratio << " in config " << config.get_id();
    }
    return p * t_ratio;
  }

  mc_weight_t move_insert_c_cdag::accept() {
	  
    time_pt tau_min = std::min(tau1,tau2);
    time_pt tau_max = std::max(tau1,tau2);
    if (tau_min < data.imp_trace.min_tau) data.imp_trace.min_tau = tau_min;
    if (tau_max > data.imp_trace.max_tau) data.imp_trace.max_tau = tau_max;
    data.updated = true;

    // insert in the tree
    data.imp_trace.confirm_insert();

    // insert in the configuration
    config.insert(tau1, op1);
    config.insert(tau2, op2);
    config.finalize();

    // insert in the determinant
    data.dets[block_index].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    if (histo_accepted) *histo_accepted << dtau;

#ifdef EXT_DEBUG
    std::cerr << "* Move move_insert_c_cdag accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif

    return data.current_sign / data.old_sign;
  }

  void move_insert_c_cdag::reject() {

    config.finalize();
    data.imp_trace.cancel_insert();
    data.dets[block_index].reject_last_try();

#ifdef EXT_DEBUG
    std::cerr << "* Move move_insert_c_cdag rejected" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
  }
}
