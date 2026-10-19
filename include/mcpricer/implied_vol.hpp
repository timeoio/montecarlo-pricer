#pragma once

#include <vector>

#include "mcpricer/option.hpp"

namespace mc {

struct ImpliedVolResult {
    double vol;
    int iterations;
    std::vector<double> history;  // sigma after each iteration, to show the convergence
};

// Starting point of Newton-Raphson: the sigma that maximises vega.
double initial_guess(const Option& option, const Market& market);

// Finds sigma such that bs_price(option, sigma) = target_price, by Newton-Raphson:
//   sigma_{n+1} = sigma_n - (bs_price(sigma_n) - target_price) / vega(sigma_n)
// market.vol is ignored. Throws if the price violates the no-arbitrage bounds.
ImpliedVolResult implied_vol(const Option& option, const Market& market, double target_price,
                             double tolerance = 1e-10, int max_iterations = 100);

}  // namespace mc
