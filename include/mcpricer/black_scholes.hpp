#pragma once

#include "mcpricer/option.hpp"

// Closed-form Black-Scholes formulas: the benchmark for Monte Carlo, and the
// function inverted by Newton-Raphson to get the implied volatility.
namespace mc {

double norm_pdf(double x);
double norm_cdf(double x);

double bs_price(const Option& option, const Market& market);

// dPrice/dSigma, identical for calls and puts. It is the derivative Newton-Raphson needs.
double bs_vega(const Option& option, const Market& market);

}  // namespace mc
