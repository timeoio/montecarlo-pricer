#pragma once

#include <algorithm>

namespace mc {

// Black-Scholes market: under the risk-neutral measure the spot follows
//   dS = (r - q) S dt + sigma S dW
struct Market {
    double spot;      // S(0)
    double rate;      // continuously compounded risk-free rate r
    double dividend;  // continuous dividend yield q
    double vol;       // volatility sigma
};

enum class OptionType { Call, Put };

// European option: can only be exercised at maturity.
struct Option {
    OptionType type;
    double strike;    // K
    double maturity;  // T, in years
};

// Cash flow received at maturity when the spot is spot_at_maturity.
inline double payoff(const Option& option, double spot_at_maturity) {
    if (option.type == OptionType::Call) return std::max(spot_at_maturity - option.strike, 0.0);
    return std::max(option.strike - spot_at_maturity, 0.0);
}

}  // namespace mc
