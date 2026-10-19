#include "mcpricer/implied_vol.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "mcpricer/black_scholes.hpp"

namespace mc {

double initial_guess(const Option& option, const Market& market) {
    // Manaster-Koehler starting point: the sigma where vega is maximal, i.e. the
    // inflection point of BS(sigma). The price is convex in sigma below it and
    // concave above, which guarantees Newton converges from there. A naive
    // guess (e.g. 20% for a far out-of-the-money option) can land where vega ~ 0:
    // the tangent is flat and Newton jumps away.
    const double forward = market.spot * std::exp((market.rate - market.dividend) * option.maturity);
    const double guess = std::sqrt(2.0 * std::abs(std::log(forward / option.strike)) / option.maturity);
    return guess > 1e-3 ? guess : 0.2;  // at the money any reasonable guess works
}

ImpliedVolResult implied_vol(const Option& option, const Market& market, double target_price, double tolerance,
                             int max_iterations) {
    // The BS price increases with sigma, from its intrinsic value (sigma -> 0)
    // to an upper bound (sigma -> infinity). Outside this range no sigma exists.
    const double spot_pv = market.spot * std::exp(-market.dividend * option.maturity);
    const double strike_pv = option.strike * std::exp(-market.rate * option.maturity);
    const double lower = option.type == OptionType::Call ? std::max(spot_pv - strike_pv, 0.0)
                                                         : std::max(strike_pv - spot_pv, 0.0);
    const double upper = option.type == OptionType::Call ? spot_pv : strike_pv;
    if (target_price <= lower || target_price >= upper)
        throw std::invalid_argument("price outside no-arbitrage bounds: no implied volatility");

    Market m = market;
    m.vol = initial_guess(option, market);
    ImpliedVolResult result{m.vol, 0, {}};

    for (int i = 1; i <= max_iterations; ++i) {
        // f(sigma) = BS(sigma) - target is the function we want to set to zero.
        // Newton replaces f by its tangent and jumps to where the tangent hits zero.
        const double f = bs_price(option, m) - target_price;
        const double vega = bs_vega(option, m);  // f'(sigma)
        if (vega < 1e-12) throw std::runtime_error("vega too small: Newton-Raphson cannot proceed");

        double next = m.vol - f / vega;
        if (next <= 0.0) next = m.vol / 2.0;  // safeguard: sigma must stay positive

        result.history.push_back(next);
        result.iterations = i;
        const bool converged = std::abs(next - m.vol) < tolerance;
        m.vol = next;
        if (converged) {
            result.vol = m.vol;
            return result;
        }
    }
    throw std::runtime_error("Newton-Raphson did not converge");
}

}  // namespace mc
