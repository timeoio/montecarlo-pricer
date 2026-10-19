#include "mcpricer/black_scholes.hpp"

#include <cmath>
#include <numbers>

namespace mc {

namespace {

// d1 = (ln(S/K) + (r - q + sigma^2 / 2) T) / (sigma sqrt(T))
double d1(const Option& option, const Market& m) {
    return (std::log(m.spot / option.strike) + (m.rate - m.dividend + 0.5 * m.vol * m.vol) * option.maturity) /
           (m.vol * std::sqrt(option.maturity));
}

}  // namespace

double norm_pdf(double x) {
    return std::exp(-0.5 * x * x) * std::numbers::inv_sqrtpi / std::numbers::sqrt2;
}

double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::numbers::sqrt2);
}

double bs_price(const Option& option, const Market& m) {
    const double d_1 = d1(option, m);
    const double d_2 = d_1 - m.vol * std::sqrt(option.maturity);
    const double spot_pv = m.spot * std::exp(-m.dividend * option.maturity);    // S e^{-qT}
    const double strike_pv = option.strike * std::exp(-m.rate * option.maturity);  // K e^{-rT}

    if (option.type == OptionType::Call) return spot_pv * norm_cdf(d_1) - strike_pv * norm_cdf(d_2);
    return strike_pv * norm_cdf(-d_2) - spot_pv * norm_cdf(-d_1);
}

double bs_vega(const Option& option, const Market& m) {
    return m.spot * std::exp(-m.dividend * option.maturity) * norm_pdf(d1(option, m)) * std::sqrt(option.maturity);
}

}  // namespace mc
