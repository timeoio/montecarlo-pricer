#include "mcpricer/monte_carlo.hpp"

#include <cmath>
#include <random>
#include <stdexcept>

namespace mc {

McResult mc_price(const Option& option, const Market& m, const McConfig& config) {
    if (option.maturity <= 0.0 || m.spot <= 0.0 || m.vol <= 0.0)
        throw std::invalid_argument("maturity, spot and vol must be positive");

    // With antithetic variates one sample = one pair of paths (Z, -Z).
    const std::size_t num_samples = config.antithetic ? config.num_paths / 2 : config.num_paths;
    if (num_samples < 2) throw std::invalid_argument("not enough paths");

    std::mt19937_64 rng(config.seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    // A European payoff only depends on S(T), and GBM gives it exactly:
    //   S(T) = S(0) * exp((r - q - sigma^2 / 2) T + sigma sqrt(T) Z),  Z ~ N(0, 1)
    const double T = option.maturity;
    const double drift = (m.rate - m.dividend - 0.5 * m.vol * m.vol) * T;
    const double diffusion = m.vol * std::sqrt(T);
    const double discount = std::exp(-m.rate * T);

    // Welford's online algorithm: updates the mean and the sum of squared
    // deviations m2 = sum (x_i - mean)^2 at each sample, in O(1) memory and a
    // single pass. Unlike sum(x^2) - n * mean^2, it never subtracts two large,
    // nearly equal numbers, so it does not lose precision when the mean is
    // large compared with the standard deviation.
    double mean = 0.0;
    double m2 = 0.0;
    for (std::size_t i = 0; i < num_samples; ++i) {
        const double z = normal(rng);
        double x = discount * payoff(option, m.spot * std::exp(drift + diffusion * z));
        if (config.antithetic) {
            // The mirror path -Z is as likely as Z. Their payoffs are negatively
            // correlated, so averaging them has a lower variance than two
            // independent paths.
            const double x_anti = discount * payoff(option, m.spot * std::exp(drift - diffusion * z));
            x = 0.5 * (x + x_anti);
        }
        const double delta = x - mean;              // deviation from the old mean
        mean += delta / static_cast<double>(i + 1);
        m2 += delta * (x - mean);                   // old deviation * new deviation
    }

    // Price = sample mean; std error = sample std deviation / sqrt(N)
    const double n = static_cast<double>(num_samples);
    const double variance = m2 / (n - 1.0);
    return {mean, std::sqrt(variance / n)};
}

}  // namespace mc
