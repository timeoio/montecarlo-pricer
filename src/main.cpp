// Demo: prices a European call by Monte Carlo, compares it with Black-Scholes,
// then recovers the implied volatility of a price with Newton-Raphson.

#include <chrono>
#include <cmath>
#include <format>
#include <iostream>
#include <string_view>

#include "mcpricer/black_scholes.hpp"
#include "mcpricer/implied_vol.hpp"
#include "mcpricer/monte_carlo.hpp"

namespace {

void section(std::string_view title) {
    std::cout << std::format("\n=== {} ===\n", title);
}

}  // namespace

int main() {
    const mc::Market market{.spot = 100.0, .rate = 0.05, .dividend = 0.0, .vol = 0.20};
    const mc::Option call{.type = mc::OptionType::Call, .strike = 100.0, .maturity = 1.0};
    const mc::Option put{.type = mc::OptionType::Put, .strike = 100.0, .maturity = 1.0};

    std::cout << std::format("Market: S0={} r={} q={} sigma={} | K={} T={}\n", market.spot, market.rate,
                             market.dividend, market.vol, call.strike, call.maturity);

    section("1. Monte Carlo vs Black-Scholes (1,000,000 paths)");
    std::cout << std::format("{:<6}{:<12}{:>10}{:>10}{:>11}{:>22}{:>8}\n", "", "Method", "BS", "MC", "Std err",
                             "95% CI", "Time");
    for (const auto& [name, option] : {std::pair{"Call", call}, std::pair{"Put", put}}) {
        for (bool antithetic : {false, true}) {
            const auto start = std::chrono::steady_clock::now();
            const auto r = mc::mc_price(option, market, {.antithetic = antithetic});
            const std::chrono::duration<double, std::milli> ms = std::chrono::steady_clock::now() - start;
            std::cout << std::format("{:<6}{:<12}{:>10.4f}{:>10.4f}{:>11.5f}   [{:.4f}, {:.4f}]{:>6.0f}ms\n", name,
                                     antithetic ? "Antithetic" : "Plain", mc::bs_price(option, market), r.price,
                                     r.std_error, r.ci_low(), r.ci_high(), ms.count());
        }
    }

    section("2. Convergence: the error shrinks like 1/sqrt(N)");
    const double bs = mc::bs_price(call, market);
    std::cout << std::format("{:>12}{:>12}{:>12}{:>12}\n", "Paths", "MC price", "|MC - BS|", "Std err");
    for (std::size_t n = 100; n <= 10'000'000; n *= 10) {
        const auto r = mc::mc_price(call, market, {.num_paths = n, .antithetic = true});
        std::cout << std::format("{:>12}{:>12.4f}{:>12.4f}{:>12.5f}\n", n, r.price, std::abs(r.price - bs),
                                 r.std_error);
    }
    std::cout << "(x100 paths -> error / 10)\n";

    section("3. Implied volatility by Newton-Raphson");
    const double market_price = 12.0;
    std::cout << std::format("Which sigma makes the BS call price equal to {}?\n", market_price);
    const auto iv = mc::implied_vol(call, market, market_price);
    mc::Market start = market;
    start.vol = mc::initial_guess(call, market);
    std::cout << std::format("{:>6}{:>14}{:>16}\n", "Iter", "sigma", "BS(sigma)");
    std::cout << std::format("{:>6}{:>14.10f}{:>16.10f}\n", 0, start.vol, mc::bs_price(call, start));
    for (std::size_t i = 0; i < iv.history.size(); ++i) {
        mc::Market m = market;
        m.vol = iv.history[i];
        std::cout << std::format("{:>6}{:>14.10f}{:>16.10f}\n", i + 1, iv.history[i], mc::bs_price(call, m));
    }
    std::cout << std::format("Implied vol = {:.4f}% after {} iterations (quadratic convergence:\n"
                             "the number of correct digits roughly doubles at each step)\n",
                             100.0 * iv.vol, iv.iterations);

    section("4. Round trip: MC price -> implied vol");
    const auto r = mc::mc_price(call, market, {.antithetic = true});
    const double vol_from_mc = mc::implied_vol(call, market, r.price).vol;
    std::cout << std::format("True sigma = {:.2f}%, implied from MC price {:.4f} -> {:.4f}%\n", 100.0 * market.vol,
                             r.price, 100.0 * vol_from_mc);
    std::cout << std::format("MC std error {:.4f} / vega {:.2f} = vol error of ~{:.3f}%\n", r.std_error,
                             mc::bs_vega(call, market), 100.0 * r.std_error / mc::bs_vega(call, market));
}
