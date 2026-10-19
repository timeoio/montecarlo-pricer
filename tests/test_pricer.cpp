// Minimal self-contained test runner (no external dependency).
// Monte Carlo checks use statistical tolerances: the MC price must lie within
// a few standard errors of the closed-form price.

#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "mcpricer/black_scholes.hpp"
#include "mcpricer/implied_vol.hpp"
#include "mcpricer/monte_carlo.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "  FAILED: " << message << '\n';
    }
}

void check_close(double actual, double expected, double tolerance, const std::string& what) {
    check(std::abs(actual - expected) <= tolerance,
          std::format("{}: got {:.8f}, expected {:.8f} (tol {:.8f})", what, actual, expected, tolerance));
}

template <typename Exception, typename F>
bool throws(F&& f) {
    try {
        f();
    } catch (const Exception&) {
        return true;
    }
    return false;
}

const mc::Market market{.spot = 100.0, .rate = 0.05, .dividend = 0.02, .vol = 0.25};
const mc::Option call{.type = mc::OptionType::Call, .strike = 105.0, .maturity = 1.5};
const mc::Option put{.type = mc::OptionType::Put, .strike = 105.0, .maturity = 1.5};
constexpr double num_std_errors = 4.0;  // false-failure probability ~ 6e-5 per check

void test_bs_reference_value() {
    // Textbook value (Hull): S=100, K=100, r=5%, sigma=20%, T=1 -> 10.4506
    const mc::Market m{.spot = 100.0, .rate = 0.05, .dividend = 0.0, .vol = 0.20};
    const mc::Option atm{.type = mc::OptionType::Call, .strike = 100.0, .maturity = 1.0};
    check_close(mc::bs_price(atm, m), 10.4506, 1e-4, "BS call");
}

void test_bs_put_call_parity() {
    // C - P = S e^{-qT} - K e^{-rT}
    const double T = call.maturity;
    const double forward = market.spot * std::exp(-market.dividend * T) - call.strike * std::exp(-market.rate * T);
    check_close(mc::bs_price(call, market) - mc::bs_price(put, market), forward, 1e-10, "put-call parity");
}

void test_mc_matches_bs() {
    for (const auto& option : {call, put}) {
        for (bool antithetic : {false, true}) {
            const auto r = mc::mc_price(option, market, {.num_paths = 400'000, .antithetic = antithetic});
            check_close(r.price, mc::bs_price(option, market), num_std_errors * r.std_error,
                        std::format("MC {} (antithetic={})", option.type == mc::OptionType::Call ? "call" : "put",
                                    antithetic));
        }
    }
}

void test_antithetic_lowers_std_error() {
    const double plain = mc::mc_price(call, market, {.num_paths = 400'000}).std_error;
    const double anti = mc::mc_price(call, market, {.num_paths = 400'000, .antithetic = true}).std_error;
    check(anti < plain, "antithetic reduces std error at equal number of paths");
}

void test_reproducible() {
    check(mc::mc_price(call, market).price == mc::mc_price(call, market).price, "same seed gives the same price");
}

void test_implied_vol_round_trip() {
    // Price with a known sigma, then check Newton-Raphson recovers it.
    for (double strike : {60.0, 90.0, 100.0, 120.0, 180.0}) {
        for (double sigma : {0.05, 0.2, 0.6, 1.2}) {
            for (auto type : {mc::OptionType::Call, mc::OptionType::Put}) {
                const mc::Option option{.type = type, .strike = strike, .maturity = 1.0};
                mc::Market m = market;
                m.vol = sigma;
                const double price = mc::bs_price(option, m);
                // Skip options whose time value is ~0 (deep in or out of the money):
                // their price carries no information about sigma.
                const double intrinsic = std::max(type == mc::OptionType::Call
                                                      ? m.spot * std::exp(-m.dividend) - strike * std::exp(-m.rate)
                                                      : strike * std::exp(-m.rate) - m.spot * std::exp(-m.dividend),
                                                  0.0);
                if (price - intrinsic < 1e-6) continue;
                const auto iv = mc::implied_vol(option, market, price);
                check_close(iv.vol, sigma, 1e-6, std::format("implied vol K={} sigma={}", strike, sigma));
            }
        }
    }
}

void test_implied_vol_rejects_arbitrage() {
    check(throws<std::invalid_argument>([] { mc::implied_vol(call, market, 1000.0); }), "price above spot throws");
    check(throws<std::invalid_argument>([] { mc::implied_vol(call, market, 0.0); }), "zero price throws");
}

void test_invalid_inputs() {
    const mc::Option bad{.type = mc::OptionType::Call, .strike = 100.0, .maturity = -1.0};
    check(throws<std::invalid_argument>([&] { mc::mc_price(bad, market); }), "negative maturity throws");
}

}  // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"bs_reference_value", test_bs_reference_value},
        {"bs_put_call_parity", test_bs_put_call_parity},
        {"mc_matches_bs", test_mc_matches_bs},
        {"antithetic_lowers_std_error", test_antithetic_lowers_std_error},
        {"reproducible", test_reproducible},
        {"implied_vol_round_trip", test_implied_vol_round_trip},
        {"implied_vol_rejects_arbitrage", test_implied_vol_rejects_arbitrage},
        {"invalid_inputs", test_invalid_inputs},
    };
    for (const auto& [name, test] : tests) {
        const int before = failures;
        test();
        std::cout << (failures == before ? "[ OK ] " : "[FAIL] ") << name << '\n';
    }
    std::cout << (failures == 0 ? "All tests passed\n" : std::format("{} check(s) failed\n", failures));
    return failures == 0 ? 0 : 1;
}
