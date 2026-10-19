#pragma once

#include <cstddef>
#include <cstdint>

#include "mcpricer/option.hpp"

namespace mc {

struct McConfig {
    std::size_t num_paths = 1'000'000;
    std::uint64_t seed = 42;   // same seed -> same result
    bool antithetic = false;   // pair each draw Z with -Z
};

struct McResult {
    double price;
    double std_error;  // standard deviation of the estimator: shrinks like 1 / sqrt(N)

    // 95% confidence interval
    double ci_low() const { return price - 1.96 * std_error; }
    double ci_high() const { return price + 1.96 * std_error; }
};

McResult mc_price(const Option& option, const Market& market, const McConfig& config = {});

}  // namespace mc
