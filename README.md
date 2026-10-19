# Monte Carlo Option Pricer

A small, dependency-free C++20 project that:

1. prices **European options** by **Monte Carlo** simulation, with **antithetic variates**,
2. checks the results against the **Black-Scholes** closed-form formula,
3. recovers the **implied volatility** from an option price with **Newton-Raphson**.

## 1. Monte Carlo pricing

Under the risk-neutral measure the stock follows a geometric Brownian motion, which has an exact solution at maturity:

$$S_T = S_0 \exp\left(\left(r - q - \tfrac{1}{2}\sigma^2\right)T + \sigma\sqrt{T}\,Z\right), \quad Z \sim \mathcal{N}(0,1)$$

The option price is the discounted expected payoff, $V_0 = e^{-rT}\,\mathbb{E}[\text{payoff}(S_T)]$. Monte Carlo estimates this expectation with an average:

1. draw $N$ independent standard normals $Z_i$;
2. compute $S_T^{(i)}$ and the discounted payoffs $X_i = e^{-rT}\,\text{payoff}(S_T^{(i)})$;
3. take the price as $\hat V = \frac{1}{N}\sum X_i$, with standard error $\hat\sigma_X / \sqrt{N}$.

By the central limit theorem, $\hat V \pm 1.96\,\hat\sigma_X/\sqrt{N}$ is a 95% confidence interval. The error decreases like $1/\sqrt{N}$, so 100× more paths give one more correct digit.

### Antithetic variates

If $Z$ is a standard normal, so is $-Z$. Each draw is therefore used twice, and the two payoffs are averaged:

$$Y_i = \tfrac{1}{2}\left(X(Z_i) + X(-Z_i)\right)$$

The estimator stays unbiased. Because the payoff is monotonic in $Z$, the two halves are negatively correlated, so $\text{Var}(Y)$ is smaller than the variance from two independent paths. The result is a lower standard error for the same number of paths, at almost no extra cost.

## 2. Implied volatility by Newton-Raphson

The market quotes option *prices*. The implied volatility is the $\sigma$ that makes the Black-Scholes price match the quoted price:

$$f(\sigma) = BS(\sigma) - P_{\text{market}} = 0$$

There is no closed-form solution, so the equation is solved with **Newton-Raphson**. The method replaces $f$ by its tangent and jumps to the point where the tangent crosses zero:

$$\sigma_{n+1} = \sigma_n - \frac{f(\sigma_n)}{f'(\sigma_n)} = \sigma_n - \frac{BS(\sigma_n) - P_{\text{market}}}{\text{vega}(\sigma_n)}$$

Black-Scholes gives the derivative exactly: it is the vega, $S e^{-qT}\,\varphi(d_1)\sqrt{T}$.

- **Quadratic convergence:** close to the solution, the number of correct digits roughly doubles at each iteration, so 3 or 4 steps reach machine precision.
- **Choice of starting point:** for options far from the money, a naive guess such as $\sigma_0 = 20\%$ can land where vega ≈ 0. The tangent is then flat and Newton jumps away. The code starts instead at $\sigma_0 = \sqrt{2\,|\ln(F/K)|/T}$ (Manaster-Koehler), the point where vega is maximal. $BS(\sigma)$ is convex below that point and concave above it, which guarantees convergence.
- **No-arbitrage bounds:** a solution exists only if the price lies strictly between the intrinsic value and its upper bound ($S e^{-qT}$ for a call). Otherwise the function throws an exception.

## Build and run

Requirements: a C++20 compiler (Clang 17+, GCC 13+, MSVC 19.30+) and CMake 3.20 or later.

```sh
cmake -S . -B build
cmake --build build
./build/mcpricer_demo                         # demo
ctest --test-dir build --output-on-failure    # tests
```

## Example output

```
Market: S0=100 r=0.05 q=0 sigma=0.2 | K=100 T=1

=== 1. Monte Carlo vs Black-Scholes (1,000,000 paths) ===
      Method              BS        MC    Std err                95% CI    Time
Call  Plain          10.4506   10.4723    0.01475   [10.4434, 10.5013]    45ms
Call  Antithetic     10.4506   10.4499    0.01043   [10.4294, 10.4703]    26ms
Put   Plain           5.5735    5.5616    0.00865   [5.5447, 5.5786]    45ms
Put   Antithetic      5.5735    5.5704    0.00664   [5.5574, 5.5834]    26ms

=== 2. Convergence: the error shrinks like 1/sqrt(N) ===
       Paths    MC price   |MC - BS|     Std err
         100     13.3083      2.8577     1.31792
        1000     10.6350      0.1845     0.34053
       10000     10.5167      0.0661     0.10409
      100000     10.4520      0.0015     0.03281
     1000000     10.4499      0.0007     0.01043
    10000000     10.4542      0.0036     0.00329

=== 3. Implied volatility by Newton-Raphson ===
Which sigma makes the BS call price equal to 12?
  Iter         sigma       BS(sigma)
     0  0.3162277660   14.8470470727
     1  0.2412039231   12.0032902377
     2  0.2411168941   12.0000000169
     3  0.2411168937   12.0000000000
     4  0.2411168937   12.0000000000
Implied vol = 24.1117% after 4 iterations

=== 4. Round trip: MC price -> implied vol ===
True sigma = 20.00%, implied from MC price 10.4499 -> 19.9981%
```

## Project structure

```
include/mcpricer/
  option.hpp         market data, option definition, payoff
  black_scholes.hpp  closed-form price and vega
  monte_carlo.hpp    Monte Carlo pricer (with antithetic variates)
  implied_vol.hpp    Newton-Raphson implied volatility
src/
  black_scholes.cpp
  monte_carlo.cpp
  implied_vol.cpp
  main.cpp           demo program
tests/
  test_pricer.cpp    tests against closed forms
```
