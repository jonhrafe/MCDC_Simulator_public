//!  RandomEngine Class =========================================================================/
/*!
*   \details   Single seeded RNG abstraction for the MC-DC simulator.
*
*   Every stochastic draw in the simulator (walker placement, step direction,
*   substrate generation and the membrane percolation draw) must go through one
*   of these so that a fixed `seed` yields bit-reproducible runs. Historically
*   those draws used a mix of fresh `std::random_device` engines and the C
*   `rand()` call, none of which honoured the user seed.
*
*   seedFrom() derives an independent, deterministic stream from a base seed plus
*   a few integer keys (worker index, walker id, purpose). The key-mixing keeps
*   this interface ready for a future counter-based RNG (Philox/Threefry) keyed
*   on (walker, step, purpose) without having to touch any call site.
*
*   \author    RNG unification
*   \version   0.1
*==============================================================================================*/

#ifndef RNG_H
#define RNG_H

#include <random>
#include <cstdint>

/*! \class RandomEngine
 *  \brief Deterministic, seedable random source wrapping std::mt19937_64.
 */
class RandomEngine
{
public:

    //! Purpose tags so distinct draws keyed on the same ids stay independent.
    enum Purpose : uint64_t {
        PLACEMENT = 1,
        CROSSING  = 2,
        SUBSTRATE = 3
    };

    RandomEngine() { seed(0); }
    explicit RandomEngine(uint64_t s) { seed(s); }

    //! Seed directly from a single 64-bit value (run through the mixer first).
    void seed(uint64_t s) { engine.seed(mix(s)); }

    //! Derive a deterministic, independent stream from a base seed and up to
    //! three integer keys (e.g. worker index, walker id, purpose).
    void seedFrom(uint64_t base, uint64_t k1, uint64_t k2 = 0, uint64_t k3 = 0)
    {
        uint64_t h = mix(base);
        h = mix(h ^ (k1 + GOLDEN));
        h = mix(h ^ (k2 + GOLDEN));
        h = mix(h ^ (k3 + GOLDEN));
        engine.seed(h);
    }

    //! Uniform double in [0,1).
    double uniform() { return dist(engine); }

    //! Uniform double in [a,b).
    double uniform(double a, double b) { return a + (b - a) * dist(engine); }

    //! Access to the underlying engine (e.g. for std::gamma_distribution).
    std::mt19937_64 &generator() { return engine; }

private:

    static const uint64_t GOLDEN = 0x9E3779B97F4A7C15ULL;

    //! splitmix64 finalizer: scrambles a 64-bit value into a well-spread seed.
    static uint64_t mix(uint64_t z)
    {
        z += GOLDEN;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    std::mt19937_64 engine;
    std::uniform_real_distribution<double> dist{0.0, 1.0};
};

#endif // RNG_H
