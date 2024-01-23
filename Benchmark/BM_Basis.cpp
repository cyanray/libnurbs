#include <random>
#include <benchmark/benchmark.h>
#include <libnurbs/Basis/BSplineBasis.hpp>
#include <libnurbs/Core/KnotVector.hpp>

using namespace libnurbs;

static void BM_Basis_Evaluate(benchmark::State& state)
{
    int degree = 3;
    KnotVector U{{0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0}};
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    for (auto _: state)
    {
        state.PauseTiming();
        double u = distr(generator);
        state.ResumeTiming();
        BSplineBasis::Evaluate(degree, U, u);
    }
}
BENCHMARK(BM_Basis_Evaluate);

static void BM_Basis_EvaluateDerivative(benchmark::State& state)
{
    int degree = 3;
    KnotVector U{{0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0}};
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    for (auto _: state)
    {
        state.PauseTiming();
        double u = distr(generator);
        state.ResumeTiming();
        BSplineBasis::EvaluateDerivative(degree, U, u, 2);
    }
}
BENCHMARK(BM_Basis_EvaluateDerivative);




// Run the benchmark
BENCHMARK_MAIN();
