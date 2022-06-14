//!  Benchmark class ============================================================================/
/*!
    \details    Implementation for benchmark and testing. The class implements performance tests
                and accuracy of the collision system.
    \author     Jonathan Rafael
    \date       November 2016
================================================================================================*/

#ifndef BENCHMARK_H
#define BENCHMARK_H
#include "simerrno.h"
#include <vector>
#include "parallelmcsimulation.h"

/*! \class Benchmark
 *  \brief Implementation for benchmark and testing. The class implements performance tests
 *         and accuracy of the collision system.
 *
 *  This class serves to compare optimization methods and accuracy checks in the dynamics.
 *
 *  All the methods assume a folder structure: test1, test2, ..., test6. Where the results are saved.
 */
class Benchmark
{
public:

    std::string output_dir;                 /*!< Output directory where the output folders are.         */

    std::vector<double> elapsed_times;      /*!< Elapsed time of all outputs.                           */

    /*! \fn  Unique and basic constructor
     *  \brief Sets the output directory
     */
    Benchmark(std::string output);

    /*! \fn  starts the Benchmarks
     */
    void start();

};

#endif // BENCHMARKCPP_H
