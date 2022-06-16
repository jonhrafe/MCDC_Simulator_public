#ifndef BENCHMARK_MCSIMULATION_H
#define BENCHMARK_MCSIMULATION_H

#include "mcsimulation.h"
#include <parameters.h>
#include <nvtx3/nvToolsExt.h> 

using namespace std;


class Benchmark
{
public:

    unsigned int benchmark_id;                  /*!< Id of the benchmark test to run */ 

    Parameters benchmark_params;                /*!< Parameters instance \see :Parameters:                              */

    std::vector <PLYObstacle> plyObstacles_list; /*!< vector with all the instances of PLYObstacles                      */
    std::vector <Sphere> spheres_list;           /*!< vector with all the instances of spheres                         */
    std::vector <Cylinder> cylinders_list;           /*!< vector with all the instances of cylinders                         */


    MCSimulation* simulation;                   /*!< Pointer to MCSimulation instance                       */

    /*! \fn  Benchmark.
     *  \brief  Default constructor.
     */
    Benchmark(); 

    /*! \fn  Benchmark.
     *  \brief Constructor with benchmark test id.
     */
    Benchmark(unsigned int id_); 


    /*! \fn  ~Benchmark.
     *  \brief Main destructor. Frees dynamicly allocated memory instances.
     */
    ~Benchmark();


    /*! \fn  startBenchmark.
     *  \brief Launch benchmark test.
     */
    void startBenchmark();


private:

        
    /*! \fn  selectBenchmark.
     *  \brief Load benchmark parameters.
    */
    void selectBenchmark();

    /*! \fn  writeBenchmark.
     *  \brief Save benchmark test results.
    */
    void writeBenchmark();


    /*! \fn  loadSpheres.
     *  \brief Load sphere list
    */
    void loadSpheres(string path_sphere_list, Parameters &params);

    /*! \fn  loadVoxel.
     *  \brief Load voxel
    */
    void loadVoxel(string path_voxel, Parameters &params);

};

#endif BENCHMARK_H