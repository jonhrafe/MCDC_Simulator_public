## Monte Carlo Diffusion and Collision Simulator

# Static building from sources:

####Assuming a bin folder:

Cretes a folder
`mkdir bin`

Compiles all the magic (statically, slowly, but without messy .o)

`g++ -O3 -std=c++11 -lpthread -std=c++0x -pthread -I. main.cpp simulablesequence.cpp vertex.cpp obstacle.cpp collision.cpp scheme.cpp voxel.cpp cylinder.cpp walker.cpp mcsimulation.cpp parallelmcsimulation.cpp trajectory.cpp triangle.cpp parameters.cpp plyobstacle.cpp pgsesequence.cpp dynamicsSimulation.cpp simerrno.cpp collisionsphere.cpp  cylindergammadistribution.cpp sentinel.cpp subdivision.cpp gradientwaveform.cpp propagator.cpp -o bin/MC-DC_Simulator`
