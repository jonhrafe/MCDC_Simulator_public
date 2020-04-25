

Clone or download the repository in your system
	``git clone https://github.com/jonhrafe/MCDC_Simulator_public.git``

## Compilation
### Linux and macOs systems (supported systems)
**Requirements:**  GCC compiler with g++ with c++11 support  (any version post 2011).

 - **Via the default command interpreter (sh)**.  In a terminal, inside the
   main directory, just type:

   `sh compile.sh`

 - **Via the Makefile**. Place the interpreter in the *src* folder and run the make file:
 
 ``cd src``
 
  ``make``

A self-contained binary file will be created in the root folder called:  
	*MC-CD_Simulator*

Voilá

### Windows

There's no official support for Windows apart from using the [Official Windows subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10). Once installed, the same installation process can be followed as above.

In addition, [MinGw](http://www.mingw.org/), and [Cygwin](https://www.cygwin.com/) can be installed to compile the source directly. In both cases, the full GCC compiler suite packages need to be installed and then use the application dependent gcc .exe file to run:

    <gcc.exe> -O3 -std=c++11 -lpthread -std=c++0x -pthread -I. main.cpp simulablesequence.cpp vertex.cpp obstacle.cpp collision.cpp scheme.cpp voxel.cpp cylinder.cpp walker.cpp mcsimulation.cpp
      parallelmcsimulation.cpp trajectory.cpp triangle.cpp parameters.cpp plyobstacle.cpp pgsesequence.cpp
       dynamicsSimulation.cpp simerrno.cpp collisionsphere.cpp  cylindergammadistribution.cpp sentinel.cpp subdivision.cpp gradientwaveform.cpp propagator.cpp -o MC-CD_Simulator

where <g+.exe> should be replaced by g++ in the Cygwin64 Terminal, or the g++ executable of mingw (*g++.exe*) .

# Getting started

To test your compiled application follow the [first simulation example in the Getting started page.](https://github.com/jonhrafe/MCDC_Simulator_public/blob/master/instructions/GettingStarted.md)

Good luck. 
