# Automatic initialization of spins in the internal or external compartment

In this example we will show several options for the automatic initialization of the spins in two geometries (obstacles):
- [Gamma distributed cylinders](#Gamma-distributed-cylinders)
- [PLY models](#PLY-Meshes)


Before starting be sure to download and compile the latest version: [Installation instructions](https://github.com/jonhrafe/MCDC_Simulator_public/blob/master/instructions/compilation.md) 

# Gamma distributed cylinders

__The Setup:__ Below we define a toy example of a substrate containing cylinders with radii samples from a Gamma distribution. The .conf file [can be found in here](https://github.com/jonhrafe/MCDC_Simulator_public/tree/master/docs/conf_file_examples)

```
N 1000
T 1000
duration 0.0359000000
diffusivity 0.6e-9

exp_prefix instructions/demos/output/cylinder_gamma_packing_test

scheme_file docs/scheme_files/PGSE_sample_scheme.scheme
scale_from_stu 1

write_txt 1
write_bin 0
write_traj_file 1

<obstacle>
<cylinder_gamma_packing>
alpha 1.5
beta 0.5
icvf 0.70
num_cylinders 1000
</cylinder_gamma_packing>
</obstacle>

ini_walkers_pos extra
num_process 1
<END>
```

###  Basic parameters description:
`` N 1000 `` Defines the number of spins to be placed to 1000 spins.

``` T 1000 ``` Defines the total number of steps placed to 1000 steps.

``` duration 0.036 ``` The total duration is set to 0.036 seconds.

``` diffusivity  0.6e-9 ``` sets the diffusion coefficient of the medium to 0.6 m^2/s (near ex-vivo diffusion)

``` exp_prefix [string] ``` output path and prefix for the experiment

``` scheme_file [string] ``` (dummy) PGSE scheme file provided in the repository

``` scale_from_stu 1 ``` 1 by the default. 1 if the PGSE is in __SI base unit__: meters, seconds. 

``` write_txt [0,1] ``` flag to indicate if the output files should be written in txt format with reduced numerical precision. __WARNING__ this will produce big files that scale quadratically with the N and T parameters.

``` write_bin [0,1] ``` flag to indicate if the output files should be written in binary format with full floating (float32) precision (recommended).

``` write_traj_file [0,1] ``` flag to indicate if the particle’s trajectories should be written to disk (warning, very big files). The trajectory file will be written in txt, or binary format depending on the write_txt and write_bin flags.

``` num_process [unsigned int] ``` short for number of processors, defines the number of processors/cores to use for the simulation. WARNING: increasing the number of processors while writing a trajectory file will split the .traj file according to the number of processors.

``` <END> ```The .conf must end with this tag. All other parameters or text after this tag are ignored.


### Defining obstacles:

In order to define a new obstacle you need to open and close the obstacles tags: ``` <obstacle> </obstacle> ```

Inside we define a new obstacle formed of gamma distributed cylinders using the tags: ``` <cylinder_gamma_packing>  <cylinder_gamma_packing> ``` which requires of the 4 parameters listed below:

``` shape 1.5 ```  Shape parameter of the gamma distribution in __micrometers__.
``` scale 0.5 ``` Scale parameter of the gamma distribution in __micrometers__.
``` icvf 0.70 ``` Desired (target) ICVF (ratio of the space filled by cylinders in the voxel).
``` num_cylinders ``` 1000 Total number of cylinders to pack inside the voxel. The voxel size will be automatically adjusted.

The parameters above will then define a gama distribution with a mean diameter (shape/scale) of 1 micrometer. The target ICVF will be 70% of the total space. Note that the packing algorithm will try to fit as close as possible to the target parameter, however, ICVF greater than 0.8 are nearly unfeasible.

### Spin initialization:

In order to sample the spins not uniformly, but only in one compartment you can use the parameter ini_walkers_pos
``` ini_walkers_pos extra ``` will then only initialize spins outside the defined cylinders (extra-axona space).

## Simulation outputs:

The .conf file is included in the [doc folder](https://github.com/jonhrafe/MCDC_Simulator_public/tree/master/docs/conf_file_examples): 

Tu run the simulation simply call the simulator executable file in the project’s root directory as:

``` 
./MC-DC_Simulator docs/conf_file_examples/gammaDistributedCylinders.conf
```
The outputs will be stored, as defined, in ``` instructions/demos/output/ ```

_Output files:_

``` cylinder_gamma_packing_test_DWI.txt ```  Real part of the un-normalized diffusion signal (in this case, issued only by the extra-axonal space).

``` cylinder_gamma_packing_test_DWI_img.txt ```  Imaginary part part of the un-normalized diffusion signal (likely 0's).

``` cylinder_gamma_packing_test_gamma_distributed_cylinder_list.txt ``` List with the position and radius of all the packed cylinders

``` cylinder_gamma_packing_test_simulation_info.txt ``` Simulation info file (important) contains all the simulation info and possible warnings.

``` cylinder_gamma_packing_test_0.hdr.txt ``` Trajectory file header.

``` cylinder_gamma_packing_test_0.traj.txt ``` All the particle's trajectories (for visualization purposes mostly).


### Visualizing the outputs

The trajfile (binary or in text) contains all the particle's positions over time. The total number of positions in the file is then ``N * (T+1)* 3``. This is: 3 positions (x,y,z), for each step from 0 to T, for all N particles. 

The format is then a list from the first position of the first particle, to the last position to the last particle.

x<sub>{1,0}</sub>

x<sub>{1,0}</sub>

y<sub>{1,0}</sub>

z<sub>{1,0}</sub>

...

x<sub>{1,T}</sub>

y<sub>{1,T}</sub>

z<sub>{1,T}</sub>

...

...

x<sub>{N,T}</sub>

x<sub>{N,T}</sub>

x<sub>{N,T}</sub>

Then, you can visualize each particle's trajectory and get something like this:
![simulation](https://user-images.githubusercontent.com/4105920/88835884-31fb9800-d1d6-11ea-8dcb-5210ae50793e.gif)

# PLY Meshes

Similarly to the previous example, we can use the automatic initialization to simulate only in the internal compartment (representing the intra-axonal space) on closed meshes of arbitrary shape.

Bellow we show a visualization of a convoluted mesh that we will use in this example. The raw mesh can be found inside the repo's examples folder [in here](https://github.com/jonhrafe/MCDC_Simulator_public/tree/master/instructions/meshes).

![simulation](https://user-images.githubusercontent.com/4105920/88836514-18a71b80-d1d7-11ea-9e2a-a43df479a889.gif)

__The Setup:__ Below we list the .conf file parameters. In this case we will only set 100 steps over 0.02 seconds with a very low diffusion coefficient of 0.6e-10 m^2/s. The .conf file [can be found in here](https://github.com/jonhrafe/MCDC_Simulator_public/tree/master/docs/conf_file_examples)

```
N 1000
T 100
duration 0.020
diffusivity 0.6e-10
exp_prefix instructions/demos/output/mesh_initialization_test

scheme_file docs/scheme_files/PGSE_sample_scheme.scheme
scale_from_stu 1

write_txt 1
write_bin 0
write_traj_file 1

<obstacle>
ply instructions/meshes/PorusMedia.ply
ply_scale 0.001
</obstacle>

<voxels>
-.010 -.010 -.050
.010 .010 .050
</voxels>

<spawning_area>
-.005 -.005 -0.005
.005 .005 .005
</spawning_area>

ini_walkers_pos intra

num_process 1

<END>
```
To run the simulation type:

``` 
./MC-DC_Simulator docs/conf_file_examples/meshIntraInitialization.conf
```
###  New parameters description:
Above we have set a new type of obstacle defined inside the tags ` <obstacle>` `</obstacle>` named `ply` obstacle.

`ply [string]` The ply obstacle descriptor requires a PLY mesh at the __milimeter's (mm) scale__. This means that the vertex units will be __assumed to be in milimeters__. The mesh model must be __completely triangulated__  and must not contain any other information but the edges and vertices of the model.

`ply_scale [float]` Optionally, you can set a scale factor to be applied to the mesh vertices. In our example, since our model scale is originally saved in micrometers for visualization purposes, we have defined a custom scale of 0.001 to scale the model back to millimeters.


`<voxels> [float] [float] [float] [float] [float] [float]</voxels>` The  voxel tags define the voxel limits (volume where the diffusion signal will be synthesized). The contained 6 numbers must define the lower and upper limits of the __voxel in milimeters__.

`<spawning_area> [float] [float] [float] [float] [float] [float] <spawning_area>` Optionally you can define a custom spawning area where the particles will be initialized uniformly. The followed 6 numbers must define the lower and upper limits of the sampling area __in milimeters__. This option is especially useful when simulating specific compartments. If this command is omitted the particles will be uniformly sampled inside the defined voxel.

`ini_walkers_pos intra` This option will forze the initialization of particles exclusively inside the provided mesh model. The algorithm first initializes __uniformely in the selected sampling area__ and then discards any particle that landed outside the mesh. 

### Visualizing the outputs

As before, you can visualize the trajectory file to check the custom initialization. 
![mesh](https://user-images.githubusercontent.com/4105920/88842439-bef71f00-d1df-11ea-9616-ff607cc2b6fe.gif)





