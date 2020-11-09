# Getting Started

## Basic usage 

    MC-DC_Simulator configuration_file.conf

 - *MCDC_Simulator* :  Application name
 - *configuration_file.conf:* file with ALL the simulation parameters listed following the simulator syntaxis. 

## Test you installation: 

Inside the folder *docs\conf_file_examples* there's several examples discussed in the tutorials.

To test your installation, run the following line placed in the simulator root folder: 

  ./MC-DC_Simulator docs/conf_file_examples/freeDiffusion.conf 

If everything goes well, simulation output will be saved inside: *MCDC_Simulator_public\instructions\demos*.

![enter image description here](https://user-images.githubusercontent.com/4105920/80267704-b2c06a80-86a2-11ea-9706-daf9ab7eae6e.gif)

Below we review the basic format and parameters included inside *docs/conf_file_examples/freeDiffusion.conf* 

## Free diffusion simulation:  

**Input:** 
 - *freeDiffusion.conf*
 
**Outputs:** 
All the outputs are stored ins separated files with the indicated "exp_prefix" name and the output appended name. 

**Output files:** 
 - real part of the DW-MRI signal (*"_DWI"*).
 - imaginary part of the DW_MRI signal (*"_DWI_img"*)
 - info file ("*_info*")
 
Below, the content of the configuration file, freeDiffusion.conf, is shown. 

    N 10000
    T 1000
    duration 0.100
    diffusivity 0.6e-9
    scale_from_stu 1

    exp_prefix instructions/demos/output/free_diffusion_test
    
    scheme_file docs/scheme_files/PGSE_sample_scheme.scheme
    
    write_txt 1
    write_bin 0
    write_traj_file 0
    
    <voxel>
    0.0 0.0 0.0
    0.5 0.5 0.5
    </voxel>
    
    num_process 5
    
    <END>

The configuration files will **read all the listed parameters until finding  the 
< END> tag**. all the parameters should be listed in one line with the name of the parameters followed by a space and the indicated value. Any path should be either a **full system path**, or **relative to the consoles execution path.**

Below the description of the parameters in the confile.

 - **`N`** [int]: total number of spin particles to diffuse.
 - **`T`** [int]: total number of steps to perform during the duration of the experiment. 
 - **`duration`** [float]: total diffusion duration in seconds, this duration should be at least as long as the longest echo time (TE) in the acquisition protocol  acquisition.
 - **`diffusivity`** [float] diffusion coefficient of the medium in m^2/s for standard units, mm^2/ms otherwise.
 - **`scale_from_stu`** [0,1] flag to indicate if the protocol's scheme_file, and the experiment's duration and diffusivity, are in standard units (SI units; m,s,etc.)   
 - **`exp_prefix`** [string] output path and prefix for the experiment 
 - scheme_file docs/scheme_files/PGSE_sample_scheme.scheme
 - **`write_txt`** [0,1]  flag to indicate if the output files should be written in txt format with reduced numerical precision.
 - **`write_bin`** [0,1] flag to indicate if the output files should be written in binary format with full floating (float32)  precision **(recommended)**.
 - **`write_traj_file`** [0,1] flag to indicate if the particles trajectories should be written to disk (**warning, very big files**). The trajectory file will be written in txt, or binary format depending on the **write_txt** and **write_bin** flags.
 - **`<voxel> <\voxel>`** Voxel tags to define the simulation voxel in MILLIMETRES. The first 3 numbers defines the **minimum voxel limit (x_min,y_min,z_min)** followed by 3 numbers defining the **maximum voxel limit (x_max, y_max, z_min)**.
 - **`num_process`** [unsigned int] short for number of processors, defines the number of processors/cores to use for  the simulation. 

**Experiment description:** 
In our free diffusion example, we have defined a simulations with 10,000 spins that will diffuse during 100 milliseconds, doing 1,000 steps (one step every 0.1 milliseconds), the output files will be stored in the folder "instructions/demos/output/free_diffusion_test" with the prefix "free_diffusion_test". The diffusion signal will be computed using the PGSE scheme file "PGSE_sample_scheme.scheme", which contains 270 acquisition with a maximum TE of 100 milliseconds.  Te results will be written in txt format, and no trajfile will be saved. We defined an arbitrary voxel from      0.0mm 0.0mm 0.0mm to 0.5mm, 0.5mm, 0.5mm. Finally we have assigned 5 processors/cores to the simulation. 

