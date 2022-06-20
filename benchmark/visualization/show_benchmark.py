import os
import numpy as np 

import matplotlib.pyplot as plt 

def give_b(exp_proto):

    gyro    = GYRO/1e3
    b       = np.around((exp_proto[:, 3]*exp_proto[:, 5]*gyro)**2 * (exp_proto[:, 4]-exp_proto[:, 5]/3), 6) *1e12 #s/mm2
    return b


def read_scheme(path_scheme):

    with open(path_scheme) as f: lines = f.readlines()

    exp_type = lines[0].split('VERSION: ')[-1].split('\n')[0]

    if (exp_type == 'STEJSKALTANNER') |  (exp_type == '1'):
        exp_proto = np.zeros((len(lines[1:]), 7))

        for i,line in enumerate(lines[1:]):            
            exp_proto[i] = [float(a) for a in line.split()]

    return exp_proto

def load_S(path_exp, bin=True):

    if bin:
        with open(path_exp+'_DWI.bfloat', 'rb') as f:
            S_real = np.fromfile(f, '<f4')
        
        try:
            with open(path_exp+'_DWI_img.bfloat', 'rb') as f:
                S_img = np.fromfile(f, '<f4')
        except FileNotFoundError:
            S_img = None
    else: 
        with open(path_exp+'_DWI.txt', 'r') as f:
            S_real = np.loadtxt(f) 

        try:
            with open(path_exp+'_DWI_img.txt', 'r') as f:
                S_img = np.loadtxt(f) 
        except FileNotFoundError:
            S_img = None

        try:
            with open(path_exp+'_DWI_intra.txt', 'r') as f:
                S_intra = np.loadtxt(f) 
            with open(path_exp+'_DWI_extra.txt', 'r') as f:
                S_extra = np.loadtxt(f) 
            
        except FileNotFoundError:
            S_intra, S_extra = None, None
        
    return S_real, S_img, S_intra, S_extra



def show_benchmark():

    # 1. Load scheme
    exp_protocol = read_scheme(PATH_SCHEME)

    # 2. Load signals 

    # 2.1 Spheres 
    exp_sphere_list     = "R_2_R_4_v_50_ICVF_0.57_gaussian_sphere_packing"
    path_sphere_list    = os.path.join(PATH_OUTPUT_FOLDER, "sphere", exp_sphere_list)

    S_real_sphere, S_img_sphere, S_intra_sphere, S_extra_sphere = load_S(path_sphere_list, bin=False)


    # 3. Visualization
    # 3.1 Signals
    
    f_sig, ax_sig = plt.subplots(1, 1, figsize=FIG_S)
    ax_sig.plot(S_real_sphere, label="Full")
    ax_sig.plot(S_intra_sphere, label="Intra")
    ax_sig.plot(S_extra_sphere, label="Extra")
    
    plt.legend()
    f_sig.savefig(os.path.join(PATH_SAVE, "signals.png"))
    

    return 

def main():
    show_benchmark()

    return 

if __name__=="__main__":

    # CONSTANT
    GYRO        = 267.51525e3
    FIG_H           = 9
    FIG_W           = 9
    FIG_S   =(FIG_H, FIG_W)


    #PATH
    PATH_SAVE           = "./benchmark/visualization"
    PATH_OUTPUT_FOLDER  = "./benchmark/output"
    PATH_SCHEME         = "./benchmark/scheme_files/10shell_fixed_DdTe.scheme"

    main()
