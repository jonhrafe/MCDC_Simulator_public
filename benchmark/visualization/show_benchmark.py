import os
import numpy as np 

import matplotlib.pyplot as plt 

def add_plot_info(ax, title='', legend='', xlabel='', ylabel='', fs=15, xlims=None, ylims=None, xscale=None, yscale=None, xticks=None, yticks=None):
    ax.set_title(title, fontsize=fs)
    ax.set_xlabel(xlabel, fontsize=int(fs*0.9))
    ax.set_ylabel(ylabel, fontsize=int(fs*0.9))
    ax.xaxis.set_tick_params(labelsize=int(fs*0.75))
    ax.yaxis.set_tick_params(labelsize=int(fs*0.75))

    if xticks is not None: ax.set_xticks(xticks)
    if yticks is not None: ax.set_yticks(yticks)
    
    
    if xlims is not None: ax.set_xlim(xlims)
    if ylims is not None: ax.set_ylim(ylims)
    

    if legend is not None:
        ax.legend(legend) if len(legend) else ax.legend(fontsize=int(fs*0.75))

    if xscale is not None: ax.set_xscale(xscale)
    if yscale is not None: ax.set_yscale(yscale)
    
    return


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
            S_real/=S_real.max()

        try:
            with open(path_exp+'_DWI_img.txt', 'r') as f:
                S_img = np.loadtxt(f)
                S_img/=S_img.max() 
        except FileNotFoundError:
            S_img = None

        try:
            with open(path_exp+'_DWI_intra.txt', 'r') as f:
                S_intra = np.loadtxt(f)
                S_intra/=S_intra.max() 
            with open(path_exp+'_DWI_extra.txt', 'r') as f:
                S_extra = np.loadtxt(f) 
                S_extra/=S_extra.max()
        except FileNotFoundError:
            S_intra, S_extra = None, None
        
    return S_real, S_img, S_intra, S_extra


def powder_averaged_S(S_, exp_proto):


    S_b_d   = np.zeros((b_vals_u.size, t_diff_u.size))

    for bb_idx, bb in enumerate(b_vals_u): 
        for tt_idx, tt in enumerate(t_diff_u):
            idx_ = np.where((b_vals==bb)&(t_diff==tt))[0]
            S_b_d[bb_idx, tt_idx] = S_[idx_].mean()

    return S_b_d 
    

def show_benchmark():


    # 2. Load signals 

    # 2.1 Spheres list
    exp_sphere_list     = "R_2_R_4_v_50_ICVF_0.57_gaussian_sphere_packing_orig"
    path_sphere_list    = os.path.join(PATH_OUTPUT_FOLDER, "sphere", exp_sphere_list)

    S_real_sphere_list, S_img_sphere_list, S_intra_sphere_list, S_extra_sphere_list = load_S(path_sphere_list, bin=False)
    S_powder_sphere_list = powder_averaged_S(S_real_sphere_list, exp_protocol)


    exp_sphere_list_cuda     = "R_2_R_4_v_50_ICVF_0.57_gaussian_sphere_packing_rep_04"
    path_sphere_list_cuda    = os.path.join(PATH_OUTPUT_FOLDER, "sphere", exp_sphere_list_cuda)

    S_real_sphere_list_cuda, _, _, _ = load_S(path_sphere_list_cuda, bin=False)
    S_powder_sphere_list_cuda = powder_averaged_S(S_real_sphere_list_cuda, exp_protocol)



    # 2.2 Spheres PLY
    exp_sphere_ply     = "R_2_R_4_v_50_ICVF_0.57_gaussian_sphere_packing"
    path_sphere_ply    = os.path.join(PATH_OUTPUT_FOLDER, "sphere", exp_sphere_list)

    S_real_sphere_ply, S_img_sphere_ply, S_intra_sphere_ply, S_extra_sphere_ply = load_S(path_sphere_ply, bin=False)
    S_powder_sphere_ply = powder_averaged_S(S_real_sphere_ply, exp_protocol)
    

    # 3. Visualization
    # 3.1 Raw signals
    
    f_sig, ax_sig = plt.subplots(1, 2, figsize=(FIG_W*2, FIG_W))
    ax_sig[0].plot(S_real_sphere_list, label="Full")
    ax_sig[0].plot(S_intra_sphere_list, label="Intra")
    ax_sig[0].plot(S_extra_sphere_list, label="Extra")
    
    ax_sig[0].plot(S_real_sphere_list_cuda, label="Full - Cuda")

    ax_sig[1].plot(S_real_sphere_ply, label="Full")
    ax_sig[1].plot(S_intra_sphere_ply, label="Intra")
    ax_sig[1].plot(S_extra_sphere_ply, label="Extra")
    

    add_plot_info(ax_sig[0], title="Sphere list", xlabel="Exp")
    add_plot_info(ax_sig[1], title="Sphere PLY", xlabel="Exp")
    

    plt.legend()
    f_sig.savefig(os.path.join(PATH_SAVE, "signals_raw.png"))

    # 3.2 Powder averaged signals
    f_sig_powder, ax_sig_powder = plt.subplots(1, 2, figsize=(FIG_W*2, FIG_W))
    ax_sig_powder[0].plot(b_vals_u, S_powder_sphere_list, label="Full")
    ax_sig_powder[0].plot(b_vals_u, S_powder_sphere_list_cuda, label="Full")
    
    ax_sig_powder[1].plot(b_vals_u, S_powder_sphere_ply, label="Full")
    
    add_plot_info(ax_sig_powder[0], title="Sphere list", xlabel=r"$b$")
    add_plot_info(ax_sig_powder[1], title="Sphere PLY", xlabel=r"$b$")

    
    plt.legend()
    f_sig_powder.savefig(os.path.join(PATH_SAVE, "signals_powder.png"))

    # 3.3 ADC
    


    return 

def main():

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


    # 1. Load scheme
    exp_protocol        = read_scheme(PATH_SCHEME)

    b_vals      = give_b(exp_protocol)
    b_vals_u    = np.unique(b_vals)
    t_diff      = exp_protocol[:, 4]
    t_diff_u    = np.unique(t_diff)


    show_benchmark()