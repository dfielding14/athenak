#!/usr/bin/env python

import numpy as np
import subprocess
import argparse
import re

parser = argparse.ArgumentParser()
args = vars(parser.parse_args())
locals().update(args)

if __name__ == '__main__':
    # resolution
    Nx_vals = [256]
    # Mach number
    Ms = [0.5]
    # density contrasts
    chis = [1e2]
    # xi = tshear/tcool
    xis = [1e1,3e1,1e2,3e2,1e3]
    # riemann solver
    rsolver = ["hllc"]
    # options for hydrodynamics
    #hyd_ops = ["hydro"]

    # choices related to resolution
    NX_MBs = {64:64, 128:64, 256:128, 512:128}
    NNODES = {64:1 , 128:1 , 256:1  , 512:8  }
    NTASKS = {64:1 , 128:8 , 256:8  , 512:64 }
    NT_PN  = {64:1 , 128:8 , 256:8  , 512:8  }
    NCPU_PT= {64:8,  128:8 , 256:8  , 512:8  }
    NGPU_PT= {64:1 , 128:1 , 256:1  , 512:1  }
    TIMES  = {64 :"00:30:00" , 128:"01:00:00",\
              256:"18:00:00" , 512:"36:00:00"}

    # specify the paths to templates and where results will be saved
    input_path = "/mnt/home/llancaster/athenak-RM/inputs/hydro/"
    slurm_path = "/mnt/home/llancaster/athenak-RM/scripts/"
    # specify templates
    input_template = "TRML_template.athinput"
    slurm_template = "rusty.trml_template"

    for Nx in Nx_vals:
        for M in Ms:
            for chi in chis:
                for xi in xis:
                    for rs in rsolver:
                        # Set the output file name
                        athinput_out = "TRML_N%d_%s_chi%1.0e_xi%1.0e.athinput"%(Nx,rs,chi,xi)
                        # Set dictionary for string substitutions
                        subs = {}
                        ## start with resolution
                        subs["nx1"] = str(Nx)
                        subs["nx2"] = str(Nx)
                        subs["nx3"] = str(int(1.5*Nx))
                        subs["nx1_mb"] = str(NX_MBs[Nx])
                        subs["nx2_mb"] = str(NX_MBs[Nx])
                        subs["nx3_mb"] = str(int(3*NX_MBs[Nx]//2))

                        ## TRML parameters
                        vrel = M*np.sqrt(5./3)
                        subs["vrel"] = str(vrel)
                        tshear = 1.0/vrel
                        # simulation run for 30 shear times
                        tsim = 30*tshear
                        subs["tlim"] = str(tsim)
                        subs["chi"] = str(chi)
                        subs["xi"] = str(xi)
                        tcool = tshear/xi
                        subs["tcool"] = str(tcool)

                        ## various other problem paramters
                        subs["dx_smooth"] = str(1./(2*Nx))
                        subs["dt_hst"] = str(tsim/1e4)
                        subs["dt_output_hydro"] = str(tsim/1e2)
                        subs["dt_slice"] = str(tsim/1e3)

                        # specify the riemann solver
                        subs["rsolver"] = rs

                        ## Read template, make substitutions, and write a file
                        ## for the input file
                        with open(input_path + input_template, 'r') as f:
                            tmp = f.read()
                        for k, v in subs.items():
                            tmp = re.sub(r'@{0}@'.format(k), v, tmp)
                        with open(input_path + athinput_out, 'w') as f:
                            f.write(tmp)

                        ## Now deal with creating the submission script #####################

                        slurm_sub = {}
                        slurm_out = "rusty.run_trml_N%d_%s_chi%1.0e_xi%1.0e.slurm"%(Nx,rs,chi,xi)

                        slurm_sub["sim_name"] = "TRML_N%d_%s_chi%1.0e_xi%1.0e"%(Nx,rs,chi,xi)
                        slurm_sub["sim_input"] = athinput_out
                        slurm_sub["nnodes"] = str(NNODES[Nx])
                        slurm_sub["ntasks"] = str(NTASKS[Nx])
                        slurm_sub["ntasks_per_node"] = str(NT_PN[Nx])
                        slurm_sub["time"] = TIMES[Nx]
                        slurm_sub["cpus_per_task"] = str(NCPU_PT[Nx])
                        slurm_sub["gpus_per_task"] = str(NGPU_PT[Nx])

                        ## Read template, make substitutions, and write a file
                        ## for the input file
                        with open(slurm_path + slurm_template, 'r') as f:
                            tmp = f.read()
                        for k, v in slurm_sub.items():
                            tmp = re.sub(r'@{0}@'.format(k), v, tmp)
                        with open(slurm_path + slurm_out, 'w') as f:
                            f.write(tmp)
