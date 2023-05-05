# 1. Setting up environment on Cedar

## 1.0 Prerequisite

All required packages can be loaded with `module`:
```bash
$ module load qt/5.15.2 gcc/9.3.0 StdEnv/2020   root/6.26.06  eigen/3.3.7
```

## 1.1 Simulation


## 1.2 Tracker

The Tracker code is located at: https://github.com/EdmondRen/MATHUSLA-Kalman-Algorithm

**Build**
```bash
$ cd tracker/build
$ cmake ..
$ make 
```

**Run**   
At this point, the tracker executable is available in the /build/ directory. Note that the /build/ directory MUST be placed in the /tracker/ directory.
```bash
$ ./tracker path_to_input_file path_to_write_output 
```

## 1.3 Analysis

We are using jupyter interface: https://jupyterhub2.cedar.computecanada.ca/

In JupyterHub, launch a terminal. Run
`/project/def-mdiamond/soft/singularity/cdmsfull_V04-07-00/cdmskernel_install.sh`

Open a Launcher, there should be a new kernel titled "Singularity V04-07". Clicking on it will open a notebook with this kernel.

# 2. Analysis outline

## 2.1 Single track study