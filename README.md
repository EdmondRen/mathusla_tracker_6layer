# 1. Setting up environment on Cedar

## 1.0 Prerequisite

Most of the required packages can be loaded with `module`. Dlib is not installed by default. You can either complile it or use the one under my home directory, by appending it to `PATH`.
```bash
$ module load qt/5.15.2 gcc/9.3.0 StdEnv/2020   root/6.26.06  eigen/3.3.7
$ PATH=$PATH:/project/def-mdiamond/tomren/mathusla/dlib-19.24/install
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
At this point, the tracker executable is available in the ./build/ directory. Note that the ./build/ directory MUST be placed in the /tracker/ directory.
```bash
$ ./tracker path_to_input_file path_to_write_output 
```

## 1.3 Analysis

We are mostly using python to do analysis. This repository contains python functions for reading/visualizing the data. There are also many jupyter notebooks that contains some analysis code. Jupyter notebook is a convenient way of running interactive data analysis. 

**Set up jupyter on cedar**   

Jupyter interface on cedar: https://jupyterhub2.cedar.computecanada.ca/   
You need to setup the python kernel before using the jupyter for the first time. 
In JupyterHub, launch a terminal. Run
`/project/def-mdiamond/soft/singularity/cdmsfull_V04-07-00/cdmskernel_install.sh`   
Open a Launcher, there should be a new kernel titled "Singularity V04-07". Clicking on it will open a notebook with this kernel.

# 2. Analysis outline

## 2.1 Single track study

# 3. File discription
## 3.1 Analysis
## 3.2 Python functions