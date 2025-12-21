# Casadi and Ipopt Installation Instructions on Ubuntu
## 1. Installation of Ipopt
Casadi's nonlinear optimization functionality relies on the Ipopt library, so Ipopt must be installed first.

### 1.1 Install Dependencies
Run the following command to install the required dependencies:
```bash
sudo apt-get install gcc g++ gfortran git patch wget pkg-config liblapack-dev libmetis-dev libblas-dev
```

### 1.2 Create a Dedicated Folder
Create a folder to store all Ipopt-related files for easy management:
```bash
mkdir ~/Ipopt_pkg
cd ~/Ipopt_pkg
```

### 1.3 Install ASL
```bash
git clone https://github.com/coin-or-tools/ThirdParty-ASL.git
cd ThirdParty-ASL
sudo ./get.ASL
sudo ./configure
sudo make
sudo make install
cd ..
```

### 1.4 Install HSL
```bash
git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
cd ThirdParty-HSL
```

Next, you need to download the coinhsl file and extract it to the ThirdParty-HSL directory. The download method for the coinhsl file is as follows:
To download the coinhsl file, you need to apply for a license. The license application URL is: [STFC IP Store](https://www.stfc.ac.uk/ip-store/). The available application types are `Coin-HSL` or `Coin-HSL Archive`.
- The `Coin-HSL` type is only available for free application by academic institutions; other institutions need to pay and go through an approval process, which takes a long time.
- The application for the `Coin-HSL Archive` type is approved immediately, so you can apply for this if you need it urgently.

After the application is approved, extract the coinhsl file in the ThirdParty-HSL directory and rename the extracted folder to `coinhsl`. Then, execute the following commands in the ThirdParty-HSL directory:
```bash
sudo ./configure
sudo make
sudo make install
cd ..
```

### 1.5 Install MUMPS
```bash
git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git
cd ThirdParty-Mumps
sudo ./get.Mumps
sudo ./configure
sudo make
sudo make install
cd ..
```

### 1.6 Install Ipopt
```bash
git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
mkdir build
cd build
sudo ../configure
sudo make
sudo make test
sudo make install
```

### 1.7 Configure Environment Variables
```bash
cd /usr/local/include
sudo cp coin-or coin -r
sudo ln -s /usr/local/lib/libcoinmumps.so.3 /usr/lib/libcoinmumps.so.3
sudo ln -s /usr/local/lib/libcoinhsl.so.2 /usr/lib/libcoinhsl.so.2
sudo ln -s /usr/local/lib/libipopt.so.3 /usr/lib/libipopt.so.3
```

If no errors occur during the above steps, Ipopt is successfully installed.

## 2. Installation of Casadi
The installation and testiting of Casadi is a simple process. The following script can be used to complete the installation with one click. For installing a higher version, the method is similar, or you can refer to the installation instructions in the [Casadi official GitHub repository](https://github.com/casadi/casadi).

```bash
#!/usr/bin/env bash

# Fail on first error.
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo -e "\033[40;32m${DIR} \033[0m"

# Download
wget https://github.com/casadi/casadi/releases/download/3.5.5/casadi-3.5.5-1.tar.gz
tar -zxvf casadi-3.5.5-1.tar.gz
echo -e "\033[40;32mdownload finish \033[0m"

cd casadi-3.5.5.1
mkdir build && cd build
cmake .. -DWITH_IPOPT=ON -DWITH_EXAMPLES=OFF
make -j4
sudo make install
sudo ldconfig

# Clean up.
sudo apt-get clean && sudo rm -rf /var/lib/apt/lists/*
sudo rm -fr casadi-3.5.5-1.tar.gz casadi-3.5.5.1
```

## 3. Casadi Installation Test
After the installation is complete, test whether Casadi is installed successfully. For convenience, the test program is packaged into a compressed package. [Click here to download the compressed package](
https://github.com/lzlbadguy/TriPlanner/raw/main/Casadi_test_package/casadi_test.tar.gz
) (Note: Replace with the actual download link of the test package) for testing.

After downloading, extract and enter the folder. The file directory is as follows:
```
.
├── CMakeLists.txt
└── src
    └── casadi_test.cpp
```

Execute the following commands to compile the test file:
```bash
mkdir build
cd build/
cmake ..
make
```

After successful compilation, execute the `./casadi_test_node` command in the build folder to run the test program. If Casadi is installed correctly, the following output will appear in the terminal:

```
casadi_test
x:[x_0, x_1, x_2]
p:[p_0, p_1]
f:((sq(x_0)+sq(x_1))+sq(x_2))
g:[((((6*x_0)+(3*x_1))+(2*x_2))-p_0), ((((p_1*x_0)+x_1)-x_2)-1)]

******************************************************************************
This program contains Ipopt, a library for large-scale nonlinear optimization.
 Ipopt is released as open source code under the Eclipse Public License (EPL).
         For more information visit https://github.com/coin-or/Ipopt
******************************************************************************

This is Ipopt version 3.14.16, running with linear solver ma27.

Number of nonzeros in equality constraint Jacobian...:        6
Number of nonzeros in inequality constraint Jacobian.:        0
Number of nonzeros in Lagrangian Hessian.............:        3

Total number of variables............................:        3
                     variables with only lower bounds:        3
                variables with lower and upper bounds:        0
                     variables with only upper bounds:        0
Total number of equality constraints.................:        2
Total number of inequality constraints...............:        0
        inequality constraints with only lower bounds:        0
   inequality constraints with lower and upper bounds:        0
        inequality constraints with only upper bounds:        0

iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
   0  4.5100000e-02 3.63e+00 4.11e-01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  5.8681488e-01 0.00e+00 1.95e+00  -1.0 3.91e-01    -  3.37e-01 1.00e+00h  1
   2  5.9327019e-01 1.11e-16 1.00e-06  -1.0 1.32e-02    -  1.00e+00 1.00e+00f  1
   3  5.6004673e-01 0.00e+00 4.92e-02  -2.5 8.93e-02    -  9.26e-01 1.00e+00f  1
   4  5.5264341e-01 0.00e+00 2.83e-08  -2.5 4.42e-02    -  1.00e+00 1.00e+00f  1
   5  5.5114453e-01 8.88e-16 1.50e-09  -3.8 2.36e-02    -  1.00e+00 1.00e+00f  1
   6  5.5102559e-01 1.78e-15 1.50e-09  -3.8 7.16e-03    -  1.00e+00 1.00e+00f  1
   7  5.5102042e-01 8.88e-16 1.84e-11  -5.7 1.77e-03    -  1.00e+00 1.00e+00f  1
   8  5.5102041e-01 8.88e-16 2.51e-14  -8.6 6.77e-05    -  1.00e+00 1.00e+00h  1
   9  5.5102041e-01 0.00e+00 9.17e-15  -9.0 9.29e-08    -  1.00e+00 1.00e+00h  1

Number of Iterations....: 9

                                   (scaled)                 (unscaled)
Objective...............:   5.5102040816326525e-01    5.5102040816326525e-01
Dual infeasibility......:   9.1685374877122625e-15    9.1685374877122625e-15
Constraint violation....:   0.0000000000000000e+00    0.0000000000000000e+00
Variable bound violation:   0.0000000000000000e+00    0.0000000000000000e+00
Complementarity.........:   9.0911698984215914e-10    9.0911698984215914e-10
Overall NLP error.......:   9.0911698984215914e-10    9.0911698984215914e-10

Number of objective function evaluations             = 10
Number of objective gradient evaluations             = 10
Number of equality constraint evaluations            = 10
Number of inequality constraint evaluations          = 0
Number of equality constraint Jacobian evaluations   = 10
Number of inequality constraint Jacobian evaluations = 0
Number of Lagrangian Hessian evaluations             = 9
Total seconds in IPOPT                               = 0.001

EXIT: Optimal Solution Found.
      solver  :   t_proc      (avg)   t_wall      (avg)    n_eval
       nlp_f  |   4.00us (400.00ns)   4.15us (415.40ns)        10
       nlp_g  |  14.00us (  1.40us)   9.50us (950.10ns)        10
    nlp_grad  |   2.00us (  2.00us)   1.40us (  1.40us)         1
  nlp_grad_f  |  14.00us (  1.27us)  11.14us (  1.01us)        11
  nlp_hess_l  |   5.00us (555.56ns)   4.56us (506.78ns)         9
   nlp_jac_g  |   7.00us (636.36ns)   5.98us (543.18ns)        11
       total  |   1.41ms (  1.41ms)   8.76ms (  8.76ms)         1
--------------------------------
Optimal solution for p = [5, 1]:
                   Objective: 0.55102
             Primal solution: [0.632653, 0.387755, 0.0204082]
           Dual solution (x): [-1.43695e-09, -2.3445e-09, -4.45467e-08]
           Dual solution (g): [-0.163265, -0.285714]
```

Congratulations! Casadi has been successfully installed.

## 4. A Common Issue During Ipopt Installation
During the configuration process, the following error may occur:
```
checking for function ma27ad_ in -L/usr/local/lib -lcoinhsl ... no checking for function ma27ad in -L/usr/local/lib -lcoinhsl ... no checking for function MA27AD_ in -L/usr/local/lib -lcoinhsl ... no checking for function MA27AD in -L/usr/local/lib -lcoinhsl ... no checking for function ma27ad__ in -L/usr/local/lib -lcoinhsl ... no checking for function ma27ad_ in -L/usr/local/lib -lcoinhsl ... no checking for function MA27AD__ in -L/usr/local/lib -lcoinhsl ... no checking for function MA27AD_ in -L/usr/local/lib -lcoinhsl ... no configure: error: Provided package HSL is not working or does not contain MA27. See config.log for details on failed checks.
```

This error indicates that HSL is not installed correctly. Perform the following operations in the HSL compilation directory:

### 4.1 Clean Up Previous Configurations
```bash
make distclean
```
This command will delete the previously generated configuration files and object files to prepare for reconfiguration.

### 4.2 Re-run Configure
```bash
./configure LIBS="-ldl"
```
In this command, we explicitly specify the `-ldl` library to resolve potential linking issues during configuration and ensure that the correct version of the dynamic link library is included during linking.

### 4.3 Check Configuration Results
Check the output of the `configure` command to ensure there are no errors such as `unrecognized command line option`.

### 4.4 Continue Compilation
```bash
make
make install
```

This document records the successful installation process during the learning phase for sharing purposes.