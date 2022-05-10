C21018 IMAP-TDIA Code Base

#####################LICENSING#######################################################
All code here in is distributed under an OSS License 
© 2022. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.

And a BSD license: 
This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
2.Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution. 
3.Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
########################End#########################################################


IMAP-TDIA Code Base
Zachary Robbins and Chonggang Xu

With this code, we hope to make openly available the tools necessary to conduct the analysis we conducted in the Sierra Nevada for the western pine beetle. These scripts provide examples as to how to calculate the necessary model inputs. This includes the median rate of development for the individual insect stages and the mortality rate of trees killed by beetles from the FIA database. It contains the scripts used for the analysis on drought, and the methods for calculating the drought driver needed to run the TDIA model as well as the code necessary to create a historic climate by subtracting the difference in historic and current monthly normal from the daily climate driver. These inputs can then be used with the FORTRAN code included for implementing the IMAP/TDIA model within python (3.x) and instructions on how to initiate it. It also contains the Markov chain Monte Carlo (MCMC) method used to fit the model to the data on tree mortality and the method for validating the model. This structure within python is how the model was implemented, and data to initiate these model runs will be included with code. When possible, notebooks are also included in HTML output, so that they can be seen running with the external data. Data sources are documented throughout, though data (such as publicly available FIA, or DAYMET data) is not included due to size links to its download are included. Other data may be omitted due to contributor priority.
 
 List of Scripts: 

1.1 Calculate_Forest_Mortality.R (Html output included):
               This script is used to calculate the FIA conditions for initiating the model and for mortality comparison. This model subsets the FIA data to latitude and elevation parameters for each model within the study area. It then calculates the percentage of mortality observed during each year. It calculates the initial conditions of the model from the mean observed ponderosa pine density for the years of 2005-2006. It then applies that percentage of mortality to the initial density for each area for each time step. This results in a .csv that was used the MCMC processes and analysis.
 
1.2 Comparing_Drought_Drivers.R (Html output included):
               This script calculates the Precipitation -Evapotranspiration, Standard Precipitation Index, and Palmers Drought Severity Index derived drought drivers for the TDIA model that could be used within the model for the four sites we modeled. It also compares those drivers to the mortality experienced by the trees for each year and calculates an R2. 

1.3 Median_Insect_Growth.IPYNB (Html output included):
               This script uses the python package Pymc3 to fit the growth curve required for each stage of insect development from published rates of development in relation to temperature. It fits each growth curve using a Hamiltonian sampler. Graphs are produced for the runs and the resulting curves. 

1.4 Historical_Climate_WPB.RMD (Html output included):
 	This script calculates a minimum and maximum temperature vector for the IMAP/TDIA model that reflects the difference in monthly means between the two time periods analyzed. To do this we will first compare the PRISM monthly historic means, with the monthly means from the observed temperatures used to run the model for each of the four sites. Then we will use the monthly difference to calculate a new daily driver that removes this difference, thus creating the historic.1.5 WPB_Submission.F90 

1.5 WPB_Submission.F90:
          The FORTRAN code to implement the IMAP/TDIA. This includes the insect growth and population model (IMAP) and the tree attack model (TDIA) and implements them together with climate drivers and produces the number of trees killed at each time step, as well as insect populations at each development stage. Included in this FORTRAN code are basic instructions as to how to implement FORTRAN code in python using the f2py module. This code can be run on a variety of FORTRAN implementation platforms, however, for our research, we primarily wrapped the FORTRAN code as a python package. This often has to be done unique to a given Python version and setup, so instructions and some links are included as to how to implement this. 

1.6 WPB_Test_Sub.cp36-win_amd64:
         An f2py wrapper for python 3.6 on windows using mgwin64 bit. This can be implemented on certain python 3.6 builds on windows, by calling import WPB_Test_Sub. The file must be in the same location as the script calling it. If you cannot implement this read the WPB_submission.F90 instructions as to how to build a python wrapper for the FORTRAN code for your given python build.

1.6 Single-Model-MCMC.py:
         This code was built to run an MCMC sampler for 4 the IMAP/TDIA models in parallel. It was used to do the parameterization of the attack model to the field data. It runs a single MCMC parameterization with 4 parallel model runs (4 sites) for both a current and proposed move. It also serves as an example of how to implement the IMAP/TDIA in python. 

1.7 Multi-Model-MCMC.py:
 	This code was built to run an MCMC sampler for 4 of the IMAP/TDIA models in parallel. It was used to do the validation runs (removing one site from each parameterization). It runs four MCMC parameterizations independent and parallel, using a multivariate_normal distribution.


