# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:02:25 2019

@author: 345578
"""

#####################LICENSING####################################################### 
#All code here in is distributed under an OSS License © 2022. Triad National 
#Security, LLC. All rights reserved. This program was produced under U.S. 
#Government contract 89233218CNA000001 for Los Alamos National Laboratory 
#(LANL), which is operated by Triad National Security, LLC for the U.S. 
#Department of Energy/National Nuclear Security Administration. All rights in 
#the program are reserved by Triad National Security, LLC, and the U.S. 
#Department of Energy/National Nuclear Security Administration. The Government
# is granted for itself and others acting on its behalf a nonexclusive, 
#paid-up, irrevocable worldwide license in this material to reproduce, prepare 
#derivative works, distribute copies to the public, perform publicly and 
#display publicly, and to permit others to do so.

#And a BSD license: This program is open source under the BSD-3 License. 
#Redistribution and use in source and binary forms, with or without modification, 
#are permitted provided that the following conditions are met:

# Redistributions of source code must retain the above copyright notice, this 
# list of conditions and the following disclaimer. 2.Redistributions in binary 
# form must reproduce the above copyright notice, this list of conditions and 
# the following disclaimer in the documentation and/or other materials provided 
# with the distribution. 3.Neither the name of the copyright holder nor the 
# names of its contributors may be used to endorse or promote products derived 
# from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
########################End#########################################################

#IMAP-TDIA Code Base Zachary Robbins and Chonggang Xu




################Single-Model-MCMC.py###############
# This code was built to run an MCMC sampler for 4 the IMAP/TDIA models in parallel. 
# It was used to do the parameterization of the attack model to the field data. 
# It runs a single mcmc parameterizations with 4 parallel model runs (4 sites) for both a current and proposed move.
# Zachary Robbins 2021


## Load in packages 
import numpy as np
from scipy.stats import beta
from scipy.stats import lognorm
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
import matplotlib
from functools import partial 
import math
from joblib import Parallel, delayed ,cpu_count,Memory
from statsmodels.graphics.tsaplots import plot_acf
import WPB_Module  as WPB ### This is the WPB module from Fortran
import os
import gc
from joblib.externals.loky import get_reusable_executor
### Looking at the major fuction
print(WPB.beetlefit.__doc__)


##############################
### Setting up variables for the model run   
WPB_LUT=pd.read_csv('WPB_Inputs/WPB_Trees_6_20.csv')
WPB_LUT=WPB_LUT[WPB_LUT['Location'].isin(['UpperMed','UpperLow','LowerLow','LowerMed'])]
###update Lookup numbers
WPB_LUT.LUN=list(range(1,21))
#print(WPB_LUT)
locations=WPB_LUT.Location.unique()
### Temperature Max
Tmax_csv=pd.read_csv('WPB_Inputs/MI_Tmax_6_20.csv')
#print(list(Tmax_csv.columns))
Tmax_csv=Tmax_csv.rename(columns={'Unnamed: 0':"Unnamed",'Dates':"Dates",
                'Low_LowElShape_Out':"LowerLow",
                'Low_MedElShape_Out':"LowerMed",
                'Low_HighElShape_Out':"LowerHigh",
                'High_LowElShape_Out':"UpperLow",
                'High_MedElShape_Out':"UpperMed",
                'High_HighElShape_Out':"UpperHigh"})
Tmax_csv.sort_values(by=['Dates'],inplace=True)
### Temperature Min 
Tmin_csv=pd.read_csv('WPB_Inputs/MI_Tmin_6_20.csv')
Tmin_csv=Tmin_csv.rename(columns={'Unnamed: 0':"Unnamed",'Dates':"Dates",
                'Low_LowElShape_Out':"LowerLow",
                'Low_MedElShape_Out':"LowerMed",
                'Low_HighElShape_Out':"LowerHigh",
                'High_LowElShape_Out':"UpperLow",
                'High_MedElShape_Out':"UpperMed",
                'High_HighElShape_Out':"UpperHigh"})
Tmin_csv.sort_values(by=['Dates'],inplace=True)
### Drought
Drought_csv=pd.read_csv('WPB_Inputs/MI_SPI_6_20.csv')
Drought_csv=Drought_csv.rename(columns={'Unnamed: 0':"Unnamed",'Dates':"Dates",
                'LE_LL':"LowerLow",
                'HE_LL':"LowerMed",
                'LE_HL':"UpperLow",
                'HE_HL':"UpperMed"})



Drought_csv.sort_values(by=['Dates'],inplace=True)
mask = (Drought_csv['Dates'] > "2001-10-01") & (Drought_csv['Dates'] <= "2018-12-31")
Drought_csv=Drought_csv.loc[mask]

CWD=Drought_csv
Dates=pd.DataFrame(Tmax_csv.Dates)


#This method of manually build an mcmc sampler derived from
# https://github.com/Joseph94m/MCMC/blob/master/MCMC.ipynb

##### Function set 
# These functions nest to run the MCMC model 
### This a formatting tool to get the Fortran to accept the input from numpy. 
def p2f(x1,dtype): 
  x1=np.array(x1,dtype=dtype,order='F')
  return x1

### This was the transition model used to fit the covaraince structure, running with this as the proposed 
### step rather than the normal covariance matrix. 

#transition_model = lambda x2: [np.random.normal(x2[0],0.002,(1,))[0].round(decimals=5), ### SD of distribution 
#                                               ##### Draw one Shape 1d take the first (Lots of things for just one )
#                              np.random.normal(x2[1],0.015,(1,))[0].round(decimals=5),
#                              np.random.normal(x2[2],0.001,(1,))[0].round(decimals=5),
#                              np.random.normal(x2[3],0.002,(1,))[0].round(decimals=5),
#                              np.random.normal(x2[4],15,(1,))[0].round(decimals=0),
#                              x2[5],
#                              np.random.normal(x2[6],.001,(1,))[0].round(decimals=5)] ##


def Bounce(xR):
   # The Bounce function takes proposed moves that are outside the realm 
   # of biological feasibility and mirrors them back on the other-side of a boundary
   ### Intercept
    if(xR[1] >=14.0):
        xR[1]=14.0-(xR[1]-14.0)    
    if(xR[1] <=10.0):
        xR[1]=10.0+(10.-xR[1])
		### Drought influence
    if(xR[2] >=1.55):
        xR[2]=1.55-(xR[2]-1.55)
		### Size influence
    if(xR[3] <=.50):
        xR[3]=.50+(.50-xR[3])
		### Initial Population
    if(xR[4] <=221):
        xR[4]=221+(221-xR[4])
    if(xR[4] >873):
        xR[4]=873-(xR[4]-873)
		### Size preference (0-1)
    if(xR[6] >1):
        xR[6]=1-(xR[6]-1)
    if(xR[6] <0):
        xR[6]=0+(0-xR[6])    
   
    return xR

def TemperatureParser(LUN,Tmin_csv,Tmax_csv,CWD,MPB_LUT,startdate,enddate):
    ### Given a look up number for sites (0-3) this function retrieves 
	### the climate variables and the number of trees for initial conditions at a site. 
    Dates=Tmax_csv['Dates']
    LU=MPB_LUT[MPB_LUT.Location==locations[LUN]]
   
    L=LU.Location
    Tmax_One=Tmax_csv[L]
    Tmin_One=Tmin_csv[L]
    CWD_One=CWD[L]
    Temp_LUT=pd.concat([Tmax_One.reset_index(drop=True),Dates.reset_index(drop=True)],axis=1)
    Temp_LUT=pd.concat([Temp_LUT.reset_index(drop=True),Tmin_One.reset_index(drop=True)],axis=1)
    Temp_LUT=pd.concat([Temp_LUT.reset_index(drop=True),CWD_One.reset_index(drop=True)],axis=1)
    Temp_LUT=Temp_LUT.sort_values('Dates')
    Temp_LUT['Dates'] = pd.to_datetime(Temp_LUT['Dates'])
    mask = (Temp_LUT['Dates'] > startdate) & (Temp_LUT['Dates'] <= enddate)
    Temp_LUT = Temp_LUT.loc[mask]
    ntgeq317=np.int(LU.Initial_above_20)
    ntgeq00=np.int(LU.Initial_below_20)
    #SubYears.iloc[:,0].isnull().values.any()
    return(Temp_LUT.iloc[:,2],Temp_LUT.iloc[:,0],Temp_LUT.iloc[:,3],ntgeq317,ntgeq00)

 
def Insect_Modeling(Tmin_csv,Tmax_csv,CWD,WPB_LUT,LUN,Xbind,startdate,enddate):
 ## This function parses out the variables from Xbind and parses out which site the 
  ## is supposed to be simulated. Then it fetches the climate variable and the number of trees (Temperature parser).
  ## It lastly runs the fortran model with the appropriate variables and stores the trees for evaluation at the level above. 
    r1=Xbind[0]
    x0=Xbind[1] 
    x1=Xbind[2]
    x2=Xbind[3]
    tparents=Xbind[4]
    size_factor=Xbind[6]
    LUN1=LUN
    ### Pull the correct Temperatures 
    SubMin,SubMax,CWDvec,ntgeq317,ntgeq00 =TemperatureParser(LUN1,
                                                             Tmin_csv,
                                                             Tmax_csv,
                                                             CWD,WPB_LUT,
                                                             startdate,
                                                             enddate)
    #### Run the model 
    NGEQ20vec,NGEQ00vec,Flight,Fec,Eggs,L1,L2,Pout,Ten,Adults= WPB.beetlefit(p2f(SubMin,"cdouble"),#TempMin
                p2f(SubMax,"cdouble"),#TempMax
                p2f(r1,"float"), #phi1
                p2f(x0,"float"), # intercept
                p2f(x1,"float"), # Drought influence
                p2f(CWDvec,"cdouble"), # Drought Vector
                p2f(x2,"float"), # Size influence 
                p2f(size_factor,"float"),# Aggregation factor
                p2f(ntgeq317,"float"), # Large Trees
                p2f(ntgeq00,"float"), # Small Trees 
                p2f(tparents,"float")) # Initial parents 
	## Set Nas to zero (Nas result when the population is zero or less than zero			
    LargeOut=NGEQ20vec
    LargeOut[np.isnan(LargeOut)]=0
    SmallOut=NGEQ00vec
    SmallOut[np.isnan(SmallOut)]=0
    ##### Sort data out to match up with the data we are testing 
    BigTreesOut=np.array([LargeOut[5*365],
                       LargeOut[6*365],LargeOut[7*365],LargeOut[8*365],
                       LargeOut[9*365],LargeOut[10*365],LargeOut[11*365],
                       LargeOut[12*365],LargeOut[13*365],LargeOut[14*365],
                       LargeOut[15*365],LargeOut[16*365],LargeOut[17*365]])
    SmallTreesOut=np.array([SmallOut[5*365],
                       SmallOut[6*365],SmallOut[7*365],SmallOut[8*365],
                       SmallOut[9*365],SmallOut[10*365],SmallOut[11*365],
                       SmallOut[12*365],SmallOut[13*365],SmallOut[14*365],
                       SmallOut[15*365],SmallOut[16*365],SmallOut[17*365]])
	### clean out parallelized 
    r1=[]
    x0=[] 
    x1=[]
    x2=[]
    tparents=[]
    size_factor=[]

    gc.collect()
    return(LUN,BigTreesOut,SmallTreesOut)
   


#Computes the likelihood of the data given a sigma (new or current) according to equation (2)
### In this version this is where all the parallelized work is happening 

def manual_log_like_normal(WPB_LUT,Tmax_csv,Tmin_csv,CWD,xbind):  
  ## This function calculates the log likelihood for each of the proposed steps for each of the four testing runs. 
  ### Set up containers 
  Outputmodel=[]
  OutputY_obs=[]
  ### Dates control which climate gets loaded in.
  startdate='10-01-2001'
  enddate='12-31-2018'
   ### This feeds the function all the stable elements. The partial function stores all the variables in a 
  ### single function to limit redundancy. Insect Modeling is prime function, and the rest is stable variables.
  ### Insect modeling sorts out the sites internally based on the look up number(LUN) that is passed in the next line. 
  PartIM=partial(Insect_Modeling,Tmin_csv=Tmin_csv,Tmax_csv=Tmax_csv,
   WPB_LUT=WPB_LUT,Xbind=xbind,CWD=CWD,
   startdate=startdate,enddate=enddate)   

  #### Each position is structured as four parallel runs, one for each site
  Run1,Run2,Run3,Run4=Parallel(n_jobs=4)(delayed(PartIM)(LUN=i) for i in  range(0,4))

  ## Create on output	
  Outputmodel1=np.concatenate([Run1[1],Run2[1],Run3[1],Run4[1],Run1[2],Run2[2],Run3[2],Run4[2]])
  OutputY_obs=np.concatenate([np.array(WPB_LUT.iloc[0,4:17]),np.array(WPB_LUT.iloc[1,4:17]),
                            np.array(WPB_LUT.iloc[2,4:17]),np.array(WPB_LUT.iloc[3,4:17]),
                            np.array(WPB_LUT.iloc[0,17:30]),np.array(WPB_LUT.iloc[1,17:30]),
                            np.array(WPB_LUT.iloc[2,17:30]),np.array(WPB_LUT.iloc[3,17:30])])
    
  
  OutputY_obs=OutputY_obs.astype(np.float32)  
  Outputmodel=Outputmodel1.astype(np.float32) 
  Outputmodel=Outputmodel.astype(np.float64) 
  Outputmodel[Outputmodel<.1]=.1
  ## Calculate the log liklihood. 
  sigma=np.array(OutputY_obs*.0983)
  mean=np.log((OutputY_obs**2)/(np.sqrt((OutputY_obs**2)+(sigma**2))) )
  sigmaprime=np.sqrt(np.log(1+((sigma**2)/(OutputY_obs**2))))
  LL = sum(
        -np.log(
          (np.exp((-(np.log(Outputmodel) - mean)**2) 
                  / (2 * sigmaprime**2)))
       / (Outputmodel*sigmaprime * np.sqrt(2 * np.pi))))

  print(LL)
  xbind=[]
  ## Clean out the parallel space 
  get_reusable_executor().shutdown(wait=True)
  del(PartIM)
  del(Run1)
  del(Run2)
  del(Run3)
  del(Run4) 
  gc.collect()
  ### return the log-liklihood score and the raw model results. 
  return(LL,Outputmodel)



def acceptance(x_likelihood, x_new_likelihood):
#	 Defines whether to accept or reject the new sample
    if x_new_likelihood<x_likelihood:
        return True
    else:
        accept=np.random.uniform(0,1)
        # Since we did a log likelihood, we need to exponentiate in order to compare to the random number
        return (accept < np.exp(x_likelihood-x_new_likelihood))

def metropolis_hastings(likelihood_computer, cov, param_init,iterations,WPB_LUT,Tmax_csv,Tmin_csv,CWD,acceptance_rule):

    #### This function manages the testing of each model using a metropolis hastings step 

    accepted = []
    rejected = []
    LLlist=[]
    alllist=[]
    Treeslist=[]
    start = time.time()
	### Set up initial position 
    x = param_init ## Proposed start
    x=Bounce(x)
    Xbind=np.array(x)
    print(Xbind)
    ### Calculate the log-liklihood of the initial position 
    x_lik,oldtrees= manual_log_like_normal(WPB_LUT,Tmax_csv,Tmin_csv,CWD,Xbind)
    print(oldtrees)
    for i in range(iterations):
	# For each iteration of testing 
        if((i/50) % 1 ==0):
            print("#################################"+str(i)+"###########################################")            
        #Clocking 
        print(i)
		# Propose a new step 
        x_new = np.random.multivariate_normal(x,cov, size=1) .tolist()[0]
        x_new=Bounce(x_new)                      
        #print(x_new)
        # Test the step and get a log-liklihood 
        Xbind=np.array(x_new)
        x_new_lik,new_trees= manual_log_like_normal(WPB_LUT,Tmax_csv,Tmin_csv,CWD,Xbind)
        print(x_new_lik,new_trees)
        # Gather all the outputs for analysis
        Treeslist.append(oldtrees)
        Treeslist.append(new_trees)
        # Check whether to accep the new move 
        if (acceptance_rule(x_lik,x_new_lik)): 
			# if true update the proposed variables 
            x = x_new
            x_lik=x_new_lik
            accepted.append(x_new)
            LLlist.append(x_new_lik)
            alllist.append(x_new)
            oldtrees=new_trees
            print("accepted")
            
        else:
			# If false store the scores and moves 
            rejected.append(x_new) 
            LLlist.append(x_lik)
            alllist.append(x_new)
            print("rejected")
    end = time.time()
    print("Run "+str(end - start)+" seconds")   
	## Once completed store a vector of the outputs. 
    return np.array(accepted), np.array(rejected),np.array(LLlist),np.array(alllist),Treeslist

### Running portion 
### For this version the proposed move is calculated by the multi-normal distribution with a covariance structure taken
### from runs in which a normal distribution was proposed with the parameter stated above. This vector of accepted moves 
### is then used to create a covariance matrix by which to propose new moves here. 

x=[]
CovSet=pd.read_csv("Outputs/Accept_Covariance.csv")
cov=np.cov(CovSet.transpose())

accepted, rejected,LLoutput,alllist,Treeslist= metropolis_hastings(manual_log_like_normal,### Computes Liklihood
                                         cov,
                                         [0.2009,10.03,1.545 ,0.506,840 ,2,0.051], ### Initial postion
                                       15000, ###Iterations
                                         WPB_LUT, ### LUT information including LUN for parproc
                                         Tmax_csv,### All Tmax
                                         Tmin_csv,### All Tmin 
                                         CWD,
                                         acceptance) ##Acceptatnce rule (Metropolis)

Trees=pd.DataFrame(Treeslist)
Trees1=np.array(Trees.iloc[::2])
np.savetxt("Accept_Practice_CovR2_1115.csv",accepted,delimiter=",")
np.savetxt("Reject_Practice_CovR2_1115.csv",rejected,delimiter=",")
np.savetxt("AllPOWPractice_CovR2_1115.csv",alllist,delimiter=",")
np.savetxt("TreesPOWPractice_CovR2_1115.csv",Treeslist,delimiter=",")
np.savetxt("LL_Practice_CovR2_1115.csv",LLoutput,delimiter=",")

