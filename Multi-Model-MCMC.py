# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:02:25 2019

@author: 345578
"""
################Multi-Model-MCMC.py###############
# This code was built to run an MCMC sampler for 4 of the IMAP/TDIA models in parallel. 
# It was used to do the validation runs (removing one site from each parameterization). 
# It runs four mcmc parameterizations independent and parallel, using a multivariate_normal distribution.
# Zachary Robbins 2021

## Load in packages 
import numpy as np
from scipy.stats import beta
from scipy.stats import lognorm
import pandas as pd
import time
import os
from functools import partial 
import math
from joblib import Parallel, delayed ,cpu_count
from joblib.externals.loky import get_reusable_executor
import gc
import WPB315  as WPB ### This is the WPB module from Fortran
print(WPB.mpbfit.__doc__)

##############################
### Setting up variables for the model run 
WPB_LUT=pd.read_csv('WPB_Inputs/WPB_Trees_6_20.csv')
#### Subset to the zones 
WPB_LUT=WPB_LUT[WPB_LUT['Location'].isin(['UpperMed','UpperLow','LowerLow','LowerMed'])]
###update Lookup numbers
WPB_LUT.LUN=list(range(1,21))
#print(WPB_LUT)
locations=WPB_LUT.Location.unique()
### Temperature Max
Tmax_csv=pd.read_csv('WPB_Inputs/MI_Tmax_6_20.csv')
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
### Ensure only the used dates are used 
mask = (Drought_csv['Dates'] > "2001-10-01") & (Drought_csv['Dates'] <= "2018-12-31")
Drought_csv=Drought_csv.loc[mask]
CWD=Drought_csv
### Set the dates for the run 
Dates=pd.DataFrame(Tmax_csv.Dates)

# This method derived from
# https://github.com/Joseph94m/MCMC/blob/master/MCMC.ipynb

##### Function set 
# These functions nest to run the MCMC model 


def p2f(x1,dtype):### This a formatting tool to get the Fortran to accept the imput.  
  x1=np.array(x1,dtype=dtype,order='F')
  return x1.copy()
###  The transition_model proposes new functions for the next step when running the model to create the 
### covariance structure needed for the covariance version of the model. 
transition_model = lambda x2: [np.random.normal(x2[0],0.0010,(1,))[0].round(decimals=5), ### SD of distribution 
                                               ##### Draw one Shape 1d take the first (Lots of things for just one )
                              np.random.normal(x2[1],0.008,(1,))[0].round(decimals=5),
                              np.random.normal(x2[2],0.0007,(1,))[0].round(decimals=5),
                              np.random.normal(x2[3],0.0010,(1,))[0].round(decimals=5),
                              np.random.normal(x2[4],25,(1,))[0].round(decimals=0),
                              x2[5],
                              np.random.normal(x2[6],.00010,(1,))[0].round(decimals=5)] ##

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
    return(Temp_LUT.iloc[:,2],Temp_LUT.iloc[:,0],Temp_LUT.iloc[:,3],ntgeq317,ntgeq00)

  
def Insect_Modeling(Tmin_csv,Tmax_csv,CWD,WPB_LUT,LUN,Xbind,startdate,enddate):
  ## This function parses out the variables from Xbind and parses out which site the 
  ## is supposed to be simulated. Then it fetches the climate variable and the number of trees (Temperature parser).
  ## It lastly runs the fortran model with the appropriate variables and stores the trees for evaluation at the level above. 
    print(LUN)
    r1=Xbind[0]
    x0=Xbind[1] 
    x1=Xbind[2]
    x2=Xbind[3]
    tparents=Xbind[4]
    size_factor=Xbind[6]
    LUN1=LUN
    if(LUN>3):
        r1=Xbind[7]
        x0=Xbind[8] 
        x1=Xbind[9]
        x2=Xbind[10]
        tparents=Xbind[11]
        size_factor=Xbind[13]
        LUN1=LUN-4       
        if(LUN>7 ):
            r1=Xbind[14]
            x0=Xbind[15] 
            x1=Xbind[16]
            x2=Xbind[17]
            tparents=Xbind[18]
            size_factor=Xbind[20]
            LUN1=LUN-8
            if(LUN>11 ):
                r1=Xbind[21]
                x0=Xbind[22] 
                x1=Xbind[23]
                x2=Xbind[24]
                tparents=Xbind[25]
                size_factor=Xbind[27]
                LUN1=LUN-12
                if(LUN> 15 ):
                    r1=Xbind[28]
                    x0=Xbind[29] 
                    x1=Xbind[30]
                    x2=Xbind[31]
                    tparents=Xbind[32]
                    size_factor=Xbind[34]
                    LUN1=LUN-16
                    if(LUN>19):
                        r1=Xbind[35]
                        x0=Xbind[36] 
                        x1=Xbind[37]
                        x2=Xbind[38]
                        tparents=Xbind[39]
                        size_factor=Xbind[41]
                        LUN1=LUN-20 
                        if(LUN>23):
                            r1=Xbind[42]
                            x0=Xbind[43] 
                            x1=Xbind[44]
                            x2=Xbind[45]
                            tparents=Xbind[46]
                            size_factor=Xbind[48]
                            LUN1=LUN-24 
                            if(LUN>27):
                                r1=Xbind[49]    
                                x0=Xbind[50] 
                                x1=Xbind[51]
                                x2=Xbind[52]
                                tparents=Xbind[53]
                                size_factor=Xbind[55]
                                LUN1=LUN-28
        
        
        
    print(LUN1)
    ### Pull the correct Temperatures 
    SubMin,SubMax,CWDvec,ntgeq317,ntgeq00 =TemperatureParser(LUN1,
                                                             Tmin_csv,
                                                             Tmax_csv,
                                                             CWD,WPB_LUT,
                                                             startdate,
                                                             enddate)

    #### Run the model 
    NGEQ20vec,NGEQ00vec,Flight,Fec,Eggs,L1,L2,Pout,Ten,Adults= WPB.mpbfit(p2f(SubMin,"cdouble"),#TempMin
                p2f(SubMax,"cdouble"),#TempMax
                p2f(r1,"float"), #phi1
                p2f(x0,"float"),
                p2f(x1,"float"),
                p2f(CWDvec,"cdouble"),
                p2f(x2,"float"),
                p2f(size_factor,"float"),
                p2f(ntgeq317,"float"),
                p2f(ntgeq00,"float"),
                p2f(tparents,"float"))
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
    gc.collect()
    #### Return the smallest number of living trees. 
    return(LUN,BigTreesOut,SmallTreesOut)
   


#Computes the likelihood of the data given a sigma (new or current) according to equation (2)
### In this version this is where all the parallelized work is happening 

def manual_log_like_normal(WPB_LUT,Tmax_csv,Tmin_csv,CWD,xbind):  
  ### This function calculates the log likelihood for each of the proposed steps for each of the four testing runs. 
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
  ## This is from the the package joblib
  #### To speed this whole deal up, the model is structured as as 24 independent runs
  # That is 4 model sets ( Each with a combination of 3 of 4 sites) repeated twice (current position, proposed step). 
  
  Run1,Run2,Run3,Run1a,Run2a,Run3a,Run1b,Run2b,Run3b,Run1c,Run2c,Run3c,Run1d,Run2d,Run3d,Run1e,Run2e,Run3e,Run1f,Run2f,Run3f,Run1g,Run2g,Run3g=Parallel(n_jobs=32)(delayed(PartIM)(LUN=i) for i in  list([0,1,2,
                                                                                4,5,6,
                                                                                8,9,11,
                                                                                12,13,15,
                                                                                16,18,19,
                                                                                20,22,23,
                                                                                25,26,27,
                                                                                29,30,31]))
  ### Clear cached memory to ensure the parallel doesn't bleed between runs. 
  get_reusable_executor().shutdown(wait=True)

## Now we parse out all those runs to evaluate the Log-liklihood comparison for each run 
   
   
####### Model set 1 
  ## Each of these sets organizes the model runs and the observed data, 
  ## then calculates the lognomral log liklihood distribution 
  ## It returns one score for the current (eg,LL1) and one for the proposed (eg LLnew1)
  Outputmodel1=np.concatenate([Run1[1],Run2[1],Run3[1],Run1[2],Run2[2],Run3[2]])
  Outputmodelnew1=np.concatenate([Run1a[1],Run2a[1],Run3a[1],Run1a[2],Run2a[2],Run3a[2]])
        
  OutputY_obs1=np.concatenate([np.array(WPB_LUT.iloc[0,4:17]),np.array(WPB_LUT.iloc[1,4:17]),
                            np.array(WPB_LUT.iloc[2,4:17]),
                            np.array(WPB_LUT.iloc[0,17:30]),np.array(WPB_LUT.iloc[1,17:30]),
                            np.array(WPB_LUT.iloc[2,17:30])])
      
  OutputY_obs1=OutputY_obs1.astype(np.float64)  
  Outputmodel1=Outputmodel1.astype(np.float64) 
  Outputmodelnew1=Outputmodelnew1.astype(np.float64) 
  Outputmodel1[Outputmodel1<.1]=.1
  Outputmodelnew1[Outputmodelnew1<.1]=.1
  
  sigma=np.array(OutputY_obs1*.0983).astype(np.float64) 
  mean=np.log((OutputY_obs1**2)/(np.sqrt((OutputY_obs1**2)+(sigma**2))) )
  sigmaprime=np.sqrt(np.log(1+((sigma**2)/(OutputY_obs1**2))))

  LL1 = sum(-np.log((np.exp(-(np.log(Outputmodel1) - mean)**2 / (2 * sigmaprime**2))
       / ((Outputmodel1*sigmaprime * np.sqrt(2 * np.pi))))))
  LLnew1 = sum(-np.log((np.exp(-(np.log(Outputmodelnew1) - mean)**2 / (2 * sigmaprime**2))
       / ((Outputmodelnew1*sigmaprime * np.sqrt(2 * np.pi))))))
  ##### Assess Model set2 
  Outputmodel2=np.concatenate([Run1b[1],Run2b[1],Run3b[1],Run1b[2],Run2b[2],Run3b[2]])
  Outputmodelnew2=np.concatenate([Run1c[1],Run2c[1],Run3c[1],Run1c[2],Run2c[2],Run3c[2]])
  

  OutputY_obs2=np.concatenate([np.array(WPB_LUT.iloc[0,4:17]),np.array(WPB_LUT.iloc[1,4:17]),
                            np.array(WPB_LUT.iloc[3,4:17]),
                            np.array(WPB_LUT.iloc[0,17:30]),np.array(WPB_LUT.iloc[1,17:30]),
                            np.array(WPB_LUT.iloc[3,17:30])])   
  OutputY_obs2=OutputY_obs2.astype(np.float64)  
  Outputmodel2=Outputmodel2.astype(np.float64) 
  Outputmodelnew2=Outputmodelnew2.astype(np.float64) 
  Outputmodel2[Outputmodel2<.1]=.1
  Outputmodelnew2[Outputmodelnew2<.1]=.1
  
  sigma=np.array(OutputY_obs2*.0983).astype(np.float64) 
  mean=np.log((OutputY_obs2**2)/(np.sqrt((OutputY_obs2**2)+(sigma**2))) )
  sigmaprime=np.sqrt(np.log(1+((sigma**2)/(OutputY_obs2**2))))
  LL2 = sum(-np.log((np.exp(-(np.log(Outputmodel2) - mean)**2 / (2 * sigmaprime**2))
       / ((Outputmodel2*sigmaprime * np.sqrt(2 * np.pi))))))
  LLnew2 = sum(-np.log(
          (np.exp(-(np.log(Outputmodelnew2) - mean)**2 / (2 * sigmaprime**2))
          
       / ((Outputmodelnew2*sigmaprime * np.sqrt(2 * np.pi))))))
 ###  Assess Model 3 
  Outputmodel3=np.concatenate([Run1d[1],Run2d[1],Run3d[1],Run1d[2],Run2d[2],Run3d[2]])
  Outputmodelnew3=np.concatenate([Run1e[1],Run2e[1],Run3e[1],Run1e[2],Run2e[2],Run3e[2]])
        
  OutputY_obs3=np.concatenate([np.array(WPB_LUT.iloc[0,4:17]),np.array(WPB_LUT.iloc[2,4:17]),
                            np.array(WPB_LUT.iloc[3,4:17]),
                            np.array(WPB_LUT.iloc[0,17:30]),np.array(WPB_LUT.iloc[2,17:30]),
                            np.array(WPB_LUT.iloc[3,17:30])])
      
  OutputY_obs3=OutputY_obs3.astype(np.float64)  
  Outputmodel3=Outputmodel3.astype(np.float64) 
  Outputmodelnew3=Outputmodelnew3.astype(np.float64) 
  Outputmodel3[Outputmodel3<.1]=.1
  Outputmodelnew3[Outputmodelnew3<.1]=.1
  
  sigma=np.array(OutputY_obs3*.0983).astype(np.float64) 
  mean=np.log((OutputY_obs3**2)/(np.sqrt((OutputY_obs3**2)+(sigma**2))) )
  sigmaprime=np.sqrt(np.log(1+((sigma**2)/(OutputY_obs3**2))))

  LL3 = sum(-np.log((np.exp(-(np.log(Outputmodel3) - mean)**2 / (2 * sigmaprime**2))
       / ((Outputmodel3*sigmaprime * np.sqrt(2 * np.pi))))))
  LLnew3 = sum(-np.log((np.exp(-(np.log(Outputmodelnew3) - mean)**2 / (2 * sigmaprime**2))
       / ((Outputmodelnew3*sigmaprime * np.sqrt(2 * np.pi)))))) 
  ###  Assess Model 4 
  Outputmodel4=np.concatenate([Run1f[1],Run2f[1],Run3f[1],Run1f[2],Run2f[2],Run3f[2]])
  Outputmodelnew4=np.concatenate([Run1g[1],Run2g[1],Run3g[1],Run1g[2],Run2g[2],Run3g[2]])
        
  OutputY_obs4=np.concatenate([np.array(WPB_LUT.iloc[1,4:17]),np.array(WPB_LUT.iloc[2,4:17]),
                            np.array(WPB_LUT.iloc[3,4:17]),
                            np.array(WPB_LUT.iloc[1,17:30]),np.array(WPB_LUT.iloc[2,17:30]),
                            np.array(WPB_LUT.iloc[3,17:30])])
      
  OutputY_obs4=OutputY_obs4.astype(np.float64)  
  Outputmodel4=Outputmodel4.astype(np.float64) 
  Outputmodelnew4=Outputmodelnew4.astype(np.float64) 
  Outputmodel4[Outputmodel4<.1]=.1
  Outputmodelnew4[Outputmodelnew4<.1]=.1
  
  sigma=np.array(OutputY_obs4*.0983).astype(np.float64) 
  mean=np.log((OutputY_obs4**2)/(np.sqrt((OutputY_obs4**2)+(sigma**2))) )
  sigmaprime=np.sqrt(np.log(1+((sigma**2)/(OutputY_obs4**2))))

  LL4 = sum(-np.log((np.exp(-(np.log(Outputmodel4) - mean)**2 / (2 * sigmaprime**2))
       / ((Outputmodel4*sigmaprime * np.sqrt(2 * np.pi))))))
  LLnew4 = sum(-np.log((np.exp(-(np.log(Outputmodelnew4) - mean)**2 / (2 * sigmaprime**2))
       / ((Outputmodelnew4*sigmaprime * np.sqrt(2 * np.pi)))))) 
  
  
  
  

  # printers for checking progress
  #print("Run1")
  #print("Old x")
  #print(LL1)
  #print("New X")
  #print(LLnew1)
  #print("Run2")
  #print("Old x")
  #print(LL2)
  #print("New X")
  #print(LLnew2)
  #print("Run3")
  #print("Old x")
  #print(LL3)
  #print("New X")
  #print(LLnew3)
  #print("Run4")
  #print("Old x")
  #print(LL4)
  #print("New X")
  #print(LLnew4)
  ### Return the LL score for each of four pairs of runs 
  return(LL1,LLnew1,LL2,LLnew2,LL3,LLnew3,LL4,LLnew4)


# Defines whether to accept or reject the new sample
def acceptance(x_likelihood, x_new_likelihood):
	### A metropolis-hastings rule acceptance 
    if x_new_likelihood<x_likelihood:
        return True
    else:
        accept=np.random.uniform(0,1)
        # Since we did a log likelihood, we need to exponentiate in order to compare to the random number
        # less likely x_new are less likely to be accepted
        return (accept < np.exp(x_likelihood-x_new_likelihood))

def metropolis_hastings(likelihood_computer, cov,iterations,WPB_LUT,Tmax_csv,Tmin_csv,CWD,acceptance_rule):

	### This function manages the testing of each model using a metropolis hastings step 
	# A) it proposes a new step for each model, and insures it is not outside the realm of the posssible (Bounce())
	# B) It binds all the current and proposed moves into a single vector which is passed to the manual_log_like_normal function 
	# resulting in 2 scores each for 4 models. 
	# C) It then runs the acceptance algorithm for each set of model. 
	# D) It returns vectors containing the all moves, accepted moves, rejected moves, and the logliklihood scores associated with each proposed move
	### A proposed starting place for each model 
	x1=[0.201323116,10.02008261,1.544441227,0.502509843,851.5005162,2,0.049387451]
    x2=[0.200345995,10.00041733,1.545792262,0.506219717,859.4871259,2,0.049914916]
    x3=[0.187506542,10.1011393,61.549765959,0.613737724,753.2828798,2,0.058028045]
    x4=[0.208986903,10.19907002,1.542914988,0.506099424,739.4243006,2,0.045617422]

    accepted1 = []
    rejected1 = []
    LLlist1=[]
    alllist1=[]
    accepted2 = []
    rejected2 = []
    LLlist2=[]
    alllist2=[]
    accepted3 = []
    rejected3 = []
    LLlist3=[]
    alllist3=[]
    accepted4 = []
    rejected4 = []
    LLlist4=[]
    alllist4=[]
    
    
    start = time.time()
    for i in range(iterations):
        if((i/50) % 1 ==0):
            print("#################################"+str(i)+"###########################################")
        print(i)
        #### A)
        x_new1 = np.random.multivariate_normal(x1,cov, size=1) .tolist()[0]
        x_new1=Bounce(x_new1)
        x_new2 =  np.random.multivariate_normal(x2,cov, size=1) .tolist()[0]
        x_new2=Bounce(x_new2)
        x_new3 =  np.random.multivariate_normal(x3,cov, size=1) .tolist()[0]
        x_new3=Bounce(x_new3)
        x_new4 = np.random.multivariate_normal(x4,cov, size=1) .tolist()[0]
        x_new4 =Bounce(x_new4)
       
	    ##### B) 
        Xbind=np.append(np.array(x1),np.array(x_new1))
        Xbind=np.append(Xbind,np.array(x2))
        Xbind=np.append(Xbind,np.array(x_new2))
        Xbind=np.append(Xbind,np.array(x3))
        Xbind=np.append(Xbind,np.array(x_new3))
        Xbind=np.append(Xbind,np.array(x4))
        Xbind=np.append(Xbind,np.array(x_new4))
        x_lik1,x_new_lik1,x_lik2,x_new_lik2,x_lik3,x_new_lik3,x_lik4,x_new_lik4= manual_log_like_normal(WPB_LUT,Tmax_csv,Tmin_csv,CWD,Xbind)
         
		##### C)
        #### Switching up the X's 
        print("FirstModel")
        if (acceptance_rule(x_lik1,x_new_lik1)):            
            x1 = x_new1
            accepted1.append(x_new1)
            LLlist1.append(x_new_lik1)
            alllist1.append(x_new1)
            print("accepted")
        else:
            rejected1.append(x_new1) 
            LLlist1.append(x_lik1)
            alllist1.append(x_new1)
            print("rejected")
        print("Second Model")
        if (acceptance_rule(x_lik2,x_new_lik2)):            
            x2 = x_new2
            accepted2.append(x_new2)
            LLlist2.append(x_new_lik2)
            alllist2.append(x_new2)
            print("accepted")
        else:
            rejected2.append(x_new2) 
            LLlist2.append(x_lik2)
            alllist2.append(x_new2)
            print("rejected")
        print("Third Model")
        if (acceptance_rule(x_lik3,x_new_lik3)):            
            x3 = x_new3
            accepted3.append(x_new3)
            LLlist3.append(x_new_lik3)
            alllist3.append(x_new3)
            print("accepted")
        else:
            rejected3.append(x_new3) 
            LLlist3.append(x_lik3)
            alllist3.append(x_new3)
            print("rejected")    
        print("Fourth Model")
        if (acceptance_rule(x_lik4,x_new_lik4)):            
            x4 = x_new4
            accepted4.append(x_new4)
            LLlist4.append(x_new_lik4)
            alllist4.append(x_new4)
            print("accepted")
        else:
            rejected4.append(x_new4) 
            LLlist4.append(x_lik4)
            alllist4.append(x_new4)
            print("rejected")        
    end = time.time()
    print("Run "+str(end - start)+" seconds") 
	# D)  	
    return np.array(accepted1), np.array(rejected1),np.array(LLlist1),np.array(alllist1),np.array(accepted2), np.array(rejected2),np.array(LLlist2),np.array(alllist2),np.array(accepted3), np.array(rejected3),np.array(LLlist3),np.array(alllist3),np.array(accepted4), np.array(rejected4),np.array(LLlist4),np.array(alllist4)

### Running portion 
### For this version the proposed move is calculated by the multi-normal distribution with a covariance structure taken
### from runs in which a normal distribution was proposed with the parameter stated above. This vector of accepted moves 
### is then used to create a covariance matrix by which to propose new moves here. 
x=[]
CovSet=pd.read_csv("Accept_Practice_TestCOV2.csv")
cov=np.cov(CovSet.transpose())
### The model is then run for the number of iterations needed. 
accepted1, rejected1,LLoutput1,alllist1,accepted2, rejected2,LLoutput2,alllist2,accepted3, rejected3,LLoutput3,alllist3,accepted4, rejected4,LLoutput4,alllist4= metropolis_hastings(manual_log_like_normal,### Computes Liklihood
                                         cov, 
                                        8500, ###Iterations
                                         WPB_LUT, ### LUT information including LUN for parproc
                                         Tmax_csv,### All Tmax
                                         Tmin_csv,### All Tmin 
                                         CWD,
                                         acceptance) ##Acceptance rule (Metropolis)
### Saving each of many outputs
# Accepted - Accepted moves
# Rejected - Rejected moves
# All - Every proposed moved
# LL - the log-liklihood score of each step. 

Run="Run_11_22"
np.savetxt("ManyTests/AcceptTesting1_"+Run+".csv",accepted1,delimiter=",")
np.savetxt("ManyTests/Reject_Testing1_"+Run+".csv",rejected1,delimiter=",")
np.savetxt("ManyTests/AllPOWTesting1_"+Run+".csv",alllist1,delimiter=",")
np.savetxt("ManyTests/LL_POWTesting1_"+Run+".csv",LLoutput1,delimiter=",")
np.savetxt("ManyTests/AcceptTesting2_"+Run+".csv",accepted2,delimiter=",")
np.savetxt("ManyTests/Reject_Testing2_"+Run+".csv",rejected2,delimiter=",")
np.savetxt("ManyTests/AllPOWTesting2_"+Run+".csv",alllist2,delimiter=",")
np.savetxt("ManyTests/LL_POWTesting2_"+Run+".csv",LLoutput2,delimiter=",")
np.savetxt("ManyTests/AcceptTesting3_"+Run+".csv",accepted3,delimiter=",")
np.savetxt("ManyTests/Reject_Testing3_"+Run+".csv",rejected3,delimiter=",")
np.savetxt("ManyTests/AllPOWTesting3_"+Run+".csv",alllist3,delimiter=",")
np.savetxt("ManyTests/LL_POWTesting3_"+Run+".csv",LLoutput3,delimiter=",")
np.savetxt("ManyTests/AcceptTesting4_"+Run+".csv",accepted4,delimiter=",")
np.savetxt("ManyTests/Reject_Testing4_"+Run+".csv",rejected4,delimiter=",")
np.savetxt("ManyTests/AllPOWTesting4_"+Run+".csv",alllist4,delimiter=",")
np.savetxt("ManyTests/LL_POWTesting4_"+Run+".csv",LLoutput4,delimiter=",")



