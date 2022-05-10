!####################LICENSING####################################################### 
!All code here in is distributed under an OSS License © 2022. 
!Triad National Security, LLC. All rights reserved. This program 
!was produced under U.S. Government contract 89233218CNA000001 
!for Los Alamos National Laboratory (LANL), which is operated by 
!Triad National Security, LLC for the U.S. Department of 
!Energy/National Nuclear Security Administration. All rights 
!in the program are reserved by Triad National Security, LLC, 
!and the U.S. Department of Energy/National Nuclear Security 
!Administration. The Government is granted for itself and others 
!acting on its behalf a nonexclusive, paid-up, irrevocable worldwide 
!license in this material to reproduce, prepare derivative works, 
!distribute copies to the public, perform publicly and display 
!publicly, and to permit others to do so.

!And a BSD license: This program is open source under the BSD-3 
!License. Redistribution and use in source and binary forms, with 
!or without modification, are permitted provided that the following 
!conditions are met:

! 1.Redistributions of source code must retain the above copyright 
!notice, this list of conditions and the following disclaimer. 
!
! 2.Redistributions in binary form must reproduce the above 
!copyright notice, this list of conditions and the following 
!disclaimer in the documentation and/or other materials provided 
! with the distribution. 
! 3.Neither the name of the copyright holder 
!nor the names of its contributors may be used to endorse or promote 
!products derived from this software without specific prior written 
!permission.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
!FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
!COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
!INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
!BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
!LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
!ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
!POSSIBILITY OF SUCH DAMAGE. 
! ########################End#########################################################

!Part of the IMAP-TDIA Code Base: Zachary Robbins and Chonggang Xu




! This is the fortran implementation of the IMAP and TDIA model for the western pine beetle
! as was implemented in the paper "Warming increased bark beetle mortality
! by thirty percent during an extreme drought in California" 
! This was authored by Zachary Robbins and Chonggang Xu.  
! This is modified from the original IMAP code authored by Devin W. Goodsman, while he was a postdoctoral researcher at
! LANL (2017-2018) except for the subroutines for fast Fourier transformations, which were
! written by Paul Swarztrauber, and the subroutines for the incomplete beta function, which
! I modified from the reference below.

!********************************************************
!* Reference:                                           *
!* "Numerical Recipes, By W.H. Press, B.P. Flannery,    *
!*  S.A. Teukolsky and T. Vetterling, Cambridge         *
!*  University Press, 1986" [BIBLI 08].                 *
!*                                                      *
!* I have refactored the incopmplete beta function code *
!* such that the functions become subroutines           *
!* (Goodsman, D.W. Oct. 2, 2018). I converted these to  *
!* subroutines to more easily call them from R and to   *
!* more consistently call them within other Fortran     *
!* programs that are subroutine-based.                  *
!********************************************************





! To host fortran code as a function set with python 3.x 
! Ensure the numpy is installed, as is mingw64 (or similar fortran compiler_)
! This might require some path setting see https://stackoverflow.com/questions/48826283/compile-fortran-module-with-f2py-and-python-3-6-on-windows-10
! Use the following function in the python console C:\Users\username >f2py -c WPB_submission.F90 -m WPB_Module
! Load into python using traditional module loading format 
! eg: 
! import WPB_Module as WPB 


! This is the main do loop for the length of the model run (Temperature vector -1). It handles the model set up in preparation to run
! BeetleSim which does the work of growing beetles and attacking trees. 

subroutine BeetleFit(Length,Tminvec, Tmaxvec,r1,x0,x1,CWDvec,x2,SizeFactor,NtGEQ317,NtGEQ00,Tparents,NGEQ20vec,NGEQ00vec,&
	FebInPopnvec,Fecoutvec, Eoutvec, L1outvec, L2outvec, Poutvec, Teoutvec, Aoutvec)

    implicit none

    Interface

        subroutine BeetleSim(Tmax, Tmin,TminNext, Parents, FA, OE, OL1, OL2, &
			NewParentstm1, OPare, Pare, ActiveParents, &
            OL3, OL4, OP, OT, NewEggstm1, NewL1tm1, &
            NewL2tm1, NewL3tm1, NewL4tm1, NewPtm1, NewTtm1, &
            Fec, E, L1, L2, L3, L4, P, Te, A, ColdestT, &
            NtGEQ317,NtGEQ00, Bt,r1,x0,x1,CWD,x2,SizeFactor,FebInPopn,&
			Fecout, Eout, L1out, L2out, Pout, Teout, Aout,Flyout,EndMPBPopn,BigTreesOut,SmallTreesOut) 
            real(kind = 8), intent(in) :: Tmax
            real(kind = 8), intent(in) :: Tmin
			real(kind = 8), intent(in) :: TminNext
            real(kind = 8), intent(inout) :: Parents
            real(kind = 8), intent(out) :: FA
            real(kind = 8), intent(in) :: CWD  ! The Current Climatic water deficit
            complex(kind = 8), intent(inout) :: OE(2**8)
            complex(kind = 8), intent(inout) :: OL1(2**8)
            complex(kind = 8), intent(inout) :: OL2(2**8)
            complex(kind = 8), intent(inout) :: OL3(2**8)
            complex(kind = 8), intent(inout) :: OL4(2**8)
            complex(kind = 8), intent(inout) :: OP(2**8)
            complex(kind = 8), intent(inout) :: OT(2**8)
			complex(kind = 8), intent(inout) :: OPare(2**8)
			real(kind = 8), intent(inout) :: NewParentstm1
            real(kind = 8), intent(inout) :: NewEggstm1
            real(kind = 8), intent(inout) :: NewL1tm1
            real(kind = 8), intent(inout) :: NewL2tm1
            real(kind = 8), intent(inout) :: NewL3tm1
            real(kind = 8), intent(inout) :: NewL4tm1
            real(kind = 8), intent(inout) :: NewPtm1
            real(kind = 8), intent(inout) :: NewTtm1
            real(kind = 8), intent(inout) :: NtGEQ317
            real(kind = 8), intent(inout) :: NtGEQ00
            real(kind = 8), intent(inout) :: Fec
			real(kind = 8), intent(inout) :: ActiveParents
			real(kind = 8), intent(inout) :: Pare
            real(kind = 8), intent(inout) :: E
            real(kind = 8), intent(inout) :: L1
            real(kind = 8), intent(inout) :: L2
            real(kind = 8), intent(inout) :: L3
            real(kind = 8), intent(inout) :: L4
            real(kind = 8), intent(inout) :: P 
            real(kind = 8), intent(inout) :: Te
            real(kind = 8), intent(inout) :: A
            real(kind = 8), intent(in) :: ColdestT
			real(kind = 8), intent(out) :: Fecout
            real(kind = 8), intent(out) :: Eout
            real(kind = 8), intent(out) :: L1out
            real(kind = 8), intent(out) :: L2out
            real(kind = 8), intent(out) :: Pout
            real(kind = 8), intent(out) :: Teout
            real(kind = 8), intent(out) :: Aout
			real(kind = 8), intent(out) :: Flyout


            real(kind = 8), intent(inout) :: Bt
            real(kind = 8), intent(in) :: r1                    ! controls beetle clustering
            real(kind = 8), intent(in) :: x0   ! The intercept on attack 
            real(kind = 8), intent(in) :: x1   ! x1 controls the impact CWD
            real(kind = 8), intent(in) :: x2   ! x2 controls the impact of size class    
			real(kind = 8), intent(in) :: SizeFactor 
            real(kind = 8), intent(in) :: FebInPopn
            real(kind = 8), intent(in) :: EndMPBPopn
            real(kind = 8), intent(out) :: BigTreesOut
            real(kind = 8), intent(out) :: SmallTreesOut
        end subroutine BeetleSim

    End Interface

    ! Here are the input and output variables and parameters for the main subroutine above.
    ! Each of the arrays is the same length as the daily temperature time series.
    integer(kind = 4), intent(in) ::Length
	real(kind = 8), intent(in) :: Tminvec(Length)
    real(kind = 8), intent(in) :: Tmaxvec(Length)
    real(kind = 8), intent(in) :: CWDvec(Length)
    real(kind = 8), intent(inout) :: NtGEQ317
    real(kind = 8), intent(inout) :: NtGEQ00
    real(kind = 8), intent(inout) :: Tparents  !Initial Conditions parents 
    real(kind = 8), intent(inout) :: r1   ! controls beetle clustering
    real(kind = 8), intent(in) :: x0   ! The intercept on attack
    real(kind = 8), intent(in) :: x1   ! x1 controls the impact CWD
    real(kind = 8), intent(in) :: x2   ! x2 controls the impact of size class    
    real(kind = 8), intent(in) :: SizeFactor   ! x2 controls the impact of size class    
    !real(kind = 8), intent(inout) ::NtGEQ20
    real(kind = 8), intent(out) :: NGEQ20vec(Length)
    real(kind =8),  intent(out) :: NGEQ00vec(Length)
    real(kind = 8), intent(out) :: FebInPopnvec(Length)
	real(kind = 8), intent(out) :: Fecoutvec (Length)
    real(kind =8),  intent(out) :: Eoutvec (Length)
    real(kind = 8), intent(out) :: L1outvec(Length)
    real(kind = 8), intent(out) :: L2outvec (Length)
    real(kind =8),  intent(out) :: Poutvec (Length)
    real(kind = 8), intent(out) :: Teoutvec(Length)			
	real(kind = 8), intent(out) :: Aoutvec (Length)	
    real(kind = 8) :: EndMPBPopn
    real(kind = 8) :: FebInPopn

    ! Here are the internal parameters and variables for the MPBSim subroutine

    real(kind = 8) :: Tmax
    real(kind = 8) :: Tmin
	real(kind = 8) :: TminNext
    real(kind = 8) :: Parents
    real(kind = 8) :: CWD  ! The Current Climatic water deficit
    real(kind = 8) :: FA

    ! Containers for the distributions of physiological age for each life stage
    complex(kind = 8) :: OE(2**8)
    complex(kind = 8) :: OL1(2**8)
    complex(kind = 8) :: OL2(2**8)
    complex(kind = 8) :: OL3(2**8)
    complex(kind = 8) :: OL4(2**8)
    complex(kind = 8) :: OP(2**8)
    complex(kind = 8) :: OT(2**8)
    complex(kind = 8) :: OPare(2**8)
    real(kind = 8) :: BigTreesOut
    real(kind = 8) :: SmallTreesOut
    real(kind = 8) :: NewParentstm1
	real(kind = 8) :: NewEggstm1
    real(kind = 8) :: NewL1tm1
    real(kind = 8) :: NewL2tm1
    real(kind = 8) :: NewL3tm1
    real(kind = 8) :: NewL4tm1
    real(kind = 8) :: NewPtm1
    real(kind = 8) :: NewTtm1
	real(kind = 8) :: ActiveParents
	real(kind = 8) :: Pare 
    real(kind = 8) :: Fec               ! the expected number of pre-eggs at each time
    real(kind = 8) :: E                 ! the expected number of eggs at each time
    real(kind = 8) :: L1                ! the expected number of L1 at each time step
    real(kind = 8) :: L2                ! the expected number of L2 at each time step
    real(kind = 8) :: L3                ! the expected number of L3 at each time step
    real(kind = 8) :: L4                ! the expected number of L4 at each time step
    real(kind = 8) :: P                 ! the expected number of pupae at each time step
    real(kind = 8) :: Te                ! the expected number of tenerals at each time step
    real(kind = 8) :: A                 ! the expected number of flying adults at each time step
	real(kind = 8) :: Fecout            ! the expected number of pre-eggs at each time
	real(kind = 8) :: Eout
    real(kind = 8) :: L1out	
    real(kind = 8):: L2out		
	real(kind = 8) :: Pout              
	real(kind = 8):: Teout
	real(kind = 8):: Aout
	real(kind = 8):: Flyout
	integer(kind = 4):: TrueLength
    ! Winter mortality is a function of the coldest temperature to date.
    real(kind = 8) :: ColdestT

    real(kind = 8) :: Bt                ! cumulative attacking beetles

    integer(kind = 4) :: i              ! iterator

 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! All of the life stages must be initialized
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Opare =0.0
	NewParentstm1 = 0.0
    ! Initializing the eggs
    OE = 0.0
    NewEggstm1 = 0.0

    ! Initializing the L1 larvae
    OL1 = 0.0
    NewL1tm1 = 0.0

    ! Initializing the L2 larvae
    OL2 = 0.0
    NewL2tm1 = 0.0

    ! Initializing the L3 larvae
    OL3 = 0.0
    NewL3tm1 = 0.0

    ! Initializing the L4 larvae
    OL4 = 0.0
    NewL4tm1 = 0.0

    ! Initializing the pupae
    OP = 0.0
    NewPtm1 = 0.0

    ! Initializing the teneral adults
    OT = 0.0
    NewTtm1 = 0.0

    ! Initializing the first element of each container
    Fec = 0.0
    E = 0.0
    L1 = 0.0
    L2 = 0.0
    L3 = 0.0
    L4 = 0.0
    P = 0.0
    Te = 0.0
    A = 0.0
    FA = 0.0
    Bt = 0.0
	
	! The number of initial parents equals the endemic reset. 
    EndMPBPopn = Tparents
    NGEQ20vec(1) = NtGEQ317
    NGEQ00vec(1) = NtGEQ00
    FebInPopn = 0.0

    ColdestT = 15.0
	! The true length of the simulation is reduced by one as the each days temperature
	! is calculated using that day and the next. 
    TrueLength=Length-1
    
    do i = 2,TrueLength
        Tmax = Tmaxvec(i)
        Tmin = Tminvec(i)
		TminNext=Tminvec(i+1)
        CWD = CWDvec(i)
        if(Tmin < ColdestT)then
            ColdestT = Tmin
        end if

        !Computing the insect population on the initialization of the current year (October 15th).
		! This resets the base population in the event of collapse
        if(i== 1 .or. i ==  1+365 .or. i == 1+2*365 .or. i ==  1+3*365 .or. &
            i ==  1+4*365 .or. i ==  1+5*365 .or. i == 1+6*365 .or. &
			i ==  1+7*365 .or. i ==  1+8*365 .or. i ==  1+9*365 .or. &
			i ==  1+10*365 .or. i ==  1+11*365 .or. i ==  1+12*365 .or.  &
			i ==  1+13*365 .or. i ==  1+14*365 .or. i ==  1+15*365 .or. i ==  1+16*365 .or. &
			i ==  1+17*365) then
        FebInPopn = Fec + E + (L1 + L2 )/(1.0 + exp(-(ColdestT + 29.22791)/2.0)) + P + Te + A
        
        end if

        ! If the total insect population on October 15th (ignoring leap years) of the year is less than
        ! the endemic population, we re-establish an endemic population level in October of the current year.
        if(i == 2 + 365 .or. i == 2 + 2*365 .or. i == 2 + 3*365 .or. &
            i == 2 + 4*365 .or. i == 2 + 5*365 .or. i == 2 + 6*365 .or. &
			i == 2 + 7*365 .or. i == 2 + 8*365 .or. i == 2 + 9*365 .or. &
			i == 2 + 10*365 .or. i == 2 + 11*365 .or. i == 2 + 12*365 .or.  &
			i ==  2+13*365 .or. i ==  2+14*365 .or. i ==  2+15*365 .or. i ==  2+16*365 .or. &
			i ==  2+17*365) then
			if(FebInPopn < EndMPBPopn)then
				Parents =EndMPBPopn
			end if
        end if

        ! On April 1st we reset the coldest temperature to date to 15.0 degrees Celsius
        if(i == 184  + 365 .or. i ==  184 + 2*365 .or. i ==  184 + 3*365 .or. &
            i == 184  + 4*365 .or. i ==  184 + 5*365 .or. i ==  184 + 6*365 .or. &
			i ==  184+ 7*365 .or. i ==  184+ 8*365 .or. i ==  184 + 9*365 .or. &
			i ==  184 + 10*365 .or. i ==  184 + 11*365 .or. i ==  184 + 12*365 .or.  &
			i ==   184+13*365 .or. i ==   184+14*365 .or. i ==   184+15*365 .or. i ==  184 +16*365 .or. &
			i ==   184+17*365) then
            ColdestT = 15.0

        end if
		! On December 1  we reset the Adults and Teneral Adults to make sure there 
		! Populations don't live longer than expected. 
        if(i == 60 + 365 .or. i == 60  + 2*365 .or. i == 60  + 3*365 .or. &
            i == 60 + 4*365 .or. i == 60  + 5*365 .or. i == 60+ 6*365 .or. &
			i == 60 + 7*365 .or. i == 60  + 8*365 .or. i == 60+ 9*365 .or. &
			i == 60 + 10*365 .or. i == 60  + 11*365 .or. i == 60+ 12*365 .or.  &
			i ==  60+13*365 .or. i ==  60 +14*365 .or. i ==  60+15*365 .or. i ==  60+16*365 .or. &
			i ==  60+17*365) then
			A=0
			TE=0
			Bt = 0.0
			FA = 0.0
        end if
        ! I use a uniform distribution of the total number of parents for seven days 
        ! to initialized the model.
        if(i <= 7)then
            ActiveParents = TParents/7
        end if

        ! Now calling the full MPB simulation for the time step.
        call BeetleSim(Tmax, Tmin,TminNext, Parents, FA, OE, OL1, OL2, &
			NewParentstm1, OPare, Pare, ActiveParents, &
            OL3, OL4, OP, OT, NewEggstm1, NewL1tm1, &
            NewL2tm1, NewL3tm1, NewL4tm1, NewPtm1, NewTtm1, &
            Fec, E, L1, L2, L3, L4, P, Te, A, ColdestT, &
            NtGEQ317,NtGEQ00, Bt,r1,x0,x1,CWD,x2,SizeFactor,FebInPopn,&
			Fecout, Eout, L1out, L2out,Pout, Teout, Aout,Flyout,EndMPBPopn,BigTreesOut,SmallTreesOut)
		!Update the number of trees. 	
		NtGEQ317= BigTreesOut
        NtGEQ00= SmallTreesOut	
        ! This will update E, L1, L2, L3, L4, P, Te, A, and FA
        ! as well as all of the physiological age arrays and
        ! scalars for individuals transitioning from one stage to
        ! another.
        ! These are the diagnostic outputs

        NGEQ20vec(i) = BigTreesOut
        NGEQ00vec(i) = SmallTreesOut
        FebInPopnvec(i)=Flyout
		Fecoutvec(i)= Fecout 
		Eoutvec(i)=Eout 
		L1outvec(i)=L1out
		L2outvec(i)=L2out
		Poutvec(i)=Pout
		Teoutvec(i)=Teout
		Aoutvec(i)=Aout
     end do

End subroutine BeetleFit

!========================================================================================
Subroutine BeetleSim(Tmax, Tmin,TminNext, Parents, FA, OE, OL1, OL2, &
			NewParentstm1, OPare, Pare, ActiveParents, &
            OL3, OL4, OP, OT, NewEggstm1, NewL1tm1, &
            NewL2tm1, NewL3tm1, NewL4tm1, NewPtm1, NewTtm1, &
            Fec, E, L1, L2, L3, L4, P, Te, A, ColdestT, &
            NtGEQ317,NtGEQ00, Bt,r1,x0,x1,CWD,x2,SizeFactor,FebInPopn,&
			Fecout, Eout, L1out, L2out, Pout, Teout, Aout,Flyout, EndMPBPopn,BigTreesOut,SmallTreesOut)
    ! This subroutine simulates the demographic processes
    ! of the mountain pine beetle for a single time step including
    ! oviposition, the egg stage, the four larval instars,
    ! the pupal stage, the teneral adult stage,
    ! the adult stage, the flying adult stage, and then attack.

    implicit none

        ! The subroutine calls the following subroutines.
    Interface

        subroutine Ovipos(Fec, Parents, med, Tmn2, NewEggs)
            real(kind = 8), intent(inout) :: Fec
            real(kind = 8), intent(in) :: Parents
            real(kind = 8), intent(in) :: med
            real(kind = 8), intent(in) :: Tmn2
            real(kind = 8), intent(out) :: NewEggs
        end subroutine Ovipos

        subroutine EPTDev(n, avec, med, mu, sigma, Tmn2, NewEPT, NewEPTtm1, OEPT, EPTcurrent, NewNext)
            integer(kind = 4), intent(in) :: n
            real(kind = 8), intent(in) :: avec(n)
            real(kind = 8), intent(in) :: med
            real(kind = 8), intent(in) :: mu
            real(kind = 8), intent(in) :: sigma
            real(kind = 8), intent(in) :: Tmn2
            real(kind = 8), intent(in) :: NewEPT
            real(kind = 8), intent(inout) :: NewEPTtm1
            complex(kind = 8), intent(inout) :: OEPT(n)
            real(kind = 8), intent(out) :: EPTcurrent
            real(kind = 8), intent(out) :: NewNext
        end subroutine

        subroutine LarvDev(n, avec, med, mu, sigma, NewL, NewLtm1, OL, Lcurrent, NewNext)
            integer(kind = 4), intent(in) :: n
            real(kind = 8), intent(in) :: avec(n)
            real(kind = 8), intent(in) :: med
            real(kind = 8), intent(in) :: mu
            real(kind = 8), intent(in) :: sigma
            real(kind = 8), intent(in) :: NewL
            real(kind = 8), intent(inout) :: NewLtm1
            complex(kind = 8), intent(inout) :: OL(n)
            real(kind = 8), intent(out) :: Lcurrent
            real(kind = 8), intent(out) :: NewNext
        end subroutine

        subroutine AdSR(NewA, Tmn2, Tmx2, Adtm1, FA)
            real(kind = 8), intent(in) :: NewA
            real(kind = 8), intent(in) :: Tmn2
            real(kind = 8), intent(in) :: Tmx2
            real(kind = 8), intent(inout) :: Adtm1
            real(kind = 8), intent(out) :: FA
        end subroutine

        Subroutine ConvolveSR(ItemA, ItemB, AconvB)
            complex(kind = 8), intent(in) :: ItemA(:)
            complex(kind = 8), intent(in) :: ItemB(:)
            complex(kind = 8), intent(out) :: AconvB(:)
        End Subroutine ConvolveSR

        Subroutine LnormPDF(x, mulog, sigmalog, PDF)
            real(kind = 8), intent(in) :: x(:)
            real(kind = 8), intent(in) :: mulog
            real(kind = 8), intent(in) :: sigmalog
            complex(kind = 8), intent(out) :: PDF(:)
        End Subroutine LnormPDF

        subroutine RegniereFunc(TC, TB, DeltaB, TM, DeltaM, omega, psi, DevR)
            real(kind = 8), intent(in) :: TC
            real(kind = 8), intent(in) :: TB
            real(kind = 8), intent(in) :: DeltaB
            real(kind = 8), intent(in) :: TM
            real(kind = 8), intent(in) :: DeltaM
            real(kind = 8), intent(in) :: omega
            real(kind = 8), intent(in) :: psi
            real(kind = 8), intent(out) :: DevR
        end subroutine RegniereFunc

        subroutine FlightFunc(TC, Flying)
            real(kind = 8), intent(in) :: TC
            real(kind = 8), intent(out) :: Flying
        end subroutine

        subroutine BeetleAttack(NtGEQ317,NtGEQ00, Bt, FA, Parents, &
            r1,x0,x1,CWD,x2, SizeFactor,FebInPopn, EndMPBPopn,BigTreesOut,SmallTreesOut)
            real(kind = 8), intent(in) :: NtGEQ317
            real(kind = 8), intent(in) :: NtGEQ00
            real(kind = 8), intent(inout) :: Bt
            real(kind = 8), intent(in) :: FA
            real(kind = 8), intent(out) :: Parents
            real(kind = 8), intent(out) :: BigTreesOut     
            real(kind = 8), intent(out) :: SmallTreesOut               
            real(kind = 8), intent(in) ::x0   ! The intercept on attack 
            real(kind = 8), intent(in) ::CWD  ! The Current Climatic water deficit
            real(kind = 8), intent(in) ::x1   ! x1 controls the impact CWD
            real(kind = 8), intent(in) ::x2   ! x2 controls the impact of size class
            real(kind = 8), intent(in) ::SizeFactor
			real(kind = 8), intent(in) :: r1
               
            real(kind = 8), intent(in) :: FebInPopn
            real(kind = 8), intent(in) :: EndMPBPopn
        end subroutine BeetleAttack

        subroutine BETAISR(A, B, X, Y)
            real(kind = 8), intent(in) :: A
            real(kind = 8), intent(in) :: B
            real(kind = 8), intent(in) :: X
            real(kind = 8), intent(out) :: Y
        end subroutine BETAISR

        subroutine BETACFSR(A, B, X, Z)
            real(kind = 8), intent(in) :: A
            real(kind = 8), intent(in) :: B
            real(kind = 8), intent(in) :: X
            real(kind = 8), intent(out) :: Z
        end subroutine BETACFSR

        subroutine GAMMLNSR(XX, YY)
            real(kind = 8), intent(in) :: XX
            real(kind = 8), intent(out) :: YY
        end subroutine GAMMLNSR

    End Interface

    ! Here are the input and output variables
    real(kind = 8), intent(in) :: Tmax
    real(kind = 8), intent(in) :: Tmin
	real(kind = 8), intent(in) :: TminNext
    real(kind = 8), intent(inout) :: Parents
    real(kind = 8), intent(out) :: FA

    complex(kind = 8), intent(inout) :: OE(2**8)
    complex(kind = 8), intent(inout) :: OL1(2**8)
    complex(kind = 8), intent(inout) :: OL2(2**8)
    complex(kind = 8), intent(inout) :: OL3(2**8)
    complex(kind = 8), intent(inout) :: OL4(2**8)
    complex(kind = 8), intent(inout) :: OP(2**8)
    complex(kind = 8), intent(inout) :: OT(2**8)
    complex(kind = 8), intent(inout) :: OPare(2**8)	
	real(kind = 8), intent(inout) :: NewParentstm1
    real(kind = 8), intent(inout) :: NewEggstm1
    real(kind = 8), intent(inout) :: NewL1tm1
    real(kind = 8), intent(inout) :: NewL2tm1
    real(kind = 8), intent(inout) :: NewL3tm1
    real(kind = 8), intent(inout) :: NewL4tm1
    real(kind = 8), intent(inout) :: NewPtm1
    real(kind = 8), intent(inout) :: NewTtm1
	
	
	real(kind = 8), intent(inout) ::Pare
	real(kind = 8), intent(inout) :: Fec                ! the expected number of pre-eggs at each time
	real(kind = 8), intent(out) :: Fecout               ! the expected number of pre-eggs at each time
    real(kind = 8), intent(inout) :: E                  ! the expected number of eggs at each time
	real(kind = 8), intent(out) :: Eout
    real(kind = 8), intent(inout) :: L1                 ! the expected number of L1 at each time step
    real(kind = 8), intent(out) :: L1out	
    real(kind = 8), intent(inout) :: L2                 ! the expected number of L2 at each time step
    real(kind = 8), intent(out) :: L2out	
	real(kind = 8), intent(inout) :: L3                 ! the expected number of L2 at each time step
	real(kind = 8), intent(inout) :: L4                 ! the expected number of L2 at each time step
    real(kind = 8), intent(inout) :: P                  ! the expected number of pupae at each time step
	real(kind = 8), intent(out) :: Pout 
    real(kind = 8), intent(inout) :: Te                 ! the expected number of tenerals at each time step
	real(kind = 8), intent(out) :: Teout
    real(kind = 8), intent(inout) :: A                  ! the expected number of flying adults at each time step
	real(kind = 8), intent(out) :: Aout
	real(kind = 8), intent(inout) :: ActiveParents 
	real(kind = 8), intent(out) ::Flyout
    real(kind = 8), intent(out) :: BigTreesOut
    real(kind = 8), intent(out) :: SmallTreesOut
    ! The smallest probability of larval winter survival as a function of the lowest temperature to date.
    real(kind = 8), intent(in) :: ColdestT
	real(kind = 8), parameter :: medA = .125
	real(kind = 8), parameter :: sigmaA = .01
	real(kind = 8) :: muA               ! mean of the log transformed random development rate in egg stage
    ! input and output variables
    real(kind = 8), intent(inout) :: NtGEQ317           ! initial susceptible host trees in the >31.7+ cm dbh size class
    real(kind = 8), intent(inout) :: NtGEQ00             ! initial susceptible host trees in the <31.7+ cm dbh size class
    real(kind = 8), intent(inout) :: Bt                 ! beetles that remain in flight from the previous step
    ! input parameters
    real(kind = 8), intent(in) :: r1                    ! controls beetle clustering
    real(kind = 8), intent(in) ::x0   ! The intercept on attack 
    real(kind = 8), intent(in) ::CWD  ! The Current Climatic water deficit
    real(kind = 8), intent(in) ::x1   ! x1 controls the impact CWD
    real(kind = 8), intent(in) ::x2   ! x2 controls the impact of size class4
    real(kind= 8), intent (in) ::SizeFactor
    real(kind = 8), intent(in) :: FebInPopn             ! February insect population
    real(kind = 8), intent(in) :: EndMPBPopn            ! endemic mountain pine beetle population

    !---------------------------------------------------------------------------------
     ! All of the parameters below are internal parameters (internal to the subroutine)

    ! iterator
    integer(kind = 4) :: j
    ! Defining the time step (1 day)

    real(kind = 8), parameter :: deltat = 1.0 ! units are days


	! These are the temperature-development parameter for each life stage
    ! parameters for oviposition rate
    real(kind = 8), parameter :: sigma0 = 0.2458         ! controls rate variability
    real(kind = 8), parameter :: TB0 = 7.750          ! base temperature in degrees C
    real(kind = 8), parameter :: DeltaB0 = 10000
    real(kind = 8), parameter :: TM0 = 34.887           ! max temperature in degree C
    real(kind = 8), parameter :: DeltaM0 = 2.0
    real(kind = 8), parameter :: omega0 = .0085
    real(kind = 8), parameter :: psi0 = 0.013
    
   ! parameters for development rate for eggs
    real(kind = 8), parameter :: sigma1=.3     
    real(kind = 8), parameter :: Rmax1=0.187          
    real(kind = 8), parameter :: Tmax1=34.741
    real(kind = 8), parameter :: Tmin1= 7.267
	real(kind = 8), parameter :: k1= 1056.233
	real(kind = 8), parameter :: dm1 = 7.588
  

    ! parameters for development rate for larvae
    real(kind = 8), parameter :: sigma2=.3     
    real(kind = 8), parameter :: Rmax2=.038            
    real(kind = 8), parameter :: VTmx2=34.701
    real(kind = 8), parameter :: VTmn2= 5.167
	real(kind = 8), parameter :: k2= 766.410
	real(kind = 8), parameter :: dm2 = 7.566
   
    ! parameters for development rate for Prepupae
    real(kind = 8), parameter :: sigma3=.3     
    real(kind = 8), parameter :: Rmax3=0.231          
    real(kind = 8), parameter :: VTmx3=37.599
    real(kind = 8), parameter :: VTmn3= 4.30
	real(kind = 8), parameter :: k3= 431.000
	real(kind = 8), parameter :: dm3 = 11.30
   
    ! parameters for development rate for Pupae
     real(kind = 8), parameter ::sigma4=.3     
    real(kind = 8), parameter :: Rmax4=.062        
    real(kind = 8), parameter :: VTmx4=31.528
    real(kind = 8), parameter :: VTmn4= 5.840
	real(kind = 8), parameter :: k4= 203.096
	real(kind = 8), parameter :: dm4 = 10.016
   
   
    ! parameters for development rate for TA
    real(kind = 8), parameter :: sigma5=.3     
    real(kind = 8), parameter :: Rmax5=.01568408 
    real(kind = 8), parameter :: VTmx5=50.8
    real(kind = 8), parameter :: VTmn5= -0.002603655
	real(kind = 8), parameter :: k5= -0.003262513
	real(kind = 8), parameter :: dm5 = 0.0005461579
   

    ! Here are variables to hold the buffered under bark temperatures
    real(kind = 8) :: Tmean     ! mean temperature at each time step in degrees Celcius
    real(kind = 8) :: Tmin2     ! the buffered under-bark minimum temperature
    real(kind = 8) :: Tmax2     ! the warmer under bark maximum temperature
	real(kind = 8) :: Tmeanall    
    ! Variables that relate to the domain of the lognormal distribution
    integer(kind = 4), parameter :: n = 2**8    ! input variable. Must be specified
    real(kind = 8), parameter :: Mx = 2.0
    real(kind = 8), parameter :: Mn = 1.0e-20
    real(kind = 8), parameter :: da = (Mx - Mn)/(2.0**8.0)

    real(kind = 8) :: avec(n) = (/(Mn + j*da, j=0,n-1)/)          ! vector defining the domain

    ! parameters that change with each iteration
	! The median parameters are calculated for each of 8 periods within the day then summed to 
	! caluclate the daily rate of development. 
	! Each stage will calculate t1-t8
    real(kind = 8) :: med0              ! median development rate in the pre-egg stage
    real(kind = 8) :: med1              ! median development rate in egg stage
    real(kind = 8) :: med2              ! median development rate in L1 stage
    real(kind = 8) :: med3              ! median development rate in L2 stage
    real(kind = 8) :: med4              ! median development rate in L3 stage
    real(kind = 8) :: med5              ! median development rate in L4 stage
    real(kind = 8) :: med6              ! median development rate in Pupa stage
    real(kind = 8) :: med7              ! median development rate in teneral adult stage
    real(kind = 8) :: mu1               ! mean of the log transformed random development rate in egg stage
    real(kind = 8) :: mu2               ! mean of the log transformed random development rate in L1 stage
    real(kind = 8) :: mu3               ! mean of the log transformed random development rate in L2 stage
    real(kind = 8) :: mu4               ! mean of the log transformed random development rate in L3 stage
    real(kind = 8) :: mu5               ! mean of the log transformed random development rate in L4 stage
    real(kind = 8) :: mu6               ! mean of the log transformed random development rate in Pupa stage
    real(kind = 8) :: mu7               ! mean of the log transformed random development rate in teneral adult stage
	real(kind = 8) :: med0t1              !Transfer median development rate in the pre-egg stage
    real(kind = 8) :: med1t1              ! Transfer median development rate in egg stage
    real(kind = 8) :: med2t1              ! Transfer median development rate in L1 stage
    real(kind = 8) :: med3t1              ! Transfer median development rate in L2 stage
    real(kind = 8) :: med4t1              ! Transfer median development rate in L3 stage
    real(kind = 8) :: med5t1              ! Transfer median development rate in L4 stage
    real(kind = 8) :: med6t1              ! Transfer median development rate in Pupa stage
    real(kind = 8) :: med7t1              !Transfer median development rate in teneral adult stage
	real(kind = 8) :: med0t2              !Transfer median development rate in the pre-egg stage
    real(kind = 8) :: med1t2              ! Transfer median development rate in egg stage
    real(kind = 8) :: med2t2              ! Transfer median development rate in L1 stage
    real(kind = 8) :: med3t2              ! Transfer median development rate in L2 stage
    real(kind = 8) :: med4t2              ! Transfer median development rate in L3 stage
    real(kind = 8) :: med5t2              ! Transfer median development rate in L4 stage
    real(kind = 8) :: med6t2              ! Transfer median development rate in Pupa stage
    real(kind = 8) :: med7t2              !Transfer median development rate in teneral adult stage
	real(kind = 8) :: med0t3              !Transfer median development rate in the pre-egg stage
    real(kind = 8) :: med1t3              ! Transfer median development rate in egg stage
    real(kind = 8) :: med2t3              ! Transfer median development rate in L1 stage
    real(kind = 8) :: med3t3              ! Transfer median development rate in L2 stage
    real(kind = 8) :: med4t3              ! Transfer median development rate in L3 stage
    real(kind = 8) :: med5t3              ! Transfer median development rate in L4 stage
    real(kind = 8) :: med6t3              ! Transfer median development rate in Pupa stage
    real(kind = 8) :: med7t3              !Transfer median development rate in teneral adult stage
	real(kind = 8) :: med0t4              !Transfer median development rate in the pre-egg stage
    real(kind = 8) :: med1t4              ! Transfer median development rate in egg stage
    real(kind = 8) :: med2t4              ! Transfer median development rate in L1 stage
    real(kind = 8) :: med3t4              ! Transfer median development rate in L2 stage
    real(kind = 8) :: med4t4              ! Transfer median development rate in L3 stage
    real(kind = 8) :: med5t4              ! Transfer median development rate in L4 stage
    real(kind = 8) :: med6t4              ! Transfer median development rate in Pupa stage
    real(kind = 8) :: med7t4              !Transfer median development rate in teneral adult stage
	real(kind = 8) :: med0t5             !Transfer median development rate in the pre-egg stage
    real(kind = 8) :: med1t5              ! Transfer median development rate in egg stage
    real(kind = 8) :: med2t5              ! Transfer median development rate in L1 stage
    real(kind = 8) :: med3t5              ! Transfer median development rate in L2 stage
    real(kind = 8) :: med4t5              ! Transfer median development rate in L3 stage
    real(kind = 8) :: med5t5              ! Transfer median development rate in L4 stage
    real(kind = 8) :: med6t5              ! Transfer median development rate in Pupa stage
    real(kind = 8) :: med7t5              !Transfer median development rate in teneral adult stage
	real(kind = 8) :: med0t6              !Transfer median development rate in the pre-egg stage
    real(kind = 8) :: med1t6              ! Transfer median development rate in egg stage
    real(kind = 8) :: med2t6              ! Transfer median development rate in L1 stage
    real(kind = 8) :: med3t6              ! Transfer median development rate in L2 stage
    real(kind = 8) :: med4t6              ! Transfer median development rate in L3 stage
    real(kind = 8) :: med5t6              ! Transfer median development rate in L4 stage
    real(kind = 8) :: med6t6              ! Transfer median development rate in Pupa stage
    real(kind = 8) :: med7t6              !Transfer median development rate in teneral adult stage
	real(kind = 8) :: med0t7              !Transfer median development rate in the pre-egg stage
    real(kind = 8) :: med1t7              ! Transfer median development rate in egg stage
    real(kind = 8) :: med2t7              ! Transfer median development rate in L1 stage
    real(kind = 8) :: med3t7              ! Transfer median development rate in L2 stage
    real(kind = 8) :: med4t7              ! Transfer median development rate in L3 stage
    real(kind = 8) :: med5t7              ! Transfer median development rate in L4 stage
    real(kind = 8) :: med6t7              ! Transfer median development rate in Pupa stage
    real(kind = 8) :: med7t7              !Transfer median development rate in teneral adult stage
	real(kind = 8) :: med0t8             !Transfer median development rate in the pre-egg stage
    real(kind = 8) :: med1t8              ! Transfer median development rate in egg stage
    real(kind = 8) :: med2t8              ! Transfer median development rate in L1 stage
    real(kind = 8) :: med3t8              ! Transfer median development rate in L2 stage
    real(kind = 8) :: med4t8              ! Transfer median development rate in L3 stage
    real(kind = 8) :: med5t8              ! Transfer median development rate in L4 stage
    real(kind = 8) :: med6t8              ! Transfer median development rate in Pupa stage
    real(kind = 8) :: med7t8              !Transfer median development rate in teneral adult stage
    ! New individuals that developed into the next life stage in the
    ! time step (these are each scalar values)
    real(kind = 8) :: NewEggs
    real(kind = 8) :: NewL1
    real(kind = 8) :: NewL2
    real(kind = 8) :: NewL3
    real(kind = 8) :: NewL4
    real(kind = 8) :: NewP
    real(kind = 8) :: NewT
    real(kind = 8) :: NewA
	real(kind = 8) :: Cmean             ! mean of current day
    real(kind = 8) :: Cmean2 			! Mean of next half
    real(kind = 8) :: CDif            	! Distance between poles
    real(kind = 8) :: CDif2             ! Distance between next pole
	real(kind = 8) :: fouram            ! 4am Temp
    real(kind = 8) :: sevenam 			! 7am Temp
    real(kind = 8) :: tenam             ! 10am Temp
    real(kind = 8) :: onepm             ! 1pm Temp
    real(kind = 8) :: fourpm            ! 4pm Temp
    real(kind = 8) :: sevenpm           ! 7pm Temp
	real(kind = 8) :: tenpm           ! 10pm Temp
    real(kind = 8) :: oneam             ! 1am Temp
	real(kind = 8) :: Tmin3
    !--------------------------------------------------------------------------------------------------

    !! We need to compute the mean phloem temperature according to the Sine 
    !! relationship we fit to phloem ~ daily air temperatures. 
	!! Each day eight periods are calculated 
    
	!!!Tmean = 0.5*(Tmax + Tmin) + 0.9 + 6.6*(Tmax - Tmin)/(2.0*24.4)
    !Tmax2 = Tmax + 6.6*(Tmax - Tmin)/24.4
    !Tmin2 = Tmin + 1.8
	!Tmeanall=0.5*(Tmax2 + Tmin2) 
	!!Tmin3=TminNext
	! The mean between Tmin(4am) and Tmax(4pm)
	CMean=(Tmax+Tmin)/2
	! The mean between Tmax(4pm) and TminNext (4am the following day)
	CMean2=(Tmax+TminNext)/2
	! The difference between Tmin and Tmax 
	CDif=(Tmax-Tmin)
	! The difference between Tmax and T min Next 
	CDif2=(Tmax-TminNext)
	! This is the implementation of the Ssine Cycle for each of the time periods in the model. 
	fouram=3.8532+(Tmin*.9677)
	sevenam=1.86978+(.93522*(CMean+(.5*CDif*-0.7071068)))
	tenam=(-.4533)+(1.00899*(CMean))
	onepm=(-1.148846)+(.985801*(CMean+(.5*CDif*0.7071068)))
	fourpm=(.0656866)+(.942395*(CMean+(.5*CDif*1)))
	sevenpm=(-0.702683)+(.979172*(CMean2+(.5*CDif2*0.7071068)))
	tenpm=.934665+(.988126*(CMean2))
	oneam=3.2294+(.9842*(CMean2+(.5*CDif2*-0.7071068)))
	Tmax2=fourpm
	Tmin2=fouram
    ! Computing the median development rate for each life stage in this time step
	! This is the sum of the eight periods. 
	! Oviposition
	call RegniereFunc(fouram, TB0, DeltaB0, TM0, DeltaM0, omega0, psi0, med0t1)   ! for pre-eggs
	call RegniereFunc(sevenam, TB0, DeltaB0, TM0, DeltaM0, omega0, psi0, med0t2)   ! for pre-eggs
	call RegniereFunc(tenam, TB0, DeltaB0, TM0, DeltaM0, omega0, psi0, med0t3)   ! for pre-eggs
	call RegniereFunc(onepm, TB0, DeltaB0, TM0, DeltaM0, omega0, psi0, med0t4)   ! for pre-eggs
	call RegniereFunc(fourpm, TB0, DeltaB0, TM0, DeltaM0, omega0, psi0, med0t5)   ! for pre-eggs
	call RegniereFunc(sevenpm, TB0, DeltaB0, TM0, DeltaM0, omega0, psi0, med0t6)   ! for pre-eggs
	call RegniereFunc(tenpm, TB0, DeltaB0, TM0, DeltaM0, omega0, psi0, med0t7)   ! for pre-eggs
	call RegniereFunc(oneam, TB0, DeltaB0, TM0, DeltaM0, omega0, psi0, med0t8)   ! for pre-eggs
    med0=med0t1+med0t2+med0t3+med0t4+med0t5+med0t6+med0t7+med0t8
    ! Eggs 
	call LoganFunc(fouram, Rmax1,Tmax1,Tmin1,k1,dm1, med1t1)   ! 4am
	call LoganFunc(sevenam, Rmax1,Tmax1,Tmin1,k1,dm1, med1t2)   ! 7am
	call LoganFunc(tenam, Rmax1,Tmax1,Tmin1,k1,dm1, med1t3)   ! 10am
	call LoganFunc(onepm, Rmax1,Tmax1,Tmin1,k1,dm1, med1t4)   ! 1 pm
	call LoganFunc(fourpm, Rmax1,Tmax1,Tmin1,k1,dm1, med1t5)   ! 4 pm
	call LoganFunc(sevenpm, Rmax1,Tmax1,Tmin1,k1,dm1, med1t6)   ! 7 pm
	call LoganFunc(tenpm, Rmax1,Tmax1,Tmin1,k1,dm1, med1t7)   ! 10am
	call LoganFunc(oneam, Rmax1,Tmax1,Tmin1,k1,dm1, med1t8)   ! 1 pm
    med1=med1t1+med1t2+med1t3+med1t4+med1t5+med1t6+med1t7+med1t8
	! Larval stage
	call LoganFunc(fouram, Rmax2,VTmx2,VTmn2,k2,dm2, med2t1)   ! for L1
	call LoganFunc(sevenam, Rmax2,VTmx2,VTmn2,k2,dm2, med2t2)   ! for L1
	call LoganFunc(tenam, Rmax2,VTmx2,VTmn2,k2,dm2, med2t3)   ! for L1
	call LoganFunc(onepm, Rmax2,VTmx2,VTmn2,k2,dm2, med2t4)   ! 1 pm
	call LoganFunc(fourpm, Rmax2,VTmx2,VTmn2,k2,dm2, med2t5)   ! for L1
	call LoganFunc(sevenpm, Rmax2,VTmx2,VTmn2,k2,dm2, med2t6)   ! for L1
	call LoganFunc(tenpm, Rmax2,VTmx2,VTmn2,k2,dm2, med2t7)   ! for L1
	call LoganFunc(oneam, Rmax2,VTmx2,VTmn2,k2,dm2, med2t8)   ! for L1
	med2=med2t1+med2t2+med2t3+med2t4+med2t5+med2t6+med2t7+med2t8
	! Pre-pupal 
	call LoganFunc(fouram, Rmax3,VTmx3,VTmn3,k3,dm3, med3t1)   ! for L2
	call LoganFunc(sevenam, Rmax3,VTmx3,VTmn3,k3,dm3, med3t2)   ! for L2
	call LoganFunc(tenam, Rmax3,VTmx3,VTmn3,k3,dm3, med3t3)   ! for L2
	call LoganFunc(onepm, Rmax3,VTmx3,VTmn3,k3,dm3, med3t4)   ! 1 pm
	call LoganFunc(fourpm, Rmax3,VTmx3,VTmn3,k3,dm3, med3t5)   ! for L2
	call LoganFunc(sevenpm, Rmax3,VTmx3,VTmn3,k3,dm3, med3t6)   ! for L2
	call LoganFunc(tenpm, Rmax3,VTmx3,VTmn3,k3,dm3, med3t7)   ! for L2
	call LoganFunc(oneam, Rmax3,VTmx3,VTmn3,k3,dm3, med3t8)   ! for L2
	med3=med3t1+med3t2+med3t3+med3t4+med3t5+med3t6+med3t7+med3t8
	! Pupal 
	call LoganFunc(fouram, Rmax4,VTmx4,VTmn4,k4,dm4, med4t1)   ! for L2
	call LoganFunc(sevenam, Rmax4,VTmx4,VTmn4,k4,dm4, med4t2)   ! for L2
	call LoganFunc(tenam, Rmax4,VTmx4,VTmn4,k4,dm4, med4t3)   ! for L2
	call LoganFunc(onepm, Rmax4,VTmx4,VTmn4,k4,dm4, med4t4)   ! 1 pm
	call LoganFunc(fourpm, Rmax4,VTmx4,VTmn4,k4,dm4, med4t5)   ! for L2
	call LoganFunc(sevenpm, Rmax4,VTmx4,VTmn4,k4,dm4, med4t6)   ! for L2
	call LoganFunc(tenpm, Rmax4,VTmx4,VTmn4,k4,dm4, med4t7)   ! for L2
	call LoganFunc(oneam, Rmax4,VTmx4,VTmn4,k4,dm4, med4t8)   ! for L2
	med4=med4t1+med4t2+med4t3+med4t4+med4t5+med4t6+med4t7+med4t8
    !!!!!! TA
	call LoganFunc(fouram, Rmax5,VTmx5,VTmn5,k5,dm5, med5t1)   ! for L2
	call LoganFunc(sevenam, Rmax5,VTmx5,VTmn5,k5,dm5, med5t2)   ! for L2
	call LoganFunc(tenam, Rmax5,VTmx5,VTmn5,k5,dm5, med5t3)   ! for L2
	call LoganFunc(onepm, Rmax5,VTmx5,VTmn5,k5,dm5, med5t4)   ! 1 pm
	call LoganFunc(fourpm, Rmax5,VTmx5,VTmn5,k5,dm5, med5t5)   ! for L2
	call LoganFunc(sevenpm, Rmax5,VTmx5,VTmn5,k5,dm5, med5t6)   ! for L2
	call LoganFunc(tenpm, Rmax5,VTmx5,VTmn5,k5,dm5, med5t7)   ! for L2
	call LoganFunc(oneam, Rmax5,VTmx5,VTmn5,k5,dm5, med5t8)   ! for L2
	med5=med5t1+med5t2+med5t3+med5t4+med5t5+med5t6+med5t7+med5t8
    
    ! The mu parameter of the lognormal distribution is given by the natural logarithm of the median
    ! development rate.
	muA = log(medA*deltat)
    mu1 = log(med1*deltat)  ! for eggs
    mu2 = log(med2*deltat)  ! for L1
    mu3 = log(med3*deltat)  ! for PP
    mu4 = log(med4*deltat)  ! for P
    mu5 = log(med5*deltat)  ! for TA
   
    !---------------------------------------------------------------------------------------------------------
    ! Now we can simulate each of the life stages by calling the appropriate subroutines
	! Parents die at temperature less than or equal to 12.2 
	
	if(Tmin2 <= -12.2)then        !pg 90 Miller and Keen
		Parents=Parents*.1
		OPare=OPare*.1
		ActiveParents=Parents*.1
	end if
	
	!!! New parents wait 8 days to oviposit.
	call EPTDev(n, avec, medA, muA, sigmaA, Tmin2, Parents, NewParentstm1, OPare, Pare, ActiveParents)
    ! Simulating oviposition:
    call Ovipos(Fec, ActiveParents, med0, Tmin2, NewEggs)
    ! The output of this subroutine (NewEggs) is input for the next subroutine below.
    ! The Ovipos subroutine also updates the scalar value for fecundity.
    ! It takes as input the number of parents which comes from an initial value
    ! or from the BeetleAttack subroutine called at the end of this sequence.
	! 50% of eggs die in the event the minimum temperature is below -20.6C
	if(Tmin2 <= -20.6)then !pg 90 Miller and Keen
		OE=OE*.5
		E=E*.5
	end if
    ! Simulating egg development:
    call EPTDev(n, avec, med1, mu1, sigma1, Tmin2, NewEggs, NewEggstm1, OE, E, NewL1)
    ! The output of this subroutine (NewL1) is input for the next subroutine below.
    ! This updates the aging distribution (OE) and NewEggstm1 and
    ! outputs a scalar for the expected number of eggs.
	
    ! Simulating development of L1 larvae:
    call LarvDev(n, avec, med2, mu2, sigma2, NewL1, NewL1tm1, OL1, L1, NewL2)
    ! The output of this subroutine (NewL2) is input for the next subroutine below.
    ! This updates the aging distribution (OL1) and  NewL1tm1 and
    ! outputs a scalar for the expected number of first instar larvae (L1).

    ! Simulating development of PP:
    call LarvDev(n, avec, med3, mu3, sigma3, NewL2, NewL2tm1, OL2, L2, NewP)
    ! The output of this subroutine (NewL3) is input for the next subroutine below.
    ! This updates the aging distribution (OL2) and  NewL2tm1 and
    ! outputs a scalar for the expected number of second instar larvae (L2).
	
	! These stages are not used in the Western pine beetle model 
    ! Simulating development of L3 larvae:
    ! call LarvDev(n, avec, med4, mu4, sigma4, NewL3, NewL3tm1, OL3, L3, NewL4)
    ! The output of this subroutine (NewL4) is input for the next subroutine below.
    ! This updates the aging distribution (OL3) and  NewL3tm1 and
    ! outputs a scalar for the expected number of third instar larvae (L3).

    ! Simulating development of L4 larvae:
    !call LarvDev(n, avec, med5, mu5, sigma5, NewL4, NewL4tm1, OL4, L4, NewP)
    ! The output of this subroutine (NewP) is input for the next subroutine below.
    ! This updates the aging distribution (OL4) and  NewL4tm1 and
    ! outputs a scalar for the expected number of fourth instar larvae (L4).

    ! Applying winter mortality  to the larval stages. 
    NewP = NewP/(1.0 + exp(-(ColdestT + 29.22791)/2.0))
    ! NewP = NewP*0.5
	
    ! Simulating pupal development:
    call EPTDev(n, avec, med4, mu4, sigma4, Tmin2, NewP, NewPtm1, OP, P, NewT)
    ! The output of this subroutine (NewT) is input for the next subroutine below.
    ! This updates the aging distribution (OP) and  NewPtm1 and
    ! outputs a scalar for the expected number of pupae (P).
	! The Pupal population is reduced to 10% in the event the temperature is below 
	! -20.6 C
	if(Tmin2 <= -20.6)then !pg 90 Miller and Keen
		OP=P*.1
		P=P*.1
	end if
    ! Simulating teneral adult development:
	
    call EPTDev(n, avec, med5, mu5, sigma5, Tmin2, NewT, NewTtm1, OT, Te, NewA)
    ! The output of this subroutine (NewA) is input for the next subroutine below.
    ! This updates the aging distribution (OT) and  NewTtm1 and
    ! outputs a scalar for the expected number of teneral adults (Te).
	! The teneral adult population is reduced to 10% in the event the population is 
	! reduced below -12.2C 
	if(Tmin2 <= -12.2)then !pg 90 Miller and Keen
		OT=OT*.1
		Te=Te*.1
	end if
    ! Simulating adult flight
	! This updates the expected number of adults (A) and flying adults (FA).
    call AdSR(NewA, Tmin2, Tmax, A, FA)

	! This updates the output variables for model monitoring. 
	Flyout=FA+Bt
	Fecout=Fec 
	Eout=E 
	L1out=L1 
	L2out=L2 
	Pout=P 
	Teout=Te
	Aout=A
	
    ! If temperatures are cold, all of the flying beetles die. This prevents them
    ! from killing trees when they should not be. The value of -10*C is from Miller and Keen 
    if(Tmin < -10.0)then
        Bt = 0.0
        FA = 0.0
    end if
	
    ! Simulating the attack of host trees	
    call BeetleAttack(NtGEQ317, NtGEQ00, Bt, FA, Parents, r1,x0,x1,CWD,x2,&
	                  SizeFactor,FebInPopn,EndMPBPopn,BigTreesOut,SmallTreesOut)
    ! This updates the density of trees in each of the size classes, and the density of beetles that remain in
    ! flight and outputs a number of parents that will start the oviposition process.

End Subroutine BeetleSim

!=================================================================================================================
subroutine BeetleAttack(NtGEQ317, NtGEQ00,Bt, FA, Parents, &
		   r1,x0,x1,CWD,x2,SizeFactor,FebInPopn,EndMPBPopn,BigTreesOut,SmallTreesOut)
    ! In this subroutine we use the TDIA model to determine the number of trees that are killed using the relationship between the 
	! the number of beetles in flight and the amount of drought experienced by the trees and their size. 
	
	
	
    implicit none

    interface

        subroutine BETAISR(A, B, X, Y)
            real(kind = 8), intent(in) :: A
            real(kind = 8), intent(in) :: B
            real(kind = 8), intent(in) :: X
            real(kind = 8), intent(out) :: Y
        end subroutine BETAISR

        subroutine BETACFSR(A, B, X, Z)
            real(kind = 8), intent(in) :: A
            real(kind = 8), intent(in) :: B
            real(kind = 8), intent(in) :: X
            real(kind = 8), intent(out) :: Z
        end subroutine BETACFSR

        subroutine GAMMLNSR(XX, YY)
            real(kind = 8), intent(in) :: XX
            real(kind = 8), intent(out) :: YY
        end subroutine GAMMLNSR

    end interface

    ! input and output variables.
    ! Tree density in size classes per ha 
    real(kind = 8), intent(in) :: NtGEQ317                ! initial susceptible host trees in the 31.7 cm dbh size class
    real(kind = 8), intent(in) :: NtGEQ00			      ! initial susceptible host trees in the 10-31.7 cm dbh size class
    real(kind = 8), intent(inout) :: Bt                   ! Cumulative attacking beetles

    ! input variable
    real(kind = 8), intent(in) :: FA                      ! Adults that just started to fly in this time step
    ! output variables
    real(kind = 8), intent(out) :: Parents                ! the density of beetles that entered trees killed in this time step
    real(kind = 8), intent(out) :: BigTreesOut
    real(kind = 8), intent(out) :: SmallTreesOut
    ! input parameters
    real(kind = 8), intent(in) ::x0                       ! The intercept on attack 
    real(kind = 8), intent(in) ::CWD                      ! The Current 4-yr lagged standard percipitaiton index
    real(kind = 8), intent(in) ::x1  					  ! x1 controls the impact CWD
    real(kind = 8), intent(in) ::x2   					  ! x2 controls the impact of size class
    real(kind = 8), intent(in) ::SizeFactor               ! The preference control for host size. 
    real(kind = 8), intent(in) :: r1                      ! controls beetle clustering
    real(kind = 8), intent(in) :: FebInPopn               ! February insect population (THIS CAN BE REMOVED 
    real(kind = 8), intent(in) :: EndMPBPopn              ! endemic Western pine beetle population threshold

    ! Here are internal variables and parameters
    real(kind = 8) :: timestep = 1.0            ! one day time step
    real(kind = 8) :: Btp1                      ! an updated value for the beetles
    real(kind = 8) :: Atp1GEQ20                 ! attacking beetles
    real(kind = 8) :: phi1                      ! controls beetle attack success
    real(kind = 8) :: phi2                      ! controls beetle attack success     
    real(kind = 8) :: Ntp1GEQ20                 ! updated susceptible host trees in the 20+ cm dbh size class
    real(kind = 8) :: Ntp1GEQ00
    real(kind = 8) :: Ptp1GEQ20                 ! first argument of the incomplete beta function for tree loss (clustering param)
    real(kind = 8) :: shape2a                   ! updated parent beetles the 20+ cm dbh size class
    real(kind = 8) :: shape2b                   ! updated parent beetles the 20+ cm dbh size class     
    real(kind = 8) :: shape1                	! second argument of the incomplete beta function for tree loss (attack thresh + 1)
    real(kind = 8) :: q1                  		! third argument of the incomplete beta function for tree loss (see below)
    real(kind = 8) :: q2        				! third argument of the incomplete beta function for tree loss (see below)    
    real(kind = 8) :: q3
    real(kind = 8) :: Ps
    real(kind = 8) :: ibeta        				! output of the incomplete beta function (subroutine)
    real(kind = 8) :: ibeta2                    ! output of the incomplete beta function (subroutine)     
    real(kind = 8) :: ibeta3 
    real(kind = 8) :: S                    		! output of the incomplete beta function (subroutine)
    real(kind = 8) :: Killed
    real(kind = 8) :: Btp2
    real(kind = 8) :: Btp3    
     ! setting the values of the arguments of the incomplete beta function
    Btp1 = Bt + FA
    ! The control on the aggregation of beetles 
	shape1 = r1
	! phi is the number of beetles necessary to kill a tree (Controlled by the size of the tree x2 and the influence
	! of drought x1. 
	! The number of beetles needed to kill the large size class (>31.6cm)
    phi1 = exp(x0+x1*CWD+x2*1)
	! The number of beetles needed to kill the smaller size class (>10cm and <31.6cm)
    phi2 = exp(x0+x1*CWD)
	
	! If there are large trees 
	if(NtGEQ317 >0.0) then 
	    !!! q1 is the mean number of beetles to effectively aggregate to a tree. 
		q1 = (r1*Btp1)/(NtGEQ317)    
		!!! The surviving number of trees is 1- the number of trees killed 
		ibeta=1-((q1)/(1+phi1))
		! The expected number of beetles to attack each the larger tree type is determined
		Ps=(1-ibeta)**SizeFactor
		! No less than 50% of beetles  will attack the preferred host (No-preference)
		if(Ps <.5) then
			Ps=.5
		end if 
		!!!! Two mortality submodels here use the expected number of beetles to attack each 
		! size class to to determine the mortality level s
		! Big Trees
		! Number of beetles to attack large trees
		Btp2 = Btp1*Ps
		! q2 the effective number of beetles to attack each tree
		q2 = (r1*Btp2)/(NtGEQ317)    
		! The surviving number of trees 
		ibeta2=1-(q2/(1+phi1)) 
		! The remaining population 
		Ntp1GEQ20 = (ibeta2 * NtGEQ317)
		
		!Little trees
		! Number of beetles to attack smaller trees
		Btp3 = Btp1 *(1-Ps)    
		q3 = (r1*Btp3)/(NtGEQ00)     
		ibeta3=1-(q3/(1+phi2))
	    ! The surviving Small Trees
		Ntp1GEQ00 = (ibeta3 * NtGEQ00)	
		
		
		! The resulting number of parents that attacked each tree
		Ptp1GEQ20 = phi1*(NtGEQ317-Ntp1GEQ20)+phi2*(NtGEQ00-Ntp1GEQ00)
		!! If trees did not die, remove attempted flight (Beetles die in attempted attacks)
		if (Ntp1GEQ20< NtGEQ317) then
			Btp2=0
		end if 	
		if(NtGEQ00< NtGEQ00) then
			Btp3=0
		end if
		Btp1=Btp2+Btp3
		!------------------------------------------------------------------------------------------------
		!! Handling nas
		!! Handling nas
		if (isnan(Btp1)) then
			Btp1=0
		end if 
		! Ensure more trees aren't killed than exist. 
		if(Ntp1GEQ00<0)then
			Ntp1GEQ00=0
		end if 
		if(Ntp1GEQ20<0)then
			Ntp1GEQ20=0
		end if
		if(isnan(Ntp1GEQ20))then
			Ntp1GEQ20=0
		end if 
		if(isnan(Ntp1GEQ00))then
			Ntp1GEQ00=0
		end if 
		
		!! Ensure that we don't end up with NA parents 
		if (isnan(phi1*(NtGEQ317-Ntp1GEQ20))) then
			Ptp1GEQ20=phi2*(NtGEQ00-Ntp1GEQ00)
		end if 
		if (isnan(phi2*(NtGEQ00-Ntp1GEQ00))) then
			Ptp1GEQ20=phi1*(NtGEQ317-Ntp1GEQ20)
		end if 
		if (isnan(Ptp1GEQ20)) then
			Ptp1GEQ20=0
		end if 
		! Now I update all of the state variables. This depends on whether the population is endemic or not.
		! when populations are in the endemic phase, they only attack weakened
		! trees that are already functionally dead from other causes.
		if(FebInPopn > EndMPBPopn)then
			Parents = Ptp1GEQ20
			Bt = Btp1
			BigTreesOut = Ntp1GEQ20
			SmallTreesOut = Ntp1GEQ00
			else
			! Under the endemic scenario beetles do not kill trees.
				Parents = Ptp1GEQ20
				Bt = Btp1
				BigTreesOut = NtGEQ317
				SmallTreesOut = NtGEQ00
		end if
	end if 	
end subroutine

!======================================================================================================
subroutine Ovipos(Fec, Parents, med, Tmn2, NewEggs)

    ! This subroutine simulate oviposition by parent beetles.
    ! The fecundity variable is the number of eggs remaining in the pool.
    ! med is median (temperature-dependent) oviposition rate.
    ! NewEggs are eggs that were laid in the time step.

    implicit none

    ! Here are the input and output variables
    real(kind = 8), intent(inout) :: Fec            ! fecundity (eggs remaining to be laid)
    real(kind = 8), intent(in) :: Parents           ! number of parents doing the ovipositing
    real(kind = 8), intent(in) :: med               ! median oviposition rate
    real(kind = 8), intent(in) :: Tmn2              ! minimum temperature under the bark
    real(kind = 8), intent(out) :: NewEggs          ! Eggs laid in the time step

    ! internal parameters
    real(kind = 8), parameter :: fmax = 48*.5        !Miller and Keen and Stephen and Dalsten

    real(kind = 8), parameter :: netp = 0.142857   ! This is 1 minus the probability of mortality of a variety of causes
    !real(kind = 8), parameter :: netp = 0.1625       ! This is 1 minus the probability of mortality of a variety of causes

    ! Aplying winter mortality to egg laying adults
     if(Tmn2 <= -12.2)then  !pg 90 Miller and Keen 
        Fec = 0.0
    end if

    ! Computing new eggs. Note this has to be done before updating the
    ! Fec variable below. I pro-actively apply mortality to the new eggs
    ! that they will experience as they develop through the juvenile stages.
    NewEggs = Fec*(1.0 - exp(-med))*netp

    ! Simulating oviposition: (Fec represents the number of eggs remaining)
    ! Below I assume that half of the individuals that fly and
    ! lay eggs are females (after tree induced colonization mortality)
    ! and each female lays an initial clutch of 82 eggs.
    Fec = Parents*fmax + Fec*exp(-med)

end subroutine Ovipos

!======================================================================================================
subroutine EPTDev(n, avec, med, mu, sigma, Tmn2, NewEPT, NewEPTtm1, OEPT, EPTcurrent, NewNext)

    ! This subroutine advances egg, pupal or teneral adult development and
    ! returns the number of individuals that move into the next stage
    ! (NewNext) in the time step as well as a count of the number
    ! of individuals currently in the stage (EPTcurrent).
    ! It can be used for the eggs, pupae or teneral adults, but it takes
    ! different life stage-dependent parameters (med, mu, sigma).

    ! This subroutine takes as input, n, which gives the size of the age domain,
    ! avec, which is the aging domain itself. The temperature-dependent
    ! median development rate (med) for the life stage, the corresponding
    ! mean of the log-normally distributed aging rate (mu),
    ! and the scale parameter of the log-normal distribution (sigma).
    ! Note that med, mu, and sigma are life stage-specific. The algorithm
    ! also takes as input the minimum daily temperature (Tmn2) which it uses
    ! to assign temperature-dependent mortality.
    ! The subroutine takes as input the number of new individuals in the
    ! current step (NewEPT), and from the previous step (NewEPTtm1).
    ! NewEPTtm1, OEPT is input and output and is the
    ! distribution of the accumulated physiological age in the life stage.

    implicit none

    ! Here are the input and output variables
    integer(kind = 4), intent(in) :: n              ! size of aging domain
    real(kind = 8), intent(in) :: avec(n)           ! The aging domain
    real(kind = 8), intent(in) :: med               ! median development rate
    real(kind = 8), intent(in) :: mu                ! mean of the log transformed random development rate
    real(kind = 8), intent(in) :: sigma             ! scale parameter of the log-normally distributed rate
    real(kind = 8), intent(in) :: Tmn2              ! the buffered under-bark minimum temperature
    real(kind = 8), intent(in) :: NewEPT            ! New individuals (just developed from previous stage)
    real(kind = 8), intent(inout) :: NewEPTtm1      ! New individuals from the previous time step
    complex(kind = 8), intent(inout) :: OEPT(n)     ! Distribution of physiological age
    real(kind = 8), intent(out) :: EPTcurrent       ! Number of individuals currently in the stage
    real(kind = 8), intent(out) :: NewNext          ! New that just developed out of the current stage into the next one

    ! Here are variables related to the convolution
    complex(kind = 8) :: OldEPT(n)                  ! copy of the distribution of physiological age
    complex(kind = 8) :: AconvB(n)                  ! An array to hold the convolution result
    complex(kind = 8) PDF(n)                        ! An array to hold the aging kernel

    ! The subroutine calls the ConvolveSR subroutine to do the
    ! convolution and the LnormPDF subroutine to create the
    ! log-normally distributed development rate.
    Interface

        Subroutine ConvolveSR(ItemA, ItemB, AconvB)
            complex(kind = 8), intent(in) :: ItemA(:)
            complex(kind = 8), intent(in) :: ItemB(:)
            complex(kind = 8), intent(out) :: AconvB(:)
        End Subroutine ConvolveSR

        Subroutine LnormPDF(x, mulog, sigmalog, PDF)
            real(kind = 8), intent(in) :: x(:)
            real(kind = 8), intent(in) :: mulog
            real(kind = 8), intent(in) :: sigmalog
            complex(kind = 8), intent(out) :: PDF(:)
        End Subroutine LnormPDF

    End Interface

    ! To compute the number of new individuals in the next stage,
    ! we need to know how many old individuals there were
    ! in the previous step.
    OldEPT = OEPT

    ! I use the -18 threshold to kill all eggs, pupae, or teneral
    ! adults as described in Regniere 2015.
    !if(Tmn2 <= -15.0)then
    !    OldEPT = 0.0
    !end if

    if(med > 0.0)then
        ! Computing the aging kernel
        call LnormPDF(avec, mu, sigma, PDF)
        PDF = PDF/sum(PDF)     ! ensuring that it sums to one

        ! Doing the convolution
        call ConvolveSR(PDF, OldEPT, AconvB)
        ! I assign all new individuals
        ! to the first age in the next stage.
        OEPT = AConvB + NewEPTtm1*PDF
        NewEPTtm1 = NewEPT

        ! Computing eggs
        EPTcurrent = sum(real(OEPT(1:128))) + NewEPT

        ! Computing NewNext: New individuals that developed into
        ! the next stage in this time step are individuals that
        ! exceeded the breakpoint of the current life stage (128).
        NewNext = sum(real(OEPT(129:n))) - sum(real(OldEPT(129:n)))

        else
            OEPT = OldEPT
            NewEPTtm1 = NewEPTtm1 + NewEPT
            EPTcurrent = sum(real(OEPT(1:128))) + NewEPTtm1

            ! Computing NewNext: if no development occurs there are no new individuals
            NewNext = 0.0
        end if

        ! Just to make sure that silly things don't happen
        if(NewNext < 0.0 .or. isnan(NewNext))then
            NewNext = 0.0
        end if

end subroutine EPTDev

!======================================================================================================
subroutine LarvDev(n, avec, med, mu, sigma, NewL, NewLtm1, OL, Lcurrent, NewNext)

    ! This subroutine advances larval development and returns
    ! the number of larvae in the current instar that advance
    ! to the next life stage (NewNext) in the time step as well as
    ! a count of the number larvae currently in the stage (Lcurrent).
    ! It can be used for any of the four mountain pine beetle larval
    ! instars but takes different instar-dependent parameters (med, mu, sigma).

    ! This subroutine takes as input, n, which gives the size of the age domain,
    ! avec, which is the aging domain itself. The temperature-dependent
    ! median development rate (med) for the larval stage
    ! (each instar has its own), the corresponding mean of the
    ! log-normally distributed aging rate (mu), and the scale parameter
    ! of the log-normal distribution (sigma). Note that mu and sigma
    ! are also larval instar-specific. The subroutine takes
    ! as input the number of new larvae in the current step (NewL)
    ! and from the previous step (NewLtm1). NewLtm1, and OL are input
    ! and output. OL is the distribution of the
    ! accumulated physiological age in the life stage.

    implicit none

    ! Here are the input and output variables
    integer(kind = 4), intent(in) :: n              ! size of aging domain
    real(kind = 8), intent(in) :: avec(n)           ! The aging domain
    real(kind = 8), intent(in) :: med               ! median development rate
    real(kind = 8), intent(in) :: mu                ! mean of the log transformed random development rate
    real(kind = 8), intent(in) :: sigma             ! scale parameter of the log-normally distributed rate
    real(kind = 8), intent(in) :: NewL              ! New larvae (just developed from the previous stage)
    real(kind = 8), intent(inout) :: NewLtm1        ! New larvae from the previous time step
    complex(kind = 8), intent(inout) :: OL(n)       ! Distribution of physiological age
    real(kind = 8), intent(out) :: Lcurrent         ! Number of larvae
    real(kind = 8), intent(out) :: NewNext          ! New individuals in the next stage (just developed out of current stage)

    ! Here are variables related to the convolution
    complex(kind = 8) :: OldL(n)                    ! copy of the distribution of physiological age
    complex(kind = 8) :: AconvB(n)                  ! An array to hold the convolution result
    complex(kind = 8) PDF(n)                        ! An array to hold the aging kernel

    ! The subroutine calls the ConvolveSR subroutine to do the
    ! convolution and the LnormPDF subroutine to create the
    ! log-normally distributed development rate.
    Interface

        Subroutine ConvolveSR(ItemA, ItemB, AconvB)
            complex(kind = 8), intent(in) :: ItemA(:)
            complex(kind = 8), intent(in) :: ItemB(:)
            complex(kind = 8), intent(out) :: AconvB(:)
        End Subroutine ConvolveSR

        Subroutine LnormPDF(x, mulog, sigmalog, PDF)
            real(kind = 8), intent(in) :: x(:)
            real(kind = 8), intent(in) :: mulog
            real(kind = 8), intent(in) :: sigmalog
            complex(kind = 8), intent(out) :: PDF(:)
        End Subroutine LnormPDF

    End Interface

    ! To compute the number of new individuals in the next stage,
    ! we need to know how many old individuals there were
    ! in the previous step.
    OldL = OL

    if(med > 0.0)then
        ! computing the aging kernel
        call LnormPDF(avec, mu, sigma, PDF)
        PDF = PDF/sum(PDF)     ! ensuring that it sums to one

        ! doing the convolution
        call ConvolveSR(PDF, OldL, AconvB)
        ! I assign all new individuals
        ! to the first age in the next stage.
        OL = AConvB + NewLtm1*PDF
        NewLtm1 = NewL

        ! computing Lcurrent
        Lcurrent = sum(real(OL(1:128))) + NewL

        ! Computing individuals that developed into the next stage
        ! in this time step: These are individuals that exceeded
        ! the breakpoint (128).
        NewNext = sum(real(OL(129:n))) - sum(real(OldL(129:n)))

    else
        OL = OldL
        NewLtm1 = NewLtm1 + NewL
        Lcurrent = sum(real(OL(1:128))) + NewLtm1

        ! Computing new individuals in the next stage
        NewNext = 0.0
    end if

    ! Just to make sure that silly things don't happen
    if(NewNext < 0.0 .or. isnan(NewNext))then
        NewNext = 0.0
    end if

end subroutine LarvDev

!======================================================================================================
subroutine AdSR(NewA, Tmn2, Tmx2, Adtm1, FA)

    ! This subroutine keeps track of the number of mature adults
    ! that have not yet flown and outputs flown adults. This subroutine
    ! calls the PropFly subroutine to compute what proportion of adults
    ! fly as a function of maximum under-bark temperature.

    ! It takes as input a scalar (NewA) that comes from individuals that
    ! have developed out of the teneral adult stage (from the EPTDev subroutine),
    ! the scalar representing the minimum under bark temperature (Tmn2),
    ! the scalar representing the maximum under bark temperature(Tmx2),
    ! and a scalar number of adults in the previous time step Adtm1.
    ! It returns an updated scalar for the number of adults (Adtm1),
    ! and a scalar of the number of adults that flew in this timestep (FA).

    ! The number of adults that flew in the time step becomes input for a
    ! differential equation representing the attack of trees by the flying
    ! adults.

    implicit none

    ! Here are the input and output variables
    real(kind = 8), intent(in) :: NewA              ! New adults
    real(kind = 8), intent(in) :: Tmn2              ! minimum temperature under the bark
    real(kind = 8), intent(in) :: Tmx2              ! maximum temperature under the bark
    real(kind = 8), intent(inout) :: Adtm1          ! Number of adult individuals
    real(kind = 8), intent(out) :: FA               ! Flying beetles

    ! An internal variable
    real(kind = 8) :: PropFly                       ! The proportion of adult beetles that fly

    ! The subroutine calls the FlightFunc subroutine.
    Interface
        Subroutine FlightFunc(TC, Flying)
            real(kind = 8), intent(in) :: TC
            real(kind = 8), intent(out) :: Flying
        End Subroutine FlightFunc
    End Interface

    ! I use the -12.2 Milller AND Keen pg 90
    if(Tmn2 <= -12.2)then
        Adtm1 = 0.0
    end if

    ! If it's warm enough, beetles fly. I'm using an empirical
    ! model based on the work of McCambridge (1971)
    call FlightFunc(Tmx2, PropFly)

    ! Computing adults that flew in this time step. Note that
    ! this needs to be done BEFORE updating the remaining adults.
    FA = Adtm1*PropFly

    ! Uptdating Adtm1
    Adtm1 = Adtm1*(1.0-PropFly) + NewA

end subroutine AdSR

!=============================================================================================================================
Subroutine ConvolveSR(ItemA, ItemB, AconvB)
    ! This subroutine convolves to arrays of the same dimensions (ItemA and ItemB)
    ! by first Fast Fourier transforming them, then multiplying their
    ! Fast Fourier transforms and then inverse Fast Fourier
    ! transforming them. However, to prevent circular creep (due to circular domain
    ! assumption for FFT) I pad the arrays on the right hand side with zeros.

    implicit none

    complex(kind = 8), intent(in) :: ItemA(:)
    complex(kind = 8), intent(in) :: ItemB(:)
    complex(kind = 8), intent(out) :: AconvB(:)

    ! n -> integer, length of sequence to be transformed--This is an input only variable
    ! inc -> integer, the distance (increment) between discretized locations in the sequence--This is an input only variable
    ! c -> complex, the sequence of data to be transformed (vector)--This is an input/output variable
    ! lenc -> integer, the dimension of the R array (typically it will be the same as n)--This is an input only variable
    ! wsave -> An array that stores the prime factorization of n along with necessary trigonometric functions
    ! the wsave array does not have to be explicitly defined as calling the initialization function rfft1i
    ! will do this for you---This is an input only variable
    ! lensav -> integer, this is the dimension of the wsave array (must be at least 2*n + int(log(real(n))) + 4)--This is an input variable only
    ! work -> the work array. This also does not have to be explicitly defined as calling the initialization function
    ! rfft1i will do this for you
    ! lenwrk -> integer, the dimension of the work array (must be at least 2*n)--This is an input only variable
    ! ier -> error flag (0-success, 1-lenr not big enough, 2-lensav not big enough, 3-lenwrk not big enough)
    !--This is an output only variable.

    ! Variables that are arguments of the fft functions
    integer(kind = 4), parameter :: m = 2**8        ! input variable. Must be specified
    integer(kind = 4), parameter :: Twom = 2*m      ! Because I use a pad of size n, we double n
    integer(kind = 4), parameter :: inc = 1         ! input variable. Must be specified
    integer(kind = 4), parameter :: lenc = Twom     ! input variable. Must be specified

    ! Computational set-up to do the Fast Fourier transforms
    ! variable related to the wsave and work arrays
    integer(kind = 4) :: lensav
    integer(kind = 4) :: lenwrk

    integer(kind = 4) :: ier                        ! output variable (not specified)

    ! Values of these two variables are obtained by initializing (calling rffti1)
    real(kind = 8), allocatable, dimension( : ) :: work
    real(kind = 8), allocatable, dimension( : ) :: wsave

    complex(kind = 8) :: Convolved(Twom)

    ! First I do the padding with zeros on the right hand side of the two input arrays
    complex(kind = 8) :: PaddedItemA(Twom)
    complex(kind = 8) :: PaddedItemB(Twom)
    PaddedItemA(1:m) = ItemA
    PaddedItemB(1:m) = ItemB
    PaddedItemA(m+1:Twom) = 0.0
    PaddedItemB(m+1:Twom) = 0.0

    ! dimension of wsave and work arrays.
    lensav = 2*Twom + int ( log ( real ( Twom, kind = 4 ) ) / log ( 2.0E+00 ) ) + 4
    lenwrk = 2*Twom

    ! Allocating memory (for the fft)
    allocate(work(1:lenwrk))
    allocate(wsave(1:lensav))

    ! initializing th fft subroutines
    ! This will assign values to the work and wsave arrays
    call cfft1i(Twom, wsave, lensav, ier)

    ! Now I fast Fourier transform the padded variables (which NEEED to be complex)
    ! Doing the fast Fourier transform (forward)
    call cfft1f(Twom, inc, PaddedItemA, lenc, wsave, lensav, work, lenwrk, ier)
    call cfft1f(Twom, inc, PaddedItemB, lenc, wsave, lensav, work, lenwrk, ier)

    ! Doing the convolution
    ! We have to multiply by 2*n to counteract the normalization
    ! that the Fortran version of the fft algorithm does automatically
    Convolved = Twom*PaddedItemA*PaddedItemB

    ! Inverse Fourier transforming the convolution (backward)
    call cfft1b(Twom, inc, Convolved, lenc, wsave, lensav, work, lenwrk, ier)

    ! To not loose individuals that developed past the threshold I add
    AconvB = Convolved(1:m)

    ! In this last funky step I take all of the individuals that developed into the
    ! padded region and add them to the last possible spot in the array so that they
    ! do not affect development in future steps.
    AconvB(m) = AconvB(m) + sum(Convolved(m+1:Twom))

    deallocate(work)
    deallocate(wsave)

End Subroutine ConvolveSR

!================================================================================================
Subroutine LnormPDF(x, mulog, sigmalog, PDF)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Here is a subroutine that produces a lognormal probability density function
    ! over a range of values.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! I use the assumed-shape convention for the 1-D array
    implicit none

    real(kind = 8), intent(in) :: x(:)
    real(kind = 8), intent(in) :: mulog
    real(kind = 8), intent(in) :: sigmalog
    complex(kind = 8), intent(out) :: PDF(:)

    real(kind = 8), parameter :: pi = 3.1415927

    PDF = 1.0/(x*sigmalog*sqrt(2.0*pi))*exp(-1.0/(2.0*sigmalog**2.0)*(log(x)-mulog)**2.0)

end Subroutine LnormPDF

!========================================================================================
subroutine RegniereFunc(TC, TB, DeltaB, TM, DeltaM, omega, psi, DevR)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The RegniereFunc represents the Regniere function for
    ! temperature-dependent insect development (Regniere et al. 2012).
    ! TC is temperature in degrees Celcius.
    ! TB, DeltaB, TM,DeltaM, omega, psi are parameters
    ! The subroutine returns DevR, the median
    ! temperature-dependent development rate.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(kind = 8), intent(in) :: TC
    real(kind = 8), intent(in) :: TB
    real(kind = 8), intent(in) :: DeltaB
    real(kind = 8), intent(in) :: TM
    real(kind = 8), intent(in) :: DeltaM
    real(kind = 8), intent(in) :: omega
    real(kind = 8), intent(in) :: psi
    real(kind = 8), intent(out) :: DevR

    ! We start by defining it as zero in case the condition below does not hold
    DevR = 0.0

    ! Computing the development rate
    if(TC >= TB .and. TC <= TM) then
        DevR = psi*(exp(omega*(TC - TB)) - (TM - TC)/(TM - TB)*exp(-omega*(TC - TB)/DeltaB) &
        - (TC - TB)/(TM - TB)*exp(omega*(TM - TB) - (TM - TC)/DeltaM))
    end if
	if(DevR <0) then
	   DevR=0
	end if 
end subroutine RegniereFunc
!========================================================================================
subroutine LoganFunc(TC,Rmax,Tmax,Tmin,k,dm,DevR)
     implicit none
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 ! The function for stage growth as defined in 
	 ! Logan, J. A. (1988). Toward an expert system for development of pest simulation models. Environmental Entomology, 17(2), 359-376.
	 ! TC is the current temperature in degrees Celsius
	 ! Rmax, Tmax, Tmin, k, dm are the the fit parameters
	 ! dm is the median rate of development. 
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     real(kind = 8), intent(in) :: TC
     real(kind = 8), intent(in) :: Rmax
     real(kind = 8), intent(in) :: Tmax
	 real(kind = 8), intent(in) :: Tmin
     real(kind = 8), intent(in) :: k
	 real(kind = 8), intent(in) :: dm
     real(kind = 8), intent(out) :: DevR
     real(kind = 8), parameter :: Euler =2.71828
     ! We start by defining it as zero in case the condition below does not hold
     DevR = 0.0
     if(TC <= TMax) then
        DevR =Rmax*(((((TC-Tmin)**2)/(((TC-Tmin)**2)+k))-(Euler**(-(Tmax-(TC-Tmin))/dm))))
     end if
	 if(DevR <0) then
	    DevR=0
	 end if 
end subroutine LoganFunc


!========================================================================================
subroutine TaylorFunc(TC,Rmax,Tmax,Tdelta,DevR)
     implicit none
     real(kind = 8), intent(in) :: TC
     real(kind = 8), intent(in) :: Rmax
     real(kind = 8), intent(in) :: Tmax
     real(kind = 8), intent(in) :: Tdelta
     real(kind = 8), intent(out) :: DevR
     real(kind = 8), parameter :: Euler =2.71828
     ! We start by defining it as zero in case the condition below does not hold
     DevR = 0.0
     if(TC <= TMax) then
        DevR =Rmax*(Euler**(-.5*(((TC-TMax)/Tdelta)**2)))
     end if
	 if(DevR <0) then
	    DevR=0
	 end if 
end subroutine TaylorFunc
!========================================================================================

subroutine FlightFunc(TC, Flying)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Here's a subroutine that outputs the proportion of beetles that fly
    ! as a function of air temperature. The subroutine is based on an
    ! empirical fit of a polynomial function to the data of
    ! McCambridge (1971).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(kind = 8), intent(in) :: TC
    real(kind = 8), intent(out) :: Flying

    ! We start by defining it as zero in case the condition below does not hold
    Flying = 0.0

    ! Calculated in the R script based on Gaylord 2014
    if(TC >= 18.6 .and. TC <= 38.9) then
        Flying = .02524971 + (-.1952827)*(TC)+(.03021792)*(TC**2.0) + (-.00167174)*(TC**3.0) +  &
         (4.096269e-05)*(TC**4.0) +(-3.632886e-07)*(TC**5.0)
    end if

    if(Flying > 1.0) then
        Flying = 1.0
    end if

    if(Flying < 0.0) then
        Flying = 0.0
    end if

end subroutine FlightFunc

!===============================================================================
subroutine BETAISR(A, B, X, Y)
    !! returns the normalized incomplete Beta function Ix(a,b) in the variable Y

    implicit none

    Interface
        ! Computes the numerator of the normalized incomplete beta function.
        subroutine BETACFSR(A, B, X, Z)
            real(kind = 8), intent(in) :: A
            real(kind = 8), intent(in) :: B
            real(kind = 8), intent(in) :: X
            real(kind = 8), intent(out) :: Z
        end subroutine BETACFSR

        ! returns the value Ln(Gamma(XX)) for XX>0.
        subroutine GAMMLNSR(XX, YY)
            real(kind = 8), intent(in) :: XX
            real(kind = 8), intent(out) :: YY
        end subroutine GAMMLNSR

    End Interface

    ! Arguments of the main subroutine
    real(kind = 8), intent(in) :: A
    real(kind = 8), intent(in) :: B
    real(kind = 8), intent(in) :: X
    real(kind = 8), intent(out) :: Y

    ! Here are internal variables used by the subroutines that are called
    ! by the main subroutine
    real(kind = 8) :: BT
    real(kind = 8) :: Z
    real(kind = 8) :: YY1
    real(kind = 8) :: YY2
    real(kind = 8) :: YY3

    If(X.LT.0.d0.OR.X.GT.1.d0) then
        print *,' BETAI: Bad argument X (must be 0<=X<=1).'
        Return
    End If

    !--------------------------------------------------------------------------------------------
    ! Here's the part of the subroutine that does the work.
    If(X.EQ.0.d0.OR.X.EQ.1.d0) then
        BT = 0.d0
        else
            ! Filling in the variables that are natural logarithms of the gamma function
            call GAMMLNSR(A+B, YY1)
            call GAMMLNSR(A, YY2)
            call GAMMLNSR(B, YY3)

            BT = exp(YY1 - YY2 - YY3 + A*log(X) + B*log(1.d0-X))
    End If

    If(X.LT.(A+1.d0)/(A+B+2.d0)) then
        call BETACFSR(A, B, X, Z)
        Y=BT*Z/A
        return
        else
            call BETACFSR(B, A, 1.d0-X, Z)
            Y=1.d0-BT*Z/B
    End If

end subroutine BETAISR

!===============================================================================
subroutine BETACFSR(A, B, X, Z)
    !! Returns the numerator of the regularized incomplete beta function in the variable Z.

    implicit none

    ! Arguments
    real(kind = 8), intent(in) :: A
    real(kind = 8), intent(in) :: B
    real(kind = 8), intent(in) :: X
    real(kind = 8), intent(out) :: Z

    ! Here are internal variables
    integer, parameter :: ITMAX=100
    integer :: M
    real(kind = 8), parameter :: EPS=3.d-7
    real(kind = 8) AM,BM,AZ,QAB,QAP,QAM,BZ,EM,TEM,D,AP,BP,APP,BPP
    real(kind = 8) AOLD

    ! Setting the values of the variables
    AM=1.d0
    BM=1.d0
    AZ=1.d0
    QAB=A+B
    QAP=A+1.d0
    QAM=A-1.d0
    BZ=1.d0-QAB*X/QAP

    ! Here's where we do the work
    Do M = 1, ITMAX
        EM=M
        TEM=EM+EM
        D=EM*(B-M)*X/((QAM+TEM)*(A+TEM))
        AP=AZ+D*AM
        BP=BZ+D*BM
        D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM))
        APP=AP+D*AZ
        BPP=BP+D*BZ
        AOLD=AZ
        AM=AP/BPP
        BM=BP/BPP
        AZ=APP/BPP
        BZ=1.d0
        If(abs(AZ-AOLD).LT.EPS*abs(AZ)) Goto 1
    End Do

    print *,' BETACF: A or B too big, or ITMAX too small.'
    Return
    1 Z = AZ
    Return

end subroutine BETACFSR

!===============================================================================
subroutine GAMMLNSR(XX, YY)

    implicit none

    ! Arguments
    real(kind = 8), intent(in) :: XX
    real(kind = 8), intent(out) :: YY

    ! Internal variables
    integer :: j
    real(kind = 8) COF(6),STP,HALF,ONE,FPF,X,TMP,SER

    ! Assigning values
    data COF,STP /76.18009173d0,-86.50532033d0,24.01409822d0, &
       -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/

    data HALF,ONE,FPF /0.5d0,1.d0,5.5d0/

    X=XX-ONE
    TMP=X+FPF
    TMP=(X+HALF)*log(TMP)-TMP
    SER=ONE

    !-------------------------------------------------------------------
    ! Here's the part of the subroutine that does the work
    Do j = 1, 6
        X=X+ONE
        SER=SER+COF(J)/X
    End Do

    YY = TMP + log(STP*SER)
    return

end subroutine GAMMLNSR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Below are the required subroutines to do complex number FFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine c1f2kb ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F2KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,2)
  real ( kind = 8 ) ch(in2,l1,2,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,1,2)

  if ( ido <= 1 .and. na /= 1 ) then

    do k = 1, l1
      chold1 = cc(1,k,1,1)+cc(1,k,1,2)
      cc(1,k,1,2) = cc(1,k,1,1)-cc(1,k,1,2)
      cc(1,k,1,1) = chold1
      chold2 = cc(2,k,1,1)+cc(2,k,1,2)
      cc(2,k,1,2) = cc(2,k,1,1)-cc(2,k,1,2)
      cc(2,k,1,1) = chold2
    end do

    return

  end if

  do k = 1, l1
    ch(1,k,1,1) = cc(1,k,1,1)+cc(1,k,1,2)
    ch(1,k,2,1) = cc(1,k,1,1)-cc(1,k,1,2)
    ch(2,k,1,1) = cc(2,k,1,1)+cc(2,k,1,2)
    ch(2,k,2,1) = cc(2,k,1,1)-cc(2,k,1,2)
  end do

  do i = 2, ido
    do k = 1, l1
      ch(1,k,1,i) = cc(1,k,i,1)+cc(1,k,i,2)
      tr2 = cc(1,k,i,1)-cc(1,k,i,2)
      ch(2,k,1,i) = cc(2,k,i,1)+cc(2,k,i,2)
      ti2 = cc(2,k,i,1)-cc(2,k,i,2)
      ch(2,k,2,i) = wa(i,1,1)*ti2+wa(i,1,2)*tr2
      ch(1,k,2,i) = wa(i,1,1)*tr2-wa(i,1,2)*ti2
    end do
  end do

  return
end

!==================================================================================
subroutine c1f2kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F2KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    13 May 2013
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,2)
  real ( kind = 8 ) ch(in2,l1,2,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,1,2)

  if ( ido <= 1 ) then

    sn = 1.0D+00 / real ( 2 * l1, kind = 8 )

    if ( na == 1 ) then

      do k = 1, l1
        ch(1,k,1,1) = sn*(cc(1,k,1,1)+cc(1,k,1,2))
        ch(1,k,2,1) = sn*(cc(1,k,1,1)-cc(1,k,1,2))
        ch(2,k,1,1) = sn*(cc(2,k,1,1)+cc(2,k,1,2))
        ch(2,k,2,1) = sn*(cc(2,k,1,1)-cc(2,k,1,2))
      end do

    else

      do k = 1, l1
        chold1 = sn*(cc(1,k,1,1)+cc(1,k,1,2))
        cc(1,k,1,2) = sn*(cc(1,k,1,1)-cc(1,k,1,2))
        cc(1,k,1,1) = chold1
        chold2 = sn*(cc(2,k,1,1)+cc(2,k,1,2))
        cc(2,k,1,2) = sn*(cc(2,k,1,1)-cc(2,k,1,2))
        cc(2,k,1,1) = chold2
      end do

    end if

  else

    do k = 1, l1
      ch(1,k,1,1) = cc(1,k,1,1)+cc(1,k,1,2)
      ch(1,k,2,1) = cc(1,k,1,1)-cc(1,k,1,2)
      ch(2,k,1,1) = cc(2,k,1,1)+cc(2,k,1,2)
      ch(2,k,2,1) = cc(2,k,1,1)-cc(2,k,1,2)
    end do

    do i = 2, ido
      do k = 1, l1
        ch(1,k,1,i) = cc(1,k,i,1)+cc(1,k,i,2)
        tr2 = cc(1,k,i,1)-cc(1,k,i,2)
        ch(2,k,1,i) = cc(2,k,i,1)+cc(2,k,i,2)
        ti2 = cc(2,k,i,1)-cc(2,k,i,2)
        ch(2,k,2,i) = wa(i,1,1)*ti2-wa(i,1,2)*tr2
        ch(1,k,2,i) = wa(i,1,1)*tr2+wa(i,1,2)*ti2
      end do
    end do

  end if

  return
end

!==================================================================================
subroutine c1f3kb ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F3KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    13 May 2013
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,3)
  real ( kind = 8 ) ch(in2,l1,3,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ), parameter :: taui =  0.866025403784439D+00
  real ( kind = 8 ), parameter :: taur = -0.5D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,2,2)

  if ( ido <= 1 .and. na /= 1 ) then

    do k = 1, l1
      tr2 = cc(1,k,1,2)+cc(1,k,1,3)
      cr2 = cc(1,k,1,1)+taur*tr2
      cc(1,k,1,1) = cc(1,k,1,1)+tr2
      ti2 = cc(2,k,1,2)+cc(2,k,1,3)
      ci2 = cc(2,k,1,1)+taur*ti2
      cc(2,k,1,1) = cc(2,k,1,1)+ti2
      cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
      ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
      cc(1,k,1,2) = cr2-ci3
      cc(1,k,1,3) = cr2+ci3
      cc(2,k,1,2) = ci2+cr3
      cc(2,k,1,3) = ci2-cr3
    end do

  else

    do k = 1, l1
      tr2 = cc(1,k,1,2)+cc(1,k,1,3)
      cr2 = cc(1,k,1,1)+taur*tr2
      ch(1,k,1,1) = cc(1,k,1,1)+tr2
      ti2 = cc(2,k,1,2)+cc(2,k,1,3)
      ci2 = cc(2,k,1,1)+taur*ti2
      ch(2,k,1,1) = cc(2,k,1,1)+ti2
      cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
      ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
      ch(1,k,2,1) = cr2-ci3
      ch(1,k,3,1) = cr2+ci3
      ch(2,k,2,1) = ci2+cr3
      ch(2,k,3,1) = ci2-cr3
    end do

    do i = 2, ido
      do k = 1, l1
        tr2 = cc(1,k,i,2)+cc(1,k,i,3)
        cr2 = cc(1,k,i,1)+taur*tr2
        ch(1,k,1,i) = cc(1,k,i,1)+tr2
        ti2 = cc(2,k,i,2)+cc(2,k,i,3)
        ci2 = cc(2,k,i,1)+taur*ti2
        ch(2,k,1,i) = cc(2,k,i,1)+ti2
        cr3 = taui*(cc(1,k,i,2)-cc(1,k,i,3))
        ci3 = taui*(cc(2,k,i,2)-cc(2,k,i,3))
        dr2 = cr2-ci3
        dr3 = cr2+ci3
        di2 = ci2+cr3
        di3 = ci2-cr3
        ch(2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
        ch(1,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
        ch(2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
        ch(1,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
      end do
    end do

  end if

  return
end

!==================================================================================
subroutine c1f3kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F3KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    13 May 2013
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,3)
  real ( kind = 8 ) ch(in2,l1,3,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ), parameter :: taui = -0.866025403784439D+00
  real ( kind = 8 ), parameter :: taur = -0.5D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,2,2)

  if ( ido <= 1 ) then

    sn = 1.0D+00 / real ( 3 * l1, kind = 8 )

    if ( na /= 1 ) then

      do k = 1, l1
        tr2 = cc(1,k,1,2)+cc(1,k,1,3)
        cr2 = cc(1,k,1,1)+taur*tr2
        cc(1,k,1,1) = sn*(cc(1,k,1,1)+tr2)
        ti2 = cc(2,k,1,2)+cc(2,k,1,3)
        ci2 = cc(2,k,1,1)+taur*ti2
        cc(2,k,1,1) = sn*(cc(2,k,1,1)+ti2)
        cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
        ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
        cc(1,k,1,2) = sn*(cr2-ci3)
        cc(1,k,1,3) = sn*(cr2+ci3)
        cc(2,k,1,2) = sn*(ci2+cr3)
        cc(2,k,1,3) = sn*(ci2-cr3)
      end do

    else

      do k = 1, l1
        tr2 = cc(1,k,1,2)+cc(1,k,1,3)
        cr2 = cc(1,k,1,1)+taur*tr2
        ch(1,k,1,1) = sn*(cc(1,k,1,1)+tr2)
        ti2 = cc(2,k,1,2)+cc(2,k,1,3)
        ci2 = cc(2,k,1,1)+taur*ti2
        ch(2,k,1,1) = sn*(cc(2,k,1,1)+ti2)
        cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
        ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
        ch(1,k,2,1) = sn*(cr2-ci3)
        ch(1,k,3,1) = sn*(cr2+ci3)
        ch(2,k,2,1) = sn*(ci2+cr3)
        ch(2,k,3,1) = sn*(ci2-cr3)
      end do

    end if

  else

    do k = 1, l1
      tr2 = cc(1,k,1,2)+cc(1,k,1,3)
      cr2 = cc(1,k,1,1)+taur*tr2
      ch(1,k,1,1) = cc(1,k,1,1)+tr2
      ti2 = cc(2,k,1,2)+cc(2,k,1,3)
      ci2 = cc(2,k,1,1)+taur*ti2
      ch(2,k,1,1) = cc(2,k,1,1)+ti2
      cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
      ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
      ch(1,k,2,1) = cr2-ci3
      ch(1,k,3,1) = cr2+ci3
      ch(2,k,2,1) = ci2+cr3
      ch(2,k,3,1) = ci2-cr3
    end do

    do i = 2, ido
      do k = 1, l1
        tr2 = cc(1,k,i,2)+cc(1,k,i,3)
        cr2 = cc(1,k,i,1)+taur*tr2
        ch(1,k,1,i) = cc(1,k,i,1)+tr2
        ti2 = cc(2,k,i,2)+cc(2,k,i,3)
        ci2 = cc(2,k,i,1)+taur*ti2
        ch(2,k,1,i) = cc(2,k,i,1)+ti2
        cr3 = taui*(cc(1,k,i,2)-cc(1,k,i,3))
        ci3 = taui*(cc(2,k,i,2)-cc(2,k,i,3))
        dr2 = cr2-ci3
        dr3 = cr2+ci3
        di2 = ci2+cr3
        di3 = ci2-cr3
        ch(2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
        ch(1,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
        ch(2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
        ch(1,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
      end do
    end do

  end if

  return
end

!==================================================================================
subroutine c1f4kb ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F4KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,4)
  real ( kind = 8 ) ch(in2,l1,4,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa(ido,3,2)

  if ( ido <= 1 .and. na /= 1 ) then

      do k = 1, l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,4)-cc(2,k,1,2)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,2)-cc(1,k,1,4)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         cc(1,k,1,1) = tr2+tr3
         cc(1,k,1,3) = tr2-tr3
         cc(2,k,1,1) = ti2+ti3
         cc(2,k,1,3) = ti2-ti3
         cc(1,k,1,2) = tr1+tr4
         cc(1,k,1,4) = tr1-tr4
         cc(2,k,1,2) = ti1+ti4
         cc(2,k,1,4) = ti1-ti4
      end do

  else

      do k = 1, l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,4)-cc(2,k,1,2)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,2)-cc(1,k,1,4)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         ch(1,k,1,1) = tr2+tr3
         ch(1,k,3,1) = tr2-tr3
         ch(2,k,1,1) = ti2+ti3
         ch(2,k,3,1) = ti2-ti3
         ch(1,k,2,1) = tr1+tr4
         ch(1,k,4,1) = tr1-tr4
         ch(2,k,2,1) = ti1+ti4
         ch(2,k,4,1) = ti1-ti4
      end do

      do i = 2, ido
         do k = 1, l1
            ti1 = cc(2,k,i,1)-cc(2,k,i,3)
            ti2 = cc(2,k,i,1)+cc(2,k,i,3)
            ti3 = cc(2,k,i,2)+cc(2,k,i,4)
            tr4 = cc(2,k,i,4)-cc(2,k,i,2)
            tr1 = cc(1,k,i,1)-cc(1,k,i,3)
            tr2 = cc(1,k,i,1)+cc(1,k,i,3)
            ti4 = cc(1,k,i,2)-cc(1,k,i,4)
            tr3 = cc(1,k,i,2)+cc(1,k,i,4)
            ch(1,k,1,i) = tr2+tr3
            cr3 = tr2-tr3
            ch(2,k,1,i) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(1,k,2,i) = wa(i,1,1)*cr2-wa(i,1,2)*ci2
            ch(2,k,2,i) = wa(i,1,1)*ci2+wa(i,1,2)*cr2
            ch(1,k,3,i) = wa(i,2,1)*cr3-wa(i,2,2)*ci3
            ch(2,k,3,i) = wa(i,2,1)*ci3+wa(i,2,2)*cr3
            ch(1,k,4,i) = wa(i,3,1)*cr4-wa(i,3,2)*ci4
            ch(2,k,4,i) = wa(i,3,1)*ci4+wa(i,3,2)*cr4
         end do
      end do

  end if

  return
end

!==================================================================================
subroutine c1f4kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F4KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,4)
  real ( kind = 8 ) ch(in2,l1,4,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa(ido,3,2)

      if ( 1 < ido ) go to 102
      sn = 1.0D+00 / real ( 4 * l1, kind = 8 )
      if (na == 1) go to 106

      do k = 1, l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,2)-cc(2,k,1,4)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,4)-cc(1,k,1,2)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         cc(1,k,1,1) = sn*(tr2+tr3)
         cc(1,k,1,3) = sn*(tr2-tr3)
         cc(2,k,1,1) = sn*(ti2+ti3)
         cc(2,k,1,3) = sn*(ti2-ti3)
         cc(1,k,1,2) = sn*(tr1+tr4)
         cc(1,k,1,4) = sn*(tr1-tr4)
         cc(2,k,1,2) = sn*(ti1+ti4)
         cc(2,k,1,4) = sn*(ti1-ti4)
      end do

      return

  106 do 107 k = 1, l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,2)-cc(2,k,1,4)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,4)-cc(1,k,1,2)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         ch(1,k,1,1) = sn*(tr2+tr3)
         ch(1,k,3,1) = sn*(tr2-tr3)
         ch(2,k,1,1) = sn*(ti2+ti3)
         ch(2,k,3,1) = sn*(ti2-ti3)
         ch(1,k,2,1) = sn*(tr1+tr4)
         ch(1,k,4,1) = sn*(tr1-tr4)
         ch(2,k,2,1) = sn*(ti1+ti4)
         ch(2,k,4,1) = sn*(ti1-ti4)
  107 continue

      return

  102 do 103 k = 1, l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,2)-cc(2,k,1,4)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,4)-cc(1,k,1,2)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         ch(1,k,1,1) = tr2+tr3
         ch(1,k,3,1) = tr2-tr3
         ch(2,k,1,1) = ti2+ti3
         ch(2,k,3,1) = ti2-ti3
         ch(1,k,2,1) = tr1+tr4
         ch(1,k,4,1) = tr1-tr4
         ch(2,k,2,1) = ti1+ti4
         ch(2,k,4,1) = ti1-ti4
  103 continue
      do 105 i = 2, ido
         do 104 k = 1, l1
            ti1 = cc(2,k,i,1)-cc(2,k,i,3)
            ti2 = cc(2,k,i,1)+cc(2,k,i,3)
            ti3 = cc(2,k,i,2)+cc(2,k,i,4)
            tr4 = cc(2,k,i,2)-cc(2,k,i,4)
            tr1 = cc(1,k,i,1)-cc(1,k,i,3)
            tr2 = cc(1,k,i,1)+cc(1,k,i,3)
            ti4 = cc(1,k,i,4)-cc(1,k,i,2)
            tr3 = cc(1,k,i,2)+cc(1,k,i,4)
            ch(1,k,1,i) = tr2+tr3
            cr3 = tr2-tr3
            ch(2,k,1,i) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(1,k,2,i) = wa(i,1,1)*cr2+wa(i,1,2)*ci2
            ch(2,k,2,i) = wa(i,1,1)*ci2-wa(i,1,2)*cr2
            ch(1,k,3,i) = wa(i,2,1)*cr3+wa(i,2,2)*ci3
            ch(2,k,3,i) = wa(i,2,1)*ci3-wa(i,2,2)*cr3
            ch(1,k,4,i) = wa(i,3,1)*cr4+wa(i,3,2)*ci4
            ch(2,k,4,i) = wa(i,3,1)*ci4-wa(i,3,2)*cr4
  104    continue
  105 continue

  return
end

!==================================================================================
subroutine c1f5kb ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F5KB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,5)
  real ( kind = 8 ) ch(in2,l1,5,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: ti11 =  0.9510565162951536D+00
  real ( kind = 8 ), parameter :: ti12 =  0.5877852522924731D+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ), parameter :: tr11 =  0.3090169943749474D+00
  real ( kind = 8 ), parameter :: tr12 = -0.8090169943749474D+00
  real ( kind = 8 ) wa(ido,4,2)

      if ( 1 < ido .or. na == 1) go to 102

      do k = 1, l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         chold1 = cc(1,k,1,1)+tr2+tr3
         chold2 = cc(2,k,1,1)+ti2+ti3
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cc(1,k,1,1) = chold1
         cc(2,k,1,1) = chold2
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         cc(1,k,1,2) = cr2-ci5
         cc(1,k,1,5) = cr2+ci5
         cc(2,k,1,2) = ci2+cr5
         cc(2,k,1,3) = ci3+cr4
         cc(1,k,1,3) = cr3-ci4
         cc(1,k,1,4) = cr3+ci4
         cc(2,k,1,4) = ci3-cr4
         cc(2,k,1,5) = ci2-cr5
      end do

      return

  102 do 103 k = 1, l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         ch(1,k,1,1) = cc(1,k,1,1)+tr2+tr3
         ch(2,k,1,1) = cc(2,k,1,1)+ti2+ti3
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2,1) = cr2-ci5
         ch(1,k,5,1) = cr2+ci5
         ch(2,k,2,1) = ci2+cr5
         ch(2,k,3,1) = ci3+cr4
         ch(1,k,3,1) = cr3-ci4
         ch(1,k,4,1) = cr3+ci4
         ch(2,k,4,1) = ci3-cr4
         ch(2,k,5,1) = ci2-cr5
  103 continue

      do 105 i = 2, ido
         do 104 k = 1, l1
            ti5 = cc(2,k,i,2)-cc(2,k,i,5)
            ti2 = cc(2,k,i,2)+cc(2,k,i,5)
            ti4 = cc(2,k,i,3)-cc(2,k,i,4)
            ti3 = cc(2,k,i,3)+cc(2,k,i,4)
            tr5 = cc(1,k,i,2)-cc(1,k,i,5)
            tr2 = cc(1,k,i,2)+cc(1,k,i,5)
            tr4 = cc(1,k,i,3)-cc(1,k,i,4)
            tr3 = cc(1,k,i,3)+cc(1,k,i,4)
            ch(1,k,1,i) = cc(1,k,i,1)+tr2+tr3
            ch(2,k,1,i) = cc(2,k,i,1)+ti2+ti3
            cr2 = cc(1,k,i,1)+tr11*tr2+tr12*tr3
            ci2 = cc(2,k,i,1)+tr11*ti2+tr12*ti3
            cr3 = cc(1,k,i,1)+tr12*tr2+tr11*tr3
            ci3 = cc(2,k,i,1)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(1,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
            ch(2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
            ch(1,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
            ch(2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
            ch(1,k,4,i) = wa(i,3,1)*dr4-wa(i,3,2)*di4
            ch(2,k,4,i) = wa(i,3,1)*di4+wa(i,3,2)*dr4
            ch(1,k,5,i) = wa(i,4,1)*dr5-wa(i,4,2)*di5
            ch(2,k,5,i) = wa(i,4,1)*di5+wa(i,4,2)*dr5
  104    continue
  105 continue

  return
end

!==================================================================================
subroutine c1f5kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F5KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,5)
  real ( kind = 8 ) ch(in2,l1,5,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: ti11 = -0.9510565162951536D+00
  real ( kind = 8 ), parameter :: ti12 = -0.5877852522924731D+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ), parameter :: tr11 =  0.3090169943749474D+00
  real ( kind = 8 ), parameter :: tr12 = -0.8090169943749474D+00
  real ( kind = 8 ) wa(ido,4,2)

      if ( 1 < ido ) go to 102
      sn = 1.0D+00 / real ( 5 * l1, kind = 8 )
      if (na == 1) go to 106

      do k = 1, l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         chold1 = sn*(cc(1,k,1,1)+tr2+tr3)
         chold2 = sn*(cc(2,k,1,1)+ti2+ti3)
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cc(1,k,1,1) = chold1
         cc(2,k,1,1) = chold2
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         cc(1,k,1,2) = sn*(cr2-ci5)
         cc(1,k,1,5) = sn*(cr2+ci5)
         cc(2,k,1,2) = sn*(ci2+cr5)
         cc(2,k,1,3) = sn*(ci3+cr4)
         cc(1,k,1,3) = sn*(cr3-ci4)
         cc(1,k,1,4) = sn*(cr3+ci4)
         cc(2,k,1,4) = sn*(ci3-cr4)
         cc(2,k,1,5) = sn*(ci2-cr5)
      end do

      return

  106 do 107 k = 1, l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         ch(1,k,1,1) = sn*(cc(1,k,1,1)+tr2+tr3)
         ch(2,k,1,1) = sn*(cc(2,k,1,1)+ti2+ti3)
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2,1) = sn*(cr2-ci5)
         ch(1,k,5,1) = sn*(cr2+ci5)
         ch(2,k,2,1) = sn*(ci2+cr5)
         ch(2,k,3,1) = sn*(ci3+cr4)
         ch(1,k,3,1) = sn*(cr3-ci4)
         ch(1,k,4,1) = sn*(cr3+ci4)
         ch(2,k,4,1) = sn*(ci3-cr4)
         ch(2,k,5,1) = sn*(ci2-cr5)
  107 continue

      return

  102 do 103 k = 1, l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         ch(1,k,1,1) = cc(1,k,1,1)+tr2+tr3
         ch(2,k,1,1) = cc(2,k,1,1)+ti2+ti3
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2,1) = cr2-ci5
         ch(1,k,5,1) = cr2+ci5
         ch(2,k,2,1) = ci2+cr5
         ch(2,k,3,1) = ci3+cr4
         ch(1,k,3,1) = cr3-ci4
         ch(1,k,4,1) = cr3+ci4
         ch(2,k,4,1) = ci3-cr4
         ch(2,k,5,1) = ci2-cr5
  103 continue

      do 105 i = 2, ido
         do 104 k = 1, l1
            ti5 = cc(2,k,i,2)-cc(2,k,i,5)
            ti2 = cc(2,k,i,2)+cc(2,k,i,5)
            ti4 = cc(2,k,i,3)-cc(2,k,i,4)
            ti3 = cc(2,k,i,3)+cc(2,k,i,4)
            tr5 = cc(1,k,i,2)-cc(1,k,i,5)
            tr2 = cc(1,k,i,2)+cc(1,k,i,5)
            tr4 = cc(1,k,i,3)-cc(1,k,i,4)
            tr3 = cc(1,k,i,3)+cc(1,k,i,4)
            ch(1,k,1,i) = cc(1,k,i,1)+tr2+tr3
            ch(2,k,1,i) = cc(2,k,i,1)+ti2+ti3
            cr2 = cc(1,k,i,1)+tr11*tr2+tr12*tr3
            ci2 = cc(2,k,i,1)+tr11*ti2+tr12*ti3
            cr3 = cc(1,k,i,1)+tr12*tr2+tr11*tr3
            ci3 = cc(2,k,i,1)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(1,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
            ch(2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
            ch(1,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
            ch(2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
            ch(1,k,4,i) = wa(i,3,1)*dr4+wa(i,3,2)*di4
            ch(2,k,4,i) = wa(i,3,1)*di4-wa(i,3,2)*dr4
            ch(1,k,5,i) = wa(i,4,1)*dr5+wa(i,4,2)*di5
            ch(2,k,5,i) = wa(i,4,1)*di5-wa(i,4,2)*dr5
  104    continue
  105 continue

  return
end

!==================================================================================
subroutine c1fgkb ( ido, ip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa )

!*****************************************************************************80
!
!! C1FGKB is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lid

  real ( kind = 8 ) cc(in1,l1,ip,ido)
  real ( kind = 8 ) cc1(in1,lid,ip)
  real ( kind = 8 ) ch(in2,l1,ido,ip)
  real ( kind = 8 ) ch1(in2,lid,ip)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ki
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) na
  real ( kind = 8 ) wa(ido,ip-1,2)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

      ipp2 = ip+2
      ipph = (ip+1)/2

      do ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
      end do

      do 111 j=2,ipph
         jc = ipp2-j
         do 112 ki=1,lid
            ch1(1,ki,j) =  cc1(1,ki,j)+cc1(1,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)-cc1(1,ki,jc)
            ch1(2,ki,j) =  cc1(2,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(2,ki,jc)
  112    continue
  111 continue

      do 118 j=2,ipph
         do 117 ki=1,lid
            cc1(1,ki,1) = cc1(1,ki,1)+ch1(1,ki,j)
            cc1(2,ki,1) = cc1(2,ki,1)+ch1(2,ki,j)
  117    continue
  118 continue

      do 116 l=2,ipph
         lc = ipp2-l
         do 113 ki=1,lid
            cc1(1,ki,l) = ch1(1,ki,1)+wa(1,l-1,1)*ch1(1,ki,2)
            cc1(1,ki,lc) = wa(1,l-1,2)*ch1(1,ki,ip)
            cc1(2,ki,l) = ch1(2,ki,1)+wa(1,l-1,1)*ch1(2,ki,2)
            cc1(2,ki,lc) = wa(1,l-1,2)*ch1(2,ki,ip)
  113    continue
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = mod((l-1)*(j-1),ip)
            war = wa(1,idlj,1)
            wai = wa(1,idlj,2)
            do 114 ki=1,lid
               cc1(1,ki,l) = cc1(1,ki,l)+war*ch1(1,ki,j)
               cc1(1,ki,lc) = cc1(1,ki,lc)+wai*ch1(1,ki,jc)
               cc1(2,ki,l) = cc1(2,ki,l)+war*ch1(2,ki,j)
               cc1(2,ki,lc) = cc1(2,ki,lc)+wai*ch1(2,ki,jc)
  114       continue
  115    continue
  116 continue

      if( 1 < ido .or. na == 1) go to 136

      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ki=1,lid
            chold1 = cc1(1,ki,j)-cc1(2,ki,jc)
            chold2 = cc1(1,ki,j)+cc1(2,ki,jc)
            cc1(1,ki,j) = chold1
            cc1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
            cc1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
            cc1(1,ki,jc) = chold2
  119    continue
  120 continue
      return

  136 do 137 ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
  137 continue

      do 135 j=2,ipph
         jc = ipp2-j
         do 134 ki=1,lid
            ch1(1,ki,j) = cc1(1,ki,j)-cc1(2,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
            ch1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
  134    continue
  135 continue

      if (ido == 1) then
        return
      end if

      do 131 i=1,ido
         do 130 k = 1, l1
            cc(1,k,1,i) = ch(1,k,i,1)
            cc(2,k,1,i) = ch(2,k,i,1)
  130    continue
  131 continue

      do 123 j=2,ip
         do 122 k = 1, l1
            cc(1,k,j,1) = ch(1,k,1,j)
            cc(2,k,j,1) = ch(2,k,1,j)
  122    continue
  123 continue

      do 126 j=2,ip
         do 125 i = 2, ido
            do 124 k = 1, l1
               cc(1,k,j,i) = wa(i,j-1,1)*ch(1,k,i,j) &
                            -wa(i,j-1,2)*ch(2,k,i,j)
               cc(2,k,j,i) = wa(i,j-1,1)*ch(2,k,i,j) &
                            +wa(i,j-1,2)*ch(1,k,i,j)
  124       continue
  125    continue
  126 continue

  return
end

!==================================================================================
subroutine c1fgkf ( ido, ip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa )

!*****************************************************************************80
!
!! C1FGKF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lid

  real ( kind = 8 ) cc(in1,l1,ip,ido)
  real ( kind = 8 ) cc1(in1,lid,ip)
  real ( kind = 8 ) ch(in2,l1,ido,ip)
  real ( kind = 8 ) ch1(in2,lid,ip)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ki
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) wa(ido,ip-1,2)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

      ipp2 = ip+2
      ipph = (ip+1)/2
      do ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
      end do

      do 111 j=2,ipph
         jc = ipp2-j
         do 112 ki=1,lid
            ch1(1,ki,j) =  cc1(1,ki,j)+cc1(1,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)-cc1(1,ki,jc)
            ch1(2,ki,j) =  cc1(2,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(2,ki,jc)
  112    continue
  111 continue

      do 118 j=2,ipph
         do 117 ki=1,lid
            cc1(1,ki,1) = cc1(1,ki,1)+ch1(1,ki,j)
            cc1(2,ki,1) = cc1(2,ki,1)+ch1(2,ki,j)
  117    continue
  118 continue

      do 116 l=2,ipph
         lc = ipp2-l
         do 113 ki=1,lid
            cc1(1,ki,l) = ch1(1,ki,1)+wa(1,l-1,1)*ch1(1,ki,2)
            cc1(1,ki,lc) = -wa(1,l-1,2)*ch1(1,ki,ip)
            cc1(2,ki,l) = ch1(2,ki,1)+wa(1,l-1,1)*ch1(2,ki,2)
            cc1(2,ki,lc) = -wa(1,l-1,2)*ch1(2,ki,ip)
  113    continue
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = mod((l-1)*(j-1),ip)
            war = wa(1,idlj,1)
            wai = -wa(1,idlj,2)
            do 114 ki=1,lid
               cc1(1,ki,l) = cc1(1,ki,l)+war*ch1(1,ki,j)
               cc1(1,ki,lc) = cc1(1,ki,lc)+wai*ch1(1,ki,jc)
               cc1(2,ki,l) = cc1(2,ki,l)+war*ch1(2,ki,j)
               cc1(2,ki,lc) = cc1(2,ki,lc)+wai*ch1(2,ki,jc)
  114       continue
  115    continue
  116 continue

      if ( 1 < ido ) go to 136
      sn = 1.0D+00 / real ( ip * l1, kind = 8 )
      if (na == 1) go to 146
      do 149 ki=1,lid
         cc1(1,ki,1) = sn*cc1(1,ki,1)
         cc1(2,ki,1) = sn*cc1(2,ki,1)
  149 continue
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ki=1,lid
            chold1 = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
            chold2 = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
            cc1(1,ki,j) = chold1
            cc1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
            cc1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
            cc1(1,ki,jc) = chold2
  119    continue
  120 continue
      return

  146 do 147 ki=1,lid
         ch1(1,ki,1) = sn*cc1(1,ki,1)
         ch1(2,ki,1) = sn*cc1(2,ki,1)
  147 continue
      do 145 j=2,ipph
         jc = ipp2-j
         do 144 ki=1,lid
            ch1(1,ki,j) = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
            ch1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
            ch1(1,ki,jc) = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
            ch1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
  144    continue
  145 continue
      return

  136 do 137 ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
  137 continue
      do 135 j=2,ipph
         jc = ipp2-j
         do 134 ki=1,lid
            ch1(1,ki,j) = cc1(1,ki,j)-cc1(2,ki,jc)
            ch1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
  134    continue
  135 continue
      do 131 i=1,ido
         do 130 k = 1, l1
            cc(1,k,1,i) = ch(1,k,i,1)
            cc(2,k,1,i) = ch(2,k,i,1)
  130    continue
  131 continue
      do 123 j=2,ip
         do 122 k = 1, l1
            cc(1,k,j,1) = ch(1,k,1,j)
            cc(2,k,j,1) = ch(2,k,1,j)
  122    continue
  123 continue
      do 126 j=2,ip
         do 125 i = 2, ido
            do 124 k = 1, l1
               cc(1,k,j,i) = wa(i,j-1,1)*ch(1,k,i,j) &
                            +wa(i,j-1,2)*ch(2,k,i,j)
               cc(2,k,j,i) = wa(i,j-1,1)*ch(2,k,i,j) &
                            -wa(i,j-1,2)*ch(1,k,i,j)
  124       continue
  125    continue
  126 continue

  return
end

!==================================================================================
subroutine c1fm1b ( n, inc, c, ch, wa, fnf, fac )

!*****************************************************************************80
!
!! C1FM1B is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  complex ( kind = 8 ) c(*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) inc2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) lid
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)

      inc2 = inc+inc
      nf = fnf
      na = 0
      l1 = 1
      iw = 1

      do k1=1,nf
         ip = fac(k1)
         l2 = ip*l1
         ido = n/l2
         lid = l1*ido
         nbr = 1+na+2*min(ip-2,4)
         go to (52,62,53,63,54,64,55,65,56,66),nbr
   52    call c1f2kb (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   62    call c1f2kb (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   53    call c1f3kb (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   63    call c1f3kb (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   54    call c1f4kb (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   64    call c1f4kb (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   55    call c1f5kb (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   65    call c1f5kb (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   56    call c1fgkb (ido,ip,l1,lid,na,c,c,inc2,ch,ch,2,wa(iw))
         go to 120
   66    call c1fgkb (ido,ip,l1,lid,na,ch,ch,2,c,c,inc2,wa(iw))
  120    l1 = l2
         iw = iw+(ip-1)*(ido+ido)
         if(ip <= 5) na = 1-na
      end do

  return
end

!==================================================================================
subroutine c1fm1f ( n, inc, c, ch, wa, fnf, fac )

!*****************************************************************************80
!
!! C1FM1F is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  complex ( kind = 8 ) c(*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) inc2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) lid
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)

      inc2 = inc+inc
      nf = fnf
      na = 0
      l1 = 1
      iw = 1

      do k1=1,nf

         ip = fac(k1)
         l2 = ip*l1
         ido = n/l2
         lid = l1*ido
         nbr = 1+na+2*min(ip-2,4)
         go to (52,62,53,63,54,64,55,65,56,66),nbr
   52    call c1f2kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   62    call c1f2kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   53    call c1f3kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   63    call c1f3kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   54    call c1f4kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   64    call c1f4kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   55    call c1f5kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   65    call c1f5kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   56    call c1fgkf (ido,ip,l1,lid,na,c,c,inc2,ch,ch,2,wa(iw))
         go to 120
   66    call c1fgkf (ido,ip,l1,lid,na,ch,ch,2,c,c,inc2,wa(iw))
  120    l1 = l2
         iw = iw+(ip-1)*(ido+ido)
         if(ip <= 5) na = 1-na
      end do

  return
end

!==================================================================================
subroutine cfft1b ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! CFFT1B: complex real(kind = 8) backward fast Fourier transform, 1D.
!
!  Discussion:
!
!    CFFT1B computes the one-dimensional Fourier transform of a single
!    periodic sequence within a complex array.  This transform is referred
!    to as the backward transform or Fourier synthesis, transforming the
!    sequence from spectral to physical space.
!
!    This transform is normalized since a call to CFFT1B followed
!    by a call to CFFT1F (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be
!    transformed.  The transform is most efficient when N is a product of
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in
!    array C, of two consecutive elements within the sequence to be transformed.
!
!    Input/output, complex ( kind = 8 ) C(LENC) containing the sequence to be
!    transformed.
!
!    Input, integer ( kind = 4 ) LENC, the dimension of the C array.
!    LENC must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be
!    initialized with a call to CFFT1I before the first call to routine CFFT1F
!    or CFFT1B for a given transform length N.  WSAVE's contents may be
!    re-used for subsequent calls to CFFT1F and CFFT1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.
!    LENWRK must be at least 2*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENC not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  complex ( kind = 8 ) c(lenc)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lenc < inc * ( n - 1 ) + 1 ) then
    ier = 1
    call xerfft ( 'cfft1b ', 4 )
  else if ( lensav < 2 * n + int ( log ( real ( n, kind = 8 ) ) &
    / log ( 2.0D+00 ) ) + 4 ) then
    ier = 2
    call xerfft ( 'cfft1b ', 6 )
  else if ( lenwrk < 2 * n ) then
    ier = 3
    call xerfft ( 'cfft1b ', 8 )
  end if

  if ( n == 1 ) then
    return
  end if

  iw1 = n + n + 1

  call c1fm1b ( n, inc, c, work, wsave, wsave(iw1), wsave(iw1+1) )

  return
end

!==================================================================================
subroutine cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! CFFT1F: complex real(kind = 8) forward fast Fourier transform, 1D.
!
!  Discussion:
!
!    CFFT1F computes the one-dimensional Fourier transform of a single
!    periodic sequence within a complex array.  This transform is referred
!    to as the forward transform or Fourier analysis, transforming the
!    sequence from physical to spectral space.
!
!    This transform is normalized since a call to CFFT1F followed
!    by a call to CFFT1B (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be
!    transformed.  The transform is most efficient when N is a product of
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in
!    array C, of two consecutive elements within the sequence to be transformed.
!
!    Input/output, complex ( kind = 8 ) C(LENC) containing the sequence to
!    be transformed.
!
!    Input, integer ( kind = 4 ) LENC, the dimension of the C array.
!    LENC must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be
!    initialized with a call to CFFT1I before the first call to routine CFFT1F
!    or CFFT1B for a given transform length N.  WSAVE's contents may be re-used
!    for subsequent calls to CFFT1F and CFFT1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.
!    LENWRK must be at least 2*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENC   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  complex ( kind = 8 ) c(lenc)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lenc < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cfft1f ', 4)
  else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0D+00 )) + 4) then
    ier = 2
    call xerfft ('cfft1f ', 6)
  else if (lenwrk < 2*n) then
    ier = 3
    call xerfft ('cfft1f ', 8)
  end if

  if (n == 1) then
    return
  end if

  iw1 = n + n + 1
  call c1fm1f ( n, inc, c, work, wsave, wsave(iw1), wsave(iw1+1) )

  return
end

!==================================================================================
subroutine cfft1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! CFFT1I: initialization for CFFT1B and CFFT1F.
!
!  Discussion:
!
!    CFFT1I initializes array WSAVE for use in its companion routines
!    CFFT1B and CFFT1F.  Routine CFFT1I must be called before the first
!    call to CFFT1B or CFFT1F, and after whenever the value of integer
!    N changes.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be
!    transformed.  The transform is most efficient when N is a product
!    of small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors
!    of N and  also containing certain trigonometric values which will be used
!    in routines CFFT1B or CFFT1F.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough.

  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) n
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < 2*n + int(log( real ( n, kind = 8 ) )/log( 2.0D+00 )) + 4) then
    ier = 2
    call xerfft ('cfftmi ', 3)
  end if

  if ( n == 1 ) then
    return
  end if

  iw1 = n+n+1

  call r8_mcfti1 ( n, wsave, wsave(iw1), wsave(iw1+1) )

  return
end

!==================================================================================
subroutine r8_factor ( n, nf, fac )

!*****************************************************************************80
!
!! R8_FACTOR factors of an integer for real real(kind = 8) computations.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number for which factorization and
!    other information is needed.
!
!    Output, integer ( kind = 4 ) NF, the number of factors.
!
!    Output, real ( kind = 8 ) FAC(*), a list of factors of N.
!
  implicit none

  real ( kind = 8 ) fac(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry

  nl = n
  nf = 0
  j = 0

  do while ( 1 < nl )

    j = j + 1

    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if

    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nf = nf + 1
      fac(nf) = real ( ntry, kind = 8 )
      nl = nq

    end do

  end do

  return
end

!==================================================================================
subroutine r8_mcfti1 ( n, wa, fnf, fac )

!*****************************************************************************80
!
!! R8_MCFTI1 sets up factors and tables, real real(kind = 8) arithmetic.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)
!
!  Get the factorization of N.
!
  call r8_factor ( n, nf, fac )
  fnf = real ( nf, kind = 8 )
  iw = 1
  l1 = 1
!
!  Set up the trigonometric tables.
!
  do k1 = 1, nf
    ip = int ( fac(k1) )
    l2 = l1 * ip
    ido = n / l2
    call r8_tables ( ido, ip, wa(iw) )
    iw = iw + ( ip - 1 ) * ( ido + ido )
    l1 = l2
  end do

  return
end

!==================================================================================
subroutine r8_tables ( ido, ip, wa )

!*****************************************************************************80
!
!! R8_TABLES computes trigonometric tables, real real(kind = 8) arithmetic.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip

  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) argz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido,ip-1,2)

  tpi = 8.0D+00 * atan ( 1.0D+00 )
  argz = tpi / real ( ip, kind = 8 )
  arg1 = tpi / real ( ido * ip, kind = 8 )

  do j = 2, ip

    arg2 = real ( j - 1, kind = 8 ) * arg1

    do i = 1, ido
      arg3 = real ( i - 1, kind = 8 ) * arg2
      wa(i,j-1,1) = cos ( arg3 )
      wa(i,j-1,2) = sin ( arg3 )
    end do

    if ( 5 < ip ) then
      arg4 = real ( j - 1, kind = 8 ) * argz
      wa(1,j-1,1) = cos ( arg4 )
      wa(1,j-1,2) = sin ( arg4 )
    end if

  end do

  return
end

!==================================================================================
subroutine xerfft ( srname, info )

!*****************************************************************************80
!
!! XERFFT is an error handler for the FFTPACK routines.
!
!  Discussion:
!
!    XERFFT is an error handler for FFTPACK version 5.1 routines.
!    It is called by an FFTPACK 5.1 routine if an input parameter has an
!    invalid value.  A message is printed and execution stops.
!
!    Installers may consider modifying the stop statement in order to
!    call system-specific exception-handling facilities.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, character ( len = * ) SRNAME, the name of the calling routine.
!
!    Input, integer ( kind = 4 ) INFO, an error code.  When a single invalid
!    parameter in the parameter list of the calling routine has been detected,
!    INFO is the position of that parameter.  In the case when an illegal
!    combination of LOT, JUMP, N, and INC has been detected, the calling
!    subprogram calls XERFFT with INFO = -1.
!
  implicit none

  integer ( kind = 4 ) info
  character ( len = * ) srname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XERFFT - Fatal error!'

  if ( 1 <= info ) then
    write ( *, '(a,a,a,i3,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter number ', info, ' had an illegal value.'
  else if ( info == -1 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameters LOT, JUMP, N and INC are inconsistent.'
  else if ( info == -2 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter L is greater than LDIM.'
  else if ( info == -3 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter M is greater than MDIM.'
  else if ( info == -5 ) then
    write ( *, '(a,a,a,a)' ) '  Within ', trim ( srname ), &
      ' input error returned by lower level routine.'
  else if ( info == -6 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter LDIM is less than 2*(L/2+1).'
  end if

  stop
end
