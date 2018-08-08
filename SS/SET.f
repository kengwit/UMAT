c
c  Single Element Test-bed
c  Prusodman Sathananthan
c  August 2018
c
       program SET
c
         implicit none
c        
c ## VARIABLE DEFINITIONS, with descriptions ##
c  
c  1. control hypo-ep integration
c  2. deformation gradients at previous and current timestep
c  3. kinematics of deformation derived from F0, F1 
c  4. stress, strain and strain increment
c  5. spin tensor and required tensor products for objective update
c
c  #1
         double precision :: dalp, alp, alpf, b
         integer i, inc
c  #2
         double precision :: F0 (9), F1(9)
c  #3
         double precision :: F(9),dF(9),l(9),D(9),W(9)
c  #4    
         double precision :: ep(9), deps(9), sig (9), hisv(11)
c  #5
         double precision :: R(9), eR(9), Re(9), sigR(9), Rsig(9)
c
c  Initialize all variables for stress/strain integration
c
         dalp = 0.001_8
         alp = 0.0_8
         alpf = 1.0_8
         b = 0.5_8
         inc = alpf/dalp      
c
         ep = 0.0_8
         deps = 0.0_8
         sig = 0.0_8
         F1 = (/1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 1.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 1.0_8 /)
c
c  Step through deformation 
c
        do i=1,inc
c
c  Get kinematics of deformation
c
           F0 = F1
           call kindefgen(F1,alp,1)
           call kindef (F0,F1,b,dF,F,l,D,W)
c
c  Apply objective update to stress and strain
c
           R = W
           call tc_2s2 (ep,R,eR)
           call tc_2s2 (R,ep,Re)
c
           deps = D
           ep = ep + deps - eR + Re
c
           call tc_2s2 (sig,R,sigR)
           call tc_2s2 (R,sig,Rsig)
           sig = sig - sigR + Rsig
c
c          call umat here
c

c
c Output results
c
           print "(21e15.6)", alp,ep,sig,hisv(1),hisv(11)
c
           alp = alp + dalp
c
        end do
c
c
       end program SET
