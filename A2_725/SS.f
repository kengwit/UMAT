        program simpleshear
c
          implicit none
c
          double precision :: dgam, gam, gf, b
          double precision :: E, v
          double precision :: F0(9), F1(9) 
          double precision :: F(9),dF(9),l(9),D(9),W(9)
          double precision :: ep(9), deps(9), eW(9), We(9)
          double precision :: sig (9), sigW(9), Wsig(9)
          double precision :: cm(2), hisv(1)
          integer i,inc
c
          dgam = 0.001_8
          gam = 0.0
          gf = 0.004_8
c
          inc = gf/dgam
          E = 205e9_8
          v = 0.25_8
c
          ep = 0.0_8
          sig = 0.0_8
          F1 = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /) 
c          
          do i=1,inc
c
c           Get kinematics of deformation
c 
            F0 = F1
            F1 = (/ 1.0_8,   gam, 0.0_8,
     &              0.0_8, 1.0_8, 0.0_8,
     &              0.0_8, 0.0_8, 1.0_8 /)
c
            call kindef (F0,F1,b, dF,F,l,D,W)
c
c           Perform objective update on strain
c
            call tc_2s2 (e,W,eW)
            call tc_2s2 (W,e,We)
            deps = D - eW + We
            ep = ep + deps
c
c           Perform objective update on stress (-sigW + Wsig, the objective update)
c
c            call tc_2s2 (sig,W,sigW)
c            call tc_2s2 (W,sig,Wsig)
c            sig = sig - sigW + Wsig           
c
c           Perform stress update on stress (sig dot hat, the stress update)
c
            call umat_elastic (cm,deps,sig,hisv)
c
            print *,           
            call t2print (sig)
            print *,
            gam = gam + dgam
c         
          end do
c                    
        end program simpleshear
