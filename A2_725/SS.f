        program simpleshear
c
          implicit none
c
          double precision :: dgam, gam, gf, inc, b
          double precision :: E, v
          double precision :: F0(9), F1(9) 
          double precision :: F(9),dF(9),l(9),D(9),W(9)
          double precision :: e(9), deps(9), eW(9), We(9)
          double precision :: sig (9)
          integer i
c
          dgam = 0.001_8
          gam = 0.0_8
          gf = 0.5_8
c
          inc = gf/dgam
          E = 205e9_8
          v = 0.25_8
c
          e = 0.0_8
          sig = 0.0_8
          F1 = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /) 
c          
          do i=1,inc
c
c           Get kinematics of deformation
c 
            F0 = F1
            F1 = (/1.0, gam, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
            call kindef (F0,F1,b, dF,F,l,D,W)
c
c           Perform objective update on strain
c
            call tc_2s2 (e,W,eW)
            call tc_2s2 (W,e,We)
            deps = D - eW + We
            e = e + deps
c
            
c            
            gam = gam + dgam
          end do
c                    
        end program simpleshear
