        program simpleshear
c
          implicit none
c
          double precision :: dgam, gam, gf, inc, b
          double precision :: E, v
          double precision :: F0(9), F1(9) 
          double precision :: F(9),dF(9),l(9),D(9),W(9)
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
          F1 = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /) 
c          
          do i=1,inc
            F0 = F1
            F1 = (/1.0, gam, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
            call kindef (F0,F1,b, dF,F,l,D,W)
            
            gam = gam + dgam
          end do
c                    
        end program simpleshear
