c
        program A3_725
c
         implicit none
c
         double precision ::  DP980 (3)
         double precision ::  ep (4), Sb, h
         double precision ::  ep_curr
         double precision ::  A(9),seq,N(9),ph(9,9) 
         double precision ::  N1(9),N2(9),Ntest,ans
         integer i
c
c         DP980 = (/ 1300.0_8, 0.002_8, 0.07_8 /)
c         ep    = (/ 0.00_8, 0.05_8, 0.20_8, 0.50_8 /)
c         print *,
c         print *,"         ep,        Sb,          h"
c         print *,"=================================="
c
c         ep_curr = 0.01_8
c         do i=1,100
c           call iso_hard (DP980(1),DP980(2),DP980(3),ep_curr,Sb,h)
c           print "(3e12.4)",ep_curr,Sb,h
c           ep_curr = ep_curr + 0.01_8
c         end do
c
c         print *,
c         print *,
c
c         do i=1,4
c           call iso_hard (DP980(1),DP980(2),DP980(3),ep(i),Sb,h)
c           print "(3e12.4)",ep(i),Sb,h
c         end do
c    
c         print *,
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8 /)
         call testvmbox (A)
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 0.5_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8 /)
         call testvmbox (A)
c
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 1.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8 /)
         call testvmbox (A)
c
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, -1._8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8 /)
         call testvmbox (A)
c
c
         A = (/ 0.0_8, 1.0_8, 0.0_8,
     &          1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8 /)
         call testvmbox (A)
c
c
         A = (/ 1.0_8, 0.1_8, -0.2_8,
     &          0.1_8, 0.5_8, 0.3_8,
     &          -0.2_8, 0.3_8, 0.25_8 /)
         call testvmbox (A)
c
        end program A3_725
