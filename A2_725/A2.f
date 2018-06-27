        subroutine testCH (A)
c
         double precision, intent (in) :: A(9)
         double precision :: B(9)
c
         call CH (A,B)
         print *,"----------------"
         print *,"Begin CH test"
         print *,"----------------"
         print *,"Input Tensor:"
         print *,
         call t2print (A)
         print *,
         print *,"Principal Values:"
         print *,
         call t2print (B)
         print *,
         print *,"----------------"
         print *,"End CH test"
         print *,"----------------"  
c
        end subroutine testCH
c
        subroutine testPD (F)
c
         double precision, intent (in) :: F(9)
         double precision :: R(9),U(9),V(9),C1(9),C2(9)
         call PD (F,R,U,V)
         call tc_2s2(R,U,C1)
         C1 = F-C1
         call tc_2s2(V,R,C2)
         C2 = F-C2
         print *,"----------------"
         print *,"Begin DECOMP test"
         print *,"----------------"
         print *,"Input Tensor:"
         print *,
         call t2print (F)
         print *,
         print *,"R:"
         print *,
         call t2print (R)
         print *,
         print *,"U:"
         print *,
         call t2print (U)
         print *,
         print *,"V:"
         print *,
         call t2print (V)
         print *,
         print *,"C1 (should be zero):"
         print *,
         call t2print (C1)
         print *,
         print *,"C2 (should be zero):"
         print *,
         call t2print (C2)
         print *,"----------------"
         print *,"End DECOMP test"
         print *,"----------------"  
c
        end subroutine testPD
c
        program A2_725
c
         double precision :: A(9)
         double precision :: f
         f = 1.0e-8_8
c
c        ZERO CHECK
c
         A = 0.0_8
         call testCH(A)
c
c        ALREADY IN PRINCIPAL SPACE
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 1.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 1.0_8 /)
         call testCH(A)
c
c       ONE ZERO PRINCIPAL VALUE
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 1.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8 /)
         call testCH(A)
c
c       TWO ZERO PRINCIPAL VALUES
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8 /)
         call testCH(A)
c
c       TOLERANCE BUFFET
c
         A = (/ 1.0_8+f,           0.0_8, 0.0_8,
     &            0.0_8, 1.0_8/(1.0_8+f), 0.0_8,
     &            0.0_8,           0.0_8,     f /)
         call testCH(A)
c
c       ALL TOLERANCE
c
         A = (/    f, 0.0_8, 0.0_8,
     &         0.0_8,     f, 0.0_8,
     &         0.0_8, 0.0_8,     f /)
         call testCH(A)
c
c       RANDOM VALUES
c
         A = (/ 300.0_8,  -20.0_8,   25.0_8,
     &          -20.0_8,  100.0_8, -125.0_8,
     &           25.0_8, -125.0_8,   50.0_8 /)
         call testCH(A)
c
        end program A2_725
