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
c
         A = (/ 4.0,1.0,6.0,1.0,0.0,5.0,6.0,5.0,9.0 /)
c         call testCH (A)
         call testPD (A) 
c
        end program A2_725
