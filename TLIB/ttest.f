c      
      subroutine testscale (A,m)
       double precision, intent (in):: A(9)
       double precision::scaled(9)
       double precision:: m
c
       print *,"----------------"
       print *,"Begin Scale Test"
       print *,"----------------"
       print *,"Input Tensor:"
       print *,
       call t2print (A)
       print *,
       print *,"Scale Factor:"
       print "(1e12.4)",m
       print *,
       print *,"Output Tensor:"
       call scale2(A,scaled,m)
       call t2print(scaled)
       print *,
       print *,"----------------"
       print *,"End Scale Test"
       print *,"----------------"
c
      end subroutine testscale
c
      subroutine testtpose (A)
       double precision, intent (in)::A(9)
       double precision::posed(9)
c
       print *,"----------------"
       print *,"Begin TPOSE Test"
       print *,"----------------"
       print *,"Input Tensor:"
       print *,
       call t2print (A)
       print *,
       print *,"Output Tensor:"
       call tpose(A,posed)
       call t2print(posed)
       print *,
       print *,"----------------"
       print *,"End TPOSE Test"
       print *,"----------------"
c
      end subroutine testtpose
c
      subroutine testdet (A)
       double precision, intent(in)::A(9)
       double precision:: det_a
c
       print *,"----------------"
       print *,"Begin DET Test"
       print *,"----------------"
       print *,"Input Tensor:"
       print *,
       call t2print (A)
       print *,
       print *,"Determinant:"
       call det(A,det_A)
       print "(1e12.4)",det_A
       print *,
       print *,"----------------"
       print *,"End DET Test"
       print *,"----------------"
c
      end subroutine testdet
c
      subroutine testinv (A)
c
       double precision, intent(in):: A(9)
       double precision::inverse(9)
c
       print *,"----------------"
       print *,"Begin Invert Test"
       print *,"----------------"
       print *,"Input Tensor:"
       print *,
       call t2print (A)
       print *,
       print *,"Inverse Tensor:"
       call invert(A,inverse)
       call t2print(inverse)
       print *,
       print *,"----------------"
       print *,"End Invert Test"
       print *,"----------------"
c
      end subroutine testinv
c
      subroutine testtc_2s2 (A,B)
c
       double precision, intent(in):: A(9),B(9)
       double precision::C(9)
c
       print *,"----------------"
       print *,"Begin tc_2s2 Test"
       print *,"----------------"
       print *,"Input Tensor A:"
       print *,
       call t2print (A)
       print *,
       print *,"Input Tensor B:"
       print *,
       call t2print (B)
       print *,
       print *,"Output Tensor:"
       call tc_2s2(A,B,C)
       call t2print(C)
       print *,
       print *,"----------------"
       print *,"End tc_2s2 Test"
       print *,"----------------"
c
      end subroutine testtc_2s2
c
      subroutine testtc_2d2 (A,B)
c
       double precision, intent(in):: A(9),B(9)
       double precision::c
c
       print *,"----------------"
       print *,"Begin tc_2d2 Test"
       print *,"----------------"
       print *,"Input Tensor A:"
       print *,
       call t2print (A)
       print *,
       print *,"Input Tensor B:"
       print *,
       call t2print (B)
       print *,
       print *,"Output:"
       call tc_2d2(A,B,c)
       print "(1e12.4)",c
       print *,
       print *,"----------------"
       print *,"End tc_2d2 Test"
       print *,"----------------"
c
      end subroutine testtc_2d2
c
      subroutine testtc_4d4 (A,B)
c
       double precision, intent(in):: A(9,9),B(9,9)
       double precision:: C(9,9)
c
       print *,"----------------"
       print *,"Begin tc_4d4 Test"
       print *,"----------------"
       print *,"Input Tensor A:"
       print *,
       call t4print (A)
       print *,
       print *,"Input Tensor B:"
       print *,
       call t4print (B)
       print *,
       print *,"Output:"
       call tc_4d4(A,B,C)
       print *,
       call t4print (C)
       print *,
       print *,"----------------"
       print *,"End tc_4d4 Test"
       print *,"----------------"
c
      end subroutine testtc_4d4
c
      subroutine testinvariants (A)
c
       double precision :: A(9)
       double precision :: x
c
       print *,"----------------"
       print *,"Begin I,J test"
       print *,"----------------"
c
       print *, "Input tensor:"
       print *,
       call t2print(A)
       print *,
       print *,"I1:"
       call I1(A,x)
       print "(1e12.4)",x
       print *,
c
       print *,
       print *,"I2:"
       call I2(A,x)
       print "(1e12.4)",x
       print *,
c
       print *,
       print *,"I3:"
       call I3(A,x)
       print "(1e12.4)",x
       print *,
c
       print *,
       print *,"J1:"
       call J1(A,x)
       print "(1e12.4)",x
       print *,
c
       print *,
       print *,"J2:"
       call J2(A,x)
       print "(1e12.4)",x
       print *,
c
       print *,
       print *,"J3:"
       call J3(A,x)
       print "(1e12.4)",x
       print *,
c
       print *,
       print *,"----------------"
       print *,"End I,J test"
       print *,"----------------"
       print *,
c
      end subroutine testinvariants
c
      subroutine testzrotation (A,de)
c
       double precision, intent (in) :: A(9)
       double precision :: B(9),Q(9)
       double precision :: de,th,pi
c
c     I don't think there's a pi constant
c     in FORTRAN...
c      
       pi = 4*atan(1.0_8)
       th = pi*de/180.0_8
c
       print *,"----------------"
       print *,"Start rot test"
       print *,"----------------"
       print *,
       print *,"Input Tensor: "
       print *, 
       call t2print (A)
       print *,
       print *,"Input Angle:"
       print "(1e12.4)",de
       print *,
       print *,"Rotated Tensor:"
       print *,
       call qzrad (th,Q)
       call transform2 (A,Q,B)
       call t2print (B)
       print *,
       print *,"----------------"
       print *,"End rot test"
       print *,"----------------"
       return       
c
      end subroutine testzrotation
c
      subroutine testIDN (A)
c
       double precision :: A(9),I4(9,9),B(9)
c
       print *,
       print *,"----------------"
       print *,"Start IDN test"
       print *,"----------------"
       print *,
       print *,"Input Tensor:"
       print *,
       call t2print(A)
       print *,
       print *,"Output Tensor (IDN4):"
       print *,
       call IDN4 (I4)
       call tc_4d2(I4,A,B)
       call t2print(B)
       print *,
       print *,"Output Tensor (SYM4):"
       print *,
       call SYM4 (I4)
       call tc_4d2(I4,A,B)
       call t2print(B)
       print *,
       print *,"Output Tensor (ANT4):"
       print *,
       call ANT4 (I4)
       call tc_4d2(I4,A,B)
       call t2print(B)
       print *,
       print *,"Output Tensor (VOL4):"
       print *,
       call VOL4 (I4)
       call tc_4d2(I4,A,B)
       call t2print(B)
       print *,
       print *,"Output Tensor (DEV4):"
       print *,
       call DEV4 (I4)
       call tc_4d2(I4,A,B)
       call t2print(B)
       print *,
       print *,"----------------"
       print *,"End IDN test"
       print *,"----------------"
c
       return
c
      end subroutine testIDN
c
      subroutine testvonMises (A)
c
       double precision, intent (in) :: A(9)
       double precision :: vm
c
       print *,
       print *,"----------------"
       print *,"Start VM test"
       print *,"----------------"
       print *,
       print *,"Input Tensor:"
       print *,
       call t2print(A)
       print *,
c       call SEQ(A,vm)
       print *,
       print *,"von Mises:"
       print *,
c       print "(1e12.4)",vm
       call VMS(A,vm)
       print *,vm
       print *,      
       print *,"----------------"
       print *,"End VM test"
       print *,"----------------"
c
       return
c
      end subroutine testvonMises
c      
c
        subroutine testCH (A)
c
         double precision, intent (in) :: A(9)
         double precision :: B(9)
         double precision :: I1A,I2A,I3A
         double precision :: I1B,I2B,I3B
c
         call CH (A,B)
         call I1 (A,I1A)
         call I2 (A,I2A)
         call I3 (A,I3A)
         call I1 (B,I1B)
         call I2 (B,I2B)
         call I3 (B,I3B)
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
         print *,"1st Invariant (should be equal)"
         print *,
         print "(2e12.4)",I1A,I1B
         print *,
         print *,
         print *,"2nd Invariant (should be equal)"
         print *,
         print "(2e12.4)",I2A,I2B
         print *,
         print *,
         print *,"3rd Invariant (should be equal)"
         print *,
         print "(2e12.4)",I3A,I3B
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
