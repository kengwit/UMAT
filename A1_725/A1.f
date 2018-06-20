      
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
      program A1_725
c
       double precision:: A(9),B(9),C(9),test_mat(9),
     &                    transform_sij(9),SA(9),SB(9),
     &                    e1(9),s1(9)
       double precision:: A4(9,9),B4(9,9),I4(9,9),L4(9,9),
     &                    C4(9,9),R4(9,9)
c
       do i=1,9
        A(i)=i
       end do
c
       do i=1,9
         do j=1,9
           A4(i,j) = i+j
           B4(i,j) = i*j
         end do
       end do
c 
       test_mat = (/ 4,1,6,9,0,5,2,8,9 /)
c
       transform_sij = (/  1.0,  0.1,  0.0,
     &                     0.1,  0.5,  0.0,
     &                     0.0,  0.0, 0.01  /)
c
       SA = (/  0.1,  1.1,  0.0,
     &          1.1, -0.1,  0.0,
     &          0.0,  0.0,  0.0  /)
c
       SB = (/  1.0,  0.1, -0.2,
     &          0.1,  0.5,  0.4,
     &         -0.2,  0.4, 0.05  /)
c
       e1 = (/  0.001_8 ,  0.0001_8,  0.0_8   ,
     &          0.0001_8, -0.0005_8,  0.0_8   ,
     &          0.0_8   ,  0.0_8   , -0.0005_8  /)
c
c
c        call testscale (A,5.0_8)
c        call testscale (test_mat,1e-20_8)
c        call testtpose (A)
c        call testtpose (test_mat)
c        call testdet (A)
c        call testdet (test_mat)
c        call testinv (A)
c        call testinv (test_mat)
c        call testtc_2s2(A,test_mat)
c         call testtc_2d2(A,test_mat)
c        call testtc_4d4(A4,B4)
c        call testinvariants(A)
c        call testzrotation(transform_sij,60.0_8)
c        call testIDN(A)
c        call DEV4(I4)
c        call t4print(I4)
        call testvonMises (SA)
        call testvonMises (SB)
c        call L4_EL(205e9_8,0.28_8,L4)
c        call tc_4d2(L4,e1,s1)
c        call t2print(s1)
c        call C4_EL(205e9_8,0.28_8,C4)
c        call tc_4d2(C4,s1,e1)
c        call t2print (e1)
c        call tc_4d4(C4,L4,R4)
c
c        call t4print(C4)  
c        call t4print(R4)    
c
      end program A1_725
