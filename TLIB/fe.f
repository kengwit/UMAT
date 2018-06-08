c
c Implemented to date
c ===================
c  I1, I2, I3  -> tensor invariants
c  A2S         -> obtain deviatoric component of tensor
c  J1, J2, J3  -> deviatoric invariants
c  L4_EL       -> determine elastic tensor using E,v
c  C4_EL       -> determine compilance tensor using E,v
c  KD          -> returns Kronecker's delta
c  IDN4        -> returns 4th order identity tensor  
c  SYM4	       -> returns symmetric identity tensor
c  ANT4        -> returns antisymmetric identity tensor
c  VOL4        -> returns volumetric identity tensor
c  DEV4        -> returns deviatoric identity tensor
c  SEQ         -> find von Mises stress of a tensor
c
       subroutine I1 (A2,x)
c
        implicit none
        double precision, intent (in) :: A2(9)
        double precision, intent (out) :: x
c
        x = A2(1) + A2(5) + A2(9)
        return
c 
       end subroutine I1
c
c
c
       subroutine I2 (A2,x)
c
        implicit none
        double precision, intent (in) :: A2(9)
        double precision, intent (out) :: x
        double precision :: trA, trASQ, ASQ(9)
c 
        call tc_2s2(A2,A2,ASQ)
        call I1(A2,trA)
        call I1(ASQ,trASQ)
c
        x = 0.5*trA*trA-0.5*trASQ
        return 
c
       end subroutine I2
c
c
       subroutine I3 (A2,x)
c
        implicit none
        double precision, intent (in) :: A2(9)
        double precision, intent (out) :: x
c
        call det(A2,x)
        return
c       
       end subroutine I3
c
c
c
       subroutine A2S (A2,S2)
c
        implicit none
        double precision, intent (in) :: A2(9)
        double precision, intent (out) :: S2(9)
        double precision :: trA, P
c
        call I1(A2,trA)
        P = trA / 3.0_8
        S2(1) = A2(1) - P
        S2(2) = A2(2)
        S2(3) = A2(3)        
        S2(4) = A2(4)
        S2(5) = A2(5) - P
        S2(6) = A2(6)        
        S2(7) = A2(7)
        S2(8) = A2(8)
        S2(9) = A2(9) - P        
        return
c
       end subroutine A2S
c
c
c      
       subroutine J1 (A2,x)
c
        implicit none
        double precision, intent (in) :: A2(9)
        double precision, intent (out) :: x
        double precision :: S2(9)
c
        call A2S (A2,S2)
        call I1 (S2,x)
        return
c
       end subroutine J1
c
c
c
       subroutine J2 (A2,x)
c
        implicit none
        double precision, intent (in) :: A2(9)
        double precision, intent (out) :: x
        double precision :: S2(9)
c
        call A2S (A2,S2)
        call I2 (S2,x)
        return
c
       end subroutine J2
c
c
c
       subroutine J3 (A2,x)
c
        implicit none
        double precision, intent (in) :: A2(9)
        double precision, intent (out) :: x
        double precision :: S2(9)
c
        call A2S (A2,S2)
        call I3 (S2,x)
        return
c
       end subroutine J3
c
c
c
       subroutine L4_EL (E,v,L4)
c
        implicit none
        double precision, intent (in) :: E,v
        double precision, intent (out) :: L4(9,9)
        double precision :: G,K,LAME
        double precision :: eye(9)
        double precision :: B4(9,9)
c
        G = 0.5_8*E/(1+v)
        K = E/3.0_8/(1-2*v)
        LAME = K-2.0_8*G/3.0_8
c
c    Lel = LAME*(IXI)+2*G*I
c
        call KD(eye)
        call tp_2dy2(eye,eye,L4)
        call scale4 (L4,L4,LAME)
        call SYM4(B4)
        call scale4(B4,B4,2*G)
        L4 = L4 + B4
c
        return
c
       end subroutine L4_EL
c
       subroutine C4_EL (E,v,C4)
c
        implicit none
        double precision, intent (in) :: E,v
        double precision, intent (out) :: C4(9,9)
        double precision :: G,K,LAME,alp
        double precision :: eye(9)
        double precision :: B4(9,9)
c
        G = 0.5_8*E/(1+v)
        K = E/3.0_8/(1-2*v)
        LAME = K-2.0_8*G/3.0_8
        alp = -LAME/2/G/(3*LAME+2*G)
c
        call KD(eye)
        call tp_2dy2(eye,eye,C4)
        call scale4 (C4,C4,alp)
        call SYM4(B4)
        call scale4(B4,B4,0.5_8/G)
        C4 = C4 + B4
c
        return
c
       end subroutine C4_EL
c
       subroutine KD (del)
c
        implicit none
        double precision, intent (out) :: del (9)
c
        del = (/ 1.0, 0.0, 0.0,
     &           0.0, 1.0, 0.0,     
     &           0.0, 0.0, 1.0  /)
c
        return
c
       end subroutine KD
c
       subroutine IDN4 (I)
c
        implicit none
        double precision, intent (out) :: I (9,9)
c
        I = 0.0_8
        I(1,1) = 1.0_8
        I(2,2) = 1.0_8
        I(3,3) = 1.0_8
        I(4,4) = 1.0_8
        I(5,5) = 1.0_8
        I(6,6) = 1.0_8
        I(7,7) = 1.0_8
        I(8,8) = 1.0_8
        I(9,9) = 1.0_8
c        
        return
c
       end subroutine IDN4
c
       subroutine SYM4 (I)
c
        implicit none
        double precision, intent (out) :: I(9,9)
        integer x,ji(9)       
c
        ji = (/ 1,4,7,2,5,8,3,6,9 /)
c
        I = 0.0_8
c
        do x=1,9
         I(x,x) = 0.5_8
         I(x,ji(x)) = 0.5_8
        end do
c
        I(1,1)= 1.0_8
        I(5,5)= 1.0_8
        I(9,9)= 1.0_8
c
        return
c
       end subroutine SYM4
c
       subroutine ANT4 (I)
c
        implicit none
        double precision, intent (out) :: I(9,9)
        integer x,ji(9)       
c
        ji = (/ 1,4,7,2,5,8,3,6,9 /)
c
        I = 0.0_8
c
        do x=1,9
         I(x,x) = 0.5_8
         I(x,ji(x)) = -0.5_8
        end do
c
        I(1,1) = 0.0_8
        I(5,5) = 0.0_8
        I(9,9) = 0.0_8
c
        return
c
       end subroutine ANT4
c
       subroutine VOL4 (I)
c
        implicit none
        double precision, intent (out) :: I(9,9)
c
        I = 0.0_8
        I(1,1) = 1.0_8/3
        I(1,5) = 1.0_8/3
        I(1,9) = 1.0_8/3
        I(5,1) = 1.0_8/3
        I(5,5) = 1.0_8/3
        I(5,9) = 1.0_8/3
        I(9,1) = 1.0_8/3
        I(9,5) = 1.0_8/3
        I(9,9) = 1.0_8/3
c
        return
c
       end subroutine VOL4
c
       subroutine DEV4 (I)
c
        implicit none
        double precision, intent (out) :: I(9,9)
        double precision :: eye4(9,9)
c   
        call VOL4(I)
        call IDN4(eye4)
        I = eye4-I
c
        return
c
       end subroutine DEV4
c
       subroutine SEQ (A,vm)
c
        implicit none
        double precision, intent (in) :: A(9)
        double precision, intent (out) :: vm
        double precision :: DEV(9,9), J2
c
c       
        call DEV4(DEV)
        call tc_4d2(DEV,A,A)
        call tc_2d2 (A,A,vm)
        vm = SQRT(1.5_8*vm)
c
        return
c
       end subroutine SEQ
c
