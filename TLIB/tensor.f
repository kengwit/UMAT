c ===============================
c
c Tensor Library for ME725
c Version 0.1 LM: 20180605
c P.Sathananthan
c
c ===============================
c
c Subroutines implemented to date
c
c  tpose (A,B)       -> Aij = Bji           (implemented & verified)
c  det (A,x)         -> det(Aij) = x        (implemented & verified)
c  scale2 (A,B,m)    -> Bij = m.*Aij        (implemented & verified)
c  scale4 (A,B,m)    -> Bijkl = m.*Aijkl    (implemented)
c  invert (A,B)      -> B = inv (A)         (implemented & verified)
c  tc_2s2 (A,B,C)    -> Cij = Aik*Bkj       (implemented)
c  tc_2d2 (A,B,c)    -> c = Aij*Bji         (implemented)
c  tc_4d4 (A4,B4,C4) -> Cijkl = Aijmn*Bmnkl (implemented)
c  tc_4d2 (A4,B,C)   -> Cij = Aijkl*Bkl     (implemented)
c  tc_2d4 (A,B4,C)   -> Cij = Akl*Bklij     (implemented)
c  transform2        -> B = QAQt            (implemented)
c  qzrad             -> rot tensor about z  (implemented)
c  qxrad             -> rot tensor about x  (implemented)
c  qyrad             -> rot tensor about y  (implemented)
c  tp_2dy2           -> Cijkl = AijXBkl
c
c -------------------------------
c
c   2nd order tensor indexing (aij)
c
c   a11 a12 a13   (1) (2) (3) 
c   a21 a22 a23   (4) (5) (6)
c   a31 a32 a33   (7) (8) (9)
c
c    index (i,j) = 3*(i-1)+j 
c
c -------------------------------
c
c   4th order tensor indexing (aijkl)
c
c   a1111 a1112 ... a1132 a1133   (1,1) (1,2) ... (1,8) (1,9)
c     .                     .       .                  .
c     .                     .       .                  .
c     .                     .       .                  .
c   a3311 a3312 ... a3332 a3333   (9,1) (9,2) ... (9,8) (9,9)
c   
c
c   index1 (i,j) = 3*(i-1)+j
c   index2 (k,l) = 3*(k-1)+l
c 
c ===============================
c
c transpose 2nd order tensor subroutine
c
c input:  Aij (a 2nd order tensor)
c output: Bji (a 2nd order tensor Aij=Bji)
c
       subroutine tpose (A,B)
c
c       variable declaration region
c       declare intent to avoid accidental changes (C.Butcher)
c     
        implicit none
        double precision, intent (in) :: A(9)
        double precision, intent (out) :: B(9)
c
c       hard-code transpose for SPEED!!!
c
        B(1)=A(1)
	B(2)=A(4)
        B(3)=A(7)
        B(4)=A(2)
        B(5)=A(5)
        B(6)=A(8)
        B(7)=A(3)
        B(8)=A(6)
        B(9)=A(9)
c
        return
c
       end subroutine tpose
c
c determinant of a 2nd order tensor subroutine
c
c
       subroutine det (A,x)
c
        implicit none
        double precision, intent (in) :: A(9)
        double precision, intent (out) :: x
c
        x = A(1)*A(5)*A(9)+A(2)*A(6)*A(7)+A(3)*A(4)*A(8)
     &     -A(1)*A(6)*A(8)-A(2)*A(4)*A(9)-A(3)*A(5)*A(7)
c
        return
       end subroutine det
c
c scale 2nd order tensor subroutine
c
       subroutine scale2 (A,B,m)
c
        implicit none
        double precision, intent (in) :: A(9)
        double precision, intent (out) :: B(9)
        double precision, intent (in) :: m
c
        B(1)=m*A(1)
	B(2)=m*A(2)
        B(3)=m*A(3)
        B(4)=m*A(4)
        B(5)=m*A(5)
        B(6)=m*A(6)
        B(7)=m*A(7)
        B(8)=m*A(8)
        B(9)=m*A(9)
c
        return
       end subroutine scale2
c
c scale 4th order tensor subroutine
c
       subroutine scale4 (A4,B4,m)
c
        implicit none
        double precision, intent (in) :: A4(9,9)
        double precision, intent (out) :: B4(9,9)
        double precision, intent (in) :: m
        integer ij
c
        do ij=1,9
         B4(ij,1)=m*A4(ij,1)
	 B4(ij,2)=m*A4(ij,2)
         B4(ij,3)=m*A4(ij,3)
         B4(ij,4)=m*A4(ij,4)
         B4(ij,5)=m*A4(ij,5)
         B4(ij,6)=m*A4(ij,6)
         B4(ij,7)=m*A4(ij,7)
         B4(ij,8)=m*A4(ij,8)
         B4(ij,9)=m*A4(ij,9)
        end do
c
        return
       end subroutine scale4
c
c invert 2nd order tensor subroutine
c        using cofactor method
c
       subroutine invert (A,B)
c
       implicit none
       double precision, intent (in) :: A(9)
       double precision, intent (out) :: B(9)
       double precision minor (9)
       double precision det_A
c
c DETERMINE WHETHER MATRIX IS SINGULAR (NOT INVERTIBLE)
c
       call det (A,det_A)
       if (det_A.EQ.0.0_8) then
         B = (/ 0,0,0,0,0,0,0,0,0 /)
         return
       end if
c
c MATRIX OF MINORS (COFACTORS FROM HERE)
c
       minor (1) =  A(5)*A(9) - A(6)*A(8)
       minor (2) = -A(4)*A(9) + A(6)*A(7)
       minor (3) =  A(4)*A(8) - A(5)*A(7)
       minor (4) = -A(2)*A(9) + A(3)*A(8)
       minor (5) =  A(1)*A(9) - A(3)*A(7)
       minor (6) = -A(1)*A(8) + A(2)*A(7)
       minor (7) =  A(2)*A(6) - A(3)*A(5)
       minor (8) = -A(1)*A(6) + A(3)*A(4)
       minor (9) =  A(1)*A(5) - A(2)*A(4)      
c
c TRANSPOSE MAT TO GET ADJUGATE AND FIND INV
c
       call tpose (minor,B)
       call scale2 (B,B,1.0_8/det_A)
c
       end subroutine invert
c
c
c Single Contraction of 2 2nd order Tensors
c
       subroutine tc_2s2 (A,B,C)
c
        implicit none
        double precision, intent (in) :: A(9),B(9)
        double precision, intent (out) :: C(9)
c
        C(1)=A(1)*B(1)+A(2)*B(4)+A(3)*B(7)
        C(2)=A(1)*B(2)+A(2)*B(5)+A(3)*B(8)
        C(3)=A(1)*B(3)+A(2)*B(6)+A(3)*B(9) 
        C(4)=A(4)*B(1)+A(5)*B(4)+A(6)*B(7)
        C(5)=A(4)*B(2)+A(5)*B(5)+A(6)*B(8)
        C(6)=A(4)*B(3)+A(5)*B(6)+A(6)*B(9)
        C(7)=A(7)*B(1)+A(8)*B(4)+A(9)*B(7)
        C(8)=A(7)*B(2)+A(8)*B(5)+A(9)*B(8)
        C(9)=A(7)*B(3)+A(8)*B(6)+A(9)*B(9)
c         
       end subroutine tc_2s2
c
c Double Contraction of 2 2nd order Tensors
c
       subroutine tc_2d2 (A,B,c)
c
        implicit none
        double precision, intent (in) :: A(9),B(9)
        double precision, intent (out) :: c
c
        c = A(1)*B(1)+A(2)*B(4)+A(3)*B(7)
     &     +A(4)*B(2)+A(5)*B(5)+A(6)*B(8)
     &     +A(7)*B(3)+A(8)*B(6)+A(9)*B(9)
c         
       end subroutine tc_2d2
c
c
c Double Contraction of 2 4th order Tensors
c
       subroutine tc_4d4 (A,B,C)
c
        implicit none
        double precision, intent (in) :: A(9,9),B(9,9)
        double precision, intent (out) :: C(9,9)
        integer ij,kl
c
        do ij=1,9
           do kl=1,9
            C(ij,kl)=A(ij,1)*B(1,kl)+A(ij,2)*B(2,kl)+A(ij,3)*B(3,kl)+
     &               A(ij,4)*B(4,kl)+A(ij,5)*B(5,kl)+A(ij,6)*B(6,kl)+
     &               A(ij,7)*B(7,kl)+A(ij,8)*B(8,kl)+A(ij,9)*B(9,kl)
           end do
        end do 
c         
       end subroutine tc_4d4
c
c
c Double Contraction of 4th and 2nd order Tensors
c
       subroutine tc_4d2 (A4,B,C)
c
        implicit none
        double precision, intent (in) :: A4(9,9),B(9)
        double precision, intent (out) :: C(9)
        integer ij
c
        do ij=1,9
            C(ij)=A4(ij,1)*B(1)+A4(ij,2)*B(2)+A4(ij,3)*B(3)+
     &            A4(ij,4)*B(4)+A4(ij,5)*B(5)+A4(ij,6)*B(6)+
     &            A4(ij,7)*B(7)+A4(ij,8)*B(8)+A4(ij,9)*B(9)
        end do 
c         
       end subroutine tc_4d2
c
c
c Double Contraction of 2nd and 4th order Tensors
c
       subroutine tc_2d4 (A,B4,C)
c
        implicit none
        double precision, intent (in) :: A(9),B4(9,9)
        double precision, intent (out) :: C(9)
        integer ij
c
        do ij=1,9
            C(ij)=A(1)*B4(1,ij)+A(2)*B4(2,ij)+A(3)*B4(3,ij)+
     &            A(4)*B4(4,ij)+A(5)*B4(5,ij)+A(6)*B4(6,ij)+
     &            A(7)*B4(7,ij)+A(8)*B4(8,ij)+A(9)*B4(9,ij)
        end do 
c         
       end subroutine tc_2d4
c
c Transformation of second order tensor
c
       subroutine transform2 (A,Q,C)
c
	implicit none
        double precision, intent (in) :: A(9), Q(9)
        double precision, intent (out) :: C(9)
        double precision :: B(9)
        double precision :: Qt(9) 
c
c       calculate qt
c
        call tpose (Q,Qt)
        print *,
        call t2print(Qt)
        print *,
c
c        let      B = QA
c        we want  B = (QA)Qt
c        then     B = BQt
c
        call tc_2s2 (A,Qt,B)
        print *,
        call t2print (B)
        print *,
        call t2print (Q)
        print *,
        call tc_2s2 (Q,B,C)
c
        return
c
       end subroutine transform2
c
c Rotation tensor for rotation about z-axis
c        IN RADIANS!!!!!!!!!!!
c
       subroutine qzrad (th,Q)
c
        implicit none
        double precision, intent (in) :: th
        double precision, intent (out) :: Q(9)
c
        Q(1) = COS(th) 
        Q(2) = SIN(th)
        Q(3) = 0.0_8
        Q(4) = -SIN(th)
        Q(5) = COS(th)
        Q(6) = 0.0_8
        Q(7) = 0.0_8
        Q(8) = 0.0_8
        Q(9) = 1.0_8
        return
c
       end subroutine qzrad
c
c Rotation tensor for rotation about x-axis
c        IN RADIANS!!!!!!!!!!!
c
       subroutine qxrad (th,Q)
c
        implicit none
        double precision, intent (in) :: th
        double precision, intent (out) :: Q(9)
c
        Q(1) = 1.0_8
        Q(2) = 0.0_8
        Q(3) = 0.0_8
        Q(4) = 0
        Q(5) = COS(th)
        Q(6) = -SIN(th)
        Q(7) = 0.0_8
        Q(8) = SIN(th)
        Q(9) = COS(th)
        return
c
       end subroutine qxrad
c
c Rotation tensor for rotation about y-axis
c        IN RADIANS!!!!!!!!!!!
c
       subroutine qyrad (th,Q)
c
        implicit none
        double precision, intent (in) :: th
        double precision, intent (out) :: Q(9)
c
        Q(1) = COS(th) 
        Q(2) = 0.0_8
        Q(3) = SIN(th)
        Q(4) = 0.0_8
        Q(5) = 1.0_8
        Q(6) = 0.0_8
        Q(7) = -SIN(th)
        Q(8) = 0.0_8
        Q(9) = COS(th)
        return
c
       end subroutine qyrad
c
c Dyadic Product of 2 2nd order tensors
c
       subroutine tp_2dy2 (A,B,C4)
c
        implicit none
        double precision, intent (in) :: A(9), B(9)
        double precision, intent (out) :: C4(9,9)
        integer ij
        
        do ij=1,9
         C4(ij,1) = A(ij)*B(1)
         C4(ij,2) = A(ij)*B(2)
         C4(ij,3) = A(ij)*B(3)
         C4(ij,4) = A(ij)*B(4)
         C4(ij,5) = A(ij)*B(5)
         C4(ij,6) = A(ij)*B(6)
         C4(ij,7) = A(ij)*B(7)
         C4(ij,8) = A(ij)*B(8)
         C4(ij,9) = A(ij)*B(9)
        end do
c
        return
c
       end subroutine tp_2dy2
c
