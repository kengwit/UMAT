       subroutine kindef (F0, F1, b, dF, Fb, lb, Db, Wb)
c
         double precision, intent(in) :: F0(9), F1(9)
         double precision, intent(in) :: b
         double precision, intent(out) :: dF(9), Fb(9),lb(9),Db(9),Wb(9)
         double precision :: Finv(9)
c
         dF = F1-F0
c
c        Fb (b=0,0.5,1 for backward, midpoint, and forward)
c
         Fb = F0-b*F0+b*F1
         call invert (Fb,Finv)
         call tc_2s2 (dF,Finv,lb)
         call SYM2 (lb,Db)
         call SKEW2 (lb,Wb)
         return
c
       end subroutine kindef
c
       subroutine umat_elastic (cm, deps, sig, hisv)
c
        double precision, intent (in) :: cm (2),deps(9)
        double precision, intent (inout) :: sig(9),hisv(1)
        double precision :: dsig(9)
        double precision :: L4(9,9)
c
c   cm(1) = E (Young's modulus) {Pa}
c   cm(2) = v (Poisson's ratio) {-}        
c
        call L4_EL (cm(1), cm(2), L4)
c
        call tc_4d2(L4,deps,dsig)
c        call t2print (dsig)
        sig = sig + dsig
c        
        return
c
       end subroutine
c
c       subroutine drate (D,W,B2,O)
cc
c        double precision, intent (in) :: D(9),W(9),B2(9)
c        double precision, intent (out) :: O(9)
c        double precision :: b(3), Bp(9)
c        double precision :: BD(9), BB (9) , BBD(9),BBDB(9)
c        double precision :: v0,v(3),alp,tol
c        double precision :: u(3)
c        integer :: i,k
cc        
c        call CH (B2,Bp)
c        b(1) = Bp(1)
c        b(2) = Bp(5)
c        b(3) = Bp(9)
c        print *,"b:"
c        print "(3e12.4)",b(1),b(2),b(3)
cc
c        tol = 1.0e-8_8
cc
c        call tc_2s2(B2,D,BD)
c        u(1) = b(2)/b(3)
c        u(2) = b(3)/b(1)
c        u(3) = b(1)/b(2)
c        print *,"u:"
c        print "(3e12.4)",u(1),u(2),u(3)
cc
c        v = 0.0_8
cc
c        if (abs(b(1)-b(2)).lt.tol .and.
c     &      abs(b(2)-b(3)).lt.tol .and.
c     &      abs(b(3)-b(1)).lt.tol) then
c          O = 0.0_8
c          print * , "I'm in case 1"
c        else if (abs(b(1)-b(2)).gt.tol .and.
c     &           abs(b(2)-b(3)).lt.tol) then
c          v0 = 1/(b(1)-b(2))*(((1+u(3))/(1-u(3)))+(2/LOG(u(3))))
c          O = v0 * BD
c          print * , "I'm in case 2"
c        else 
c          call tc_2s2(B2,B2,BB)
c          call tc_2s2(BB,D,BBD)
c          call tc_2s2(BB,DB,BBDB)
c          alp = -1/((b(1)-b(2))*(b(2)-b(3))*(b(3)-b(1)))
c          do k=1,3
c             do i=1,3
cc                 v(k) = v(k) + 1
c                v(k) = v(k) + ((-1*b(i)**(3-k))
c     &                 *(((1+u(i))/(1-u(i)))+(2/LOG(u(i)))))
cc                print *,((1+u(i))/(1-u(i)))
c                print *,(2/LOG(u(i)))
c             end do
c          end do
c          O = v(1)*BD+v(2)*BBD+v(3)*BBDB
cc          print * , "I'm in case 3"
c          print *,"v:"
c          print "(3e12.4)",v(1),v(2),v(3)
c        end if
cc
c        O = W + O
c        return
cc
c       end subroutine drate
c
