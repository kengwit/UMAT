       subroutine kindef (F0, F1, b, dF, Fb, lb, Db, Wb)
c
         double precision, intent(in) :: F0(9), F1(9)
         double precision, intent(in) :: b
         double precision, intent(out) :: dF(9), Fb(9),lb(9),Db(9),Wb(9)
         double precision :: Ft (9)
c
         dF = F1-F0
c
c        Fb (b=0,0.5,1 for backward, midpoint, and forward)
c
         Fb = b*F0+F1-b*F1
         call tpose (Fb,Ft)
         call tc_2s2 (dF,Ft,lb)
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
        print *,cm(1),cm(2)
        call L4_EL (cm(1), cm(2), L4)
c
c        call t4print (L4)
        dsig = tc_4d2(L,deps)
        sig = sig + dsig
c        
        return
c
       end subroutine
c
