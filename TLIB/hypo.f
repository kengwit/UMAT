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
       subroutine drate (D,W,B,O)
c
        double precision, intent (in) :: W(9),D(9),B(9)
c
        return
c
       end subroutine drate
