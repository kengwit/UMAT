c
c       simple Swift's hardening law      
c
c       
        subroutine iso_hard (K,eo,n,ep,Sb,h)
c
         double precision, intent (in) :: K,eo,n,ep
         double precision, intent (out) :: Sb,h
c
         Sb = K * (eo+ep) ** n
         h = K*n*(eo+ep)**(n-1.0_8)
c
         return
c
        end subroutine iso_hard
c
c       COMMENT THIS!
c
        subroutine vm_box (A, cflag, seq, s1d, s2d)
c
         implicit none
         double precision, intent (in) :: A(9)
         double precision, intent (out) :: seq, s1d(9), s2d(9,9)
         double precision :: S(9) 
         integer, intent(in) :: cflag
c
c
         seq = SQRT(A(1)*A(1)+A(5)*A(5)+A(9)*A(9)
     &        -A(1)*A(5)-A(5)*A(9)-A(9)*A(1)
     &        +3.0_8*A(6)*A(6)+3.0_8*A(7)*A(7)+3.0_8*A(2)*A(2))
c
         s1d = 0.0_8
         s2d = 0.0_8

         if (cflag.eq.2) then
            call A2S (A,S)
            s1d = 1.5_8/seq*S
         else if (cflag.eq.3) then
         else
         end if
c
         return
c
        end subroutine vm_box
c
c
c
        subroutine umat_43 (cm, sig, deps, hisv, tt, dt)
c
         implicit none
         double precision, intent (in) :: cm(5), deps (9), tt, dt
         double precision, intent (inout) :: sig (9), hisv(11)
         double precision :: L4 (9,9), strial (9), dsig(9)
         double precision :: sbar, seq, h, sflow, N(9), ph (9,9)
         double precision :: dlam, ddlam, dep,tol,sflow_t, pl_corr(9)
         double precision :: LN(9), NLN
         integer :: cnt
         tol = cm(3)*1e-8_8
c
c        cm(1) = E (Young's Modulus) {Pa}
c        cm(2) = v (Poisson's ratio) {-}
c        cm(3) = K (Hardening coeff) {Pa)
c        cm(4) = eo (Onset strain? ) {m/m}
c        cm(5) = n (Hardening exp. ) {-}
c
c ==========================================
c           1. Calculate Trial Stress
c ==========================================
         call L4_EL (cm(1), cm(2), L4)
         call tc_4d2 (L4,deps,dsig)
c         print "(3e12.4)",dsig
c         print *,
         strial = sig + dsig
c ==========================================
c           2. Find Eqv. and Yield Stress
c ==========================================
         call vm_box (sig,2,seq,N,ph)
         call iso_hard (cm(3),cm(4),cm(5),hisv(1),sbar,h)
         sflow = seq - sbar
c ==========================================
c           3. Check for Yielding
c ==========================================
         if (sflow.lt.0.0_8) then
            sig = strial
c ==========================================
c           4. Initiate Plastic Loop
c ==========================================
         else
c           
c  prime the loop, set all counters to 0
c  set sflow = 100 to ensure this runs once
c
            dlam = 0.0_8
            dep = 0.0_8
            cnt = 1
            sflow = 100.0_8
c
c           plastic loop time!
c
            do while (sflow.gt.tol.OR.cnt.lt.5)
c
c            +  use normal and elastic tensor to calculate
c            +  plastic correction
c
               call tc_4d2(L4,N,LN)
               pl_corr = dlam * LN
               sig = strial - pl_corr
c
c            +  recalculate seq and sflow to determine whether
c            +  material is still yielding
c
               call VMS (sig,seq)
               sflow_t = sflow + h*dep
               sflow = seq - sflow_t
c
c            +  update plastic multiplier increment, plastic multiplier
c            +  and plastic strain (the same in assoc. flow) to reflect
c            +  plastic correction
c               
               call tc_2d2 (N,LN,NLN)
               ddlam = sflow / (NLN+h)
               dlam = dlam + ddlam
               dep = dlam
               cnt = cnt + 1
c
            end do
            hisv(1) = hisv(1) + dep
            hisv(11) = seq
         end if
         return
c
        end subroutine umat_43
c
