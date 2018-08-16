c
c   Subroutine to find the JH2 yield stress, which is dependent on 
c   pressure, strain rate and damage accumulated. All values are
c   normalized with respect to the HEL (stress wrt SHEL, pressure
c   wrt PHEL). 
c  
c   inputs
c   ======
c   epd = equivalent strain rate         (history var)
c   P = current hydrostatic pressure     (history var)
c   D = current damage                   (history var)
c   A = intact strength coefficient      (JH2 param.)
c   B = fractured strength coefficient   (JH2 param.)
c   c = strain rate coefficient          (JH2 param.)
c   n = intact strength exponent         (JH2 param.)
c   m = fractured strength exponent      (JH2 param.)
c   EPSI = quasi-static threshold s.rate (JH2 param.)
c   T = hydrostatic pressure to failure  (JH2 param.)
c   SFmax = maximum failure strength     (JH2 param.)
c   PHEL = Pressure comp. of HEL         (JH2 param.)
c   SHEL = Deviatoric comp. of HEL       (JH2 param.)
c
c   outputs
c   =======
c   syield = yield stress
c
        subroutine JH2_yield (
     &             epd,P,D,A,B,c,n,m,EPSI,T,SFmax,HEL,PHEL,syield)
c
         implicit none
         double precision, intent (in) :: epd,P,D,A,B,c,n,m,T,
     &                                    SFmax,HEL,PHEL,EPSI
         double precision, intent (inout) :: syield
         double precision :: si, sf, Ps, Ts, es, SHEL
c
c        Find normalized pressure and tensile strength constant (wrt PHEL)
c        =================================================================
         Ps = P / PHEL
         Ts = T / PHEL
         SHEL = 1.5_8 * (HEL - PHEL)
c
c        Find normalized strain rate (wrt EPSI) and find strain rate term
c        ================================================================
         es = 1.0_8+c*log(epd/EPSI)
c
c        Find normalized strength at D=0 (si) , D=1 (sf)
c        ===============================================
         si = A * (Ps+Ts) ** n * es
         sf = B * Ps ** m * es
c
c        Check if normalized fractured strength exceeds SFmax, if so sf=SFmax
c        ====================================================================
         if (sf.GT.SFmax) then
            sf = SFmax
         end if
c
c        Calculate flow stress based on si, sf and damage
c        ================================================
         syield = SHEL * (si - D*si + D*sf)
c 
         return
c
        end subroutine JH2_yield
c
c  Subroutine to update damage in the JH2 model 
c
c  inputs
c  ======
c  eps = effective plastic strain increment
c  P = current hydrostatic pressure
c  T = hydrostatic tensile stress to failure
c  PHEL = pressure component of HEL
c  D1 = damage coefficient
c  D2 = damage exponent
c   
c  outputs
c  =======
c  D = damage parameter (to be updated in subroutine)
c
c
        subroutine JH2_damage (eps,P,T,PHEL,D1,D2,D)
c
         double precision, intent (in) :: eps,P,T,PHEL,D1,D2
         double precision, intent (inout) :: D
         double precision :: epf
c
         Ps = P / PHEL
         Ts = T / PHEL
c
         epf = D1 * (Ps + Ts) ** D2
         D = D + eps / epf
c
         return
c
        end subroutine JH2_damage
c
c  Subroutine to calculate equivalent stress as well as normals
c  required for the plastic loop
c
c  inputs
c  ======
c  A = input tensor
c  cflag = changes outputs of the function (note only cflag=2 is implemented to date)
c
c  outputs
c  =======
c  seq = von Mises stress as output
c  s1d = 1st derivatives of Seq wrt S (is the nomal tensor in case of von Mises)
c  s2d = 2nd derivatives of Seq	wrt S (apparently used for implicit formulations, not implemented here)
c
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
         call VMS (A,seq)
c
         s1d = 0.0_8
         s2d = 0.0_8

         if (cflag.eq.2) then
            call A2S (A,S)
c            call t2print (S)
            s1d = (1.5_8/seq)*S
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
c  inputs
c  ======
c  P = hydrostatic pressure
c  D = damage
c  A = intact strength coefficent 
c  B = damaged strength coefficient
c  M = damaged strength exponent
c  N = intact strength exponent
c  HEL = hugoniot elastic limit
c  PHEL = hydrostatic pressure at HEL
c  Lel = elastic tensor
c  N = normal vector (from vm_box)
c  D1 = damage accumulation coefficient
c  D2 = damage accumulation exponent
c
c  outputs
c  =======
c  h = hardening rate (ds_bar/dep)
c
        subroutine JH2_hardening 
     &  (epd,P,D,A,B,c,M,N,T,HEL,PHEL,SFmax,EPSI,Lel,NOR,D1,D2,h)

c
          implicit none
          double precision, intent (in) :: epd,P,D,A,B,c,M,N,T
          double precision, intent (in) :: HEL,PHEL,SFmax,EPSI
          double precision, intent (in) :: Lel(9,9), NOR(9), D1, D2
          double precision, intent (inout) :: h
          double precision :: es, si, sf, ef, Ps, Ts, SHEL
          double precision :: alp, bet, kap, KIJ(9), LN(9), KLN
c
c
c        Find normalized pressure and tensile strength constant (wrt PHEL)
c        =================================================================
         Ps = P / PHEL
         Ts = T / PHEL
         SHEL = 1.5 * (HEL - PHEL)
c
c        Find normalized strain rate (wrt EPSI) and find strain rate term
c        ================================================================
         es = 1.0_8+c*log(epd/EPSI)
c
c        Find normalized strength at D=0 (si) , D=1 (sf)
c        ===============================================
         si = A * (Ps+Ts) ** n * es
         sf = B * Ps ** m * es
c
c        Check if normalized fractured strength exceeds SFmax, if so sf=SFmax
c        ====================================================================
         if (sf.GT.SFmax) then
            sf = SFmax
         end if
c
c        Find failure strain at given pressure
c        =====================================
         ef = D1 ** (Ps+Ts) ** D2
         if (ef.lt.1e-20) then
            ef = 1e-20
         end if
c
c        Find intermediates to find total derivative
c        ===========================================
         call KD (KIJ)
         call tc_4d2 (Lel,NOR,LN)
         call tc_2d2 (KIJ,LN,KLN)
         kap = Ps * KLN / PHEL / 3.0_8
         alp = A * N * (Ps-Ts)**(n-1.0_8)
         bet = B * m * (Ps)**(m-1.0_8)
c
         h = SHEL * ((1.0_8-D)*alp*kap+D*bet*kap+(sf-si)/ef) 
c
         return
c
        end subroutine JH2_hardening
c
c      
        subroutine umat_JH2 (cm, sig, deps, hisv, tt, dt)
c
          implicit none
          double precision, intent (in) :: cm (20), deps (9), tt, dt
          double precision, intent (inout) :: sig (9), hisv (12)
          double precision :: L4 (9,9), strial (9), dsig(9), E, v
          double precision :: seq, sbar, h, sflow, NOR(9), ph (9,9)
          double precision :: LN (9), NLN, pl_corr(9), SHEL, syield
          double precision :: MID,RO,G,A,B,C,M,N,EPSI,T,SFmax,HEL,PHEL
          double precision :: BETA, D1, D2, K1, K2, K3, FS
          double precision :: P,D, dlam, ddlam, dep, ef, Ps, Ts, sbar_t  
          double precision :: tol
          integer cnt
c
c         cm(1) = MID   , material id          
c         cm(2) = RO    , density
c         cm(3) = G     , shear modulus
c         cm(4) = A     , intact material strength coefficient
c         cm(5) = B     , fractured material strength exponent
c         cm(6) = C     , strain rate coefficent
c         cm(7) = M     , fractured material strength exponent
c         cm(8) = N     , intact material strength coefficient
c         cm(9) = EPSI  , quasi-static strain rate threshold
c         cm(10) = T    , hydrostatic tensile strength
c         cm(11) = SFmax, maximum normalized fractured strength
c         cm(12) = HEL  , hugoniot elastic limit 
c         cm(13) = PHEL , pressure component of HEL
c         cm(14) = BETA , fraction of elastic energy converted to pressure 
c         cm(15) = D1   , damage accumulation coefficient
c         cm(16) = D2   , damage accumulation exponent
c         cm(17) = K1   , EOS term 1 (if K2,K3=0, K1 is the bulk mod.)
c         cm(18) = K2   , EOS term 2
c         cm(19) = K3   , EOS term 3
c         cm(20) = FS   , failure criteria
c
         MID = cm(1)
         RO = cm(2)
         G = cm(3)
         A = cm(4)
         B = cm(5)
         C = cm(6)
         M = cm(7)
         N = cm(8)
         EPSI = cm(9)
         T = cm(10)
         SFmax = cm(11)
         HEL = cm(12) 
         PHEL = cm(13)
         BETA = cm(14)
         D1 = cm(15)
         D2 = cm(16)
         K1 = cm(17)
         K2 = cm(18)
         K3 = cm(19)
         FS = cm(20)
         D = hisv(12)
         call I1 (sig, P)
         P = P / 3.0_8
         syield = 0.0_8
c
c ================================
c    1. Calculate Trial Stress
c ================================                   
c
c         Get elastic constants from G,K
c         ==============================
          E = 9.0_8*G*K1/(3.0_8*K1+G)
          v = (3.0_8*G-2.0_8*K1)/(2.0_8*(3.0_8*K1+G))
          SHEL = 1.5_8*(HEL-PHEL)
c
c         Get elastic fourth order tensor L4
c         ==================================
          call L4_EL (E,v,L4)
c
c         Calculate stress increment (dsig = L:D)
c         =======================================
          call tc_4d2 (L4,deps,dsig)
c
c         Update strial with dsig (sn+1,trial = sn + dsig)
c         ================================================
          strial = sig + dsig
c
c ================================
c      2. Find Seq and Syield
c ================================
c
c         Find equivalent stress and normal using vm_box
c         ==============================================
          call vm_box (sig,2,seq,NOR,ph)
c
c         Find yield stress using JH2_yield
c         =================================
          call JH2_yield 
     &    (0.0_8,P,D,A,B,c,M,N,EPSI,T,SFmax,HEL,PHEL,syield)
c
c         Check for yielding
c         ==================
          sflow = seq - syield
c
          if (sflow.lt.0.0_8) then
             sig = strial
          else
c         =====================
c         Initiate PLASTIC LOOP
c         =====================
          call JH2_hardening 
     &    (0.0_8,P,D,A,B,c,M,N,T,HEL,PHEL,SFmax,EPSI,L4,NOR,D1,D2,h)
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
            do while (DABS(sflow).gt.tol.AND.cnt.lt.5)
c
c            +  use normal and elastic tensor to calculate
c            +  plastic correction
c
               call tc_4d2(L4,NOR,LN)
               pl_corr = dlam * LN
               sig = strial - pl_corr
c
c            +  recalculate seq and sflow to determine whether
c            +  material is still yielding
c
               call VMS (sig,seq)
               sbar_t = sbar + h*dep
               sflow = seq - sbar_t
c
c            +  update plastic multiplier increment, plastic multiplier
c            +  and plastic strain (the same in assoc. flow) to reflect
c            +  plastic correction
c               
               call tc_2d2 (NOR,LN,NLN)
               ddlam = sflow / (NLN+h)
               dlam = dlam + ddlam
               dep = dlam
               cnt = cnt + 1
c
            end do
            hisv(1) = hisv(1) + dep
            hisv(11) = seq
c
c           +  update damage variable based on plastic strain increment
c
            Ps = P / PHEL
            Ts = T / PHEL
            ef = D1 ** (Ps+Ts) ** D2
            D = D + dep / ef
            hisv(12) = D
          end if
c         return
c
        end subroutine umat_JH2
c
