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
     &             epd,P,D,A,B,c,m,n,EPSI,T,SFmax,HEL,PHEL,syield)
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
c         es = 1.0_8+c*log(epd/EPSI)
          es = 1.0_8
c
c        Find normalized strength at D=0 (si) , D=1 (sf)
c        ===============================================
         if ((Ps+Ts).gt.0.0_8) then
             si = A * (Ps+Ts) ** n
         else
             si = 0.0_8
         end if
c
         if (Ps.gt.0.0_8) then
             sf = B * Ps ** m
         else
             sf = 0.0_8
         end if
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
         double precision :: S(9),vm
         integer, intent(in) :: cflag
c
c
          call VMS (A,seq)
          vm = seq
          if (seq.lt.1e-20) then
             vm = 1e-20
          end if
c
         s1d = 0.0_8
         s2d = 0.0_8

         if (cflag.eq.2) then
            call A2S (A,S)
            s1d = (1.5_8/vm)*S
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
          double precision :: alp, bet, kap, KIJ(9), LK(9), NLK
          double precision :: x,y,z
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
         if (Ps+Ts.gt.0.0_8) then
             si = A * (Ps+Ts) ** N
         else
             si = 0.0_8
         end if
c
         if (Ps.gt.0.0_8) then
             sf = B * Ps ** M
         else 
             sf = 0.0_8
         end if
c
c        Check if normalized fractured strength exceeds SFmax, if so sf=SFmax
c        ====================================================================
         if (sf.GT.SFmax) then
             sf = SFmax
         end if
c
c        Find failure strain at given pressure
c        =====================================
         if ((Ps+Ts).gt.0.0) then
             ef = D1 * (Ps+Ts) ** D2
         else
             ef = 1e-20_8
         end if
c
c        Find intermediates to find total derivative
c        ===========================================
         call KD (KIJ)
         call tc_4d2 (Lel,KIJ,LK)
         call tc_2d2 (N,LK,NLK)
         kap =  NLK / 3.0_8
         if ((Ps+Ts).gt.0.0_8) then
             alp = SHEL/PHEL * A * N * (Ps+Ts)**(N-1.0_8)
         else
             alp = 0.0_8
         end if
c
         if (Ps.gt.0.0_8) then
             bet = SHEL/PHEL * B * M * (Ps)**(M-1.0_8)
         else
             bet = 0.0_8
         end if
c     
         x = (1.0_8-D)*alp*kap
         y = D*bet*kap
         z = SHEL*(sf-si)/ef
         if ((D.lt.0.001).or.(D.gt.0.9999)) then
            z = 0.0_8
         end if
         h = x+y+z 
         if (DABS(h).lt.1e-20_8) then
             h = 1e-20
         end if
c
         return
c
        end subroutine JH2_hardening
c
c      
        subroutine umat_JH2 (cm, Sio, Dio, hisv, tt, dt, efail)
c
          implicit none
          double precision, intent (in) :: cm (*), Dio (*), tt, dt
          integer, intent (inout) :: efail
          double precision, intent (inout) :: Sio (*), hisv (*)
          double precision :: L4 (9,9), strial (9), dsig(9), E, v
          double precision :: seq, sbar, h, sflow, NOR(9), ph (9,9)
          double precision :: LN (9), NLN, pl_corr(9), SHEL, syield
          double precision :: MID,RO,G,A,B,C,M,N,EPSI,T,SFmax,HEL,PHEL
          double precision :: BETA, D1, D2, K1, K2, K3, FS
          double precision :: P,D, dlam, ddlam, dep, ef, Ps, Ts, sbar_t  
          double precision :: tol, epd, deps(9), sig(9)
          integer cnt
c
          sig(1) = Sio (1)
          sig(5) = Sio (2)
          sig(9) = Sio (3)
          sig(2) = Sio (4)
          sig(3) = Sio (5)
          sig(6) = Sio (6)
          sig(4) = sig (2)
          sig(7) = sig (3)
          sig(8) = sig (6)
c
          deps(1) = Dio (1)
          deps(5) = Dio (2)
          deps(9) = Dio (3)
          deps(2) = 0.5_8*Dio (4)
          deps(3) = 0.5_8*Dio (5)
          deps(6) = 0.5_8*Dio (6)
          deps(4) = deps (2)
          deps(7) = deps (3)
          deps(8) = deps (6)

c
c         cm(1) = K     , bulk modulus (for timestep)          
c         cm(2) = G     , shear modulus (for timestep)
c         cm(3) = DAM   , element deletion (not used)
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
c         cm(17) = D1   , damage accumulation coefficient
c         cm(18) = D2   , damage accumulation exponent
c         cm(19) = K1   , EOS term 1 (if K2,K3=0, K1 is the bulk mod.)
c         cm(20) = K2   , EOS term 2
c         cm(21) = K3   , EOS term 3
c         cm(22) = FS   , failure criteria
c
c         print *,"JH2"
c         print "(8e12.4)",cm(1),cm(2),cm(3),cm(4),
c     &                    cm(5),cm(6),cm(7),cm(8)
c         print "(8e12.4)",cm(9),cm(10),cm(11),cm(12),
c     &                    cm(13),cm(14),cm(15),cm(16)
c         print "(8e12.4)",cm(17),cm(18),cm(19),cm(20),
c     &                    cm(21),cm(22),cm(23),cm(24) 
c         print *,""
c
c         print *,"HISV"
c         print "(4e12.4)",hisv(1),hisv(2),hisv(3),hisv(4),
c     &                    hisv(5),hisv(6),hisv(7),hisv(8),
c     &                    hisv(9),hisv(10),hisv(11),hisv(12)
c         print *,""
c
         K1 = cm(1)
         G = cm(2)
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
         D1 = cm(17)
         D2 = cm(18)
         K1 = cm(19)
         K2 = cm(20)
         K3 = cm(21)
         FS = cm(22)
         D = 0.0_8
         call I1 (sig, P)
         P = -P / 3.0_8
         syield = 0.0_8
         tol = G * 1.0e-7_8
         epd = 0.0_8
c
         if (tt.eq.0.0) then
            hisv(1) = 0.0
            hisv(2) = 0.0
            hisv(3) = 0.0
            hisv(4) = 0.0
            print *,"Material Constants:"
            print *,"G:",G
            print *,"K:",K1
            print *,""
         end if
c
         if (ABS(epd).lt.1e-20) then
            epd = 1e-20
         end if
c
c ================================
c    1. Calculate Trial Stress
c ================================                   
c
c         Get elastic constants from G,K
c         ==============================
          E = 9.0_8*G*K1/(3.0_8*K1+G)
          v = (3.0_8*K1-2.0_8*G)/(2.0_8*(3.0_8*K1+G))
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
     &    (epd,P,D,A,B,c,M,N,EPSI,T,SFmax,HEL,PHEL,syield)
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
     &    (epd,P,hisv(2),A,B,c,M,N,T,HEL,PHEL,SFmax,EPSI,L4,NOR,D1,D2,h)
c           
c  prime the loop, set all counters to 0
c  set sflow = 100 to ensure this runs once
c
          if ((P+T).gt.0.0_8) then
            dlam = 0.0_8
            dep = 0.0_8
            cnt = 1
            sflow = tol*100
            call tc_4d2(L4,NOR,LN)
c
c           plastic loop time!
c
            do while (DABS(sflow).gt.tol.AND.cnt.lt.5)
c
c            +  use normal and elastic tensor to calculate
c            +  plastic correction
c
               pl_corr = dlam * LN
               sig = strial - pl_corr
c
c            +  recalculate seq and sflow to determine whether
c            +  material is still yielding
c
               call VMS (sig,seq)
               sbar_t = syield + h*dep
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
            end do
c
            hisv(1) = hisv(1) + dep
            hisv(3) = seq            
c
c           +  update damage variable based on plastic strain increment
c
            Ps = P / PHEL
            Ts = T / PHEL
c
            ef = D1 * (Ps+Ts) ** D2
            D = dep/ef
c           
            hisv(2) = hisv (2) + ABS(D)
c            
            if (hisv(2).gt.1.0) then
               hisv (2) = 1.0_8
            end if
c
c
          else
            hisv(2) = 1.0_8
            hisv(4) = 1.0_8
            call EEQ (deps,dep)
            hisv(1) = hisv(1) + dep
            sig = 0.0_8
            efail = 1
c (P+T)>0
         end if
c (Sflow)>0
         end if
c
          Sio (1) = sig(1)
          Sio (2) = sig(5)
          Sio (3) = sig(9)
          Sio (4) = 0.5_8*(sig(2)+sig(4))
          Sio (5) = 0.5_8*(sig(3)+sig(7))
          Sio (6) = 0.5_8*(sig(6)+sig(8))
c
          return
c
        end subroutine umat_JH2
c
