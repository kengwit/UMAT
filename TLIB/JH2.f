c
c   Subroutine to find the JH2-flow stress, which is dependent on 
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
c   sflow = flow stress
c
c
        subroutine JH2_flow (
     &             epd,P,D,A,B,c,n,m,EPSI,T,SFmax,PHEL,SHEL,sflow)
c
         implicit none
         double precision, intent (in) :: epd,P,D,A,B,c,n,m,T,
     &                                    SFmax,PHEL,SHEL,EPSI
         double precision, intent (inout) :: sflow
         double precision :: si, sf, Ps, Ts, es
c
         Ps = P / PHEL
         Ts = T / PHEL
         es = 1.0_8+c*log(epd/EPSI)
         si = A * (Ps+Ts) ** n * es
         sf = B * Ps ** m * es
         if (sf.GT.SFmax) then
            sf = SFmax
         end if
         sflow = SHEL * (si - D*si + D*sf)
c
         return
c
        end subroutine JH2_flow
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
c
c
c      
        subroutine umat_JH2 (cm, sig, deps, hisv, tt, dt)
c
          implicit none
          double precision, intnet (in) :: cm (20), deps (9), tt, dt
          double precision, intent (inout) :: sig (9), hisv (12)
          double precision :: L4 (9,9), strial (9), dsig(9)
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
          return
c
        end subroutine umat_JH2 ()
c
