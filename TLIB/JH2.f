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
