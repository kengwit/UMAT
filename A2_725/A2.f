c
        program A2_725
c
         double precision :: F(9),R(9),U(9),V(9),Vp(9),Up(9)
         double precision :: VP1,VP2,VP3
         double precision :: UP1,UP2,UP3
         double precision :: eu1,eu2,eu3
         double precision :: la1,la2,la3
         double precision :: gam, lam
         gam = 2.0_8
         lam = 0.5_8
c
c    Simple Shear Def. Grad.
c
c         F = (/ 1.0_8, gam, 0.0_8,
c     &          0.0_8, 1.0_8, 0.0_8,
c     &          0.0_8, 0.0_8, 1.0_8 /)
c
c   Plain Strain Stretching Def. Grad.
c
c         F = (/ 1.0_8, gam, 0.0_8,
c     &          0.0_8, lam, 0.0_8,
c     &          0.0_8, 0.0_8, 1.0_8/lam /)
c
c   Equal Biaxial Stretching Def. Grad.
c
         F = (/ lam, gam, 0.0_8,
     &          0.0_8, lam, 0.0_8,
     &          0.0_8, 0.0_8, 1/(lam*lam) /)
c
        call PD (F,R,U,V)
c
        call CH (V,Vp)
        call CH (U,Up)
        VP1 = Vp(1)
        VP2 = Vp(5)
        VP3 = Vp(9)
        eu1 = LOG(VP1)
        eu2 = LOG(VP2)
        eu3 = LOG(VP3)
c
c        UP1 = Up(1)
c        UP2 = Up(5)
c        UP3 = Up(9)
c        la1 = LOG(UP1)
c        la2 = LOG(UP2)
c        la3 = LOG(UP3)
c
        print "(3e12.4)",eu1,eu2,eu3
c
c        print "(3e12.4)",la1,la2,la3
c
        end program A2_725
