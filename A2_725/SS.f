        program simpleshear
c
          implicit none
c
          double precision :: dgam, gam, gf, b
          double precision :: F0(9), F1(9) 
          double precision :: F(9),dF(9),l(9),D(9),W(9)
          double precision :: ep(9), deps(9), eW(9), We(9)
          double precision :: sig (9), sigW(9), Wsig(9)
		  double precision :: sigp(9), epp(9)
          double precision :: cm(2), hisv(1)
          integer i,inc
c
          dgam = 0.001_8
          gam = 0.0_8
          gf = 3.0_8
          b = 0.5_8
c
          inc = gf/dgam
          cm(1) = 205e9_8
          cm(2) = 0.25_8
c
          ep = 0.0_8
          deps = 0.0_8
          sig = 0.0_8
          D = 0.0_8
          F1 = (/ 1.0_8, 0.0_8, 0.0_8,
     &            0.0_8, 1.0_8, 0.0_8,
     &            0.0_8, 0.0_8, 1.0_8 /) 
c          
          do i=1,inc
c
c           Get kinematics of deformation
c 
            F0 = F1
            F1 = (/ 1.0_8,   gam, 0.0_8,
     &              0.0_8, 1.0_8, 0.0_8,
     &              0.0_8, 0.0_8, 1.0_8 /)
c
            call kindef (F0,F1,b, dF,F,l,D,W)
c
c           Perform objective update on strain
c
            call tc_2s2 (ep,W,eW)
            call tc_2s2 (W,ep,We)
            deps = D - eW + We
            ep = ep + deps
c
c           Perform objective update on stress (-sigW + Wsig, the objective update)
c
            call tc_2s2 (sig,W,sigW)
            call tc_2s2 (W,sig,Wsig)
            sig = sig - sigW + Wsig           
c
c           Perform stress update on stress (sig dot hat, the stress update)
c
            call CH (ep, epp)
			call CH (sig, sigp)
            call umat_elastic (cm,deps,sig,hisv)
            print "(25e15.6)",gam,ep,sig,
     &            epp(1),epp(5),epp(9),sigp(1),sigp(5),sigp(9)
c
c            print *,
c            print *,"gamma:"
c            print "(1e12.4)",gam
c            print *,
c            print *,"ep"           
c            call t2print (ep)
c            print *,
c            print *,"deps"
c            call t2print (deps)
c            print *,
c            print *,"W"
c            call t2print (W)
c            print *,
c            print *,"D"
c            call t2print (D)
c            print *,
c            print *,"sig"
c            call t2print (sig)
c            print *,
            gam = gam + dgam
c         
          end do
c                    
        end program simpleshear
c
