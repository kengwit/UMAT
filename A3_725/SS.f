        program simpleshear
c
          implicit none
c
          double precision :: dgam, gam, gf, b
          double precision :: F0(9), F1(9) 
          double precision :: F(9),dF(9),l(9),D(9),W(9),R(9),logc
          double precision :: ep(9), deps(9), eR(9), Re(9)
          double precision :: sig (9), sigR(9), Rsig(9)
          double precision :: sigp(9), epp(9)
          double precision :: cm(5), hisv(11), gamout(6)
          integer incout(6),place
          integer i,inc
c
          dgam = 0.0005_8
          gam = 0.0_8
          gf = 6.001_8
          b = 0.5_8
          gamout = (/ 0.005_8,0.1_8,0.5_8,1.0_8,3.0_8,6.0_8 /)
          incout = gamout / dgam
          incout = incout + 1
          place = 1
c
          inc = gf/dgam
          cm(1) = 205000_8
          cm(2) = 0.25_8
          cm(3) = 1300.0_8
          cm(4) = 0.002_8
          cm(5) = 0.07_8
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
             if (i.ne.1) then
                logc = dgam/(4.0_8+gam*gam)
                logc = logc + dgam*gam*0.25_8/(SQRT(4+gam*gam)
     &                        *ASINH(gam/2))
                R = logc*(/ 0,1,0,-1,0,0,0,0,0 /)
             else
	        R = 0.0_8
             end if			 
c            R = W
            call tc_2s2 (ep,R,eR)
            call tc_2s2 (R,ep,Re)
            deps = D - eR + Re
            ep = ep + deps
c
c           Perform objective update on stress (-sigW + Wsig, the objective update)
c
            call tc_2s2 (sig,R,sigR)
            call tc_2s2 (R,sig,Rsig)
            sig = sig - sigR + Rsig           
c
c           Perform stress update on stress (sig dot hat, the stress update)
c
            call CH (ep, epp)
	    call CH (sig, sigp)
c            call umat_elastic (cm,deps,sig,hisv)
            call umat_43 (cm,sig,deps,hisv,0.0_8,0.0_8)
            if (i.eq.incout(place)) then
               print "(27e15.6)",gam,ep,sig,
     &               epp(1),epp(5),epp(9),sigp(1),sigp(5),sigp(9),
     &               hisv(1),hisv(11)
               place = place + 1
c
            end if
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
c            print *,"R"
c            call t2print (R)
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
