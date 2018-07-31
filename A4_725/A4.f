       program A4_725
c
         implicit none
c
c        parameters of interest: stress, strain, eq. plastic strain
c        THESE ARE OUTPUTS
c
         double precision :: sig(9), ep (9), eps
         double precision :: dtrial (9), strial (9)
         double precision :: d33trial (7), s33trial (7)
         double precision :: deps (9), cm (5), hisv (11)
         double precision :: xnm1, xnm2, fnm1, fnm2, deltas
         double precision :: tol
         integer :: i,final,cnt
c
c        cm(1): Young's modulus
c        cm(2): Poisson's Ratio
c        cm(3): Swift's K
c        cm(4): Swift's e0
c        cm(5): Swift's n
c
         cm(1) = 205000_8
         cm(2) = 0.25_8
         cm(3) = 1300.0_8
         cm(4) = 0.002_8
         cm(5) = 0.07_8
c
         sig = 0.0_8
         ep  = 0.0_8
         deps = 0.0_8
         eps = 0.0_8
         hisv = 0.0_8
	 tol = cm(3) * 1e-8_8
c
         dtrial = 0.0_8
         strial = 0.0_8
         d33trial = 0.0_8
         s33trial = 0.0_8
         cnt = 1
         final = 1000
c
         do i = 1,final 
c
              deps = 0.0_8
              deps (1) = 0.001_8
              eps = hisv(1)
c
c             initial guess for D33
c
              deps (9) = 0.000_8
c
c  The following code is an implementation of the Secant
c  method as documented on the wiki page for the Secant 
c  method... It is pretty clunky as I am iterating off a
c  umat call. Because of this, I have to make sure that I am
c  not affecting the stress and strain state of the previous 
c  timestep. All states should be saved before the iteration,
c  reset at each iteration and ONLY UPDATED after the solution
c  has converged (at least for plain strain)
c
c  Step 1: 
c
c  Make a guess at where the root is (upper and lower bound)
c  I used +/- D11. Since D22 = 0, D33 should be the shinkage
c  due to D11 so it seems reasonable to assume it is somewhere
c  between (-D11,+D11)...
c
c
              d33trial(1) = -deps(1)
              d33trial(2) = deps(1)
c
c Step 2: Find stress at guesses so that we can find the derivative
c         for the secant method
c
              dtrial = deps
              strial = sig
              hisv(1) = eps
              dtrial(9) = d33trial(1)
              call umat_43 (cm,strial,dtrial,hisv,0.0_8,0.0_8)
              s33trial(1) = strial (9)
              dtrial = deps
              strial = sig
              hisv(1) = eps
              dtrial(9) = d33trial(2)
              call umat_43 (cm,strial,dtrial,hisv,0.0_8,0.0_8)
              s33trial(2) = strial (9)
              deltas = 100
              cnt = 3
c
c Step 3: Iterate on D33 s.t. S33 = 0
c
              do while ((DABS(deltas).GT.tol*100).AND.(cnt.LT.7))
                    hisv(1) = eps
                    dtrial = deps
                    strial = sig
                    xnm1 = d33trial (cnt-1)
                    xnm2 = d33trial (cnt-2)
                    fnm1 = s33trial (cnt-1)
                    fnm2 = s33trial (cnt-2)
                    d33trial(cnt) = 
     &                      xnm1 - fnm1 * ((xnm1-xnm2)/(fnm1-fnm2))
                    dtrial(9) = d33trial (cnt)
                    call umat_43(cm,strial,dtrial,hisv,0.0_8,0.0_8)
                    s33trial(cnt) = strial (9)
c                    print "(9e15.6)",s33trial(cnt)
                    deltas = s33trial (cnt) - s33trial(cnt-1)
                    cnt = cnt + 1
              end do
c
c Step 4: Update stress and strain accordingly
c
              sig = strial
              deps = dtrial
              eps = hisv(1)
c
c  No spin in deformation so no objective update?
c  Regardless, I do not have F so I can't find W
c
              ep = ep + deps
c
              print "(20e15.6)",ep,sig,hisv(1),hisv(11)
c
         end do
c
       end program A4_725
c
