       program A4_725
c
         implicit none
c
c        parameters of interest: stress, strain, eq. plastic strain
c        THESE ARE OUTPUTS
c
         double precision :: sig(9), ep (9), eps
         double precision :: dtrial (9), strial (9)
         double precision :: d22trial (7), s22trial (7)
         double precision :: d33trial (7), s33trial (7)
         double precision :: deps (9), cm (5), hisv (11)
         double precision :: x1nm1, x1nm2, f1nm1, deltas22
         double precision :: x2nm1, x2nm2, f2nm1, deltas33
         double precision :: f1x1m2, f1x2m2, f2x1m2, f2x2m2
         double precision :: df1dx1, df1dx2, df2dx1, df2dx2
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
         d22trial = 0.0_8
         s22trial = 0.0_8
         d33trial = 0.0_8
         s33trial = 0.0_8
         cnt = 1
         final = 2
c
         do i = 1,final 
c
              deps = 0.0_8
              deps (1) = 0.001_8
              eps = hisv(1)
c
c   Step 1: Guess bounds of d22, d33
c
              d22trial(1) = -deps(1)
              d22trial(2) = deps(1)
              d33trial(1) = -deps(1)
              d33trial(2) = deps(1)
c
c   Step 2: Find cooresponding s22,s33
c
              dtrial = deps
              strial = sig
              hisv(1) = eps
              dtrial(5) = d22trial(1)
              dtrial(9) = d33trial(1)
              call umat_43(cm,strial,dtrial,hisv,0.0_8,0.0_8)
              s22trial(1) = strial(5)
              s33trial(1) = strial(9)             
c
              dtrial = deps
              strial = sig
              hisv(1) = eps
              dtrial(5) = d22trial(2)
              dtrial(9) = d33trial(2)
              call umat_43(cm,strial,dtrial,hisv,0.0_8,0.0_8)
              s22trial(2) = strial(5)
              s33trial(2) = strial(9)
              deltas22 = 100
              deltas33 = 100
              cnt = 3             
c             
              do while (    (DABS(deltas22).GT.tol*100)
     &                 .AND.(DABS(deltas33).GT.tol*100)
     &                 .AND.(cnt.LT.200)   )
c
                  dtrial = deps
                  strial = sig
                  hisv(1) = eps
c
                  x1nm1 = d22trial(cnt-1)
                  x1nm2 = d22trial(cnt-2)
                  x2nm1 = d33trial(cnt-1)
                  x2nm2 = d33trial(cnt-2)
c
                  f1nm1 = s22trial(cnt-1)
                  f2nm1 = s33trial(cnt-1)
c
                  dtrial = deps
                  strial = sig
                  hisv(1) = eps
                  dtrial(5) = d22trial(cnt-2)
                  dtrial(9) = d33trial(cnt-1)
                  call umat_43(cm,strial,dtrial,hisv,0.0_8,0.0_8)
                  f1x1m2 = strial(5)
                  f2x1m2 = strial(9)             
c
                  dtrial = deps
                  strial = sig
                  hisv(1) = eps
                  dtrial(5) = d22trial(cnt-1)
                  dtrial(9) = d33trial(cnt-2)
                  call umat_43(cm,strial,dtrial,hisv,0.0_8,0.0_8)
                  f1x2m2 = strial(5)
                  f2x2m2 = strial(9)             
c
                  df1dx1 = (f1nm1-f1x1m2)/(x1nm1-x1nm2)
                  df1dx2 = (f1nm1-f1x2m2)/(x2nm1-x2nm2)
                  df2dx1 = (f2nm1-f2x1m2)/(x1nm1-x1nm2)
                  df2dx2 = (f2nm1-f2x2m2)/(x2nm1-x2nm2)
c
                  d22trial(cnt) = x1nm1
     &                          - f1nm1 / df1dx1 - f2nm1 / df2dx1
c
                  d33trial(cnt) = x2nm1 
     &                          - f1nm1 / df1dx2 - f2nm1 / df2dx2

c
                  hisv(1) = eps
                  dtrial = deps
                  strial = sig
                  dtrial(5) = d22trial(cnt)
                  dtrial(9) = d33trial(cnt)
                  call umat_43(cm,strial,dtrial,hisv,0.0_8,0.0_8)
                  s22trial(cnt) = strial(5)
                  s33trial(cnt) = strial(9)
                  deltas22 = s22trial(cnt) - s22trial(cnt-1)
                  deltas33 = s33trial(cnt) - s33trial(cnt-1)
                  print "(2e15.6)",d22trial(cnt-1),d22trial(cnt)
c                  print "(4e15.6)",f1nm1,df1dx2,f2nm1,df2dx2
                  cnt = cnt + 1                 
c
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
               print "(9e15.6)",sig
c              print "(20e15.6)",ep,sig,hisv(1),hisv(11)
c
         end do
c
       end program A4_725
c
