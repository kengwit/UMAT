c
        program A2_725
c
         double precision :: A(9)
         double precision :: f,r
         f = 1.0e-8_8
         r = ABSPOW(-4.0_8,0.5_8)
         print *,r
c
c        ZERO CHECK
c
         A = 0.0_8
         call testCH(A)
c
c        ALREADY IN PRINCIPAL SPACE
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 1.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 1.0_8 /)
         call testCH(A)
c
c       ONE ZERO PRINCIPAL VALUE
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 1.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8 /)
         call testCH(A)
c
c       TWO ZERO PRINCIPAL VALUES
c
         A = (/ 1.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8,
     &          0.0_8, 0.0_8, 0.0_8 /)
         call testCH(A)
c
c       TOLERANCE BUFFET
c
         A = (/ 1.0_8+f,           0.0_8, 0.0_8,
     &            0.0_8, 1.0_8/(1.0_8+f), 0.0_8,
     &            0.0_8,           0.0_8,     f /)
         call testCH(A)
c
c       ALL TOLERANCE
c
         A = (/    f, 0.0_8, 0.0_8,
     &         0.0_8,     f, 0.0_8,
     &         0.0_8, 0.0_8,     f /)
         call testCH(A)
c
c       RANDOM VALUES
c
         A = (/ 300.0_8,  -20.0_8,   25.0_8,
     &          -20.0_8,  100.0_8, -125.0_8,
     &           25.0_8, -125.0_8,   50.0_8 /)
         call testPD(A)
c
        end program A2_725
