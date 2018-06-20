        subroutine testCH (A)
c
         double precision, intent (in) :: A(9)
         double precision :: B(9)
c
         call CH (A,B)
         print *,"----------------"
         print *,"Begin CH test"
         print *,"----------------"
         print *,"Input Tensor:"
         print *,
         call t2print (A)
         print *,
         print *,"Principal Values:"
         print *,
         call t2print (B)
         print *,
         print *,"----------------"
         print *,"End CH test"
         print *,"----------------"  
c
        end subroutine testCH
c
        program A2_725
c
         double precision :: A(9)
c
         A = (/ 4,1,6,1,0,5,6,5,9 /)
         call testCH (A)
c
        end program A2_725
