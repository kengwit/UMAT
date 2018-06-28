      subroutine t2print (A)
c
       implicit none
       double precision, intent (in) :: A(9)
       integer i,base
       do i=0,2
        base = 3*i
        print '(3e12.4)',A(base+1),A(base+2),A(base+3)
       end do
c
      end subroutine t2print

      subroutine t4print (A)
c
       implicit none
       double precision, intent (in) :: A(9,9)
       integer i
       do i=1,9
        print '(9e10.2)',A(i,1),A(i,2),A(i,3)
     &                 ,A(i,4),A(i,5),A(i,6)
     &                 ,A(i,7),A(i,8),A(i,9)
       end do
c
      end subroutine t4print
c
      subroutine PI (p)
c
       implicit none
       double precision, intent (out) :: p
c
       p = 4.0_8*atan(1.0_8)
c
      end subroutine pi
c
      double precision function ABSPOW (a,b)
c
       implicit none
       double precision, intent (in) :: a,b
c
       ABSPOW = ABS(a) ** b
c
       return
c
      end function ABSPOW
c
      double precision function ABSQRT (a)
c
       implicit none
       double precision, intent (in) :: a 
c
       ABSQRT = SQRT(ABS(a))
       print *,a
       print *,
       print *,ABSQRT
c
      end function ABSQRT
c
