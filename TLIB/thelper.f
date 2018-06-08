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
