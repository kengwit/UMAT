c
c       simple Swift's hardening law      
c
c       
        subroutine iso_hard (K,eo,n,ep,Sb,h)
c
         double precision, intent (in) :: K,eo,n,ep
         double precision, intent (out) :: Sb,h
c
         Sb = K * (eo+ep) ** n
         h = K*n*(eo+ep)**(n-1.0_8)
c
         return
c
        end subroutine iso_hard
