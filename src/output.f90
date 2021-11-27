!==========================================================
subroutine output( ) 
!==========================================================

use param
use globvar
implicit none     
integer:: i,d,npart,n
character(len=4)::number

n = itimestep/save_step
write(number,1000) n
1000 format(I4.4)

open(1,file="C:\Users\edwar\Documents\outputdata/f_xv"//number//".dat")
open(2,file="C:\Users\edwar\Documents\outputdata/f_state"//number//".dat")
open(3,file="C:\Users\edwar\Documents\outputdata/f_other"//number//".dat")
open(4,file="C:\Users\edwar\Documents\outputdata/f_force"//number//".dat")

if (n .eq. 1.d0) then
  open(70,file="C:\Users\edwar\Documents\outputdata/time.dat")
else
  open(70,file="C:\Users\edwar\Documents\outputdata/time.dat",status="old",position="append")
endif

write(70,1070) itimestep, time, cputime+s2-s1

do i = 1, ntotal
  write(1,1001) i, parts(i)%x(:), parts(i)%vx(:)
  write(2,1002) i, mass, parts(i)%rho, parts(i)%p, c
  write(3,1003) i, parts(i)%itype, hsml
end do

1001  format(1x, I6, 6(2x, e24.17))
1002  format(1x, I6, 7(2x, e24.17))
1003  format(1x, I6, 2x, I4, 2x, e24.17)
1004  format(1x, I6, 8(2x, e24.17))
1070  format(1x, I10, 2(2x, e24.17))
                                        
close(1)
close(2)
close(3)
close(4)
close(70)
      
end