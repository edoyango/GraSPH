subroutine flink_list( )
! subroutine to search for particle interactions using cell-linked list.

use param
implicit none
type particleincellarray
  type(particles),pointer:: p
end type particleincellarray

integer,parameter:: maxpcell=12
integer:: ngridx(dim),i,j,k,d,icell,jcell,xi,yi
real(8):: maxgridx(dim),mingridx(dim),dcell
integer,allocatable:: pincell(:,:)
type(particleincellarray),allocatable:: cells(:,:,:)
integer,parameter:: sweep(2,4) = reshape((/ 1,-1,&
                                            1, 0,&
							                0, 1,&
							                1, 1/),(/2,4/))

!Generate grid geometry and number of grids in x, y directions
mingridx(:) = parts(1)%x(:)
maxgridx(:) = parts(1)%x(:)
do i = 2,ntotal+nvirt
  do d = 1,dim
    if (parts(i)%x(d)<mingridx(d)) mingridx(d) = parts(i)%x(d)
	if (parts(i)%x(d)>maxgridx(d)) maxgridx(d) = parts(i)%x(d)
  end do
end do

dcell = scale_k*hsml
maxgridx(:) = maxgridx(:) + 2d0*dcell
mingridx(:) = mingridx(:) - 2d0*dcell
ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
maxgridx(:) = mingridx(:) + ngridx(:)*dcell

allocate( pincell(ngridx(1),ngridx(2)),&
			cells(maxpcell,ngridx(1),ngridx(2)) )
			
!Mapping particles to grid cells
pincell(:,:) = 0

do i=1,ntotal+nvirt
  icell = int((parts(i)%x(1) - mingridx(1))/dcell) + 1
  jcell = int((parts(i)%x(2) - mingridx(2))/dcell) + 1 
  pincell(icell,jcell) = pincell(icell,jcell) + 1
  cells(pincell(icell,jcell),icell,jcell)%p => parts(i)
enddo

! Searching for interactions. Algorithm loops through every cell with indices icell,jcell.
! First checks particles within cell are interacting, then checks particles interacting with
! particles in adjacent cells. 
niac = 0
do icell = 2,ngridx(1)-2
	do jcell = 2,ngridx(2)-2
	
		! finding pairs within cell icell,jcell
		do i = 1,pincell(icell,jcell)-1
			do j = i+1,pincell(icell,jcell)
				call check_if_interact(cells(i,icell,jcell)%p,cells(j,icell,jcell)%p)
			end do
		end do
		
		! finding pairs between particles in cell icell,jcell and particles in cell xi,yi
		do k = 1,4
			xi = icell + sweep(1,k)
			yi = jcell + sweep(2,k)
			do i = 1,pincell(icell,jcell)
				do j = 1,pincell(xi,yi)
					call check_if_interact(cells(i,icell,jcell)%p,cells(j,xi,yi)%p)
				end do
			end do
		end do
		
	end do
end do

!Statistics for the interaction, Print information to screen
if (mod(itimestep,print_step)==0) then
	call print_interact_info
endif    

end subroutine flink_list

!--------------------------------------------------------------------------------------------------
subroutine check_if_interact(p_i,p_j)
! subroutine to chekc if two particles are interacting and consequently adding to pair list
use param

implicit none
type(particles),intent(in),target:: p_i,p_j
real(8):: dxiac(dim),r
dxiac(:) = p_i%x(:) - p_j%x(:)
r = SQRT(SUM(dxiac*dxiac))
if (r < hsml*scale_k) then
	niac = niac + 1
	if (niac < maxinter) then
		pairs(niac)%i => p_i
		pairs(niac)%j => p_j
		call kernel(r,dxiac,hsml,pairs(niac)%w,pairs(niac)%dwdx(:))
	else
		print *,' >>> Error <<< : Too many interactions'
		stop
	end if
end if
end subroutine check_if_interact

!--------------------------------------------------------------------------------------------------
subroutine print_interact_info
! subroutine to print interaction summary data
use param

implicit none
integer,allocatable:: ns(:)
integer:: maxp,minp,sumiac,maxiac,miniac,noiac,i,j,k,d

allocate(ns(ntotal+nvirt))
ns(1:ntotal+nvirt) = 0
do k = 1,niac
	ns(pairs(k)%i%ind) = ns(pairs(k)%i%ind) + 1
	ns(pairs(k)%j%ind) = ns(pairs(k)%j%ind) + 1
end do
   sumiac = 0
   maxiac = 0
   noiac  = 0
   miniac = 1000
   do i=1,ntotal+nvirt
     sumiac = sumiac + ns(parts(i)%ind)
     if (ns(parts(i)%ind) > maxiac) then
       maxiac = ns(parts(i)%ind)
   	maxp = i
     endif
     if (ns(parts(i)%ind) < miniac) then 
       miniac = ns(parts(i)%ind)
       minp = i
     endif
     if (ns(parts(i)%ind).eq.0) noiac = noiac + 1
   enddo
   print *,' >> Statistics: interactions per particle:'
   print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
   print *,'**** Particle:',minp, ' minimal interactions:',miniac
   print *,'**** Average :',real(sumiac)/real(ntotal)
   print *,'**** Total pairs : ',niac
   print *,'**** Particles with no interactions:',noiac
deallocate(ns)
	
end subroutine print_interact_info