module globvar
	
	use derived_types
	implicit none
	
	! particle array
	type(particles),allocatable,target:: parts(:)
	
	! interaction array
	type(interactions),allocatable:: pairs(:)
	
	!global 1D variables
	integer:: ntotal,nvirt
	integer:: mp,np,op,pp,maxn,maxinter,maxnsend
	integer:: niac
	integer:: nstart,current_ts,maxtimestep,itimestep,yesorno
	
	real(8):: scale_k
	
	!diagnostics
	integer:: nvalid
	
	!timing
	real(8):: cputime,s1,s2,time,testtime(20)
	
end module globvar