c	Extract the dgrav-component of displacement from 'getg.in'.
	open(2,file='getg.in')
	rewind(2)
	open(8,file='getg.out')
	i=0
5	read(2,76,end=99) xc,yc,dispx1,dispy1,dispz1,exx1,eyy1,exy1,exz1,
     &	eyz1,ezz1,omxy1,dgrav1
	i=i+1
c	Write out the gravity anomaly in migrogals.
	write(8,*)  dgrav1*1000.
	go to 5
99	continue
	close(2)
	close(8)
76	format(2f10.3,11e13.6e2)
	end
