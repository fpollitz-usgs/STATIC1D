c	Extract the z-component of displacement from 'getz.in'.
	open(2,file='getz.in')
	rewind(2)
	open(8,file='getz.out')
	i=0
5	read(2,76,end=99) xc,yc,dispx1,dispy1,dispz1,exx1,eyy1,exy1,exz1,eyz1,dgrav
	i=i+1
c	Write out the uplift in cm.
	write(8,*)  dispz1
	go to 5
99	continue
	close(2)
	close(8)
76	format(2f10.3,9e13.6e2)
	end
