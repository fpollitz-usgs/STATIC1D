c	Write out the magnitude of horizontal motion.
	open(2,file='disph.in')
	rewind(2)
	open(4,file='disph.gmt')
	open(8,file='disph.out')
5	read(2,76,end=99) xc,yc,dispx1,dispy1,dispz1,exx1,eyy1,exy1,exz1,eyz1,ezz1
c	Mult by 10. to convert from cm to mm.
	write(4,*) dispx1,dispy1
	write(8,76) 10.*sqrt(dispx1**2+dispy1**2)
c	write(8,76) 10.*dispz1
	go to 5
99	continue
	close(2)
	close(4)
	close(8)
c	write(6,*)'end of disph'
76	format(2f10.3,9e13.6e2)
	end
