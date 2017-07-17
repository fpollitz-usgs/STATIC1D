c	Multiply the contents of 'multw2.in' by the value in 'multw.val'
c	These are assumed to be in the format of output from 'strainw.f'
	open(2,file='multw2.in')
	rewind(2)
	read(5,*) val
	open(8,file='multw2.out')
5	read(2,76,end=99) xc,yc,dispx1,dispy1,dispz1,exx1,eyy1,exy1,exz1,eyz1,
     &	ezz1,omxy1,dgrav1
76	format(2f10.3,11e13.6e2)
	dispx=dispx1*val
	dispy=dispy1*val
	dispz=dispz1*val
	exx=exx1*val
	eyy=eyy1*val
	exy=exy1*val
	exz=exz1*val
	eyz=eyz1*val
	ezz=ezz1*val
	omxy=omxy1*val
	dgrav=dgrav1*val
	write(8,76) xc,yc,dispx,dispy,dispz,exx,eyy,exy,exz,eyz,ezz,omxy,dgrav
	go to 5
99	continue
	close(2)
	close(4)
	close(8)
c	write(6,*)'end of multw2'
	end
