c	Add the contents of 'addw2.in1' and 'addw2.in2' together.
c	These are assumed to be in the format of output from 'strainw.f'
	open(2,file='addw2.in1')
	rewind(2)
	open(4,file='addw2.in2')
	rewind(4)
	open(8,file='addw2.out')
5	read(2,76,end=99) xc,yc,dispx1,dispy1,dispz1,exx1,eyy1,exy1,exz1,
     &	eyz1,ezz1,omxy1,dgrav1
76	format(2f10.3,11e13.6e2)
	read(4,76) xc,yc,dispx2,dispy2,dispz2,exx2,eyy2,exy2,exz2,
     &	eyz2,ezz2,omxy2,dgrav2
	dispx=dispx1+dispx2
	dispy=dispy1+dispy2
	dispz=dispz1+dispz2
	exx=exx1+exx2
	eyy=eyy1+eyy2
	exy=exy1+exy2
	exz=exz1+exz2
	eyz=eyz1+eyz2
	ezz=ezz1+ezz2
	omxy=omxy1+omxy2
	dgrav=dgrav1+dgrav2
	write(8,76) xc,yc,dispx,dispy,dispz,exx,eyy,exy,exz,eyz,ezz,omxy,dgrav  
	go to 5
99	continue
	close(2)
	close(4)
	close(8)
c	write(6,*)'end of addw2'
	end
