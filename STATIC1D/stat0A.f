	program STAT0
c***
c	Version of STAT0 which calculates Green's functions for
c	[ndep] different depths and writes them all onto STAT0.OUT.
c	Thus [ndep] times more space is required as for STAT0B.
c	Anticipating degree numbers
c	as high as 25000 use a step of [lste] in the l-loop
	parameter (ndep=33)
c***
c	Find weight coefficients for static displacement field
c	due to buried seismic source 
c	at origin of spherical coordinate
c	system.  m=0, 1, 2 harmonics contribute.
c	Do depths dmin thru dmax km at intervals of (dmax-dmin)/(ndep-1) km.
c	Find coefficients which control
c	geometric dependence of displacement fields.
c	Input file   
c	'earth.model' has format
c	N  Nd 	(Number of layers; Number of different layers)
c	rb(1)	     |  rt(1)   |dens(1)|  kappa(1)  | mu(1)  | eta(1)
c	...		...	 ...	   ...	       ...      ...
c	rb(N)	     |  rt(N)   |dens(N)|  kappa(N)  | mu(N)  | eta(N)
c	bottom radius|top radius|density|bulk modulus|rigidity|viscosity
c	(km)          (km)      (g/cm**3) (10**11 dyne/cm**2)  10**19
c							       dyne-s/cm**2
c		                                               10**18 Pa-s
	character*80 aread
	real*4 kappa,mu,dens(200)
	real*8 ub,ua,sb,sa,vb,rx,va,ry
	real*8 z11,z12,z13,z21,z23
	real*8 aj(4,4),d1(4),d2(4),d3(4),d4(4),d5(4),d6(4)
	real*8 ajt(2,2),dt1(2),dt2(2),dt3(2)
	real*8 x11,x12,x13,x14,x15,x16
	real*8 x21,x22,x23,x24,x25,x26
	real*8 us1,vs1,us2,vs2,up1,vp1
	real*8 up2,vp2,us3,vs3,us4,vs4,up3,vp3,up4,vp4
	real*8 us5,vs5,us6,vs6,up5,vp5,up6,vp6
	real*8 dur3,dur4,dvr3,dvr4,ur3,ur4,vr3,vr4
	real*8 dur5,dur6,dvr5,dvr6,ur5,ur6,vr5,vr6
	real*8 ur1,ur2,vr1,vr2,dur1,dur2,dvr1,dvr2
	real*8 wr1,dwr1,wr2,dwr2,wr3,dwr3
	real*8 r,dr,r2,biga,fl21,flp1,odep
	real*8 y7,y8,ddens,r0
	real*8 lams,mus
	real*8 c(4),b1(4),b2(4),b3(4),b4(4),b5(4),b6(4)
	real*8 b1s(4),b2s(4),b3s(4),b4s(4),b5s(4),b6s(4)
	real*8 ct(2),det,t1(2),t2(2),t3(2),t1s(2),t2s(2),t3s(2)
	real*8 am1(4,4)
	real*8 ylim(3)
	real*8 a00(4,4),bc0(6,6),c0(6),a0(6)
	real*8 bc1(3,3),c1(3),a1(3)
	dimension eta(200)
	common/mat/rb(200),rt(200),kappa(200),mu(200),l,r0
	pi=3.14159265
	twopi=2.*pi
	write(6,*)'lmin,lmax?'
	read(5,*) lmin,lmax
c *
		lste=1
c *
	flste=real(lste)
c	Round up lmin,lmax to the nearest (corresponding) multiple of [lste].
	lmin=lste*((lmin-1)/lste)+lste
	lmax=lste*((lmax-1)/lste)+lste
	write(6,*)'Max, Min depth [km] of Greens functions?'
	read(5,*) dmax,dmin
	write(6,*)'observation depth (km)?'
	read(5,*) odepr
	open(2,file='stat0.out',form='unformatted')
21	format(i2)
c 
cc	source depth=-0.5+1.0*ide km.
c	source depth=dmin+(dmax-dmin)*(ide-1)/(ndep-1) km
	ide=0
44	ide=ide+1
c		write(0,*)'ide=',ide
	if(ide.gt.ndep) go to 94
c ** *	
c	Read in earth model
	open(4,file='earth.model_stat')
	rewind(4)
	read(4,5) n,nd,bigr
	write(6,5) n,nd,bigr
5	format(2i2,f10.3)
c	bigr=earth radius in km
c	dep=1.0-0.717*real(ide)/bigr
c	dep=1.0-(-0.5+real(ide))/bigr
	dep=1.0-(dmin+(dmax-dmin)*real(ide-1)/real(ndep-1))/bigr
	  write(6,*)'dep=',dep
	if(ide.eq.1) write(2) bigr,dmax,dmin
	if(ide.eq.1) write(2) lmin,lmax
	if(ide.eq.1) odep=dble(1.-odepr/bigr)
	if(ide.eq.1) write(2) odep
	do 10 j=1,n
	i=2*j-1
c	following designed to read 'earth.model' files in either the old
c	or new format (which may include mu' for every layer).
	read(4,13) aread 
	read(aread,17,err=11) rb(i),rt(i),dens(i),kappa(i),mu(i),adum,eta(i)
	go to 22
11	read(aread,15) rb(i),rt(i),dens(i),kappa(i),mu(i),eta(i)
	adum=0.
22	continue
	write(6,17) rb(i),rt(i),dens(i),kappa(i),mu(i),adum,eta(i)
	dens(i+1)=dens(i)
c		write(6,*)'dens(',i,')=',dens(i)
c		write(6,*)'dens(',i+1,')=',dens(i+1)
	mu(i+1)=mu(i)
	eta(i+1)=eta(i)
	kappa(i+1)=kappa(i)
	rb(i+1)=(rb(i)+rt(i))/2.
	rt(i+1)=rt(i)
	rt(i)=rb(i+1)
	rb(i+1)=rb(i+1)/bigr + (1.-6371./bigr)
	rt(i+1)=rt(i+1)/bigr + (1.-6371./bigr)
	rb(i)=rb(i)/bigr + (1.-6371./bigr)
	rt(i)=rt(i)/bigr + (1.-6371./bigr)
10	continue
13	format(a80)
15	format(5f9.3,e13.6e2)
17	format(6f9.3,e13.6e2)
	n=2*n 
	close(4)
c	Put source at boundary between layer [isdep-1] and [isdep].
c	This involves putting in an extra layer within a pre-existing
c	layer.
	if(dep.lt.rb(n)) go to 12
	rt(n)=dep
	rb(n+1)=dep
	rt(n+1)=1.0
	dens(n+1)=dens(n)
	eta(n+1)=eta(n)
	kappa(n+1)=kappa(n)
	mu(n+1)=mu(n)
	isdep=n+1
	go to 18
12	isdep=0
	do 16 j=2,n
	i=n-j+2
	if(rb(i).lt.dep.or.rb(i-1).gt.dep) go to 14
	write(6,*)'put in source layer at depth',bigr*(1.-dep),'km'
	write(6,*)'rb(',i,')=',rb(i),'rt(',i,')=',rt(i)
	isdep=i
	rt(i-1)=dep
	rb(i+1)=rb(i)
	rt(i+1)=rt(i)
	rt(i)=rb(i)
	rb(i)=dep
	eta(i+1)=eta(i)
	kappa(i+1)=kappa(i)
	mu(i+1)=mu(i)
	eta(i)=eta(i-1)
	kappa(i)=kappa(i-1)
	mu(i)=mu(i-1)
14	if(isdep.gt.0) go to 16
	rt(i+1)=rt(i)
	rb(i+1)=rb(i)
	eta(i+1)=eta(i)
	kappa(i+1)=kappa(i)
	mu(i+1)=mu(i)
16	continue
18	continue
	n=n+1
		write(6,*)'number of layers=',n
cc	check the model for correct interpolation.
c	write(6,*)'BBB'
c	do j=1,n
c	write(6,*) rb(j),rt(j),dens(j),kappa(j),mu(j),eta(j)
c	enddo
	write(2) kappa(isdep),mu(isdep)
	nrad=10
c ** *
	do 50 l=lmin,lmax,lste
	write(6,*)'l=',l
	depmax=dep-2.5*(2.*pi/(real(l)+0.5))
	depmin=dep+2.5*(2.*pi/(real(l)+0.5))
		write(6,*)'depmax,depmin,odep=',depmax,depmin,odep
	  awr=(1.-real(depmax))*bigr
	  write(6,*)'max integration depth=',awr
c	depmax=dimensionless radius 2 wavelengths into the earth.
	x11=0.d0
	x12=0.d0
	x13=0.d0
	x14=0.d0
	x15=0.d0
	x16=0.d0
	x21=0.d0
	x22=0.d0
	x23=0.d0
	x24=0.d0
	x25=0.d0
	x26=0.d0
	fl21=dble(l*(l+1))
	flp1=dble(l+1)
c	SPHEROIDAL MODES
c	Start matrix propagation upward for solutions b1,b2.  Stop
c	at layer [isdep-1].
c	Give values at CMB those appropriate for homogeneous solid
c	in 0<r<rb(1).
	is=0
	do 30 i=1,isdep-1
	if(i.lt.2) go to 30
c	Start integration at most 2.5 wavelengths into the earth.
	if(rt(i).lt.depmax.or.rb(i).gt.depmin) go to 30
	if(dble(depmax).ge.odep.or.dble(depmin).le.odep) go to 30
	mus=dble(mu(i))
	lams=dble(kappa(i)-2.*mu(i)/3.)
	biga=lams+2.d0*mus
	r0=dble(rb(i))
c	Give boundary conditions at top of layer (i-1) appropriate for
c	homogeneous solid in 0 < r < rt(i-1).
	if(is.ne.0) go to 28
	is=1
	j=i-1
	call matra1(aj,j,1)
	do 27 k=1,4
	b1(k)=aj(k,1)
	b2(k)=aj(k,2)
27	continue
c--	If following lines are not commented out, they
c	give boundary conditions appropriate for an incompressible
c	fluid in the sphere 0 < r < rt(i-1).
	b1(1)=dble(l)
	b1(2)=0.d0
	b1(3)=1.d0
	b1(4)=0.d0
	b2(1)=0.d0
	b2(2)=0.d0
	b2(3)=1.d0
	b2(4)=0.d0
c--	
28	continue
c	Determine coefficient matrices d1 and d2 for layer i.
	do 33 m=1,2
	if(m.eq.1) call equal(b1,c,4)
	if(m.eq.2) call equal(b2,c,4)
	call matra1(aj,i,0)
	if(m.eq.1) call ainver(d1,aj,c,4)
	if(m.eq.2) call ainver(d2,aj,c,4)
33	continue
29	dr=dble(rt(i)-rb(i))/dble(nrad)
c	Integrate equations of motion to propagate solutions up to
c	the surface.
	r=dble(rb(i))-dr/2.d0
	do 31 j=1,nrad
	r2=r*r 
	r=r+dr
	r0=dble(r)  
	j0=j0+1
	call matra1(aj,i,2)
	do 32 m=1,2
	if(m.eq.1) call prodr(aj,d1,4,c)
	if(m.eq.2) call prodr(aj,d2,4,c)
	if(m.eq.1) call equal(c,b1,4)
	if(m.eq.2) call equal(c,b2,4)
32	continue
c	Calculate the contribution of this depth to the gravity anomaly.
	us1=b1(1)/r
	vs1=b1(3)/r
	us2=b2(1)/r
	vs2=b2(3)/r
	up1=(b1(2)-(lams/r)*(2.d0*b1(1)-fl21*b1(3)))/biga 
	up2=(b2(2)-(lams/r)*(2.d0*b2(1)-fl21*b2(3)))/biga 
c	up is equal to dy1/dr.
	vp1=(1.d0/r)*(b1(3)-b1(1))+b1(4)/mus 
	vp2=(1.d0/r)*(b2(3)-b2(1))+b2(4)/mus
c	vp is equal to dy3/dr.
	if(r.ge.odep) go to 19
	ur1=us1
	ur2=us2
	vr1=vs1
	vr2=vs2
	dur1=up1
	dur2=up2
	dvr1=vp1
	dvr2=vp2
c	Determine contribution to gravity anomaly coeff.
19	y7=flp1*dlog(r)
	y8=0.d0
	if(y7.gt.-20.d0) y8=dexp(y7)
	x11=x11-(2.d0*us1+up1-fl21*vs1)*
     &	y8*r*dble(dens(i))*dr
	x12=x12-(2.d0*us2+up2-fl21*vs2)*
     &	y8*r*dble(dens(i))*dr
c	x11, x12 represent the gravity 
c	anomaly due to intrinsic density
c	changes following faulting.
c	At top of layer, evaluate d(dens)/dr contribution.
	if(j.ne.nrad) go to 31
	if(i.ne.n) ddens=dble(dens(i+1)-dens(i))
	if(i.eq.n) ddens=1.000d0-dble(dens(i))
c	Use water density of 1.000 g-cm^3 in above line
cBOUGUER	  if(i.eq.n) go to 31
c	  Previous line: want Bouguer anomaly.
	x21=x21-ddens*b1(1)*y8*r
	x22=x22-ddens*b2(1)*y8*r
c	x21, x22 represent the gravity 
c	anomaly due to vertical motions
c	of a volume of stratified material.
c*
c	If Bouguer anomaly is desired, then omit the x21,x22
c	contribution from the i=n (free surface) term above.
c	i.e., uncomment the cBOUGUER line above.
c*
31	continue
c		write(6,*)'TOP layer',i,'b1=',b1,'b2=',b2
c	Make sure that matrix propagation  is done accurately.
	do 53 m=1,2
	if(m.eq.1) call equal(b1,c,4)
	if(m.eq.2) call equal(b2,c,4)
	call matra1(aj,i,1)
	if(m.eq.1) call prodr(aj,d1,4,c)
	if(m.eq.2) call prodr(aj,d2,4,c)
	if(m.eq.1) call equal(c,b1,4)
	if(m.eq.2) call equal(c,b2,4)
53	continue
46	bmag=abs(real(b1(1)))+real(l)*abs(real(b1(3)))
	if(bmag.lt.1.e+10) go to 30
c	  write(6,*)'bmag=',bmag
	do 36 j1=1,4
	b1(j1)=b1(j1)*1.d-10
36	b2(j1)=b2(j1)*1.d-10
	ur1=ur1*1.d-10
	ur2=ur2*1.d-10
	vr1=vr1*1.d-10
	vr2=vr2*1.d-10
	dur1=dur1*1.d-10
	dur2=dur2*1.d-10
	dvr1=dvr1*1.d-10
	dvr2=dvr2*1.d-10
	x11=x11*1.d-10
	x12=x12*1.d-10
	x21=x21*1.d-10
	x22=x22*1.d-10
30	continue
	if(is.eq.1) go to 26
c	Write out zero displacements if source is too deep.
	y1=0.
	y2=0.
	y3=0.
	y4=0.
	bfacl=0.
	bfacm=0.
	write(2) l,y1,y2,y3,y4
	write(2) bfacl,bfacm
	write(2) l,y1,y2,y3,y4
	write(2) bfacl,bfacm
	write(2) l,y1,y2,y3,y4
	write(2) bfacl,bfacm
	write(2) l,y1,y2,y3,y4
	write(2) bfacl,bfacm
	write(2) l,y1,y2,y3,y4
	write(2) bfacl,bfacm
	write(2) l,y1,y2
	write(2) l,y1,y2
	go to 50
26	continue
c	Now add in contribution to gravity anomaly from source radius.
	x11=x11-dble(dens(isdep-1))*(-b1(1))*y8*rt(isdep-1)
	x12=x12-dble(dens(isdep-1))*(-b2(1))*y8*rt(isdep-1)
c	Store b1 - b2 at source radius.
	do 40 j1=1,4
	b1s(j1)=b1(j1)
	b2s(j1)=b2(j1)
40	continue
c	  write(6,*)'finished b1,b2'
c *
c	Start matrix propagation upward for solutions b3,b4,b5,b6.  Stop
c	at layer [n].
c	Give boundary conditions at top of layer [isdep-1] appropriate for
c	homogeneous solid in this layer.
	j=isdep-1
	call matra1(aj,j,1)
	do 127 k=1,4
	b3(k)=aj(k,1)
	b4(k)=aj(k,2)
	b5(k)=aj(k,3)
	b6(k)=aj(k,4)
c	Now add in contribution to gravity anomaly from source radius.
	x13=-dble(dens(isdep-1))*b3(1)*y8*rt(isdep-1)
	x14=-dble(dens(isdep-1))*b4(1)*y8*rt(isdep-1)
	x15=-dble(dens(isdep-1))*b5(1)*y8*rt(isdep-1)
	x16=-dble(dens(isdep-1))*b6(1)*y8*rt(isdep-1)
c	Store these solutions at source radius.
	b3s(k)=b3(k)
	b4s(k)=b4(k)
	b5s(k)=b5(k)
	b6s(k)=b6(k)
127	continue
	do 130 i=isdep,n
	if(rb(i).gt.depmin) go to 130
	mus=dble(mu(i))
	lams=dble(kappa(i)-2.*mu(i)/3.)
	biga=lams+2.d0*mus
	r0=dble(rb(i))
c	Determine coefficient matrices d3 - d6 for layer i.
	do 133 m=1,4
	if(m.eq.1) call equal(b3,c,4)
	if(m.eq.2) call equal(b4,c,4)
	if(m.eq.3) call equal(b5,c,4)
	if(m.eq.4) call equal(b6,c,4)
	call matra1(aj,i,0)
	if(m.eq.1) call ainver(d3,aj,c,4)
	if(m.eq.2) call ainver(d4,aj,c,4)
	if(m.eq.3) call ainver(d5,aj,c,4)
	if(m.eq.4) call ainver(d6,aj,c,4)
133	continue
	dr=dble(rt(i)-rb(i))/dble(nrad)
c	Integrate equations of motion to propagate solutions up to
c	the surface.
	r=dble(rb(i))-dr/2.d0
	do 131 j=1,nrad
	r2=r*r 
	r=r+dr
	r0=dble(r)  
	j0=j0+1
	call matra1(aj,i,2)
	do 132 m=1,4
	if(m.eq.1) call prodr(aj,d3,4,c)
	if(m.eq.2) call prodr(aj,d4,4,c)
	if(m.eq.3) call prodr(aj,d5,4,c)
	if(m.eq.4) call prodr(aj,d6,4,c)
	if(m.eq.1) call equal(c,b3,4)
	if(m.eq.2) call equal(c,b4,4)
	if(m.eq.3) call equal(c,b5,4)
	if(m.eq.4) call equal(c,b6,4)
132	continue
c	Determine contribution to gravity anomaly coeff.
	us3=b3(1)/r
	vs3=b3(3)/r
	us4=b4(1)/r
	vs4=b4(3)/r
	us5=b5(1)/r
	vs5=b5(3)/r
	us6=b6(1)/r
	vs6=b6(3)/r
	up3=(b3(2)-(lams/r)*(2.d0*b3(1)-fl21*b3(3)))/biga 
	up4=(b4(2)-(lams/r)*(2.d0*b4(1)-fl21*b4(3)))/biga 
	up5=(b5(2)-(lams/r)*(2.d0*b5(1)-fl21*b5(3)))/biga 
	up6=(b6(2)-(lams/r)*(2.d0*b6(1)-fl21*b6(3)))/biga 
c	up is equal to dy1/dr.
	vp3=(1.d0/r)*(b3(3)-b3(1))+b3(4)/mus 
	vp4=(1.d0/r)*(b4(3)-b4(1))+b4(4)/mus
	vp5=(1.d0/r)*(b5(3)-b5(1))+b5(4)/mus 
	vp6=(1.d0/r)*(b6(3)-b6(1))+b6(4)/mus
c	vp is equal to dy3/dr.
	if(r.ge.odep) go to 114
	ur3=us3
	ur4=us4
	ur5=us5
	ur6=us6
	vr3=vs3
	vr4=vs4
	vr5=vs5
	vr6=vs6
	dur3=up3
	dur4=up4
	dur5=up5
	dur6=up6
	dvr3=vp3
	dvr4=vp4
	dvr5=vp5
	dvr6=vp6
c	Determine contribution to gravity anomaly coeff.
114	y7=flp1*dlog(r)
	y8=0.d0
	if(y7.gt.-20.d0) y8=dexp(y7)
	x13=x13-(2.d0*us3+up3-fl21*vs3)*
     &	y8*r*dble(dens(i))*dr
	x14=x14-(2.d0*us4+up4-fl21*vs4)*
     &	y8*r*dble(dens(i))*dr
	x15=x15-(2.d0*us5+up5-fl21*vs5)*
     &	y8*r*dble(dens(i))*dr
	x16=x16-(2.d0*us6+up6-fl21*vs6)*
     &	y8*r*dble(dens(i))*dr
c	x15, x16 represent the gravity 
c	anomaly due to intrinsic density
c	changes following faulting.
c	At top of layer, evaluate d(dens)/dr contribution.
	if(j.ne.nrad) go to 131
	if(i.ne.n) ddens=dble(dens(i+1)-dens(i))
	if(i.eq.n) ddens=1.000d0-dble(dens(i))
c	Use water density of 1.000 g-cm^3 in above line
c	  if(i.eq.n) go to 131
c	  Previous line: want Bouger anomaly.
	x23=x23-ddens*b3(1)*y8*r
	x24=x24-ddens*b4(1)*y8*r
	x25=x25-ddens*b5(1)*y8*r
	x26=x26-ddens*b6(1)*y8*r
c	x26, x26 represent the gravity 
c	anomaly due to vertical motions
c	of a volume of stratified material.
c	If Bouger anomaly is desired, then omit the x23,x24,x25,x26
c	contribution from the i=n (free surface) term above.
131	continue
c		write(6,*)'TOP layer',i,'b3=',b3,'b4=',b4
c	Make sure that matrix propagation  is done accurately.
	do 153 m=1,4
	if(m.eq.1) call equal(b3,c,4)
	if(m.eq.2) call equal(b4,c,4)
	if(m.eq.3) call equal(b5,c,4)
	if(m.eq.4) call equal(b6,c,4)
	call matra1(aj,i,1)
	if(m.eq.1) call prodr(aj,d3,4,c)
	if(m.eq.2) call prodr(aj,d4,4,c)
	if(m.eq.3) call prodr(aj,d5,4,c)
	if(m.eq.4) call prodr(aj,d6,4,c)
	if(m.eq.1) call equal(c,b3,4)
	if(m.eq.2) call equal(c,b4,4)
	if(m.eq.3) call equal(c,b5,4)
	if(m.eq.4) call equal(c,b6,4)
153	continue
	bmag=abs(real(b3(1)))+real(l)*abs(real(b3(3)))
	if(bmag.lt.1.e+10) go to 130
	  write(6,*)'bmag=',bmag
	do 136 j1=1,4
	b3(j1)=b3(j1)*1.d-10
	b4(j1)=b4(j1)*1.d-10
	b5(j1)=b5(j1)*1.d-10
	b6(j1)=b6(j1)*1.d-10
	b3s(j1)=b3s(j1)*1.d-10
	b4s(j1)=b4s(j1)*1.d-10
	b5s(j1)=b5s(j1)*1.d-10
	b6s(j1)=b6s(j1)*1.d-10
136	continue
	ur3=ur3*1.d-10
	ur4=ur4*1.d-10
	ur5=ur5*1.d-10
	ur6=ur6*1.d-10
	vr3=vr3*1.d-10
	vr4=vr4*1.d-10
	vr5=vr5*1.d-10
	vr6=vr6*1.d-10
	dur3=dur3*1.d-10
	dur4=dur4*1.d-10
	dur5=dur5*1.d-10
	dur6=dur6*1.d-10
	dvr3=dvr3*1.d-10
	dvr4=dvr4*1.d-10
	dvr5=dvr5*1.d-10
	dvr6=dvr6*1.d-10
	x13=x13*1.d-10
	x14=x14*1.d-10
	x15=x15*1.d-10
	x16=x16*1.d-10
	x23=x23*1.d-10
	x24=x24*1.d-10
	x25=x25*1.d-10
	x26=x26*1.d-10
130	continue
c
c	  write(6,*)'finished b3 - b6'
c *	Now have available solutions b1 - b6 at source radius and
c *	b3 - b6 at surface.
cc	Next require inverse of A at source radius, where dy/dr=Ay.
	call matra2(am1,isdep)
	ylim(1)=dsqrt((dble(2*l+1))/(4.d0*dble(pi)))
	ylim(2)=-0.5d0*ylim(1)*dsqrt(fl21)
	ylim(3)=0.25d0*ylim(1)*dsqrt(dble((l-1)*(l+2))*fl21)
	  write(6,*)'ylim=',ylim
c	Do m=0, 1, 2 in that order.
c	June 29, 1999.  Then do vertical force, which is a m=0 deformation
c	represented by m1=4 below.
	m1=-1
38 	m1=m1+1
	if(m1.gt.4) go to 138
c	  write(6,*)'m1=',m1
c	Determine matrix elements of 6 X 6 boundary condition matrix bc0.
cc	First determine 4 X 4 a00 from eqn.'s (13)-(15) of notes.
	mus=dble(mu(isdep))
	lams=dble(kappa(isdep)-2.*mu(isdep)/3.)
	biga=lams+2.d0*mus
	r=dble(rb(isdep))
	  if(m1.eq.1) write(6,*) isdep,r,dble(rt(isdep-1))
	if(m1.gt.0) go to 163
	a00(1,1)=r*r
	a00(1,2)=0.d0
	a00(1,3)=0.d0
	a00(1,4)=0.d0
	a00(2,1)=0.d0
	a00(2,2)=r*r
	a00(2,3)=0.d0
	a00(2,4)=0.d0
	a00(3,1)=0.d0
	a00(3,2)=0.d0
	a00(3,3)=-2.d0*r*r
	a00(3,4)=0.d0
	a00(4,1)=0.d0
	a00(4,2)=0.d0
	a00(4,3)=0.d0
	a00(4,4)=r*r
	c0(1)=ylim(1)
	c0(2)=0.d0
	c0(3)=0.d0
	c0(4)=0.d0
	c0(5)=0.d0
	c0(6)=0.d0
	go to 134
163	if(m1.eq.2) go to 151
cc	For m1=1,3 have a00 always equal to am1.
c	do 142 i=1,4
c	do 141 j=1,4
c141	a00(i,j)=am1(i,j)
c142	continue
c
	if(m1.eq.3) go to 152
	if(m1.eq.4) go to 154
c	below lines are then for m1=1 (second type of m=0 deformation).
	a00(1,1)=1.d0
	a00(1,2)=0.d0
	a00(1,3)=0.d0
	a00(1,4)=0.d0
	a00(2,1)=0.d0
	a00(2,1)=0.d0
	a00(2,2)=1.d0
	a00(2,3)=0.d0
	a00(2,4)=0.d0
	a00(3,1)=0.d0
	a00(3,2)=0.d0
	a00(3,3)=1.d0
	a00(3,4)=0.d0
	a00(4,1)=0.d0
	a00(4,2)=0.d0
	a00(4,3)=0.d0
	a00(4,4)=1.d0
	c0(1)=-2.d0*ylim(1)*(1.d0/(r*r))*(mus/biga-0.5d0)
	c0(2)=(-2.d0*mus/(r*r*r))*ylim(1)*(3.d0-4.d0*mus/biga)
	c0(3)=0.d0
	c0(4)=-0.5d0*c0(2)
	c0(5)=0.d0
	c0(6)=0.d0
	go to 134
151	a00(2,1)=r*r
	a00(2,2)=0.d0
	a00(2,3)=0.d0
	a00(2,4)=0.d0
	a00(3,1)=0.d0
	a00(3,2)=r*r
	a00(3,3)=0.d0
	a00(3,4)=0.d0
	a00(1,1)=0.d0
	a00(1,2)=0.d0
	a00(1,3)=r*r*fl21/2.d0
	a00(1,4)=0.d0
	a00(4,1)=0.d0
	a00(4,2)=0.d0
	a00(4,3)=0.d0
	a00(4,4)=r*r
	c0(1)=ylim(2)
	c0(2)=0.d0
	c0(3)=0.d0
	c0(4)=0.d0
	c0(5)=0.d0
	c0(6)=0.d0
	go to 134
152	continue
	a00(1,1)=r*r
	a00(1,2)=0.d0
	a00(1,3)=0.d0
	a00(1,4)=0.d0
	a00(2,1)=0.d0
	a00(2,2)=r*r*r
	a00(2,3)=0.d0
	a00(2,4)=0.d0
	a00(3,1)=0.d0
	a00(3,2)=0.d0
	a00(3,3)=-2.d0*r*r
	a00(3,4)=0.d0
	a00(4,1)=0.d0
	a00(4,2)=0.d0
	a00(4,3)=0.d0
	a00(4,4)=1.d0
	c0(1)=0.d0
	c0(2)=0.d0
	c0(3)=0.d0
	c0(4)=-0.5d0*ylim(1)*
     &	dsqrt(dble((l+2)*(l-1))/dble(l*(l+1)))*mus/(r*r*r)
	c0(5)=0.d0
	c0(6)=0.d0
	go to 134
154	continue
	a00(1,1)=1.d0
	a00(1,2)=0.d0
	a00(1,3)=0.d0
	a00(1,4)=0.d0
	a00(2,1)=0.d0
	a00(2,2)=1.d0
	a00(2,3)=0.d0
	a00(2,4)=0.d0
	a00(3,1)=0.d0
	a00(3,2)=0.d0
	a00(3,3)=1.d0
	a00(3,4)=0.d0
	a00(4,1)=0.d0
	a00(4,2)=0.d0
	a00(4,3)=0.d0
	a00(4,4)=1.d0
	c0(1)=0.d0
	c0(2)=-ylim(1)/(r*r)
	c0(3)=0.d0
	c0(4)=0.d0
	c0(5)=0.d0
	c0(6)=0.d0
c	First four rows of bc0 correspond to boundary conditions
c	at source.
134	do 135 i=1,6
	do 137 j=1,4
	if(i.eq.1) c(j)=-b1s(j)
	if(i.eq.2) c(j)=-b2s(j)
	if(i.eq.3) c(j)=b3s(j)
	if(i.eq.4) c(j)=b4s(j)
	if(i.eq.5) c(j)=b5s(j)
	if(i.eq.6) c(j)=b6s(j)
137	continue
	do 140 j=1,4
	bc0(j,i)=0.d0
	do 145 k=1,4
	bc0(j,i)=bc0(j,i)+a00(j,k)*c(k)
145	continue
140	continue
135	continue
c	5th row of bc0 corresponds to zero shear stress at surface.
c	6th row of bc0 corresponds to zero normal stress at surface.
	bc0(5,1)=0.d0
	bc0(5,2)=0.d0
	bc0(5,3)=b3(2)
	bc0(5,4)=b4(2)
	bc0(5,5)=b5(2)
	bc0(5,6)=b6(2)
	bc0(6,1)=0.d0
	bc0(6,2)=0.d0
	bc0(6,3)=b3(4)
	bc0(6,4)=b4(4)
	bc0(6,5)=b5(4)
	bc0(6,6)=b6(4)
c	err type strain (m1=0), ett+epp (m1=1),
c	ert,erp (m1=2), or ett-epp,etp (m1=3).
c	  if(m1.eq.1) write(6,*)'c0=',c0
cTE
	if(l.eq.1) then
	bc0(1,1)=1.e+6
	bc0(2,1)=1.e+6
	bc0(3,1)=1.e+6
	bc0(4,1)=1.e+6
	endif
c--
	call ainver(a0,bc0,c0,6)
c	  if(m1.eq.1) write(6,*)'a0=',a0
c	Determine displacements and their derivatives at surface.
	if(odep.gt.dble(dep)) go to 24
	y1=real(a0(1)*dur1+a0(2)*dur2)
	y2=real(a0(1)*ur1+a0(2)*ur2)
	y3=real(a0(1)*dvr1+a0(2)*dvr2)
	y4=real(a0(1)*vr1+a0(2)*vr2)
	go to 25
24	y1=real(a0(3)*dur3+a0(4)*dur4+a0(5)*dur5+a0(6)*dur6)
	y2=real(a0(3)*ur3+a0(4)*ur4+a0(5)*ur5+a0(6)*ur6)
	y3=real(a0(3)*dvr3+a0(4)*dvr4+a0(5)*dvr5+a0(6)*dvr6)
	y4=real(a0(3)*vr3+a0(4)*vr4+a0(5)*vr5+a0(6)*vr6)
c		write(6,*)'odep',odep,'gt dep',dep
c		write(6,*)'a0(3),a0(34),a0(5),a0(6)=',a0(3),a0(4),a0(5),a0(6)
c		write(6,*) 'dur3,dur4,dur5,dur6=',dur3,dur4,dur5,dur6
c		write(6,*) 'LOOK','l=',l,'m1=',m1,y1,y2,y3,y4
25	write(2) l,y1,y2,y3,y4
c	TEST the surface BCs.
	c2cof=real(a0(3)*b3(2)+a0(4)*b4(2)+a0(5)*b5(2)+a0(6)*b6(2))
	c4cof=real(a0(3)*b3(4)+a0(4)*b4(4)+a0(5)*b5(4)+a0(6)*b6(4))
	write(6,*)'Normal stress at surface=',c2cof
	write(6,*)'Shear  stress at surface=',c4cof
c	write(6,*)' '
c	write(6,*) 'L,Y1,Y2,Y3,Y4=',l,y1,y2,y3,y4
c	write(6,*)' '
c	  shs=real(y1+2.d0*lams*y2/biga-lams*fl21*y4/biga)*real(biga)
c	  shs=real(a0(3)*b3(2)+a0(4)*b4(2)+a0(5)*b5(2)+a0(6)*b6(2))
c	  write(6,*)'shs=',shs
	  ub=(a0(1)*b1s(1)+a0(2)*b2s(1))
	  sb=(a0(1)*b1s(2)+a0(2)*b2s(2))
	  ua=(a0(3)*b3s(1)+a0(4)*b4s(1)+
     &	  a0(5)*b5s(1)+a0(6)*b6s(1))
	  sa=(a0(3)*b3s(2)+a0(4)*b4s(2)+
     &	  a0(5)*b5s(2)+a0(6)*b6s(2))
	  vb=(a0(1)*b1s(3)+a0(2)*b2s(3))
	  rx=(a0(1)*b1s(4)+a0(2)*b2s(4))
	  va=(a0(3)*b3s(3)+a0(4)*b4s(3)+
     &	  a0(5)*b5s(3)+a0(6)*b6s(3))
	  ry=(a0(3)*b3s(4)+a0(4)*b4s(4)+
     &	  a0(5)*b5s(4)+a0(6)*b6s(4))
	  write(6,*)'just above source:'
	  write(6,*)'u,nor,v,sh=',ua,sa,va,ry
	  write(6,*)'just below source:'
	  write(6,*)'u,nor,v,sh=',ub,sb,vb,rx
c**16AUG2006	Use following lines to eliminate dilatation,
c	for use in isolating uplift effect --> stat0Aupl
c	x11=0.
c	x12=0.
c	x13=0.
c	x14=0.
c	x15=0.
c	x16=0.
c**--
c**16AUG2006	Use following lines to eliminate uplift effect,
c	for use in isolating dilatation effect --> stat0Adil
c	x21=0.
c	x22=0.
c	x23=0.
c	x24=0.
c	x25=0.
c	x26=0.
c**--
	bfacl=real((x11+x21)*a0(1)+(x12+x22)*a0(2)+(x13+x23)*a0(3)
     &	+(x14+x24)*a0(4)+(x15+x25)*a0(5)+(x16+x26)*a0(6))
	bfacm=real((x21)*a0(1)+(x22)*a0(2)+(x23)*a0(3)
     &	+(x24)*a0(4)+(x25)*a0(5)+(x26)*a0(6))
	write(2) bfacl,bfacm
	go to 38
138	continue
	  write(6,*)'starting TOROIDAL modes'
c	TOROIDAL MODES
c	Start matrix propagation upward for solution t1.  Stop
c	at layer [isdep-1].
	is=0
	depmax=dep-2.5*(2.*pi/(real(l)+0.5))
	depmin=dep+2.5*(2.*pi/(real(l)+0.5))
	do 230 i=1,isdep-1
	if(i.lt.2) go to 230
c	Start integration at most 2 wavelengths into the earth.
	if(rt(i).lt.depmax.or.rb(i).gt.depmin) go to 230
	if(dble(depmax).ge.odep.or.dble(depmin).le.odep) go to 230
	if(rt(i).lt.depmax) go to 230
c	Give boundary conditions at top of layer (i-1) appropriate for
c	homogeneous solid in 0 < r < rt(i-1).
	if(is.ne.0) go to 228
	is=1
c	  if(l.gt.397) write(6,*)'l=',l,'start at bot of layer',i
c	  if(l.gt.397) write(6,*)'bot radius=',rb(i)
	j=i-1
	call matra(ajt,j,1)
	do 229 k=1,2
	t1(k)=ajt(k,1)
229	continue
c--	If following lines are not commented out, they
c	give boundary conditions appropriate for an incompressible
c	fluid in the sphere 0 < r < rt(i-1).
	t1(1)=1.d0
	t1(2)=0.d0
c--
228	continue
	dr=dble(rt(i)-rb(i))/dble(nrad)
	r=dble(rb(i))-dr/2.d0
	mus=dble(mu(i))
	call matra(ajt,i,0)
	call ainver(dt1,ajt,t1,2)
	do 231 j=1,nrad
	r=r+dr 
	r0=r
	call matra(ajt,i,2)
	call prodr(ajt,dt1,2,t1) 
	if(r.ge.odep) go to 231
	wr1=t1(1)
	dwr1=t1(1)/r+t1(2)/mus
231	continue
c	Make sure that matrix propagation is done accurately.
	call matra(ajt,i,1)
	call prodr(ajt,dt1,2,t1)
	  bmag=abs(real(t1(1))/r)+abs(real(t1(2)/mus))
	if(bmag.lt.1.e+10) go to 230
	t1(1)=t1(1)*1.d-10
	t1(2)=t1(2)*1.d-10
	wr1=wr1*1.d-10
	dwr1=dwr1*1.d-10
230	continue
c	Store t1 at source radius.
	do 250 j1=1,2
	t1s(j1)=t1(j1)
250	continue
c	  write(6,*)'finished solution t1'
c	Start matrix propagation upward for solutions t2,t3.  Stop
c	at layer [n].
c	Give boundary conditions at top of layer [isdep-1] appropriate for
c	homogeneous solid in this layer.
	j=isdep-1
	call matra(ajt,j,1)
	do 227 k=1,2
	t2(k)=ajt(k,1)
	t3(k)=ajt(k,2)
c	Store these solutions at source radius.
	t2s(k)=t2(k)
	t3s(k)=t3(k)
227	continue
	do 330 i=isdep,n
	if(rb(i).gt.depmin) go to 330
	mus=dble(mu(i))
c	Determine coefficient matrices dt2 - dt3 for layer i.
	do 333 m=1,2
	if(m.eq.1) call equal(t2,ct,2)
	if(m.eq.2) call equal(t3,ct,2)
	call matra(ajt,i,0)
	if(m.eq.1) call ainver(dt2,ajt,ct,2)
	if(m.eq.2) call ainver(dt3,ajt,ct,2)
333	continue
	dr=dble(rt(i)-rb(i))/dble(nrad)
	r=dble(rb(i))-dr/2.d0
	do 331 j=1,nrad
	r=r+dr 
	r0=r
	call matra(ajt,i,2)
	do 232 m=1,2
	if(m.eq.1) call prodr(ajt,dt2,2,ct) 
	if(m.eq.2) call prodr(ajt,dt3,2,ct) 
	if(m.eq.1) call equal(ct,t2,2)
	if(m.eq.2) call equal(ct,t3,2)
232	continue
	if(r.ge.odep) go to 331
	wr2=t2(1)
	dwr2=t2(1)/r+t2(2)/mus
	wr3=t3(1)
	dwr3=t3(1)/r+t3(2)/mus
331	continue
c	Make sure that matrix propagation is done accurately.
	do 253 m=1,2
	if(m.eq.1) call equal(t2,ct,2)
	if(m.eq.2) call equal(t3,ct,2)
	call matra(ajt,i,1)
	if(m.eq.1) call prodr(ajt,dt2,2,ct)
	if(m.eq.2) call prodr(ajt,dt3,2,ct)
	if(m.eq.1) call equal(ct,t2,2)
	if(m.eq.2) call equal(ct,t3,2)
253	continue
	  bmag=abs(real(t2(1))/r)+abs(real(t2(2)/mus))
	if(bmag.lt.1.e+10) go to 330
	t2(1)=t2(1)*1.d-10
	t2(2)=t2(2)*1.d-10
	t3(1)=t3(1)*1.d-10
	t3(2)=t3(2)*1.d-10
	t2s(1)=t2s(1)*1.d-10
	t2s(2)=t2s(2)*1.d-10
	t3s(1)=t3s(1)*1.d-10
	t3s(2)=t3s(2)*1.d-10
	wr2=wr2*1.d-10
	dwr2=dwr2*1.d-10
	wr3=wr3*1.d-10
	dwr3=dwr3*1.d-10
330	continue
c
c	Now have available t1 - t3 at source radius, and t2 - t3 at surface.
cc	Next require inverse of A, where dy/dr=Ay.
	mus=dble(mu(isdep))
	r=dble(rb(isdep))
	det=-dble(l*l+l+1)/(r*r)
	am1(1,1)=(-3.d0/r)/det
	am1(1,2)=-(1.d0/mus)/det
	am1(2,1)=-mus*dble((l-1)*(l+2))/(r*r)
	am1(2,2)=(1.d0/r)/det
c	Do m=1, 2 in that order.
	m1=0
238	m1=m1+1
	if(m1.gt.2) go to 50
c	  write(6,*)'m1=',m1
c	Determine matrix elements of 4 X 4 boundary condition matrix bc1.
c	First determine 2 X 2 a00 from eqn's (25)-(26) of notes.
	if(m1.gt.1) go to 252
	a00(1,1)=r*r*fl21/2.d0
	a00(1,2)=0.d0
	a00(2,1)=0.d0
	a00(2,2)=r*r
	c1(1)=ylim(2)
	c1(2)=0.d0
	c1(3)=0.d0
	go to 234
252	a00(1,1)=r*r
	a00(1,2)=0.d0
	a00(2,1)=0.d0
	a00(2,2)=1.d0
	c1(1)=0.d0
	c1(2)=-0.5d0*ylim(1)*
     &	dsqrt(dble((l+2)*(l-1))/dble(l*(l+1)))*mus/(r*r*r)
	c1(3)=0.d0
c	First two rows of bc1 correspond to boundary conditions at source.
234	do 235 i=1,3
	do 237 j=1,2
	if(i.eq.1) c(j)=-t1s(j)
	if(i.eq.2) c(j)=t2s(j)
	if(i.eq.3) c(j)=t3s(j)
237	continue
	do 240 j=1,2
	bc1(j,i)=0.d0
	do 245 k=1,2
	bc1(j,i)=bc1(j,i)+a00(j,k)*c(k)
245	continue
240	continue
235	continue
c	3rd row of bc1 corresponds to zero shear stress at surface.
	bc1(3,1)=0.d0
	bc1(3,2)=t2(2)
	bc1(3,3)=t3(2)
c	ert, erp type strain (m1=1) or ett-epp, etp (m1=2).
cTE
	if(l.eq.1) then
	bc1(1,1)=1.e+6
	bc1(2,1)=1.e+6
	endif
c--
	call ainver(a1,bc1,c1,3)
c	Determine displacement and their derivatives at surface.
	  write(6,*)'odep,dep=',odep,dep
	if(odep.gt.dble(dep)) go to 124
	y1=real(a1(1)*(dwr1-wr1))
	y2=real(a1(1)*wr1)
	go to 125
124	y1=real(a1(2)*(dwr2-wr2)+a1(3)*(dwr3-wr3))
cc	ABOVE LINE is zero when evaluated at surface.  So to avoid underflow,
cc	just set it equal to zero.
c	y1=0.
	y2=real(a1(2)*wr2+a1(3)*wr3)
125	write(2) l,y1,y2
	  wb=real(a1(1)*t1s(1))
	  sb=real(a1(1)*t1s(2))
	  wa=real(a1(2)*t2s(1)+a1(3)*t3s(1))
	  sa=real(a1(2)*t2s(2)+a1(3)*t3s(2))
c	  write(6,*)'t1s=',t1s
c	  write(6,*)'t2s=',t2s
c	  write(6,*)'t3s=',t3s
c	  write(6,*)'t2=',t2
c	  write(6,*)'t3=',t3
c	  write(6,*)'a1=',a1
	  write(6,*)'w,sh just above=',wa,sa
	  write(6,*)'w,sh just below=',wb,sb
c	  ssur=real(a1(2)*t2(2)+a1(3)*t3(2))
c	  write(6,*)'surface traction=',ssur
	go to 238
50	continue
	go to 44
94	continue
	close(2)
	end 
	subroutine matra1(aj,i,n)
c	Compute matrix elements for spheroidal modes.
c	n=0: evaluate at bottom of layer i
c	n=1: evaluate at top of layer i
c	Oherwise evaluate at r0.
	real*8 r0
	real*4 mu,kappa
	real*8 mus,lams,biga,r,flp,flm,r2,flp3,flm3
	real*8 a1p,a1m,a2p,a2m,a3p,a3m,a4p,a4m
	real*8 b1p,b1m,b2p,b2m,b3p,b3m,b4p,b4m,rlp,rlm  
	common/mat/rb(200),rt(200),kappa(200),mu(200),l,r0
	real*8 aj(4,4)
	r=r0
	mus=dble(mu(i))
	lams=dble(kappa(i)-2.*mu(i)/3.)
	biga=lams+2.d0*mus   
10	if(n.eq.1) r=dble(rt(i))
	if(n.eq.0) r=dble(rb(i))
c	Determine matrix elements (Takeuchi and Saito, 1972, p. 243-244).
	flp=dble(l)
	flm=-flp-1.d0
	flp3=2.d0*(2.d0*flp+3.d0)
	flm3=2.d0*(2.d0*flm+3.d0)
	r2=r*r
	aj(1,1)=flp
	aj(1,3)=flm
	aj(2,1)=2.d0*mus*(flp*(flp-1.d0))
	aj(2,3)=2.d0*mus*(flm*(flm-1.d0))
	aj(3,1)=(1.d0)
	aj(3,3)=(1.d0)
	aj(4,1)=2.d0*mus*(flp-1.d0)
	aj(4,3)=2.d0*mus*(flm-1.d0)
	a1p=-(flp+2.d0)*r2/flp3
	a1m=-(flm+2.d0)*r2/flm3
	a2p=-(lams+2.d0*mus*(flp+2.d0)*(flp+1.d0)/flp3)*r2
	a2m=-(lams+2.d0*mus*(flm+2.d0)*(flm+1.d0)/flm3)*r2
	a3p=-r2/flp3
	a3m=-r2/flm3
	a4p=-2.d0*mus*(flp+1.d0)*r2/flp3
	a4m=-2.d0*mus*(flm+1.d0)*r2/flm3
	b1p=flp*(flp+1.d0)*r2/flp3
	b1m=flm*(flm+1.d0)*r2/flm3
	b2p=2.d0*mus*flp*(flp+1.d0)*(flp+1.d0)*r2/flp3
	b2m=2.d0*mus*flm*(flm+1.d0)*(flm+1.d0)*r2/flm3
	b3p=(flp+3.d0)*r2/flp3
	b3m=(flm+3.d0)*r2/flm3
	b4p=mus*2.d0*flp*(flp+2.d0)*r2/flp3
	b4m=mus*2.d0*flm*(flm+2.d0)*r2/flm3
	aj(1,2)=mus*(flp+1.d0)*a1p+biga*b1p
	aj(1,4)=mus*(flm+1.d0)*a1m+biga*b1m
	aj(2,2)=mus*(flp+1.d0)*a2p+biga*b2p
	aj(2,4)=mus*(flm+1.d0)*a2m+biga*b2m
	aj(3,2)=mus*(flp+1.d0)*a3p+biga*b3p
	aj(3,4)=mus*(flm+1.d0)*a3m+biga*b3m
	aj(4,2)=mus*(flp+1.d0)*a4p+biga*b4p
	aj(4,4)=mus*(flm+1.d0)*a4m+biga*b4m 
c	Normalize matrix elements.
	rlp=(r/dble(rb(i)))**l
	rlm=(dble(rb(i))/r)**(l+1)
	do 20 m=1,2
	aj(1,m)=aj(1,m)*rlp/r
	aj(1,m+2)=aj(1,m+2)*rlm/r 
	aj(2,m)=aj(2,m)*rlp/r2
	aj(2,m+2)=aj(2,m+2)*rlm/r2
	aj(3,m)=aj(3,m)*rlp/r
	aj(3,m+2)=aj(3,m+2)*rlm/r 
	aj(4,m)=aj(4,m)*rlp/r2
	aj(4,m+2)=aj(4,m+2)*rlm/r2
20	continue
	return
	end 
	subroutine matra(aj,i,n)
c	Compute matrix elements for toroidal modes.
c	n=0: evaluate at bottom of layer i
c	n=1: evaluate at top of layer i
c	Oherwise evaluate at r0.
	real*8 r0
	real*4 mu,kappa
	real*8 mus
	common/mat/rb(200),rt(200),kappa(200),mu(200),l,r0
	real*8 aj(2,2)
	r=real(r0)
	mus=dble(mu(i))
	if(n.eq.1) r=rt(i)
	if(n.eq.0) r=rb(i)
	aj(1,1)=dble((r/rb(i))**l)
	aj(1,2)=dble((rb(i)/r)**(l+1))
	aj(2,1)=mus*dble(l-1)*dble((r/rb(i))**(l-1)*(1./rb(i)))
	aj(2,2)=-mus*dble(l+2)*dble((rb(i)/r)**(l+2)*(1./rb(i)))
	return
	end
	subroutine matra2(am1,i)
c	Compute matrix A which belongs with dy/dr = Ay (spheroidal modes)
c	in bottom of layer i.
c	Then return its inverse in matrix am1.
	real*8 r0
	real*4 mu,kappa
	real*8 mus,lams,gams,biga,r,fl21
	common/mat/rb(200),rt(200),kappa(200),mu(200),l,r0
	real*8 a(4,4),am1(4,4)
	r=dble(rb(i))
	mus=dble(mu(i))
	lams=dble(kappa(i)-2.*mu(i)/3.)
	biga=lams+2.d0*mus 
	gams=lams+mus-lams*lams/biga
	fl21=dble(l*(l+1))
c	Specify matrix A.
	a(1,1)=-2.d0*lams/(biga*r)
	a(1,2)=1.d0/biga
	a(1,3)=lams*fl21/(biga*r)
	a(1,4)=0.d0
	a(2,1)=4.d0*gams/(r*r)
	a(2,2)=2.d0*(lams/biga-1.d0)/r
	a(2,3)=-2.d0*gams*fl21/(r*r)
	a(2,4)=fl21/r
	a(3,1)=-1.d0/r
	a(3,2)=0.d0
	a(3,3)=1.d0/r
	a(3,4)=1.d0/mus
	a(4,1)=-2.d0*gams/(r*r)
	a(4,2)=-lams/(biga*r)
	a(4,3)=(-2.d0*mus+(gams+mus)*fl21)/(r*r)
	a(4,4)=-3.d0/r
	call matinv(a,am1,4)
	return
	end
	subroutine matinv(b,bm1,n)
c       Find inverse bm1 of n X n matrix b.
	real*8 bs(4,4),bm1(4,4),b(4,4),a(4),c(4)
c	Determine inverse one row at a time.
	do 5 i=1,n
	do 10 j=1,n
	c(j)=0.d0
	if(i.eq.j) c(j)=1.d0
10	continue
c	Write b into bs before calling ainver.
	do 15 j=1,n
	do 20 k=1,n
	bs(j,k)=b(j,k)
20	continue
15	continue
c
	call ainver(a,bs,c,n)
	do 25 j=1,n
	bm1(j,i)=a(j)
25	continue
5	continue
	return
	end
	subroutine prodr(a,c,n,b)
c	Form matrix product b=(A)c
	real*8 b(n) 
	real*8 a(n,n),c(n)
	do 5 i=1,n
	b(i)=0.d0
	do 4 j=1,n
4	b(i)=b(i)+a(i,j)*c(j)
5	continue
	return
	end 
	subroutine ainver(a,b,c,n)
c       Find solution (a) of the matrix equation (B)(a)=(c).
c       Method:  Gaussian elimination.
	integer*2 perm 
	dimension perm(10)
	real*8 a(n),b(n,n),c(n),bsave,fac
c		if(n.eq.6) write(6,*) 'AINVER: b=',b
        do 5 i=1,n
5       perm(i)=i
	i=0
10      i=i+1
	if (i.gt.n) go to 35
c       Find maximum in row i.
	amax=0.
	imax=i
	do 15 j=i,n
	t=real(b(i,j))
	if (t.lt.0.) t=-t
	if (t.lt.amax) go to 15
	amax=t
	imax=j
15      continue
	j=imax
c       Switch columns i and j.
	do 20 m=1,n
	bsave=b(m,i)
	b(m,i)=b(m,j)
20      b(m,j)=bsave
	iperm=perm(i)
	perm(i)=perm(j)
	perm(j)=iperm
c       Eliminate ith column.
c	  if(n.eq.6) write(6,*)'pivot row',i,'=',b(i,i)
c	  if(i.eq.n) pause  
	do 25 j=1,n
	if (j.eq.i) go to 25
	fac=b(j,i)/b(i,i)
	do 30 k=i,n
30      b(j,k)=b(j,k)-fac*b(i,k)
	c(j)=c(j)-fac*c(i)
25      continue
	go to 10
35      do 40 i=1,n
	k=perm(i)
40      a(k)=c(i)/b(i,i)
	return
	end
	subroutine equal(a,b,n)
c	Set the n x 1 matrix (b) equal to (a).
	real*8 a(n),b(n)
	do 10 j=1,n
10	b(j)=a(j)
	return
	end 

