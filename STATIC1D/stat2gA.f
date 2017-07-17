c	Program STAT2gA
c***	JULY 1, 2005
c	Version of stat2A which allows number of points in downdip integration
c	to be specified with [ndep].
	parameter (ndep=33)
c***
c	Version of STAT1 which uses Green's functions for
c	ndep different depths and reads them all from STAT0.OUT.
c	Thus ndep times more space is required as for STAT1B.
c	Use LEGL subroutines which
c	calulate spherical harmonics in steps of [lste].  All quantities
c	must be then multiplied by [lste] before output.
c***	THIS VERSION: 9 Nov 2001.  Response Greens functions
c	at appropriate depths are evaluated with spherical harmonic sum
c	ahead of time, then stored, leading to extremely rapid computation!
c***
c	Determine displacement 
c	and horizontal strain fields for spheroidal motion (coseismic).
c	Both the vertical and longitudinal components are written out.
c	Read in displacement coefficients from
c	input file 'stat0.out'.  
c	Also determine gravity anomaly in units of mGal.
c	This program handles finite faults.  The dip, strike, and rake follow
c	the convention of Ben Menahem and Singh...
c	If the rake is >181 (<-181) deg., a tensile (compressional) dislocation
c	is assumed.
c**
c	NOVEMBER 19, 1996: Furthermore, if |rake| is > 210, then compensate
c	the tensile (compressional) mechanism with an isotropic component
c	such that tr(M)=0.
c**
c	The input location of each fault segment 
c	is at the lowermost corner of the fault segment
c	rectangle closest to its strike direction.
c	*** NOTE ***
c	Handles dipping faults with any shear dislocation.
c	This program computes displacements and strains at any number
c	of observation points (up to 3600), and handles faults composed
c	of several segments of varying length and strike (up to
c	2000 segments).  The dip and depth range of each fault segment
c	are fixed.  Surface observations assumed.
c**
c	Output strains in units of 10**(-6) .
	real*4 mrr,mtt,mpp,mrt,mrp,mtp  
	real*4 kappa,mu,lam
	common/source/mrr,mtt,mpp,mrt,mrp,mtp,depth,slat,slon
	real*8 cl,sl,xl0(25001),dxl0(25001),xl1(25001),dxl1(25001)
	real*8 xl2(25001),dxl2(25001)
	dimension y1(7*ndep,25000)
	dimension y2(7*ndep,25000)
	dimension y3(5*ndep,25000)
	dimension y4(5*ndep,25000)
	dimension fm(6)
	dimension x1s(25000),x2s(25000),x3s(25000)
	dimension dx1s(25000),dx2s(25000),dx3s(25000)
	dimension ddx1s(25000),ddx2s(25000),ddx3s(25000)
	dimension olat(3600),olon(3600),flat(2000),flon(2000)
     	dimension fbigl(2000),fstr(2000)
	dimension fwt(2000),frake(2000)
	dimension bfacl(5*ndep,25000),bfacm(5*ndep,25000)
	dimension deltr(40)
	dimension rd1g(6,ndep,40)
	dimension rd1S(6,ndep,40),r1hS(6,ndep,40),r1tS(6,ndep,40)
	dimension rttS(6,ndep,40),rppS(6,ndep,40),rtrS(6,ndep,40)
	dimension rtpS(6,ndep,40),rprS(6,ndep,40),rrrS(6,ndep,40)
	dimension romT(6,ndep,40)
	dimension r1hT(6,ndep,40),r1tT(6,ndep,40)
	dimension rttT(6,ndep,40),rppT(6,ndep,40),rtrT(6,ndep,40)
	dimension rtpT(6,ndep,40),rprT(6,ndep,40)
	dimension vd1g(6)
	dimension vd1S(6),v1hS(6),v1tS(6),vttS(6),vppS(6)
	dimension vtrS(6),vtpS(6),vprS(6),vrrS(6)
	dimension vomT(6)
	dimension v1hT(6),v1tT(6),vttT(6),vppT(6)
	dimension vtrT(6),vtpT(6),vprT(6)
	real*8 yd1g(40)
	real*8 yd1S(40),y1hS(40),y1tS(40),yttS(40),yppS(40)
	real*8 ytrS(40),ytpS(40),yprS(40),yrrS(40)
	real*8 yomT(40)
	real*8 y1hT(40),y1tT(40),yttT(40),yppT(40)
	real*8 ytrT(40),ytpT(40),yprT(40)
	real*8 sarr(40)
	real*8 sd1g(6,ndep,40)
	real*8 sd1S(6,ndep,40),s1hS(6,ndep,40),s1tS(6,ndep,40),sttS(6,ndep,40),
     &	sppS(6,ndep,40)
	real*8 strS(6,ndep,40),stpS(6,ndep,40),sprS(6,ndep,40),srrS(6,ndep,40)
	real*8 somT(6,ndep,40)
	real*8 s1hT(6,ndep,40),s1tT(6,ndep,40),sttT(6,ndep,40),sppT(6,ndep,40)
	real*8 strT(6,ndep,40),stpT(6,ndep,40),sprT(6,ndep,40)
c	integer*4 sid(10)
	real*8 yi1(ndep),yi2(ndep),yi3(ndep),yi4(ndep),bi1(ndep),bi2(ndep)
	real*8 sl1(ndep),sl2(ndep),sl3(ndep),sl4(ndep),sb1(ndep),sb2(ndep)
	real*8 dep,odep,rad1,rad10
	real*8 ypoint,evaleq 
	real*8 wkv1(7),wkv2(7),wkv3(5),wkv4(5),wb1(5),wb2(5)
	real*8 ymu(ndep),ylam(ndep),smu(ndep),slam(ndep)
	real*8 rmin,rmax,deltp
	real*4 len,kmin
	parameter (dstmax=21000.)
c	dstmax has the maximum fault-observation pt. distance. 
	common/maxl/lmax
c 
c	ethr=earth's radius in km.
c	bigr=earth radius in 10**6 cm.
	pi=3.1415926  
	twopi=2.*3.1415926
	rad=180./3.1415926
	open(2,file='stat0.out',form='unformatted')
	rewind(2)
	read(2) ethr,dmax,dmin
c	The displacement coefficients are evaluated at depths
c	dmin through dmax km in increments of (dmax-dmin)/(ndep-1) km.
	  write(6,*)'ethr=',ethr
	bigr=ethr/10.
	read(2) lmin,lmax
	  write(6,*)'lmin,lmax=',lmin,lmax
	kmin=twopi/(real(lmax)+0.5)
	read(2) odep
	do 55 j=1,ndep
	  write(6,*)'j=',j
	read(2) kappa,mu
	  write(6,*)'kappa,mu=',kappa,mu
	lam=kappa-2.*mu/3.
	ymu(j)=dble(mu)
	ylam(j)=dble(lam)
	  write(6,*)'ymu,ylam=',ymu(j),ylam(j)
c *
		lste=1
c *
	flste=real(lste)
	do 50 l=lmin,lmax,lste
	lr=l/lste
c		write(6,*)'l,lr=',l,lr
	do 69 idec=1,5
	k=ndep*idec-ndep+j
	read(2) ldum,y1(k,lr),y2(k,lr),y3(k,lr),y4(k,lr)
	read(2) bfacl(k,lr),bfacm(k,lr)
c		if(j.eq.ndep) write(6,*) ldum,y1(k,lr),y2(k,lr),y3(k,lr),y4(k,lr)
c		if(j.eq.ndep) write(6,*) bfacl(k,lr),bfacm(k,lr)
c		if(lr.eq.3900) write(6,*)'y1(',k,lr,')=',y1(k,lr)
69	continue
	do 58 idec=6,7
	k=ndep*idec-ndep+j
	read(2) ldum,y1(k,lr),y2(k,lr)
c		if(j.eq.ndep) write(6,*) ldum,y1(k,lr),y2(k,lr)
58	continue
c	  pause
c	Note: y1 - y4(1 thru 5*ndep,l) are Spheroidal motion coeff.
c	      y1 - y2(5*ndep+1 thru 7*ndep,l) are Toroidal motion coeff.
50	continue
55	continue
	close(2)
	call splneq(ndep,ymu,smu)
	call splneq(ndep,ylam,slam)
	rad10=dble(1.-dmax/ethr)
	rad1=dble(1.-dmin/ethr)
	  write(6,*)'rad1,rad10=',rad1,rad10
	write(6,*)'max depth, min depth (km), dip(deg.)?'
	read(5,*) dmax,dmin,dip
	cdip=cos(dip/rad)
	sdip=sin(dip/rad)
	c2dip=cdip*cdip-sdip*sdip
	dip=cos(dip/rad)/sin(dip/rad)
	rmin=1.d0-dble(dmax/ethr)
	rmax=1.d0-dble(dmin/ethr)
	  write(6,*)'rmin,rmax=',ethr,rmin,rmax
c 
	write(6,*)'finite fault with iseg # segments.  iseg=?'
	read(5,*) iseg
c	iobs=0
	fleng=0.
	do 39 i=1,iseg
	write(6,*)'segment #',i,'lat,lon(deg.),length(km)'
	write(6,*) 'strike(deg.),rake(deg.),slip(cm)?'
	read(5,*) flat(i),flon(i),fbigl(i),fstr(i),frake(i),fwt(i)
	fleng=fleng+fbigl(i)
	flat(i)=(pi/2.-flat(i)/rad)
	flon(i)=flon(i)/rad
	fstr(i)=fstr(i)/rad 
	frake(i)=frake(i)/rad 
39	continue
c	write(6,*)'total fault length =',fleng,'km'
c 
	write(6,*)'number of observation points?'
	read(5,*) ipts
	write(6,*)'ipts=',ipts
	do 40 i=1,ipts 
c	write(6,*)'point #',i,'lat,lon(deg.)?'
	read(5,*) olat(i),olon(i)
	olat(i)=(pi/2.-olat(i)/rad)
	olon(i)=olon(i)/rad
40	continue
c-----------
c	open(2,file='latlon.inDEF')
c	rewind(2)
c	read(2,*) ipts
c	write(6,*)'ipts=',ipts
c	do i=1,ipts 
c	write(6,*)'point #',i,'lat,lon(deg.)?'
c	read(2,*) olat(i),olon(i)
c	olat(i)=(pi/2.-olat(i)/rad)
c	olon(i)=olon(i)/rad
c	enddo
c	close(2)
c-----------
	write(6,*)'finished reading displacement coefficients'
c	Divide moment tensor by nmesh (# elements in fault subdivision).
c	NOTE: nmesh1 must be <= 17 (important for storage space
c	of interpolated spheroidal and toroidal coeff. later).
c	nmesh=21
c	nmesh=8
c	nmesh=2
c	nmesh1=11
c	Always use nmesh1=ndep in this program !!
	nmesh1=ndep
c	nmesh1=9
c	nmesh1=ndep
c	nmesh2=36
c	nmesh2=17
c	nmesh2=3
c	nmesh2=25 (used for 1838, 1857 etc modeling)
c	nmesh2=2
c	nmesh2=50
	nmesh2=25
c	  nmesh2=3*int(fleng)+2
c	nmesh2=2
c	nmesh2=21
c	nmesh1-1=# depth points
c	nmesh2-1=# horizontal length points 
	nmesh0=nmesh1
c	if(nmesh1.eq.2) nmesh0=17
	dh=(rmax-rmin)*dip/real(nmesh0-1)
c*****	
	write(6,*)'calculating spline interpolation coefficient'
	write(6,*)'and evaluation functions ahead of time'
c	Calculate the required spheroidal and toroidal motion
c	coefficients at interpolated depths.  It saves a lot of time to
c	do this ahead of time and retrieve these values later.
	do 150 l=lmin,lmax,lste
	lr=l/lste
	do 275 idec=1,7
c	Form arrays yi1--yi4  and bi1--bi2 containing eigenfunction
c	values at discrete radius points for use in spline interpolation.
	do 274 j=1,ndep
	k=ndep*idec-ndep+j
	yi1(j)=dble(y1(k,lr))
	yi2(j)=dble(y2(k,lr))
	if(idec.gt.5) go to 274
	yi3(j)=dble(y3(k,lr))
	yi4(j)=dble(y4(k,lr))
	bi1(j)=dble(bfacl(k,lr))
	bi2(j)=dble(bfacm(k,lr))
274	continue 
c		if(l.eq.3900) write(6,*)'l=',l,'IDEC=',idec
c		if(l.eq.3900) write(6,*)'yi1=',yi1
c		if(l.eq.3900) write(6,*)'yi2=',yi2
c		if(l.eq.3900) write(6,*)'yi1(1)=',yi1(1)
c		if(l.eq.3900) write(6,*)'yi1(2)=',yi1(2)
c		if(l.eq.3900) write(6,*)'yi1(3)=',yi1(3)
	call splneq(ndep,yi1,sl1)
	call splneq(ndep,yi2,sl2)
	if(idec.gt.5) go to 180
	call splneq(ndep,yi3,sl3)
	call splneq(ndep,yi4,sl4)
	call splneq(ndep,bi1,sb1)
	call splneq(ndep,bi2,sb2)
180	bigh=-dh/2.
c		write(6,*)'doing spline coeff., nmesh0=',nmesh0
	idip=0
181	idip=idip+1
	if(idip.eq.nmesh0) go to 275
	kdip=(nmesh0-1)*idec-(nmesh0-1)+idip
	bigh=bigh+dh
	dep=dble(rmin+bigh/dip)
c	if(dep.lt.rad10.or.dep.gt.rad1) write(6,*)'A.idip=',idip,'going to 272'
	if(dep.lt.rad10.or.dep.gt.rad1) go to 272
c	if(abs(dep-odep).lt.0.25*kmin)  write(6,*)'B.idip=',idip,'going to 272'
c		write(6,*)'dep,odep=',dep,odep
c		write(6,*)'kmin=',kmin
c		pause
cNOTE	if(abs(dep-odep).lt.0.25*kmin) go to 272
	ypoint=1.d0+(dep-rad1)*dble(ndep-1)/(rad10-rad1)
c	write over existing y1,y2,y3,y4,bfacl,bfacm arrays.
	y1(kdip,lr)=real(evaleq(ypoint,ndep,yi1,sl1))
c		a=real(evaleq(ypoint,ndep,yi1,sl1))
c		write(6,*)'kdip,lr=',kdip,lr
c		write(6,*)'a,y1=',a,y1(kdip,lr)
	y2(kdip,lr)=real(evaleq(ypoint,ndep,yi2,sl2))
c		if(l.eq.3900) write(6,*)'181 loop: y2(',kdip,',',lr,')=',y2(kdip,lr)
	if(idec.gt.5) go to 181
	y3(kdip,lr)=real(evaleq(ypoint,ndep,yi3,sl3))
	y4(kdip,lr)=real(evaleq(ypoint,ndep,yi4,sl4))	
	bfacl(kdip,lr)=real(evaleq(ypoint,ndep,bi1,sb1))
	bfacm(kdip,lr)=real(evaleq(ypoint,ndep,bi2,sb2))
	go to 181
c	NOTE--with 272--
c	If part of the fault depth range lies outside of the depth range
c	of the displacement eigenfunctions calculated by stat0A,
c	then do not include that fault depth.
c	Also, if the observation depth is within 0.25*kmin of the running
c	fault depth (kmin=minimum wavelength of expansion),
c	then do not include that fault depth.
272	continue
	y1(kdip,lr)=0.
	y2(kdip,lr)=0.
	if(idec.gt.5) go to 181
	y3(kdip,lr)=0.
	y4(kdip,lr)=0.
	bfacl(kdip,lr)=0.
	bfacm(kdip,lr)=0.
	go to 181
c	
275	continue
150	continue
c		pause
	write(6,*)'done calculating spline interpolation coefficients'
c *****
c	Next, calculate response Greens functions for m=0, 1, and 2 at 
c	specific values of DELTA.
	write(6,*)'calculate response Greens functions for m=0, 1, and 2 at' 
	write(6,*)'specific values of DELTA'
	do i=1,40
	deltr(i)=dstmax**(0.025*real(i))/ethr
	enddo
	do 105 i=1,40
	deltpf=deltr(i)
	write(6,*)'i=',i,' DELTA=',deltpf
	cotd=cos(deltpf)/sin(deltpf)
	slf=sin(deltpf)
	call lgndr0(deltpf,xl0,dxl0)
	call lgndr1(deltpf,xl1,dxl1)
	call lgndr2(deltpf,xl2,dxl2) 
c	Store Legendre functions for later use.
	do 364 lk=lmin,lmax+1
	faclk=sqrt(real(lk*lste*(lk*lste+1)))
c	Note multiplication by [lste] because of the stepped l-summation.
	x1=real(xl0(lk))
	x2=real(xl1(lk))
	x3=real(xl2(lk))
	dx1=real(dxl0(lk))
	dx2=real(dxl1(lk))
	dx3=real(dxl2(lk))
	ddx1=(-faclk*faclk)*x1-cotd*dx1
	ddx2=(1./(slf*slf)-faclk*faclk)*x2-cotd*dx2 
	ddx3=(4./(slf*slf)-faclk*faclk)*x3-cotd*dx3
	x1s(lk)=x1
	x2s(lk)=x2
	x3s(lk)=x3
	dx1s(lk)=dx1
	dx2s(lk)=dx2
	dx3s(lk)=dx3
	ddx1s(lk)=ddx1
	ddx2s(lk)=ddx2
	ddx3s(lk)=ddx3
364	continue 
	do 107 im=1,6
	fm(1)=0.
	fm(2)=0.
	fm(3)=0.
	fm(4)=0.
	fm(5)=0.
	fm(6)=0.
	fm(im)=1.	
	do 106 idip=1,nmesh0-1
	rd1g(im,idip,i)=0.
	rd1S(im,idip,i)=0.
	r1hS(im,idip,i)=0.
	r1tS(im,idip,i)=0.
	rttS(im,idip,i)=0.
	rppS(im,idip,i)=0.
	rtrS(im,idip,i)=0.
	rtpS(im,idip,i)=0.
	rprS(im,idip,i)=0.
	rrrS(im,idip,i)=0.
	romT(im,idip,i)=0.
	r1hT(im,idip,i)=0.
	r1tT(im,idip,i)=0.
	rttT(im,idip,i)=0.
	rppT(im,idip,i)=0.
	rtrT(im,idip,i)=0.
	rtpT(im,idip,i)=0.
	rprT(im,idip,i)=0.
c		if(im.ne.6) go to 106
	cp=1.
	sp=1.
	cp2=1.
	sp2=1.
	do 371 l=lmin,lmax,lste
c		if(l.gt.250) go to 371
	lr=l/lste
c	With facf apply tapering over last half of wavenumber range.
	fl=real(l)
	flmax=real(lmax)
	iw0m=lmax/2
	iw0x=lmax-iw0m+1
	fiw0m=real(iw0m)
	facf=1.0
	if(l.ge.iw0x) facf=1.0-cos((flmax-fl)/(fiw0m-1.) * real(pi)/2.)
c
	facl=sqrt(real(l*(l+1)))
c
	dfac=facf*(1.e+6)/ethr**2
	dfac1=dfac 
	dfacg=dfac*(real(l+1)/real(2*l+1))*8.3818e-4
c	Now account for the fact that the observation radius is generally
c	different from the surface radius.  This affects displ calculations
c	(for example, y2 and y4-arrays have functions
c	y1(r)/r and y3(r)/r, where r is dimensionless radius,
c	which equals 1 only at the surface), but not strain calculations.
	dfac=dfac*real(odep)
c	SPHEROIDAL MODES
	do 375 idec=1,5
	kdip=(nmesh0-1)*idec-(nmesh0-1)+idip
	wkv1(idec)=dble(y1(kdip,lr))
	wkv2(idec)=dble(y2(kdip,lr))
c		write(6,*)'kdip,lr,wkv2=',kdip,lr,wkv2(idec)
	wkv3(idec)=dble(y3(kdip,lr))
	wkv4(idec)=dble(y4(kdip,lr))
	wb1(idec)=dble(bfacl(kdip,lr))
	wb2(idec)=dble(bfacm(kdip,lr))
375	continue
c		write(6,*)'wkv4=',wkv4
c	Note extra factor of 2 in z2 - z5 components in order to account
c	for m=-1 and m=-2 contributions.
c * *	Moment tensor excitation
	z1=fm(1)*real(wkv2(1))+fm(2)*real(wkv2(2))
	z2=2.*fm(3)*real(wkv2(3))
	z3=-2.*fm(4)*real(wkv2(3))
	z4=2.*fm(5)*real(wkv2(4))
	z5=-2.*fm(6)*real(wkv2(4))
c		write(6,*)'fm(5),wkv2(4),z4=',fm(5),wkv2(4),z4
c		write(6,*)'fm(6),wkv2(4),z5=',fm(6),wkv2(4),z5
	dz1=fm(1)*real(wkv1(1))+fm(2)*real(wkv1(2))
	dz2=2.*fm(3)*real(wkv1(3))
	dz3=-2.*fm(4)*real(wkv1(3))
	dz4=2.*fm(5)*real(wkv1(4))
	dz5=-2.*fm(6)*real(wkv1(4))
	dz1h=fm(1)*real(wkv3(1))+fm(2)*real(wkv3(2))
	dz2h=2.*fm(3)*real(wkv3(3))
	dz3h=-2.*fm(4)*real(wkv3(3))
	dz4h=2.*fm(5)*real(wkv3(4))
	dz5h=-2.*fm(6)*real(wkv3(4))
	z1h=fm(1)*real(wkv4(1))+fm(2)*real(wkv4(2))
	z2h=2.*fm(3)*real(wkv4(3))
	z3h=-2.*fm(4)*real(wkv4(3))
	z4h=2.*fm(5)*real(wkv4(4))
	z5h=-2.*fm(6)*real(wkv4(4))
	z1o=fm(1)*real(wb1(1))+fm(2)*real(wb1(2))
	z2o=2.*fm(3)*real(wb1(3))
	z3o=-2.*fm(4)*real(wb1(3))
	z4o=2.*fm(5)*real(wb1(4))
	z5o=-2.*fm(6)*real(wb1(4))
c * *
c	  write(6,*)'fm(6),wkv4(4),z5h=',fm(6),wkv4(4),z5h
	x1=x1s(lr)
	x2=x2s(lr)
	x3=x3s(lr)
	dx1=dx1s(lr)
	dx2=dx2s(lr)
	dx3=dx3s(lr)
	ddx1=ddx1s(lr)
	ddx2=ddx2s(lr)
	ddx3=ddx3s(lr)
c 
	rd1g(im,idip,i)=rd1g(im,idip,i)+(z1o*x1 +z2o*x2*cp
     &	+z3o*x2*sp +z4o*x3*cp2 +z5o*x3*sp2)*dfacg
	rd1S(im,idip,i)=rd1S(im,idip,i)+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac
	r1hS(im,idip,i)=r1hS(im,idip,i)+(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac
c		write(6,*)'r1hS(',im,idip,i,')=',r1hS(im,idip,i)
	r1tS(im,idip,i)=r1tS(im,idip,i)+(-z2h*x2*sp
     &	+z3h*x2*cp -2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac/slf 
	rttS(im,idip,i)=rttS(im,idip,i)+(z1h*ddx1 +z2h*ddx2*cp
     &	+z3h*ddx2*sp +z4h*ddx3*cp2 +z5h*ddx3*sp2)*dfac1/bigr 
	rttS(im,idip,i)=rttS(im,idip,i)+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac1/bigr 
	rppS(im,idip,i)=rppS(im,idip,i)+(-z2h*x2*cp -z3h*x2*sp 
     &	-4.*z4h*x3*cp2 -4.*z5h*x3*sp2)*(dfac1/slf)/(slf*bigr) 
	rppS(im,idip,i)=rppS(im,idip,i)+(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac1*cotd/bigr 
	rppS(im,idip,i)=rppS(im,idip,i)+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac1/bigr 
	rtrS(im,idip,i)=rtrS(im,idip,i)+(dz1h*dx1 +dz2h*dx2*cp
     &	+dz3h*dx2*sp +dz4h*dx3*cp2 +dz5h*dx3*sp2)*dfac1/(2.*bigr) 
	rtrS(im,idip,i)=rtrS(im,idip,i)-(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac1/(2.*bigr) 
	rtrS(im,idip,i)=rtrS(im,idip,i)+(z1*dx1 +z2*dx2*cp
     &	+z3*dx2*sp +z4*dx3*cp2 +z5*dx3*sp2)*dfac1/(2.*bigr) 
	rtpS(im,idip,i)=rtpS(im,idip,i)+(-z2h*dx2*sp
     &	+z3h*dx2*cp -2.*z4h*dx3*sp2 +2.*z5h*dx3*cp2)*dfac1/(2.*bigr*slf) 
	rtpS(im,idip,i)=rtpS(im,idip,i)-(-z2h*x2*sp +z3h*x2*cp 
     &	-2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac1*cotd/(bigr*slf) 
	rtpS(im,idip,i)=rtpS(im,idip,i)+(-z2h*dx2*sp +z3h*dx2*cp 
     &	-2.*z4h*dx3*sp2 +2.*z5h*dx3*cp2)*dfac1/(2.*slf*bigr)
	rprS(im,idip,i)=rprS(im,idip,i)+(-z2*x2*sp
     &	+z3*x2*cp -2.*z4*x3*sp2 +2.*z5*x3*cp2)*dfac1/(2.*slf*bigr) 
	rprS(im,idip,i)=rprS(im,idip,i)+(-dz2h*x2*sp +dz3h*x2*cp 
     &	-2.*dz4h*x3*sp2 +2.*dz5h*x3*cp2)*dfac1/(2.*bigr*slf) 
	rprS(im,idip,i)=rprS(im,idip,i)-(-z2h*x2*sp +z3h*x2*cp 
     &	-2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac1/(2.*bigr*slf) 
c		if(l.gt.7490) write(6,*)'rprS(',im,idip,i,')=',rprS(im,idip,i)
	rrrS(im,idip,i)=rrrS(im,idip,i)+(dz1*x1 +dz2*x2*cp
     &	+dz3*x2*sp +dz4*x3*cp2 +dz5*x3*sp2)*dfac1/bigr 
c	Done with Spheroidal motion component of degree l.
c	TOROIDAL MODES
	do 475 idec=6,7
	kdip=(nmesh0-1)*idec-(nmesh0-1)+idip
	wkv1(idec)=dble(y1(kdip,lr))
	wkv2(idec)=dble(y2(kdip,lr))
475	continue
	z1=-2.*fm(4)*real(wkv2(6))
	z2=-2.*fm(3)*real(wkv2(6))
	z3=-2.*fm(6)*real(wkv2(7))
	z4=-2.*fm(5)*real(wkv2(7))
	dz1=-2.*fm(4)*real(wkv1(6)+wkv2(6))
	dz2=-2.*fm(3)*real(wkv1(6)+wkv2(6))
	dz3=-2.*fm(6)*real(wkv1(7)+wkv2(7))
	dz4=-2.*fm(5)*real(wkv1(7)+wkv2(7))
	romT(im,idip,i)=romT(im,idip,i)-facl*facl*((-z1*cp-z2*sp)*x2
     &	-(z3*cp2+z4*sp2)*x3)*dfac/(2.*bigr)
	r1tT(im,idip,i)=r1tT(im,idip,i)+(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac
	r1hT(im,idip,i)=r1hT(im,idip,i)-(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac/slf 
	rttT(im,idip,i)=rttT(im,idip,i)-(-z1*dx2*sp +z2*dx2*cp
     &	-2.*z3*dx3*sp2 +2.*z4*dx3*cp2)*dfac1/(bigr*slf) 
	rttT(im,idip,i)=rttT(im,idip,i)+(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1*cotd/(bigr*slf) 
	rppT(im,idip,i)=rppT(im,idip,i)+(-z1*dx2*sp +z2*dx2*cp 
     &	-2.*z3*dx3*sp2 +2.*z4*dx3*cp2)*dfac1/(bigr*slf)
	rppT(im,idip,i)=rppT(im,idip,i)-(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1*cotd/(bigr*slf) 
	rtpT(im,idip,i)=rtpT(im,idip,i)+(z1*ddx2*cp +z2*ddx2*sp
     &	+z3*ddx3*cp2 +z4*ddx3*sp2)*dfac1/(2.*bigr)
	rtpT(im,idip,i)=rtpT(im,idip,i)-(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac1*cotd/(2.*bigr)
	rtpT(im,idip,i)=rtpT(im,idip,i)-(-z1*x2*cp -z2*x2*sp
     &	-4.*z3*x3*cp2 -4.*z4*x3*sp2)*(dfac1/slf)/(2.*bigr*slf)
	rtrT(im,idip,i)=rtrT(im,idip,i)-(-dz1*x2*sp +dz2*x2*cp
     &	-2.*dz3*x3*sp2 +2.*dz4*x3*cp2)*dfac1/(2.*bigr*slf) 
	rtrT(im,idip,i)=rtrT(im,idip,i)+(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1/(2.*bigr*slf) 
	rprT(im,idip,i)=rprT(im,idip,i)+(dz1*dx2*cp +dz2*dx2*sp
     &	+dz3*dx3*cp2 +dz4*dx3*sp2)*dfac1/(2.*bigr)
	rprT(im,idip,i)=rprT(im,idip,i)-(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac1/(2.*bigr)
c	Done with Toroidal motion component of degree l.
371	continue
106	continue
107	continue
105	continue
	write(6,*)'response Greens functions determined'
	write(6,*)'rd1S(1,1,35)=',rd1S(1,1,35)
	write(6,*)'rd1S(1,8,35)=',rd1S(1,8,35)
	write(6,*)'rd1S(1,15,35)=',rd1S(1,15,35)
c	pause
c *****
c	Set up spline interpolation arrays for the response functions.
	do 111 im=1,6
	write(6,*)'setting up spline interpolation arrays for response functions'
	write(6,*)'im=',im
	do 112 idip=1,nmesh0-1
	write(6,*)'idip=',idip
	do 113 i=1,40
	yd1g(i)=dble(rd1g(im,idip,i))
	yd1S(i)=dble(rd1S(im,idip,i))
	y1hS(i)=dble(r1hS(im,idip,i))
	y1tS(i)=dble(r1tS(im,idip,i))
	yttS(i)=dble(rttS(im,idip,i))
	yppS(i)=dble(rppS(im,idip,i))
	ytrS(i)=dble(rtrS(im,idip,i))
	ytpS(i)=dble(rtpS(im,idip,i))
	yprS(i)=dble(rprS(im,idip,i))
	yrrS(i)=dble(rrrS(im,idip,i))
	yomT(i)=dble(romT(im,idip,i))
	y1hT(i)=dble(r1hT(im,idip,i))
	y1tT(i)=dble(r1tT(im,idip,i))
	yttT(i)=dble(rttT(im,idip,i))
	yppT(i)=dble(rppT(im,idip,i))
	ytrT(i)=dble(rtrT(im,idip,i))
	ytpT(i)=dble(rtpT(im,idip,i))
	yprT(i)=dble(rprT(im,idip,i))
113	continue
c		if(im.eq.1) write(6,*)'y1hS=',y1hS
c
	call splneq(40,yd1g,sarr)
	do i=1,40
	sd1g(im,idip,i)=sarr(i)
	enddo
	call splneq(40,yd1S,sarr)
	do i=1,40
	sd1S(im,idip,i)=sarr(i)
	enddo
	mu=real(evaleq(ypoint,ndep,ymu,smu))
	call splneq(40,y1hS,sarr)
	do i=1,40
	s1hS(im,idip,i)=sarr(i)
c		if(im.eq.1) write(6,*)'i=',i,'s1hS=',s1hS(im,idip,i)
	enddo
	call splneq(40,y1tS,sarr)
	do i=1,40
	s1tS(im,idip,i)=sarr(i)
	enddo
	call splneq(40,yttS,sarr)
	do i=1,40
	sttS(im,idip,i)=sarr(i)
	enddo
	call splneq(40,yppS,sarr)
	do i=1,40
	sppS(im,idip,i)=sarr(i)
	enddo
	call splneq(40,ytrS,sarr)
	do i=1,40
	strS(im,idip,i)=sarr(i)
	enddo
	call splneq(40,ytpS,sarr)
	do i=1,40
	stpS(im,idip,i)=sarr(i)
	enddo
	call splneq(40,yprS,sarr)
	do i=1,40
	sprS(im,idip,i)=sarr(i)
	enddo
	call splneq(40,yrrS,sarr)
	do i=1,40
	srrS(im,idip,i)=sarr(i)
	enddo
c
	call splneq(40,yomT,sarr)
	do i=1,40
	somT(im,idip,i)=sarr(i)
	enddo
	call splneq(40,y1hT,sarr)
	do i=1,40
	s1hT(im,idip,i)=sarr(i)
	enddo
	call splneq(40,y1tT,sarr)
	do i=1,40
	s1tT(im,idip,i)=sarr(i)
	enddo
	call splneq(40,yttT,sarr)
	do i=1,40
	sttT(im,idip,i)=sarr(i)
	enddo
	call splneq(40,yppT,sarr)
	do i=1,40
	sppT(im,idip,i)=sarr(i)
	enddo
	call splneq(40,ytrT,sarr)
	do i=1,40
	strT(im,idip,i)=sarr(i)
	enddo
	call splneq(40,ytpT,sarr)
	do i=1,40
	stpT(im,idip,i)=sarr(i)
	enddo
	call splneq(40,yprT,sarr)
	do i=1,40
	sprT(im,idip,i)=sarr(i)
	enddo
c
112	continue
111	continue
	write(6,*)'done determining spline interpolations of response functions'

	open(2,file='stat2g.out')
c	Begin sweep over observation points.
	do 60 ip=1,ipts 
	dispx=0.
	dispy=0.
	dispz=0.
	exx=0.
	eyy=0.
	exy=0.
	exz=0.
	eyz=0.
	ezz=0.
	omxy=0.
	dgrav=0.
c	Compute displacements, strains and gravity anomaly 
c	at points specified by
c	arrays olat and olon.  Begin summation over fault segments.
	jf=0
65	jf=jf+1
	if(jf.gt.iseg) go to 84
c		write(6,*)'jf=',jf
	write(6,*) flat(jf),flon(jf),fbigl(jf),fstr(jf)
	sstr=sin(fstr(jf))
	cstr=cos(fstr(jf))
	s2str=2.*sstr*cstr
	c2str=cstr*cstr-sstr*sstr 
	srak=sin(frake(jf))
	crak=cos(frake(jf))
	dlen=(fbigl(jf)/ethr)/real(nmesh2-1)
	dlon=olon(ip)-flon(jf)
c	Find angular distance and azimuth to observation point from
c	the initial point on the fault segment.
c	Note: delta is the angular distance from the earthquake.
c	Note: phi is the azimuth meazsured positive counterclockwise
c	from south.
	cdelt=cos(flat(jf))*cos(olat(ip))+sin(flat(jf))*
     &	sin(olat(ip))*cos(dlon)
	if (cdelt.le.0.9999) delta=acos(cdelt)
	if(cdelt.gt.0.9999) delta=sqrt((flat(jf)-olat(ip))**2+
     &	(dlon*sin(flat(jf)))**2)
	spsi=sin(dlon)*sin(olat(ip))/sin(delta)
	cpsi=(cos(olat(ip))-cos(flat(jf))*cdelt)/(sin(flat(jf))*sin(delta))
	phi=pi-atan2(spsi,cpsi) 
	dwrite=rad*delta
	pwrite=rad*phi
	  write(6,*)'delta,phi=',dwrite,'deg.',pwrite 
	cphi=cos(phi)
	sphi=sin(phi)
cNEW
	spsi=sin(dlon)*sin(flat(jf))/sin(delta)
	cpsi=(cos(flat(jf))-cos(olat(ip))*cdelt)/(sin(olat(ip))*sin(delta))
	peps=atan2(spsi,cpsi)
	cphio=cos(peps)
	sphio=sin(peps)
c--
	if(jf.gt.1) go to 59
	cphi1=cphi
	sphi1=sphi
	delt1=delta 
59	delta=delta*(6371./ethr)
c	Angular distances are a factor of (6371./ethr) larger on the
c	earth with radius ethr km.
	u1=fwt(jf)
	  write(6,*)'after 59, delta=',delta
c	u1=magnitude of slip.
c	Integrate displacements over fault elements.
	bigh=-dh/2.
	idip=0
81	idip=idip+1
	if(idip.eq.nmesh0) go to 65
c	DO DEPTH INTEGRAL.
	  write(6,*)'idip=',idip
	bigh=bigh+dh
	dep=dble(rmin+bigh/dip)
	ypoint=1.d0+(dep-rad1)*dble(ndep-1)/(rad10-rad1)
	  write(6,*)'dep,rad1,rad10=',dep,rad1,rad10
	  write(6,*)'ypoint=',ypoint
c ***	Determine moment tensor elements for fault element (ilen,idip).
c	Units of moment tensor are 10**20 N-m.
c	First get [mu] and [lam] for fault element at radius [dep].
	mu=real(evaleq(ypoint,ndep,ymu,smu))
	lam=real(evaleq(ypoint,ndep,ylam,slam))
	  write(6,*)'mu,lam=',mu,lam
	amesh=real((nmesh1-1)*(nmesh2-1))
c	Note that shear modulus [mu] is already input in units of 10**10 Pa.
c	Moment tensor from Ben Menahem and Singh, eqn. 4.115b for
c	shear dislocation, and derived from eqn. (4.101), (4.110), and
c	(4.113) for a tensile dislocation.
	if(abs(frake(jf)).gt.3.159046) go to 49
c	Next line is shear moment.
	shrm1=u1*fbigl(jf)*((dmax-dmin)/sdip)*mu*(1.e-6)/amesh
	p1=srak*sdip*cdip*s2str+crak*sdip*c2str
	p2=crak*sdip*s2str-srak*sdip*cdip*c2str
	p3=-crak*cdip*sstr+srak*c2dip*cstr
	p4=srak*c2dip*sstr+crak*cdip*cstr
	p5=srak*sdip*cdip
	mrr=shrm1*2.*p5
	mtt=-shrm1*(p2+p5)
	mpp=shrm1*(p2-p5)
	mrt=-shrm1*p4
	mrp=-shrm1*p3
	mtp=-shrm1*p1
	go to 51
49	srak=1.
	if(frake(jf).lt.-3.159046) srak=-1.
	  if(srak.eq.1.0) write(6,*)'tensile dislocation'
	  if(srak.eq.-1.0) write(6,*)'compressional dislocation'
        shrm1=2.*u1*fbigl(jf)*((dmax-dmin)/sdip)*srak*mu*(1.e-6)/amesh
	mrr=shrm1*cdip**2
	mtt=shrm1*(sdip*sstr)**2
	mpp=shrm1*(sdip*cstr)**2
	mrt=shrm1*sdip*cdip*sstr
	mrp=shrm1*sdip*cdip*cstr
	mtp=shrm1*sdip**2*sstr*cstr
	shrm2=0.5*shrm1*(lam/mu)
	mrr=mrr+shrm2
	mtt=mtt+shrm2
	mpp=mpp+shrm2
	if(abs(frake(jf)).lt.3.665191) go to 51
c	Now put in an isotropic component such that tr(M)=0.
	trm=mrr+mtt+mpp
	mrr=mrr-trm/3.
	mtt=mtt-trm/3.
	mpp=mpp-trm/3.
51	continue
	  write(6,*)'mrr,mtt,mpp,mrt,mrp,mtp='
	  write(6,*) mrr,mtt,mpp,mrt,mrp,mtp
c ***
c	Determine multiplying factors for the cases m=0,1,2 
	fm(1)=(1./(2.*mu*mu+3.*lam*mu))*((lam+mu)*mrr-(lam/2.)*(mtt+mpp))
	fm(2)=(1./(2.*mu*mu+3.*lam*mu))*(-lam*mrr+(lam/2.+mu)*(mtt+mpp))
	fm(3)=(1./(2.*mu))*mrt
	fm(4)=-(1./(2.*mu))*mrp
	fm(5)=(1./(2.*mu))*(mtt-mpp)
	fm(6)=-(1./mu)*mtp
c		fm(1)=0.
c		fm(2)=0.
c		fm(3)=0.
c		fm(4)=0.
c		fm(5)=0.
c		fm(6)=0.
	len=-dlen/2.
	ilen=0
83	ilen=ilen+1
	if(ilen.eq.nmesh2) go to 81
c	  write(6,*)'ilen=',ilen
	len=len+dlen
	bigx=delta*sphi+len*sstr+bigh*cstr
	bigy=-delta*cphi+len*cstr-bigh*sstr
	deltp=dble(sqrt(bigx*bigx+bigy*bigy))
	phip=atan2(bigy,bigx)+pi/2.	
	  deltpw=real(deltp)*rad
	  phipw=real(phip)*rad 
c	  write(6,*)'bigx,bigy,deltp,phip=',bigx,bigy,deltp,phipw 
	cp=cos(phip)
	sp=sin(phip)
	cp2=cos(2.*phip)
	sp2=sin(2.*phip)
	cl=dcos(deltp)
	sl=dsin(deltp)
	cotd=real(cl/sl)
	slf=real(sl)
	deltpf=real(deltp)
c*****
c	Now interpolate response functions at this deltpf.
	ypoint=dble(40.*log(deltpf*ethr)/log(dstmax))
	do 114 im=1,6
	do i=1,40
	sarr(i)=sd1g(im,idip,i)
	yd1g(i)=dble(rd1g(im,idip,i))
	enddo
	vd1g(im)=real(evaleq(ypoint,40,yd1g,sarr))
	do i=1,40
	sarr(i)=sd1S(im,idip,i)
	yd1S(i)=dble(rd1S(im,idip,i))
	enddo
	vd1S(im)=real(evaleq(ypoint,40,yd1S,sarr))
	do i=1,40
	sarr(i)=s1hS(im,idip,i)
	y1hS(i)=dble(r1hS(im,idip,i))
	enddo
	v1hS(im)=real(evaleq(ypoint,40,y1hS,sarr))
	do i=1,40
	sarr(i)=s1tS(im,idip,i)
	y1tS(i)=dble(r1tS(im,idip,i))
	enddo
	v1tS(im)=real(evaleq(ypoint,40,y1tS,sarr))
	do i=1,40
	sarr(i)=sttS(im,idip,i)
	yttS(i)=dble(rttS(im,idip,i))
	enddo
	vttS(im)=real(evaleq(ypoint,40,yttS,sarr))
	do i=1,40
	sarr(i)=sppS(im,idip,i)
	yppS(i)=dble(rppS(im,idip,i))
	enddo
	vppS(im)=real(evaleq(ypoint,40,yppS,sarr))
	do i=1,40
	sarr(i)=strS(im,idip,i)
	ytrS(i)=dble(rtrS(im,idip,i))
	enddo
	vtrS(im)=real(evaleq(ypoint,40,ytrS,sarr))
	do i=1,40
	sarr(i)=stpS(im,idip,i)
	ytpS(i)=dble(rtpS(im,idip,i))
	enddo
	vtpS(im)=real(evaleq(ypoint,40,ytpS,sarr))
	do i=1,40
	sarr(i)=sprS(im,idip,i)
	yprS(i)=dble(rprS(im,idip,i))
	enddo
	vprS(im)=real(evaleq(ypoint,40,yprS,sarr))
	do i=1,40
	sarr(i)=srrS(im,idip,i)
	yrrS(i)=dble(rrrS(im,idip,i))
	enddo
	vrrS(im)=real(evaleq(ypoint,40,yrrS,sarr))
c
	do i=1,40
	sarr(i)=somT(im,idip,i)
	yomT(i)=dble(romT(im,idip,i))
	enddo
	vomT(im)=real(evaleq(ypoint,40,yomT,sarr))
	do i=1,40
	sarr(i)=s1hT(im,idip,i)
	y1hT(i)=dble(r1hT(im,idip,i))
	enddo
	v1hT(im)=real(evaleq(ypoint,40,y1hT,sarr))
	do i=1,40
	sarr(i)=s1tT(im,idip,i)
	y1tT(i)=dble(r1tT(im,idip,i))
	enddo
	v1tT(im)=real(evaleq(ypoint,40,y1tT,sarr))
	do i=1,40
	sarr(i)=sttT(im,idip,i)
	yttT(i)=dble(rttT(im,idip,i))
	enddo
	vttT(im)=real(evaleq(ypoint,40,yttT,sarr))
	do i=1,40
	sarr(i)=sppT(im,idip,i)
	yppT(i)=dble(rppT(im,idip,i))
	enddo
	vppT(im)=real(evaleq(ypoint,40,yppT,sarr))
	do i=1,40
	sarr(i)=strT(im,idip,i)
	ytrT(i)=dble(rtrT(im,idip,i))
	enddo
	vtrT(im)=real(evaleq(ypoint,40,ytrT,sarr))
	do i=1,40
	sarr(i)=stpT(im,idip,i)
	ytpT(i)=dble(rtpT(im,idip,i))
	enddo
	vtpT(im)=real(evaleq(ypoint,40,ytpT,sarr))
	do i=1,40
	sarr(i)=sprT(im,idip,i)
	yprT(i)=dble(rprT(im,idip,i))
	enddo
	vprT(im)=real(evaleq(ypoint,40,yprT,sarr))
114	continue
c*****
cOLD	  gamma1=sp
cOLD	  gamma2=-cp 
cNEW
	bigx=delta*sphio+len*sstr+bigh*cstr
	bigy=-delta*cphio+len*cstr-bigh*sstr
	deltp=dble(sqrt(bigx*bigx+bigy*bigy))
	phip=atan2(bigy,bigx)+pi/2.	
	  gamma1=sin(phip)
	  gamma2=-cos(phip)
c--
	dgrav1=vd1g(1)*fm(1)+vd1g(2)*fm(2)+vd1g(3)*fm(3)*cp
     &	+vd1g(4)*fm(4)*sp+vd1g(5)*fm(5)*cp2+vd1g(6)*fm(6)*sp2
	omxy1=vomT(4)*fm(4)*cp+vomT(3)*fm(3)*sp+vomT(6)*fm(6)*cp2
     &	+vomT(5)*fm(5)*sp2
	disp1=vd1S(1)*fm(1)+vd1S(2)*fm(2)+vd1S(3)*fm(3)*cp
     &	+vd1S(4)*fm(4)*sp+vd1S(5)*fm(5)*cp2+vd1S(6)*fm(6)*sp2
	disp1h=v1hS(1)*fm(1)+v1hS(2)*fm(2)+(v1hS(3)+v1hT(3))*fm(3)*cp
     &	+(v1hS(4)+v1hT(4))*fm(4)*sp+(v1hS(5)+v1hT(5))*fm(5)*cp2
     &	+(v1hS(6)+v1hT(6))*fm(6)*sp2
	disp1t=(v1tS(3)+v1tT(3))*fm(3)*sp
     &	+(v1tS(4)+v1tT(4))*fm(4)*cp+(v1tS(5)+v1tT(5))*fm(5)*sp2
     &	+(v1tS(6)+v1tT(6))*fm(6)*cp2
	ett=vttS(1)*fm(1)+vttS(2)*fm(2)+(vttS(3)+vttT(3))*fm(3)*cp
     &	+(vttS(4)+vttT(4))*fm(4)*sp+(vttS(5)+vttT(5))*fm(5)*cp2
     &	+(vttS(6)+vttT(6))*fm(6)*sp2
	epp=vppS(1)*fm(1)+vppS(2)*fm(2)+(vppS(3)+vppT(3))*fm(3)*cp
     &	+(vppS(4)+vppT(4))*fm(4)*sp+(vppS(5)+vppT(5))*fm(5)*cp2
     &	+(vppS(6)+vppT(6))*fm(6)*sp2
	etr=vtrS(1)*fm(1)+vtrS(2)*fm(2)+(vtrS(3)+vtrT(3))*fm(3)*cp
     &	+(vtrS(4)+vtrT(4))*fm(4)*sp+(vtrS(5)+vtrT(5))*fm(5)*cp2
     &	+(vtrS(6)+vtrT(6))*fm(6)*sp2
	etp=(vtpS(3)+vtpT(3))*fm(3)*sp
     &	+(vtpS(4)+vtpT(4))*fm(4)*cp+(vtpS(5)+vtpT(5))*fm(5)*sp2
     &	+(vtpS(6)+vtpT(6))*fm(6)*cp2
	epr=(vprS(3)+vprT(3))*fm(3)*sp
     &	+(vprS(4)+vprT(4))*fm(4)*cp+(vprS(5)+vprT(5))*fm(5)*sp2
     &	+(vprS(6)+vprT(6))*fm(6)*cp2
	err=vrrS(1)*fm(1)+vrrS(2)*fm(2)+vrrS(3)*fm(3)*cp
     &	+vrrS(4)*fm(4)*sp+vrrS(5)*fm(5)*cp2+vrrS(6)*fm(6)*sp2
c
c	Rotate strain tensor into x-y-z Cartesian coordinates.
c	Also rotate displacements into x-y coordinates.
	dispx1=gamma1*disp1h-gamma2*disp1t
	dispy1=gamma2*disp1h+gamma1*disp1t
	dispz1=disp1
	exx1=gamma1**2*ett+gamma2**2*epp-2.*gamma1*gamma2*etp
	exy1=gamma1*gamma2*(ett-epp)+(gamma1**2-gamma2**2)*etp
	eyy1=gamma2**2*ett+gamma1**2*epp+2.*gamma1*gamma2*etp
	exz1=etr*gamma1-epr*gamma2
	eyz1=etr*gamma2+epr*gamma1
	ezz1=err
70	dispx=dispx+dispx1
	dispy=dispy+dispy1
	dispz=dispz+dispz1
	exx=exx+exx1
	eyy=eyy+eyy1
	exy=exy+exy1
	exz=exz+exz1
	eyz=eyz+eyz1
	ezz=ezz+ezz1
	omxy=omxy+omxy1
	dgrav=dgrav+dgrav1
	go to 83
84	continue
c	Multiply everything by [lste] because of the stepped l-summation.
	dispx=dispx*flste
	dispy=dispy*flste
	dispz=dispz*flste
	exx=exx*flste
	eyy=eyy*flste
	ezz=ezz*flste
	exz=exz*flste
	eyz=eyz*flste
	exy=exy*flste
	dgrav=dgrav*flste
c	Find cartesian corrdinates of observation point relative
c	to fault segment #1.
	xc=(twopi*6371./360.)*rad*(delt1)*sphi1  
	yc=-(twopi*6371./360.)*rad*(delt1)*cphi1  
c	xc measured in east direction; yc measured in north direction (km).
c	displacements in cm.  gravity anomaly in mGal. 
	write(2,76) xc,yc,dispx,dispy,dispz,exx,eyy,exy,exz,eyz,ezz,omxy,dgrav  
	write(6,76) xc,yc,dispx,dispy,dispz,exx,eyy,exy,exz,eyz,ezz,omxy,dgrav
76	format(2f10.3,11e13.6e2)
77	format(2f10.3,2e13.6e2)
60	continue
	close(2)
	write(6,*)'end of stat1'
	end  
	 
