c	NOTE: version of subroutines lgndr0, lgndr1, lgndr2 in which
c	upward recurrence is carried out only up to lmax+1 rather
c	than automatically to 25001.
c
	subroutine lgndr0(ang,x0,dx0)
c	Returns the Legendre function X sub n sup 0 (cos(ang))
c	for n=1 through lmax+1. (ang in radians).
c	Array x0 contains X sub n sup 0
c	Array dx0 contains (d/d theta) X sub n sup 0 (x=cos(theta)).
	real*8 x,y,fac(25001),tox,bj,bjm,bjp,dang,x0(1),dx0(1)
	common/maxl/lmax
	dang=dble(ang)
	x=dcos(dang)
	y=dsin(dang)
	do 5 n=1,lmax+1  
	fac(n)=dsqrt((dble(n)+0.5d0)/(2.d0*3.1415926535d0))
5	continue 
	if(x.lt.-1.0d0.or.x.gt.1.0d0) pause 
     &	'x must be nonnegative in lgend'
		tox=2.d0*x
		bjm=1.d0
		bj=x
		x0(1)=x*fac(1)  
		n=lmax+1
		do 11 j=1,n-1
		  bjp=((dble(j)+0.5d0)*tox*bj-dble(j)*bjm)/dble(j+1)
		  bjm=bj
		  bj=bjp
		  x0(j+1)=bj*fac(j+1)
	  dx0(j)=-(dble(j+1)*x*bjm-dble(j+1)*bj)*fac(j)/y  
11		continue
	return
	end  
	subroutine lgndr1(ang,x1,dx1)
c	Returns the Legendre function X sub n sup 1 (cos(ang))
c	for n=1 through lmax+1. (ang in radians).
c	Array x1 contains X sub n sup 1
c	Array dx1 contains (d/d theta) X sub n sup 1 (x=cos(theta)).
	real*8 x,y,fac(25001),tox,bj,bjm,bjp,dang,x1(1),dx1(1)
	common/maxl/lmax
	dang=dble(ang)
	x=dcos(dang)
	y=dsin(dang)
	do 5 n=1,lmax+1  
	fac(n)=dsqrt((dble(n)+0.5d0)/(2.d0*3.1415926535d0))
	fac(n)=-fac(n)*dsqrt(1.d0-x*x)/dsqrt(dble(n*(n+1)))
5	continue 
	if(x.lt.-1.0d0.or.x.gt.1.0d0) pause 
     &	'x must be nonnegative in lgend'
		tox=2.d0*x
		bjm=1.d0
		bj=3.d0*x
		x1(1)=bjm*fac(1)
		x1(2)=bj*fac(2)  
		n=lmax+1
		do 11 j=1,n-1
	        if(j.eq.1) go to 12
		  bjp=((dble(j)+0.5d0)*tox*bj-dble(j+1)*bjm)/dble(j)
		  bjm=bj
		  bj=bjp
		  x1(j+1)=bj*fac(j+1)
12	  dx1(j)=-(dble(j+1)*x*bjm-dble(j)*bj)*fac(j)/y  
11		continue
	return
	end  
	subroutine lgndr2(ang,x2,dx2)
c	Returns the Legendre function X sub n sup 2 (cos(ang))
c	for n=1 through lmax+1. (ang in radians).
c	Array x1 contains X sub n sup 2
c	Array dx1 contains (d/d theta) X sub n sup 2 (x=cos(theta)).
	real*8 x,y,fac(25001),tox,bj,bjm,bjp,dang,x2(1),dx2(1)
	common/maxl/lmax
	dang=dble(ang)
	x=dcos(dang)
	y=dsin(dang)
	do 5 n=2,lmax+1  
	fac(n)=dsqrt((dble(n)+0.5d0)/(2.d0*3.1415926535d0))
	fac(n)=fac(n)*(1.d0-x*x)/dsqrt(dble(n*(n+1)))
	fac(n)=fac(n)/dsqrt(dble((n-1)*(n+2)))
5	continue 
	if(x.lt.-1.0d0.or.x.gt.1.0d0) pause 
     &	'x must be nonnegative in lgend'
		tox=2.d0*x
		bjm=3.d0
		bj=15.d0*x
		x2(1)=0.d0
		x2(2)=bjm*fac(2)
		x2(3)=bj*fac(3)  
		n=lmax+1
		do 11 j=1,n-1
		if(j.eq.1) go to 11
	        if(j.eq.2) go to 12
		  bjp=((dble(j)+0.5d0)*tox*bj-dble(j+2)*bjm)/dble(j-1)
		  bjm=bj
		  bj=bjp
		  x2(j+1)=bj*fac(j+1)
12	  dx2(j)=-(dble(j+1)*x*bjm-dble(j-1)*bj)*fac(j)/y  
11		continue
	return
	end  
