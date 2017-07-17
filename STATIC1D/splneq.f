      subroutine splneq(nn, u, s)
c$$$$$ calls no other routines
c  finds coeffs for a spline interpolation of equally spaced data
c  based on 'spline interpolation ona  digital computer' by r.f.thompson
c  nn  number of data points (may be negative - see d1,d2 below)
c  u  array of function values to be interpolated, assumed to samples at equal i
c     intervals of a twice differentiable function.
c  s  array to hold computed values of 2nd derivative of spline fit at sample
c     points.  these values are required by evaleq  to interpolate
c  if the user wishes to force specific values on the derivatives at the end
c  points, he should put h*du(1)/dx  nad  h*du(n)/dx  in  s(1),s(2), then call
c  splneq  with nn=-number of terms in series. h = sample spacing in x.
c  normally the derivatives are found by fitting a parabola through the
c  1st and last 3 points.
c  if the number of terms is between 1 and 3, straight-line interpolation is don
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension u(1),s(1),a(13)
c
      n=iabs(nn)
      if (n.le.3) go to 5000
      d1=-0.5d0*u(3)  +2.d0*u(2)  -1.5d0*u(1)
      dn= 0.5d0*u(n-2)-2.d0*u(n-1)+1.5d0*u(n)
      if (nn.gt.0) go to 1000
      d1=s(1)
      dn=s(2)
 1000 a(1)=2.0d0
      a(2)=3.5d0
      s(1)=u(2)-u(1)-d1
      s(2)=u(1)-2.d0*u(2)+u(3)-0.5d0*s(1)
      n1=n-1
      do 3000 i=3,n1
      if (i.gt.13) go to 3000
      k=i
      a(k)=4.d0-1.d0/a(k-1)
 3000 s(i)=u(i-1)-2.d0*u(i)+u(i+1)-s(i-1)/a(k-1)
      s(n)=u(n1)-u(n)+dn-s(n1)/a(k)
      s(n)=6.d0*s(n)/(2.d0-1.d0/a(k))
      n2=n-2
c  compute 2nd derivatives by back-substitution
c  the array  a  tends to a constant (2+sqrt(3)) so only 13 elements are needed
      do 4000 j=1,n2
      i=n-j
      k=min0(i,k)
 4000 s(i)=(6.d0*s(i)-s(i+1))/a(k)
      s(1)=3.d0*s(1)-0.5d0*s(2)
      return
c  series too short for cubic spline.  fit straight lines.
 5000 do 5500 i=1,n
 5500 s(i)=0.d0
      return
      end  
      function evaleq(y, nn, u, s)
c$$$$$ calls no other routines
c  performs spline interpolation of equally spaced data.
c  based on 'spline interpolation on a digital computer'( by r.f.thompson.
c  evaluates a spline interpolate in a set of equally spaced samples.
c  the routine  splneq  should be called first, to establish the array  s .
c  y  the  coordinate at which interpolate is required, with y=1 for 1st
c     sample point, y=2 for the 2nd, etc.  if actual spacing is  h  and  x1 is
c     the 1st sample coordinate use  y = 1.0 + (x-x1)/h
c  nn  number of samples of function in original set.
c  u  array containing function samples.
c  s  array of normalized 2nd derivatives, computed by  splneq.  the derivatives
c     have been multiplied by h**2, where h is the sample spacing.
c  if  y  is out of the range (1,nn), the 1st or last sample value is used.
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension u(1),s(1)
c
      if (y.le.1.d0) go to 1500
      if (y.ge.dble(real(nn))) go to 2000
      k1=y
      k=k1+1
      dk=k-y
      dk1=y-k1
      ff1=s(k1)*dk*dk*dk
      ff2=s(k)*dk1*dk1*dk1
      evaleq=(dk*(6.d0*u(k1)-s(k1))+dk1*(u(k)*6.d0-s(k))+ff1+ff2)/6.d0
      return
c  out of range.  supply constant values
 1500 evaleq=u(1)
      return
 2000 evaleq=u(nn)
      return
      end
