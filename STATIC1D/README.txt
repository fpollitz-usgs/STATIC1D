		Program Package STATIC1D

		by Fred F. Pollitz

These programs are an implementation of the Direct Green's Function method
for static deformation and free-air (or Bouguer) gravity anomaly
described by Pollitz (1996, 1997).  They solve the equations of
static equilibrium in a spherically layered isotropic medium using
a decomposition into spheroidal and toroidal motions.  For each spherical
harmonic degree l and azimuthal order number m, the (l,m) response function is
deternined subject to jumps in the displacement-stress vector at the source
radius, a zero-traction boundary condition at Earth's surface, and a homogeneous
isotropic elastic solid at the base of the specified Earth model.
The programs are very flexible, being suitable for calculations 
of the static displacement field ranging from local to global.

There are two main programs:

1) stat0A:      Computes response functions for spherical harmonic degrees
                from input l=(lmin) up to input l=[lmax], and azimuthal order 
                numbers m=0, 1, and 2.  The observation-pt depth and a range of
		source depths are fixed according to input.

2) stat2gA:     Computes static displacement at specified observation points
                and source geometry, assuming a distribution of point sources
                along a plane specified by input parameters.
                The observation-pt depth for static displacements
		is that specified in stat0A, and the gravity anomaly
		is evaluated at Earth's surface.

---------------------------------------------------------------------
---------------------------------------------------------------------

		INPUT AND OUTPUT CONVENTIONS


STAT0A evaluates Greens functions for the static deformation at a certain depth
resulting from shear dislocations that may exist in a certain depth interval

---Input format to STAT0A---

nice stat0A << ! > /dev/null
1 150
50.5 0.0
0.
!

Use spherical harmonic degrees l=1 through l=150.

Anticipate fault planes (to be specified more precisely in the
input to STAT2gA) that have a range of lower/upper edge depths
between 50.5 and 0.0.
Note that the minimum wavelength represented in the calculations is
approximately the circumference of the Earth divided by the maximum
spherical harmonic degree.  In this case it is 40000km/150 ~ 267 km.
It is possible to go up to much higher maximum spherical harmonic degree --
I routinely go up to a maximum of 5000 or higher.

Evaluate the resulting displacement and free-air gravity fields at 0. km depth.
This can be changed to other depth values (e.g. use 25. if deformation
at depth 25 km is desired)

---------------------------------------------------------------------

STAT2gA evauates the static deformation and free-air (or Bouguer) 
gravity anomaly resulting
from slip on an input fault plane, at a set of points specified
in the input file to STAT2gA.  The observation depth is the last number 
specified in the input to STAT0A

---Input format to STAT2gA---

For example, the first few lines of stat2g.in_SUM1Dgr are

50.  30.0  35.
1
 9.7938041 92.9885300  162.5 350.0 105. 1000.
1600
     0.000000    87.000000
     0.000000    87.461540
     0.000000    87.923080
     0.000000    88.384613

The lower and upper edge depths are 50 km and 30 km, and the fault dip is
35 deg.  The fault length and strike are 162.5 km and 350 deg. (measured
clockwise from due North).  The rake of slip is 105 deg. (mostly dip slip
with small component of right-lateral slip).  The slip is 1000 cm.
The point on the lower edge of the fault plane closest to the strike
direction is lat=9.7938041 deg.N, long=92.9885300 deg.E., i.e. since
the strike is 350 deg. this is the
northernmost point on the lower edge of the fault.

The static deformation is evaluated at 1600 points at the lat,lon's
given in the list in this input file.  The output file 'stat2g.out'
has lines in the format

	write(2,76) xc,yc,dispx,dispy,dispz,exx,eyy,exy,exz,eyz,ezz,omxy,dgrav  
76	format(2f10.3,11e13.6e2)

xc,yc are approximations to the Cartesian coordinates of the observation pt.,
relative to the lower edge point of the fault.  These numbers are not
important.  The actual calculation is done in a spherical geometry, and
the lat,lon-locations of the observation points are already known
from the input file to STAT2gA

dispx,dispy,dispz are static displacement in cm measured in the local
due East, North, and Up directions, respectively.

With x=local due East, y=local due North, z=Up, there follows
six strain components in dimensionless units of microstrain.

omxy is the rotation = (1/2)*(du_x/dy - du_y/dx) in dimensionless
units of microstrain

dgrav is the free-air gravity anomaly in milligals.

---------------------------------------------------------------------

		NOTES

1) If the Bouguer gravity anomaly is desired, then uncomment
the cBOUGUER line in STAT0A.f and re-run the example. The density
jump at Earth's surface assumes that the solid Earth is overlain by
water.  If that is not the case (and this is an issue only for the
Bouguer gravity anomaly), then replace the line
	if(i.eq.n) ddens=1.000d0-dble(dens(i))
with
	if(i.eq.n) ddens=-dble(dens(i))

2) A fluid-solid boundary condition is implemented at the base
of the input earth model, i.e. the sphere just below the lowest layer
in 'earth.model_stat' is assumed to be an incompressible fluid.  
It is also possible to implement a solid-solid boundary condition.
In this case, the following lines in STAT0A.f should be commented out:

	b1(1)=dble(l)
	b1(2)=0.d0
	b1(3)=1.d0
	b1(4)=0.d0
	b2(1)=0.d0
	b2(2)=0.d0
	b2(3)=1.d0
	b2(4)=0.d0

	t1(1)=1.d0
	t1(2)=0.d0

3) To avoid ringing effects related to the cutoff at spherical harmonic
degree lmax, a cosine taper is applied over the l-interval from
lmax/2 to lmax.  To turn off this tapering, the following line in STAT2gA.f
should be commented out:

	if(l.ge.iw0x) facf=1.0-cos((flmax-fl)/(fiw0m-1.) * real(pi)/2.)

---------------------------------------------------------------------
---------------------------------------------------------------------

		EXAMPLE 1

o Compile with 'make all'

o Run stat2g.xEX1 to compute coseismic deformation and 
	gravity anomaly, with maximum spherical harmonic degree of 
	150, 100, 200, and 75 in separate computations.
	The last number in every line of the output files 
	stat2g.outCOSEISMIC150
	stat2g.outCOSEISMIC100
	stat2g.outCOSEISMIC200
	stat2g.outCOSEISMIC75
	has gravity anomaly in millgals.

o Run getg.xCO to extract vertical displacement and gravity anomaly from
	stat2g.outCOSEISMIC150
	stat2g.outCOSEISMIC100
	stat2g.outCOSEISMIC200
	stat2g.outCOSEISMIC75
	Since GETG multiplies by 1000,
	the resulting gravity anomaly is in microgals.
	These values -- in
	getg.outCOSEISMIC200
	getg.outCOSEISMIC150
	getg.outCOSEISMIC100
	getg.outCOSEISMIC75
	-- are in GMT format.
	The resulting vertical displacement is in cm.
	These values -- in
	getz.outCOSEISMIC200
	getz.outCOSEISMIC150
	getz.outCOSEISMIC100
	getz.outCOSEISMIC75
	-- are in GMT format.

	'getzCOSEISMIC.eps' is a plot of the vertical uplift and
	free-air gravity patterns for l=1 up to various lmax.
	It is identical to Figure 1 of Pollitz (2006).

---------------------------------------------------------------

		EXAMPLE 2

We wish to re-run the previous calculation of 2004 Sumatra coseismic
displacement field and gravity anomaly from l=1 to 200, at Earth's surface,
at a different set of observation points in a global calculation.  To avoid  
revising the observation-pt lat,lon list of every STAT2gA-input file, we allow 
STAT2gA to read in the observation-pt lat,lons from a file 'latlon.inDEF'

o In STAT2gA, uncomment the lines

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

	and re-compile STAT2gA

o Run stat2g.xEX2

	This accomplishes the same task as EXAMPLE 1 for the 
	maximum spherical harmonic degree of 200, but it first
	copies 'latlon.inworld' into 'latlon.inDEF', from which
	stag2gA takes its input of observation-pt lat,lons.
	After computing the global displacement field
	(moving the output file to 'stat2g.out_world'), program
	DISPH extracts the horizontal displacement field, 
	multiplies by 10 to convert from cm to mm,
	and the resulting output file 'disph.gmt_sumatra_co_l=1-200'
	is in GMT format (lon,lat,u_E(mm),u_N(mm)).

	A plot of the global coseismic displacement field
	is in 'disphCOSEISMIC.eps'

------------------------------------------------------------------------

		REFERENCES

Pollitz, F.F., Coseismic deformation from earthquake faulting on a layered 
spherical Earth, Geophys. J. Int., 125, 1-14 (1996). 

Pollitz, F.F., Gravity anomaly from faulting on a layered spherical Earth with 
application to central Japan, Phys. Earth Planet. Int., 99, 259-271 (1997). 

Pollitz, F.F., A New Class of Earthquake Observations, Science, 313, 619-620 (2006).

