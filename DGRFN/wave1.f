c	Program WAVE1
c***
c	Determine seismic velocity time series using the
c	response functions computed by WAVE0S and stored in 'wave0.out'.
c	The source depth and observation-pt depth are previously specified
c	in the input to WAVE0S.  WAVE1 essentially supplies the source
c	information and writes out the time series at a number of observation pts.
c***	INPUT
c	The corner period is [corper] sec.  The seismic displacement spectrum
c	is cosine-tapered in the frequency domain between f and 2f,
c	where f=1/corper.
c	This program assumes a source extending along a horizontal line segment.  
c	It handles dipping faults with any shear dislocation.
c	The source location (more precisely, the endpoint of the line segment
c	closest to the strike direction) is [flat],[flon].
c	The seismic moment tensor is computed for a line-segment souce with
c	input parameters 
c	fault length ([fbigl]), dip ([dip]), rake ([frake]), strike [fstr],
c	and scalar moment ([fmom]), then collapsed onto a single source depth.
c	That is, this moment tensor is uniformly distributed along a horizontal 
c	line segment of a given length ([fbigl]).
c	A number of such line segments [iseg], each with different strike, etc.,
c	may be read in.  
c	The dip of each fault segment is fixed, i.e.
c	all fault segments share the same [dip].
c	The dip, strike, and rake follow
c	the convention of Ben Menahem and Singh...
c	There are [ipts] observation points with locations [olat(i)],[olon(i)],
c	for i=1,...,ipts.
c***	OUTPUT
c	Three-component velocity time series written out to 
c	'wave1.outx', 'wave1.outy', and 'wave1.outz'
c	where x=East, y=North, z=Up.
c	Units are cm/sec.
c***
	character*3 comp
c	Output strains in units of 10**(-6) .
	real*4 mrr,mtt,mpp,mrt,mrp,mtp  
	parameter (maxsrc=30)
	real*4 mrr0(maxsrc),mtt0(maxsrc),mpp0(maxsrc),mrt0(maxsrc),mrp0(maxsrc),mtp0(maxsrc)  
	real*4 kappa,mu,lam
cOLD	common/source/mrr,mtt,mpp,mrt,mrp,mtp,depth,slat,slon
	real*8 cl,sl,xl0(25001),dxl0(25001),xl1(25001),dxl1(25001)
	real*8 xl2(25001),dxl2(25001)
c	parameter (maxlen=2048/3)
	parameter (maxlen=2*256/3)
c	parameter (maxlen=2*700/3)
c	parameter (maxlen=1300)
c	parameter (maxlen=590)
	parameter (iommax=4096)
cUSE	double complex y1(7*maxlen,25000)
cUSE	double complex y2(7*maxlen,25000)
cUSE	double complex y3(5*maxlen,25000)
cUSE	double complex y4(5*maxlen,25000)
	double complex y1(7*maxlen,11500)
	double complex y2(7*maxlen,11500)
	double complex y3(5*maxlen,11500)
	double complex y4(5*maxlen,11500)
	double complex bfacl,bfacm 
	double complex dfac,dfac1
	real*8 fm(6)
	real*8 x1s(25000),x2s(25000),x3s(25000)
	real*8 dx1s(25000),dx2s(25000),dx3s(25000)
	real*8 ddx1s(25000),ddx2s(25000),ddx3s(25000)
	real*8 x1,x2,x3,dx1,dx2,dx3,ddx1,ddx2,ddx3,faclk,cotd
	parameter (maxpts=3000)
	dimension olat(maxpts),olon(maxpts)
	dimension flat(maxsrc),flon(maxsrc)
     	dimension fbigl(maxsrc),fstr(maxsrc),fdip(maxsrc)
	dimension fmom(maxsrc),frake(maxsrc)
c	integer*4 sid(10)
	real*8 dep,odep
cDONT NEED	real*8 rads
	double complex wkv1(7),wkv2(7),wkv3(5),wkv4(5)
	double complex z1,z2,z3,z4,z5,dz1,dz2,dz3,dz4,dz5
	double complex z1h,z2h,z3h,z4h,z5h,dz1h,dz2h,dz3h,dz4h,dz5h
	double complex disp1,disp1h,disp1t
	double complex ett,epp,err,etp,etr,epr
	double complex dispx1,dispy1,dispz1,exx1,eyy1,ezz1,exy1,exz1,eyz1
	double complex dispx(maxlen),dispy(maxlen),dispz(maxlen)
	double complex exx(maxlen),eyy(maxlen),ezz(maxlen)
	double complex exy(maxlen),exz(maxlen),eyz(maxlen)
	dimension uexx(4*iommax,maxpts),uexy(4*iommax,maxpts),uexz(4*iommax,maxpts)
	dimension ueyy(4*iommax,maxpts),ueyz(4*iommax,maxpts),uezz(4*iommax,maxpts)
cOLD	real*8 ymu,ylam
	real*8 deltp
	real*4 len,kmin
	common/maxl/lmax
	complex*8 xo,oome,ome,ui
	common/spec/xo(iommax),xt(4*iommax)
	common/invtyp/ityp
c 
c	ityp=0 specifies displacement seismograms in INVFT subroutine
c	ityp=1 specifies velocity seismograms in INVFT subroutine
	ityp=1
cfor DGRFN-M8.6 wave1- calculations	ityp=0
c	ethr=earth's radius in km.
c	bigr=earth radius in 10**6 cm.
	pi=3.1415926535  
	twopi=2.*pi
	ui=cmplx(0.,1.0)
	rad=180./3.1415926
c	component.txt, if it exists, will have one line
c	rtz
c	if output seismograms as Radial, Transverse, Vertical
c	is desired instead of the default convention East, North, Vertical
	open(2,file='component.txt')
	rewind(2)
	read(2,5,end=10) comp
10	continue
	close(2)
5	format(a3)
	write(6,5) comp
c
	open(2,file='wave0.out',form='unformatted')
	rewind(2)
	read(2) ethr,sdep
	  write(6,*)'ethr=',ethr
	bigr=ethr/10.
	read(2) lmin,lmax
	  write(6,*)'lmin,lmax=',0,lmax
	kmin=twopi/(real(lmax)+0.5)
	read(2) odep
	read(2) kappa,mu
	  write(6,*)'kappa,mu=',kappa,mu
	lam=kappa-2.*mu/3.
c	ymu=dble(mu)
c	ylam=dble(lam)
c	  write(6,*)'ymu,ylam=',ymu,ylam
c *	dt and nleno must be the same as in WAVE0S
	read(2) dt,nleno
	ommin=(2.*pi/dt)/real(nleno)
	obeta=2.0*ommin
c--
	write(6,*)'Reading in spheroidal and toroidal motion depth functions'
	do 55 iom=1,nleno/3
c		write(6,*)'reading iom=',iom,'out of',nleno
	do 50 lr=lmin,lmax+1
	l=lr
c	Do l=0 case last.
	if(lr.eq.(lmax+1)) l=0
c		write(6,*)'l,lr=',l,lr
	do 69 idec=1,5
	k=(nleno/3)*idec-(nleno/3)+iom
	read(2) ldum,y1(k,lr),y2(k,lr),y3(k,lr),y4(k,lr)
	read(2) bfacl,bfacm
69	continue
	if(lr.eq.(lmax+1)) go to 50
	do 58 idec=6,7
	k=(nleno/3)*idec-(nleno/3)+iom
	read(2) ldum,y1(k,lr),y2(k,lr)
cTE
c	if(ldum.eq.883.or.ldum.eq.884) then
c		write(6,*)'TOR:', ldum,y1(k,lr),y2(k,lr)
c	endif
c--
58	continue
c	  pause
c	Note: y1 - y4(1 thru 80,l) are Spheroidal motion coeff.
c	      y1 - y2(81 thru 112,l) are Toroidal motion coeff.
50	continue
55	continue
	close(2)
cDONT NEED	rads=dble(1.-sdep/ethr)
	write(6,*)'Done reading in spheroidal and toroidal motion depth functions'
	write(6,*)'[corper] corner period (sec)?'
	read(5,*) corper
cTEST
c		lmax=2000
	write(6,*)'  '
	write(6,*)'  '
	write(6,*)'NOTE: The source depth is fixed at the value'
	write(6,*)'read in from WAVE0.OUT.'
	write(6,*)'The source plane is disstributed onto'
	write(6,*)'a line segment at a single depth'
c 
	write(6,*)'finite fault with [iseg] # segments.  iseg=?'
	read(5,*) iseg
c	iobs=0
	do 39 i=1,iseg
	write(6,*)'segment #',i,'lat,lon(deg.),length(km)'
	write(6,*) 'strike(deg.),dip(deg.),rake(deg.),moment (10^20 N m)?'
	read(5,*) flat(i),flon(i),fbigl(i),fstr(i),fdip(i),frake(i),fmom(i)
	flat(i)=(pi/2.-flat(i)/rad)
	flon(i)=flon(i)/rad
	fstr(i)=fstr(i)/rad 
	frake(i)=frake(i)/rad 
	if(fmom(i).lt.0.0) then
	write(6,*)'Read in moment tensor components'
	write(6,*)'mrr,mtt,mpp,mrt,mrp,mtp (10^20 N m)?'
	write(6,*)'where r=r hat=Up, t=theta hat=South, p=phi hat=East'
	read(5,*) mrr0(i),mtt0(i),mpp0(i),mrt0(i),mrp0(i),mtp0(i)
	endif
39	continue
c 
	write(6,*)'number of observation points [ipts]?'
	read(5,*) ipts
	write(6,*)'ipts=',ipts
	do 40 i=1,ipts 
c	write(6,*)'point #',i,'lat,lon(deg.)?'
	read(5,*) olat(i),olon(i)
	olat(i)=(pi/2.-olat(i)/rad)
	olon(i)=olon(i)/rad
40	continue
	write(6,*)'finished reading displacement coefficients'
c	Divide moment tensor by nmesh (# elements in fault subdivision).
	nmesh2=2
c	nmesh2-1=# horizontal length points 
c*****	
	open(2,file='wave1.outx')
	open(4,file='wave1.outy')
	open(8,file='wave1.outz')
c	Begin sweep over observation points.
	do 60 i=1,ipts 
		write(6,*)'Doing obs pt #',i,' out of ',ipts
	do 84 iom=1,nleno/3
		write(6,*)'Doing freq # ',iom,' out of ',nleno
	dispx(iom)=0.
	dispy(iom)=0.
	dispz(iom)=0.
	exx(iom)=0.
	eyy(iom)=0.
	exy(iom)=0.
	exz(iom)=0.
	eyz(iom)=0.
	ezz(iom)=0.
c	Compute displacements and strains
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
	cdip=cos(fdip(jf)/rad)
	sdip=sin(fdip(jf)/rad)
	c2dip=cdip*cdip-sdip*sdip
	srak=sin(frake(jf))
	crak=cos(frake(jf))
	dlen=(fbigl(jf)/ethr)/real(nmesh2-1)
	dlon=olon(i)-flon(jf)
c	Find angular distance and azimuth to observation point from
c	the initial point on the fault segment.
c	Note: delta is the angular distance from the earthquake.
c	Note: phi is the azimuth meazsured positive counterclockwise
c	from south.
	cdelt=cos(flat(jf))*cos(olat(i))+sin(flat(jf))*
     &	sin(olat(i))*cos(dlon)
	if (cdelt.le.0.9999) delta=acos(cdelt)
	if(cdelt.gt.0.9999) delta=sqrt((flat(jf)-olat(i))**2+
     &	(dlon*sin(flat(jf)))**2)
	spsi=sin(dlon)*sin(olat(i))/sin(delta)
	cpsi=(cos(olat(i))-cos(flat(jf))*cdelt)/(sin(flat(jf))*sin(delta))
	phi=pi-atan2(spsi,cpsi) 
	dwrite=rad*delta
	pwrite=rad*phi
	  write(6,*)'olat(',i,')=',olat(i),'olon=',olon(i)
	  write(6,*)'flat(',jf,')=',flat(jf),'flon=',flon(jf)
	  write(6,*)'delta,phi=',dwrite,'deg.',pwrite 
	cphi=cos(phi)
	sphi=sin(phi)
cNEW
	spsi=sin(dlon)*sin(flat(jf))/sin(delta)
	cpsi=(cos(flat(jf))-cos(olat(i))*cdelt)/(sin(olat(i))*sin(delta))
	peps=atan2(spsi,cpsi)
	cphio=cos(peps)
	sphio=sin(peps)
c--
	if(jf.gt.1) go to 59
	cphi1=cphi
	sphi1=sphi
	delt1=delta 
59	delta=delta*(6371./ethr)
	  write(6,*)'ethr=',ethr
	  write(6,*)'delta (after mult)=',delta
c	Angular distances are a factor of (6371./ethr) larger on the
c	earth with radius ethr km.
	  write(6,*)'after 59, delta=',delta
c	Integrate displacements over fault elements.
c ***	Determine moment tensor elements for fault element (ilen).
c	Units of moment tensor are 10**20 N-m.
cDONT NEED	mu=real(ymu)
cDONT NEED	lam=real(ylam)
	  write(6,*)'mu,lam=',mu,lam
	amesh=real(nmesh2-1)
c	Use input moment tensor elements if scalar moment < 0
	if(fmom(jf).lt.0.0) then
	mrr=mrr0(jf)/amesh
	mtt=mtt0(jf)/amesh
	mpp=mpp0(jf)/amesh
	mrt=mrt0(jf)/amesh
	mrp=mrp0(jf)/amesh
	mtp=mtp0(jf)/amesh
	go to 51
	endif
c	Note that shear modulus [mu] is already input in units of 10**10 Pa.
c	Moment tensor from Ben Menahem and Singh, eqn. 4.115b for
c	shear dislocation, and derived from eqn. (4.101), (4.110), and
c	(4.113) for a tensile dislocation.
c	Next line is shear moment.
	shrm1=fmom(jf)/amesh
c--
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
51	continue
	  write(6,*)'mrr,mtt,mpp,mrt,mrp,mtp='
	  write(6,*) mrr,mtt,mpp,mrt,mrp,mtp
c * *
c	Determine multiplying factors for the cases m=0,1,2 
	fm(1)=dble((1./(2.*mu*mu+3.*lam*mu))*((lam+mu)*mrr-(lam/2.)*(mtt+mpp)))
	fm(2)=dble((1./(2.*mu*mu+3.*lam*mu))*(-lam*mrr+(lam/2.+mu)*(mtt+mpp)))
	fm(3)=dble((1./(2.*mu))*mrt)
	fm(4)=dble(-(1./(2.*mu))*mrp)
	fm(5)=dble((1./(2.*mu))*(mtt-mpp))
	fm(6)=dble(-(1./mu)*mtp)
	len=-dlen/2.
	ilen=0
83	ilen=ilen+1
	if(ilen.eq.nmesh2) go to 65
c	  write(6,*)'ilen=',ilen
	len=len+dlen
	bigx=delta*sphi+len*sstr
	bigy=-delta*cphi+len*cstr
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
	cotd=cl/sl
	slf=real(sl)
	deltpf=real(deltp)
c	If fault segment is within 0.2 km of observation point, ignore
c	the fault segment
cOLD	if(deltpf.lt.0.00003) go to 83
	call lgndr0(deltpf,xl0,dxl0)
	call lgndr1(deltpf,xl1,dxl1)
	call lgndr2(deltpf,xl2,dxl2) 
c	Store Legendre functions for later use.
	do 64 lk=lmin,lmax+1
	faclk=dsqrt(dble(lk*(lk+1)))
	x1=xl0(lk)
	x2=xl1(lk)
	x3=xl2(lk)
	dx1=dxl0(lk)
	dx2=dxl1(lk)
	dx3=dxl2(lk)
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
64	continue 
	do 71 lr=lmin,lmax+1
	l=lr
c	Do l=0 case last.
	if(lr.eq.(lmax+1)) l=0
	dispx1=0.d0
	dispy1=0.d0
	dispz1=0.d0
	exx1=0.d0
	eyy1=0.d0
	exy1=0.d0
	exz1=0.d0
	eyz1=0.d0
	ezz1=0.d0
c	With facf apply tapering over last quarter of frequency range.
	fl=real(l)
	flmax=real(lmax)
	iw0m=lmax/5
	iw0x=lmax-iw0m+1
	fiw0m=real(iw0m)
	facf=1.0
	if(l.ge.iw0x) facf=0.5-0.5*cos((flmax-fl)/(fiw0m-1.) * real(pi))
c		write(6,*)'A: l=',l,'facf=',facf,'num=',flmax-fl,'den=',fiw0m-1.
c
	dfac=dble(facf*(1.e+6)/ethr**2)
	dfac1=dfac 
c	Now account for the fact that the observation radius is generally
c	different from the surface radius.  This affects displ calculations
c	(for example, y2 and y4-arrays have functions
c	y1(r)/r and y3(r)/r, where r is dimensionless radius,
c	which equals 1 only at the surface), but not strain calculations.
	dfac=dfac*real(odep)
c	SPHEROIDAL MODES
	do 75 idec=1,5
	k=(nleno/3)*idec-(nleno/3)+iom
	wkv1(idec)=y1(k,lr)
	wkv2(idec)=y2(k,lr)
	wkv3(idec)=y3(k,lr)
	wkv4(idec)=y4(k,lr)
75	continue
c		write(6,*)'wkv4=',wkv4
c	Note extra factor of 2 in z2 - z5 components in order to account
c	for m=-1 and m=-2 contributions.
c * *	Moment tensor excitation
	z1=fm(1)*(wkv2(1))+fm(2)*(wkv2(2))
	z2=2.*fm(3)*(wkv2(3))
	z3=-2.*fm(4)*(wkv2(3))
	z4=2.*fm(5)*(wkv2(4))
	z5=-2.*fm(6)*(wkv2(4))
c		write(6,*)'fm(5),wkv2(4),z4=',fm(5),wkv2(4),z4
c		write(6,*)'fm(6),wkv2(4),z5=',fm(6),wkv2(4),z5
	dz1=fm(1)*(wkv1(1))+fm(2)*(wkv1(2))
	dz2=2.*fm(3)*(wkv1(3))
c		write(6,*)'z2,dz2=',z2,dz2
	dz3=-2.*fm(4)*(wkv1(3))
	dz4=2.*fm(5)*(wkv1(4))
	dz5=-2.*fm(6)*(wkv1(4))
	dz1h=fm(1)*(wkv3(1))+fm(2)*(wkv3(2))
	dz2h=2.*fm(3)*(wkv3(3))
	dz3h=-2.*fm(4)*(wkv3(3))
	dz4h=2.*fm(5)*(wkv3(4))
	dz5h=-2.*fm(6)*(wkv3(4))
	z1h=fm(1)*(wkv4(1))+fm(2)*(wkv4(2))
	z2h=2.*fm(3)*(wkv4(3))
	z3h=-2.*fm(4)*(wkv4(3))
	z4h=2.*fm(5)*(wkv4(4))
	z5h=-2.*fm(6)*(wkv4(4))
c * *
c	  write(6,*)'fm(6),wkv4(4),z5h=',fm(6),wkv4(4),z5h
	if(l.eq.0) then
	x1=dble(sqrt(0.25/pi))
	x2=0.d0
	x3=0.d0
	dx1=0.d0
	dx2=0.d0
	dx3=0.d0
	ddx1=0.d0
	ddx2=0.d0
	ddx3=0.d0
	else
	x1=x1s(lr)
	x2=x2s(lr)
	x3=x3s(lr)
	dx1=dx1s(lr)
	dx2=dx2s(lr)
	dx3=dx3s(lr)
	ddx1=ddx1s(lr)
	ddx2=ddx2s(lr)
	ddx3=ddx3s(lr)
		endif
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
c 
	disp1=(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac
	disp1h=(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac
	disp1t=(-z2h*x2*sp
     &	+z3h*x2*cp -2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac/slf 
	ett=(z1h*ddx1 +z2h*ddx2*cp
     &	+z3h*ddx2*sp +z4h*ddx3*cp2 +z5h*ddx3*sp2)*dfac1/bigr 
	ett=ett+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac1/bigr 
	epp=(-z2h*x2*cp -z3h*x2*sp 
     &	-4.*z4h*x3*cp2 -4.*z5h*x3*sp2)*(dfac1/slf)/(slf*bigr) 
	epp=epp+(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac1*cotd/bigr 
	epp=epp+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac1/bigr 
	etr=(dz1h*dx1 +dz2h*dx2*cp
     &	+dz3h*dx2*sp +dz4h*dx3*cp2 +dz5h*dx3*sp2)*dfac1/(2.*bigr) 
	etr=etr-(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac1/(2.*bigr) 
	etr=etr+(z1*dx1 +z2*dx2*cp
     &	+z3*dx2*sp +z4*dx3*cp2 +z5*dx3*sp2)*dfac1/(2.*bigr) 
	etp=(-z2h*dx2*sp
     &	+z3h*dx2*cp -2.*z4h*dx3*sp2 +2.*z5h*dx3*cp2)*dfac1/(2.*bigr*slf) 
	etp=etp-(-z2h*x2*sp +z3h*x2*cp 
     &	-2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac1*cotd/(bigr*slf) 
	etp=etp+(-z2h*dx2*sp +z3h*dx2*cp 
     &	-2.*z4h*dx3*sp2 +2.*z5h*dx3*cp2)*dfac1/(2.*slf*bigr)
	epr=(-z2*x2*sp
     &	+z3*x2*cp -2.*z4*x3*sp2 +2.*z5*x3*cp2)*dfac1/(2.*slf*bigr) 
	epr=epr+(-dz2h*x2*sp +dz3h*x2*cp 
     &	-2.*dz4h*x3*sp2 +2.*dz5h*x3*cp2)*dfac1/(2.*bigr*slf) 
	epr=epr-(-z2h*x2*sp +z3h*x2*cp 
     &	-2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac1/(2.*bigr*slf) 
	err=(dz1*x1 +dz2*x2*cp
     &	+dz3*x2*sp +dz4*x3*cp2 +dz5*x3*sp2)*dfac1/bigr 
c		write(6,*)'etp,err=',etp,err
c		write(6,*)'x1,x2,x3=',x1,x2,x3
c	Rotate strain tensor into x-y-z Cartesian coordinates.
c	Also rotate displacements into x-y-z coordinates.
	dispx1=dispx1+gamma1*disp1h-gamma2*disp1t
		if(iom.eq.1) write(6,*)'SPH: l,dispx1=',l,dispx1
	dispy1=dispy1+gamma2*disp1h+gamma1*disp1t 
c**
	  if(iom.eq.2) write(6,*) 'SPH: l,disp1h,disp1t,dispy1=',l,disp1h,disp1t,dispy1
	dispz1=dispz1+disp1
c		if(l.eq.500) then
c		write(6,*)'l,disp1=',l,disp1
c		write(6,*)'disp1t=',disp1t,'gamma1=',gamma1
c		write(6,*)'disp1h=',disp1h,'gamma2=',gamma2
c		endif
	exx1=exx1+gamma1**2*ett+gamma2**2*epp-2.*gamma1*gamma2*etp
	exy1=exy1+gamma1*gamma2*(ett-epp)+(gamma1**2-gamma2**2)*etp
	eyy1=eyy1+gamma2**2*ett+gamma1**2*epp+2.*gamma1*gamma2*etp
	exz1=exz1+etr*gamma1-epr*gamma2
	eyz1=eyz1+etr*gamma2+epr*gamma1
	ezz1=ezz1+err
c	Done with Spheroidal motion component of degree l.
c	  write(6,*)'l=',l,'SPH dispz1=',dispz1
	if(l.eq.0) go to 71
c	TOROIDAL MODES
	do 175 idec=6,7
	k=(nleno/3)*idec-(nleno/3)+iom
	wkv1(idec)=y1(k,lr)
	wkv2(idec)=y2(k,lr)
175	continue
	z1=-2.*fm(4)*(wkv2(6))
	z2=-2.*fm(3)*(wkv2(6))
	z3=-2.*fm(6)*(wkv2(7))
	z4=-2.*fm(5)*(wkv2(7))
	dz1=-2.*fm(4)*(wkv1(6)+wkv2(6))
	dz2=-2.*fm(3)*(wkv1(6)+wkv2(6))
	dz3=-2.*fm(6)*(wkv1(7)+wkv2(7))
	dz4=-2.*fm(5)*(wkv1(7)+wkv2(7))
	disp1t=(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac
c		write(6,*)'TOR: z1,z2,z3,,z4=',z1,z2,z3,z4
	disp1h=-(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac/slf 
	ett=-(-z1*dx2*sp +z2*dx2*cp
     &	-2.*z3*dx3*sp2 +2.*z4*dx3*cp2)*dfac1/(bigr*slf) 
	ett=ett+(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1*cotd/(bigr*slf) 
	epp=(-z1*dx2*sp +z2*dx2*cp 
     &	-2.*z3*dx3*sp2 +2.*z4*dx3*cp2)*dfac1/(bigr*slf)
	epp=epp-(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1*cotd/(bigr*slf) 
	etp=(z1*ddx2*cp +z2*ddx2*sp
     &	+z3*ddx3*cp2 +z4*ddx3*sp2)*dfac1/(2.*bigr)
	etp=etp-(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac1*cotd/(2.*bigr)
	etp=etp-(-z1*x2*cp -z2*x2*sp
     &	-4.*z3*x3*cp2 -4.*z4*x3*sp2)*(dfac1/slf)/(2.*bigr*slf)
	etr=-(-dz1*x2*sp +dz2*x2*cp
     &	-2.*dz3*x3*sp2 +2.*dz4*x3*cp2)*dfac1/(2.*bigr*slf) 
	etr=etr+(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1/(2.*bigr*slf) 
	epr=(dz1*dx2*cp +dz2*dx2*sp
     &	+dz3*dx3*cp2 +dz4*dx3*sp2)*dfac1/(2.*bigr)
	epr=epr-(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac1/(2.*bigr)
c	Rotate strain tensor into x-y-z Cartesian coordinates.
c	Also rotate displacements into x-y coordinates.
	dispx1=dispx1+gamma1*disp1h-gamma2*disp1t
	dispy1=dispy1+gamma2*disp1h+gamma1*disp1t
	  if(iom.eq.2) write(6,*) 'TOR: l,disp1h,disp1t,dispy1=',l,disp1h,disp1t,dispy1
	exx1=exx1+gamma1**2*ett+gamma2**2*epp-2.*gamma1*gamma2*etp
	exy1=exy1+gamma1*gamma2*(ett-epp)+(gamma1**2-gamma2**2)*etp
	eyy1=eyy1+gamma2**2*ett+gamma1**2*epp+2.*gamma1*gamma2*etp
	exz1=exz1+etr*gamma1-epr*gamma2
	eyz1=eyz1+etr*gamma2+epr*gamma1
c	Done with Toroidal motion component of degree l.
		if(iom.eq.1) write(6,*)'TOR: l,dispx1=',l,dispx1
	dispx(iom)=dispx(iom)+dispx1 
	dispy(iom)=dispy(iom)+dispy1 
		if(iom.eq.2) write(6,*)'LINEA: dispy(2)=',dispy(iom)
	dispz(iom)=dispz(iom)+dispz1 
c		write(6,*)'l,ilen=',l,ilen,'latest dispz(',iom,')=',dispz(iom)
	exx(iom)=exx(iom)+exx1
	eyy(iom)=eyy(iom)+eyy1
	exy(iom)=exy(iom)+exy1
	exz(iom)=exz(iom)+exz1
	eyz(iom)=eyz(iom)+eyz1
	ezz(iom)=ezz(iom)+ezz1
c	  write(6,*) 'l=',l,'SPH+TOR dispz1=',dispz1
71	continue
c	  write(6,*)'latest dispx,dispy,dispz=',dispx(iom),dispy(iom),dispz(iom)
	go to 83
84	continue
c	Use the needed frequency samples for the inverse FT.
c	jmax1=int((2.*pi*2./corper)/real(ommin))
	jmax1=nleno/2
c--
	do iom=1,nleno/3
c	Multiply everything by 1/(ui*ome) to account for step-function
c	character of source.
	ome=real(iom-1)*ommin - ui*obeta
	oome=-ui/ome
	dispx(iom)=dispx(iom)*oome
	dispy(iom)=dispy(iom)*oome
	dispz(iom)=dispz(iom)*oome
	exx(iom)=exx(iom)*oome
	eyy(iom)=eyy(iom)*oome
	ezz(iom)=ezz(iom)*oome
	exz(iom)=exz(iom)*oome
	eyz(iom)=eyz(iom)*oome
	exy(iom)=exy(iom)*oome
	enddo
c
	  do iom=1,jmax1
	  xo(iom)=dispx(iom)
c**	Use following line for radial displacement **
	if(comp.eq.'rtz') xo(iom)=dispx(iom)*gamma1+dispy(iom)*gamma2
c**
		write(6,*)'DISPX: xo(',iom,')=',xo(iom)
	  enddo
		write(6,*)'entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,4*jmax1
	t=0.25*dt*real(it)
	write(2,*) t,xt(it)
	enddo
c
	  do iom=1,jmax1
	  xo(iom)=dispy(iom)
c**	Use following line for transverse displacement **
	if(comp.eq.'rtz') xo(iom)=-dispx(iom)*gamma2+dispy(iom)*gamma1
c**
		write(6,*)'DISPY: xo(',iom,')=',xo(iom)
	  enddo
		write(6,*)'entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,4*jmax1
	t=0.25*dt*real(it)
	write(4,*) t,xt(it)
	enddo
c
	  do iom=1,jmax1
	  xo(iom)=dispz(iom)
		write(6,*)'DISPZ: xo(',iom,')=',xo(iom)
	  enddo
		write(6,*)'entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,4*jmax1
	t=0.25*dt*real(it)
	write(8,*) t,xt(it)
	enddo
c
	  do iom=1,jmax1
	  xo(iom)=exx(iom)
		write(6,*)'EXX: xo(',iom,')=',xo(iom)
	  enddo
		write(6,*)'entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,4*jmax1
	uexx(it,i)=xt(it)
	enddo
c
	  do iom=1,jmax1
	  xo(iom)=exy(iom)
		write(6,*)'EXY: xo(',iom,')=',xo(iom)
	  enddo
		write(6,*)'entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,4*jmax1
	uexy(it,i)=xt(it)
	enddo
c
	  do iom=1,jmax1
	  xo(iom)=exz(iom)
		write(6,*)'EXZ: xo(',iom,')=',xo(iom)
	  enddo
		write(6,*)'entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,4*jmax1
	uexz(it,i)=xt(it)
	enddo
c
	  do iom=1,jmax1
	  xo(iom)=eyy(iom)
		write(6,*)'EYY: xo(',iom,')=',xo(iom)
	  enddo
		write(6,*)'entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,4*jmax1
	ueyy(it,i)=xt(it)
	enddo
c
	  do iom=1,jmax1
	  xo(iom)=eyz(iom)
		write(6,*)'EYZ: xo(',iom,')=',xo(iom)
	  enddo
		write(6,*)'entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,4*jmax1
	ueyz(it,i)=xt(it)
	enddo
c
	  do iom=1,jmax1
	  xo(iom)=ezz(iom)
		write(6,*)'EZZ: xo(',iom,')=',xo(iom)
	  enddo
		write(6,*)'entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,4*jmax1
	uezz(it,i)=xt(it)
	enddo
c
60	continue
	close(2)
	close(4)
	close(8)
c	Write out strain seismograms in units of microstrain/sec
		write(6,*)'will write out strains'
	open(2,file='wave1.out_strains')
		write(6,*)'line A'
		write(6,*)'ipts,jmax1=',ipts,jmax1
	do i=1,ipts
		write(6,*)'doing i=',i,'out of',ipts
	do it=1,4*jmax1
	t=0.25*dt*real(it)
	write(2,70) t,uexx(it,i),ueyy(it,i),uexy(it,i),uexz(it,i),ueyz(it,i),uezz(it,i)
	enddo
	enddo
	close(2)
70	format(7e14.6e2)
c
	write(6,*)'end of wave1'
	end  
	 
