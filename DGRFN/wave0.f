	program WAVE0
c***
c	Find Green's functions for seismic deformation field
c	due to buried seismic source at the origin of a spherical coordinate 
c	system.  Deformation field is calculated as a sum of spherical
c	harmonics from l=0 to [lmax].  Only m=0, 1, 2 azimuthal order numbers 
c	contribute.  Green's functions are calculated for
c	one source depth and one observation depth
c	Spherical harmonic degree numbers up to 25000 may be used.
c***	INPUT
c	maximum spherical harmonic degree [lmax]
c	source depth [sdep] (km)
c	observation depth [odepr] (km)
c	time increment [dt] (sec), used to calculate Nyquist freq. [ommin]
c	# samples in frequency domain [nleno]
c	'earth.model' has Earth's radius, # layers, elastic parameters and density
c	and bulk and shear attenuation parameters in each layer.
c	Attenuative dispersion assumes a reference frequency of 1 Hz.
c	Input file 'earth.model' has format
c	N  Nd 	(Number of layers; Number of different layers)
c	rb(1)	     |  rt(1)   |dens(1)|  kappa(1)  | mu(1)  | qbet(1) | qkap(1)
c	...		...	 ...	   ...	       ...      ...     | ...
c	rb(N)	     |  rt(N)   |dens(N)|  kappa(N)  | mu(N)  | qbet(N) | qkap(N)
c	bottom depth|top depth|density|bulk modulus|rigidity| Q_beta  | Q_kappa
c	(km)          (km)      (g/cm**3) (10**11 dyne/cm**2)  
c***	OUTPUT
c	Green's functions are in wave0.out, which will be used by WAVE1.
c***
c							       
c		                                               
	character*1 csys
	character*80 aread
	real*4 kappa,mu
	double complex ub,ua,wb,wa,sb,sa,vb,rx,va,ry
	double complex us1,vs1,us2,vs2,up1,vp1
	double complex up2,vp2,us3,vs3,us4,vs4,up3,vp3,up4,vp4
	double complex dur3,dur4,dvr3,dvr4,ur3,ur4,vr3,vr4
	double complex ur1,ur2,vr1,vr2,dur1,dur2,dvr1,dvr2
	double complex wr1,dwr1,wr2,dwr2
	double complex usupp,vsupp,upupp,vpupp,urupp,vrupp,durupp,dvrupp
	double complex uslow,vslow,uplow,vplow,urlow,vrlow,durlow,dvrlow
	double complex bmf1(201,51),bmf2(201,51),bmf3(201,51),bmf4(201,51)
	double complex bmf5(201,51),bm(5)
	real*8 valt
	real*8 r,dr,r2,fl21,flp1,odep
	real*8 denss
	double complex biga,lams,mus,kaps
	double complex c(4,4),b1(4),b2(4),b3(4),b4(4),bupp(4),blow(4)
	double complex bm1(5),bm2(5),dydrm(5)
c	double complex t1try(2),bm1try(5)
	double complex tmp1,tmp2,tmp3,tmp4,tmp5
	double complex b1s(4),b2s(4),b3s(4),b4s(4)
	double complex ct(2),t1(2),t2(2),t1s(2),t2s(2)
	double complex bfacl,bfacm,y1,y2,y3,y4
	real*8 aden1,aden2
	real*8 ylim(3)
	real*8 a00(4,4)
	double complex bc0(4,4),c0(4),a0(4)
	double complex bc0s(4,4),c0s(4)
	double complex bc1(2,2),c1(2),a1(2)
	double complex row1,row2,row3,row4
	dimension eta(200)
c--
	double complex dydr(4),dtdr(2)
	real*8 htot,eps,rsav,hdid,hnext,yscalm(5),yscal(4),yscalt(2),r0,rintp,dl
	common/radlev/rintp,bigr
c--
	double complex w0,ommin,obeta,ui
	real*8 w0r,w0i,wmag,rterm,sterm
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),w0,r0
	common/qval/qbet(200),qkap(200)
	common/aterm/rterm,sterm
	common/matl/l
	pi=3.14159265
	twopi=2.*pi
	ui=dcmplx(0.,1.d0)
	lmin=1
	write(6,*)'lmax?'
	read(5,*) lmax
c *
	write(6,*)'source depth (km)?'
	read(5,*) sdep
	write(6,*)'observation depth (km)?'
	read(5,*) odepr
	open(2,file='wave0.out',form='unformatted')
c21	format(i2)
c 
c ** *	
c	Read in earth model
	open(4,file='earth.model')
	rewind(4)
	read(4,5) n,bigr
	write(6,5) n,bigr
5	format(i4,f10.3)
c	bigr=earth radius in km
	dep=1.0-sdep/bigr
	  write(6,*)'dep=',dep
	write(2) bigr,sdep
	write(2) lmin,lmax
	odep=dble(1.-odepr/bigr)
	write(2) odep
	do 10 j=1,n
	i=2*j-1
c--
	read(4,15) rb(i),rt(i),dens(i),kappa(i),mu(i),qbet(i),qkap(i)
	write(6,15) rb(i),rt(i),dens(i),kappa(i),mu(i),qbet(i),qkap(i)
c	Replace qbet and qkap with their inverses.
	qbet(i)=1./qbet(i)
	qkap(i)=1./qkap(i)
c--
	dens(i+1)=dens(i)
	mu(i+1)=mu(i)
	eta(i+1)=eta(i)
	kappa(i+1)=kappa(i)
	qbet(i+1)=qbet(i)
	qkap(i+1)=qkap(i)
	rb(i+1)=(rb(i)+rt(i))/2.
	rt(i+1)=rt(i)
	rt(i)=rb(i+1)
cOLD	rb(i+1)=rb(i+1)/bigr + (1.-6371./bigr)
cOLD	rt(i+1)=rt(i+1)/bigr + (1.-6371./bigr)
cOLD	rb(i)=rb(i)/bigr + (1.-6371./bigr)
cOLD	rt(i)=rt(i)/bigr + (1.-6371./bigr)
	rb(i+1)=1. - rb(i+1)/bigr 
	rt(i+1)=1. - rt(i+1)/bigr 
	rb(i)=1. - rb(i)/bigr 
	rt(i)=1. - rt(i)/bigr
10	continue
13	format(a80)
15	format(7f9.3)
	n=2*n 
	close(4)
c	Put obs pt at boundary between layer [iodep-1] and [iodep].
c	This involves putting in an extra layer within a pre-existing
c	layer.
	if(odep.eq.1.d0) then
	iodep=n+1
	go to 120
	endif
	if(odep.lt.rb(n)) go to 112
	rt(n)=odep
	rb(n+1)=odep
	rt(n+1)=1.0
	dens(n+1)=dens(n)
	qbet(n+1)=qbet(n)
	qkap(n+1)=qkap(n)
	eta(n+1)=eta(n)
	eta(n+1)=eta(n)
	kappa(n+1)=kappa(n)
	mu(n+1)=mu(n)
	iodep=n+1
	go to 118
112	iodep=0
	do 116 j=2,n
	i=n-j+2
	if(rb(i).lt.odep.or.rb(i-1).gt.odep) go to 114
	write(6,*)'put in obs pt at depth',bigr*(1.-odep),'km'
	write(6,*)'rb(',i,')=',rb(i),'rt(',i,')=',rt(i)
	iodep=i
	rt(i-1)=odep
	rb(i+1)=rb(i)
	rt(i+1)=rt(i)
	rt(i)=rb(i)
	rb(i)=odep
	  dens(i+1)=dens(i)
	qbet(i+1)=qbet(i)
	qkap(i+1)=qkap(i)
	eta(i+1)=eta(i)
	kappa(i+1)=kappa(i)
	mu(i+1)=mu(i)
	  dens(i)=dens(i-1)
	qbet(i)=qbet(i-1)
	qkap(i)=qkap(i-1)
	eta(i)=eta(i-1)
	kappa(i)=kappa(i-1)
	mu(i)=mu(i-1)
114	if(iodep.gt.0) go to 116
	rt(i+1)=rt(i)
	rb(i+1)=rb(i)
	  dens(i+1)=dens(i)
	qbet(i+1)=qbet(i)
	qkap(i+1)=qkap(i)
	eta(i+1)=eta(i)
	kappa(i+1)=kappa(i)
	mu(i+1)=mu(i)
116	continue
118	continue
	n=n+1
120	continue
		write(6,*)'After inserting obs depth: number of layers=',n
c	Put source at boundary between layer [isdep-1] and [isdep].
c	This involves putting in an extra layer within a pre-existing
c	layer.
	if(dep.lt.rb(n)) go to 12
	rt(n)=dep
	rb(n+1)=dep
	rt(n+1)=1.0
	dens(n+1)=dens(n)
	qbet(n+1)=qbet(n)
	qkap(n+1)=qkap(n)
	eta(n+1)=eta(n)
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
	  dens(i+1)=dens(i)
	qbet(i+1)=qbet(i)
	qkap(i+1)=qkap(i)
	eta(i+1)=eta(i)
	kappa(i+1)=kappa(i)
	mu(i+1)=mu(i)
	  dens(i)=dens(i-1)
	qbet(i)=qbet(i-1)
	qkap(i)=qkap(i-1)
	eta(i)=eta(i-1)
	kappa(i)=kappa(i-1)
	mu(i)=mu(i-1)
14	if(isdep.gt.0) go to 16
	rt(i+1)=rt(i)
	rb(i+1)=rb(i)
	  dens(i+1)=dens(i)
	qbet(i+1)=qbet(i)
	qkap(i+1)=qkap(i)
	eta(i+1)=eta(i)
	kappa(i+1)=kappa(i)
	mu(i+1)=mu(i)
16	continue
18	continue
	n=n+1
		write(6,*)'After inserting source: number of layers=',n
	if(dep.lt.odep) then
	write(6,*)'Must revise layer # of obs depth'
	iodep=iodep+1
	endif
	write(6,*)'final layer # of obs pt=',iodep
	write(6,*)'final layer # of source=',isdep
c	check the model for correct interpolation.
	do j=1,n
	write(6,15) bigr*(1.-rb(j)),bigr*(1.-rt(j)),dens(j),kappa(j),mu(j),
     &	1./qbet(j),1./qkap(j)
	enddo
	write(2) kappa(isdep),mu(isdep)
	nrad=5
c ** **
	write(6,*)'time increment?'
	read(5,*) dt
	write(6,*)'# time samples?'
	read(5,*) nleno
	nleno=nleno*2
	write(2) dt,nleno
	ommin=dble(2.*pi/dt)/dble(nleno)
	obeta=2.0d0*ommin
c	Anticipating corner periods in WAVE1 that are >= 7 x (1/dt),
c	do only nleno/3 frequencies (thus anticipated value of jmax1 in
c	INVFT, called by WAVE1, will not exceed nleno/3).
	do 100 iom=1,nleno/3
	w0=dble(iom-1)*ommin - ui*obeta
c	Assume a reference frequency of twopi rad/sec (i.e. 1 Hz).
	w0r=dreal(w0)/dble(twopi)
	w0i=dimag(w0)/dble(twopi)
	wmag=dsqrt(w0r*w0r+w0i*w0i)
	rterm=dlog(wmag)
	sterm=datan2(w0i,w0r)
c ** *
	do 50 lr=lmin,lmax+1
	l=lr
c	Do l=0 case last.
	if(lr.eq.(lmax+1)) l=0
	write(6,*)'l=',l,'iom=',iom
c	depmax=dep-2.5*(2.*pi/(real(l)+0.5))
c	depmin=dep+2.5*(2.*pi/(real(l)+0.5))
	depmax=rb(1)
	depmin=1.01
		write(6,*)'depmax,depmin,odep=',depmax,depmin,odep
	  awr=(1.-real(depmax))*bigr
	  write(6,*)'max integration depth=',awr
c	depmax=dimensionless radius 2 wavelengths into the earth.
	dl=dble(l)
	fl21=dble(l*(l+1))
	flp1=dble(l+1)
c	SPHEROIDAL MODES
c	Start matrix propagation upward for solution bm1.  Stop
c	at layer [isdep-1].
c	Give values at CMB those appropriate for homogeneous solid
c	in 0<r<rb(1).
	is=0
	do 30 i=1,isdep-1
c	if(i.lt.2) go to 30
c	Start integration at most 2.5 wavelengths into the earth.
	if(rt(i).lt.depmax.or.rb(i).gt.depmin) go to 30
cOLD	if(dble(depmax).ge.odep.or.dble(depmin).le.odep) go to 30
cOLD	mus=dble(mu(i))
cOLD	lams=dble(kappa(i)-2.*mu(i)/3.)
cOLD	biga=lams+2.d0*mus
c		r0=0.5d0*dble(rb(i)+rt(i))
c		call interp(mus,lams,denss)
c		biga=lams+2.d0*mus
	if(is.ne.0) go to 28
	is=1
cOLD	Give boundary conditions appropriate for an incompressible
cOLD	fluid.
cOLD	b1(1)=dble(l)
cOLD	b1(2)=0.d0
cOLD	b1(3)=1.d0
cOLD	b1(4)=0.d0
cOLD	b2(1)=0.d0
cOLD	b2(2)=0.d0
cOLD	b2(3)=1.d0
cOLD	b2(4)=0.d0
cOLD	bm1(1)=b1(1)*b2(2)-b1(2)*b2(1)
cOLD	bm1(2)=b1(1)*b2(3)-b1(3)*b2(1)
cOLD	bm1(3)=b1(1)*b2(4)-b1(4)*b2(1)
cOLD	bm1(4)=b1(2)*b2(3)-b1(3)*b2(2)
cOLD	bm1(5)=b1(2)*b2(4)-b1(4)*b2(2)
	r0=dble(rb(i))
	csys='s'
		write(6,*)'entering lbound - SPH'
	call lbound(bm1,5,csys)
		write(6,*)'out of lbound - SPH: bm1=',bm1
c	
28	continue

29	dr=dble(rt(i)-rb(i))/dble(nrad)
c	Integrate equations of motion to propagate solutions up to
c	the surface.
	r=dble(rb(i))
c	Do bm1 using bsstep.
c		write(6,*)'nrad=',nrad
	iwk=0
c	Have index ,1 correspond to bottom of layer.
	bmf1(i,1)=bm1(1)
	bmf2(i,1)=bm1(2)
	bmf3(i,1)=bm1(3)
	bmf4(i,1)=bm1(4)
	bmf5(i,1)=bm1(5)
	do 31 j=1,nrad
	rsav=r
	rintp=r+dr/2.d0
	call derivs5(r,bm1,dydrm)
	htot=dr
	eps=3.d-6
	yscalm(1)=zabs(bm1(1))+dl*zabs(bm1(3))+zabs(bm1(2))*dl**2+zabs(bm1(4))*dl
     &		 +zabs(bm1(5))
	if(l.gt.0) then
	yscalm(2)=yscalm(1)/dl**2
	yscalm(3)=yscalm(1)/dl
	else
	yscalm(1)=zabs(bm1(2))
	yscalm(2)=zabs(bm1(2))
	yscalm(3)=zabs(bm1(2))
	endif
	yscalm(4)=yscalm(3)
	yscalm(5)=yscalm(1)
c		write(6,*)'entering bsstep A'
	call bsstep(bm1,dydrm,5,r,htot,eps,yscalm,hdid,hnext)
c		write(6,*)'out of bsstep A'
	bmf1(i,j+1)=bm1(1)
	bmf2(i,j+1)=bm1(2)
	bmf3(i,j+1)=bm1(3)
	bmf4(i,j+1)=bm1(4)
	bmf5(i,j+1)=bm1(5)
	if((hdid-htot)/htot.ne.0.d0) iwk=1
c		write(6,*)'hdid,htot=',hdid,htot,'iwk=',iwk
c
32	r=rsav+dr    
31	continue 
c		write(6,*)'TOP layer',i,'bm1=',bm1
c		r0=r
c		call lbound(bm1try,5,csys)
c		write(6,*)'lbound: scaled bm1='
c		do jtry=1,5
c		write(6,*) bm1try(jtry)*bm1(5)/bm1try(5)
c		enddo
46	bmag=zabs(bm1(1))+real(l)*zabs(bm1(3))
	if(bmag.lt.1.e+10) go to 30
	do j1=1,5
	bm1(j1)=bm1(j1)*1.d-20
	enddo
	go to 46
30	continue
c	  write(6,*)'finished bm1'
c *
c	Start matrix propagation downward for solution bm2.  Stop
c	at layer [isdep].
	b3(1)=dble(l+1)
	b3(2)=0.d0
	b3(3)=1.d0
	b3(4)=0.d0
	b4(1)=0.d0
	b4(2)=0.d0
	b4(3)=1.d0
	b4(4)=0.d0
	bm2(1)=b3(1)*b4(2)-b3(2)*b4(1)
	bm2(2)=b3(1)*b4(3)-b3(3)*b4(1)
	bm2(3)=b3(1)*b4(4)-b3(4)*b4(1)
	bm2(4)=b3(2)*b4(3)-b3(3)*b4(2)
	bm2(5)=b3(2)*b4(4)-b3(4)*b4(2)
		write(6,*)'Initial bm2=',bm2
		write(6,*)'bm2(2)=',bm2(2)
	do 130 i=n,isdep,-1
	if(rb(i).gt.depmin) go to 130
cOLD	mus=dble(mu(i))
cOLD	lams=dble(kappa(i)-2.*mu(i)/3.)
cOLD	biga=lams+2.d0*mus
c		mus=dble(mu(i)) * dcmplx(1.0,qbet(i))
c		kaps=dble(kappa(i)) * dcmplx(1.0,qkap(i))
c		lams=kaps-(2.d0/3.d0)*mus
c		biga=lams+2.d0*mus
	dr=-dble(rt(i)-rb(i))/dble(nrad)
c	Integrate equations of motion to propagate solutions downward from
c	the surface.
cOLD	r=dble(rb(i))
	r=dble(rt(i))
c	Do bm2 using bsstep.
c		write(6,*)'nrad=',nrad
	iwk=0
c	Have index ,1 correspond to top of layer.
	bmf1(i,1)=bm2(1)
	bmf2(i,1)=bm2(2)
	bmf3(i,1)=bm2(3)
	bmf4(i,1)=bm2(4)
	bmf5(i,1)=bm2(5)
	do 131 j=1,nrad
	rsav=r
	rintp=r+dr/2.d0
	call derivs5(r,bm2,dydrm)
	htot=dr
	eps=3.d-6
	yscalm(1)=zabs(bm2(1))+dl*zabs(bm2(3))+zabs(bm2(2))*dl**2+zabs(bm2(4))*dl
     &		 +zabs(bm2(5))
	if(l.gt.0) then
	yscalm(2)=yscalm(1)/dl**2
	yscalm(3)=yscalm(1)/dl
	else
	yscalm(1)=zabs(bm1(2))
	yscalm(2)=zabs(bm1(2))
	yscalm(3)=zabs(bm1(3))
	endif
	yscalm(4)=yscalm(3)
	yscalm(5)=yscalm(1)
c		write(6,*)'entering bsstep B'
c		write(6,*)'bm2=',bm2
c		write(6,*)'dydrm=',dydrm
c		write(6,*)'yscalm=',yscalm
	call bsstep(bm2,dydrm,5,r,htot,eps,yscalm,hdid,hnext)
c		write(6,*)'out of bsstep B'
	bmf1(i,j+1)=bm2(1)
	bmf2(i,j+1)=bm2(2)
	bmf3(i,j+1)=bm2(3)
	bmf4(i,j+1)=bm2(4)
	bmf5(i,j+1)=bm2(5)
	if((hdid-htot)/htot.ne.0.d0) iwk=1
c		write(6,*)'hdid,htot=',hdid,htot,'iwk=',iwk
c
132	r=rsav+dr    
131	continue
		write(6,*)'TOP layer',i,'bm2=',bm2
146	bmag=zabs(bm2(1))+real(l)*zabs(bm2(3))
	if(bmag.lt.1.e+10) go to 130
	do j1=1,5
	bm2(j1)=bm2(j1)*1.d-20
	enddo
	go to 146
130	continue
c
c	  write(6,*)'finished bm2'
c *	Now have available solutions b1 - b4 at source radius.
cc	Next require inverse of A at source radius, where dy/dr=Ay.
	ylim(1)=dsqrt((dble(2*l+1))/(4.d0*dble(pi)))
	ylim(2)=-0.5d0*ylim(1)*dsqrt(fl21)
	if(l.ne.0) ylim(3)=0.25d0*ylim(1)*dsqrt(dble((l-1)*(l+2))*fl21)
	  write(6,*)'ylim=',ylim
c	Do m=0, 1, 2 in that order.
c	June 29, 1999.  Then do vertical force, which is a m=0 deformation
c	represented by m1=4 below.
	m1=-1
38 	m1=m1+1
	if(m1.gt.4) go to 138
c	Avoid m1=2 and m1=3 cases for l=0
	if(l.eq.0.and.m1.eq.2) then
	write(6,*)'m1=',m1,'Skipping inversion because def=0'
	y1=0.d0
	y2=0.d0
	y3=0.d0
	y4=0.d0
	go to 25
	endif
	if(l.eq.0.and.m1.eq.3) then
	write(6,*)'m1=',m1,'Skipping inversion because def=0'
	y1=0.d0
	y2=0.d0
	y3=0.d0
	y4=0.d0
	go to 25
	endif
c	  write(6,*)'m1=',m1
c	Determine matrix elements of 6 X 6 boundary condition matrix bc0.
cc	First determine 4 X 4 a00 from eqn.'s (13)-(15) of notes.
cOLD	mus=dble(mu(isdep))
cOLD	lams=dble(kappa(isdep)-2.*mu(isdep)/3.)
cOLD	biga=lams+2.d0*mus
		r0=0.5d0*dble(rb(isdep)+rt(isdep))
		call interp(mus,lams,denss)
		biga=lams+2.d0*mus
	r=dble(rb(isdep))
	  if(m1.eq.0) write(6,*) isdep,r,dble(rt(isdep-1))
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
	go to 134
163	if(m1.eq.2) go to 151
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
		write(6,*)'line A: ylim(1),mus,biga=',ylim(1),mus,biga
	c0(1)=-2.d0*ylim(1)*(1.d0/(r*r))*(mus/biga-0.5d0)
		write(6,*)'line B: c0(1)=',c0(1)
	c0(2)=(-2.d0*mus/(r*r*r))*ylim(1)*(3.d0-4.d0*mus/biga)
	c0(3)=0.d0
	c0(4)=-0.5d0*c0(2)
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
c	Four rows of bc0 correspond to boundary conditions
c	at source.
134	continue

	c(1,2)=-bm1(3)
	c(2,2)=-bm1(5)
	c(1,4)=bm2(3)
	c(2,4)=bm2(5)
	if(l.gt.0) then
	c(1,1)=-bm1(1)
	c(2,1)=0.d0
	c(3,1)=bm1(4)
	c(4,1)=bm1(5)
	c(3,2)=(1.d0/fl21)*bm1(1)
	c(4,2)=0.d0
	c(1,3)=bm2(1)
	c(2,3)=0.d0
	c(3,3)=-bm2(4)
	c(4,3)=-bm2(5)
	c(3,4)=-(1.d0/fl21)*bm2(1)
	c(4,4)=0.d0
	endif

	do 135 i=1,4
	do 140 j=1,4
	bc0(j,i)=0.d0
	do 145 k=1,4
	bc0(j,i)=bc0(j,i)+a00(j,k)*c(k,i)
145	continue
140	continue
135	continue
c		Save the bc0 array before it's written over by ainver
		do ic1=1,4
		do jc1=1,4
		bc0s(ic1,jc1)=bc0(ic1,jc1)
		enddo
		c0s(ic1)=c0(ic1)
		enddo
c	err type strain (m1=0), ett+epp (m1=1),
c	ert,erp (m1=2), or ett-epp,etp (m1=3).
	  write(6,*)'m1=',m1,'entering ainver with c0=',c0
	if(l.eq.0) then
	bc1(1,1)=bc0(1,2)
	bc1(1,2)=bc0(1,4)
	bc1(2,1)=bc0(2,2)
	bc1(2,2)=bc0(2,4)
	c1(1)=c0(1)
	c1(2)=c0(2)
c		write(6,*)'bc0=',bc0
c		write(6,*)'bc1=',bc1
c		write(6,*)'c1=',c1
	call ainver(a1,bc1,c1,2)
	a0(1)=0.d0
	a0(2)=a1(1)
	a0(3)=0.d0
	a0(4)=a1(2)
	else	
	call ainver(a0,bc0,c0,4)
	endif
	  write(6,*)'lr,l=',lr,l,'a0=',a0
c - - - - - -
c	Determine displacements and their derivatives at surface.
	if(odep.gt.dble(dep)) go to 24
c	Do downward integration of `lower' displ-stress vector starting at source.
	do j=1,4
	blow(j)=-(a0(1)*c(j,1)+a0(2)*c(j,2))
	enddo
	write(6,*)'INITIAL blow=',blow
	do 530 i=isdep-1,1,-1
c	if(rb(i).gt.depmin) go to 530
c	if(i.lt.2) go to 530
cOLD	mus=dble(mu(i))
cOLD	lams=dble(kappa(i)-2.*mu(i)/3.)
cOLD	biga=lams+2.d0*mus
		r0=0.5d0*dble(rb(i)+rt(i))
		call interp(mus,lams,denss)
		biga=lams+2.d0*mus
	dr=-dble(rt(i)-rb(i))/dble(nrad)
c	Integrate equations of motion to propagate solutions downward from
c	the source.
	r=dble(rt(i))
c	Do blow using bsstep.
	iwk=0
	do 531 j=1,nrad
	rsav=r
	rintp=r+dr/2.d0
	call derivs4(r,blow,dydr)
	htot=dr
	eps=3.d-6
	if(l.gt.0) then
	yscal(1)=zabs(blow(1))+dl*zabs(blow(3))+zabs(blow(2))/dl+zabs(blow(4))
	yscal(2)=yscal(1)*dl
	yscal(3)=yscal(1)/dl
	yscal(4)=yscal(1)
	endif
	if(l.eq.0) then
	yscal(1)=zabs(blow(1))+zabs(blow(2))
	yscal(2)=yscal(1)
	yscal(3)=yscal(1)
	yscal(4)=yscal(1)
	endif
		if(yscal(1).eq.0.d0) then
		yscal(1)=1.d-10
		yscal(2)=1.d-10
		yscal(3)=1.d-10
		yscal(4)=1.d-10
		endif
c		write(6,*)'entering bsstep blow'
	call bsstep(blow,dydr,4,r,htot,eps,yscal,hdid,hnext)
	if((hdid-htot)/htot.ne.0.d0) iwk=1
	bm(1)=bmf1(i,nrad+1-j)
	bm(2)=bmf2(i,nrad+1-j)
	bm(3)=bmf3(i,nrad+1-j)
	bm(4)=bmf4(i,nrad+1-j)
	bm(5)=bmf5(i,nrad+1-j)
	valt=dlog(zabs(bm(4))+(1.d0/dl**2)*zabs(bm(1)))
     &	-dlog(zabs(bm(2))+(1.d0/dl**2)*zabs(bm(4)))-dlog(dl)
c			write(6,*)'before blow-bm relation: blow(1)=',blow(1)
c			write(6,*)'before blow-bm relation: blow(2)=',blow(2)
	if(valt.lt.1.5d0) then
c		write(6,*)'use 1st relation'
c	Update values of blow(2) and blow(4) using the relation between
c	the minors and the displacement-stress vector components.
	blow(2)=(bm(4)/bm(2))*blow(1)+(bm(1)/bm(2))*blow(3)
	blow(4)=(1.d0/(dl*(dl+1.d0)))*(bm(1)/bm(2))*blow(1)+
     &	(bm(3)/bm(2))*blow(3)
c		write(6,*)'used 1st relation'
	endif
	if(valt.ge.1.5d0) then
c		write(6,*)'use 2nd relation'
c	Update values of blow(1) and blow(3) using the relation between
c	the minors and the displacement-stress vector components.
	blow(1)=(bm(3)/bm(5))*blow(2)-(bm(1)/bm(5))*blow(4)
	blow(3)=-(1.d0/(dl*(dl+1.d0)))*(bm(1)/bm(5))*blow(2)+
     &	(bm(4)/bm(5))*blow(4)
c		write(6,*)'used 2nd relation'
	endif
c			write(6,*)'after blow-bm relation: blow(1)=',blow(1)
c			write(6,*)'after blow-bm relation: blow(2)=',blow(2)
c			write(6,*)'= = = = = = ='
c
	uslow=blow(1)/r
	vslow=blow(3)/r
	uplow=(blow(2)-(lams/r)*(2.d0*blow(1)-fl21*blow(3)))/biga 
c	up is equal to dy1/dr.
	vplow=(1.d0/r)*(blow(3)-blow(1))+blow(4)/mus 
c	vp is equal to dy3/dr.
cOLD	if(r.le.odep) go to 532
	if(i.eq.iodep.and.j.eq.nrad) then
	write(6,*)'IODEP4: store eigenfn at depth',bigr*(1.-r)
	urlow=uslow
	vrlow=vslow
	durlow=uplow
	dvrlow=vplow
	endif
532	r=rsav+dr    
531	continue
c	bmag=zabs(blow(1))+real(l)*zabs(blow(3))
c	if(bmag.lt.1.e+10) go to 530
cc	  write(6,*)'bmag=',bmag
c	do j1=1,4
c	blow(j1)=blow(j1)*1.d-10
c	enddo
c	urlow=urlow*1.d-10
c	vrlow=vrlow*1.d-10
c	durlow=durlow*1.d-10
c	dvrlow=dvrlow*1.d-10
530	continue
	write(6,*)'NEW blow AT LOWER BOUNDARY=',blow
	y1=durlow
	y2=urlow
	y3=dvrlow
	y4=vrlow
	go to 25
c - - - - - -
24	continue
c	Do upward integration of `upper' displ-stress vector starting at source.
	do j=1,4
	bupp(j)=a0(3)*c(j,3)+a0(4)*c(j,4)
	enddo
	write(6,*)'INITIAL bupp=',bupp
	do 430 i=isdep,n
c	if(rb(i).gt.depmin) go to 430
cOLD	mus=dble(mu(i))
cOLD	lams=dble(kappa(i)-2.*mu(i)/3.)
cOLD	biga=lams+2.d0*mus
		r0=0.5d0*dble(rb(i)+rt(i))
		call interp(mus,lams,denss)
		biga=lams+2.d0*mus
	dr=dble(rt(i)-rb(i))/dble(nrad)
c	Integrate equations of motion to propagate solutions upward from
c	the source.
	r=dble(rb(i))
c	Do bupp using bsstep.
	iwk=0
	do 431 j=1,nrad
	rsav=r
	rintp=r+dr/2.d0
	call derivs4(r,bupp,dydr)
	htot=dr
	eps=3.d-6
	if(l.gt.0) then
	yscal(1)=zabs(bupp(1))+dl*zabs(bupp(3))+zabs(bupp(2))/dl+zabs(bupp(4))
	yscal(2)=yscal(1)*dl
	yscal(3)=yscal(1)/dl
	yscal(4)=yscal(1)
	endif
	if(l.eq.0) then
	yscal(1)=zabs(bupp(1))+zabs(bupp(2))
	yscal(2)=yscal(1)
	yscal(3)=yscal(1)
	yscal(4)=yscal(1)
	endif
		if(yscal(1).eq.0.d0) then
		yscal(1)=1.d-10
		yscal(2)=1.d-10
		yscal(3)=1.d-10
		yscal(4)=1.d-10
		endif
c		write(6,*)'B calling bsstep; bupp=',bupp
c		write(6,*)'dydr=',dydr
	call bsstep(bupp,dydr,4,r,htot,eps,yscal,hdid,hnext)
c		write(6,*)'B out of bsstep'
	if((hdid-htot)/htot.ne.0.d0) iwk=1
	bm(1)=bmf1(i,nrad+1-j)
	bm(2)=bmf2(i,nrad+1-j)
	bm(3)=bmf3(i,nrad+1-j)
	bm(4)=bmf4(i,nrad+1-j)
	bm(5)=bmf5(i,nrad+1-j)
	valt=dlog(zabs(bm(4))+(1.d0/dl**2)*zabs(bm(1)))
     &	-dlog(zabs(bm(2))+(1.d0/dl**2)*zabs(bm(4)))-dlog(dl)
c			write(6,*)'before bupp-bm relation: bupp(1)=',bupp(1)
c			write(6,*)'before bupp-bm relation: bupp(2)=',bupp(2)
	if(valt.lt.1.5d0) then
c		write(6,*)'use 1st relation'
c	Update values of bupp(2) and bupp(4) using the relation between
c	the minors and the displacement-stress vector components.
	bupp(2)=(bm(4)/bm(2))*bupp(1)+(bm(1)/bm(2))*bupp(3)
	bupp(4)=(1.d0/(dl*(dl+1.d0)))*(bm(1)/bm(2))*bupp(1)+
     &	(bm(3)/bm(2))*bupp(3)
c		write(6,*)'used 1st relation'
	endif
	if(valt.ge.1.5d0) then
c		write(6,*)'use 2nd relation'
c	Update values of bupp(1) and bupp(3) using the relation between
c	the minors and the displacement-stress vector components.
	bupp(1)=(bm(3)/bm(5))*bupp(2)-(bm(1)/bm(5))*bupp(4)
	bupp(3)=-(1.d0/(dl*(dl+1.d0)))*(bm(1)/bm(5))*bupp(2)+
     &	(bm(4)/bm(5))*bupp(4)
c		write(6,*)'used 2nd relation'
	endif
c			write(6,*)'after bupp-bm relation: bupp(1)=',bupp(1)
c			write(6,*)'after bupp-bm relation: bupp(2)=',bupp(2)
c			write(6,*)'= = = = = = ='
c
	usupp=bupp(1)/r
	vsupp=bupp(3)/r
	upupp=(bupp(2)-(lams/r)*(2.d0*bupp(1)-fl21*bupp(3)))/biga 
c	up is equal to dy1/dr.
	vpupp=(1.d0/r)*(bupp(3)-bupp(1))+bupp(4)/mus 
c	vp is equal to dy3/dr.
cOLD	if(r.ge.odep) go to 432
	if(i.eq.(iodep-1).and.j.eq.nrad) then
	write(6,*)'IODEP3: store eigenfn at depth',bigr*(1.-r)
	urupp=usupp
	vrupp=vsupp
	durupp=upupp
	dvrupp=vpupp
	endif
432	r=rsav+dr    
431	continue
c	bmag=zabs(bupp(1))+real(l)*zabs(bupp(3))
c	if(bmag.lt.1.e+10) go to 430
cc	  write(6,*)'bmag=',bmag
c	do j1=1,4
c	bupp(j1)=bupp(j1)*1.d-10
c	enddo
c	urupp=urupp*1.d-10
c	vrupp=vrupp*1.d-10
c	durupp=durupp*1.d-10
c	dvrupp=dvrupp*1.d-10
430	continue
	write(6,*)'NEW bupp AT SURFACE=',bupp
	y1=durupp
	y2=urupp
	y3=dvrupp
	y4=vrupp
c- - - - - -
25	write(2) l,y1,y2,y3,y4
c	TEST the solution of the inverse problem.
	row1=bc0s(1,1)*a0(1)+bc0s(1,2)*a0(2)+bc0s(1,3)*a0(3)+bc0s(1,4)*a0(4)
	write(6,*)'row1,c0s(1)=',row1,c0s(1)
	row2=bc0s(2,1)*a0(1)+bc0s(2,2)*a0(2)+bc0s(2,3)*a0(3)+bc0s(2,4)*a0(4)
	write(6,*)'row2,c0s(2)=',row2,c0s(2)
	row3=bc0s(3,1)*a0(1)+bc0s(3,2)*a0(2)+bc0s(3,3)*a0(3)+bc0s(3,4)*a0(4)
	write(6,*)'row3,c0s(3)=',row3,c0s(3)
	row4=bc0s(4,1)*a0(1)+bc0s(4,2)*a0(2)+bc0s(4,3)*a0(3)+bc0s(4,4)*a0(4)
	write(6,*)'row4,c0s(4)=',row4,c0s(4)
c	write(6,*)' '
	bfacl=0.d0
	bfacm=0.d0
	write(2) bfacl,bfacm
		write(6,*)'bfacl,bfacm=',bfacl,bfacm
	go to 38
138	continue
	  write(6,*)'starting TOROIDAL modes'
c	Don't do l=0 case for toroidal modes.
	if(lr.eq.(lmax+1)) go to 50
		write(6,*)'A lr,l=',lr,l,'iom=',iom
c	TOROIDAL MODES
c	Start matrix propagation upward for solution t1.  Stop
c	at layer [isdep-1].
	is=0
c	depmax=dep-2.5*(2.*pi/(real(l)+0.5))
c	depmin=dep+2.5*(2.*pi/(real(l)+0.5))
	depmax=rb(1)
	depmin=1.01
	do 230 i=1,isdep-1
c	if(i.lt.2) go to 230
c	Start integration at most 2 wavelengths into the earth.
	if(rt(i).lt.depmax.or.rb(i).gt.depmin) go to 230
cOLD	if(dble(depmax).ge.odep.or.dble(depmin).le.odep) go to 230
	if(rt(i).lt.depmax) go to 230
c	Give boundary conditions at top of layer (i-1) appropriate for
c	homogeneous solid in 0 < r < rt(i-1).
	if(is.ne.0) go to 228
	is=1
cOLD	  t1(1)=1.d0
cOLD	  t1(2)=0.d0
	r0=dble(rb(i))
	csys='t'
c		write(6,*)'entering lbound - TOR'
	call lbound(t1,2,csys)
c		write(6,*)'out of lbound - TOR: t1=',t1
228	continue
	dr=dble(rt(i)-rb(i))/dble(nrad)
	r=dble(rb(i))
		r0=0.5d0*dble(rb(i)+rt(i))
		call interp(mus,lams,denss)
	iwk=0
	do 231 j=1,nrad
	rsav=r
c		write(6,*)'entering derivs2'
	rintp=r+dr/2.d0
	call derivs2(r,t1,dtdr)
	htot=dr
	eps=3.d-6
	yscalt(1)=zabs(t1(1))+zabs(t1(2))/dl
	yscalt(2)=yscalt(1)*dl
c		write(6,*)'CC TOR'
	call bsstep(t1,dtdr,2,r,htot,eps,yscalt,hdid,hnext)
c		write(6,*)'DD TOR'
	if((hdid-htot)/htot.ne.0.d0) iwk=1
		if(iwk.eq.1) then
		write(6,*)'iwk=1 TOR l=',l
		stop
		endif
c
cOLD	if(r.ge.odep) go to 232
	if(i.eq.(iodep-1).and.j.eq.nrad) then
	write(6,*)'IODEP1: store eigenfn at depth',bigr*(1.-r)
	wr1=t1(1)/r
	dwr1=t1(1)/r+t1(2)/mus
	endif
232	r=rsav+dr    
231	continue
c		write(6,*)'TOP layer',i,'t1=',t1
c		r0=r
c		call lbound(t1try,2,csys)
c		write(6,*)'lbound: scaled t1='
c		do jtry=1,2
c		write(6,*) t1try(jtry)*t1(2)/t1try(2)
c		enddo
	  bmag=zabs(t1(1)/r)+zabs(t1(2)/mus)
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
c	Start matrix propagation downward for solutions t.  Stop
c	at layer [n].
	t2(1)=1.d0
	t2(2)=0.d0
	wr2=t2(1)
	dwr2=wr2
	do 330 i=n,isdep,-1
	if(rb(i).gt.depmin) go to 330
		r0=0.5d0*dble(rb(i)+rt(i))
		call interp(mus,lams,denss)
	dr=-dble(rt(i)-rb(i))/dble(nrad)
	r=dble(rt(i))
	iwk=0
	do 331 j=1,nrad
	rsav=r
c		write(6,*)'entering derivs2'
	rintp=r+dr/2.d0
	call derivs2(r,t2,dtdr)
	htot=dr
	eps=3.d-6
	yscalt(1)=zabs(t2(1))+zabs(t2(2))/dl
	yscalt(2)=yscalt(1)*dl
	call bsstep(t2,dtdr,2,r,htot,eps,yscalt,hdid,hnext)
	if((hdid-htot)/htot.ne.0.d0) iwk=1
c
		if(iwk.eq.1) then
		write(6,*)'iwk=1 TOR l=',l
		stop
		endif
cOLD	if(r.lt.odep) go to 332
	if(i.eq.iodep.and.j.eq.nrad) then
	write(6,*)'IODEP2: store eigenfn at depth',bigr*(1.-r)
	wr2=t2(1)/r
	dwr2=t2(1)/r+t2(2)/mus
	endif
332	r=rsav+dr    
331	continue
	  bmag=zabs(t2(1)/r)+zabs(t2(2)/mus)
	if(bmag.lt.1.e+10) go to 330
	t2(1)=t2(1)*1.d-10
	t2(2)=t2(2)*1.d-10
	wr2=wr2*1.d-10
	dwr2=dwr2*1.d-10
330	continue
c	Store these solutions at source radius.
	do 227 k=1,2
	t2s(k)=t2(k)
227	continue
c
c	Now have available t1 at source radius, and t2 at surface.
cc	Next require inverse of A, where dy/dr=Ay.
cOLD	mus=dble(mu(isdep))
		r0=0.5d0*dble(rb(isdep)+rt(isdep))
		call interp(mus,lams,denss)
	r=dble(rb(isdep))
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
	go to 234
252	a00(1,1)=r*r
	a00(1,2)=0.d0
	a00(2,1)=0.d0
	a00(2,2)=1.d0
	c1(1)=0.d0
	c1(2)=-0.5d0*ylim(1)*
     &	dsqrt(dble((l+2)*(l-1))/dble(l*(l+1)))*mus/(r*r*r)
c	First two rows of bc1 correspond to boundary conditions at source.
234	continue
	do j=1,2
	c(j,1)=-t1s(j)
	c(j,2)=t2s(j)
	enddo
	do 235 i=1,2
	do 240 j=1,2
	bc1(j,i)=0.d0
	do 245 k=1,2
	bc1(j,i)=bc1(j,i)+a00(j,k)*c(k,i)
245	continue
240	continue
235	continue
c	ert, erp type strain (m1=1) or ett-epp, etp (m1=2).
		write(6,*)'TOR: entering ainver: c1=',c1
		write(6,*)'TOR: entering ainver: bc1='
		write(6,*) bc1(1,1),bc1(1,2)
		write(6,*) bc1(2,1),bc1(2,2)
	call ainver(a1,bc1,c1,2)
		write(6,*)'TOR: out of ainver: a1=',a1
c	Determine displacement and their derivatives at surface.
	  write(6,*)'odep,dep=',odep,dep
	if(odep.gt.dble(dep)) go to 124
	y1=(a1(1)*(dwr1-wr1))
	y2=(a1(1)*wr1)
	go to 125
124	y1=(a1(2)*(dwr2-wr2))
cc	ABOVE LINE is zero when evaluated at surface.  So to avoid underflow,
cc	just set it equal to zero.
c	y1=0.
	y2=(a1(2)*wr2)
125	write(2) l,y1,y2
	  wb=(a1(1)*t1s(1))
	  sb=(a1(1)*t1s(2))
	  wa=(a1(2)*t2s(1))
	  sa=(a1(2)*t2s(2))
c	  write(6,*)'t1s=',t1s
c	  write(6,*)'t2s=',t2s
c	  write(6,*)'t2=',t2
c	  write(6,*)'a1=',a1
	  write(6,*)'w,sh just above=',wa,sa
	  write(6,*)'w,sh just below=',wb,sb
c	  ssur=real(a1(2)*t2(2)+a1(3)*t3(2))
c	  write(6,*)'surface traction=',ssur
	go to 238
50	continue
100	continue
	close(2)
	end 

	subroutine ainver(a,b,c,n)
c       Find solution (a) of the matrix equation (B)(a)=(c).
c       Method:  Gaussian elimination.
	integer*2 perm 
	dimension perm(10)
	double complex a(n),b(n,n),c(n),bsave,fac
c		if(n.eq.6) write(6,*) 'AINVER: b=',b
        do 5 i=1,n
5       perm(i)=i
	i=0
10      i=i+1
	if (i.gt.n) go to 35
c       Find maximum in row i.
	amax=0.d0
	imax=i
	do 15 j=i,n
	t=zabs(b(i,j))
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

