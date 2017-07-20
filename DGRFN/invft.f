	subroutine invft(nleno,dt,obeta,corper)
c	Obtain inverse Fourier transform of spectrum in xo-array with
c	time increment dt.  The
c	Filter the time series xo of length nleno and time step dt
c	using corner period corper.
c	Result output into array xt.
c
c	The input time series is specified in the frequency domain
c	at freq. samples that are a distance obeta below the real omega axis.
c	The input spectrum is assumed to be sampled in intervals
c	of angular frequency equal to ommin, as specified below.
c
	parameter (iommax=4096)
	complex*8 xo
	complex*8 xval,ui
	real*4 xt,ommin,cost0,cosdt0,sint0,sindt0
	real*4 facf(iommax)
	common/spec/xo(iommax),xt(4*iommax)
	common/invtyp/ityp
	ui=cmplx(0.,1.)
	pi=3.14159265
c**	take inverse FT with corner freq. at period [corper] sec.
c	(low pass filter).
	ommin=2.*pi/(dt*real(nleno))
	jmax1=int((2.*pi*2./corper)/real(ommin))
c		write(6,*)'invft: jmax1=',jmax1
	jmax0=jmax1/2
c	Use a high-pass filter at period 400 sec
c	jmax2=int((2.*pi*2./400.)/real(ommin))
c		write(6,*)'jmax2=',jmax2
	do 43 j=1,jmax1
	facf(j)=1.0
	if(j.eq.1) facf(j)=0.5
c	Above line: we pick up zero frequency only once, so
c	want it multiplied by (1/(2*pi)) instead of (1/pi) in j-loop below.
c---------low-pass filter
	if(j.gt.jmax0) facf(j)=0.5-0.5*cos((real(j-jmax1)/real(
     &	jmax1-jmax0))*(pi))
c---------high-pass filter
c	if(j.lt.jmax2) facf(j)=0.5+0.5*cos((real(jmax2-j)/real(jmax2-1))
c     &	*(pi))
c---------
43	continue
	de0=exp(obeta*0.25*dt)
	e0=1./de0
	do 46 k=1,4*nleno
	xt(k)=0.0
	cosdt0=cos(0.25*dt*real(k)*ommin)
	sindt0=sin(0.25*dt*real(k)*ommin)
	cost0=cosdt0
	sint0=-sindt0
	do 40 j=1,jmax1
	xval=xo(j)
	cossav=cost0
	cost0=cost0*cosdt0-sint0*sindt0
	sint0=cossav*sindt0+sint0*cosdt0
cDISPL	Displacement response.
	if(ityp.eq.0) then
	xt(k)=xt(k)+(1./pi)*(real(xval)*cost0-aimag(xval)*sint0)*
     &	ommin*facf(j)
	endif
cVEL	Velocity response.  Mult bu ui*(real(j-1)*ommin - ui*obeta) 
	if(ityp.eq.1) then
	xt(k)=xt(k)+(1./pi)*((real(xval)*cost0-aimag(xval)*sint0)*obeta
     &	-(aimag(xval)*cost0+real(xval)*sint0)*real(j-1)*ommin)*
     &	ommin*facf(j)
	endif
40	continue
c	Restore the exponential factor arising from integrating
c	below the real omega axis.
	e0=e0*de0
	xt(k)=xt(k)*e0
46	continue
c* * *
	return
	end
