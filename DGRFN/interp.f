	subroutine interp(mus,lams,denss)
c	Give back mus ,lams, denss at radius r0 interpolated between layers.
	double complex w0
	double complex mus,lams,kaps
	real*8 denss,r0,f,rterm,sterm
	real*8 pi,twopi
	real*4 mu,kappa
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),w0,r0
	common/qval/qbet(200),qkap(200)
	common/aterm/rterm,sterm
	common/matl/l
	pi=3.14159265358979d0
	j=0
5	j=j+1
c		if(j.eq.34) write(6,*)'INTERP34: j=',j,'r0=',r0,'rt(',j,')=',rt(j)
c		if(j.eq.35) write(6,*)'INTERP35: j=',j,'r0=',r0,'rt(',j,')=',rt(j)
	if(rt(j).lt.real(r0)) go to 5
c		write(6,*)'INTERP: j=',j,'r0=',r0,'rt(',j,')=',rt(j)
cOLD	f=(r0-dble(rb(j)))/dble(rt(j)-rb(j))
c		write(6,*)'INTERP: qbet(',j,')=',qbet(j)
c	Assume a reference frequency of twopi rad/sec (i.e. 1 Hz).
c	This is done when computing the real and imaginary parts of log(w0/twopi)
c	in wave0s.f, which are in rterm and sterm, respectively.
	mus=dble(mu(j)) * 
     &	dcmplx(1.d0+(2.d0/pi)*dble(qbet(j))*rterm,dble(qbet(j))*(1.d0+(2.d0/pi)*sterm))
c		write(6,*)'w0=',w0,'rterm,sterm=',rterm,sterm,'mus=',mus
c		write(6,*)'qbet(',j,')=',qbet(j)
c		write(6,*)'(1.d0+(2.d0/pi)*sterm)=',(1.d0+(2.d0/pi)*sterm)
	kaps=dble(kappa(j)) * 
     &	dcmplx(1.d0+(2.d0/pi)*dble(qkap(j))*rterm,dble(qkap(j))*(1.d0+(2.d0/pi)*sterm))
	lams=kaps-(2.d0/3.d0)*mus
	denss=dble(dens(j))
	return
	end
