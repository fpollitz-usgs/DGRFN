	subroutine lbound(zyst,j,csys)
c
* Calculates STARTING VALUES for integration of the spheroidal
* and toroidal motion upwards to the source level.
* The values are derived from a homogeneous isotropic Earth model
* assumed beneath the starting radius r1. The results are
* spherical Bessel functions. To avoid explicit calculation of
* those functions, we use ratios of them. These ratios can 
* be computed from a recurrence relation. See subroutine ZETL().
c Author: J.R. Dalkolmo
c
c	This subroutine is modified from that written by J.R. Dalkolmo
c	(Fred Pollitz, October, 2010)
c
	character*1 csys
	real*4 kappa,mu 
      double complex zyst(j),zom,zxa2,zeta,zxb2,zetb,zl2m,zmue,za,zb
	double complex w0
	double complex mus,lams
	real*8 denss,rmus
	real*8 rterm,sterm
      real*8 r0,elp1,xfi,ro1,el,rsq
	real*8 rintp,fac
	real*8 tiny
	parameter (tiny=1.e-20)
	common/radlev/rintp,bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),w0,r0
	common/qval/qbet(200),qkap(200)
	common/aterm/rterm,sterm
	common/matl/l
c	Obtain complex-valued elastic parameters at radius r0
	call interp(mus,lams,denss)
cTE	To crudely implement fluid outer core, just reset mu=0 here.
c	lams=65.0d0
c	denss=10.d0
c	mus=0.d0
c--
	rmus=dreal(mus)
	fac=dble(sqrt(10.)/bigr)
c * * *
c	Re-cast needed values in terms of Dalkolmo's nomenclature.
	zmue=mus
	zl2m=lams+2.d0*mus
	xfi=r0
	ro1=denss
	zom=w0/fac
	elp1 = dble(l*(l+1))
c		write(6,*)'TEST in lbound'
c		write(6,*)'zmue=',zmue
c		write(6,*)'zl2m=',zl2m
c		write(6,*)'xfi=',r0
c		write(6,*)'ro1=',ro1
c		write(6,*)'zom=',w0
c		write(6,*)'elp1=',elp1
c		write(6,*)'csys=',csys
c * * *
      rsq = xfi*xfi
      el = dble(l)

	

      if(csys.eq.'s') then            !   Spheroidal motion

       zxa2=ro1*zom*zom*rsq/zl2m
c		write(6,*)'about to enter zetl for zeta'
c		write(6,*)'ro1,zom,rsq,zl2m=',ro1,zom,rsq,zl2m
c      print *, 'zxa2=',zxa2
c		write(6,*)'doing zeta'
       call zetl(zxa2,l,zeta)   ! Ratio of spherical Bessel functions.
c		write(6,*)'zeta=',zeta
       if(rmus.gt.tiny) then
         zxb2=ro1*zom*zom*rsq/zmue
c		write(6,*)'doing zetb'
         call zetl(zxb2,l,zetb)
c		write(6,*)'zetb=',zetb
       endif

c Starting values for integration UPWARDS in liquid or solid layer.

      if (l.eq.0) then
        zyst(1)=-zeta
        zyst(2)=-ro1*xfi*zom*zom + 4.d0*zmue/xfi*zeta
	zyst(3)=0.d0
	zyst(4)=0.d0
	zyst(5)=0.d0
        return
      else
        if(rmus.lt.tiny) then
          zyst(1) = -(1.d0-zeta/el)/(zom*zom*ro1*xfi)
          zyst(2) = 1.d0/el
	  zyst(3)=0.d0
	  zyst(4)=0.d0
	  zyst(5)=0.d0
          return 
        else
          za=zeta/el
          zb=zetb/(el+1.d0)
          zyst(2)=(-za+zb*(za-1.d0))/el
          zyst(1)=zmue/xfi*(zxb2/el+2.d0*elp1*zyst(2))
          zyst(3)=zmue/xfi*(-2.d0*zyst(2)+zxb2/elp1*(za-1.d0))
          zyst(4)=-zmue/xfi*(zxb2/(el*el)*(1.d0-zb)+4.d0*zyst(2))
          zyst(5)=zmue*zmue/rsq*(-4.d0*(el-1.d0)*(el+2.d0)*zyst(2)+
     1       zxb2/el*(zxb2/elp1-2.d0*(el-1.d0)*(2.d0*el+1.d0)/elp1
     2                       -4.d0/(el+1.d0)*za-2.d0/el*zb))
c          print *,(zyst(i),i=1,5)
        endif
       endif

      else if(csys.eq.'t') then   !    Toroidal motion
      
       if (rmus.gt.tiny) then
          zxb2=ro1*zom*zom*rsq/zmue
          call zetl(zxb2,l,zetb)
          zyst(1) = 1.d0/el
c         print *, zetb
          zyst(2) = zmue/(xfi*el)*(el-1.d0-zetb)
       else         ! liquid layer
          zyst(1) = 1.d0/el
          zyst(2) = 0.d0
       endif

      else
       stop 'Error in <stavani>: No motion specified!'
      endif

      return
      end 

      Subroutine   Zetl (zx2,l,z)
c
c Computes ratio of (complex) spherical Bessel-functions 
c    z(n):= x*j(n+1)/j(n)
c
c by calculation of the recurrence formula
c    z(n)=x**2/(2*n+3 -z(n+1)).
c (Takeuchi and Saito (1972), Methods in computional Physics, p.243)
c
c The recurrence formula is treated as a continued fraction
c and evaluated with the MODIFIED LENTZ'S METHOD. See Press, W.H. et
c al, 1992, Numerical Recipes, 2nd Edition, Cambridge, Chapter 5.2,
c p. 169 for further information.
c
c Author: J.R. Dalkolmo


	real*8 tiny,eps
	integer maxit,l,i
      parameter(maxit=18000,tiny=1.d-30)
	double complex zx2,z,zd,zc,zdelta
      common /acc/ eps
	save
c	Next line from Fred Pollitz; we don't use common/acc/ 
	eps=3.d-6
c First iteration
      z      = tiny
      zd     = dcmplx(2*l+3)
      zc     = zd+zx2/z
      zd     = 1.d0/zd
      zdelta = zc*zd
      z      = z*zdelta

c Remaining iterations
      do 10 i=2,maxit
		zd=dcmplx(2*(l+i)+1)-zx2*zd
		if (zd.eq.0.d0) zd=tiny
		zc=dcmplx(2*(l+i)+1)-zx2/zc
		if (zc.eq.0.d0) zc=tiny
		zd=1.d0/zd
		zdelta=zc*zd
		z=z*zdelta
c	write(6,*)'iteration #',i,'zdelta=',zdelta
cORIG		if ((zabs(zdelta)-1.d0).lt.eps) then
c	The above line is the original line in the gemini-2.2 release.
c	It is a bug and should be replaced with the following line.
c	(Fred Pollitz, October, 2010)
		if ((zabs(zdelta-1.d0)).lt.eps) then
c	write(6,*)'Zetl: zabs(zdelta)=',zabs(zdelta),'eps=',eps
c	write(6,*)'zdelta=',zdelta
	write(6,*)'i=',i,'out of max of ',maxit
c	  print *, zsqrt(zx2),l,i
		  return
		endif
 10   continue
      print *, 'zx2,l,eps,z,maxit,tiny=',zx2,l,eps,z,maxit,tiny
      pause '*** Iteration failed to converge in ZETL ! ***'

      end

	subroutine derivs2(x,yn,yout)
	real*8 r0,r2,facl,fac
	real*8 x,rintp
	double complex w0,om2,yn(2),yout(2),aj(2,2)
	double complex mus,lams
	real*8 denss,rterm,sterm
	real*4 kappa,mu
	common/radlev/rintp,bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),w0,r0
	common/qval/qbet(200),qkap(200)
	common/aterm/rterm,sterm
	common/matl/l
	fac=dble(sqrt(10.)/bigr)
	om2=(w0/fac)*(w0/fac)
c       fac = sqrt (10**11/(R**2)), where R = bigr*10**5 cm.
c       this converts w0 from rad/sec to "internal" units.
	r0=rintp
c		write(6,*)'entering DERIVS2: r0=',r0
	call interp(mus,lams,denss)
c		write(6,*)'out of DERIVS2: r0=',r0,'mus=',mus,'l=',l
	r0=x
	r2=r0*r0
	facl=dble((l-1)*(l+2))
c	facl=(l-1.d0)*(l+2.d0)
	  aj(1,1)=1.d0/r0
	  aj(1,2)=1.d0/mus
  	  aj(2,1)=facl*mus/r2-denss*om2
	  aj(2,2)=-3.d0/r0
	call prod(aj,yn,2,yout) 
	return
	end

	subroutine derivs4(x,yn,yout)
	real*8 r0,rterm,sterm
	real*8 x,rintp
	double complex w0,yn(4),yout(4),aj(4,4)
	real*4 kappa,mu
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),w0,r0
	common/qval/qbet(200),qkap(200)
	common/aterm/rterm,sterm
	common/matl/l
	common/radlev/rintp,bigr
c		write(6,*)'derivs4: x=',x
c		write(6,*)'derivs4: yn=',yn
	r0=x
c		write(6,*)'entering amat4: r0=',r0
	call amat4(aj,1)
c		write(6,*)'derivs4: out of amat4: aj=',aj
	call prod(aj,yn,4,yout) 
c		write(6,*)'derivs4: out of prod: dydr=',yout
	return
	end

	subroutine derivs5(x,yn,yout)
	real*8 r0,rterm,sterm
	real*8 x,rintp
	double complex w0,yn(5),yout(5),aj(5,5)
	real*4 kappa,mu
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),w0,r0
	common/qval/qbet(200),qkap(200)
	common/aterm/rterm,sterm
	common/matl/l
	common/radlev/rintp,bigr
c		write(6,*)'derivs5: x[input]=',x
	r0=x
c		write(6,*)'derivs5: yn=',yn
c		write(6,*)'derivs5: r0=',r0,'rintp=',rintp
c		write(6,*)'derivs5: qbet(1),qbet(2)=',qbet(1),qbet(2)
	call amat5(aj,1)
c		write(6,*)'out of amat5: aj=',aj
	call prod(aj,yn,5,yout) 
c		write(6,*)'derivs5: yout=',yout
c		pause
	return
	end

    	subroutine amat4(aj,i)
c 	Use to determine A(r) (returned in aj) in layer at radius r=r0.
c	dY(r)/dr = A(r)Y.
	real*4 kappa,mu 
	double complex w0,om2,aj(4,4)
	real*8 fac
	double complex mus,lams,biga
	real*8 rintp
	real*8 fl21,r,r0,denss,r2,rterm,sterm
	common/radlev/rintp,bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),w0,r0
	common/qval/qbet(200),qkap(200)
	common/aterm/rterm,sterm
	common/matl/l
	fac=dble(sqrt(10.)/bigr)
	om2=(w0/fac)*(w0/fac)
cOLD	om2=(w0/4.96355d-4)*(w0/4.96355d-4)
c       4.96355e-4 = sqrt (10**11/(R**2)), where R = 6371*10**5 cm.
c       this converts w0 from rad/sec to "internal" units.
c	mus=dble(mu(i))
c	lams=dble(kappa(i)-2.*mu(i)/3.)
c	denss=dble(dens(i)) 
c		write(6,*)'amat4: r0=',r0,'rintp=',rintp
	r=r0
	r0=rintp
	call interp(mus,lams,denss)
	r0=r
	biga=lams+2.d0*mus
	fl21=dble(l*(l+1))
	r=r0
	r2=r*r
	aj(1,1)=-2.d0*lams/(r*biga)
	aj(1,2)=1.d0/biga
	aj(1,3)=fl21*lams/(r*biga)
	aj(1,4)=0.d0
	aj(2,1)=(4.d0*(lams+mus)/r)/r+(2.d0/r)*lams*aj(1,1)
     &	-om2*denss  
	aj(2,2)=-2.d0/r+(2.d0/r)*lams*aj(1,2)
	aj(2,3)=-fl21*(2.d0*(lams+mus)/r)/r+(2.d0/r)*lams*aj(1,3)
	aj(2,4)=fl21/r
	aj(3,1)=-1.d0/r
	aj(3,2)=0.d0
	aj(3,3)=1.d0/r
	aj(3,4)=1.d0/mus
	aj(4,1)=-(lams/r)*aj(1,1)+2.d0*(mus-biga)/r2 
	aj(4,2)=-(lams/r)*aj(1,2)
	aj(4,3)=-(lams/r)*aj(1,3)+(fl21*biga-2.d0*mus)/r2-om2*denss 
	aj(4,4)=-3.d0/r
c	aj contains A(r).  
	return
	end 

    	subroutine amat5(aj,i0)
c 	Use to determine A(r) (returned in aj) in layer at radius r=r0.
c	dY(r)/dr = A(r)Y.
	real*4 kappa,mu 
	double complex w0,om2,aj(5,5)
	real*8 fac
	double complex mus,lams,biga
	real*8 rintp
	real*8 fl21,r,r0,denss,r2,rterm,sterm
	common/radlev/rintp,bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),w0,r0
	common/qval/qbet(200),qkap(200)
	common/aterm/rterm,sterm
	common/matl/l
	fac=dble(sqrt(10.)/bigr)
	om2=(w0/fac)*(w0/fac)
cOLD	om2=(w0/4.96355d-4)*(w0/4.96355d-4)
c	4.96355e-4 = sqrt (10**11/(R**2)), where R = 6371*10**5 cm.
c	this converts w0 from rad/sec to "internal" units.
c	mus=dble(mu(i0))
c	lams=dble(kappa(i0)-2.*mu(i0)/3.)
c	denss=dble(dens(i0)) 
	r=r0
	r0=rintp
	call interp(mus,lams,denss)
c		write(6,*)'amat5: out of interp: mus,lams=',mus,lams
c		if(r0.ge.0.999572) write(6,*)'out of interp: r0=',r0,'denss=',denss
	r0=r
c		write(6,*)'amat5: mus,lams,denss=',mus,lams,denss
c		pause
	biga=lams+2.d0*mus
c	fl21=l*(l+1.d0)
	fl21=dble(l*(l+1))
	r=r0
	r2=r*r
	aj(1,1)=-2.d0/r
	aj(1,2)=-2.d0*fl21*(biga-lams*lams/biga-mus)/r2
	aj(1,3)=fl21/r
	aj(1,4)=-fl21*(lams/biga)/r
	aj(1,5)=0.d0
	aj(2,1)=0.d0
	aj(2,2)=(1.d0-2.d0*lams/biga)/r
	aj(2,3)=1.d0/mus
	aj(2,4)=1.d0/biga
	aj(2,5)=0.d0
	aj(3,1)=-2.d0*(lams/biga)/r
	aj(3,2)=-om2*denss+fl21*(biga-lams*lams/biga)/r2-2.d0*mus/r2
	aj(3,3)=-(3.d0+2.d0*lams/biga)/r
	aj(3,4)=0.d0
	aj(3,5)=1.d0/biga
	aj(4,1)=2.d0/r
	aj(4,2)=-om2*denss+4.d0*(biga-lams*lams/biga-mus)/r2
	aj(4,3)=0.d0
	aj(4,4)=-(1.d0-2.d0*lams/biga)/r
	aj(4,5)=1.d0/mus
	aj(5,1)=4.*(biga-lams*lams/biga-mus)/r2
	aj(5,2)=0.d0
	aj(5,3)=-om2*denss+4.*(biga-lams*lams/biga-mus)/r2
	aj(5,4)=-om2*denss+fl21*(biga-lams*lams/biga)/r2-2.d0*mus/r2
	aj(5,5)=(2.d0*lams/biga-5.d0)/r
c	aj contains A(r).  
c	  write(6,*)'B: aj=',aj
c	  pause
90	return
	end 

	subroutine prod(a,c,n,b)
c	Form matrix product b=(A)c
	double complex a,b,c 
	dimension a(n,n),c(n),b(n)
	do 5 i=1,n
	b(i)=0.d0
	do 4 j=1,n
4	b(i)=b(i)+a(i,j)*c(j)
5	continue
	return
	end 

