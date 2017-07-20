	subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout)
cOLD	implicit real*8 (a-h,o-z)
	parameter (nmax=10)
cOLD	dimension y(nvar),dydx(nvar),yout(nvar),ym(nmax),yn(nmax)
c--
      double complex y(nvar),dydx(nvar),yout(nvar),ym(nmax),yn(nmax),swap
      real*8 xs,htot,h,h2,x
c--
c		write(6,*)'mmid: nvar=',nvar
c		write(6,*)'dydx=',dydx
	h=htot/dble(nstep)
	do 11 i=1,nvar
	  ym(i)=y(i)
	  yn(i)=y(i)+h*dydx(i)
11	continue
	x=xs+h
	if(nvar.eq.2) call derivs2(x,yn,yout)
c		if(nvar.eq.4) then
c		write(6,*)'mmid: will call derivs4 with x=',x
c		write(6,*)'yn=',yn
c		endif
	if(nvar.eq.4) call derivs4(x,yn,yout)
	if(nvar.eq.5) call derivs5(x,yn,yout)
c		write(6,*)'mmid: yout=',yout
	h2=2.d0*h
	do 13 n=2,nstep
	do 12 i=1,nvar
	  swap=ym(i)+h2*yout(i)
	  ym(i)=yn(i)
	  yn(i)=swap
12	continue
	  x=x+h
	if(nvar.eq.2) call derivs2(x,yn,yout)
	if(nvar.eq.4) call derivs4(x,yn,yout)
	if(nvar.eq.5) call derivs5(x,yn,yout)
13	continue
	do 14 i=1,nvar
	yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14	continue
	return
	end

