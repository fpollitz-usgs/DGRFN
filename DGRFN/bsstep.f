	subroutine bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext)
cOLD	implicit real*8 (a-h,o-z)
	real*8 shrink,grow
	parameter (nmax=10,imax=14,nuse=7,one=1.d0,shrink=0.95d0,grow=1.2d0)
	double complex y(nv),dydx(nv),yerr(nmax),
     &	ysav(nv),dysav(nv),yseq(nv)
	real*8 h,htry,x,xsav,eps,yscal,hdid,hnext,xest,errmax
	real*8 rintp
	dimension yscal(nv)
	dimension nseq(imax)
	common/radlev/rintp,bigr
	data nseq/2,4,6,8,12,16,24,32,48,64,96,128,192,256/
	h=htry
	xsav=x
	do 11 i=1,nv
	  ysav(i)=y(i)
	  dysav(i)=dydx(i)
11	continue
1	do 10 i=1,imax
c		write(6,*)'BSSTEP: entering mmid -- rintp=',rintp
c		write(6,*)'BSSTEP: entering mmid - dysav=',dysav
	if(nv.eq.2) call mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq)
	if(nv.eq.4) call mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq)
	if(nv.eq.5) call mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq)
c		write(6,*)'bsstep: out of mmid: yseq=',yseq
	  xest=(h/(dble(nseq(i))))**2
	  call rzextr(i,xest,yseq,y,yerr,nv,nuse)
c		write(6,*)'out of rzextr: yerr=',yerr
c		write(6,*)'yscal=',yscal
	  errmax=0.d0
	do 12 j=1,nv
cOLD	  t=zabs(yerr(j)/yscal(j))
cOLD	  if(errmax.lt.t) errmax=t
	  errmax=max(errmax,zabs(yerr(j)/yscal(j)))
c		write(6,*)'j=',j,yerr(j)/yscal(j)
12	continue
c		write(6,*)'** i=',i,' ** errmax=',errmax
c		write(6,*)'eps=',eps
c		write(6,*)'  '
	  errmax=errmax/eps
c		write(6,*)'NEW errmax=',errmax,'one=',one
	  if(errmax.lt.one) then
c		write(6,*)'yerr=',yerr,'yscal=',yscal
	  	x=x+h
		hdid=h
		if(i.eq.nuse) then
		  hnext=h*shrink
		else if(i.eq.nuse-1) then
		  hnext=h*grow
		else
		  hnext=(h*dble(nseq(nuse-1)))/dble(nseq(i))
		endif
		return
	  endif
10	continue
	h=0.25d0*h/dble(2**((imax-nuse)/2))
	if(x+h.eq.x) pause 'step size underflow'
	go to 1
	end	
