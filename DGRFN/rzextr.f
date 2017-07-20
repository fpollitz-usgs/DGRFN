	subroutine rzextr(iest,xest,yest,yz,dy,nv,nuse)
cOLD	implicit real*8 (a-h,o-z)
	parameter (imax=14,nmax=10,ncol=7)
	double complex yest(nv),yz(nv),dy(nmax),d(nmax,ncol),v,c,b,b1,yy,ddy
	real*8 x(imax),fx(ncol),xest
	save d,x
	x(iest)=xest
	if(iest.eq.1) then
	do 11 j=1,nv
		  yz(j)=yest(j)
	 	  d(j,1)=yest(j)
		  dy(j)=yest(j)
11	continue
	else
		m1=min(iest,nuse)
	do 12 k=1,m1-1
		fx(k+1)=x(iest-k)/xest
12	continue
	do 14 j=1,nv
		yy=yest(j)
		v=d(j,1)
		c=yy
		d(j,1)=yy
		do 13 k=2,m1
		  b1=fx(k)*v
		  b=b1-c
		  if(b.ne.0.d0) then
			b=(c-v)/b
			ddy=c*b
			c=b1*b
		  else
			ddy=v
		  endif
		  if(k.ne.m1) v=d(j,k)
		  d(j,k)=ddy
		  yy=yy+ddy
13	continue
		dy(j)=ddy
		yz(j)=yy
14	continue
	endif
	return
	end
