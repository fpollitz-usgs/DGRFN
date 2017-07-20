# Makefile for DECAY, DECAY4, VTORDEP, VSPHDEP, STRAINA, INVERTV,
# STAT0A, STAT2A
#
#FC=ifort
#FFLAGS=-O2 -ftz -ip -ipo -axP -align all -extend_source
FC=gfortran
FFLAGS=-O2 -ffixed-form -ffixed-line-length-132 

all : wave0 wave1

wave0SAVE  : wave0SAVE.f msubs.f bsstep.f rzextr.f mmid.f interp.f
	${FC} -o $@ ${FFLAGS} wave0SAVE.f msubs.f bsstep.f rzextr.f mmid.f interp.f 

wave0  : wave0.f msubs.f bsstep.f rzextr.f mmid.f interp.f
	${FC} -o $@ ${FFLAGS} wave0.f msubs.f bsstep.f rzextr.f mmid.f interp.f 
#	gfortran -ffixed-form -ffixed-line-length-132 -o wave0 wave0.f msubs.f bsstep.f rzextr.f mmid.f interp.f

wave0-ndep  : wave0-ndep.f msubs.f bsstep.f rzextr.f mmid.f interp.f
	${FC} -o $@ ${FFLAGS} wave0-ndep.f msubs.f bsstep.f rzextr.f mmid.f interp.f 

wave0-2dep  : wave0-2dep.f msubs.f bsstep.f rzextr.f mmid.f interp.f
	${FC} -o $@ ${FFLAGS} wave0-2dep.f msubs.f bsstep.f rzextr.f mmid.f interp.f 

wave0-layer  : wave0-layer.f msubs.f bsstep.f rzextr.f mmid.f interp.f
	${FC} -o $@ ${FFLAGS} wave0-layer.f msubs.f bsstep.f rzextr.f mmid.f interp.f 

wave1  : wave1.f lgndrl.f invft.f 
#	${FC} -o $@ ${FFLAGS} wave1.f lgndrl.f invft.f 
	gfortran -ffixed-form -ffixed-line-length-132 -o wave1 wave1.f lgndrl.f invft.f

wave1-displ  : wave1-displ.f lgndrl.f splneq.f invft.f 
	${FC} -o $@ ${FFLAGS} wave1-displ.f lgndrl.f splneq.f invft.f 

wave1-ndep  : wave1-ndep.f lgndrl.f splneq.f invft.f 
	${FC} -o $@ ${FFLAGS} wave1-ndep.f lgndrl.f splneq.f invft.f 

wave1-layer  : wave1-layer.f lgndrl.f invft.f 
	${FC} -o $@ ${FFLAGS} wave1-layer.f lgndrl.f invft.f 

wave1mte  : wave1mte.f lgndrl.f invft.f 
	${FC} -o $@ ${FFLAGS} wave1mte.f lgndrl.f invft.f 

wave1mtr  : wave1mtr.f lgndrl.f hfn.f invft.f 
	${FC} -o $@ ${FFLAGS} wave1mtr.f lgndrl.f hfn.f invft.f 

hglook  : hglook.f hfn.f 
	${FC} -o $@ ${FFLAGS} hglook.f hfn.f

hfs  : hfs.f invft.f 
	${FC} -o $@ ${FFLAGS} hfs.f invft.f

