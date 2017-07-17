# Makefile for DECAY, DECAY4, DECAY4M, VSPHDEP,
# VSPHM, VTORDEP, and STRAINA
# Options below are for Absoft compiler on MAC OSX10.3
#FC=f77
#FFLAGS=-O -q -W -s
# Options below are for Intel compiler on MAC OSX10.6
FC=ifort
FFLAGS=-O2 -ftz -ip -ipo -axP -mp -align all -extend_source

all : multw2 stat0A stat2gA addw2 getg

multw2 : multw2.f
	${FC} -o $@ ${FFLAGS} multw2.f

stat0A : stat0A.f
	${FC} -o $@ ${FFLAGS} stat0A.f

stat2gA : stat2gA.f lgndrl.f splneq.f
	${FC} -o $@ ${FFLAGS} stat2gA.f lgndrl.f splneq.f

addw2 : addw2.f
	${FC} -o $@ ${FFLAGS} addw2.f

getg : getg.f
	${FC} -o $@ ${FFLAGS} getg.f

getz : getz.f
	${FC} -o $@ ${FFLAGS} getz.f
disph : disph.f
	${FC} -o $@ ${FFLAGS} disph.f

