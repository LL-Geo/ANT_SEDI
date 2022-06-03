###################################
# -*- Makefile -*-
# Huian LI, created on June 13, 2008
####################################

# FC	= mpif90
# LD	= mpif90

FC      = af90
LD      = af77

FCFLAGS	= -s -O2 -c
LDFLAGS =

RM	= rm

rift2d_fail3: rift2d_fail3.o
	$(LD) $(LDFLAGS) -o $@ $<

rift2d_fail3.o: rift2d_fail3.f
	$(FC) $(FCFLAGS) -c -o $@ $<


clean: 
	$(RM) rift2d_fail3.o rift2d_fail3
