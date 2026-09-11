#--------------------------------------------------------------
#                 Makefile for quantum module
#--------------------------------------------------------------

# List module object filenames
quantum_objects=\
bath.o \
bath_fft.o \
bath_log.o \
bath_ou.o \
data.o \
export.o \
initialize.o \
interface.o \
llg.o \
llg_mpi.o \
llg_serial.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/quantum/,$(quantum_objects))
