# Makefile for the solver 
# construted 19/04/2016. 
# updated --/--/----
# Typing 'make' will create the executable file.

#  -g    compile with debug mode, (-O2 optimaze more, -O3 optimaze most).
#  -I.   look in the current directory (.) for the include file example.h.
#  -Wall turns on most, but not all, compiler warnings
#  for single process define  CC = gcc, while for MPI define CC = mpicc.

CC      = mpicc
CFLAGS  = -I. -g
LIBS    = -lm 

TARGET = solver.exe
default: $(TARGET)

# To create the executable file we need the object files main.o, startm.o, etc. :
#
$(TARGET):   main.o start.o read.o set.o loop.o tmarch.o bound.o flux.o vflux.o thermochem.o save.o clean.o  
	$(CC) $(CFLAGS) -o $(TARGET) main.o start.o read.o set.o loop.o tmarch.o bound.o flux.o vflux.o thermochem.o save.o clean.o $(LIBS)

main.o:  main.c start.c read.c set.c loop.c clean.c
	$(CC) $(CFLAGS) -c main.c 

start.o:  start.c comm.h mpi.h
	$(CC) $(CFLAGS) -c start.c 

read.o:  read.c comm.h chemdata.h mpi.h
	$(CC) $(CFLAGS) -c read.c

thermochem.o:  thermochem.c clean.c comm.h chemdata.h mpi.h
	$(CC) $(CFLAGS) -c thermochem.c

set.o:  set.c thermochem.c clean.c comm.h chemdata.h mpi.h 
	$(CC) $(CFLAGS) -c set.c

loop.o:  loop.c thermochem.c flux.c clean.c save.c comm.h mpi.h
	$(CC) $(CFLAGS) -c loop.c

tmarch.o:  tmarch.c comm.h
	$(CC) $(CFLAGS) -c tmarch.c

bound.o:  bound.c comm.h
	$(CC) $(CFLAGS) -c bound.c

flux.o:  flux.c vflux.c bound.c thermochem.c comm.h
	$(CC) $(CFLAGS) -c flux.c

vflux.o:  vflux.c comm.h chemdata.h 
	$(CC) $(CFLAGS) -c vflux.c

save.o:  save.c clean.c comm.h chemdata.h 
	$(CC) $(CFLAGS) -c save.c

clean.o:  clean.c comm.h mpi.h 
	$(CC) $(CFLAGS) -c clean.c

# To start over from scratch, type 'make clean'.  This
# removes the executable file, as well as old .o object
# files and *~ backup files:
#
clean: 
	$(RM) solver.exe *.o *~
