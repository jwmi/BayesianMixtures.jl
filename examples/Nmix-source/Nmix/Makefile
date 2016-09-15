FC = f77
CC = cc
FFLAGS = -O2

Nmix: Nmix.o pnorm.o algama.o rgamma.o gauss4.o sd.o
	f77 -o Nmix Nmix.o pnorm.o algama.o rgamma.o gauss4.o sd.o
sd.o: sd.c
	cc -c -o sd.o -DRETS -DSUNF sd.c



