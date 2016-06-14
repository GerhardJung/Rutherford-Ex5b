#
CFLAGS=-O2
LFLAGS=-O2

all:Rutherford

Rutherford : rutherford.c
	g++ rutherford.c -o Rutherford $(CFLAGS)
clean:
	rm Rutherford
