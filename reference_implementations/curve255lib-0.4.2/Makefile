all: curve25519.o curve25519 curve25519test
curve25519: curve25519.o curve25519cmd.o base32.o
curve25519test: curve25519.o curve25519test.o base32.o

CFLAGS=-O2 -Wall
LDLIBS=-lgmp

clean:
	rm -f *.o curve25519test curve25519 curve25519cmd
