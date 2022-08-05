CC      = gcc
CFLAGS  = -g
RM      = rm -f

default: all

all: decode

decode: decode.c
	$(CC) $(CFLAGS) -o decode decode.c
	
del:
	$(RM) decode