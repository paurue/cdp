PROG = cdp

OBJECTS = sptk.o 
MAINOBJ = $(PROG).o
MAINSRC = $(PROG).c

CC = gcc 
CFLAGS = -Wall -O3 
LIBS = -lm  -lfftw3 -lgsl 

all: $(MAINOBJ) $(OBJECTS)
	$(CC) $(CFLAGS) $(OPTS) $(MAINOBJ) $(OBJECTS) -o $(PROG) $(LIBS) 

$(MAINOBJ):	$(MAINSRC) $(OBJECTS)
	$(CC) $(CFLAGS) $(OPTS) -c $(MAINSRC)

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o
