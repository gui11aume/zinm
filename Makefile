OBJECTS= zinm.o

CFLAGS= -std=c99 -g -Wall -O3
LDLIBS= -lm
CC= gcc

all: do

do: $(OBJECTS) main.c
	$(CC) $(CFLAGS) $(INCLUDES) main.c $(OBJECTS) $(LDLIBS) -o do

clean:
	rm -f $(OBJECTS) do

%.o: %.c %.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

