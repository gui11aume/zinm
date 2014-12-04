OBJECTS= zinm.o gen.o pso.o

CFLAGS= -std=gnu99 -g -Wall -O3
LDLIBS= -lpthread -lm
CC= gcc

all: do

do: $(OBJECTS) main.c 
	$(CC) $(CFLAGS) $(INCLUDES) main.c $(OBJECTS) $(LDLIBS) -o do

clean:
	rm -f $(OBJECTS) do

%.o: %.c %.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

