OBJECTS= zinm.o

CFLAGS= -std=gnu99 -g -Wall -O3
LDLIBS= -lm
CC= gcc

all: $(OBJECTS)

%.o: %.c %.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS)

