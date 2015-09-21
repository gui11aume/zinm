RLIB= Rzinm.so

SRC_DIR= ../..
INC_DIR= ../..

INCLUDES= $(addprefix -I, $(INC_DIR))
OBJECTS= zinm.o
SOURCES= Rzinm.c

CC= gcc
CFLAGS= -std=gnu99 -fPIC -Wall -g -O0

all: $(RLIB)

%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

Rzinm.so: $(OBJECTS) $(SOURCES)
	R CMD SHLIB $(INCLUDES) $(SOURCES) $(OBJECTS)

clean:
	rm -f $(OBJECTS) $(RLIB)
