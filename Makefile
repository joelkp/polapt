.POSIX:
CC=cc
CFLAGS=-std=c99 -W -Wall -O3
LFLAGS=-lm
BIN=polapt
OBJ=polapt.o

all: $(BIN)
clean:
	rm -f $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) $(LFLAGS) -o $(BIN)

polapt.o: polapt.c
	$(CC) -c $(CFLAGS) polapt.c
