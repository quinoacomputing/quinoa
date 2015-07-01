
CC=gcc

CFLAGS=-O2 -Wall

SRC =	*.c

OBJ =	$(SRC:.c=.o)

OUT =	../lib/librngsse.a

$(OUT): $(OBJ)
	ar rcs $(OUT) $(OBJ)

$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c $(SRC)

clean: 
	rm $(OBJ)
