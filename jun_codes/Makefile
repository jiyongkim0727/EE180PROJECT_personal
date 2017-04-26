CC = gcc
CFLAGS = -lmraa -lfann

default: extract

extract: main.c
	$(CC) $(CFLAGS)  main.c -o $@

clean:
	rm -f extract
