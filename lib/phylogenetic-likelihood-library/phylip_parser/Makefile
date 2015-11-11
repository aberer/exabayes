CC           = gcc
CFLAGS       = -g -Wall 

all: parser parser-mf

parser: lexer.c phylip.c ssort.c main.c
	@echo "===> COMPILE $@"
	$(CC) $(CFLAGS) $+ -o $@

parser-mf: lexer.c phylip.c ssort.c main.c
	@echo "===> COMPILE $@"
	$(CC) $+ -o $@ -fmudflap -lmudflap

clean:
	-rm parser parser-mf

.PHONY: all clean
