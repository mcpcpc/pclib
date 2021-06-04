.POSIX:
ALL_WARNING = -Wall -Wextra -pedantic
ALL_LDFLAGS = $(LDFLAGS)
ALL_CFLAGS = $(CPPFLAGS) $(CFLAGS) -std=c99 $(ALL_WARNING)
PREFIX = /usr/local
LDLIBS = -lm
BINDIR = $(PREFIX)/bin

all: mat
install: all
	mkdir -p $(DESTDIR)$(BINDIR)
	cp -f pca $(DESTDIR)$(BINDIR)
	chmod 755 $(DESTDIR)$(BINDIR)/mat
mat: mat.o
	$(CC) $(ALL_LDFLAGS) -o mat mat.o $(LDLIBS)
mat.o: mat.c
	$(CC) $(ALL_CFLAGS) -c mat.c
clean:
	rm -f mat *.o
uninstall:
	rm -f $(DESTDIR)$(BINDIR)/mat
.PHONY: all install uninstall clean