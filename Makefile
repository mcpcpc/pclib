.POSIX:
ALL_WARNING = -Wall -Wextra -pedantic
ALL_LDFLAGS = $(LDFLAGS)
ALL_CFLAGS = $(CPPFLAGS) $(CFLAGS) -std=c99 $(ALL_WARNING)
PREFIX = /usr/local
LDLIBS = -lm
BINDIR = $(PREFIX)/bin

all: pca
install: all
	mkdir -p $(DESTDIR)$(BINDIR)
	cp -f pca $(DESTDIR)$(BINDIR)
	chmod 755 $(DESTDIR)$(BINDIR)/pca
pca: pca.o
	$(CC) $(ALL_LDFLAGS) -o pca pca.o $(LDLIBS)
pca.o: pca.c
	$(CC) $(ALL_CFLAGS) -c pca.c
clean:
	rm -f pca *.o
uninstall:
	rm -f $(DESTDIR)$(BINDIR)/pca
.PHONY: all install uninstall clean
