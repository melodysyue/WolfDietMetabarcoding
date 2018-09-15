
SOURCES = goodtaxon.c \
          readdnadb.c \
          smothsort.c \
          sortword.c \
          hashsequence.c \
          strictprimers.c \
          aproxpattern.c \
          merge.c \
          queue.c \
          libstki.c \
          sortmatch.c \
          pairtree.c \
          pairs.c \
          taxstats.c \
          apat_search.c \
		  filtering.c \
		  PrimerSets.c \
		  ahocorasick.c

SRCS=$(SOURCES)
         
OBJECTS= $(patsubst %.c,%.o,$(SOURCES))

LIBFILE= libecoprimer.a
RANLIB= ranlib


include ../global.mk


all: $(LIBFILE)

clean:
	rm -rf $(OBJECTS) $(LIBFILE)

$(LIBFILE): $(OBJECTS)
	ar -cr $@ $?
	$(RANLIB) $@
