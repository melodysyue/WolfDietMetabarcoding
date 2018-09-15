MACHINE=MAC_OS_X
LIBPATH= -LlibecoPCR -Llibecoprimer -Llibthermo
MAKEDEPEND = gcc -D$(MACHINE) -M $(CPPFLAGS) -o $*.d $<

CC=gcc
CFLAGS= -W -Wall -m64 -g
#CFLAGS= -W -Wall -O5 -m64 -g
#CFLAGS= -W -Wall -O0 -m64  -g
#CFLAGS= -W -Wall -O5 -fast -g

default: all

%.o: %.c
	$(CC) -D$(MACHINE) $(CFLAGS) -c -o $@ $<

%.P : %.c
	$(MAKEDEPEND)
	@sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' < $*.d > $@; \
	rm -f $*.d; [ -s $@ ] || rm -f $@

include $(SRCS:.c=.P)
