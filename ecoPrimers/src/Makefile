EXEC=ecoPrimers

PRIMER_SRC= ecoprimer.c
PRIMER_OBJ= $(patsubst %.c,%.o,$(PRIMER_SRC))


SRCS= $(PRIMER_SRC) 

LIB= -lecoprimer -lecoPCR -lthermo -lz -lm 

LIBFILE= libecoPCR/libecoPCR.a \
         libecoprimer/libecoprimer.a \
         libthermo/libthermo.a \
	 


include global.mk

all: $(EXEC)


########
#
# ecoPrimer compilation
#
########

# executable compilation and link

ecoPrimers: $(PRIMER_OBJ) $(LIBFILE) 
	$(CC) -g $(LDFLAGS) -O5 -m64 -o $@ $< $(LIBPATH) $(LIB)
	

########
#
# library compilation
#
########

libecoPCR/libecoPCR.a:
	$(MAKE) -C libecoPCR
	
libecoprimer/libecoprimer.a:
	$(MAKE) -C libecoprimer
	
libthermo/libthermo.a:
	$(MAKE) -C libthermo
	
########
#
# project management
#
########

clean:
	rm -f *.o	
	rm -f *.P
	rm -f $(EXEC)
	$(MAKE) -C libecoPCR clean
	$(MAKE) -C libecoprimer clean
	$(MAKE) -C libthermo clean
		

		
########
#
# clean for k2 to remove .o and .P files
#
########

k2clean:
	rm -f *.o	
	rm -f *.P
	rm -f libecoPCR/*.o
	rm -f libecoPCR/*.P
	rm -f libecoprimer/*.o
	rm -f libecoprimer/*.P
	rm -f libthermo/*.o
	rm -f libthermo/*.P
