CC = af95 
CFLAGS = -Ofast -m64 
COPTIONS = 
FILE = beam
OBJS = constant.o params.o service.o nonlin.o solver.o reader.o writer.o initbeam.o plasma.o fmain.o

prog: $(OBJS)
	@echo "linking..."
	$(CC) $(CFLAGS) -o$(FILE) $(OBJS) 
	make clean


$(OBJS): 
	$(CC) -I $(CFLAGS) $(OPTIONS) -c $*.f90

clean:
	@echo "cleaning up..."
	rm -f *.o
	rm -f *~
	rm -f *.mod



