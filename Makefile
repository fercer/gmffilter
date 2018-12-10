LIBVERSION:=1.0

SRCDIR:=./include
LIBDIR:=./lib
ODIR:=./obj

EXAMPLESDIR:=./examples
EXAMPLESSRCDIR:=$(EXAMPLESDIR)/src
EXAMPLESODIR:=$(EXAMPLESDIR)/obj
EXAMPLESBINDIR:=$(EXAMPLESDIR)/bin

SRCS:=$(wildcard $(SRCDIR)/*.c)
HDRS:=$(wildcard $(SRCDIR)/*.h)
OBJS:=$(patsubst $(SRCDIR)/%,$(ODIR)/%, $(patsubst %.c,%.o, $(SRCS)))

LIBS:=libgmf.so.1

make: $(OBJS) $(LIBDIR)
	gcc -shared -L/usr/lib/x86_64-linux-gnu -lfftw3 -o $(LIBDIR)/$(LIBS).$(LIBVERSION) $(OBJS) -lc

.PHONY: $(OBJS)
$(OBJS): $(ODIR)
	gcc -c -fpic -Wl,-soname,$LIBS$ -I/usr/include/ $(SRCS) -o $(OBJS) -std=c99

$(ODIR):
	mkdir $(ODIR)

$(LIBDIR):
	mkdir $(LIBDIR)

install:
	cp $(LIBDIR)/$(LIBS).$(LIBVERSION) /usr/local/lib
	
	
example_gmf: $(EXAMPLESODIR)/test_gmf.o $(EXAMPLESBINDIR)
	gcc -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu $(EXAMPLESODIR)/test_gmf.o -o $(EXAMPLESBINDIR)/test_gmf.bin -lgmf -lfftw3 -lm
	
.PHONY: $(EXAMPLESODIR)/test_gmf.o
$(EXAMPLESODIR)/test_gmf.o: $(EXAMPLESODIR)
	@echo "Compiling the examples into " $(EXAMPLESODIR)
	gcc -c -I$(SRCDIR) -I/usr/include $(EXAMPLESSRCDIR)/test_gmf.c -o $(EXAMPLESODIR)/test_gmf.o -std=c99

$(EXAMPLESODIR):
	@echo "Creating directory to store the compiled objects of the examples ..."
	@echo $(EXAMPLESODIR)
	mkdir $(EXAMPLESODIR)
	
$(EXAMPLESBINDIR):
	mkdir $(EXAMPLESBINDIR)

$(EXAMPLESDIR):
	mkdir $(EXAMPLESDIR)