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

REQUIRED_INCLUDES_PATH:=-I/home/fercer/Apps/gmffilter/include -I/usr/include
REQUIRED_LIBS_PATH:=-L/usr/lib/x86_64-linux-gnu
REQUIRED_LIBS_LINK:=-lfftw3

LIBPREFIX:=lib
LIBNAME:=gmf
LIBSUFFIX:=.so
LIBVERSION:=.1.0

make: $(ODIR) $(LIBDIR)
	gcc -c -fpic -Wl,-soname,$(LIBPREFIX)$(LIBNAME) -DNDEBUG $(REQUIRED_INCLUDES_PATH) $(SRCS) -o $(OBJS) -std=c99 -O2
	gcc -shared $(REQUIRED_LIBS_PATH) -o $(LIBDIR)/$(LIBPREFIX)$(LIBNAME)$(LIBSUFFIX)$(LIBVERSION) $(OBJS) $(REQUIRED_LIBS_LINK)

$(ODIR):
	mkdir $(ODIR)

$(LIBDIR):
	mkdir $(LIBDIR)

install:
	cp $(LIBDIR)/$(LIBPREFIX)$(LIBNAME)$(LIBSUFFIX)$(LIBVERSION) /usr/local/lib

debug: $(ODIR) $(LIBDIR)
	gcc -c -fpic -Wl,-soname,$(LIBPREFIX)$(LIBNAME) $(REQUIRED_INCLUDES_PATH) $(SRCS) -o $(OBJS) -std=c99 -g
	gcc -shared $(REQUIRED_LIBS_PATH) -o $(LIBDIR)/$(LIBPREFIX)$(LIBNAME)$(LIBSUFFIX)$(LIBVERSION) $(OBJS) $(REQUIRED_LIBS_LINK)

python: $(ODIR) $(LIBDIR)
	gcc -c -fpic -Wl,-soname,$(LIBPREFIX)$(LIBNAME) -DNDEBUG -DBUILDING_PYTHON_MODULE $(REQUIRED_INCLUDES_PATH) -I/usr/include/python2.7 $(SRCS) -o $(OBJS) -std=c99 -O2
	gcc -shared $(REQUIRED_LIBS_PATH) -o $(LIBDIR)/$(LIBNAME).pyd $(OBJS) $(REQUIRED_LIBS_LINK) -lpython2.7
	
install_python:
	cp $(LIBDIR)/$(LIBNAME).pyd /home/fercer/anaconda2/envs/sigproc_env/lib/python2.7/site-packages
	
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