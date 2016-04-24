export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig

# Set up basic variables
GFORTRAN = gfortran
CFLAGS = -c
LDFLAGS = 

# For fortran order does matter. 
# Cannot use $(wildcard *.f90) here
SOURCES = math.f90 \
			handlers.f90 \
			main.f90

OBJECTS = $(SOURCES:.f90=.o)
EXECUTABLE = prog

# Add flags for gtk
CFLAGS += `pkg-config --cflags gtk-2-fortran`
LDFLAGS += `pkg-config --libs gtk-2-fortran`

# Flags for LAPACK
LDFLAGS += -llapack -lblas

$(EXECUTABLE): $(OBJECTS)
	$(GFORTRAN) $(OBJECTS) -o $@ $(LDFLAGS)

$(OBJECTS) : $(SOURCES)
	$(GFORTRAN) $(CFLAGS) $(patsubst %.o, %.f90, $@)

clean:
	rm $(OBJECTS) $(EXECUTABLE) *.mod
