#---------------------------------------------------------------------#
# Makefile for Monte Carlo B-W code                                   #
#                                                                     #
# C. D. Woodgate, Warwick                                        2023 #
#---------------------------------------------------------------------#

# Specify a particular compiler with "make compiler=intel", etc.
# Alternatively, these options can be specified as environment variables
# e.g. "export compiler=gfortran" can be added to $HOME/.bashrc

SYSTEM=$(shell uname -s)

# Compiler specific flags
# Note: you MUST specify a compiler option. None are specified by default.

ifeq ($(strip $(compiler)),)
  MAKECMDGOALS = error
error:
	@echo 'You need to set a value for the compiler variable'
	@echo ' eg. "make compiler=gfortran"'
	@echo 'Alternatively, you can add "export compiler=gfortran" to $$HOME/.bashrc'
endif

# Intel
ifeq ($(strip $(compiler)),intel)
  FC = mpiifort
  FFLAGS = -O3 
  FFLAGS = -O0 -warn all -check all 
  FFLAGS += -module $(OBJDIR)
  CC = icc
  CFLAGS = -O3
endif

# gfortran
ifeq ($(strip $(compiler)),gfortran)
  FC = mpifort
#  FFLAGS = -g -Wall -Wextra -fcheck=bounds
#  FFLAGS = -O3 -g -Wall -Wextra -fcheck=all
  FFLAGS = -O3 -Wall -Wextra -fimplicit-none
  FFLAGS += -I/usr/local/include -I$(OBJDIR) -J$(OBJDIR)
  LDFLAGS=-lgcc -lopenblas
  CC=gcc -I$(INCDIR)
  CFLAGS=-O3
endif

SRCDIR=src
BINDIR=bin
DATDIR=data
OBJDIR=obj
INCDIR=include

ifeq ($(SYSTEM),Darwin)
         FFLAGS += $(shell nf-config --fflags)
         LDFLAGS += $(shell nf-config --flibs) \
                    -lnetcdf -lnetcdff
else
         FFLAGS +=$(shell nf-config --fflags)
         LDFLAGS += $(shell nf-config --flibs) 
endif

# Command to use for linking and executable
LD=$(FC)
EXE=brawl.run

MODFILES=mt19937ar.c kinds.f90 shared_data.f90 io.f90 comms.f90 write_netcdf.f90 \
         write_xyz.f90 metropolis_output.f90 command_line.f90 c_functions.f90 \
         display.f90 bw_hamiltonian.f90 analytics.f90 random_site.f90 \
         metropolis.f90 nested_sampling.f90 tmmc.f90 wang-landau.f90 \
         energy_spectrum.f90 config_output.f90 initialise.f90 constants.f90 \
         derived_types.o

SRCFILES=$(MODFILES) main.f90

EXFILES=$(MODFILES) howto_examples.f90 example.f90

OBJFILES:=$(SRCFILES:.f90=.o)
OBJFILES:=$(OBJFILES:.c=.o)

EXOBJFILES:=$(EXFILES:.f90=.o)
EXOBJFILES:=$(EXOBJFILES:.c=.o)

VPATH = $(SRCDIR):$(SRCDIR)/core:$(OBJDIR):$(INCDIR)

alloy: $(OBJFILES)
	$(FC) $(FFLAGS) -o $(EXE) $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS)

example: $(EXOBJFILES)
	$(FC) $(FFLAGS) -o $(EXE) $(addprefix $(OBJDIR)/,$(EXOBJFILES)) $(LDFLAGS)

# Purge build files and executable
clean :
	@rm -rf $(OBJDIR) $(BINDIR) $(EXE)

# Rules for building object files
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $(OBJDIR)/$@ $<

# Rules for building object files
%.o: %.c
	$(CC) $(CFLAGS) -c -o $(OBJDIR)/$@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(OBJFILES): | $(OBJDIR)

################
# Dependencies #
################
shared_data.o: kinds.o
derived_types.o: kinds.o
comms.o: kinds.o shared_data.o
io.o: kinds.o shared_data.o command_line.o display.o comms.o
write_netcdf.o: kinds.o shared_data.o
write_xyz.o: shared_data.o kinds.o analytics.o
metropolis_output.o: shared_data.o kinds.o analytics.o
command_line.o: kinds.o
c_functions.o: mt19937ar.o
display.o: kinds.o shared_data.o constants.o derived_types.o
bw_hamiltonian.o: kinds.o shared_data.o c_functions.o io.o
analytics.o: shared_data.o kinds.o display.o io.o
random_site.o: shared_data.o kinds.o c_functions.o analytics.o
metropolis.o: kinds.o shared_data.o c_functions.o bw_hamiltonian.o random_site.o analytics.o initialise.o
nested_sampling.o: kinds.o shared_data.o c_functions.o bw_hamiltonian.o random_site.o analytics.o initialise.o metropolis.o
initialise.o: kinds.o shared_data.o c_functions.o bw_hamiltonian.o random_site.o comms.o
main.o: initialise.o shared_data.o kinds.o c_functions.o write_netcdf.o\
	write_xyz.o metropolis_output.o command_line.o display.o metropolis.o
