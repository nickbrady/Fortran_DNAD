## Copyright (c) 2019 alberto Otero de la Roza <aoterodelaroza@gmail.com>
## This file is frere software; distributed under GNU/GPL version 3.

## FC=compiler
## FCSYNTAX=syntax-only compilation flag
## FCMODDIR=flag to set read/write .mod and .smod directory
## FCMODREADDIR=flag to set write-only .mod and .smod directory
## MODDIR=.mod and .smod directory (leave blank for root)
FC:=gfortran
FCSYNTAX:=-fsyntax-only
FCMODDIR:=-J
MODDIR:=.mod
# FC:=ifort
# FCSYNTAX:=-syntax-only
# FCMODDIR:=-module
# FCMODREADDIR:=-I
# MODDIR:=.mod

#### user input ends here ####

## some tricks for text manipulation
null:=
space:=$(null) $(null)
$(space):=$(space)
define \n


endef


## no implicit rules
.SUFFIXES:

## auxiliary programs
AWK:=gawk
SED:=sed
RM:=rm -f
MKDIR:=mkdir -p
TEST:=test

## known fortran extensions
FORTEXT:=f fpp for ftn f90 f95 f03 f08

## locate the source files (-iregex case-insensitive)
# FILEPATH:=/Users/nicholasbrady/Downloads/example-final-1/Test/
SOURCES:=$(shell gfind . -regextype posix-awk -iregex '.*\.($(subst $( ),|,$(FORTEXT)))$$')
$(info $$SOURCES is [${SOURCES}])

COREMODULES:=/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/CoreFortran_Subroutines_Modules
SOURCES +=$(shell gfind $(COREMODULES) -maxdepth 1 -regextype posix-awk -iregex '.*\.($(subst $( ),|,$(FORTEXT)))$$')
$(info $$SOURCES is [${SOURCES}])

# SOURCES:=$(foreach ext,$(FORTEXT), $(wildcard *.$(ext)) $(wildcard **/*.$(ext)))
# $(info $$SOURCES is [${SOURCES}])

FORTEXT:=$(suffix $(SOURCES))
FORTEXT:=$(sort $(FORTEXT))						# produces the unique extensions
FORTEXT:=$(subst ., ,$(FORTEXT))			# remove the .'s'
$(info $$FORTEXT is [${FORTEXT}])

## compilation and syntax-compilation commands
COMPILE.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c
MAKEMOD.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) $(FCSYNTAX) -c

## create the mod and smod directory; define slashed version of MODDIR
ifneq ($(MODDIR),)
  $(shell $(TEST) -d $(MODDIR) || $(MKDIR) -p $(MODDIR))
  MODDIRSLSH:=$(MODDIR)/
else
  MODDIRSLSH:=./
endif

## create the temporary mod and smod directory
ifneq ($(FCMODREADDIR),)
  MODDIRTMP:=.tmp$(MODDIR)
  $(shell $(TEST) -d $(MODDIRTMP) || $(MKDIR) -p $(MODDIRTMP))
  MAKEMOD.f08+= $(FCMODDIR) $(MODDIR)
  COMPILE.f08+= $(FCMODREADDIR) $(MODDIR) $(FCMODDIR) $(MODDIRTMP)
else
  MAKEMOD.f08+= $(FCMODDIR) $(MODDIR)
  COMPILE.f08+= $(FCMODDIR) $(MODDIR)
endif

## define the anchors and the objects variables
# $(call source-to-extension,source-file-list,new-extension)
define source-to-extension
  $(strip \
    $(foreach ext,$(FORTEXT),\
      $(subst .$(ext),.$2,$(filter %.$(ext),$1))))
endef
OBJECTS:=$(call source-to-extension,$(SOURCES),o)
ANCHORS:=$(call source-to-extension,$(SOURCES),anc)

## define the anchors and the objects variables
BASESOURCE:=$(basename $(SOURCES))
OBJECTS:=$(addsuffix .o, $(BASESOURCE))
ANCHORS:=$(addsuffix .anc, $(BASESOURCE))

## default target, main and clean targets
all: main

main: $(OBJECTS)
	$(FC) -o $@ $+

.PHONY: clean
clean:
	-$(RM) *.mod *.smod $(OBJECTS) $(ANCHORS) main
	-$(TEST) -d $(MODDIR) && $(RM) -r $(MODDIR)
	-$(TEST) -d $(MODDIRTMP) && $(RM) -r $(MODDIRTMP)

# syntax-only compilation rule: all anchor files depend on their source
# $(call modsource-pattern-rule,extension)
define modsource-pattern-rule
%.anc: %.$1
	$$(MAKEMOD.f08) $$<
	@touch $$@
endef
$(foreach ext,$(FORTEXT),$(eval $(call modsource-pattern-rule,$(ext))))

# %.anc: $1 -- syntax-only compilation rule: all anchor files depend on their source
# $(call modsource-pattern-rule,source_file)
# define modsource-pattern-rule
# $(addsuffix .anc, $(basename $1)): $1
# 	$$(MAKEMOD.f08) $$<
# 	@touch $$@
# endef
# $(foreach src,$(SOURCES),$(eval $(call modsource-pattern-rule,$(src))))

# compilation rule: objects depend on their anchor file
%.o: %.anc
	$(COMPILE.f08) $(OUTPUT_OPTION) $(wildcard $(addprefix $*.,$(FORTEXT)))
ifdef MODDIRTMP
	-@$(RM) $(MODDIRTMP)/*.mod $(MODDIRTMP)/*.smod
endif
	@touch $@

# define compilation-rule
# $(addsuffix .o, $(basename $1)): $(addsuffix .anc, $(basename $1))
# 	$$(COMPILE.f08) $$(OUTPUT_OPTION) $$1
# 	@echo $$(COMPILE.f08) $$(OUTPUT_OPTION) $($1)
# 	@touch $$@
# endef
# $(foreach src,$(SOURCES),$(eval $(call compilation-rule,$(src))))
#
# ifdef MODDIRTMP
# 	-@$(RM) $(MODDIRTMP)/*.mod $(MODDIRTMP)/*.smod
# endif

## automatically generate the dependency rules
$(eval $(subst $( ),$(\n),$(shell $(AWK) --traditional -f $(COREMODULES)/makedepf08.awk $(SOURCES) | sort | uniq | $(SED) -e 's!^.mod/!$(MODDIRSLSH)!' -e 's!:.mod/!:$(MODDIRSLSH)!')))
