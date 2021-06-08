#****************************************************************************#
#SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001
#$Id: StupidTest.mak 1014 2009-03-01 16:50:36Z favorov $
#****************************************************************************#

name=StupidTest
exename=StupidTest
md=.
td=./routine-tests
od=./obj

srcdirlist=$(md):$(td)

empty=
space=$(empty) $(empty)
includeflags = $(foreach dir,$(subst :,$(space),$(srcdirlist)),-I$(dir))
#this strange invocation is just preparing -I flag from srcdirlist.

OPTIMISE=YES

include ~/include/ccvars

.PHONY: all clean

vpath %.c $(srcdirlist)
vpath %.cpp $(srcdirlist)

all: $(exename)$(EXEEXT) 

OBJS=$(od)/$(name).o \
$(od)/SymbolsCounter.o \
$(od)/Random.o \
$(od)/Sequences.o 


$(od)/%.o: %.c
	$(CC) $(CCFLAGS) $< -o $@
	
$(od)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $< -o $@

$(exename)$(EXEEXT): $(OBJS)
	$(CPP) -o $(exename)$(EXEEXT) $(OBJS) $(LINKFLAGS)
	chmod 755 $(exename)$(EXEEXT) 

AtgcHeaders= $(md)/Atgc.hpp $(md)/Exception.hpp
SequencesHeaders=$(AtgcHeaders) $(md)/Sequences.hpp 
SymbolsHeaders=$(md)/SymbolsCounter.hpp $(md)/MarkovChainState.hpp \
	$(md)/Exception.hpp
KullbakHeaders=$(md)/SymbolsCounter.hpp $(md)/KullbakCounter.hpp \
	$(md)/Exception.hpp
MCMCStHeaders=$(md)/MarkovChainState.hpp $(md)/Exception.hpp
RandomHeaders= $(md)/Random.h

$(od)/SymbolsCounter.o $(od)/Sequences.o: $(AtgcHeaders)
$(od)/SymbolsCounter.o $(od)/$(name).o: $(SymbolsHeaders)
$(od)/Random.o $(od)/$(name).o: $(RandomHeaders)
$(od)/$(name).o: $(KullbakHeaders)
$(od)/Sequences.o $(od)/SymbolsCounter.o $(od)/$(name).o: $(SequencesHeaders)


$(od)/$(name).o : $(td)/$(name).cpp
$(od)/Sequences.o : $(md)/Sequences.cpp
$(od)/SymbolsCounter.o : $(md)/SymbolsCounter.cpp
$(od)/Random.o : $(md)/Random.c

clean:
	rm -f $(OBJS)
	rm -r -f *~
