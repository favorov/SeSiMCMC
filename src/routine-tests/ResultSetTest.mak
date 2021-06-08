#this file can be started only from parent directory
name=ResultSetTest
#-pedantic -Wall -ffor-scope
md=.
td=./routine-tests
od=./obj

srcdirlist=$(md):$(td)

empty=
space=$(empty) $(empty)
includeflags = $(foreach dir,$(subst :,$(space),$(srcdirlist)),-I$(dir))
#this strange invocation is just preparing -I flag from srcdirlist.

include ~/include/ccvars

.PHONY: all clean

vpath %.c $(srcdirlist)
vpath %.cpp $(srcdirlist)

all: $(name)$(EXEEXT) 

OBJS=$(od)/$(name).o \
$(od)/SymbolsCounter.o \
$(od)/Sequences.o \
$(od)/ResultsSet.o \
$(od)/MCMC.o \
$(od)/Random.o \

$(od)/%.o: %.c
	$(CC) $(CCFLAGS) $< -o $@
	
$(od)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $< -o $@

$(name)$(EXEEXT): $(OBJS)
	$(CPP) -o $(name)$(EXEEXT) $(OBJS) $(LINKFLAGS) 

AtgcHeaders= $(md)/Atgc.hpp $(md)/Exception.hpp
SequencesHeaders=$(AtgcHeaders) $(md)/Sequences.hpp 
ResultsHeaders= $(md)/MarkovChainState.hpp $(md)/ResultsSet.hpp $(md)/MCMC.hpp\
	$(md)/Exception.hpp
SymbolsHeaders=$(md)/SymbolsCounter.hpp $(md)/MarkovChainState.hpp \
	$(md)/Exception.hpp
MCMCHeaders=$(md)/MCMC.hpp $(md)/MarkovChainState.hpp $(md)/Exception.hpp
RandomHeaders= $(md)/Random.h

$(od)/SymbolsCounter.o $(od)/MCMC.o $(od)/Sequences.o: $(AtgcHeaders)
$(od)/MCMC.o $(od)/$(name).o:  $(MCMCHeaders)
$(od)/SymbolsCounter.o $(od)/MCMC.o $(od)/$(name).o: $(SymbolsHeaders)
$(od)/Sequences.o $(od)/SymbolsCounter.o $(od)/MCMC.o $(od)/$(name).o: $(SequencesHeaders)
$(od)/Random.o $(od)/MCMC.o $(od)/$(name).o: $(RandomHeaders)
$(od)/$(name).o: $(ResultsHeaders)

$(od)/$(name).o : $(td)/$(name).cpp
$(od)/Sequences.o : $(md)/Sequences.cpp
$(od)/SymbolsCounter.o : $(md)/SymbolsCounter.cpp
$(od)/MCMC.o : $(md)/MCMC.cpp
$(od)/Random.o : $(md)/Random.c
$(od)/ResultsSet.o: $(md)/ResultsSet.cpp 



clean:
	rm -f $(OBJS)
	rm -r -f *~

