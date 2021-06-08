#this file can be started only from parent directory
name=SymbolsCounterTest
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
$(od)/Random.o \
$(od)/Sequences.o


$(od)/%.o: %.c
	$(CC) $(CCFLAGS) $< -o $@
	
$(od)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $< -o $@

$(name)$(EXEEXT): $(OBJS)
	$(CPP) -o $(name)$(EXEEXT) $(OBJS) $(LINKFLAGS) 

AtgcHeaders= $(md)/Atgc.hpp $(md)/Exception.hpp

SequencesHeaders=$(AtgcHeaders) $(md)/Sequences.hpp 

SymbolsHeaders=$(md)/SymbolsCounter.hpp $(md)/MarkovChainState.hpp \
  $(md)/Exception.hpp

RandomHeaders= $(md)/Random.h

$(od)/$(name).o $(od)/SymbolsCounter.o : $(SequencesHeaders) $(SymbolsHeaders)
$(od)/Sequences.o : $(SequencesHeaders)
$(od)/Sequences.o $(od)/Random.o : $(RandomHeaders)

$(od)/$(name).o : $(td)/$(name).cpp
$(od)/SymbolsCounter.o : $(md)/SymbolsCounter.cpp
$(od)/Sequences.o : $(md)/Sequences.cpp
$(od)/Random.o : $(md)/Random.c

clean:
	rm -f $(OBJS)
	rm -r -f *~
