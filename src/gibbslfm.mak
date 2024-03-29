#****************************************************************************#
#SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
#$Id$
#****************************************************************************#

name=gibbslfm
exename=SeSiMCMC
md=.
td=.
od=./obj

srcdirlist=$(md):$(td)

empty=
space=$(empty) $(empty)
includeflags = $(foreach dir,$(subst :,$(space),$(srcdirlist)),$(INCLUDEKEY)$(dir)) $(INCLUDECLOSETERM)
#this strange invocation is just preparing -I flag from srcdirlist.

OPTIMISE=YES

include ccvars

.PHONY: all objs clean tgz 

vpath %.c $(srcdirlist)
vpath %.cpp $(srcdirlist)

all: $(exename)$(EXEEXT) 

OBJS=$(od)/$(name).o \
$(od)/SymbolsCounter.o \
$(od)/Sequences.o \
$(od)/MCMC.o \
$(od)/ResultsSet.o \
$(od)/Random.o \
$(od)/confread.o


objs:$(OBJS)

$(od)/%.o: %.c
	$(CC) $(CCFLAGS) $< -o $@
	
$(od)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $< -o $@

$(exename)$(EXEEXT): $(OBJS)
	$(CPP) -o $(exename)$(EXEEXT) $(OBJS) $(LINKFLAGS)
	chmod 755 $(exename)$(EXEEXT) 

#generated by g++ -MM *.c*
$(od)/MCMC.o: MCMC.cpp Exception.hpp Sequences.hpp Atgc.hpp Random.h \
  MarkovChainState.hpp SymbolsCounter.hpp KullbakCounter.hpp MCMC.hpp \
  Logger.hpp Diagnostics.hpp RandomMappings.hpp
$(od)/Random.o: Random.c Random.h
$(od)/ResultsSet.o: ResultsSet.cpp ResultsSet.hpp Sequences.hpp Exception.hpp \
  Atgc.hpp Random.h SymbolsCounter.hpp MarkovChainState.hpp \
  KullbakCounter.hpp MCMC.hpp Logger.hpp Diagnostics.hpp
$(od)/Sequences.o: Sequences.cpp Sequences.hpp Exception.hpp Atgc.hpp Random.h \
  ResultsSet.hpp SymbolsCounter.hpp MarkovChainState.hpp \
  KullbakCounter.hpp MCMC.hpp Logger.hpp Diagnostics.hpp
$(od)/SymbolsCounter.o: SymbolsCounter.cpp Atgc.hpp Exception.hpp Random.h \
  SymbolsCounter.hpp Sequences.hpp MarkovChainState.hpp
$(od)/confread.o: confread.c confread.h
$(od)/gibbslfm.o: gibbslfm.cpp Sequences.hpp Exception.hpp Atgc.hpp Random.h \
  MarkovChainState.hpp ResultsSet.hpp SymbolsCounter.hpp \
  KullbakCounter.hpp MCMC.hpp Logger.hpp Diagnostics.hpp confread.h

clean:
	rm -f $(OBJS)
	rm -r -f *~

tgz:
	rm -f ../www/src/SeSiMCMC.tgz
	tar --exclude *~ --exclude *tgz --exclude *o --exclude *exe -czvf ../www/src/SeSiMCMC.tgz ./* obj
