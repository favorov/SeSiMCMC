#this file can be started only from parent directory
name=DiagnosticsTest
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

OBJS=$(od)/$(name).o 

$(od)/%.o: %.c
	$(CC) $(CCFLAGS) $< -o $@
	
$(od)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $< -o $@

$(name)$(EXEEXT): $(OBJS)
	$(CPP) -o $(name)$(EXEEXT) $(OBJS) $(LINKFLAGS) 

$(od)/$(name).o: $(md)/Diagnostics.hpp

$(od)/$(name).o : $(td)/$(name).cpp

clean:
	rm -f $(OBJS)
	rm -r -f *~

