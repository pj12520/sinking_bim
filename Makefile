OBJS=main.o inout.o object.o testing.o
EXE=impact
HEAD=inout.h object.h const.h testing.h
all: $(EXE)

# This rule uses the automatic variables, '$^' and '$@'
$(EXE): $(OBJS)
	g++ $^ -o $@

# This is a 'static pattern rule'
$(OBJS): %.o : %.cc $(HEAD)
	g++ -c $< -o $@

.PHONY: clean spotless

clean:
	\rm -f $(OBJS)

spotless:
	\rm -f $(OBJS) $(EXE)
