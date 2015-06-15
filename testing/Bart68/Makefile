OBJS=main.o # inout.o object.o testing.o ellip.o axisym.o build.o build_sup.o solve.o
EXE=Bart_test
HEAD='# inout.h object.h const.h testing.h ellip.h axisym.h build.h build_sup.h solve.h
all: $(EXE)

# This rule uses the automatic variables, '$^' and '$@'
$(EXE): $(OBJS)
	g++ $^ -lgsl -lgslcblas -o $@

# This is a 'static pattern rule'
$(OBJS): %.o : %.cc 
	g++ -c $< -o $@

.PHONY: clean spotless

clean:
	\rm -f $(OBJS)

spotless:
	\rm -f $(OBJS) $(EXE)
