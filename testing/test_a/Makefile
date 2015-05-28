OBJS=test_a.o ../../axisym.o ../../ellip.o
EXE=test_a
HEAD=../../const.h ../../axisym.h ../../ellip.h
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
