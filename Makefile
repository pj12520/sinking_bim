OBJS=main.o geo.o inout.o object.o testing.o ellip.o axisym.o build.o build_sup.o solve.o vel.o 
EXE=impact
HEAD=geo.h inout.h object.h const.h testing.h ellip.h axisym.h build.h build_sup.h solve.h vel.h  interp_1d.h dfridr.h interp_curve.h
all: $(EXE)

# This rule uses the automatic variables, '$^' and '$@'
$(EXE): $(OBJS)
	g++ $^ -lgsl -lgslcblas -o $@

# This is a 'static pattern rule'
$(OBJS): %.o : %.cc $(HEAD)
	g++ -c $< -o $@

.PHONY: clean spotless

clean:
	\rm -f $(OBJS)

spotless:
	\rm -f $(OBJS) $(EXE)
