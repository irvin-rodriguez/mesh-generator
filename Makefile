CPPFLAGS=-Wall -Werror -pedantic -std=gnu++98 -ggdb3
RMOBJS=$(filter-out $(PROVIDED_OBJS), $(wildcard *.o))

all: mesh-step1 mesh-step2 mesh-step3 mesh-step4
	@echo "Built mesh-step1, mesh-step2, mesh-step3, and mesh-step4"

mesh-step1: mesh-step1.o mesh.o triangle.o
	g++ -o $@ $^

mesh-step2: mesh-step2.o mesh.o triangle.o
	g++ -o $@ $^

mesh-step3: mesh-step3.o mesh.o triangle.o
	g++ -o $@ $^

mesh-step4: mesh-step4.o mesh.o triangle.o
	g++ -o $@ $^

mesh-step1.o: mesh-step1.cpp mesh.hpp node.hpp triangle.hpp
	g++ $(CPPFLAGS) -c $<

mesh-step2.o: mesh-step2.cpp mesh.hpp node.hpp triangle.hpp
	g++ $(CPPFLAGS) -c $<

mesh-step3.o: mesh-step3.cpp mesh.hpp node.hpp triangle.hpp
	g++ $(CPPFLAGS) -c $<

mesh-step4.o: mesh-step4.cpp mesh.hpp node.hpp triangle.hpp
	g++ $(CPPFLAGS) -c $<

mesh.o: mesh.cpp mesh.hpp triangle.hpp node.hpp
	g++ $(CPPFLAGS) -c $<

triangle.o: triangle.cpp triangle.hpp node.hpp
	g++ $(CPPFLAGS) -c $<

.PHONY: clean format

clean:
	rm -f $(RMOBJS) *~ mesh-step1 mesh-step2 mesh-step3 mesh-step4

format:
	clang-format -i *.cpp *.hpp
