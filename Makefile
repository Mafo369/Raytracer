.PHONY:all clean
DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
#$(shell cp ./lodepng-master/lodepng.cpp ./ >/dev/null)
DEPFLAGS = -MT $@ -MMD -MF $(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

TARGETS=mrt unit-test

all: $(TARGETS)

clean:
	@rm -rf *.o ./lodepng-master/lodepng.o *~ $(TARGETS) .d

CC=g++
CFLAGS=-Wall -std=c++14 -g -I./glm-master/ -fopenmp -I./lodepng-master/ -O3
SRCS=main.cpp image.cpp raytracer.cpp scene.cpp kdtree.cpp ./lodepng-master/lodepng.cpp unit-test.cpp

OBJ=main.o

%.o : %.cpp

%.o: %.cpp $(DEPDIR)/%.d
	$(CC) -c $(CFLAGS) $(DEPFLAGS) $*.cpp -o $*.o
	$(POSTCOMPILE)

./lodepng-master/%.o : ./lodepng-master/%.cpp

./lodepng-master/%.o: ./lodepng-master/%.cpp $(DEPDIR)/%.d
	$(CC) -c $(CFLAGS) $(DEPFLAGS) ./lodepng-master/$*.cpp -o ./lodepng-master/$*.o
	$(POSTCOMPILE)

mrt: main.o image.o scene.o raytracer.o kdtree.o ./lodepng-master/lodepng.o
	$(CC) $(CFLAGS) $^ -o $@

unit-test: unit-test.o image.o raytracer.o scene.o raytracer.o kdtree.o ./lodepng-master/lodepng.o
	$(CC) $(CFLAGS) $^ -o $@

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS))))
