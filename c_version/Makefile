Q=@

SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp,%.o,$(SRC))


.PHONY: default run clean

default: test.exe

test.exe: $(OBJ)
	@echo Linking $@
	$(Q)g++ -o $@ $? -lm -lnlopt

%.o: %.cpp
	@echo Compiling $@
	$(Q)g++ -c $< -o $@ -I/usr/include/eigen3

run: test.exe
	@echo Running $<
	$(Q)./$<

clean:
	$(Q)rm -rf $(OBJ) test.exe
