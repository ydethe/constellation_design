Q=@

SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp,%.o,$(SRC))
CPPFLAGS := -Os

.PHONY: default run clean

default: test.exe

test.exe: $(OBJ)
	@echo Linking $@
	$(Q)g++ -o $@ $^ -lm -lnlopt -lgsl

%.o: %.cpp
	@echo Compiling $@
	$(Q)g++ -c $< -o $@ -I/usr/local/include/eigen3 $(CPPFLAGS)

run: test.exe
	@echo Running $<
	$(Q)./$<

clean:
	$(Q)rm -rf $(OBJ) test.exe
