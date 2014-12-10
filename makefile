CXX ?= g++
CFLAGS = -O3 -fPIC
#CFLAGS = -Wall -Wconversion -O3 -fPIC


SHVER = 2
OS = $(shell uname)
LIBSVM_DIR = ./libsvm-weights-3.18/
all: bonsai

lib: svm.o
	if [ "$(OS)" = "Darwin" ]; then \
		SHARED_LIB_FLAG="-dynamiclib -Wl,-install_name,libsvm.so.$(SHVER)"; \
	else \
		SHARED_LIB_FLAG="-shared -Wl,-soname,libsvm.so.$(SHVER)"; \
	fi; \
	$(CXX) $${SHARED_LIB_FLAG} svm.o -o libsvm.so.$(SHVER)

svm.o: $(LIBSVM_DIR)/svm.cpp $(LIBSVM_DIR)/svm.h
	$(CXX) $(CFLAGS) -c $(LIBSVM_DIR)/svm.cpp
bonsai: main.cpp bonsai.cpp bonsai.h fom.cpp fom.h libsvm_container.h svm_interface.h svm.o
	$(CXX) $(CFLAGS) -fopenmp main.cpp svm.o bonsai.cpp fom.cpp svm_interface.cpp -I$(ROOTSYS)/include -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic  -o bonsai_svm -lm 
clean:
	rm -f *~ svm.o svm-train svm-predict svm-scale bonsai_svm libsvm.so.$(SHVER)