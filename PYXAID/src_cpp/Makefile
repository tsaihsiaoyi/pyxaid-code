# FLAGS= -fno-for-scope -O2 -fPIC
FLAGS= -fPIC -std=c++11
CPP= ${CXX}

all: pyxaid

random.o: random.cpp random.h
	${CPP} ${FLAGS} -c random.cpp

matrix.o: matrix.cpp matrix.h
	${CPP} ${FLAGS} -c matrix.cpp

state.o: state.cpp state.h
	${CPP} ${FLAGS} -c state.cpp

io.o: io.cpp io.h
	${CPP} ${FLAGS} -c io.cpp

InputStructure.o: InputStructure.cpp InputStructure.h
	${CPP} ${FLAGS} -c InputStructure.cpp

ElectronicStructure.o: ElectronicStructure.cpp ElectronicStructure.h
	${CPP} ${FLAGS} -c ElectronicStructure.cpp

namd.o: namd.cpp namd.h
	${CPP} ${FLAGS} -c namd.cpp

aux.o: aux.cpp aux.h
	${CPP} ${FLAGS} -c aux.cpp

pyxaid_core.o: pyxaid_core.cpp 
	${CPP} ${FLAGS} -c pyxaid_core.cpp

pyxaid: pyxaid_core.o aux.o matrix.o state.o ElectronicStructure.o namd.o InputStructure.o io.o random.o
	${CPP} ${FLAGS} -o pyxaid pyxaid_core.o aux.o matrix.o state.o ElectronicStructure.o namd.o InputStructure.o io.o random.o

clean:
	rm *.o
	rm pyxaid