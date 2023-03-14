# FLAGS= -fno-for-scope -O2 -fPIC
FLAGS= -fPIC -std=c++11
CPP= icpc
# IDENTIFY YOUR INCLUDE PATH WHICH CONTAIN BOOST AND PYTHON
p1=${CONDA_PREFIX}/include
p2=${p1}/python*
pl1=${CONDA_PREFIX}/lib

I=-I ${p1} -I ${p2}
L=-L ${pl1}

all: pyxaid_core.so

random.o: random.cpp random.h
	${CPP} ${FLAGS} ${I} -c random.cpp

matrix.o: matrix.cpp matrix.h
	${CPP} ${FLAGS} ${I} -c matrix.cpp

state.o: state.cpp state.h
	${CPP} ${FLAGS} ${I} -c state.cpp

io.o: io.cpp io.h
	${CPP} ${FLAGS} ${I} -c io.cpp

InputStructure.o: InputStructure.cpp InputStructure.h
	${CPP} ${FLAGS} ${I} -c InputStructure.cpp

ElectronicStructure.o: ElectronicStructure.cpp ElectronicStructure.h
	${CPP} ${FLAGS} ${I} -c ElectronicStructure.cpp

namd.o: namd.cpp namd.h
	${CPP} ${FLAGS} ${I} -c namd.cpp

namd_export.o: namd_export.cpp namd_export.h
	${CPP} ${FLAGS} ${I} -c namd_export.cpp

aux.o: aux.cpp aux.h
	${CPP} ${FLAGS} ${I} -c aux.cpp

wfc_basic_methods.o: wfc_basic_methods.cpp wfc.h
	${CPP} ${FLAGS} ${I} -c wfc_basic_methods.cpp

wfc_QE_methods.o: wfc_QE_methods.cpp wfc.h aux.cpp
	${CPP} ${FLAGS} ${I} -c wfc_QE_methods.cpp

wfc_functions.o: wfc_functions.cpp wfc.h
	${CPP} ${FLAGS} ${I} -c wfc_functions.cpp

wfc_export.o: wfc_basic_methods.o wfc_QE_methods.o wfc_functions.o wfc_export.cpp
	${CPP} ${FLAGS} ${I} -c wfc_export.cpp

pyxaid_core.o: pyxaid_core.cpp 
	${CPP} ${FLAGS} ${I} -c pyxaid_core.cpp



pyxaid_core.so: pyxaid_core.o wfc_export.o wfc_functions.o wfc_QE_methods.o wfc_basic_methods.o \
        aux.o matrix.o state.o ElectronicStructure.o namd.o namd_export.o InputStructure.o io.o random.o
	${CPP} ${FLAGS} ${I} -shared -o pyxaid_core.so pyxaid_core.o wfc_export.o wfc_functions.o \
        wfc_QE_methods.o wfc_basic_methods.o aux.o matrix.o state.o ElectronicStructure.o namd.o \
        namd_export.o InputStructure.o io.o random.o ${L} -lboost_python
	cp pyxaid_core.so ../.

clean:
	rm *.o
	rm pyxaid_core.so
	if [ -f ../pyxaid_core.so ]; then rm ../pyxaid_core.so; fi

