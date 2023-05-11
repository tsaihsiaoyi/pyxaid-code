FLAGS= -fPIC -std=c++11
CPP= ${CXX}

pyxaid: obj_dir pyxaid.o aux.o matrix.o state.o ElectronicStructure.o namd.o InputStructure.o io.o random.o
	${CPP} ${FLAGS} -o ../pyxaid objs/*.o

obj_dir:
	mkdir -p objs

random.o: random.cpp random.h
	${CPP} ${FLAGS} -c random.cpp -o objs/random.o

matrix.o: matrix.cpp matrix.h
	${CPP} ${FLAGS} -c matrix.cpp -o objs/matrix.o

state.o: state.cpp state.h
	${CPP} ${FLAGS} -c state.cpp -o objs/state.o

io.o: io.cpp io.h
	${CPP} ${FLAGS} -c io.cpp -o objs/io.o

InputStructure.o: InputStructure.cpp InputStructure.h
	${CPP} ${FLAGS} -c InputStructure.cpp -o objs/InputStructure.o

ElectronicStructure.o: ElectronicStructure.cpp ElectronicStructure.h
	${CPP} ${FLAGS} -c ElectronicStructure.cpp -o objs/ElectronicStructure.o

namd.o: namd.cpp namd.h
	${CPP} ${FLAGS} -c namd.cpp -o objs/namd.o

aux.o: aux.cpp aux.h
	${CPP} ${FLAGS} -c aux.cpp -o objs/aux.o

pyxaid.o: pyxaid.cpp 
	${CPP} ${FLAGS} -c pyxaid.cpp -o objs/pyxaid.o

clean:
	if [ -d objs ]; then rm -rf objs; fi

distclean: clean
	if [ -f ../pyxaid ]; then rm ../pyxaid; fi