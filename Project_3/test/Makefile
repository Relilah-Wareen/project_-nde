build:
	g++ test.cpp ../src/ODESolver.cpp ../src/IVPFactory.cpp ../src/LinearMultistepMethod.cpp ../src/EmbeddedRKMethod.cpp ../src/RKMethod.cpp -o test
	g++ convtest.cpp ../src/IVPFactory.cpp ../src/LinearMultistepMethod.cpp ../src/ODESolver.cpp ../src/EmbeddedRKMethod.cpp ../src/RKMethod.cpp -o convtest

run:
	./test
	./convtest
	python convana.py
	
plot:
	python plot.py

clean:
	rm -f test
	rm -f convtest