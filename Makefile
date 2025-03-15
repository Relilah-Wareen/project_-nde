.PHONY: build run plot report clean

build:
	cd test && make build

run:
	cd test && make run 

plot:
	cd test && make plot 

clean:
	cd test && make clean
	cd report && make clean

report:
	cd report && make report
