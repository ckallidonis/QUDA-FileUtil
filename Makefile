include ./make.inc

DIRS = extract convert

all: 
	-for d in $(DIRS); do (cd $$d && $(MAKE)); done

convert: 
	cd $@ && $(MAKE)

extract:
	cd $@ && $(MAKE)

clean:
	-for d in $(DIRS); do (cd $$d && $(MAKE) clean); done
