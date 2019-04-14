#toplevel make file

install:
	cd ./src && $(MAKE) all
	chmod +x bin/genA
debug:
	cd ./src && $(MAKE) debug

clean: 
	cd ./src && $(MAKE) clean	
