#!/bin/sh
	if [ ! -d "bin" ]; then 
		mkdir bin
	fi
	cd src/
	make
	mv -f asvbm ../bin
	cd ../script
	chmod a+x callset_ensemble.sh
	ln -sf $PWD/callset_ensemble.sh ../bin
	cd ..

