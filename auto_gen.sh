#!/bin/sh
	if [ ! -d "bin" ]; then 
		mkdir bin
	fi
	cd src/
	make
	mv -f sv_stat ../bin
	cd ../

