#!/bin/sh
	mkdir bin
	cd src/
	make
	mv -f sv_stat ../bin
	cd ../

