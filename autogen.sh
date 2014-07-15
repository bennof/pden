#!/bin/sh
set -ex


lib_src="\
fft.c \
gradient.c \
pden.c \
pden_math.c \
pden_mrc.c \
pden_xplor.c \
powerspectra.c \
render.c \
sfrefine.c \
babinet.c \
split2fft.c"

bin_src=" \
pdcalcps.c \
pdrefinesf.c \
pdapplysf.c \
pdnormalize.c"


fbuild () {
	printf "Fast build ...\n"
	for i in src/*.c; do 
		gcc -Wall -std=c99 -O2  -c -pedantic -DVERBOE=1 -DDOUBLE=double -Dreal=double -o ${i%*.c}.o $i  || return 1; 
	done

	for i in $lib_src
	do
		lib_obj="$lib_obj src/${i%*.c}.o"
	done
	ar crs libpden.a $lib_obj

	for i in $bin_src
	do
		gcc -Wall -O2 -o ${i%*.c} src/${i%*.c}.o libpden.a -lfftw3 -lm
	done
}


dbuild () {
	printf "Debug build ...\n"
	for i in src/*.c; do 
		gcc -Wall -std=c99 -g  -c -pedantic -DVERBOSE=1 -DDOUBLE=double -Dreal=double -o ${i%*.c}.o $i  || return 1; 
	done

	for i in $lib_src
	do
		lib_obj="$lib_obj src/${i%*.c}.o"
	done
	ar crs libpden.a $lib_obj

	for i in $bin_src
	do
		gcc -Wall -g -o ${i%*.c} src/${i%*.c}.o libpden.a -lfftw3 -lm
	done
}

clean_build () {
	rm src/*.o
	rm libpden.a
	for i in $bin_src
	do
		rm ${i%*.c}
	done
}

init () {
	printf "Init ...\n"
	autoscan
}

autotools () {
	printf "AutoGen ...\n"
	aclocal &&\
	autoheader &&\
	automake --add-missing &&\
	autoconf
}

autogen () {
	autoreconf --force --install -m
}


default () {
	autogen;
	cp src/libpden.a .
}

clean () {
	rm -rf config.log config.status configure depcomp missing install-sh autom4te.cache compile Makefile.in aclocal.m4 Makefile src/Makefile src/Makefile.in src/stamp-h1 src/config.h src/config.h.in src/config.h.in~ src/.deps 
	clean_build
}

case $1 in
init) shift; init $@;;
fbuild) shift; fbuild $@;;
dbuild) shift; dbuild $@;;
auto) shift; autotools $@;;
gen) shift; autogen $@;;
clean) shift; clean $@;;
*) default $@;;
esac
