#!/bin/sh
set -ex


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
	autoreconf --force --install || echo "autoreconf failed - try configure"
}

debug () {
	autogen;
	./configure  --enable-debug $@
	make
}

single () {
	autogen;
	./configure --enable-single $@
	make
}

default () {
	autoreconf --force --install || echo "autoreconf failed - try configure"
	./configure $@
	make
}

doxygen () {
	echo "In the country of the blind the one-eyed man is king."
	cat <<EOF > doxygen.conf
DOXYFILE_ENCODING      = UTF-8
PROJECT_NAME           = "PDEN - Protein Density"
PROJECT_NUMBER         =
PROJECT_BRIEF          =
PROJECT_LOGO           =
OUTPUT_DIRECTORY       = doc
OUTPUT_LANGUAGE        = English
OPTIMIZE_OUTPUT_FOR_C  = YES
INPUT                  = src
INPUT_ENCODING         = UTF-8
FILE_PATTERNS          = *.c *.h
GENERATE_LATEX         = NO
EOF
	doxygen doxygen.conf
}


clean () {
	make clean
	rm -rf config.log config.status configure depcomp missing install-sh autom4te.cache compile Makefile.in aclocal.m4 Makefile src/Makefile src/Makefile.in src/stamp-h1 src/config.h src/config.h.in src/config.h.in~ src/.deps 
	rm -rf doc doxygen.conf

}

help () {
set +x
echo "HELP:" 
echo "" 
echo "USAGE: autogen.sh <mode> <args>" 
echo "" 
echo "MODE:" 
echo "init     init autotools (autoscan) (dev only)"
echo "auto     run typical autotools (dev only)"
echo "reconf   run autoreconf and stop (create a configure script)"
echo "debug    run autoreconf and build with debugging flags"
echo "single   run autoreconf and build with single precision (float)"
echo "clean    remove all generated files"
echo "doxygen  create documentation using doxygen" 
echo "*        run autoreconf and build" 
echo "" 
echo "ARGS:" 
echo "all args are passed to configure if configure is executed" 
exit 0;
}


case $1 in
--help|-h|help) shift; help $@;;
init) shift; init $@;;
auto) shift; autotools $@;;
reconf) shift; autogen $@;;
debug) shift; debug $@;;
single) shift; single $@;;
clean) shift; clean $@;;
doxygen) shift; doxygen $@;;
*) default $@;;
esac
