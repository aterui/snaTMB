all: snglmm.so

snglmm.so: snglmm.cpp
	Rscript --vanilla -e "TMB::compile('snglmm.cpp','-O0 -g')"

clean:
	rm -rf *o
