all: snaTMB.so

snaTMB.so: snaTMB.cpp
	Rscript --vanilla -e "TMB::compile('snaTMB.cpp','-O0 -g')"

clean:
	rm -rf *o
