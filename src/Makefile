all: semtmb.so

semtmb.so: semtmb.cpp
	Rscript --vanilla -e "TMB::compile('semtmb.cpp','-O0 -g')"

clean:
	rm -rf *o
