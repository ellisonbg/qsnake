all: rdirac.so

rdirac.so: rdirac.f
	f2py --fcompiler=gnu95 --f77flags="-W -Wall -fdefault-real-8" -c -m rdirac rdirac.f

clean:
	py.cleanup
	rm -f rdirac.so
