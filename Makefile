PYTHON_INCLUDE_DIR = `python-config --includes`

all:
	swig -python -c++ potentials.i
	g++ -fPIC -std=c++11 -shared -o _potentials.so potentials_wrap.cxx -I $(PYTHON_INCLUDE_DIR) -I.

clean:
	rm -rf *.so *.py *.pyc *.cxx
