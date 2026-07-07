GSmulti		:	DF.o ExtMath.o FF.o GS.o IDLinterface.o Plasma.o MWmain.o
			g++ $^ -shared -o GSmulti.so
DF.o		:	DF.cpp DF.h ExtMath.h FF.h GS.h IDLinterface.h MWmain.h Plasma.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX DF.cpp
ExtMath.o	:	ExtMath.cpp DF.h ExtMath.h FF.h GS.h IDLinterface.h MWmain.h Plasma.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX ExtMath.cpp
FF.o		:	FF.cpp DF.h ExtMath.h FF.h GS.h IDLinterface.h MWmain.h Plasma.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX FF.cpp
GS.o		:	GS.cpp DF.h ExtMath.h FF.h GS.h IDLinterface.h MWmain.h Plasma.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX GS.cpp
IDLinterface.o	:	IDLinterface.cpp DF.h ExtMath.h FF.h GS.h IDLinterface.h MWmain.h Plasma.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX IDLinterface.cpp
Plasma.o	:	Plasma.cpp DF.h ExtMath.h FF.h GS.h IDLinterface.h MWmain.h Plasma.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX Plasma.cpp
MWmain.o	:	MWmain.cpp DF.h ExtMath.h FF.h GS.h IDLinterface.h MWmain.h Plasma.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX MWmain.cpp