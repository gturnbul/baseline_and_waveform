TARGET=OffsetCalc.exe
BOOST_INCLUDE=/opt/homebrew/Cellar/boost/1.87.0/include
BOOST_LIB=/opt/homebrew/Cellar/boost/1.87.0/lib

all: $(TARGET)

OffsetCalc.o: OffsetCalc.C
	g++ -c -I$(BOOST_INCLUDE) `root-config --cflags` OffsetCalc.C -o OffsetCalc.o

myDictionary.o: myDictionary.cxx
	g++ -c -I$(BOOST_INCLUDE) `root-config --cflags` myDictionary.cxx -o myDictionary.o

$(TARGET): OffsetCalc.o myDictionary.o
	g++ -o $(TARGET) OffsetCalc.o myDictionary.o -I$(BOOST_INCLUDE) -L$(BOOST_LIB) `root-config --libs` -lboost_system -lSpectrum

clean:
	rm -f $(TARGET) OffsetCalc.o myDictionary.o
