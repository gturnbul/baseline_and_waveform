TARGET=waveform_type_cuts.exe
BOOST_INCLUDE=/opt/homebrew/Cellar/boost/1.87.0/include
BOOST_LIB=/opt/homebrew/Cellar/boost/1.87.0/lib

all: $(TARGET)

waveform_type_cuts.o: waveform_type_cuts.C
	g++ -c -I$(BOOST_INCLUDE) `root-config --cflags` waveform_type_cuts.C -o waveform_type_cuts.o

myDictionary.o: myDictionary.cxx
	g++ -c -I$(BOOST_INCLUDE) `root-config --cflags` myDictionary.cxx -o myDictionary.o

$(TARGET): waveform_type_cuts.o myDictionary.o
	g++ -o $(TARGET) waveform_type_cuts.o myDictionary.o -I$(BOOST_INCLUDE) -L$(BOOST_LIB) `root-config --libs` -lboost_system -lSpectrum

clean:
	rm -f $(TARGET) waveform_type_cuts.o myDictionary.o
