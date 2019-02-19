CXXFLAGS = -O0 -g -Wall -fmessage-length=0 -std=c++11

OBJS = sv_stat.o

LIBS +=

TARGET = sv_stat

all: $(TARGET) clean

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) 

clean:
	rm -f $(OBJS)
	
clean-all: clean
	rm -f $(TARGET)
