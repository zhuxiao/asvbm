CXXFLAGS = -O0 -g -Wall -fmessage-length=0 -std=c++11

OBJS = sv_stat.o dataLoader.o convert.o ref_reg_size_stat.o num_stat.o \
	   size_dif_stat.o type_num_stat.o util.o

LIBS +=

TARGET = sv_stat

all: $(TARGET) clean

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) 

clean:
	rm -f $(OBJS)
	
clean-all: clean
	rm -f $(TARGET)
