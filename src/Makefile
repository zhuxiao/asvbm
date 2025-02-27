CXXFLAGS = -O0 -g -Wall -fmessage-length=0 -std=c++11

OBJS = asvbm_main.o paras.o dataLoader.o convert.o ref_reg_size_stat.o \
	   num_stat.o size_dif_stat.o size_num_stat.o util.o meminfo.o gnuplotcall.o htmlGenerate.o shared_fp.o benchmark_integration.o edlib.o

LIBS += -lhts -lpthread 

TARGET = asvbm

all: $(TARGET) clean

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) 

clean:
	rm -f $(OBJS)
	
clean-all: clean
	rm -f $(TARGET)
