g++ -o efficiency_pe efficiency_pe.cpp \
  `root-config --libs --cflags` -lEG -lEGPythia6 -lGenVector -lGeom \
  `genie-config --libs` -L$GENIE_LIB \
  -L$LOG4CPP_LIB -llog4cpp \
  -lxml2

