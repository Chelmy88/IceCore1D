#! /bin/bash

mkdir prof
CPUPROFILE=prof/program_name.prof DYLD_INSERT_LIBRARIES=/usr/local/Cellar/gperftools/2.7/lib/libprofiler.dylib bin/ice_model
pprof --pdf bin/ice_model prof/program_name.prof > prof/program_name.pdf
