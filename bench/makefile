# Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# Use the MFEM build directory
MFEM_ROOT ?= mfem_dir_is_not_defined
SRC = ./
CONFIG_MK = $(MFEM_ROOT)/build_make/config/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

SEQ_EXAMPLES = test_periodic test_full_periodic test_inclusion
PAR_EXAMPLES = 

ifeq ($(MFEM_USE_MPI),NO)
   EXAMPLES = $(SEQ_EXAMPLES)
else
   EXAMPLES = $(PAR_EXAMPLES) $(SEQ_EXAMPLES)
endif

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec

# Remove built-in rule
%: %.cpp

# Replace the default implicit rule for *.cpp files
%: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) -g $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

all: $(EXAMPLES) 

.PHONY: 

# Testing: Parallel vs. serial runs
RUN_MPI = $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(MFEM_MPI_NP)
%-test-par: %
	@$(call mfem-test,$<, $(RUN_MPI), Parallel example, "", "")
%-test-seq: %
	@$(call mfem-test,$<,, Serial example, -r 4, "")

# Testing: Specific execution options
xxx-test-seq: xxx
	@$(call mfem-test,$<,, Serial example, -m $(MFEM_DIR)/data/xxx.mesh -r 4, "")


# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built or MFEM_ROOT is not defined )

clean: clean-build clean-exec $(SUBDIRS_CLEAN)

clean-build:
	rm -f *.o *~ $(SEQ_EXAMPLES) $(PAR_EXAMPLES)
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -f *~
	rm -rf ex1dir L_shape_laplace_dir *vtk
