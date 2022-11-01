# Makefile for code_implementation

CC = gcc
PATH_INCLUDE = ../include_file
PATH_SOURCE = ../source_file
OBJECT = $(PATH_SOURCE)/main_implement.o \
	 $(PATH_SOURCE)/linear_system.o \
	 $(PATH_SOURCE)/linear_system_solver.o \
	 $(PATH_SOURCE)/pss_solver.o \
	 $(PATH_SOURCE)/pss_solver_iteration.o \
	 $(PATH_SOURCE)/pss_solver_iterative_scheme.o \
	 $(PATH_SOURCE)/pss_cg_solver.o \
	 $(PATH_SOURCE)/pss_bicgstab_solver.o \
	 $(PATH_SOURCE)/pss_gmres_solver.o \
	 $(PATH_SOURCE)/gmres_restart_iteration.o

matrix_split_solver: $(OBJECT)
	$(CC) -lm -o matrix_split_solver $(OBJECT)

$(PATH_SOURCE)/main_implement.o: $(PATH_SOURCE)/main_implement.c \
		$(PATH_INCLUDE)/main_implement.h

$(PATH_SOURCE)/linear_system.o: $(PATH_SOURCE)/linear_system.c \
		$(PATH_INCLUDE)/linear_system.h

$(PATH_SOURCE)/linear_system_solver.o: $(PATH_SOURCE)/linear_system_solver.c \
		$(PATH_INCLUDE)/linear_system.h \
		$(PATH_INCLUDE)/linear_system_solver.h

$(PATH_SOURCE)/pss_solver.o: $(PATH_SOURCE)/pss_solver.c \
		$(PATH_INCLUDE)/pss_solver.h \
		$(PATH_INCLUDE)/linear_system.h

$(PATH_SOURCE)/pss_solver_iteration.o: $(PATH_SOURCE)/pss_solver_iteration.c \
		$(PATH_INCLUDE)/pss_solver_iteration.h \
		$(PATH_INCLUDE)/linear_system.h

$(PATH_SOURCE)/pss_solver_iterative_scheme.o: $(PATH_SOURCE)/pss_solver_iterative_scheme.c \
		$(PATH_INCLUDE)/pss_solver_iterative_scheme.h \
		$(PATH_INCLUDE)/linear_system.h

$(PATH_SOURCE)/pss_cg_solver.o: $(PATH_SOURCE)/pss_cg_solver.c \
		$(PATH_INCLUDE)/pss_cg_solver.h \
		$(PATH_INCLUDE)/pss_solver_iteration.h

$(PATH_SOURCE)/pss_bicgstab_solver.o: $(PATH_SOURCE)/pss_bicgstab_solver.c \
		$(PATH_INCLUDE)/pss_bicgstab_solver.h \
		$(PATH_INCLUDE)/pss_solver_iteration.h

$(PATH_SOURCE)/pss_gmres_solver.o: $(PATH_SOURCE)/pss_gmres_solver.c \
		$(PATH_INCLUDE)/pss_gmres_solver.h \
		$(PATH_INCLUDE)/pss_solver_iteration.h

$(PAHT_SOURCE)/gmres_restart_iteration.o: $(PATH_SOURCE)/gmres_restart_iteration.c \
		$(PATH_INCLUDE)/gmres_restart_iteration.h \
		$(PATH_INCLUDE)/pss_gmres_solver.h

clean:
	rm matrix_split_solver $(OBJECT)
