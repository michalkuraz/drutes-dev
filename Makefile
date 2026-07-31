#                             _____________________  _______________________
#                             ___  __ \__  __ \_  / / /_  /___  ____/_  ___/
#                             __  / / /_  /_/ /  / / /_  __/_  __/  _____ \ 
#                             _  /_/ /_  _, _// /_/ / / /_ _  /___  ____/ / 
#                             /_____/ /_/ |_| \____/  \__/ /_____/  /____/  
#                                                                           
# 
#---------------------------------------------D R U t E S-----------------------------------------
#                             (Dual Richards' Unsaturated Equation Solver)
#                                           M a k e f i l e 

.DEFAULT_GOAL := all

SHELL := /bin/bash
#build directories
BUILD = build
OBJDIR = $(BUILD)/objs
MODDIR = $(BUILD)/mods
BINDIR = bin

LOGDIR = $(BUILD)/logs
LOGFILE = $(LOGDIR)/compile.log

$(BUILD):
	mkdir -p $(BUILD)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

$(BINDIR):
	mkdir -p $(BINDIR)
	
$(LOGDIR):
	mkdir -p $(LOGDIR)
	
	
# ---------- NetCDF detection ----------
#HAVE_NETCDF := $(shell command -v nf-config >/dev/null 2>&1 && echo yes || echo no)
# ---------- NetCDF detection ----------
HAVE_NETCDF := $(shell command -v nf-config >/dev/null 2>&1 && \
                        command -v nc-config >/dev/null 2>&1 && echo yes || echo no)

ifeq ($(HAVE_NETCDF),yes)
  NETCDF_FFLAGS := $(shell nf-config --fflags)
  NETCDF_FLIBS := -L$(shell nc-config --libdir) $(shell nf-config --flibs)
  CPPFLAGS_NETCDF := -DHAVE_NETCDF
  NETCDF_MSG    := compiled with NetCDF support
else
  NETCDF_FFLAGS :=
  NETCDF_FLIBS  :=
  CPPFLAGS_NETCDF :=
  NETCDF_MSG    := compiled without NetCDF support
endif

# -------- compiler --------
FC = gfortran

# -------- debugging flags (development) --------
# FFLAGS = $(CPPFLAGS_NETCDF) -fimplicit-none -fcoarray=single -fbounds-check -fbacktrace -g \
        -fdefault-real-8 -O0 -finit-real=nan -Wsurprising -J$(MODDIR) $(NETCDF_FFLAGS)

# -------- optimized flags (production) --------
FFLAGS = $(CPPFLAGS_NETCDF) -fimplicit-none -fcoarray=single -fdefault-real-8 -O3 \
         -finit-real=nan -ffpe-summary=none -fno-backtrace \
         -J$(MODDIR) $(NETCDF_FFLAGS)
         
         
d=drutes_obj-`date -I`

all: | $(LOGDIR)
	@set -o pipefail; \
	start=$$(date +%s.%N); \
	echo "==== DRUtES compilation log ====" | tee $(LOGFILE); \
	echo "Compiler: $(FC)" | tee -a $(LOGFILE); \
	echo "Flags:    $(FFLAGS)" | tee -a $(LOGFILE); \
	echo "" | tee -a $(LOGFILE); \
	if $(MAKE) --no-print-directory build_target 2>&1 | tee -a $(LOGFILE); then \
		status="SUCCESS"; \
	else \
		status="FAILED"; \
	fi; \
	end=$$(date +%s.%N); \
	elapsed=$$(awk "BEGIN {print $$end - $$start}"); \
	echo "" | tee -a $(LOGFILE); \
	echo "Elapsed: $$elapsed seconds" | tee -a $(LOGFILE); \
	echo "Status:  $$status" | tee -a $(LOGFILE); \
	echo "===================================" | tee -a $(LOGFILE); \
	if [ "$$status" = "FAILED" ]; then exit 1; fi

build_target: $(BINDIR)/drutes

$(BINDIR)/drutes: $(OBJDIR)/main.o $(ALL_objs) | $(BINDIR)
	$(FC) $(FFLAGS) -o $@ $(OBJDIR)/main.o $(ALL_objs) $(NETCDF_FLIBS)

#----------------objects definitions-------------------------------
CORE_obj := $(OBJDIR)/typy.o $(OBJDIR)/global_objs.o $(OBJDIR)/globals.o $(OBJDIR)/globals1D.o $(OBJDIR)/globals2D.o  $(OBJDIR)/debug_tools.o $(OBJDIR)/core_tools.o $(OBJDIR)/pde_objs.o $(OBJDIR)/dummy_procs.o $(OBJDIR)/global4solver.o
POINTERMAN_obj := $(OBJDIR)/manage_pointers.o
RE_obj := $(OBJDIR)/re_constitutive.o $(OBJDIR)/re_reader.o $(OBJDIR)/re_globals.o $(OBJDIR)/re_total.o $(OBJDIR)/re_pointers.o $(OBJDIR)/re_analytical.o $(OBJDIR)/re_evap_methods.o
MATHTOOLS_obj :=  $(OBJDIR)/linalg.o $(OBJDIR)/integral.o $(OBJDIR)/solver_interfaces.o $(OBJDIR)/simplelinalg.o $(OBJDIR)/gmres_solver.o
TOOLS_obj := $(OBJDIR)/printtools.o $(OBJDIR)/simegen.o $(OBJDIR)/read_inputs.o $(OBJDIR)/drutes_init.o $(OBJDIR)/geom_tools.o $(OBJDIR)/postpro.o $(OBJDIR)/readtools.o $(OBJDIR)/objfnc.o $(OBJDIR)/datetime.o
FEMTOOLS_obj := $(OBJDIR)/feminittools.o $(OBJDIR)/capmat.o $(OBJDIR)/stiffmat.o $(OBJDIR)/fem.o $(OBJDIR)/fem_tools.o $(OBJDIR)/femmat.o
DECOMPO_obj :=  $(OBJDIR)/decomp_tools.o $(OBJDIR)/schwarz_dd.o  $(OBJDIR)/decomp_vars.o $(OBJDIR)/decomposer.o $(OBJDIR)/schwarz_dd2subcyc.o
PMAoo_obj := $(OBJDIR)/fullmatrix.o $(OBJDIR)/mtx.o $(OBJDIR)/mtx_int.o $(OBJDIR)/mtxiotools.o $(OBJDIR)/pmatools.o $(OBJDIR)/solvers.o $(OBJDIR)/sparsematrix.o $(OBJDIR)/sparsematrix_int.o $(OBJDIR)/matmod.o $(OBJDIR)/reorder.o
BOUSSINESQ_obj := $(OBJDIR)/boussglob.o $(OBJDIR)/boussread.o $(OBJDIR)/boussfnc.o $(OBJDIR)/bousspointers.o
ADE_obj := $(OBJDIR)/ADE_fnc.o $(OBJDIR)/ADE_reader.o $(OBJDIR)/ADE_globals.o $(OBJDIR)/ADE_pointers.o
REDUAL_obj := $(OBJDIR)/Re_dual_totH.o $(OBJDIR)/Re_dual_globals.o $(OBJDIR)/Re_dual_pointers.o $(OBJDIR)/Re_dual_reader.o $(OBJDIR)/Re_dual_tab.o $(OBJDIR)/Re_dual_coupling.o $(OBJDIR)/Re_dual_bc.o
HEAT_obj := $(OBJDIR)/heat_fnc.o $(OBJDIR)/heat_pointers.o $(OBJDIR)/heat_globals.o $(OBJDIR)/heat_reader.o
KINWAVE_obj := $(OBJDIR)/kinreader.o $(OBJDIR)/kinglobs.o $(OBJDIR)/kinfnc.o $(OBJDIR)/kinpointer.o
FROZEN_obj := $(OBJDIR)/freeze_globs.o $(OBJDIR)/freeze_helper.o $(OBJDIR)/freeze_fnc.o $(OBJDIR)/freeze_reader.o $(OBJDIR)/freeze_pointers.o 
REevap_obj :=  $(OBJDIR)/evapglob.o $(OBJDIR)/evappointers.o $(OBJDIR)/evap_RE_constitutive.o $(OBJDIR)/evap_heat_constitutive.o $(OBJDIR)/evapreader.o $(OBJDIR)/evapbc4heat.o

ifeq ($(HAVE_NETCDF),yes)

	NETCDF_obj := $(OBJDIR)/init_netcdf.o $(OBJDIR)/ncglobvars.o $(OBJDIR)/netcdfflux.o $(OBJDIR)/ncpointers.o $(OBJDIR)/nctools.o $(OBJDIR)/ncdem.o  $(OBJDIR)/ncfluxarea.o $(OBJDIR)/lsconstitutive.o $(OBJDIR)/ncmesh.o $(OBJDIR)/ncmap.o

else

	NETCDF_obj :=

endif

MODEL_objs := $(RE_obj) $(BOUSSINESQ_obj) $(ADE_obj) $(REDUAL_obj) $(HEAT_obj) $(FROZEN_obj) $(KINWAVE_obj) $(REevap_obj)  $(NETCDF_obj)

ALL_objs := $(CORE_obj) $(TOOLS_obj) $(POINTERMAN_obj) $(MATHTOOLS_obj) $(FEMTOOLS_obj) $(DECOMPO_obj)  $(PMAoo_obj) $(MODEL_objs) 
#-----------------------------------------------------------------

#-------begin CORE_obj--------------------------------
$(OBJDIR)/typy.o: src/core/typy.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/typy.f90 -o $@
	
$(OBJDIR)/global_objs.o: $(OBJDIR)/typy.o $(PMAoo_obj) src/core/global_objs.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/global_objs.f90 -o $@
	
$(OBJDIR)/global4solver.o: $(OBJDIR)/typy.o src/core/global4solver.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/global4solver.f90 -o $@
	
$(OBJDIR)/pde_objs.o: $(OBJDIR)/typy.o $(OBJDIR)/global_objs.o $(PMAoo_obj) $(OBJDIR)/globals.o $(OBJDIR)/decomp_vars.o src/core/pde_objs.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/pde_objs.f90 -o $@
	
$(OBJDIR)/globals.o: $(OBJDIR)/typy.o $(OBJDIR)/global_objs.o src/core/globals.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/globals.f90 -o $@

$(OBJDIR)/globals1D.o: $(OBJDIR)/typy.o $(OBJDIR)/global_objs.o src/core/globals1D.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/globals1D.f90 -o $@
	
$(OBJDIR)/globals2D.o: $(OBJDIR)/typy.o $(OBJDIR)/global_objs.o src/core/globals2D.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/globals2D.f90 -o $@
	
$(OBJDIR)/core_tools.o: $(OBJDIR)/typy.o $(OBJDIR)/global_objs.o $(OBJDIR)/globals.o src/core/core_tools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/core_tools.f90 -o $@
	
$(OBJDIR)/dummy_procs.o: $(OBJDIR)/typy.o $(OBJDIR)/global_objs.o $(OBJDIR)/globals.o $(OBJDIR)/pde_objs.o src/core/dummy_procs.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/dummy_procs.f90 -o $@

$(OBJDIR)/debug_tools.o: $(OBJDIR)/typy.o $(OBJDIR)/core_tools.o src/core/debug_tools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/debug_tools.f90 -o $@

#---------end CORE_obj------------------------------


#------begin MATHTOOLS_obj-----------------------------
$(OBJDIR)/linalg.o: $(CORE_obj) src/mathtools/linalg.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/mathtools/linalg.f90 -o $@

$(OBJDIR)/gmres_solver.o: $(CORE_obj) src/mathtools/gmres_solver.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/mathtools/gmres_solver.f90 -o $@
	
$(OBJDIR)/integral.o: $(CORE_obj) $(OBJDIR)/linalg.o src/mathtools/integral.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/mathtools/integral.f90 -o $@

$(OBJDIR)/simplelinalg.o: $(CORE_obj) $(PMAoo_obj) $(OBJDIR)/re_globals.o $(OBJDIR)/linalg.o src/mathtools/simplelinalg.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/mathtools/simplelinalg.f90 -o $@

$(OBJDIR)/solver_interfaces.o: $(CORE_obj) $(PMAoo_obj) $(OBJDIR)/readtools.o $(OBJDIR)/simplelinalg.o $(OBJDIR)/gmres_solver.o src/mathtools/solver_interfaces.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/mathtools/solver_interfaces.f90 -o $@

#------end MATHTOOLS_obj---------------------------------


#--------begin PMAoo_obj------------------------
$(OBJDIR)/pmatools.o: $(OBJDIR)/typy.o src/pma++/pmatools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/pmatools.f90 -o $@

$(OBJDIR)/mtx.o: $(OBJDIR)/typy.o $(OBJDIR)/pmatools.o src/pma++/mtx.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/mtx.f90 -o $@

$(OBJDIR)/mtx_int.o: $(OBJDIR)/typy.o $(OBJDIR)/pmatools.o src/pma++/mtx_int.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/mtx_int.f90 -o $@

$(OBJDIR)/mtxiotools.o: $(OBJDIR)/typy.o src/pma++/mtxiotools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/mtxiotools.f90 -o $@

$(OBJDIR)/fullmatrix.o: $(OBJDIR)/typy.o $(OBJDIR)/mtx.o src/pma++/fullmatrix.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/fullmatrix.f90 -o $@
	
$(OBJDIR)/sparsematrix.o: $(OBJDIR)/typy.o $(OBJDIR)/mtx.o src/pma++/sparsematrix.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/sparsematrix.f90 -o $@	

$(OBJDIR)/sparsematrix_int.o: $(OBJDIR)/typy.o $(OBJDIR)/mtx.o src/pma++/sparsematrix_int.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/sparsematrix_int.f90 -o $@
	
$(OBJDIR)/solvers.o: $(OBJDIR)/global4solver.o $(OBJDIR)/typy.o $(OBJDIR)/mtx.o src/pma++/solvers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/solvers.f90 -o $@	
	
$(OBJDIR)/matmod.o: $(OBJDIR)/typy.o $(OBJDIR)/mtx.o src/pma++/matmod.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/matmod.f90 -o $@

$(OBJDIR)/datasetup.o: $(OBJDIR)/typy.o $(OBJDIR)/mtx.o src/pma++/datasetup.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/datasetup.f90 -o $@

$(OBJDIR)/reorder.o: $(OBJDIR)/typy.o $(OBJDIR)/mtx.o $(OBJDIR)/datasetup.o $(OBJDIR)/solvers.o src/pma++/reorder.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/pma++/reorder.f90 -o $@
#-------end PMA++_obj---------------------------



#-------begin TOOLS_obj----------------------------------
$(OBJDIR)/readtools.o: $(CORE_obj) src/tools/readtools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/tools/readtools.f90 -o $@

$(OBJDIR)/printtools.o: $(CORE_obj) src/tools/printtools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/tools/printtools.f90 -o $@
	
$(OBJDIR)/geom_tools.o: $(CORE_obj) $(MATHTOOLS_obj) $(OBJDIR)/core_tools.o $(OBJDIR)/readtools.o src/tools/geom_tools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/tools/geom_tools.f90 -o $@

$(OBJDIR)/simegen.o: $(CORE_obj) $(OBJDIR)/core_tools.o $(OBJDIR)/geom_tools.o src/tools/simegen.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/tools/simegen.f90 -o $@

$(OBJDIR)/read_inputs.o: $(OBJDIR)/simegen.o $(OBJDIR)/objfnc.o $(CORE_obj) $(OBJDIR)/readtools.o src/tools/read_inputs.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/tools/read_inputs.f90 -o $@

$(OBJDIR)/drutes_init.o: $(OBJDIR)/read_inputs.o $(OBJDIR)/readtools.o $(OBJDIR)/core_tools.o $(CORE_obj) src/tools/drutes_init.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/tools/drutes_init.f90 -o $@

$(OBJDIR)/postpro.o: $(CORE_obj) $(MATHTOOLS_obj) $(OBJDIR)/geom_tools.o src/tools/postpro.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/tools/postpro.f90 -o $@

$(OBJDIR)/objfnc.o: $(CORE_obj) $(OBJDIR)/readtools.o src/tools/objfnc.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/tools/objfnc.f90 -o $@
	
$(OBJDIR)/datetime.o: $(CORE_obj)  src/tools/datetime.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/tools/datetime.f90 -o $@
#-------end TOOLS_obj------------------------------------



#-------begin RE_obj--------------------------------
$(OBJDIR)/re_globals.o: $(CORE_obj) src/models/RE/re_globals.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE/re_globals.f90 -o $@

$(OBJDIR)/re_constitutive.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/re_globals.o src/models/RE/re_constitutive.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE/re_constitutive.f90 -o $@

$(OBJDIR)/re_total.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/re_globals.o $(OBJDIR)/re_constitutive.o src/models/RE/re_total.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE/re_total.f90 -o $@

$(OBJDIR)/re_reader.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/re_globals.o src/models/RE/re_reader.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE/re_reader.f90 -o $@

$(OBJDIR)/re_pointers.o: $(CORE_obj) $(OBJDIR)/re_globals.o $(OBJDIR)/re_constitutive.o $(OBJDIR)/re_total.o $(OBJDIR)/re_reader.o $(OBJDIR)/re_evap_methods.o src/models/RE/re_pointers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE/re_pointers.f90 -o $@
	
$(OBJDIR)/re_analytical.o: $(CORE_obj) $(OBJDIR)/re_globals.o $(OBJDIR)/re_constitutive.o src/models/RE/re_analytical.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE/re_analytical.f90 -o $@

$(OBJDIR)/re_evap_methods.o: $(CORE_obj) $(OBJDIR)/re_globals.o $(OBJDIR)/re_constitutive.o src/models/RE/re_evap_methods.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE/re_evap_methods.f90 -o $@
	
#-------end RE_obj--------------------------------

#------begin HEAT_obj -----------------------------------
$(OBJDIR)/heat_globals.o: $(CORE_obj) src/models/heat/heat_globals.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/heat/heat_globals.f90 -o $@

$(OBJDIR)/heat_fnc.o: $(CORE_obj) $(OBJDIR)/heat_globals.o src/models/heat/heat_fnc.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/heat/heat_fnc.f90 -o $@

$(OBJDIR)/heat_reader.o: $(CORE_obj) $(OBJDIR)/heat_globals.o $(OBJDIR)/heat_fnc.o src/models/heat/heat_reader.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/heat/heat_reader.f90 -o $@

$(OBJDIR)/heat_pointers.o: $(CORE_obj) $(RE_obj) $(OBJDIR)/heat_globals.o $(OBJDIR)/heat_fnc.o $(OBJDIR)/heat_reader.o src/models/heat/heat_pointers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/heat/heat_pointers.f90 -o $@
#------end HEAT_obj-------------------------------------


#------begin frozen_obj -----------------------------------
$(OBJDIR)/freeze_globs.o: $(CORE_obj) src/models/soilfreeze/freeze_globs.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/soilfreeze/freeze_globs.f90 -o $@

$(OBJDIR)/freeze_helper.o: $(CORE_obj) $(RE_obj) $(OBJDIR)/freeze_globs.o src/models/soilfreeze/freeze_helper.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/soilfreeze/freeze_helper.f90 -o $@

$(OBJDIR)/freeze_fnc.o: $(CORE_obj) $(OBJDIR)/freeze_helper.o $(OBJDIR)/freeze_globs.o src/models/soilfreeze/freeze_fnc.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/soilfreeze/freeze_fnc.f90 -o $@

$(OBJDIR)/freeze_reader.o: $(CORE_obj) $(OBJDIR)/freeze_globs.o src/models/soilfreeze/freeze_reader.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/soilfreeze/freeze_reader.f90 -o $@

$(OBJDIR)/freeze_pointers.o: $(CORE_obj) $(RE_obj) $(HEAT_obj) $(OBJDIR)/freeze_globs.o $(OBJDIR)/freeze_reader.o src/models/soilfreeze/freeze_pointers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/soilfreeze/freeze_pointers.f90 -o $@
#------end frozen_obj -----------------------------------

#-------begin ADE_obj-------------------------------
$(OBJDIR)/ADE_globals.o: $(CORE_obj) src/models/ADE/ADE_globals.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/ADE/ADE_globals.f90 -o $@

$(OBJDIR)/ADE_fnc.o: $(CORE_obj) $(OBJDIR)/ADE_globals.o src/models/ADE/ADE_fnc.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/ADE/ADE_fnc.f90 -o $@

$(OBJDIR)/ADE_reader.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/ADE_globals.o src/models/ADE/ADE_reader.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/ADE/ADE_reader.f90 -o $@

$(OBJDIR)/ADE_pointers.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/ADE_globals.o $(OBJDIR)/ADE_reader.o $(RE_obj) src/models/ADE/ADE_pointers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/ADE/ADE_pointers.f90 -o $@
#------end ADE_obj---------------------------------



#-------begin REDUAL_obj-----------------------------
$(OBJDIR)/Re_dual_globals.o: $(CORE_obj) src/models/RE_dual/Re_dual_globals.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE_dual/Re_dual_globals.f90 -o $@

$(OBJDIR)/Re_dual_reader.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/Re_dual_globals.o src/models/RE_dual/Re_dual_reader.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE_dual/Re_dual_reader.f90 -o $@

$(OBJDIR)/Re_dual_totH.o: $(CORE_obj) $(TOOLS_obj) $(RE_obj) $(OBJDIR)/Re_dual_globals.o $(OBJDIR)/Re_dual_reader.o src/models/RE_dual/Re_dual_totH.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE_dual/Re_dual_totH.f90 -o $@

$(OBJDIR)/Re_dual_coupling.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/Re_dual_globals.o $(OBJDIR)/Re_dual_reader.o $(OBJDIR)/Re_dual_totH.o src/models/RE_dual/Re_dual_coupling.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE_dual/Re_dual_coupling.f90 -o $@

$(OBJDIR)/Re_dual_tab.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/Re_dual_globals.o $(OBJDIR)/Re_dual_reader.o $(OBJDIR)/Re_dual_totH.o $(OBJDIR)/Re_dual_coupling.o src/models/RE_dual/Re_dual_tab.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE_dual/Re_dual_tab.f90 -o $@

$(OBJDIR)/Re_dual_bc.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/Re_dual_globals.o src/models/RE_dual/Re_dual_bc.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE_dual/Re_dual_bc.f90 -o $@

$(OBJDIR)/Re_dual_pointers.o: $(CORE_obj) $(RE_obj) $(OBJDIR)/Re_dual_reader.o $(OBJDIR)/Re_dual_totH.o $(OBJDIR)/Re_dual_tab.o $(OBJDIR)/Re_dual_bc.o src/models/RE_dual/Re_dual_pointers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/RE_dual/Re_dual_pointers.f90 -o $@
#-------end REDUAL_obj-------------------------------


#-------begin BOUSSINESQ-----------------------------
$(OBJDIR)/boussglob.o: $(CORE_obj) $(TOOLS_obj) src/models/boussinesq/boussglob.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/boussinesq/boussglob.f90 -o $@

$(OBJDIR)/boussread.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/boussglob.o src/models/boussinesq/boussread.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/boussinesq/boussread.f90 -o $@

$(OBJDIR)/boussfnc.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/boussglob.o src/models/boussinesq/boussfnc.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/boussinesq/boussfnc.f90 -o $@

$(OBJDIR)/bousspointers.o: $(CORE_obj) $(OBJDIR)/boussfnc.o $(OBJDIR)/boussglob.o $(OBJDIR)/boussread.o src/models/boussinesq/bousspointers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/boussinesq/bousspointers.f90 -o $@
#-------end BOUSSINESQ-------------------------------

#------begin FEMTOOLS_obj-----------------------------
$(OBJDIR)/fem_tools.o: $(CORE_obj) $(MATHTOOLS_obj) $(TOOLS_obj) $(PMAoo_obj) src/femtools/fem_tools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/femtools/fem_tools.f90 -o $@

$(OBJDIR)/feminittools.o: $(CORE_obj) $(MATHTOOLS_obj) $(TOOLS_obj) $(RE_obj) $(PMAoo_obj) src/femtools/feminittools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/femtools/feminittools.f90 -o $@

$(OBJDIR)/capmat.o: $(CORE_obj) src/femtools/capmat.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/femtools/capmat.f90 -o $@

$(OBJDIR)/stiffmat.o: $(CORE_obj) $(MATHTOOLS_obj) $(OBJDIR)/fem_tools.o src/femtools/stiffmat.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/femtools/stiffmat.f90 -o $@

$(OBJDIR)/femmat.o: $(CORE_obj) $(PMAoo_obj) $(OBJDIR)/fem_tools.o $(OBJDIR)/stiffmat.o $(OBJDIR)/capmat.o $(OBJDIR)/decomp_vars.o src/femtools/femmat.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/femtools/femmat.f90 -o $@

$(OBJDIR)/fem.o: $(CORE_obj) $(MATHTOOLS_obj) $(DECOMPO_obj) $(TOOLS_obj) $(OBJDIR)/femmat.o src/femtools/fem.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/femtools/fem.f90 -o $@
#------end FEMTOOLS_obj------------------------------


#------begin KINWAVE_obj-----------------------------
$(OBJDIR)/kinglobs.o: $(CORE_obj) src/models/kinwave/kinglobs.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/kinwave/kinglobs.f90 -o $@

$(OBJDIR)/kinfnc.o: $(CORE_obj) $(OBJDIR)/kinglobs.o src/models/kinwave/kinfnc.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/kinwave/kinfnc.f90 -o $@

$(OBJDIR)/kinreader.o: $(CORE_obj) $(OBJDIR)/kinglobs.o src/models/kinwave/kinreader.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/kinwave/kinreader.f90 -o $@

$(OBJDIR)/kinpointer.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/kinglobs.o $(OBJDIR)/kinreader.o src/models/kinwave/kinpointer.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/kinwave/kinpointer.f90 -o $@
#------end KINWAVE_obj-------------------------------


#------begin evaporation_obj-------------------------
$(OBJDIR)/evapglob.o: $(CORE_obj) src/models/REevap/evapglob.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/REevap/evapglob.f90 -o $@

$(OBJDIR)/evapreader.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/evapglob.o src/models/REevap/evapreader.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/REevap/evapreader.f90 -o $@

$(OBJDIR)/evap_RE_constitutive.o: $(CORE_obj) $(RE_obj) $(OBJDIR)/evapglob.o src/models/REevap/evap_RE_constitutive.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/REevap/evap_RE_constitutive.f90 -o $@

$(OBJDIR)/evap_heat_constitutive.o: $(CORE_obj) $(HEAT_obj) $(OBJDIR)/evap_RE_constitutive.o src/models/REevap/evap_heat_constitutive.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/REevap/evap_heat_constitutive.f90 -o $@

$(OBJDIR)/evappointers.o: $(CORE_obj) $(HEAT_obj) $(OBJDIR)/evapbc4heat.o $(OBJDIR)/evapreader.o $(OBJDIR)/evapglob.o $(OBJDIR)/evap_RE_constitutive.o $(OBJDIR)/evap_heat_constitutive.o src/models/REevap/evappointers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/REevap/evappointers.f90 -o $@

$(OBJDIR)/evapbc4heat.o: $(CORE_obj) $(RE_obj) $(OBJDIR)/evap_RE_constitutive.o $(OBJDIR)/evap_heat_constitutive.o src/models/REevap/evapbc4heat.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/REevap/evapbc4heat.f90 -o $@
#------end evaporation_obj-------------------------


#------begin netcdf_obj----------------------------

$(OBJDIR)/ncglobvars.o: $(CORE_obj) $(TOOLS_obj) src/models/fluxLS/ncglobvars.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/ncglobvars.f90 -o $@


$(OBJDIR)/nctools.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/ncglobvars.o src/models/fluxLS/nctools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/nctools.f90 -o $@


$(OBJDIR)/ncdem.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/nctools.o $(OBJDIR)/ncglobvars.o src/models/fluxLS/ncdem.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/ncdem.f90 -o $@


$(OBJDIR)/ncmesh.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/nctools.o $(OBJDIR)/ncglobvars.o $(OBJDIR)/ncdem.o src/models/fluxLS/ncmesh.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/ncmesh.f90 -o $@


$(OBJDIR)/ncmap.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/nctools.o $(OBJDIR)/ncglobvars.o $(OBJDIR)/ncmesh.o src/models/fluxLS/ncmap.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/ncmap.f90 -o $@


$(OBJDIR)/netcdfflux.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/ncglobvars.o src/models/fluxLS/netcdfflux.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/netcdfflux.f90 -o $@


$(OBJDIR)/ncfluxarea.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/ncglobvars.o src/models/fluxLS/ncfluxarea.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/ncfluxarea.f90 -o $@


$(OBJDIR)/lsconstitutive.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/ncglobvars.o $(OBJDIR)/netcdfflux.o src/models/fluxLS/lsconstitutive.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/lsconstitutive.f90 -o $@


$(OBJDIR)/init_netcdf.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/nctools.o $(OBJDIR)/netcdfflux.o $(OBJDIR)/ncglobvars.o $(OBJDIR)/ncdem.o $(OBJDIR)/ncmesh.o $(OBJDIR)/ncmap.o $(OBJDIR)/ncfluxarea.o src/models/fluxLS/init_netcdf.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/init_netcdf.f90 -o $@


$(OBJDIR)/ncpointers.o: $(CORE_obj) $(TOOLS_obj) $(OBJDIR)/ncglobvars.o $(OBJDIR)/init_netcdf.o $(OBJDIR)/lsconstitutive.o src/models/fluxLS/ncpointers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/models/fluxLS/ncpointers.f90 -o $@



#-------begin POINTERS_obj--------------------------------
$(OBJDIR)/manage_pointers.o: $(CORE_obj) $(TOOLS_obj) $(FEMTOOLS_obj) $(MATHTOOLS_obj) $(DECOMPO_obj) $(MODEL_objs) src/pointerman/manage_pointers.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -cpp -c src/pointerman/manage_pointers.f90 -o $@
#-------end pointers_obj--------------------------------


#-------begin DECOMPO_obj--------------------------------
$(OBJDIR)/decomp_vars.o: $(PMAoo_obj) src/decompo/decomp_vars.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/decompo/decomp_vars.f90 -o $@

$(OBJDIR)/decomposer.o: $(CORE_obj) $(TOOLS_obj) $(PMAoo_obj) $(OBJDIR)/decomp_vars.o $(OBJDIR)/decomp_tools.o src/decompo/decomposer.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/decompo/decomposer.f90 -o $@

$(OBJDIR)/decomp_tools.o: $(CORE_obj) $(MATHTOOLS_obj) $(OBJDIR)/decomp_vars.o src/decompo/decomp_tools.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/decompo/decomp_tools.f90 -o $@

$(OBJDIR)/schwarz_dd.o: $(CORE_obj) $(MATHTOOLS_obj) $(OBJDIR)/femmat.o $(OBJDIR)/decomp_vars.o $(OBJDIR)/decomposer.o $(OBJDIR)/decomp_tools.o src/decompo/schwarz_dd.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/decompo/schwarz_dd.f90 -o $@

$(OBJDIR)/schwarz_dd2subcyc.o: $(CORE_obj) $(MATHTOOLS_obj) $(OBJDIR)/femmat.o $(OBJDIR)/decomp_vars.o $(OBJDIR)/decomposer.o $(OBJDIR)/decomp_tools.o src/decompo/schwarz_dd2subcyc.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/decompo/schwarz_dd2subcyc.f90 -o $@
#-------end DECOMPO_obj--------------------------------



#----build main---------
$(OBJDIR)/main.o: $(ALL_objs) src/core/main.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c src/core/main.f90 -o $@
#-----------------------

cleanall:
	rm -rf ./build bin/*
	
clean:
	rm -rf ./build
	
git:
	cat /etc/hostname > sync.stamp && date >> sync.stamp & rm -rf *.o *.mod bin/* && git commit -a

push: 
	git push

tar :
	 tar -czf $d.tgz src Makefile drutes.conf 


