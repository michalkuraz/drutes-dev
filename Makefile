##compiler options

in := cd objs
out := cd ..

#options for debugging, use for development  
c= gfortran -fimplicit-none  -fcoarray=single -fbounds-check -fbacktrace -g -g3 -fdefault-real-8 -O0 -finit-real=nan -Wsurprising

#options for optimized compilation, use for production purposes on well debugged versions
#c=gfortran -fimplicit-none  -fcoarray=single -fdefault-real-8 -O3 -finit-real=nan -ffpe-summary=none -fno-backtrace  

d=drutes_obj-`date -I`

all : main.o $(ALL_objs)
	    [ -d bin ] || mkdir bin && $c -g -o bin/drutes main.o $(ALL_objs)
install : main.o $(ALL_objs)
	    [ -d bin ] || mkdir bin && $c -g -o bin/drutes main.o $(ALL_objs) && [ -d obj ] || mkdir obj && mv *.o *.mod obj
	
dir="obj"

servers="miguel@neptun01.fsv.cvut.cz:~  miguel@matsrv-lin01.fsv.cvut.cz:~ miguel@cml.fsv.cvut.cz:~"


#----------------objects definitions-------------------------------
CORE_obj := typy.o global_objs.o globals.o globals1D.o globals2D.o  debug_tools.o core_tools.o pde_objs.o dummy_procs.o global4solver.o
POINTERMAN_obj := manage_pointers.o
RE_obj := re_constitutive.o re_reader.o re_globals.o re_total.o re_pointers.o re_analytical.o re_evap_methods.o
MATHTOOLS_obj :=  linalg.o integral.o solver_interfaces.o simplelinalg.o gmres_solver.o
TOOLS_obj := printtools.o simegen.o read_inputs.o drutes_init.o geom_tools.o postpro.o readtools.o objfnc.o
FEMTOOLS_obj := feminittools.o capmat.o stiffmat.o fem.o fem_tools.o femmat.o
DECOMPO_obj :=  decomp_tools.o schwarz_dd.o  decomp_vars.o decomposer.o schwarz_dd2subcyc.o
PMAoo_obj := fullmatrix.o mtx.o mtx_int.o mtxiotools.o pmatools.o solvers.o sparsematrix.o sparsematrix_int.o matmod.o reorder.o
BOUSSINESQ_obj := boussglob.o boussread.o boussfnc.o bousspointers.o
ADE_obj := ADE_fnc.o ADE_reader.o ADE_globals.o ADE_pointers.o
REDUAL_obj := Re_dual_totH.o Re_dual_globals.o Re_dual_pointers.o Re_dual_reader.o Re_dual_tab.o Re_dual_coupling.o Re_dual_bc.o
HEAT_obj := heat_fnc.o heat_pointers.o heat_globals.o heat_reader.o
KINWAVE_obj := kinreader.o kinglobs.o kinfnc.o kinpointer.o
FROZEN_obj := freeze_globs.o freeze_helper.o freeze_fnc.o freeze_reader.o freeze_pointers.o 
REevap_obj :=  evapglob.o evappointers.o evap_RE_constitutive.o evap_heat_constitutive.o evapreader.o evapbc4heat.o

MODEL_objs := $(RE_obj)  $(BOUSSINESQ_obj) $(ADE_obj) $(REDUAL_obj)  $(HEAT_obj) $(LTNE_obj) $(FROZEN_obj) $(KINWAVE_obj) $(REevap_obj)

ALL_objs := $(CORE_obj) $(TOOLS_obj) $(POINTERMAN_obj) $(MATHTOOLS_obj) $(FEMTOOLS_obj) $(DECOMPO_obj)  $(PMAoo_obj) $(MODEL_objs)
#-----------------------------------------------------------------

#-------begin CORE_obj--------------------------------
typy.o: src/core/typy.f90
	$c -c src/core/typy.f90 
global_objs.o: typy.o $(PMAoo_obj) src/core/global_objs.f90
	$c -c src/core/global_objs.f90
global4solver.o: typy.o src/core/global4solver.f90
	$c -c src/core/global4solver.f90
pde_objs.o: typy.o global_objs.o $(PMAoo_obj) globals.o decomp_vars.o  src/core/pde_objs.f90
	$c -c src/core/pde_objs.f90
globals.o: typy.o global_objs.o src/core/globals.f90
	$c -c src/core/globals.f90
globals1D.o: typy.o global_objs.o src/core/globals1D.f90
	$c -c src/core/globals1D.f90
globals2D.o: typy.o global_objs.o src/core/globals2D.f90
	$c -c src/core/globals2D.f90
core_tools.o: typy.o global_objs.o globals.o   src/core/core_tools.f90
	$c -c src/core/core_tools.f90
dummy_procs.o: typy.o global_objs.o globals.o pde_objs.o src/core/dummy_procs.f90
	$c -c src/core/dummy_procs.f90
debug_tools.o: typy.o core_tools.o src/core/debug_tools.f90
	$c -c src/core/debug_tools.f90
#---------end CORE_obj------------------------------


#------begin MATHTOOLS_obj-----------------------------
linalg.o: $(CORE_obj) src/mathtools/linalg.f90
	$c -c src/mathtools/linalg.f90
gmres_solver.o: $(CORE_obj) src/mathtools/gmres_solver.f90
	$c -c src/mathtools/gmres_solver.f90
integral.o: $(CORE_obj) linalg.o src/mathtools/integral.f90
	$c -c src/mathtools/integral.f90
simplelinalg.o:  $(CORE_obj) $(PMAoo_obj) re_globals.o linalg.o src/mathtools/simplelinalg.f90
	$c -c src/mathtools/simplelinalg.f90
solver_interfaces.o:  $(CORE_obj) $(PMAoo_obj) readtools.o simplelinalg.o gmres_solver.o src/mathtools/solver_interfaces.f90
	$c -c src/mathtools/solver_interfaces.f90

#------end MATHTOOLS_obj---------------------------------


#--------begin PMAoo_obj------------------------
pmatools.o: typy.o src/pma++/pmatools.f90
	$c -c  src/pma++/pmatools.f90
mtx.o: typy.o pmatools.o  src/pma++/mtx.f90
	$c -c src/pma++/mtx.f90
mtx_int.o: typy.o pmatools.o  src/pma++/mtx_int.f90
	$c -c src/pma++/mtx_int.f90	
mtxiotools.o: typy.o src/pma++/mtxiotools.f90
	$c -c src/pma++/mtxiotools.f90 
fullmatrix.o: typy.o mtx.o src/pma++/fullmatrix.f90
	$c -c src/pma++/fullmatrix.f90
sparsematrix.o: typy.o mtx.o src/pma++/sparsematrix.f90
	$c -c src/pma++/sparsematrix.f90
sparsematrix_int.o: typy.o mtx.o src/pma++/sparsematrix_int.f90
	$c -c src/pma++/sparsematrix_int.f90	
solvers.o: global4solver.o  typy.o mtx.o src/pma++/solvers.f90
	$c -c src/pma++/solvers.f90
matmod.o: typy.o mtx.o src/pma++/matmod.f90
	$c -c  src/pma++/matmod.f90
datasetup.o: typy.o mtx.o src/pma++/datasetup.f90
	$c -c src/pma++/datasetup.f90
reorder.o: typy.o mtx.o datasetup.o solvers.o src/pma++/reorder.f90
	$c -c  src/pma++/reorder.f90
#-------end PMA++_obj---------------------------

#-------begin TOOLS_obj----------------------------------
readtools.o: $(CORE_obj) src/tools/readtools.f90
	$c -c src/tools/readtools.f90
printtools.o: $(CORE_obj) src/tools/printtools.f90
	$c -c src/tools/printtools.f90
geom_tools.o: $(CORE_obj) $(MATHTOOLS_obj) core_tools.o readtools.o src/tools/geom_tools.f90
	$c -c src/tools/geom_tools.f90
simegen.o:  $(CORE_obj) core_tools.o geom_tools.o src/tools/simegen.f90
	$c -c src/tools/simegen.f90
read_inputs.o:  simegen.o objfnc.o $(CORE_obj) readtools.o src/tools/read_inputs.f90
	$c -c src/tools/read_inputs.f90
drutes_init.o: read_inputs.o readtools.o core_tools.o $(CORE_obj) src/tools/drutes_init.f90
	$c -c src/tools/drutes_init.f90
postpro.o: $(CORE_obj) $(MATHTOOLS_obj) geom_tools.o src/tools/postpro.f90
	$c -c src/tools/postpro.f90
objfnc.o: $(CORE_obj) readtools.o src/tools/objfnc.f90
	$c -c  src/tools/objfnc.f90
#-------end TOOLS_obj------------------------------------

#-------begin RE_obj--------------------------------
re_globals.o: $(CORE_obj) src/models/RE/re_globals.f90
	$c -c  src/models/RE/re_globals.f90
re_constitutive.o: $(CORE_obj) $(TOOLS_obj) re_globals.o src/models/RE/re_constitutive.f90
	$c -c src/models/RE/re_constitutive.f90
re_total.o: $(CORE_obj) $(TOOLS_obj) re_globals.o re_constitutive.o src/models/RE/re_total.f90
	$c -c src/models/RE/re_total.f90
re_reader.o:  $(CORE_obj) $(TOOLS_obj) re_globals.o  src/models/RE/re_reader.f90
	$c -c src/models/RE/re_reader.f90	
re_pointers.o:  $(CORE_obj) re_globals.o re_constitutive.o re_total.o re_reader.o re_evap_methods.o src/models/RE/re_pointers.f90
	$c -c src/models/RE/re_pointers.f90
re_analytical.o:  $(CORE_obj) re_globals.o re_constitutive.o src/models/RE/re_analytical.f90
	$c -c src/models/RE/re_analytical.f90
re_evap_methods.o: $(CORE_obj) re_globals.o re_constitutive.o src/models/RE/re_evap_methods.f90
	$c -c src/models/RE/re_evap_methods.f90
#-------end CONSTITUTIVE_obj--------------------------------

#------begin HEAT_obj -----------------------------------
heat_globals.o: $(CORE_obj) src/models/heat/heat_globals.f90
	$c -c src/models/heat/heat_globals.f90
heat_fnc.o: $(CORE_obj) heat_globals.o src/models/heat/heat_fnc.f90
	$c -c src/models/heat/heat_fnc.f90
heat_reader.o: $(CORE_obj) heat_globals.o heat_fnc.o src/models/heat/heat_reader.f90
	$c -c src/models/heat/heat_reader.f90
heat_pointers.o: $(CORE_obj) $(RE_obj) heat_globals.o heat_fnc.o heat_reader.o src/models/heat/heat_pointers.f90
	$c -c src/models/heat/heat_pointers.f90
#------end HEAT_obj-------------------------------------


#------begin frozen_obj -----------------------------------
freeze_globs.o: $(CORE_obj) src/models/soilfreeze/freeze_globs.f90
	$c -c src/models/soilfreeze/freeze_globs.f90
freeze_helper.o: $(CORE_obj) $(RE_obj) freeze_globs.o src/models/soilfreeze/freeze_helper.f90
	$c -c src/models/soilfreeze/freeze_helper.f90
freeze_fnc.o: $(CORE_obj) freeze_helper.o freeze_globs.o src/models/soilfreeze/freeze_fnc.f90
	$c -c src/models/soilfreeze/freeze_fnc.f90
freeze_reader.o: $(CORE_obj) freeze_globs.o src/models/soilfreeze/freeze_reader.f90
	$c -c src/models/soilfreeze/freeze_reader.f90	
freeze_pointers.o: $(CORE_obj) $(RE_obj) $(HEAT_obj) freeze_globs.o freeze_reader.o src/models/soilfreeze/freeze_pointers.f90
	$c -c src/models/soilfreeze/freeze_pointers.f90
#------end frozen_obj -----------------------------------

#-------begin ADE_obj-------------------------------
ADE_globals.o: $(CORE_obj) src/models/ADE/ADE_globals.f90
	$c -c src/models/ADE/ADE_globals.f90
ADE_fnc.o: $(CORE_obj) ADE_globals.o src/models/ADE/ADE_fnc.f90
	$c -c src/models/ADE/ADE_fnc.f90
ADE_reader.o: $(CORE_obj) $(TOOLS_obj) ADE_globals.o src/models/ADE/ADE_reader.f90
	$c -c src/models/ADE/ADE_reader.f90
ADE_pointers.o: $(CORE_obj) $(TOOLS_obj) ADE_globals.o  ADE_reader.o  $(RE_obj) src/models/ADE/ADE_pointers.f90
	$c -c src/models/ADE/ADE_pointers.f90
#------end ADE_obj---------------------------------


#-------begin REDUAL_obj-----------------------------
Re_dual_globals.o: $(CORE_obj) src/models/RE_dual/Re_dual_globals.f90
	$c -c src/models/RE_dual/Re_dual_globals.f90
Re_dual_reader.o: $(CORE_obj) $(TOOLS_obj) Re_dual_globals.o src/models/RE_dual/Re_dual_reader.f90
	$c -c src/models/RE_dual/Re_dual_reader.f90
Re_dual_totH.o: $(CORE_obj) $(TOOLS_obj) $(RE_obj) Re_dual_globals.o Re_dual_reader.o src/models/RE_dual/Re_dual_totH.f90
	$c -c src/models/RE_dual/Re_dual_totH.f90
Re_dual_coupling.o: $(CORE_obj) $(TOOLS_obj) Re_dual_globals.o Re_dual_reader.o Re_dual_totH.o src/models/RE_dual/Re_dual_coupling.f90
	$c -c src/models/RE_dual/Re_dual_coupling.f90
Re_dual_tab.o: $(CORE_obj) $(TOOLS_obj) Re_dual_globals.o Re_dual_reader.o Re_dual_totH.o Re_dual_coupling.o src/models/RE_dual/Re_dual_tab.f90
	$c -c src/models/RE_dual/Re_dual_tab.f90	
Re_dual_bc.o: $(CORE_obj) $(TOOLS_obj) Re_dual_globals.o src/models/RE_dual/Re_dual_bc.f90
	$c -c src/models/RE_dual/Re_dual_bc.f90
Re_dual_pointers.o: $(CORE_obj) $(RE_obj) Re_dual_reader.o Re_dual_totH.o Re_dual_tab.o Re_dual_bc.o src/models/RE_dual/Re_dual_pointers.f90
	$c -c src/models/RE_dual/Re_dual_pointers.f90
#-------end REDUAL_obj-------------------------------


#-------begin BOUSSINESQ-----------------------------
boussglob.o:  $(CORE_obj) $(TOOLS_obj) src/models/boussinesq/boussglob.f90
	$c -c src/models/boussinesq/boussglob.f90
boussread.o: $(CORE_obj) $(TOOLS_obj) boussglob.o src/models/boussinesq/boussread.f90
	$c -c src/models/boussinesq/boussread.f90
boussfnc.o:  $(CORE_obj) $(TOOLS_obj) boussglob.o src/models/boussinesq/boussfnc.f90
	$c -c src/models/boussinesq/boussfnc.f90
bousspointers.o: $(CORE_obj)  boussfnc.o boussglob.o boussread.o src/models/boussinesq/bousspointers.f90
	$c -c src/models/boussinesq/bousspointers.f90
#-------end BOUSSINESQ-------------------------------



#------begin FEMTOOLS_obj-----------------------------
fem_tools.o:  $(CORE_obj) $(MATHTOOLS_obj) $(TOOLS_obj) $(PMA++_obj) src/femtools/fem_tools.f90
	$c -c  src/femtools/fem_tools.f90
feminittools.o:  $(CORE_obj) $(MATHTOOLS_obj) $(TOOLS_obj) $(RE_obj) $(PMAoo_obj) src/femtools/feminittools.f90
	$c -c src/femtools/feminittools.f90
capmat.o:$(CORE_obj) src/femtools/capmat.f90
	$c -c src/femtools/capmat.f90
stiffmat.o: $(CORE_obj) $(LINALG_obj) fem_tools.o src/femtools/stiffmat.f90
	$c -c src/femtools/stiffmat.f90
femmat.o: $(CORE_obj) $(PMA++_obj) fem_tools.o stiffmat.o capmat.o decomp_vars.o src/femtools/femmat.f90
	$c -c src/femtools/femmat.f90
fem.o: $(CORE_obj) $(LINALG_obj) $(DECOMPO_obj) $(TOOLS_obj) femmat.o src/femtools/fem.f90
	$c -c src/femtools/fem.f90
#------end FEMTOOLS_obj------------------------------


#------begin KINWAVE_obj-----------------------------
kinglobs.o: $(CORE_obj) src/models/kinwave/kinglobs.f90
	$c -c src/models/kinwave/kinglobs.f90
kinfnc.o: $(CORE_obj) kinglobs.o src/models/kinwave/kinfnc.f90
	$c -c src/models/kinwave/kinfnc.f90
kinreader.o: $(CORE_obj) kinglobs.o src/models/kinwave/kinreader.f90
	$c -c src/models/kinwave/kinreader.f90
kinpointer.o: $(CORE_obj) $(TOOLS_obj) kinglobs.o kinreader.o src/models/kinwave/kinpointer.f90
	$c -c src/models/kinwave/kinpointer.f90
#------end KINWAVE_obj-------------------------------




#------begin evaporation_obj-------------------------
evapglob.o: $(CORE_obj) src/models/REevap/evapglob.f90
	$c -c src/models/REevap/evapglob.f90
evapreader.o: $(CORE_obj) $(TOOLS_obj) evapglob.o src/models/REevap/evapreader.f90
	$c -c src/models/REevap/evapreader.f90
evap_RE_constitutive.o: $(CORE_obj) $(RE_obj) evapglob.o src/models/REevap/evap_RE_constitutive.f90
	$c -c src/models/REevap/evap_RE_constitutive.f90
evap_heat_constitutive.o: $(CORE_obj) $(HEAT_obj) evap_RE_constitutive.o src/models/REevap/evap_heat_constitutive.f90
	$c -c src/models/REevap/evap_heat_constitutive.f90	
evappointers.o: $(CORE_obj) $(HEAT_obj) evapbc4heat.o evapreader.o evapglob.o evap_RE_constitutive.o  evap_heat_constitutive.o  src/models/REevap/evappointers.f90
	$c -c src/models/REevap/evappointers.f90
evapbc4heat.o: $(CORE_obj) $(RE_obj) evap_RE_constitutive.o evap_heat_constitutive.o src/models/REevap/evapbc4heat.f90
	$c -c src/models/REevap/evapbc4heat.f90
#------end evaporation_obj-------------------------




#-------begin POINTERS_obj--------------------------------
manage_pointers.o: $(CORE_obj) $(TOOLS_obj) $(CORE_obj) $(FEMTOOLS_obj) $(LINALG_obj) $(DECOMPO_obj) $(MODEL_objs) $(EVAPORATION_obj) src/pointerman/manage_pointers.f90 
	$c -c src/pointerman/manage_pointers.f90
#-------end pointers_obj--------------------------------


#-------begin DECOMPO_obj--------------------------------
decomp_vars.o:  $(PMAoo_obj) src/decompo/decomp_vars.f90
	$c -c src/decompo/decomp_vars.f90
decomposer.o: $(CORE_obj) $(TOOLS_obj) $(PMAoo_obj) decomp_vars.o decomp_tools.o src/decompo/decomposer.f90
	$c -c src/decompo/decomposer.f90
decomp_tools.o: $(CORE_obj) $(MATHTOOLS_obj) decomp_vars.o  src/decompo/decomp_tools.f90
	$c -c src/decompo/decomp_tools.f90
schwarz_dd.o:  $(CORE_obj) $(MATHTOOLS_obj)  femmat.o decomp_vars.o decomposer.o decomp_tools.o src/decompo/schwarz_dd.f90
	$c -c src/decompo/schwarz_dd.f90
schwarz_dd2subcyc.o: $(CORE_obj) $(MATHTOOLS_obj)  femmat.o decomp_vars.o decomposer.o decomp_tools.o src/decompo/schwarz_dd2subcyc.f90
	$c -c src/decompo/schwarz_dd2subcyc.f90
#-------end DECOMPO_obj--------------------------------






#----build main---------
main.o:  $(ALL_objs) src/core/main.f90
	$c -c src/core/main.f90 
#-----------------------


cleanall:
	rm -rf *.o *.mod bin/*
	
clean:
	rm -rf *.o *.mod
	
git:
	cat /etc/hostname > sync.stamp && date >> sync.stamp & rm -rf *.o *.mod bin/* && git commit -a

push: 
	git push

syncup:
	cat /etc/hostname > sync.stamp && date >> sync.stamp && rsync -avztu -e ssh --delete --exclude 'out' --exclude '*.o' --exclude '*.mod' --exclude 'bin'  --exclude '*~' --exclude '*attr' --exclude '.git' ./ miguel@cml.fsv.cvut.cz:~/drutes-obj/

syncdown:
	tar -czf /tmp/git.tgz .git && rsync -avztu -e ssh --delete --exclude 'out/*' miguel@cml.fsv.cvut.cz:~/drutes-obj/ ./ && echo "last sync:" && cat sync.stamp && rm -rf .git && tar -xzf /tmp/git.tgz
		
pull-dd:
	git checkout dev-dd && git pull https://www.github.com/michalkuraz/drutes-dev dev-dd



save:
	 tar -czf $d.tgz src Makefile drutes.conf  ; for i in `echo $(servers)` ; do scp -P 22  $d.tgz $$i; done

tar :
	 tar -czf $d.tgz src Makefile drutes.conf 


