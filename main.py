from inputs import *
from band_structure import *
from real_structure import *



inputs=inputs()
inputs.run_calculations()
mode=input("Should I color a fermi surface with orbital character or velocity o/[v]?")



degauss=0.1
band_str=band_structure(degauss)

band_str.read_ef_and_dos()
band_str.read_bvec()
band_str.read_all_kpoints()
band_str.read_ene()
band_str.read_qtl()
#band_str.read_almblm()
band_str.which_bands_cross_ef()
band_str.calc_kweights()

#real_str=real_structure(band_str.E_f)


if mode=='o':
 band_str.write_to_qtl_fermisurfer()
else: 
 band_str.calc_fermi_velocity()
 band_str.write_to_fermisurfer()
 
