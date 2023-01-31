from inputs import *
from band_structure import *
from real_structure import *




#inputs.run_calculations()
degauss=0.1
band_str=band_structure(degauss)
band_str.read_qtl()

band_str.read_ef_and_dos()
band_str.read_bvec()
band_str.read_all_kpoints()
band_str.read_ene()
band_str.read_almblm()
band_str.which_bands_cross_ef()
band_str.calc_kweights()
#print(len(band_str.ENE[0]))
#band_str.calc_fermi_velocity()

real_str=real_structure(band_str.E_f) 
#real_str.read_volume()
#real_str.read_potential()
#real_str.read_radwf()
#real_str.calc_overlap()
print('real strcture ends')


#band_str.write_to_fermisurfer()
#band_str.write_to_orbital_fermisurfer(real_str.overlap)

band_str.write_to_qtl_fermisurfer()
