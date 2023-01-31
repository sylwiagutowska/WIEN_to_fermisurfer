import numpy as np
import os
import glob
#import hopfield_tetra
#from math import factorial as sil
from scipy.special import sph_harm
from scipy.integrate import quad
from scipy.misc import derivative
from multiprocessing import Process,Pool
Ry_to_eV=13.605693122994
ab_to_ang=0.529177210903

#j_to_ev=1/(1.602176634)*1e19
#hbar=6.582119569e-16 #ev/s
#hbar=hbar/13.6056980659 #ry/s
pi=3.141592653589793238462643
#Na=6.02214076e23
sqrtpm1=(np.pi)**(-0.5)



SOC=0

#if os.stat("DOS/DOS.inso").st_size != 0:  
# CIN=1./(137.0359895*2) #for treating big and small component of wave function
# SOC=1
def checkso(prefix):
  so=''
  if os.stat(prefix+'.inso').st_size!=0: so=' -so '
  return so


class inputs:
  prefix=os.getcwd().split('/')[-1]
  CIN=(1/137.0359895)
  qtl_file=prefix+'.qtl'
  dos_file=prefix+'.outputt'
  almblm_file=prefix+'.almblm'
  ene_file=prefix+'.energy'
  pot_file=prefix+'.vsp'
  radwf_file=prefix+'.radwf'
  scf_file=prefix+'.scf0'
  so=checkso(prefix) 
  n_l=5 #n_l nieograniczone daje dokladnie te same wyniki co n_l=5; n_l=4 zmienia wynik o 0.01%
  nk=48
  def run_calculations():
   yn='n'
#   yn=input('should I run calc? y/[n]: ')
   if yn=='y': 
    run_calc(inputs)




def run_calc(inputs):
 so=''
 if os.stat(inputs.prefix+'.inso').st_size!=0: so=' -so '
 os.system('rm -r hopfield_calc_old; mv hopfield_calc hopfield_calc_old; cp -r run_lapw  hopfield_calc')
 os.system('cp .machines hopfield_calc')
 os.chdir('hopfield_calc')
# os.system('pwd')
# exit()
 os.system('rename_files '+inputs.prefix+' hopfield_calc')
 os.system("x kgen <<'EOF'\n0\n+"+str(inputs.nk)+" "+str(inputs.nk)+" "+str(inputs.nk)+" \n1\nEOF")
 os.system("sed -i 's/NR2V/R2V/' hopfield_calc.in0")
 os.system("sed -i 's/CONT 1/CONT 0/' hopfield_calc.in1")
 os.system("sed -i 's/STOP 1/STOP 0/' hopfield_calc.in1")
 os.system("sed -i 's/1      (GLOBAL/0      (GLOBAL/' hopfield_calc.in1")
 os.system("run_lapw -ec 0.00001 -cc 0.001")
 os.system("rm -r run_lapw "+so+"-p; save_lapw -d run_lapw")
 os.system("x lapw1 -p")
 if len(so): os.system("x lapwso -p")
 os.system("x lapw2 "+so+" -alm -p")
 os.system("x lapw2 "+so+" -qtl -p")
 os.system("x tetra ; x tetra "+so+"-p")
 os.chdir('..')

#ins=inputs()
#print ins.prefix
#run_calc(ins)
