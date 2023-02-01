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
  qtl_file='for_qtl_fs/for_qtl_fs.qtl'
  dos_file=prefix+'.outputt'
  almblm_file=prefix+'.almblm'
  ene_file=prefix+'.energy'
  pot_file=prefix+'.vsp'
  radwf_file=prefix+'.radwf'
  scf_file=prefix+'.scf0'
  so=checkso(prefix) 
  n_l=5 #n_l nieograniczone daje dokladnie te same wyniki co n_l=5; n_l=4 zmienia wynik o 0.01%
  def run_calculations(self):
    run_calc(self.prefix)
 #   os.system("cd hopfield_calc")
 
    


def run_calc(prefix):
 name='for_qtl_fs'
 so=''
 if os.stat(prefix+'.inso').st_size!=0: so=' -so '
 os.system('rm -r '+name+'; save_lapw -d '+name)
 os.system('cp .machines '+name+'/.')
 os.chdir(name)
# os.system('pwd')
# exit()
 os.system('rename_files '+prefix+' '+name)
 
 h=open(name+'.in2')
 tmp=h.readlines()
 h.close()
 line1=tmp[1]
 noe=float(line1.split()[1])
 noe=round(1.5*noe)
 line=line1[:8]+"{:6d}".format(noe)+line1[14:]
 tmp[1]=line
 h=open(name+'.in2','w')
 for l in tmp:
  h.write(l)
 h.close()
 if os.stat(name+'.in2c').st_size!=0:
  h=open(name+'.in2c')
  tmp=h.readlines()
  h.close()
  tmp[1]=line
  h=open(name+'.in2c','w')
  for l in tmp:
   h.write(l)
  h.close()
 os.system("x lapw1 -qtl -p")
 if len(so):  os.system("x lapwso -qtl -p")
 os.system("x lapw2 "+so+" -qtl -p")
 os.chdir('..')

#ins=inputs()
#print ins.prefix
#run_calc(ins)
