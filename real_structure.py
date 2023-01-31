from inputs import *

class real_structure(inputs):
 def __init__(self,almblm):
  self.Vtot=[]
  self.Vtot_interstitial=[]
  self.kpoints_interstitial=[]
  self.LM=[]
  self.RADWF=[]
  self.RADWFsmall=[]
  self.RMESH=[]
  self.DR=[]
  self.DVDR=[]
  self.n_r=0
  self.n_at=0
  self.volume=0
 def read_potential(self):
  self.Vtot,self.LM,self.n_r,self.n_at,self.Vtot_interstitial,self.kpoints_interstitial=read_potential(inputs.pot_file)
 def read_radwf(self):
  self.RADWF,self.RADWFsmall,self.RMESH,self.DR=read_radwf(inputs.radwf_file,inputs.CIN)
 def read_volume(self):
  self.volume=read_volume(inputs.scf_file)
 def calc_norm_of_radwf(self):
  self.radwf_norm=calc_norm_of_radwf(self.RADWF,self.DR)
 def calc_overlap(self):
  self.overlap=calc_overlap(self.RADWF,self.RADWFsmall,self.DR)
def read_volume(scf_file):
 print('Read volume...')
 h=open(scf_file)
 tmp=h.readlines()
 h.close()
 for i in tmp: 
  if ':VOL' in i: 
   volume=float(i.split()[-1])
   break
 return volume 

def read_potential(pot_file):
 print('Read intrasitial total potential...')
 h=open(pot_file)
 tmp=h.readlines()
 h.close()
 Vtot=[] #[atom][lm][r]
 LM=[] #[atom][list of lm numbers]
 Vtot_interstitial=[]
 kpoints_interstitial=[]
 for ni,i in enumerate(tmp):
  if 'IN INTERSTITIAL' in i: break
  if 'ATOMNUMBER' in i: 
   Vtot.append([])
   LM.append([])
  if 'VLM(R)' in i:
   Vtot[-1].append([])
   LM[-1].append([int(i.split()[3]),int(i.split()[5])])
  elif len(i.split())==0: continue
  elif len(Vtot)!=0:
    m=i[3:-1]
    for j in range(int(len(m)/19)):
     Vtot[-1][-1].append(float(m[j*19:(j+1)*19]))
 n_at=len(Vtot) #number of atoms
 n_r=len(Vtot[0][0]) #number of r-points
 print (n_at,n_r, LM)

 no_of_vk=0
 

 print (n_at,n_r, LM)
 return Vtot,LM,n_r,n_at,Vtot_interstitial,kpoints_interstitial

def read_radwf(radwf_file,CIN):
 print("CIN",CIN)
 print('Reading radial wave funtions...')
#     IF(MODUS.EQ.'ALM  ') then
#        write(23,4645)jatom,jri(jatom),r0(jatom),dx(jatom),rmt(jatom)
#        do l=0,lmax2
#           write(23,*) l
#           if (l.le.lomax) then
#              write(23,4646) (RRAD1(jrj,l),RRAD2(jrj,l),RADE1(jrj,l),RADE2(jrj,l), &
#                   a1lo(jrj,1,l),b1lo(jrj,1,l),a1lo(jrj,2,l),b1lo(jrj,2,l), &
#                   a1lo(jrj,3,l),b1lo(jrj,3,l),jrj=1,jri(jatom))
#           else
#              write(23,4647) (RRAD1(jrj,l),RRAD2(jrj,l),RADE1(jrj,l),RADE2(jrj,l), &
#                    jrj=1,jri(jatom))
#           endif
#        enddo
#     endif

 h=open(radwf_file,'r')
 tmp=h.readlines()
 h.close()

 RADWF=[]   #RADWF[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-2]-large  component of radwf and udot, and lo if exists
 RADWFsmall=[] #RADWFsmall[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-  small component of radwf and udot, and lo if exists
 mesh_info=[]


 #read radwf
 #factor i^l/sqrt(4pi) comes from the formula from main.F in SRC_lapw7. by RMT we will multiply later
 #! psi(r) = e^ikR Sum(lm) w_lm,a(|r-R-R_a|) Y*(*)_lm(T_a^-1(r-R-R_a))   with
 #!
 #!  w_lm,a(r) = 4pi*i^l [ A_l,m,a *      u_l,a(r,E _l,a) +
 #!                       B_l,m,a * d/dE u_l,a(r,E _l,a) +
 #!                       C_l,m,a *      u_l,a(r,E'_l,a) ] * Rmt_a^2
 #where R_a=pos of atom, R=pos of considered unit cell (here R=0, because we consider only 1 unit cell 
 #--not needed, lapw2 already does it
 n_r,n_at=[],0
 ni=0
 while ni<len(tmp):
  n_at+=1
  n_r.append(int(tmp[ni].split()[1]))
  RADWF.append([])
  RADWFsmall.append([])
  mesh_info.append([float(m) for m in tmp[ni].split()])
  ni+=1
  while ni<len(tmp) and len(tmp[ni].split())==1 :
   l=int(tmp[ni].split()[0])
   RADWF[-1].append([])
   RADWFsmall[-1].append([])
   ni+=1
   for i in tmp[ni:ni+n_r[-1]]:
    try: 
     RADWF[-1][-1].append(np.array([float(i.split()[0]),float(i.split()[2]),float(i.split()[4]),float(i.split()[6]),float(i.split()[8])    ]))  #with RLO
     RADWFsmall[-1][-1].append(np.array([CIN*float(i.split()[1]),CIN*float(i.split()[3]),CIN*float(i.split()[5]),CIN*float(i.split()[7]) ,CIN*float(i.split()[9])    ]))
    except:
     try: 
      RADWF[-1][-1].append(np.array([float(i.split()[0]),float(i.split()[2]),0.,0.,0.]))  #w/o RLO
      RADWFsmall[-1][-1].append(np.array([CIN*float(i.split()[1]),CIN*float(i.split()[3]),0.,0.,0.]))
     except:continue
   ni=ni+n_r[-1]


 for i in range(len(RADWF)):
  print('For atom no. '+str(i)+' I found '+str(len(RADWF[i]))+' wave functions u_l(r) at the r-meshes as follows:')
#  if not np.any([sum(j) for j in
  for j in range(len(RADWF[i])):
   print(len(RADWF[i][j])),
  print(' ')

 print('Making r-mesh...')
 RMESH=[] #RMESH[i][j] i -atom, j -r-point
#(r_i=r0*exp((i-1)*dx)
#   1 781   0.0000100000   0.0159348926   2.5000000000
# write(23,4645)jatom,jri(jatom),r0(jatom),dx(jatom),rmt(jatom)
 for i in mesh_info:
  [at,nr,r0,dx,rmt]=i
  RMESH.append([])
  for j in range(int(nr)):
   RMESH[-1].append(r0*np.exp(j*dx))
  print('RMT='+str(rmt)+' should ='+str(RMESH[-1][-1]))

 RMESH=[ np.array([m for m in i]) for i in RMESH]
 DR=[np.diff(RMESH[at]) for at in range(n_at)] #len=n_r-1 (element 0 is excluded in diff) 
 return RADWF, RADWFsmall, RMESH, DR

def calc_norm_of_radwf(RADWF,DR):
  norm_of_radwf=[[ [ sum( 
   [(RADWF[at][l][r+1][m]**2)*DR[at][r] for r in range(len(RADWF[at][l])-1)]
                     )**0.5 for m in range (len(RADWF[at][l][0])) 
                ] for l in range(len(RADWF[at])) ] for at in range(len(RADWF))]
  return norm_of_radwf

def calc_overlap(RADWF,RADWFsmall, DR):
  norma=calc_norm_of_radwf(RADWF,DR)
  print(norma)
  overlap=[[[[ [ sum( \
   [(RADWF[at][l][r][m]*RADWF[at][l2][r][m2])*DR[at][r] for r in range(len(RADWF[at][l])-1)]\
                     ) if (norma[at][l][m]>0 and norma[at][l2][m2]>0) else 0 for m2 in range (len(RADWF[at][l2][0])) \
                ] for m in range (len(RADWF[at][l][0])) \
                ] for l2 in range(len(RADWF[at])) ]  for l in range(len(RADWF[at])) ] for at in range(len(RADWF))]
  norma=calc_norm_of_radwf(RADWFsmall,DR)
  print(norma)
  overlapsmall=[[[[ [ sum( \
   [(RADWFsmall[at][l][r][m]*RADWFsmall[at][l2][r][m2])*DR[at][r] for r in range(len(RADWFsmall[at][l])-1)]\
                     ) if (norma[at][l][m]>0 and norma[at][l2][m2]>0) else 0 for m2 in range (len(RADWF[at][l2][0])) \
                ] for m in range (len(RADWFsmall[at][l][0])) \
                ] for l2 in range(len(RADWFsmall[at])) ]  for l in range(len(RADWFsmall[at])) ] for at in range(len(RADWFsmall))]
  overlap=np.array(overlap)+np.array(overlapsmall)
  return overlap



