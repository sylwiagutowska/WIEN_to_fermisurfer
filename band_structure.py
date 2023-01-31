from inputs import *
from operator import itemgetter
class band_structure(inputs):
 def __init__(self,degauss):
  self.E_f=0
  self.dos=0
  self.n_k=0 #nonequivalents
  self.n_k_total=0 #all
  self.n_at=0
  self.nk3=[0,0,0]
  self.ALMBLM=[]
  self.ENE=[] 
  self.ENE_kweights=[] 
  self.degauss=degauss
  self.all_kpoints_indexes=[]
  self.bvec=[] #reciprocal lattice vectors
  self.nat=0
  self.No_of_bands=[]
  self.No_of_fermi_bands=0
  self.fermi_vec=[]
  self.fermi_vel=[]
  self.qtl=[]
  self.avec=[] #direct lattice becotrs
 def read_ef_and_dos(self):
  self.E_f,self.dos,self.nat=read_ef_and_dos(inputs.dos_file,inputs.prefix)
  self.n_at=self.nat

 def read_ene(self):
  try: self.ENE,self.kweights,self.n_k_total,self.n_k=read_ene(inputs.ene_file+'so',self.n_at,self.degauss,self.E_f)
  except: 
   try: self.ENE,self.kweights,self.n_k_total,self.n_k=read_ene(inputs.ene_file,self.n_at,self.degauss,self.E_f)
   except: self.ENE,self.kweights,self.n_k_total,self.n_k=read_ene(inputs.ene_file+'soup',self.n_at,self.degauss,self.E_f)

 def calc_kweights(self):
  self.degauss=adjust_degauss(self.kweights,self.ENE,self.E_f,self.n_k,self.n_k_total,self.dos)
  print('degauss equal ',self.degauss,' gives the best estimation of dos')
  self.ENE_kweights,self.only_ENE_weights=calc_kweights(self.kweights,self.ENE,self.E_f,self.n_k,self.degauss)
 def read_almblm(self):
  self.ALMBLM,self.No_of_bands,self.n_k,self.n_at=read_almblm(inputs.almblm_file,inputs.n_l,self.n_k)
 def which_bands_cross_ef(self):
  self.ENE,self.ALMBLM,self.qtl,self.No_of_fermi_bands=which_bands_cross_ef(self.ENE,self.ALMBLM,self.qtl,self.E_f) #this ene is transposed from [nk][ib] to [ib][nk] 
  #self.ENE=self.ENE #*13.606
  #self.E_f=self.E_f #*13.606
 def read_all_kpoints(self):
  self.all_kpoints_indexes,self.nk3,self.n_k=read_all_kpoints(inputs.prefix+'.outputkgen', self.bvec)
  self.n_k_total=self.nk3[0]*self.nk3[1]*self.nk3[2]
 def read_bvec(self):
  self.bvec,self.avec=read_bvec(inputs.prefix+'.outputkgen')
 def calc_fermi_velocity(self):
  self.fermi_vel=calc_fermi_velocity(self.ENE,self.all_kpoints_indexes,self.avec,self.nk3,self.No_of_fermi_bands)
 def write_to_fermisurfer(self):
  if len(self.fermi_vel):
   write_to_fermisurfer(self.bvec,self.ENE,self.all_kpoints_indexes,self.E_f,self.fermi_vel)
  else:
   write_to_fermisurfer(self.bvec,self.ENE,self.all_kpoints_indexes,self.E_f)
 def write_to_orbital_fermisurfer(self,RADWF):
  write_to_orbital_fermisurfer(self.bvec,self.ENE,self.all_kpoints_indexes,self.E_f,self.ALMBLM,RADWF,self.only_ENE_weights,inputs.dos_file,self.qtl)
 def write_to_qtl_fermisurfer(self):
  write_to_qtl_fermisurfer(self.bvec,self.ENE,self.all_kpoints_indexes,self.E_f,self.only_ENE_weights,inputs.dos_file,self.qtl)
 def read_qtl(self):
  self.qtl=read_qtl(inputs.qtl_file,inputs.so)

def calc_fermi_velocity(ENE,all_kpoints_index,avec,nk3,No_of_fermi_bands):
    de=[0 for i in range(3)]
    vf=[ [0 for i in range(nk3[0]*nk3[1]*nk3[2])] for ib in range(No_of_fermi_bands)]
    vf_vector=[ [[] for i in range(nk3[0]*nk3[1]*nk3[2])] for ib in range(No_of_fermi_bands)]
    kp_no=0
    for ib in range(No_of_fermi_bands): 
     for i0 in range(nk3[0]):
        i0p = (i0 + 1)% (nk3[0])
        i0m = (i0 - 1)% (nk3[0])
        for i1 in range(nk3[1]):
          i1p= (i1 + 1)% (nk3[1])
          i1m = (i1 - 1)% (nk3[1])
          for i2 in range(nk3[2]):
            i2p = (i2 + 1)% (nk3[2])
            i2m = (i2 - 1)% (nk3[2])
            kp_no=i0*nk3[1]*nk3[2]+i1*nk3[2]+i2
            kpi0p=all_kpoints_index[i0p*nk3[1]*nk3[2]+i1*nk3[2]+i2]
            kpi0m=all_kpoints_index[i0m*nk3[1]*nk3[2]+i1*nk3[2]+i2]
            kpi1p=all_kpoints_index[i0*nk3[1]*nk3[2]+i1p*nk3[2]+i2]
            kpi1m=all_kpoints_index[i0*nk3[1]*nk3[2]+i1m*nk3[2]+i2]
            kpi2p=all_kpoints_index[i0*nk3[1]*nk3[2]+i1*nk3[2]+i2p]
            kpi2m=all_kpoints_index[i0*nk3[1]*nk3[2]+i1*nk3[2]+i2m]
            de[0] = ENE[ib][kpi0p] - ENE[ib][kpi0m]
            de[1] = ENE[ib][kpi1p] - ENE[ib][kpi1m]
            de[2] = ENE[ib][kpi2p] - ENE[ib][kpi2m]
            de=[de[i]*nk3[i]*0.5 for i in range(3)]
            vf_vector[ib][kp_no]=[\
              (avec[0][ii] * de[0] + avec[1][ii] * de[1] + avec[2][ii] * de[2]) for ii in range(3)]
            vf[ib][kp_no]=sum([\
              (avec[0][ii] * de[0] + avec[1][ii] * de[1] + avec[2][ii] * de[2])**2 for ii in range(3)])**0.5
    h=open('fermi_velocity.dat','w')
    h.write(str(nk3)+str(kp_no)+'\n')
    for i in vf_vector: 
     for j in i: h.write(str(j)+'\n')
    h.close()
    return vf

def write_to_fermisurfer(bvec,ENE,all_kpoints_index,Ef,vf=1):
#ENE[i][j] i-kpoint, j- band
 n_k2=int(round(len(all_kpoints_index)**(1/3.)))
 h=open('FS.frmsf','w')
 for i in range(3): h.write(str(n_k2)+' ')
 h.write('\n1\n')
 h.write(str(len(ENE))+'\n')
 for i in bvec:
  for j in i: h.write(str(j)+' ')
  h.write('\n')
 for i in ENE:
  for k in all_kpoints_index:
   h.write(str((i[k]-Ef))+'\n')
 if vf!=1:
  for i in range(len(ENE)):
   for k in range(len(all_kpoints_index)):
     h.write(str(vf[i][k])+'\n')
 h.close()



def write_to_orbital_fermisurfer(bvec,ENE,all_kpoints_index,Ef,ALMBLM,overlap,ene_weights,dosfile,qtl):
#ENE[i][j] i-kpoint, j- band
#ALM[k][at][jat][l][m][x]
 h=open(dosfile)
 dosy=[float(i) for i in h.readlines()[-4].split()[3:]]
 h.close()
 n_k2=int(round(len(all_kpoints_index)**(1/3.)))
 no_at=0
 for at in range(len(ALMBLM[0])): #over atoms
  for jat in range(len(ALMBLM[0][at])): #over atoms
   no_at+=1
   for l in range(len(ALMBLM[0][at][jat][0])): #over l
    suma=0.
 #  for m in range(len(ALMBLM[0][at][0][l])): #over atoms
    h=open('FS_atom_'+str(no_at)+'_orb_'+str(l)+'.frmsf','w')
    for i in range(3): h.write(str(n_k2)+' ')
    h.write('\n1\n')
    h.write(str(len(ENE))+'\n')
    for i in bvec:
     for j in i: h.write(str(j)+' ')
     h.write('\n')
    for i in ENE:
     for k in all_kpoints_index:
      h.write(str((i[k]-Ef))+'\n')
    for i in range(len(ENE)):
     wsp=[]
     wsp2=[]
     for k in range(len(ENE[i])):
      wsps=[ sum([sum([sum([sum([ ALMBLM[k][at][jat][i][l][m][x]*np.conjugate(ALMBLM[k][at][jat][i][l2][m2][x2])*overlap[at][l2][l][x2][x] for x in range(len(ALMBLM[k][at][jat][i][l][m]))]) for x2 in range(len(ALMBLM[k][at][jat][i][l2][m2]))])  for m2 in range(len(ALMBLM[0][at][jat][0][l2])) ]) for l2 in range(len(ALMBLM[0][at][jat][0])) ])  for m in range(len(ALMBLM[0][at][jat][0][l])) ] 
    #  wsps=[ sum([sum([ ALMBLM[k][at][jat][i][l][m][x]*np.conjugate(ALMBLM[k][at][jat][i][l][m][x2])*overlap[at][l][l][x2][x] for x in range(len(ALMBLM[k][at][jat][i][l][m]))]) for x2 in range(len(ALMBLM[k][at][jat][i][l][m]))])    for m in range(2*l+1) ]  #we chooose only l=l2,m=m2 because of the orthogonality of spherical harmonics 
      wsp.append(abs(sum(wsps)))
      wsp2.append(wsp[-1]*ene_weights[i][k])
     for k in all_kpoints_index: #radwf_norm[at][l][x]
#       wsps=[sum( [abs(ALMBLM[k][at][i][l][m][x])**2 for x in range(len(ALMBLM[k][at][i][l][m]))])  for m in range(len(ALMBLM[0][at][0][l])) ] 
      
#      wsp=sum( [abs(ALMBLM[k][at][i][l][m][0])**2  for m in range(len(ALMBLM[0][at][0][l])) ]  ) *ene_weights[i][k] 
       h.write(str(wsp[k])+'\n')
       suma+=wsp2[k]
    print(l,':dos=',suma/len(all_kpoints_index)),
    try: print(dosy[l]/(suma/len(all_kpoints_index)))
    except: continue
    h.close()

def write_to_qtl_fermisurfer(bvec,ENE,all_kpoints_index,Ef, ene_weights,dosfile,qtl):
#ENE[i][j] i-kpoint, j- band
#ALMBLM[k][at][j][l][m][0-3]
 #print(ene_weights)
 h=open(dosfile)
 dosy=[float(i) for i in h.readlines()[-4].split()[3:]]
 h.close()
 n_k2=int(round(len(all_kpoints_index)**(1/3.)))
 for at in range(len(qtl)): #over atoms
#  for l in range(len(ALMBLM[0][at][0])): #over l
  for l in range(len(qtl[at][0][0])): #over l
    suma=0
 #  for m in range(len(ALMBLM[0][at][0][l])): #over atoms
    h=open('FS_atom_'+str(at)+'_orb_'+str(l)+'.frmsf','w')
    for i in range(3): h.write(str(n_k2)+' ')
    h.write('\n1\n')
    h.write(str(len(ENE))+'\n')
    for i in bvec:
     for j in i: h.write(str(j)+' ')
     h.write('\n')
    for i in ENE:
     for k in all_kpoints_index:
      h.write(str((i[k]-Ef))+'\n')
    for i in range(len(ENE)):
     wsp=[]
     for k in range(len(ENE[i])):
      wsp.append(qtl[at][i][k][l] )
     for k in all_kpoints_index: 
       h.write(str(wsp[k])+'\n')
       suma+=wsp[k]*ene_weights[i][k]
    print(l,suma/len(all_kpoints_index)),
    try: print(dosy[l]/(suma/len(all_kpoints_index)))
    except: continue
    h.close()

def rot_z(kv):
 t=np.pi/4
 Rz=np.array([ [np.cos(t),-np.sin(t),0],[np.sin(t),np.cos(t),0] ,[0,0,1]])
 return np.matmul(Rz,kv)

def kpoints_from_cartesian_to_crystal(n_k_total,b_vec):
 nk=round(n_k_total**(1/3))

# print (nk)
 KPOINTS=[]
 m=0
 for i in range(nk): 
  for j in range(nk):
   for k in range(nk):
    kv=np.transpose(np.array([i,j,k]))
#    kv=rot_z(kv)
    KPOINTS.append([m for m in np.matmul(np.transpose(b_vec),kv)]+[m])
    m+=1
# for i in KPOINTS: 

# for i in KPOINTS:
#  for j in range(3):
#   if i[j]<0: i[j]+=30.288209999999996
 return sorting(KPOINTS)

def read_all_kpoints(kgenfile,b_vec):
        #read all kpoints 
        f=open(kgenfile,'r')
        tmp=f.readlines()
        f.close()
        ALL=[]
 #       print(kgenfile)
        for ni,i in enumerate(tmp):
         if 'DIVISION OF RECIPROCAL LATTICE VECTORS' in i:
          nk3=[ int(m) for m in i.split()[6:]]
#          print (nk3)
         if 'point' in i:
          for j in (tmp[ni+1:ni+1+int((nk3[0]+1)*(nk3[1]+1)*(nk3[2]+1))]):
                k=j.split()
                ALL.append( [int(k[0]),int(k[4])]) #number of point and number of point, to which he is equivalent
#                if k[0]==k[4]: print(np.matmul(np.transpose(b_vec),np.array([int(k[1]),int(k[2]),int(k[3])])))
          break
 #       print (len(ALL),nk3)
        #matrix: equiv[i] inform, that i-th kpoint is equiv to equiv[i]-th kpoint
        equiv=[]
        noneq=[]
        noneq_index=[]
        for ni,i in enumerate(ALL):
          if i[0]==i[1]: #if is nonequivalent
           equiv.append(len(noneq))
           noneq.append(i[0]-1)
          else: equiv.append(noneq.index(i[1]-1)) #is is equivalent to i[1]-th point
        print('number of non-equivalent k-points: '+str(len(noneq)))

        equiv2=[]
        for i in range(nk3[0]):
         for j in range(nk3[1]):
          for k in range(nk3[2]):
           equiv2.append(equiv[i*(nk3[0]+1)*(nk3[1]+1)+j*(nk3[2]+1)+k])
#        print(len(equiv2))
        nk_noneq=len(noneq)
        return equiv2,nk3,nk_noneq #2

def read_bvec(fil):
 b_vec=[]
 f=open(fil,'r')
 tmp=f.readlines()
 f.close()
 for ni,i in enumerate(tmp):
  if 'R1 = ' in i:
   a_vec=np.array([ [ float(m) for m in j.split()[2:]] for j in tmp[ni:ni+3]])
  if 'G1        G2        G3' in i:
   b_vec=np.array([ [ float(m) for m in j.split()] for j in tmp[ni+1:ni+4]]).transpose()
   break
 return b_vec,a_vec


def read_ef_and_dos(dos_file,prefix):
 print('Read EF and total DOS from '+dos_file+'...')
 E_f,dos=0,0
 try:
  h=open(dos_file,'r')
  tmp=h.readlines()
  h.close()
  nat=0
  for nj,j in enumerate(tmp):
   if 'JATOM' in j: nat+=1
   if 'EF and DOS at fermi level' in j:
     E_f,dos=[round(float(m),5) for m in tmp[nj+1].split()[0:2]]
     break
  if E_f==0: raise ValueError('Fermi Energy not found in files')
 except:
  dos,nat=0,0
  try: h=open(prefix+'.scf2')
  except: h=open(prefix+'.scf2up')
  for i in h.readlines():
#   if ':NAT0' in i: nat=int(i.split()[2])
   if ':FER' in i: E_f=float(i.split()[-1])
  h.close()
 if E_f==0: raise ValueError('Fermi Energy not found in files')
 print (E_f,dos,nat)
 return E_f,dos,nat


#ENE[i][j] i-kpoint, j- band
#ALM[i][l][m][k][j]
#ENE_weights[k][iband]
def which_bands_cross_ef(ENE,ALMBLM,qtl,E_f):
 chosen_bands=[]
 nbnd=min([len(k) for k in ENE])
 band_ranges=[[999,-999] for i in range(nbnd)]
 for k in range(1,len(ENE)):
  for i in range(nbnd):
   if ENE[k][i]<band_ranges[i][0]: band_ranges[i][0]=ENE[k][i]
   elif ENE[k][i]>band_ranges[i][1]: band_ranges[i][1]=ENE[k][i]
 for ni,i in enumerate(band_ranges):
  print ni,i
  if i[0]<E_f and i[1]>E_f: chosen_bands.append(ni)

 print ('bands crossing ef: ',chosen_bands,'EF=',E_f)
 mmin=min(chosen_bands)
 mmax=max(chosen_bands)+1
 #print(np.array(ENE[:][mmin:mmax]).shape)
 ENE= np.transpose(np.array([ k[mmin:mmax] for k in ENE]))
# for i in only_ENE_weights:
#  i=i/sum(i)
 No_of_fermi_bands=len(ENE)
 print(np.array(ALMBLM).shape)
 ALMBLM0=list(ALMBLM)
 ALMBLM=[[[ jat[mmin:mmax] for jat in at] for at in k] for k in ALMBLM] 
 qtl=[ at[mmin:mmax] for at in qtl] 
 for k in range(len(ALMBLM)):
  for at in range(len(ALMBLM[k])): 
   for jat in range(len(ALMBLM[k][at])): 
    for m in range(No_of_fermi_bands-len(ALMBLM[k][at][jat])): ALMBLM[k][at][jat].append([[[0 for i in j] for j in z] for z in ALMBLM0[k][at][jat][0]])

 return ENE,ALMBLM,qtl,No_of_fermi_bands

def read_qtl(qtl_file,so):
 print('Reading alm coefficients and kpoints from '+qtl_file+'...')
 tmp=[]
 for i in range(1,64):
  try: 
   h=open(qtl_file+'_'+str(i),'r')
   tmp.extend([m.split()  for m in h.readlines() if len(m.split())!=0])
   h.close()
  except:
   break
 if len(tmp)==0:
   h=open(qtl_file,'r')
   tmp.extend([m.split()  for m in h.readlines() if len(m.split())!=0])
   h.close()

 qtl=[]
 for i in tmp:
  if 'JATOM' in i:
   qtl.append([])
   orbitals=i[-1].split(',')
   chosen_orb=[]
   for ni,i in enumerate(orbitals):
    if i in [str(m) for m in range(9)]: chosen_orb.append(ni)
   print(orbitals,chosen_orb)
  elif 'BAND' in i: qtl[-1].append([])
  elif len(i)>4 and i[1]=='1': 
   qtl[-1][-1].append([])
   for no in chosen_orb: qtl[-1][-1][-1].append(float(i[no+2]))
 print('qtl here',len(qtl),len(qtl[0]),len(qtl[0][0]),len(qtl[0][0][0]))

 if so:
  for at in qtl:
   for nb in range(0,len(at),2):
    for nk in range(len(at[nb])):
     for nl in range(len(at[nb][nk])):
       at[nb][nk][nl]=(at[nb][nk][nl]+at[nb+1][nk][nl])/2.
       at[nb+1][nk][nl]=at[nb][nk][nl]
    
 return qtl #qtl[at][band][k][l]




def read_almblm(almblm_file,n_l,n_k2):
 #write(24,4893)l,m,index,alm(INDEX,NUM),blm(INDEX,NUM),(clm(INDEX,NUM,jlo),jlo=1,3)
 print('Reading alm coefficients and kpoints from '+almblm_file+'...')
 tmp=[]
 for i in range(1,64):
  try: 
   h=open(almblm_file+'_'+str(i),'r')
   tmp.extend([m.split()  for m in h.readlines() if len(m.split())!=0])
   h.close()
  except:
   break
 if len(tmp)==0:
   h=open(almblm_file,'r')
   tmp.extend([m.split()  for m in h.readlines() if len(m.split())!=0])
   h.close()
 #ALM are in file only for occupied (or partially occupied) bands, so the number of bands is smaller here than in *energy files
 #in fact, ALM are equal to ! X_l,m,a = Sum(K) c_K/sqrt(V) exp(i(K+k)R_a) Y(*)_lm(T_a^-1(K+k)) X_l,a(|K+k|), (c_k = gaunt coeff.), to liczy augpw.frc, ale w almblm mamy to juz gotowe
 ALMBLM=[] 
 for ni,i in enumerate(tmp):
  if 'K-POINT' in i[0]: 
   kp=int(i[-1])-1
   if tmp[ni+1][0]=='1': #if jatom=1
    ALMBLM.append([[]])
   else:
    ALMBLM[kp].append([])
  elif 'ATOM' in i[1]: 
   ALMBLM[kp][-1].append([])
  elif 'weight' in i[-1]:  
   ALMBLM[kp][-1][-1].append( [ [] for k in range(n_l) ] )
#   weight=float(i[1])
  elif len(i)>0 and 'jatom,nemin,nemax' not in i:
   if int(i[0])<n_l:   #ALM[k][jat][at][l][m][x]
    ALMBLM[kp][-1][-1][-1][int(i[0])].append([complex(float(i[3]),float(i[4])),complex(float(i[5]),float(i[6])),complex(float(i[7]),float(i[8])),complex(float(i[9]),float(i[10])),complex(float(i[11]),float(i[12])) ]) 
#    ALMBLM[-1][-1][-1][int(i[0])].append([weight*complex(float(i[3]),float(i[4])),weight*complex(float(i[5]),float(i[6])),weight*complex(float(i[7]),float(i[8])),weight* complex(float(i[9]),float(i[10])),weight* complex(float(i[11]),float(i[12])) ]) 

 n_k=len(ALMBLM)
 n_at=len(ALMBLM[0])

 print('nk,nat=',n_k,n_at)
 No_of_bands=[ len(i[0]) for i in ALMBLM]
 if n_k2!=n_k: raise ValueError('no of k_points of energies:',n_k2,' and of almblm:',str(n_k),' is  not the same')

 '''
 for i in ALMBLM: 
  if len(i)<n_at: print '1'
 #transpose to get: from ALMBLM[kp][at] to ALMBLM[at][kp]
 ALMBLM0=[[ALMBLM[i][j] for i in range(len(ALMBLM)) if len(ALMBLM[i])-1<j] for j in range(n_at)]
 #ALMBLM[i][k][j][l][m][0-3]  i-atoms,k-kpoint, j-band, l (orbital No), m, [0: Re[Alm]+j*Im[Alm], 1:  Re[Blm]+j*Im[Blm]]
 

#alm and blm as numpy arrays
#ALMBLM[k][i][j][l][m][0-3] -> ALM[i][l][m][k][j]
 ALMBLM=[[[[np.array([ALMBLM0[i][k][j][l][m] for j in range(len(ALMBLM0[i][k]))]) for k in range(n_k) ] for m in range(2*l+1) ] for l in range(n_l) ] for i in range(n_at)]
 '''

 return ALMBLM,No_of_bands,n_k,n_at

def read_ene(ene_file,n_at,degauss,E_f):
 print('Reading band energies...'),
 ENE=[]
 tmp=[]
 print('Files:'),
 kweights=[]
# print nat
 for i in range(1,64):
  try: 
   h=open(ene_file+'_'+str(i),'r')
   print(ene_file+'_'+str(i)),
   tmp.extend([m.split() for m in h.readlines()[n_at*2:]])
   h.close()
  except:
   break
 if i==1:
   h=open(ene_file,'r')
   print(ene_file),
   tmp.extend([m.split() for m in h.readlines()[n_at*2:]])
   h.close()
 
#ENE[i][j] i-kpoint, j- band
 for ni,i in enumerate(tmp): 
   if len(i)>3  and len(tmp[ni+1])==2: 
     kweights.append(float(i[-1]))
     ENE.append([])
   elif len(i)==2: 
    try: ENE[-1].append(float(i[1]))
    except: continue
 n_k2,n_band=len(kweights),min([ len(i) for i in ENE])
 n_k_total=sum(kweights)
 print("Irreducible number of kpoints="+str(n_k2)+"; Total number of kpoints="+str(n_k_total))

 return ENE, kweights,n_k_total,n_k2

def calc_kweights(kweights,ENE,E_f,n_k,degauss):
 nb=len(ENE)
 ENE_weights=[ [] for j in range(nb)]
 only_ENE_weights=[ [] for j in range(nb)]
# weights for integrals over k (cold smearing)
#  ! cold smearing  (Marzari-Vanderbilt-DeVita-Payne)
#  if (n.eq. - 1) then
#     arg = min (200.d0, (x - 1.0d0 / sqrt (2.0d0) ) **2)
#     w0gauss = sqrtpm1 * exp ( - arg) * (2.0d0 - sqrt ( 2.0d0) * x)
#     return
#           w0g1 = w0gauss ( (ef1 - et (ibnd, ikk) ) / degauss1, ngauss1) &
#                / degauss1
 for ib in range(nb): 
  for k in range(n_k):
          eef=(ENE[ib][k]-E_f)/degauss
          arg=min(200., (eef - 2.**(-0.5) ) **2)
          only_ENE_weights[ib].append(max((sqrtpm1*np.exp(-arg)*(2.-(2**0.5*eef)))/degauss,0))
#          only_ENE_weights[k].append( 1/((2*np.pi)**0.5*degauss) *np.exp(-(eef)**2/2))
          ENE_weights[ib].append(kweights[k]*only_ENE_weights[ib][-1])
 return ENE_weights,only_ENE_weights

def adjust_degauss(kweights,ENE,E_f,n_k,n_k_total,dos):
 min_diff_dos=[1000,1000]
 for i in range(1,200):
  degauss=0.001*i
  weight=calc_kweights(kweights,ENE,E_f,n_k,degauss)[0]
  dos2=sum([ sum(m) for m in weight])/n_k_total
  if abs(dos2-dos)<min_diff_dos[1]: min_diff_dos=[degauss,abs(dos2-dos)]
 return min_diff_dos[0]
