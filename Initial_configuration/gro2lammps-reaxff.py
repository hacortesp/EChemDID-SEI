import numpy as np

def readgro(filename):
    atype = []; index = []; box = []; xyz=[]; diff_atype=[];mass=[]
    with open(filename) as f:
         l=0
         flg = 0
         for line in f:
             data = line.split()
             if len(data)==1:
                 natoms=int(data[0])
             if len(data) == 6 and l<natoms+2:
                 l+=1
                 index.append(int(data[2]))
                 atype.append(data[1])
                 xyz.append([float(s) for s in data[3:6]])
             if len(data) == 9:
                 for i in range(3):
                     box.append(float(data[i]))
         if 'P' in atype:
             diff_atype.append('P')
             mass.append(float(30.974))
         if 'F' in atype:
             diff_atype.append('F')
             mass.append(float(18.998))
         if 'Ox' in atype:
             diff_atype.append('Ox')
             mass.append(float(16.000))
         if 'Oc' in atype:
             diff_atype.append('Oc')
             mass.append(float(16.000))
         if 'Ca' in atype:
             diff_atype.append('Ca')
             mass.append(float(12.00))
         if 'Cc' in atype:
             diff_atype.append('Cc')
             mass.append(float(12.00))
         if 'Ce' in atype:
             diff_atype.append('Ce')
             mass.append(float(12.00))
         if 'LIa' in atype:
             diff_atype.append('LIa')
             mass.append(float(7.00))
         if 'LIc' in atype:
             diff_atype.append('LIc')
             mass.append(float(7.00))
         if 'LIi' in atype:
             diff_atype.append('LIi')
             mass.append(float(7.00))
         if 'LIe' in atype:
             diff_atype.append('LIe')
             mass.append(float(7.00))
         if 'H' in atype:
             diff_atype.append('H')
             mass.append(float(1))

    return box, index, atype, xyz, natoms, diff_atype, mass

def write_file(outfile,tot_natoms,natoms,a,b,c,mass,atype,ntypes,coord):
    o = open(outfile,'w')
    o.write('New lammps configuration\n')
    o.write(' {0} {1}\n'.format(tot_natoms,'atoms'))
    o.write(' {0} {1}\n'.format(0,'bonds'))
    o.write(' {0} {1}\n'.format(0,'angles'))
    o.write(' {0} {1}\n'.format(0,'dihedrals'))
    o.write(' {0} {1}\n'.format(0,'impropers'))
    o.write(' {0} {1}\n'.format(len(mass),'atom types'))
    o.write(' {0} {1}\n'.format(0,'bond types'))
    o.write(' {0} {1}\n'.format(0,'angle types'))
    o.write(' {0} {1}\n'.format(0,'dihedral types'))
    o.write(' {0} {1}\n'.format(0,'improper types'))
    o.write(' {0} {1} {2}\n'.format(0.00000000,a,'xlo xhi'))
    o.write(' {0} {1} {2}\n'.format(0.00000000,b,'ylo yhi'))
    o.write(' {0} {1} {2}\n\n'.format(0.00000000,c,'zlo zhi'))
    o.write(' Masses\n\n')
    for i in range(len(mass)):
        o.write(' {0} {1}\n'.format(i+1,mass[i]))
    o.write('\n{0}\n\n'.format(' Atoms # Full'))
    count = 0
    for i in range(ntypes):
        for j in range(natoms[i]):
            count += 1
            o.write('{:<} {:<} {:<} {:<8.4f}{:^12.6f}{:^12.6f}{:^12.6f}\n'.format(count,1,i+1,0.0,(coord[count-1][0])*a,(coord[count-1][1])*b,(coord[count-1][2])*c))

#the type names in conf.gro must have the element symbol
in_name='conf.gro'
outfile='data.lammps'

box, index, atype, xyz, natoms, diff_atype, mass =readgro(in_name)
#print(natoms)
#P F O OC C LI H elements, OC=O in carboxyl group
P=[];F=[];Ox=[];Oc=[];Ca=[];Ce=[];Cc=[];LIa=[];LIc=[];LIe=[];LIi=[];H=[]
for i in range(natoms):
    xyz[i][0]=xyz[i][0]/box[0]
    xyz[i][1]=xyz[i][1]/box[1]
    xyz[i][2]=xyz[i][2]/box[2]
    if 'P' in atype[i]:
        P.append(xyz[i])
    if 'F' in atype[i]:
        F.append(xyz[i])
    if 'Ox' in atype[i]:
        Ox.append(xyz[i])
    if 'Oc' in atype[i]:
        Oc.append(xyz[i])
    if 'Ca' in atype[i]:
        Ca.append(xyz[i])
    if 'Cc' in atype[i]:
        Cc.append(xyz[i])
    if 'Ce' in atype[i]:
        Ce.append(xyz[i])
    if 'LIa' in atype[i]:
        LIa.append(xyz[i])
    if 'LIc' in atype[i]:
        LIc.append(xyz[i])
    if 'LIi' in atype[i]:
        LIi.append(xyz[i])
    if 'LIe' in atype[i]:
        LIe.append(xyz[i])
    if 'H' in atype[i]:
        H.append(xyz[i])

atoms=P+F+Ox+Oc+Ca+Cc+Ce+LIa+LIc+LIi+LIe+H

a=10*box[0]; b=10*box[1]; c=10*box[2]
na_type=[len(P),len(F),len(Ox),len(Oc),len(Ca),len(Cc),len(Ce),len(LIa),len(LIc),len(LIi),len(LIe),len(H)]
natype_aux=[]
for i in na_type:
    if i!=0:
        natype_aux.append(i)

tot_natoms=len(atype)
tot_types=len(diff_atype)

print('number of atoms:', tot_natoms)
print('number of atoms of each type:')
print('P',na_type[0],'F',na_type[1],'Ox',na_type[2],'Oc',na_type[3],'Ca',na_type[4],'Cc',na_type[5],'Ce',na_type[6],'LIa',na_type[7],'LIc',na_type[8],'LIi',na_type[9],'LIe',na_type[10],'H',na_type[11])
#print(natype_aux)
#print(a,b,c)
#print(mass)
#print(diff_atype)
#print(tot_types)
#print(len(atoms))
#print(atype)

write_file(outfile,tot_natoms,natype_aux,a,b,c,mass,diff_atype,tot_types,atoms)
