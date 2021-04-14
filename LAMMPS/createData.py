#Lammps data file for simulation of polyelectrolyte coplexation
#creates initial positions of monomers
#Created on Tue Oct 24 2019
#@author: sabin
import numpy as np
import random
#import math
import copy

l = 1 #bead separation
chain_length = 60#no of beads in a chain
bond_length = chain_length-1 #no of bonds in a chain
n_chain = 160           #no. of chains



b1 = 31.0
b2 = 31.0
b3 = 31.0

xMin = -b1
yMin = -b2
zMin = -b3
xMax = b1
yMax = b2
zMax = b3
n_atoms = n_chain*chain_length
n_bonds = n_chain*(bond_length)
bondLength=1     #bond length
threshold=0.5    #threshold for overlap


alpha = 0.6
charge_mag = 14.97
Ones = int(alpha*chain_length)#no of charged beads --make sure it is an integer
Zeros = int(chain_length - Ones)
lis_a = list(charge_mag*np.ones(Ones))
lis_b = list(np.zeros(Zeros))
lis_c = list(-charge_mag*np.ones(Ones))
lis_ab = lis_a+lis_b
lis_cb = lis_c+lis_b
random.shuffle(lis_ab)
random.shuffle(lis_cb)
n_count = Ones*(n_chain/2) #no of counterions of one type
ns = 1121   #number of salt ions of one charge type 
n_salt=2*ns   #
#print(lis_ab)
n_atoms = n_atoms+2*n_count+n_salt
###########
outfile = open('data.complexation','w')
outfile.write('LAMMPS Description\n\n{0} atoms\n{1} bonds\n\n'.format(n_atoms,n_bonds))
outfile.write('4 atom types\n1 bond types\n\n')
outfile.write('{0} {1} xlo xhi\n'.format(-b1, b1))
outfile.write('{0} {1} ylo yhi\n'.format(-b2, b2))
outfile.write('{0} {1} zlo zhi\n\n'.format(-b3, b3))
outfile.write('Masses\n\n1 1.0\n2 1.0\n3 1.0\n4 1.0\n\n')
outfile.write('Atoms\n\n')

def distance(x1,y1,z1,x2,y2,z2):
	return ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5

 
def minDistance1(x1,y1,z1,lst):
	return min([distance(x1,y1,z1, i[0], i[1], i[2]) for i in lst])


xy_list =[[np.random.uniform(-b1,b1), np.random.uniform(-b2,b2), np.random.uniform(-b3,b3)]]

for i in range(1,n_chain):
    x_temp = np.random.uniform(-b1,b1)
    y_temp = np.random.uniform(-b2,b2)
    z_temp = np.random.uniform(-b3,b3)
    #rint(x_temp)
    while minDistance1(x_temp,y_temp, z_temp, xy_list)<threshold:
        x_temp = np.random.uniform(-b1,b1)
        y_temp = np.random.uniform(-b2,b2)
        z_temp = np.random.uniform(-b3,b3)
    x = x_temp
    y = y_temp
    z = z_temp
    xy_list.append([x,y,z])
    
#print(xy_list)
	

atom_no =0

globalLst = copy.copy(xy_list)


for i in xy_list:
    mol_id = xy_list.index(i)+1
    if xy_list.index(i)%2==0:
        charge_lis = lis_ab
        atom_type = 1
    else:
        charge_lis = lis_cb
        atom_type = 2
        
    for j in range(chain_length):
        atom_no = atom_no +1
        charge = charge_lis[j]
        if j == 0:
            x = i[0]
            y = i[1]
            z = i[2]
        
        if j>0:
            theta = np.pi*np.random.random()
            phi = 2*np.pi*np.random.random()
            xTemp = x+bondLength*np.sin(theta)*np.cos(phi)
            yTemp = y+bondLength*np.sin(theta)*np.sin(phi)
            zTemp = z+bondLength*np.cos(theta)
            while abs(xTemp)>b1 or abs(yTemp)>b2 or abs(zTemp)>b3 or minDistance1(xTemp, yTemp, zTemp, globalLst)<threshold:
                theta = np.pi*np.random.random()
                phi = 2*np.pi*np.random.random()
                xTemp = x+bondLength*np.sin(theta)*np.cos(phi)
                yTemp = y+bondLength*np.sin(theta)*np.sin(phi)
                zTemp = z+bondLength*np.cos(theta)
            x = xTemp
            y = yTemp
            z = zTemp
            globalLst.append([xTemp, yTemp, zTemp])
            print(atom_no)
        
        outfile.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(int(atom_no),int(mol_id),int(atom_type),charge,x,y,z))


####Create counterions here
#positive

nIons = n_count+ns

mol_id = mol_id+1
print(mol_id)
for i in range(nIons):
	atom_no = atom_no+1
	atom_type = 3
	charge = charge_mag
	x = b1*np.random.random()
	y = b2*np.random.random()
	z = b3*np.random.random()
      
        
	outfile.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(int(atom_no),int(mol_id),int(atom_type),charge,x,y,z))

#negative
mol_id = mol_id+1
print(mol_id)
for i in range(nIons):
	atom_no = atom_no+1
	atom_type = 4
	charge = -charge_mag
	x = b1*np.random.random()
	y = b2*np.random.random()
	z = b3*np.random.random()
        
	outfile.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(int(atom_no),int(mol_id),int(atom_type),charge,x,y,z))

##Create Bonds
##Bonds

outfile.write('\nBonds\n\n')
i1=1
while i1<=n_chain:
	bond_count = 1
	while bond_count <= bond_length:
		bond_no = (i1-1)*bond_length + bond_count 
		atom1 = (i1-1)*chain_length+bond_count
		atom2 = atom1+1
		bond_type = 1
		outfile.write('{0} {1} {2} {3}\n'.format(bond_no,bond_type,atom1,atom2))
		bond_count = bond_count+1

	i1=i1+1
		

		
outfile.close()

