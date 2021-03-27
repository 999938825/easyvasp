import os,sys,shutil,re
from easyvasp.incar import InCar
from pymatgen.io.vasp.inputs import *
from pymatgen.io.vasp.outputs import *
from pymatgen.io.cif import CifParser
from pymatgen.io.xyz import XYZ
from pymatgen.core.periodic_table import Element
from pymatgen.core import Structure, Lattice, Molecule,SymmOp
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.operations import SymmOp
from pymatgen.analysis.local_env import VoronoiNN, JmolNN
from scipy.spatial import Voronoi

class easyjob(object):
    def __init__(self, structure, foldname):
        self.foldname = foldname
        self.spath = os.getcwd()
        self.rootpath = os.path.split(os.path.realpath(__file__))[0]
        self.poscar = Poscar(structure)
        self.potcar_file = None
        self.potcar_list = []
        self.__dict__['potcar_style'] = 'pbe'
        self.kpoint_file = None
        self.kpoint_mode='G'
        self.kpoint_list=[1,1,1]
        self.incar = None
        self.nelect = 0
        self.eledata = []
        self.load()

    def __setattr__(self, key, value):
        if key=='foldname':
            self.__dict__[key] = value
            self.__dict__['path'] = os.path.join(os.getcwd(),self.foldname)
        elif key == 'kpoint':
            if value=='line':
                self.kpoint_mode ='line'
            elif isinstance(value,float) or isinstance(value,int):
                self.kpoint_mode=='G'
                self.automatic_set_kpoint(value)
            elif isinstance(value,list) and len(value)==3:
                self.kpoint_mode == 'G'
                self.kpoint_list=value
            else:
                raise Exception("The KPoint format is incorrect!")
        elif key == 'potcar_style':
            self.__dict__[key] = value
            self.refresh_potcar()
        elif key == 'cn':
            return self.get_ave_cn()
        else:
            self.__dict__[key] = value

    def __getitem__(self, index):
        return self.poscar.structure[index]

    def __len__(self):
        return len(self.poscar.structure)

    def load(self):
        self.reset_loc(only_jd=True)
        self.get_incar('optim')
        self.kpoint = 0.2
        self.refresh_potcar()

    @staticmethod
    def read_file(filename):
        '''
        根据path读取各个文件，不论是vasp还是cif还是xyz
        :param filename:  文件的路径
        :return:   如果是vasp或者cif 则返回structure对象，如果是xyz则返回mol对象
        '''
        suffix=os.path.split(filename)[-1]
        suffix=suffix.split('.')[-1]
        if suffix=='vasp' or suffix=='POSCAR' or suffix=='CONTCAR':
            structure = Poscar.from_file(filename).structure
        elif suffix=='cif':
            structure =CifParser(filename,occupancy_tolerance=0.7).get_structures(primitive=False)[0]
        elif suffix=='xyz':
            structure = XYZ.from_file(filename).molecule
        elif suffix == 'CHGCAR' or suffix == 'CHG':
            structure = Chgcar.from_file(filename)
        return structure

    def get_script(self):
        #复制脚本到当前目录
        shutil.copy(os.path.join(self.rootpath, 'file', 'sub.sh'), os.path.dirname(self.path))
        shutil.copy(os.path.join(self.rootpath, 'file', 'vaspjob.pbs'), os.path.dirname(self.path))

    def get_distance(self,index1,index2):
        return self.poscar.structure.get_distance(index1,index2)

    def get_angle(self,index1,index2,index3):
        return self.poscar.structure.get_angle(index1,index2,index3)

    def reset_loc(self,distance=2,only_jd=False):
        k = self.poscar.structure.frac_coords[:, 2]
        k[k < 0] = k[k < 0] + 1
        k.sort()
        k1 = np.append(k[1:], k[0])
        dk = k1 - k
        dk[-1] = dk[-1] + 1
        if self.poscar.structure.lattice.c*dk.max() >= 6:
            if not only_jd:
                dmax = k[dk.argmax()]
                if len(k) > dk.argmax() + 1:
                    dmin = k[dk.argmax() + 1]
                else:
                    dmin = k[0]
                zz = self.poscar.structure.lattice._matrix[2] * (-dmin + (distance / self.poscar.structure.lattice.c))
                self.move(zz)
            self.issuf = True
        else:
            self.issuf = False
        minnum = np.min(self.poscar.structure.cart_coords, axis=0)[-1]
        maxnum = np.max(self.poscar.structure.cart_coords, axis=0)[-1]
        self.vacuum_thick = self.poscar.structure.lattice.c - (maxnum - minnum)

    def c_to_f(self,loc):
        return self.poscar.structure.lattice.get_fractional_coords(loc)

    def f_to_c(self,loc):
        return self.poscar.structure.lattice.get_cartesian_coords(loc)

    @property
    def frac_coords(self):
        f = self.poscar.structure.frac_coords
        f = np.insert(f, 0, [0,0,0],axis=0)
        return f

    @property
    def cart_coords(self):
        f = self.poscar.structure.cart_coords
        f = np.insert(f, 0, [0,0,0],axis=0)
        return f

    def refresh_potcar(self,recommend=True):
        self.potcar_file=None
        self.potcar_list=[]
        for i in self.poscar.site_symbols:
            if recommend and self.potcar_style in ['paw','pbe']:
                if i in ['Li', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Cs', 'Ba', 'Fr', 'Ra']:
                    self.potcar_list.append(str(i) + '_sv')
                elif i in ['Na', 'Cr', 'Mn', 'Tc', 'Ru', 'Rh', 'Hf', 'Ta', 'W']:
                    self.potcar_list.append(str(i) + '_pv')
                elif i in ['Ga', 'Ge', 'In', 'Sn', 'Tl', 'Pb', 'Bi', 'Po', 'At']:
                    self.potcar_list.append(str(i) + '_d')
                elif i in ['Pr', 'Nd', 'Pm', 'Sm', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Lu']:
                    self.potcar_list.append(str(i) + '_3')
                elif i in ['Eu', 'Yb']:
                    self.potcar_list.append(str(i) + '_2')
                else:
                    self.potcar_list.append(i)
            else:
                pot =  self.potcar_list.append(i)
        self.get_potcar()

    def get_potcar(self):
        self.potcar_file = ''
        filepath = os.path.join(self.rootpath, 'Pseudopotentials', self.potcar_style)
        POTCAR=''
        for i in self.potcar_list:
            tmppath = os.path.join(filepath, str(i), 'POTCAR')
            with open(tmppath, 'r') as f:
                tmptext = f.read()
                POTCAR = POTCAR + tmptext
        self.potcar_file = POTCAR
        self.get_nelect()
        return POTCAR

    def automatic_set_kpoint(self,kppa):
        latt = self.poscar.structure.lattice
        a=40/(latt.abc[0]*6.2832*kppa)
        b=40/(latt.abc[1]*6.2832*kppa)
        c=40/(latt.abc[2]*6.2832*kppa)
        a = int(round(max(1, a)))
        b = int(round(max(1, b)))
        c = int(round(max(1, c)))
        if self.issuf:
            c=1
        self.kpoint_list[0]=a
        self.kpoint_list[1]=b
        self.kpoint_list[2]=c

    def get_kpoints(self):
        x = self.kpoint_list[0]
        y = self.kpoint_list[1]
        z = self.kpoint_list[2]
        model = self.kpoint_mode
        filepath = os.path.join(self.rootpath, 'file', 'kpoints_line.txt')
        if model == 'line':
            with open(filepath, 'r') as f:
                self.kpoint_file = f.read()
        elif model == 'G':
            text = 'Auto' + '\n' + '0' + '\n' + 'G' + '\n' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n'
            self.kpoint_file = text

    def get_incar(self, incarname):
        filepath = os.path.join(self.rootpath, 'file', 'INCAR_' + incarname + '.txt')
        if os.path.isfile(filepath):
            incar = InCar(filepath)
            self.incar = incar
            if self.issuf:
                self.incar.ISIF = 2
            else:
                self.incar.ISIF = 3
            isip2list = [Element.Fe, Element.Co, Element.Ni, Element.Mn, Element.Cr]
            for i in isip2list:
                if i in set(self.poscar.structure.species):
                    self.incar.set('ISPIN', 2)
            return incar
        else:
            errorinfo = 'INCAR-%s文件不存在' % (incarname)
            raise RuntimeError(errorinfo)

    def fix(self,splitpoint,reverse=False,abs_coords=True,axis='z',v='xyz'):
        axisdict={'x':0,'y':1,'z':2}
        axis=axisdict[axis]
        if v=='xyz':
            vv=[False, False, False]
        elif v=='xy':
            vv = [False, False, True]
        elif v=='x':
            vv=[False, True, True]
        elif v=='y':
            vv = [True, False, True]
        elif v=='z':
            vv = [True, True, False]
        else:
            vv= [True, True, True]

        if not abs_coords:
            minz=self.poscar.structure.frac_coords.min(axis=0)[axis]
            maxz=self.poscar.structure.frac_coords.max(axis=0)[axis]
            splitz=(maxz-minz)*splitpoint+minz
        else:
            splitz=splitpoint
        fixmat=[]
        fixnum=0
        for i in self.poscar.structure.frac_coords:
            if reverse:
                if i[axis] >= splitz:
                    fixnum += 1
                    fixmat.append(vv)
                elif i[axis] < splitz:
                    fixmat.append([True, True, True])
            else:
                if i[axis] <= splitz:
                    fixnum += 1
                    fixmat.append(vv)
                elif i[axis] > splitz:
                    fixmat.append([True, True, True])
        print('fix '+str(fixnum)+' atoms')
        self.poscar.selective_dynamics=fixmat

    def set_vacuum(self,thick=20):
        self.reset_loc()
        new_lattice=self.poscar.structure.lattice
        new_c=new_lattice.c-self.vacuum_thick+thick
        new_lattice=Lattice.from_parameters(new_lattice.a,new_lattice.b,new_c,new_lattice.alpha,new_lattice.beta,new_lattice.gamma)
        min=np.mean(self.poscar.structure.frac_coords[:,2])
        newcoor=self.poscar.structure.frac_coords-[0,0,min]
        newcoor[:,2]=newcoor[:,2]*self.poscar.structure.lattice.c/new_c
        new_structure=Structure(new_lattice,self.poscar.structure.species,newcoor,coords_are_cartesian=False)
        self.poscar=Poscar(new_structure)
        self.reset_loc()

    def get_cn(self,i):
        vnn = VoronoiNN()
        cn = vnn.get_cn(self.poscar.structure, i, use_weights=True)
        return cn

    def get_ave_cn(self):
        coordination_numbers = {}
        vnn = VoronoiNN()
        for spec in self.poscar.structure.composition.as_dict().keys():
            coordination_numbers[spec] = 0.0
        for j, atom in enumerate(self.poscar.structure):
            cn = vnn.get_cn(self.poscar.structure, j, use_weights=True)
            coordination_numbers[atom.species_string] += cn
        elements = self.poscar.structure.composition.as_dict()
        for el in coordination_numbers:
            coordination_numbers[el] = coordination_numbers[el] / elements[el]
        return coordination_numbers

    def fix_element(self,eledata,v='xyz'):
        kk=[]
        if v=='xyz':
            vv=[False, False, False]
        elif v=='xy':
            vv = [False, False, True]
        elif v=='x':
            vv=[False, True, True]
        elif v=='y':
            vv = [True, False, True]
        elif v=='z':
            vv = [True, True, False]
        for i in self.poscar.structure:
            bds=str(i.species)[:-1]
            if bds in eledata:
                kk.append(vv)
            else:
                kk.append([True,True,True])
        self.poscar.selective_dynamics = kk

    def get_nelect(self):
        nelect = 0
        k = 0
        self.eledata = []
        potlist = self.potcar_file.split('\n')
        for index, item in enumerate(potlist):
            if len(item) > 2:
                if item.split()[0] in ['PAW_PBE','PAW_GGA','US']:
                    nelect += float(potlist[index + 1]) * self.poscar.natoms[k]
                    self.eledata.append(float(potlist[index + 1]))
                    k += 1
        self.nelect = nelect
        return int(nelect)

    def reverse(self,frac=True):
        oloc=np.mean(self.poscar.structure.frac_coords[:,2])
        opp = SymmOp.from_rotation_and_translation(rotation_matrix=((1, 0, 0), (0, -1, 0), (0, 0, -1)))
        self.poscar.structure.apply_operation(opp,fractional=frac)
        nloc=np.mean(self.poscar.structure.frac_coords[:,2])
        self.move(self.poscar.structure.lattice.get_cartesian_coords([0,0,oloc-nloc]))

    def move(self,v,frac=False):
        if frac:
            v0=self.f_to_c([0,0,0])
            v1=self.f_to_c(v)
            v=v1-v0
        opp = SymmOp.from_rotation_and_translation(rotation_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),translation_vec=v)
        self.poscar.structure.apply_operation(opp)

    def move_to(self, i, v, frac=True):
        if type(i) is int:
            v=v-self.poscar.structure.frac_coords[i]
            self.move_by(i,v,frac=frac)
        else:
            print('error,Unknown type')

    def move_by(self, i, v, frac=True):
        self.poscar.structure.translate_sites(i, v, frac_coords=frac)

    def find_path(self,end_structure,images,sites=[],charge=None):
        from pymatgen.analysis.path_finder import NEBPathfinder
        from pymatgen.analysis.path_finder import FreeVolumePotential
        from pymatgen.analysis.path_finder import ChgcarPotential
        init_struct = self.poscar.structure
        final_struct = end_structure
        if charge is None:
            dim=[50,50,120]
            v=FreeVolumePotential(final_struct, dim)
            obj=NEBPathfinder(init_struct,final_struct,sites,v.get_v(), n_images=images)
            obj.interpolate()
            new_path=obj.images
        else:
            v=ChgcarPotential(charge)
            obj = NEBPathfinder(init_struct, final_struct, [80,81,82], v.get_v(), n_images=images)
            obj.interpolate()
            new_path = obj.images

        joblist=[]
        for index,item in enumerate(new_path):
            filepath=os.path.join(self.foldname,'{0:02d}'.format(index))
            joblist.append(easyjob(item,filepath))
        return joblist

    def make_surface(self, miller_index=[0, 0, 1],bonds=None, min_slab_size=8, min_vacuum_size=20, max_normal_search=None,primitive=False):
        slabs = SlabGenerator(self.poscar.structure.copy(), miller_index=miller_index, min_slab_size=min_slab_size,
                              min_vacuum_size=min_vacuum_size,
                              max_normal_search=max_normal_search, primitive=primitive, lll_reduce=True)
        if bonds is not None:
            slabs = slabs.get_slabs(tol=0.1,bonds=bonds)
        else:
            slabs = slabs.get_slabs(tol=0.1)
        suflist = []
        i = 0
        for slab in slabs:
            i += 1
            newname = str(miller_index[0]) + str(miller_index[1]) + str(miller_index[2])+'-' + str(i)
            suflist.append(easyjob(slab, os.path.join(self.path, newname)))
        return suflist

    def try_orthogonal(self, maxnum=2):
        #试图通过扩胞使得xy方向正交
        sita = self.poscar.structure.lattice.angles[2]
        xa = self.poscar.structure.lattice.a * math.cos(sita / 180 * math.pi)
        ya = self.poscar.structure.lattice.a * math.sin(sita / 180 * math.pi)
        xb = self.poscar.structure.lattice.b
        minx = 10000
        for m1 in range(maxnum * -1, maxnum + 1):
            for m2 in range(maxnum * -1, maxnum + 1):
                for m3 in range(maxnum * -1, maxnum + 1):
                    for m4 in range(maxnum * -1, maxnum + 1):
                        p1 = math.sqrt(m1 * m1 + m2 * m2) * math.sqrt(m3 * m3 + m4 * m4)
                        if p1 > 0 and m1 * m4 - m2 * m3 != 0:
                            p2 = m1 * m3 * (xa * xa + ya * ya) + m2 * m4 * xb * xb + (m1 * m4 + m2 * m3) * xa * xb
                            sss = p2 / p1
                            if abs(sss) < minx:
                                minx = sss
                                nm1 = m1
                                nm2 = m2
                                nm3 = m3
                                nm4 = m4
        numzz = 0
        if nm1 > 0: numzz += 1
        if nm2 > 0: numzz += 1
        if nm3 > 0: numzz += 1
        if nm4 > 0: numzz += 1
        if numzz <= 1:
            nm1 = -1 * nm1
            nm2 = -1 * nm2
            nm3 = -1 * nm3
            nm4 = -1 * nm4
        print('set supercell: (', nm1, nm2, 0, ') (', nm3, nm4, 0, ')')
        self.poscar.structure.make_supercell([[nm1, nm2, 0], [nm3, nm4, 0], [0, 0, 1]])

    def get_primitive(self,tolerance=0.25):
        self.poscar.structure=self.poscar.structure.get_primitive_structure(tolerance=tolerance)

    def supercell(self, ii=1, jj=1, kk=1, amin=None, bmin=None, cmin=None):
        if amin == None:
            self.poscar.structure.make_supercell([ii, 1, 1])
        else:
            ii = int(amin / self.poscar.structure.lattice.a)+1
            self.poscar.structure.make_supercell([ii, 1, 1])
        if bmin == None:
            self.poscar.structure.make_supercell([1, jj, 1])
        else:
            jj = int(bmin / self.poscar.structure.lattice.b)+1
            self.poscar.structure.make_supercell([1, jj, 1])
        if cmin == None:
            self.poscar.structure.make_supercell([1, 1, kk])
        else:
            kk = int(cmin / self.poscar.structure.lattice.c)+1
            self.poscar.structure.make_supercell([1, 1, kk])
        print('set supercell: (', ii, jj, kk, ')')

    def adsorb(self,mol,loc=[0.5,0.5,0.5],spec_z=False, theta=[0, 0, 0], v=[0, 0, 0], distance=2, ifpca=True,frac=True):
        from sklearn.decomposition import PCA
        if ifpca:
            pca = PCA(n_components=3)
            pca.fit(mol.cart_coords)
            new_coords = pca.transform(mol.cart_coords)
            newmol = Molecule(coords=new_coords, species=mol.species)
        else:
            newmol = mol
        newmol.rotate_sites(axis=[1, 0, 0], theta=theta[0] / 180 * math.pi)
        newmol.rotate_sites(axis=[0, 1, 0], theta=theta[1] / 180 * math.pi)
        newmol.rotate_sites(axis=[0, 0, 1], theta=theta[2] / 180 * math.pi)
        if not spec_z:
            if frac:
                ztop = np.max(self.poscar.structure.frac_coords[:,2], axis=0)
                zz = ztop + distance/self.poscar.structure.lattice.c
                loc = [loc[0] + v[0], loc[1] + v[1], zz + v[2]]
            else:
                ztop = np.max(self.poscar.structure.cart_coords, axis=0)[2]
                loc[2]=distance+ztop
        else:
            loc = [loc[0] + v[0], loc[1] + v[1], loc[2] + v[2]]
        if frac:
            loc=self.f_to_c(loc)
        sufsite = AdsorbateSiteFinder(self.poscar.structure)
        ads_struct = sufsite.add_adsorbate(newmol, loc, reorient=True,translate=False)
        return ads_struct

    def add(self,element,local,frac=True):
        if frac:
            self.poscar.structure.append(Element(element),local)
        else:
            self.poscar.structure.append(Element(element), self.f_to_c(local))

    def find_adsorb_site(self,how='all',distance=2.5,frac=False):
        sufsite = AdsorbateSiteFinder(self.poscar.structure)
        ads_sites = sufsite.find_adsorption_sites(distance)
        result=ads_sites[how]
        if frac:
            for index,item in enumerate(result):
                result[index]=self.c_to_f(item)
        return result

    def copy(self):
        return easyjob(self.poscar.structure.copy(),self.foldname)

    def remove(self,i):
        s=True
        if str(i).isdigit():
            del self.poscar.structure[i]
        else:
            while s:
                for index,item in enumerate(self.poscar.structure):
                    if item.species.elements[0]==Element(i):
                        del self.poscar.structure[index]
                        break
                    s=False
        self.refresh_potcar()

    def replace(self,atom1,atom2):
        if str(atom1).isdigit():
            atom1=int(atom1)
        self.poscar.structure[atom1]=atom2

    def apply_strain(self,strain):
        self.poscar.structure.apply_strain(strain)

    def split_from_element(self,atomselect=['Na']):
        bottop=-1000
        atomselect=[Element(i) for i in atomselect]
        for site in self.poscar.structure:
            if site.specie in atomselect:
                bottop=max(site.z,bottop)
        sufmodel=self.poscar.structure.copy()
        sufdl=[]
        molmodel=self.poscar.structure.copy()
        moldl=[]
        for index,item in enumerate(sufmodel):
            if item.z > bottop+0.001:
                sufdl.append(index)
        for index,item in enumerate(molmodel):
            if item.z < bottop+0.001:
                moldl.append(index)
        sufmodel.remove_sites(sufdl)
        molmodel.remove_sites(moldl)
        suf = easyjob(sufmodel, os.path.join(self.path,'suf'))
        mol = easyjob(molmodel, os.path.join(self.path, 'mol'))
        suf.getincar('scf')
        mol.getincar('scf')
        return suf,mol

    def split_from_z(self,z,frac=True):
        """
        介绍：按照z坐标分割结构成为两个新结构
        参数：z：z是z方向的坐标  frac：指定是否使用分数坐标形式,默认为使用
        返回值：suf是z坐标以下的easyjob对象，mol是z坐标以上的easyjob对象
        """
        sufmodel = self.poscar.structure.copy()
        molmodel = self.poscar.structure.copy()
        sufdl = []
        moldl = []
        for index,item in enumerate(sufmodel):
            if frac:
                if self.c_to_f([item.x,item.y,item.z])[2] > z+0.0001:
                    sufdl.append(index)
            else:
                if item.z > z + 0.01:
                    moldl.append(index)
        for index,item in enumerate(molmodel):
            if frac:
                if self.c_to_f([item.x,item.y,item.z])[2] < z+0.0001:
                    moldl.append(index)
            else:
                if item.z < z+0.01:
                    moldl.append(index)
        sufmodel.remove_sites(sufdl)
        molmodel.remove_sites(moldl)
        suf = easyjob(sufmodel, os.path.join(self.path,'suf'))
        mol = easyjob(molmodel, os.path.join(self.path, 'mol'))
        suf.get_incar('scf')
        mol.get_incar('scf')
        return suf,mol

    def split_from_index(self,v):
        sufmodel = self.poscar.structure.copy()
        molmodel = self.poscar.structure.copy()
        sufdl = np.array([v])
        moldl = np.array([i for i in np.arange(len(self.poscar.structure)-1) if i not in sufdl])
        sufmodel.remove_sites(sufdl)
        molmodel.remove_sites(moldl)
        suf = easyjob(sufmodel, os.path.join(self.path, 'suf'))
        mol = easyjob(molmodel, os.path.join(self.path, 'mol'))
        suf.get_incar('scf')
        mol.get_incar('scf')
        return suf,mol

    def write_job(self):
        self.write_potcar()
        self.write_poscar()
        self.write_incar()
        self.write_kpoints()
        self.get_script()
        os.chdir(self.spath)

    def write_potcar(self):
        self.get_potcar()
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.chdir(self.path)
        self.refresh_potcar()
        with open('POTCAR', 'w') as f:
            f.write(self.potcar_file)
        os.chdir(self.rootpath)

    def write_kpoints(self):
        self.get_kpoints()
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.chdir(self.path)
        with open('KPOINTS', 'w') as f:
            f.write(self.kpoint_file)
        os.chdir(self.rootpath)

    def write_incar(self):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.chdir(self.path)
        self.incar.tofile('INCAR')
        os.chdir(self.rootpath)

    def write_poscar(self, direct=True):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.chdir(self.path)
        self.poscar.write_file('POSCAR', direct=direct)
        os.chdir(self.rootpath)

if __name__ == '__main__':
    rootpath = os.path.split(os.path.realpath(__file__))[0]
    print(rootpath)
