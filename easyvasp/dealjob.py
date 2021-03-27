from pymatgen.io.vasp.outputs import *
from pymatgen.io.vasp.inputs import Poscar, Potcar
import pandas as pd
import numpy as np
import sys, os
import shutil
from scipy import constants as con
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType
from pymatgen.electronic_structure.plotter import BSDOSPlotter, BSPlotter, BSPlotterProjected, DosPlotter, CohpPlotter
from pymatgen.electronic_structure.cohp import CompleteCohp
import matplotlib.pyplot

class dealjob(object):
    def __init__(self, path):
        self.path = path
        self.name = os.path.split(path)[-1]
        self.poscar = None
        self.potcar = None
        self.potcartype = []
        self.flist = []
        self.filist = []
        self.load()

    def __getattr__(self, name):
        if name == 'energy' or name == 'E':
            return self.get_energy()
        if name == 'fermi':
            return self.get_fermi_energy()
        if name == 'vbm':
            return self.get_gap()[2]
        if name == 'cbm':
            return self.get_gap()[1]
        if name == 'gap':
            return self.get_gap()[0]
        if name == 'zpe':
            return self.get_zpe()[0]
        if name == 'ifreq':
            return self.get_zpe()[1]
        if name == 'mag':
            return self.get_mag()
        if name == 'work_function' or name == 'wf' or name == 'workfunction':
            return self.get_work_function()[1]
        if name == 'vacuum_level' or name == 'vl' or name == 'vacuumlevel':
            return self.get_work_function()[2]
        if name == 'H':
            return self.get_therm(298)[1]
        if name == 'G':
            return self.get_therm(298)[2]
        if name == 'cvdt':
            return self.get_therm(298)[3]
        if name == 'TS' or name == 'ts':
            return self.get_therm(298)[5]
        if name == 'time':
            return self.get_time()
        if name == 'nelect':
            self.get_nelect()
            return self.nelect
        if name == 'POSCAR' or 'poscar':
            return self.poscar.structure
        if name == 'd_band_center':
            return self.get_d_band_center()
        return None

    def load(self):
        self.poscar = Poscar.from_file(os.path.join(self.path, 'CONTCAR'))

    def get_potcar(self):
        self.potcar = ''
        file = os.path.join(self.path, 'POTCAR')
        with open(file, 'r') as f:
            self.potcar = f.read()

    def get_nelect(self):
        self.get_potcar()
        nelect = 0
        k = 0
        self.eledata = []
        potlist = self.potcar.split('\n')
        for index, item in enumerate(potlist):
            if len(item) > 2:
                if item.split()[0] == 'PAW_PBE':
                    nelect += float(potlist[index + 1]) * self.poscar.natoms[k]
                    self.eledata.append(float(potlist[index + 1]))
                    self.potcartype.append(item.split()[1].split('_')[0])
                    k += 1
        self.nelect = nelect

    def grep(self, file, s):
        filepath = os.path.join(self.path, file)
        f = open(filepath)
        content = f.read()
        str = s + ' *' + '(.*?)' + ' ' + '.*'
        list = re.findall(str, content)
        f.close()
        return list

    def grep_tail(self, file, s):
        filepath = os.path.join(self.path, file)
        f = open(filepath)
        content = f.read()
        str = s + ' *' + '(.*?)' + '\\n'
        list = re.findall(str, content)
        f.close()
        return list

    def get_energy(self):
        strlist = self.grep('OUTCAR', 'energy  without entropy=')
        return float(strlist[-1])

    def get_fermi_energy(self):
        strlist = self.grep('OUTCAR', 'E-fermi :')
        return float(strlist[-1])

    def get_energy_list(self):
        strlist = self.grep('OUTCAR', 'energy  without entropy=')
        return [float(i) for i in strlist]

    def get_time(self):
        strlist = self.grep_tail('OUTCAR', 'User time \(sec\):')
        return float(strlist[0])

    def get_mag(self):
        filepath = os.path.join(self.path, 'OSZICAR')
        f = open(filepath)
        content = f.read()
        str = 'mag=.*'
        list = re.findall(str, content)
        f.close()
        if len(list) != 0:
            return round(float(list[-1].split()[-1]), 2)
        else:
            return 0

    def get_freq(self):
        if len(self.flist)==0:
            filepath = os.path.join(self.path, 'OUTCAR')
            with open(filepath, 'r') as f:
                for line in f.readlines():
                    sp = line.split()
                    if len(sp) >= 4 and sp[-1] == 'meV':
                        freq = round(float(sp[-2]), 4)
                        if sp[1] == 'f':
                            self.flist.append(freq)
                        else:
                            self.filist.append(freq)
    def get_zpe(self):
        self.get_freq()
        if len(self.filist) == 0:
            zpe = round(np.sum(np.array(self.flist)) / 2000, 4)
            fp = 0
        else:
            zpe = round(np.sum(np.array(self.flist)) / 2000, 4)
            fp = round(np.sum(np.array(self.filist)) / 2000, 4)
        return zpe, fp

    def therm_info(self, temp):
        a = self.get_therm(temp)
        print(
            f'E is {a[0]:.2f} eV, H is {a[1]:.2f},G is {a[2]:.2f}, cvdt is {a[3]:.2f} , zpe is {a[4]:.2f} ,TS is {a[5]:.2f}')

    def get_therm(self, temp):
        temp = int(temp)
        E = self.energy
        zpe, fp = self.get_zpe()
        cvdt = self.getre(temp)
        TS = self.get_tol_ts(temp)
        H = round(E + zpe + cvdt, 4)
        G = round(H - TS, 4)
        return E, H, G, cvdt, zpe, TS

    def getre(self, temp=300):
        tolre = 0
        for i in self.flist:
            tolre += self.getcvdt(i, temp)
        return tolre

    def getcvdt(self, fre, Temp):
        cvdt = 0
        for i in range(1, Temp*100, 1):
            cv = self.getcv(fre, i/100)
            cvdt += cv * 0.01
        cvdt = cvdt / 1000 / 96.485  # 单位eV
        return cvdt

    def getcv(self, fre, Temp):
        fre=max(50,fre)
        fre = float(fre) * 100
        sitav = fre * con.h * con.c / con.k
        cvf = sitav / Temp * (math.exp(-1 * sitav / Temp) / (1 - math.exp(-1 * sitav / Temp)))
        Cv = con.R * cvf * cvf
        return Cv

    def get_tol_ts(self, temp=300):
        tolts = 0
        for i in self.flist:
            i=max(50,i)
            tolts += self.get_ts(i, temp)
        return tolts

    def get_ts(self, nu, temp=300):
        nu = float(nu) * 100
        temp = float(temp)
        h_p = con.h
        k_b = con.k
        R_gas = con.R
        l_s = con.c
        beta = 1 / (k_b * temp)
        x_i = h_p * nu * l_s * beta
        if x_i <= 100:
            pf_l = x_i / (math.exp(x_i) - 1)
            pf_r = math.log(1 - math.exp(-x_i))
            pf = pf_l - pf_r
            entropy = R_gas * pf
            TS = entropy * temp / 1000 / 96.485 / 3
        else:
            TS = 0
        return TS

    def get_bader(self):
        self.get_nelect()
        eledict = dict(zip(self.potcartype, self.eledata))
        elelist = []
        elenum = []
        for i in self.poscar.structure:
            element = str(i.species)[:-1]
            elelist.append(element)
            elenum.append(eledict[element])
        chglist = []
        with open(os.path.join(self.path, 'ACF.dat')) as f:
            data = f.readlines()
            for i in data:
                if len(i.split()) == 7:
                    chglist.append(i.split()[4])
        chglist = np.array(chglist).astype('float64')
        elenum = np.array(elenum).astype('float64')
        valchg = np.around(elenum - chglist, 5)
        tab = np.vstack((np.array(elelist), np.array(valchg))).T
        return tab

    def get_diffchg(self, sufpath, molpath):
        tolpath = self.path
        totalchg = Chgcar.from_file(os.path.join(self.path, 'CHGCAR'))
        sufchg = Chgcar.from_file(os.path.join(sufpath, 'CHGCAR'))
        molchg = Chgcar.from_file(os.path.join(molpath, 'CHGCAR'))
        diffchg = totalchg - sufchg - molchg
        diffchg.write_file(os.path.join(tolpath, 'DIFFCHGCAR'))
        return diffchg

    def get_edos(self, stack=True, xlim=[-10, 10], sigma=0.15):
        dos_vasprun = Vasprun(os.path.join(self.path, "vasprun.xml"), parse_projected_eigen=True)
        dos_data = dos_vasprun.complete_dos
        plotter = DosPlotter(stack=stack, sigma=sigma)
        element_dos = dos_data.get_element_dos()
        plotter.add_dos_dict(element_dos)
        plotter.save_plot(f'{self.path}-edos.png', img_format=u'png', xlim=xlim, ylim=[0, 100])

    def get_sdos(self, dos_dict, stack=True, xlim=[-10, 10], ylim=None, sigma=0.15):
        dos_vasprun = Vasprun(os.path.join(self.path, "vasprun.xml"), parse_projected_eigen=True)
        dos_data = dos_vasprun.complete_dos
        plotter = DosPlotter(stack=stack, sigma=sigma)
        edos = dos_data.get_element_dos()
        for key in dos_dict:
            if type(dos_dict[key]) is list:
                for index, item in enumerate(dos_dict[key]):
                    if index == 0:
                        temp = dos_data.get_site_dos(dos_data.structure[item])
                    else:
                        temp = temp + dos_data.get_site_dos(dos_data.structure[item])
                plotter.add_dos(key, temp)
            elif type(dos_dict[key]) is str:
                temp = edos[Element(dos_dict[key])]
                plotter.add_dos(key, temp)
            elif type(dos_dict[key]) is int:
                temp = dos_data.get_site_dos(dos_data.structure[dos_dict[key]])
                plotter.add_dos(key, temp)
            else:
                print('Warning: invalid value')

        if ylim is not None:
            plotter.save_plot(f'{self.path}-dos.png', img_format=u'png', xlim=xlim, ylim=ylim)
        else:
            plotter.save_plot(f'{self.path}-dos.png', img_format=u'png', xlim=xlim)

    def get_odos(self, stack=True, xlim=[-10, 10], ylim=[0, 50], sigma=0.15):
        dos_vasprun = Vasprun(os.path.join(self.path, "vasprun.xml"), parse_projected_eigen=True)
        dos_data = dos_vasprun.complete_dos
        plotter = DosPlotter(stack=stack, sigma=sigma)
        o_dos = dos_data.get_spd_dos()
        plotter.add_dos_dict(o_dos)
        plotter.save_plot(f'{self.path}-odos.png', img_format=u'png', xlim=xlim, ylim=ylim)

    def get_gap(self):
        dos_vasprun = Vasprun(os.path.join(self.path, "vasprun.xml"), parse_projected_eigen=True)
        dos_data = dos_vasprun.complete_dos
        gap = dos_data.get_gap()
        cbm, vbm = dos_data.get_cbm_vbm()
        efermi = dos_data.efermi
        return (gap, cbm, vbm ,efermi)

    def get_cohp(self, xlim=[-4, 4], ylim=[-10, 6]):
        cohp = CompleteCohp.from_file('LOBSTER', filename=os.path.join(self.path, "COHPCAR.lobster"),
                                      structure_file=os.path.join(self.path, "CONTCAR"))
        plotter = CohpPlotter()
        plotter.add_cohp(self.name, cohp)
        plotter.save_plot(f'{self.path}.png', xlim=xlim, ylim=ylim, img_format='png')

    def get_work_function(self, axis='z'):
        axisdict = {'x': 0, 'y': 1, 'z': 2}
        locpot = Locpot.from_file(os.path.join(self.path, "LOCPOT"))
        line = locpot.get_average_along_axis(axisdict[axis])
        vacuum_level = np.max(line)
        ax = np.argmax(line)
        klist = line[ax - 10:ax + 10]
        if np.var(klist) >= 0.1:
            print('There may be an error in the calculation, please check the result carefully')
        work_function = vacuum_level - self.fermi
        return line, work_function, vacuum_level

    def get_d_band_center(self, spin=None, atom=None):
        dos_vasprun = Vasprun(os.path.join(self.path, "vasprun.xml"), parse_projected_eigen=True)
        dos_data = dos_vasprun.complete_dos
        if not atom:
            o_dos = dos_data.get_spd_dos()
        else:
            o_dos = dos_data.get_site_spd_dos(self.poscar.structure[atom])
        ddos = o_dos[OrbitalType.d]
        spindict = {'up': Spin.up, 'down': Spin.down}
        if not spin:
            dos = ddos.get_densities(spin=None)
        else:
            dos = ddos.get_densities(spin=spindict[spin])
        e = ddos.energies
        dsum = 0
        esum = 0
        for index in range(len(e) - 1):
            dx = e[index + 1] - e[index]
            dsum += e[index] * dos[index] * dx
            esum += dx * dos[index]
        ed = dsum / esum

        dsum = 0
        esum = 0
        for index in range(len(e) - 1):
            dx = e[index + 1] - e[index]
            dsum += (e[index] - ed) * dos[index] * dx
            esum += dx * dos[index]
        Wd = dsum / esum
        d_center = ed + Wd / 2 - self.fermi
        return round(d_center, 3)
