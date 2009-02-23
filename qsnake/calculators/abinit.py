from struct import calcsize, unpack
from collections import namedtuple

import os
from tempfile import mkdtemp

from numpy import array, reshape

class Abinit(object):

    def __init__(self, sage_path=None, atoms=None):
        if atoms is not None:
            self.init(atoms)
        self._sage_path = sage_path

    def init(self, atoms):
        self._atoms = atoms
        print "Abinit initialized with the configuration:"
        print atoms

    def calculate(self):
        tmp = mkdtemp()
        files = file(tmp+"/a.files", "w")
        files.write("a.in\na.out\nai\nao\ntmp\n01h.pspgth")
        files.close()
        file(tmp+"/01h.pspgth", "w").write("""\
Goedecker-Teter-Hutter  Wed May  8 14:27:44 EDT 1996
1   1   960508                     zatom,zion,pspdat
2   1   0    0    2001    0.       pspcod,pspxc,lmax,lloc,mmax,r2well
0.2000000 -4.0663326  0.6778322 0 0     rloc, c1, c2, c3, c4
0 0 0                              rs, h1s, h2s
0 0                                rp, h1p
  1.36 .2   0.6                    rcutoff, rloc
""")
        a_in = file(tmp+"/a.in", "w")
        a_in.write("""\
acell 10 10 10
ntypat 1
znucl 1
natom 2
typat 1 1
xcart
-0.7 0.0 0.0
0.7 0.0 0.0
ecut 10.0
nkpt 1
nstep 10
toldfe 1.0d-6
diemac 2.0
""")
        a_in.close()
        print "calling abinis in '%s'" % tmp
        r = os.system("cd "+tmp+"; "+self._sage_path+"/local/bin/abinis < a.files > log")
        if r != 0:
            raise RuntimeError("Abinis returned: %d" % r)

        density = self.parse_density(tmp+"/ao_DEN")
        log = file(tmp+"/log").readlines()
        for l in log:
            if l.find("etotal") != -1:
                print l

    def parse_density(self, filename):
        def read(f, fmt):
            return unpack(fmt, f.read(calcsize(fmt)))

        f = file(filename, "rb")
        codvsn = f.read(18)
        (headform, fform, bantot, date, intxc, ixc, natom, ngfft1, ngfft2,
                ngfft3, nkpt, nspden, nspinor, nsppol, nsym, npsp, ntypat,
                occopt, pertcase, usepaw) = read(f, "20i")
        ecut, ecutdg, ecutsm, ecut_eff = read(f, "4d")
        qptn = array(read(f, "3d"))
        rprimd = reshape(array(read(f, "9d")), (3, 3))
        stmbias, tphysel, tsmear = read(f, "3d")
        istwfk = array(read(f, "%di" % nkpt))
        nband = array(read(f, "%di" % nkpt*nsppol))
        npwarr = array(read(f, "%di" % nkpt))
        so_psp = array(read(f, "%di" % npsp))
        symafm = array(read(f, "%di" % nsym))
        symrel = reshape(array(read(f, "%di" % 3*3*nsym)), (3, 3, nsym))
        typat = array(read(f, "%di" % natom))
        kpt = reshape(array(read(f, "%dd" % 3*nkpt)), (3, nkpt))
        occ = array(read(f, "%dd" % bantot))
        tnons = reshape(array(read(f, "%dd" % 3*nsym)), (3, nsym))
        znucltypat = array(read(f, "%dd" % ntypat))
        wtk = array(read(f, "%dd" % nkpt))
        for ipsp in range(npsp):
            title = f.read(132)
            znuclpsp, zionpsp = read(f, "2d")
            pspso, pspdat, pspcod, pspxc, lmn_size = read(f, "5i")
        f.read(26)
        f.read(6)
        if usepaw == 0:
            residm, = read(f, "d")
            xred = reshape(array(read(f, "%dd" % 3*natom)), (3, natom))
            etotal, fermie = read(f, "2d")
            print residm, xred, etotal, fermie
        else:
            raise NotImplementedError("Only reading usepaw == 0 implemented so far")

        cplex = 1
        for ispden in range(nspden):
            rhor = array(read(f, "%dd" % cplex*ngfft1*ngfft2*ngfft3))

        density = []
        i = 0
        for i3 in range(ngfft3):
            for i2 in range(ngfft2):
                for i1 in range(ngfft1):
                    density.append([i1+1, i2+1, i3+1, rhor[i]])
                    i += 1
        return density
