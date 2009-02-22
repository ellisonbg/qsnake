import os
from tempfile import mkdtemp

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
        log = file(tmp+"/log").readlines()
        for l in log:
            if l.find("etotal") != -1:
                print l
