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

        result = self.parse_density(tmp+"/ao_DEN")
        log = file(tmp+"/log").readlines()
        for l in log:
            if l.find("etotal") != -1:
                print l
        print result["header"]["etotal"]

    def parse_header(self, f):
        """
        Reads the header from an open file "f".

        It's probably fortran implementation dependent. This method assumes,
        that the following fortran statement:

        write(unit=header) codvsn,headform,fform

        Produces:

        <int: size of the block> <codvsn> <headform> <fform> <int: size of the
        block again, needs to match the first int>

        """
        def read(f, fmt):
            return unpack(fmt, f.read(calcsize(fmt)))

        # the size of the next section:
        head_size, = read(f, "i")
        assert head_size == calcsize("2i") + 6
        codvsn = f.read(6)
        headform, fform = read(f, "2i")
        head_size_end, = read(f, "i")
        assert head_size == head_size_end
        # the size of the next section:
        head_size, = read(f, "i")
        assert head_size == calcsize("18i4d3d9d3d")
        (bantot, date, intxc, ixc, natom, ngfft1, ngfft2, ngfft3, nkpt, nspden,
            nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase, usepaw) = \
                read(f, "18i")
        ecut, ecutdg, ecutsm, ecut_eff = read(f, "4d")
        qptn = array(read(f, "3d"))
        rprimd = reshape(array(read(f, "9d")), (3, 3))
        stmbias, tphysel, tsmear = read(f, "3d")
        head_size_end, = read(f, "i")
        assert head_size == head_size_end
        # the size of the next section:
        head_size, = read(f, "i")
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
        head_size_end, = read(f, "i")
        assert head_size == head_size_end
        # the size of the next block:
        psp_size, = read(f, "i")
        assert psp_size == calcsize("2d5i")+132
        for ipsp in range(npsp):
            title = (f.read(132)).strip()
            znuclpsp, zionpsp = read(f, "2d")
            pspso, pspdat, pspcod, pspxc, lmn_size = read(f, "5i")
        psp_size_end, = read(f, "i")
        assert psp_size == psp_size_end
        if usepaw == 0:
            final_size, = read(f, "i")
            # the size of the next block:
            assert final_size == calcsize("d"+"%d" % (3*natom) +"d2d")
            residm, = read(f, "d")
            xred = reshape(array(read(f, "%dd" % 3*natom)), (3, natom))
            etotal, fermie = read(f, "2d")
            final_size_end, = read(f, "i")
            assert final_size == final_size_end
        else:
            raise NotImplementedError("Only reading usepaw == 0 implemented so far")

        # Let's return some interesting information. We can of course return
        # more variables above if needed.
        result = {
                "etotal": etotal,
                "fermie": fermie,
                "ngfft": array([ngfft1, ngfft2, ngfft3]),
                "nspden": nspden
                }

        return result

    def parse_density(self, filename):
        """
        Reads the density file from abinit and decodes everything there is in
        there.

        It's probably fortran implementation dependent. This method assumes,
        that the following fortran statement:

        write(unit=header) codvsn,headform,fform

        Produces:

        <int: size of the block> <codvsn> <headform> <fform> <int: size of the
        block again, needs to match the first int>

        """
        def read(f, fmt):
            return unpack(fmt, f.read(calcsize(fmt)))

        f = file(filename, "rb")
        header = self.parse_header(f)

        rhor_size, = read(f, "i")
        cplex = 1
        ngfft1, ngfft2, ngfft3 = header["ngfft"]
        nspden = header["nspden"]
        assert rhor_size == cplex*ngfft1*ngfft2*ngfft3*calcsize("d")
        for ispden in range(nspden):
            rhor = array(read(f, "%dd" % cplex*ngfft1*ngfft2*ngfft3))
        rhor_size_end, = read(f, "i")
        assert rhor_size == rhor_size_end

        density = []
        i = 0
        for i3 in range(ngfft3):
            for i2 in range(ngfft2):
                for i1 in range(ngfft1):
                    density.append([i1+1, i2+1, i3+1, rhor[i]])
                    i += 1

        # Let's return some interesting information. We can of course return
        # more variables above if needed.
        result = {
                "density": density,
                "header": header,
                }

        return result
