from dftatom.atom import atom, get_config
from dftatom.radial import create_log_grid2

def check(a, b, eps = 1.0):
    assert len(a) == len(b)
    for x,y in zip(a,b):
        assert x[:4]==y[:4]
        #print abs(x[4]-y[4]), x[4], y[4]
        assert abs(x[4]-y[4]) < eps

def test_boron_m1():
    "a test for wave4h_original, nonrelativistic atom"
    state = atom(5, iter=10, relat=-1)
    correct_state = [
            (1, 0, 1, 2, -6.564347), 
            (2, 0, 1, 2, -0.344701), 
            (2, 1, 1, 1, -0.136603)
            ]

    check(get_config(state), correct_state, 0.004)

def test_boron0():
    "a test for a schrodinger solver, nonrelativistic atom"
    state = atom(5, iter=10, relat=0)
    correct_state = [
            (1, 0, 1, 2, -6.564347), 
            (2, 0, 1, 2, -0.344701), 
            (2, 1, 1, 1, -0.136603)
            ]

    check(get_config(state), correct_state, 0.0002)

def test_boron2():
    "a test for a dirac solver, nonrelativistic atom"
    state = atom(5, iter=10, relat=2)
    correct_state = [
            (1, 0, 1, 2, -6.562952), 
            (2, 0, 1, 2, -0.344764), 
            (2, 1, -1, 1, -0.136646)
            ]

    check(get_config(state), correct_state, 0.002)

def test_lead2():
    "a test for a dirac solver, relativistic atom"
    state = atom(82, iter=10, relat=2)
    correct_state = [
            (1, 0, 1, 2, -3209.559491),
            (2, 0, 1, 2, -574.220773),
            (2, 1, -1, 6, -551.761782),
            (2, 1, 1, 6, -472.409787),
            (3, 0, 1, 2, -137.902276),
            (3, 1, -1, 6, -127.717013),
            (3, 1, 1, 6, -109.992083),
            (3, 2, -1, 10, -93.196164),
            (3, 2, 1, 10, -89.401978),
            (4, 0, 1, 2, -31.188199),
            (4, 1, -1, 6, -26.770858),
            (4, 1, 1, 6, -22.420348),
            (4, 2, -1, 10, -15.202782),
            (4, 2, 1, 10, -14.386516),
            (5, 0, 1, 2, -5.264036),
            (4, 3, -1, 14, -4.998502),
            (4, 3, 1, 14, -4.813673),
            (5, 1, -1, 6, -3.748596),
            (5, 1, 1, 6, -2.927281),
            (5, 2, -1, 10, -0.839121),
            (5, 2, 1, 10, -0.743864),
            (6, 0, 1, 2, -0.448677),
            (6, 1, -1, 2, -0.176692),
            ]

    check(get_config(state), correct_state, 0.16)

def test_lead_m1():
    "a test for wave4h_original, relativistic atom"
    state = atom(82, iter=10, relat=-1)
    correct_state = [
            (1, 0, 1, 2, -2901.078061),
            (2, 0, 1, 2, -488.843335),
            (2, 1, 1, 6, -470.877785),
            (3, 0, 1, 2, -116.526852),
            (3, 1, 1, 6, -107.950391),
            (3, 2, 1, 10, -91.889924),
            (4, 0, 1, 2, -25.753330),
            (4, 1, 1, 6, -21.990564),
            (4, 2, 1, 10, -15.030026),
            (4, 3, 1, 14, -5.592532),
            (5, 0, 1, 2, -4.206797),
            (5, 1, 1, 6, -2.941657),
            (5, 2, 1, 10, -0.902393),
            (6, 0, 1, 2, -0.357187),
            (6, 1, 1, 2, -0.141831),
            ]

    check(get_config(state), correct_state, 0.02)

def test_boron0_NIST():
    """a test for a schrodinger solver, nonrelativistic atom
    But using more than 10000 grid points and a NIST grid.
    """
    r = create_log_grid2(rmin=1./(1600*5), rmax = 20., N = 10010)
    state = atom(5, iter=20, relat=0, grid = r,eps=1e-10)
    correct_state = [
            (1, 0, 1, 2, -6.564347), 
            (2, 0, 1, 2, -0.344701), 
            (2, 1, 1, 1, -0.136603)
            ]

    check(get_config(state), correct_state, 0.99*1e-6)

def test_pb0_NIST():
    """a test for a schrodinger solver, relativistic atom
    """
    r = create_log_grid2(rmin=1./(1600*202), rmax = 20., N = 10000)
    state = atom(82, iter=20, relat=0, grid = r, eps = 1e-10)
    correct_state = [
            (1, 0, 1, 2, -2901.078061),
            (2, 0, 1, 2, -488.843335),
            (2, 1, 1, 6, -470.877785),
            (3, 0, 1, 2, -116.526852),
            (3, 1, 1, 6, -107.950391),
            (3, 2, 1, 10, -91.889924),
            (4, 0, 1, 2, -25.753330),
            (4, 1, 1, 6, -21.990564),
            (4, 2, 1, 10, -15.030026),
            (4, 3, 1, 14, -5.592532),
            (5, 0, 1, 2, -4.206797),
            (5, 1, 1, 6, -2.941657),
            (5, 2, 1, 10, -0.902393),
            (6, 0, 1, 2, -0.357187),
            (6, 1, 1, 2, -0.141831),
            ]

    check(get_config(state), correct_state, 0.99*1e-6)
