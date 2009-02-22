import math

import numpy

from dftatom.solvers import (broyden1, broyden2, broyden3, broyden_generalized,
        broyden_modified, broyden1_modified, linearmixing, norm, vackar,
        anderson, anderson2, excitingmixing)

def F(x):
    def p3(y):
        return float(y.T*y)*y
    try:
        x=tuple(x.flat)
    except:
        pass
    x=numpy.matrix(x).T

    d=numpy.matrix(numpy.diag([3,2,1.5,1,0.5]))
    c=0.01
    f=-d*x-c*p3(x)

    return tuple(f.flat)

xin=[1,1,1,1,1]


def test_broyden1():
    x= broyden1(F,xin,iter=11,alpha=1)
    assert norm(x)<1e-9
    assert norm(F(x))<1e-9

def test_broyden2():
    x= broyden2(F,xin,iter=12,alpha=1)
    assert norm(x)<1e-9
    assert norm(F(x))<1e-9

def test_broyden3():
    x= broyden3(F,xin,iter=12,alpha=1)
    assert norm(x)<1e-9
    assert norm(F(x))<1e-9

def test_exciting():
    x= excitingmixing(F,xin,iter=20,alpha=0.5)
    assert norm(x)<1e-5
    assert norm(F(x))<1e-5

def test_anderson():
    x= anderson(F,xin,iter=12,alpha=0.03,M=5)
    assert norm(x)<0.33

def test_anderson2():
    x= anderson2(F,xin,iter=12,alpha=0.6,M=5)
    assert norm(x)<0.2

def test_broydengeneralized():
    x= broyden_generalized(F,xin,iter=60,alpha=0.5,M=0)
    assert norm(x)<1e-7
    assert norm(F(x))<1e-7
    x= broyden_generalized(F,xin,iter=61,alpha=0.1,M=1)
    assert norm(x)<2e-4
    assert norm(F(x))<2e-4

def xtest_broydenmodified():
    x= broyden_modified(F,xin,iter=12,alpha=1)
    assert norm(x)<1e-9
    assert norm(F(x))<1e-9

def test_broyden1modified():
    x= broyden1_modified(F,xin,iter=35,alpha=1)
    assert norm(x)<1e-9
    assert norm(F(x))<1e-9

def test_vackar():
    x= vackar(F,xin,iter=11,alpha=1)
    assert norm(x)<1e-9
    assert norm(F(x))<1e-9

def test_linearmixing():
    x= linearmixing(F,xin,iter=60,alpha=0.5)
    assert norm(x)<1e-7
    assert norm(F(x))<1e-7
