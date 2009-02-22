#! /usr/bin/env python

#import cgitb;cgitb.enable(format="text")
import math

from numpy import array, ones

import solvers
import radial
import utils

def totalenergy(r,V,Vext,Vh,Vxc,n,energies):
    Ts=0
    for e in energies:
        Ts+=e
    Ts-=integrate(mul2(r,mul2(r,mul2(V,n))),r)
    Epot=0
    Epot+=integrate(mul2(r,mul2(r,mul2(Vext,n))),r)
    Epot+=integrate(mul2(r,mul2(r,mul2(Vh,n))),r)
#    Epot+=integrate(mul2(r,mul2(r,mul2(Vxc,n))),r)
    return Ts+Epot

def KSiteration_potential(potential,r,Z,relat=2,eps=1e-5):
    v=array([-Z/x for x in r])
    u=v + potential
    density=radial.KS_construct_density(radial.KS_solve(u,r,Z,relat=relat,eps=eps),r,Z)
    Vxc=radial.getVxc(density,r,relat)
    Vh=radial.getVh(density,r,Z)
#    print totalenergy(r,u,v,Vh,Vxc,density,energies)
    return Vxc+Vh

def show(eigs, correct=None):
    spec={0:"s",1:"p",2:"d",3:"f",4:"g",5:"h"}
    c = get_config(eigs)
    if correct:
        assert len(c) == len(correct)
        print "Results:"
        print 'configuration: energy | NIST energy ("exact")'
        for X, Y in zip(c, correct):
            n,l,p,num,E = X
            if p==1: s="+"
            else: s="-"
            print "%d%s(%2d) j=l%s1/2: %.10g | %.10g"%(n,spec[l],num,s,E,Y[-1])
    else:
        for n,l,p,num,E in c:
            if p==1: s="+"
            else: s="-"
            print "%d%s(%2d) j=l%s1/2: %.10g"%(n,spec[l],num,s,E)

def get_config(eigs,n=None):
    Es,data=eigs
    if n == None:
        n = len(Es)
    done={}
    conf=[]
    for E,i in Es[:n]:
        y,n,l,p=data[i]
        if done.has_key((n,l,p)): 
            done[(n,l,p)]+=1
            continue
        done[(n,l,p)]=1
        conf.append((n,l,p,E))
    c = {}
    for d in done:
        n,l,p = d
        if not c.has_key((n,l)):
            c[(n,l)] = 0
        c[(n,l)]+=done[d]
    return [(n,l,p,c[n,l],E) for n,l,p,E in conf]

def atom(Z=82,alpha=0.3,iter=20,relat=2,rmin=0.00625,rmax=20.,inc=1.0247,
        grid=None,eps=1e-5):
    if grid == None:
        r = radial.create_log_grid(Z,rmin,rmax,inc)
    else:
        r = grid
#def atom(Z=82,alpha=0.3,iter=20,relat=2,rmin=0.00625,rmax=20.,k0=1.0247):
#    r=radial.create_hyperbolic_grid(Z,rmin,rmax,k0)
    print "radial grid points: %d"%len(r)
    if relat == -1:
        r=radial.create_hyperbolic_grid(Z,k0=2000)

    #start with a uniform density
    density = ones(len(r),"d")
    density = Z/radial.integrate(density*r*r, r) * density
    Vxc=radial.getVxc(density,r,relat)
    Vh=radial.getVh(density,r,Z)
    pot = Vxc+Vh

    def F(pot):
        return KSiteration_potential(pot,r,Z,relat,eps=eps)-pot

    pot=solvers.excitingmixing(F,pot,iter=iter,alpha=alpha)
    #pot=solvers.broyden2(F,pot,iter=iter,alpha=alpha)
    #pot=solvers.broyden3(F,pot,iter=iter,alpha=alpha)

    v=array([-Z/x for x in r])
    u=v+pot
    x=radial.KS_solve(u,r,Z,relat=relat,eps=eps)
    return (x[0][:Z], x[1])
    #return (x[0][:Z], x[1], r)
