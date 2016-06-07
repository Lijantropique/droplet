# -*- coding: utf-8 -*-
"""
droplet Test!
May 2016  - Canada
"""

import pytest
import droplet
import numpy as np

TOL=1E-6 # tolerance for calculations

def test_boundaries_extreme():
    #TODO: check limits <= values or just <
    """DOC"""
    # temp = Kelvin
    # pressure = MPa
    with pytest.raises(ValueError):
        droplet.Droplet(temperature=273.14999, pressure=100)
    with pytest.raises(ValueError):
        droplet.Droplet(temperature=1073.15001, pressure=100)
    with pytest.raises(ValueError):
        droplet.Droplet(temperature=2273.15001, pressure=50)

def test_boundaries_temperature():
    """DOC"""
    n = 10 # number of points in each dimension

    temperature = np.linspace(0, 273.14999, n)
    pressure = np.linspace(0, 100, n)
    tv, pv = np.meshgrid(temperature, pressure)
    for pair in zip(tv.ravel(),pv.ravel()):
        with pytest.raises(ValueError):
            droplet.Droplet(*pair)

    temperature = np.linspace(1073.15001, 2273.15, n)
    pressure = np.linspace(50.00001, 100, n)
    tv, pv = np.meshgrid(temperature, pressure)
    for pair in zip(tv.ravel(),pv.ravel()):
        with pytest.raises(ValueError):
            droplet.Droplet(*pair)

    temperature = np.linspace(2273.15001, 3000, n)
    pressure = np.linspace(0, 100, n)
    tv, pv = np.meshgrid(temperature, pressure)
    for pair in zip(tv.ravel(),pv.ravel()):
        with pytest.raises(ValueError):
            droplet.Droplet(*pair)

def test_boundaries_pressure():
    """DOC"""
    n = 10 # number of points in each dimension

    temperature = np.linspace(273.15, 1073.15, n)
    pressure = np.linspace(100.00001, 200, n)
    tv, pv = np.meshgrid(temperature, pressure)
    for pair in zip(tv.ravel(),pv.ravel()):
        with pytest.raises(ValueError):
            droplet.Droplet(*pair)

    temperature = np.linspace(1073.15001, 2273.15, n)
    pressure = np.linspace(50.00001, 100, n)
    tv, pv = np.meshgrid(temperature, pressure)
    for pair in zip(tv.ravel(),pv.ravel()):
        with pytest.raises(ValueError):
            droplet.Droplet(*pair)

def test_boundary2_3():
    drop = droplet.dispatcher['region23_boundary'](temperature=0.623150000E3)
    assert (drop[1] - 0.165291643E2) <= TOL
    drop = droplet.dispatcher['region23_boundary'](pressure=0.165291643E2)
    assert (drop[0] - 0.623150000E3) <= TOL

def test_region1_1():
    drop = droplet.Droplet(temperature=300, pressure=3)
    assert drop.region == 'Region1'
    assert (drop.v - 0.100215168E-2) <= TOL
    assert (drop.h - 0.115331273E3) <= TOL
    assert (drop.u - 0.112324818E3) <= TOL
    assert (drop.s - 0.392294792) <= TOL
    assert (drop.cp -0.417301218E1) <= TOL
    assert (drop.w - 0.150773921E4) <= TOL

def test_region1_2():
    drop = droplet.Droplet(temperature=300, pressure=80)
    assert drop.region == 'Region1'
    assert (drop.v - 0.971180894E-3) <= TOL
    assert (drop.h - 0.184142828E3) <= TOL
    assert (drop.u - 0.106448356E3) <= TOL
    assert (drop.s - 0.368563852) <= TOL
    assert (drop.cp -0.401008987E1) <= TOL
    assert (drop.w - 0.163469054E4) <= TOL

def test_region1_3():
    drop = droplet.Droplet(temperature=500, pressure=3)
    assert drop.region == 'Region1'
    assert (drop.v - 0.120241800E-2) <= TOL
    assert (drop.h - 0.975542239E3) <= TOL
    assert (drop.u - 0.971934985E3) <= TOL
    assert (drop.s - 0.258041912E1) <= TOL
    assert (drop.cp -0.465580682E1) <= TOL
    assert (drop.w - 0.124071337E4) <= TOL

def test_region2_1():
    drop = droplet.Droplet(temperature=300, pressure=0.0035)
    assert drop.region == 'Region2'
    assert (drop.v - 0.394913866E2) <= TOL
    assert (drop.h - 0.254991145E4) <= TOL
    assert (drop.u - 0.241169160E4) <= TOL
    assert (drop.s - 0.852238967E1) <= TOL
    assert (drop.cp -0.191300162E1) <= TOL
    assert (drop.w - 0.427920172E3) <= TOL

def test_region2_2():
    drop = droplet.Droplet(temperature=700, pressure=0.0035)
    assert drop.region == 'Region2'
    assert (drop.v - 0.923015898E2) <= TOL
    assert (drop.h - 0.333568375E4) <= TOL
    assert (drop.u - 0.301262819E4) <= TOL
    assert (drop.s - 0.101749996E2) <= TOL
    assert (drop.cp -0.208181274E1 ) <= TOL
    assert (drop.w - 0.644289068E3) <= TOL

def test_region2_3():
    drop = droplet.Droplet(temperature=700, pressure=30)
    assert drop.region == 'Region2'
    assert (drop.v - 0.542946619E-2) <= TOL
    assert (drop.h - 0.263149474E4) <= TOL
    assert (drop.u - 0.246861076E4) <= TOL
    assert (drop.s - 0.517540298E1) <= TOL
    assert (drop.cp -0.103505092E2) <= TOL
    assert (drop.w - 0.480386523E3) <= TOL

def test_region2Meta_1():
    drop = droplet.Droplet(temperature=450, pressure=1)
    assert drop.region == 'Region2-Meta'
    assert (drop.v - 0.192516540) <= TOL
    assert (drop.h - 0.276881115E4) <= TOL
    assert (drop.u - 0.257629461E4) <= TOL
    assert (drop.s - 0.656660377E1) <= TOL
    assert (drop.cp -0.276349265E1) <= TOL
    assert (drop.w - 0.498408101E3) <= TOL

def test_region2Meta_2():
    drop = droplet.Droplet(temperature=440, pressure=1)
    assert drop.region == 'Region2-Meta'
    assert (drop.v - 0.186212297) <= TOL
    assert (drop.h - 0.274015123E4) <= TOL
    assert (drop.u - 0.255393894E4) <= TOL
    assert (drop.s - 0.650218759E1) <= TOL
    assert (drop.cp -0.298166443E1) <= TOL
    assert (drop.w - 0.489363295E3) <= TOL

def test_region2Meta_3():
    drop = droplet.Droplet(temperature=450, pressure=1.5)
    assert drop.region == 'Region2-Meta'
    assert (drop.v - 0.121685206) <= TOL
    assert (drop.h - 0.272134539E4) <= TOL
    assert (drop.u - 0.253881758E4) <= TOL
    assert (drop.s - 0.629170440E1) <= TOL
    assert (drop.cp -0.362795578E1) <= TOL
    assert (drop.w - 0.481941819E3) <= TOL

#pag 32 Table 33 same input as normal, the program should solve the equation
# check python for goalseek.
def test_region3_1():
    drop = droplet.Droplet(temperature=650, pressure=0.255837018e2)
    assert drop.region == 'Region3'
    assert (drop.v - 1/500) <= TOL
    assert (drop.h - 0.186343019e4) <= TOL
    assert (drop.u - 0.181226279e4) <= TOL
    assert (drop.s - 0.405427273e1) <= TOL
    assert (drop.cp -0.138935717e2) <= TOL
    assert (drop.w - 0.502005554e3) <= TOL

def test_region3_2():
    drop = droplet.Droplet(temperature=650, pressure=0.222930643e2)
    assert drop.region == 'Region3'
    assert (drop.v - 1/200) <= TOL
    assert (drop.h - 0.237512401e4) <= TOL
    assert (drop.u - 0.226365868e4) <= TOL
    assert (drop.s - 0.485438792e1) <= TOL
    assert (drop.cp -0.446579342e2) <= TOL
    assert (drop.w - 0.383444594e3) <= TOL

def test_region3_3():
    drop = droplet.Droplet(temperature=750, pressure=0.783095639e2)
    assert drop.region == 'Region3'
    assert (drop.v - 1/500) <= TOL
    assert (drop.h - 0.225868845e4) <= TOL
    assert (drop.u - 0.210206932e4) <= TOL
    assert (drop.s - 0.446971906e1) <= TOL
    assert (drop.cp -0.634165359e1) <= TOL
    assert (drop.w - 0.760696041e3) <= TOL

def test_region4_1():
    drop = droplet.Droplet(temperature=300)
    assert drop.region == 'Region4'
    assert (drop.pressure - 0.353658941e-2) <= TOL
    drop = droplet(temperature=500)
    assert drop.region == 'Region4'
    assert (drop.pressure - 0.263889776e1) <= TOL
    drop = droplet(temperature=600)
    assert drop.region == 'Region4'
    assert (drop.pressure - 0.123443146e2) <= TOL

def test_region4_2():
    drop = droplet.Droplet(pressure=0.1)
    assert drop.region == 'Region4'
    assert (drop.temperature - 0.372755919e3) <= TOL
    drop = droplet(pressure=1)
    assert drop.region == 'Region4'
    assert (drop.temperature - 0.453035632e3) <= TOL
    drop = droplet(pressure=10)
    assert drop.region == 'Region4'
    assert (drop.temperature - 0.584149488e3) <= TOL

def test_region5_1():
    drop = droplet.Droplet(temperature=1500, pressure=0.5)
    assert drop.region == 'Region5'
    assert (drop.v - 0.138455090e1) <= TOL
    assert (drop.h - 0.521976855e4) <= TOL
    assert (drop.u - 0.452749310e4) <= TOL
    assert (drop.s - 0.965408875e1) <= TOL
    assert (drop.cp -0.261609445e1) <= TOL
    assert (drop.w - 0.917068690e3) <= TOL

def test_region5_2():
    drop = droplet.Droplet(temperature=1500, pressure=30.0)
    assert drop.region == 'Region5'
    assert (drop.v - 0.230761299e-1) <= TOL
    assert (drop.h - 0.516723514e4) <= TOL
    assert (drop.u - 0.447495124e4) <= TOL
    assert (drop.s - 0.772970133e1) <= TOL
    assert (drop.cp -0.272724317e1) <= TOL
    assert (drop.w - 0.928548002e3) <= TOL

def test_region5_3():
    drop = droplet.Droplet(temperature=2000, pressure=30.0)
    assert drop.region == 'Region5'
    assert (drop.v - 0.311385219e-1) <= TOL
    assert (drop.h - 0.657122604e4) <= TOL
    assert (drop.u - 0.563707038e4) <= TOL
    assert (drop.s - 0.853640523e1) <= TOL
    assert (drop.cp -0.288569882e1) <= TOL
    assert (drop.w - 0.106736948e4) <= TOL
