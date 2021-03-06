# -*- coding: utf-8 -*-
"""
droplet!
Water properties base on IAPWS Industrial 1997
http://www.iapws.org/relguide/IF97-Rev.html
"""
#TODO: Complete description including variables, and regions
#TODO: Create Objects: put the regions as static methods @staticmethod?
#TODO: Specify usage
#TODO: GUI

import os
import json
import math
from scipy.optimize import fsolve

# Reference Constants
R = 0.461526    # [kJ kg⁻¹ K⁻¹]
Tc = 647.096    # [K]
Pc = 22.064     # [MPa]
rhoc = 322      # [kg m⁻³]

folder_path = os.path.dirname(os.path.abspath(__file__)) + os.sep

def region1(temperature,pressure):
    P_star = 16.53
    T_star = 1386
    pi = pressure/P_star
    tao = T_star/temperature
    file_name = 'region1.json'
    with open(folder_path + 'data/'+ file_name) as data_file:
        data = json.load(data_file)

    # Gibbs Free Energy
    gibbs = lambda item: item['ni']*((7.1-pi)**item['Ii']) * ((tao-1.222)**item['Ji'])
    gamma = sum([gibbs(item) for item in data])

    # Partial derivative with respect to pi
    dgibbs_pi = lambda item: -item['ni']*item['Ii']*((7.1-pi)**(item['Ii']-1)) * ((tao-1.222)**item['Ji'])
    gamma_pi = sum([dgibbs_pi(item) for item in data])

    # Second partial derivative with respect to pi
    dgibbs_pi2 = lambda item: item['ni']*item['Ii']*(item['Ii']-1)*((7.1-pi)**(item['Ii']-2)) * ((tao-1.222)**item['Ji'])
    gamma_pi2 = sum([dgibbs_pi2(item) for item in data])

    # Partial derivative with respect to tao
    dgibbs_tao = lambda item: item['ni']*((7.1-pi)**item['Ii']) * item['Ji'] * ((tao-1.222)**(item['Ji']-1))
    gamma_tao = sum([dgibbs_tao(item) for item in data])

    # Second partial derivative with respect to tao
    dgibbs_dtao2 = lambda item: item['ni']*((7.1-pi)**item['Ii']) * item['Ji'] * (item['Ji']-1)* ((tao-1.222)**(item['Ji']-2))
    gamma_tao2 = sum([dgibbs_dtao2(item) for item in data])

    # Second partial derivative with respect to pi and tao
    dgibbs_pi_tao = lambda item: -item['ni']*(item['Ii'])*((7.1-pi)**(item['Ii']-1)) * item['Ji'] * ((tao-1.222)**(item['Ji']-1))
    gamma_pi_tao = sum([dgibbs_pi_tao(item) for item in data])

    # Specific Volume:
    # 0.001 is a factor to report the specific volume in [m³ kg⁻¹]
    v = gamma_pi*pi*R*temperature/pressure * 0.001

    # Specific Density:
    # [kg m⁻³]
    rho = 1/v

    # Specific Internal Energy:
    # [kJ kg⁻¹]
    u = (gamma_tao*tao-gamma_pi*pi)*R*temperature

    # Specific Entropy:
    # [kJ kg⁻¹ K⁻¹]
    s = (gamma_tao*tao-gamma)*R

    # Specific Enthalpy:
    # [kJ kg⁻¹]
    h = gamma_tao*tao*R*temperature

    # Specific Isobaric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    cp = -gamma_tao2*tao*tao*R

    # Specific Isochoric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    aux = -gamma_tao2*tao*tao + ((gamma_pi-tao*gamma_pi_tao)**2/gamma_pi2)
    cv = aux*R

    # Speed of Sound:
    # 31.6227766016838 is a factor to report the speed of sound in [m s⁻¹]
    aux = gamma_pi**2 / ((((gamma_pi - tao*gamma_pi_tao)**2) /(tao**2 * gamma_tao2))-gamma_pi2)
    w = (aux*R*temperature)**0.5 * 31.6227766016838

    return v, rho, u, s, h, cp, cv, w

def region2(temperature, pressure):
    P_star = 1
    T_star = 540
    pi = pressure/P_star
    tao = T_star/temperature

    # *********IDEAL*********
    file_name = 'region2_ideal.json'
    with open(folder_path + 'data/'+ file_name) as data_file:
        data = json.load(data_file)

    # Ideal Gibbs Free Energy
    gibbs = lambda item: item['n0']*(tao**item['J0'])
    gamma_ideal = math.log(pi) + sum([gibbs(item) for item in data])

    # Partial derivative with respect to pi
    gammai_pi = 1/pi

    # Second partial derivative with respect to pi
    gammai_pi2 = -1/(pi**2)

    # Partial derivative with respect to tao
    dgibbs_tao = lambda item: item['n0']*item['J0']*(tao**(item['J0']-1))
    gammai_tao = sum([dgibbs_tao(item) for item in data])

    #Second partial derivative with respect to tao
    dgibbs_dtao2 = lambda item: item['n0']*item['J0']*(item['J0']-1)*(tao**(item['J0']-2))
    gammai_tao2 = sum([dgibbs_dtao2(item) for item in data])

    # Second partial derivative with respect to pi and tao
    gammai_pi_tao = 0

    # *********RESIDUAL*********
    file_name = 'region2_residual.json'
    with open(folder_path + 'data/'+ file_name) as data_file:
        data = json.load(data_file)

    # Residual Gibbs Free Energy
    gibbs = lambda item: item['ni']*(pi**item['Ii'])*((tao-0.5)**item['Ji'])
    gamma_residual=  sum([gibbs(item) for item in data])

    # Partial derivative with respect to pi
    dgibbs_pi = lambda item: item['ni']*item['Ii']*(pi**(item['Ii']-1))*((tao-0.5)**item['Ji'])
    gammar_pi = sum([dgibbs_pi(item) for item in data])

    # Second partial derivative with respect to pi
    dgibbs_pi2 = lambda item: item['ni']*item['Ii']*(item['Ii']-1)*(pi**(item['Ii']-2))*((tao-0.5)**item['Ji'])
    gammar_pi2 = sum([dgibbs_pi2(item) for item in data])

    # Partial derivative with respect to tao
    dgibbs_tao = lambda item: item['ni']*(pi**item['Ii'])*(item['Ji'])*((tao-0.5)**(item['Ji']-1))
    gammar_tao = sum([dgibbs_tao(item) for item in data])

    # Second partial derivative with respect to tao
    dgibbs_dtao2 = lambda item: item['ni']*(pi**item['Ii'])*(item['Ji'])*(item['Ji']-1)*((tao-0.5)**(item['Ji']-2))
    gammar_tao2 = sum([dgibbs_dtao2(item) for item in data])

    #Second partial derivative with respect to pi and tao
    dgibbs_pi_tao = lambda item: item['ni']*(item['Ii'])*(pi**(item['Ii']-1))*(item['Ji'])*((tao-0.5)**(item['Ji']-1))
    gammar_pi_tao = sum([dgibbs_pi_tao(item) for item in data])

    # Specific Volume:
    # 0.001 is a factor to report the specific volume in [m³ kg⁻¹]
    aux = gammai_pi+gammar_pi
    v = pi*aux*R*temperature/pressure * 0.001

    # Specific Density:
    # [kg m⁻³]
    rho = 1/v

    # Specific Internal Energy:
    # [kJ kg⁻¹]
    aux = tao*(gammai_tao+gammar_tao)-pi*(gammai_pi+gammar_pi)
    u = aux*R*temperature

    # Specific Entropy:
    # [kJ kg⁻¹ K⁻¹]
    aux = tao*(gammai_tao + gammar_tao)-(gamma_ideal+gamma_residual)
    s = aux*R

    # Specific Enthalpy:
    # [kJ kg⁻¹]
    aux = gammai_tao + gammar_tao
    h = aux*tao*R*temperature

    # Specific Isobaric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    aux = gammai_tao2 + gammar_tao2
    cp = -aux*tao*tao*R

    # Specific Isochoric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    aux1 = (1 + pi*gammar_pi - tao*pi*gammar_pi_tao)**2
    aux2 = 1 - pi*pi*gammar_pi2
    aux = -(gammai_tao2 + gammar_tao2)*tao*tao - (aux1/aux2)
    cv = aux*R

    #Speed of Sound:
    # 31.6227766016838 is a factor to report the speed of sound in [m s⁻¹]
    aux1 = 1 + 2*pi*gammar_pi + (pi*gammar_pi)**2
    aux2 = (1 + pi*gammar_pi - tao*pi*gammar_pi_tao)**2
    aux3 = tao**2*(gammai_tao2 + gammar_tao2)
    aux = aux1 / (1 - pi**2*gammar_pi2 + (aux2/aux3))
    w = (aux*R*temperature)**0.5 * 31.6227766016838

    return v, rho, u, s, h, cp, cv, w

def region3(temperature, pressure):
    #TODO: Look for a better option to find density

    def check_pressure(density):
        delta = density/rhoc
        file_name = 'region3.json'
        with open(folder_path + 'data/'+ file_name) as data_file:
            data = json.load(data_file)

        # Helmholtz Free Energy
        helmholtz = lambda item: item['ni']*(delta**item['Ii'])*(tao**item['Ji'])
        phi = data[0]['ni']*math.log(delta) + sum([helmholtz(item) for item in data[1:]])

        # Partial derivative with respect to delta
        dhelmholtz_delta = lambda item: item['ni']*item['Ii']*(delta**(item['Ii']-1))*(tao**item['Ji'])
        phi_delta = data[0]['ni']/delta + sum([dhelmholtz_delta(item) for item in data[1:]])

        # Second partial derivative with respect to delta
        dhelmholtz_delta2 = lambda item: item['ni']*item['Ii']*(item['Ii']-1)*(delta**(item['Ii']-2))*(tao**item['Ji'])
        phi_delta2 = -data[0]['ni']/delta**2 + sum([dhelmholtz_delta2(item) for item in data[1:]])

        # Partial derivative with respect to tao
        dhelmholtz_tao = lambda item: item['ni']*(delta**item['Ii'])*item['Ji']*(tao**(item['Ji']-1))
        phi_tao = sum([dhelmholtz_tao(item) for item in data[1:]])

        # Second partial derivative with respect to tao
        dhelmholtz_tao2 = lambda item: item['ni']*(delta**item['Ii'])*item['Ji']*(item['Ji']-1)*(tao**(item['Ji']-2))
        phi_tao2 = sum([dhelmholtz_tao2(item) for item in data[1:]])

        # Second partial derivative with respect delta and tao
        dhelmholtz_delta_tao = lambda item: item['ni']*item['Ii']*(delta**(item['Ii']-1))*item['Ji']*(tao**(item['Ji']-1))
        phi_delta_tao = sum([dhelmholtz_delta_tao(item) for item in data[1:]])

        # Pressure Calculated:
        # 1000 is a factor to report the pressure in [MPa]
        pressure_calc = phi_delta*delta*density*R*temperature/1000

        return pressure_calc, phi, phi_delta, phi_delta2, phi_tao, phi_tao2, phi_delta_tao

    tao = Tc/temperature
    ########### Density LookUp ###########
    cont = 0
    density = 0.1
    while ((abs(check_pressure(density)[0]-pressure)>=1E-3) and (cont<10000)):
        density += 0.1
        cont +=1

    _, phi, phi_delta, phi_delta2, phi_tao, phi_tao2, phi_delta_tao = find_pressure(density)
    delta = density/rhoc
    ########### Density LookUp ###########

    # Specific Volume:
    v = 1/density

    # Specific Internal Energy:
    # [kJ kg⁻¹]
    u = tao*phi_tao*R*temperature

    # Specific Entropy:
    # [kJ kg⁻¹ K⁻¹]
    s = (tao*phi_tao-phi)*R

    # Specific Enthalpy:
    # [kJ kg^-1]
    aux = tao*phi_tao + delta*phi_delta
    h = aux*R*temperature

    # Specific Isobaric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    aux1 = (delta*phi_delta - delta*tao*phi_delta_tao)**2
    aux2 = (2*delta*phi_delta + delta*delta*phi_delta2)
    aux = -tao*tao*phi_tao2 + (aux1/aux2)
    cp = aux*R

    # Specific Isochoric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    aux = -tao*tao*phi_tao2
    cv = aux*R

    # Speed of Sound:
    # 31.6227766016838 is a factor to report the speed of sound in [m s⁻¹]
    aux1 = (delta*phi_delta - delta*tao*phi_delta_tao)**2
    aux2 = tao*tao*phi_tao2
    aux = 2*delta*phi_delta + delta*delta*phi_delta2 - (aux1/aux2)
    w = (aux*R*temperature)**0.5 * 31.6227766016838

    # Joule-Thomson coefficient:
    # 1000 is a factor to report the Joule-Thomson coefficient in [K MPa⁻¹]
    aux1 = -(delta*phi_delta + delta*delta*phi_delta2 + delta*tao*phi_delta_tao)
    aux2 = (delta*phi_delta - delta*tao*phi_delta_tao)**2
    aux3 = tao*tao*phi_tao2*(2*delta*phi_delta + delta*delta*phi_delta2)
    aux = aux1 /(aux2-aux3)
    jtc = aux/R/density*1000

    # Isothermal throttling coefficient:
    aux1 = delta*phi_delta - delta*tao*phi_delta_tao
    aux2 = 2*delta*phi_delta + delta*delta*phi_delta2
    aux = 1-(aux1/aux2)
    itc = aux/density

    # Isentropic temperature-pressure coefficient:
    # 1000 is a  factor to report the Isentropic temperature-pressure coefficient in [K MPa⁻¹]
    aux1 = delta*phi_delta - delta*tao*phi_delta_tao
    aux2 = (delta*phi_delta - delta*tao*phi_delta_tao)**2
    aux3 = tao*tao*phi_tao2*(2*delta*phi_delta + delta*delta*phi_delta2)
    aux = aux1 /(aux2-aux3)
    itpc= aux/R/density*1000

    return v, density, u, s, h, cp, cv, w

def region4(temperature=-1, pressure=-1):
    #TODO at this point either reg1/reg2 or reg3/reg2 shuould be used
    # and the properties should be the same (use average?)
    P_star = 1
    T_star = 1

    file_name = 'region4.json'
    with open(folder_path + 'data/'+ file_name) as data_file:
        data = json.load(data_file)

    if pressure==-1:
        theta = (temperature/T_star) + data[8]['ni']/((temperature/T_star)-data[9]['ni'])
        A = theta**2 + data[0]['ni']*theta + data[1]['ni']
        B = data[2]['ni']*(theta**2) + data[3]['ni']*theta + data[4]['ni']
        C = data[5]['ni']*(theta**2) + data[6]['ni']*theta + data[7]['ni']
        # Saturation Pressure
        # [MPa]
        pressure_sat = P_star * ((2*C / (-B + (B*B - 4*A*C)**0.5))**4)
        beta = (pressure_sat/P_star)**(1/4)
        temperature_sat = temperature

    if temperature==-1:
        beta = (pressure/P_star)**(1/4)
        E = beta**2 + data[2]['ni']*beta + data[5]['ni']
        F = data[0]['ni']*(beta**2) + data[3]['ni']*beta + data[6]['ni']
        G = data[1]['ni']*(beta**2) + data[4]['ni']*beta + data[7]['ni']
        D = 2*G / (-F-(F*F - 4*E*G)**0.5)
        # Saturation Temperature
        # [K]
        temperature_sat = T_star * (data[9]['ni'] + D - (((data[9]['ni']+D)**2 - 4*(data[8]['ni']+data[9]['ni']*D))**0.5))/2
        theta = (temperature_sat/T_star) + data[8]['ni']/((temperature/T_star)-data[9]['ni'])
        pressure_sat=pressure

    return temperature_sat, pressure_sat

def region5(temperature, pressure):
    P_star = 1
    T_star = 1000
    pi = pressure/P_star
    tao = T_star/temperature

    # *********IDEAL*********
    file_name = 'region5_ideal.json'
    with open(folder_path + 'data/'+ file_name) as data_file:
        data = json.load(data_file)

    # Ideal Gibbs Free Energy
    gibbs = lambda item: item['n0']*(tao**item['J0'])
    gamma_ideal = math.log(pi) + sum([gibbs(item) for item in data])

    # Partial derivative with respect to pi
    gammai_pi = 1/pi

    # Second partial derivative with respect to pi
    gammai_pi2 = -1/(pi**2)

    # Partial derivative with respect to tao
    dgibbs_tao = lambda item: item['n0']*item['J0']*(tao**(item['J0']-1))
    gammai_tao = sum([dgibbs_tao(item) for item in data])

    # Second partial derivative with respect to tao
    dgibbs_dtao2 = lambda item: item['n0']*item['J0']*(item['J0']-1)*(tao**(item['J0']-2))
    gammai_tao2 = sum([dgibbs_dtao2(item) for item in data])

    # Second partial derivative with respect to pi and tao
    gamma_pi_tao = 0

    # *********RESIDUAL*********
    file_name = 'region5_residual.json'
    with open(folder_path + 'data/'+ file_name) as data_file:
        data = json.load(data_file)

    # Residual Gibbs Free Energy
    gibbs = lambda item: item['ni']*(pi**item['Ii'])*(tao**item['Ji'])
    gamma_residual=  sum([gibbs(item) for item in data])

    # Partial derivative with respect to pi
    dgibbs_pi = lambda item: item['ni']*item['Ii']*(pi**(item['Ii']-1))*(tao**item['Ji'])
    gammar_pi = sum([dgibbs_pi(item) for item in data])

    # Second partial derivative with respect to pi
    dgibbs_pi2 = lambda item: item['ni']*item['Ii']*(item['Ii']-1)*(pi**(item['Ii']-2))*(tao**item['Ji'])
    gammar_pi2 = sum([dgibbs_pi2(item) for item in data])

    # Partial derivative with respect to tao
    dgibbs_tao = lambda item: item['ni']*(pi**item['Ii'])*(item['Ji'])*(tao**(item['Ji']-1))
    gammar_tao = sum([dgibbs_tao(item) for item in data])

    # Second partial derivative with respect to tao
    dgibbs_dtao2 = lambda item: item['ni']*(pi**item['Ii'])*(item['Ji'])*(item['Ji']-1)*(tao**(item['Ji']-2))
    gammar_tao2 = sum([dgibbs_dtao2(item) for item in data])

    # Second partial derivative with respect to pi and tao
    dgibbs_pi_tao = lambda item: item['ni']*(item['Ii'])*(pi**(item['Ii']-1))*(item['Ji'])*(tao**(item['Ji']-1))
    gammar_pi_tao = sum([dgibbs_pi_tao(item) for item in data])

    # Specific Volume:
    # 0.001 is a  factor to report the specific volume in [m³ kg⁻¹]
    aux = gammai_pi+gammar_pi
    v = pi*aux*R*temperature/pressure * 0.001

    # Specific Density:
    # [kg m⁻³]
    rho = 1/v

    # Specific Internal Energy:
    # [kJ kg⁻¹]
    aux = tao*(gammai_tao+gammar_tao)-pi*(gammai_pi+gammar_pi)
    u = aux*R*temperature

    # Specific Entropy:
    # [kJ kg⁻¹ K⁻¹]
    aux = tao*(gammai_tao + gammar_tao)-(gamma_ideal+gamma_residual)
    s = aux*R

    # Specific Enthalpy:
    # [kJ kg⁻¹]
    aux = gammai_tao + gammar_tao
    h = aux*tao*R*temperature

    # Specific Isobaric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    aux = gammai_tao2 + gammar_tao2
    cp = -aux*tao*tao*R

    # Specific Isochoric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    aux1 = (1 + pi*gammar_pi - tao*pi*gammar_pi_tao)**2
    aux2 = 1 - pi*pi*gammar_pi2
    aux = -(gammai_tao2 + gammar_tao2)*tao*tao - (aux1/aux2)
    cv = aux*R

    # Speed of Sound:
    # 31.6227766016838 is a  factor to report the speed of sound in [m s⁻¹]
    aux1 = 1 + 2*pi*gammar_pi + (pi*gammar_pi)**2
    aux2 = (1 + pi*gammar_pi - tao*pi*gammar_pi_tao)**2
    aux3 = tao**2*(gammai_tao2 + gammar_tao2)
    aux = aux1 / (1 - pi**2*gammar_pi2 + (aux2/aux3))
    w = (aux*R*temperature)**0.5 * 31.6227766016838

    return v, rho, u, s, h, cp, cv, w

def region2_meta(temperature, pressure):
    P_star = 1
    T_star = 540
    pi = pressure/P_star
    tao = T_star/temperature

    # *********IDEAL*********
    file_name = 'region2_meta_ideal.json'
    with open(folder_path + 'data/'+ file_name) as data_file:
        data = json.load(data_file)

    # Ideal Gibbs Free Energy
    gibbs = lambda item: item['n0']*(tao**item['J0'])
    gamma_ideal = math.log(pi) + sum([gibbs(item) for item in data])

    # Partial derivative with respect to pi
    gammai_pi = 1/pi

    # Second partial derivative with respect to pi
    gammai_pi2 = -1/(pi**2)

    # Partial derivative with respect to tao
    dgibbs_tao = lambda item: item['n0']*item['J0']*(tao**(item['J0']-1))
    gammai_tao = sum([dgibbs_tao(item) for item in data])

    # Second partial derivative with respect to tao
    dgibbs_dtao2 = lambda item: item['n0']*item['J0']*(item['J0']-1)*(tao**(item['J0']-2))
    gammai_tao2 = sum([dgibbs_dtao2(item) for item in data])

    # Second partial derivative with respect to pi and tao
    gammai_pi_tao = 0

    # *********RESIDUAL*********
    file_name = 'region2_meta_residual.json'
    with open(folder_path + 'data/'+ file_name) as data_file:
        data = json.load(data_file)

    # Residual Gibbs Free Energy
    gibbs = lambda item: item['ni']*(pi**item['Ii'])*((tao-0.5)**item['Ji'])
    gamma_residual=  sum([gibbs(item) for item in data])

    # Partial derivative with respect to pi
    dgibbs_pi = lambda item: item['ni']*item['Ii']*(pi**(item['Ii']-1))*((tao-0.5)**item['Ji'])
    gammar_pi = sum([dgibbs_pi(item) for item in data])

    # Second partial derivative with respect to pi
    dgibbs_pi2 = lambda item: item['ni']*item['Ii']*(item['Ii']-1)*(pi**(item['Ii']-2))*((tao-0.5)**item['Ji'])
    gammar_pi2 = sum([dgibbs_pi2(item) for item in data])

    # Partial derivative with respect to tao
    dgibbs_tao = lambda item: item['ni']*(pi**item['Ii'])*(item['Ji'])*((tao-0.5)**(item['Ji']-1))
    gammar_tao = sum([dgibbs_tao(item) for item in data])

    # Second partial derivative with respect to tao
    dgibbs_dtao2 = lambda item: item['ni']*(pi**item['Ii'])*(item['Ji'])*(item['Ji']-1)*((tao-0.5)**(item['Ji']-2))
    gammar_tao2 = sum([dgibbs_dtao2(item) for item in data])

    # Second partial derivative with respect to pi and tao
    dgibbs_pi_tao = lambda item: item['ni']*(item['Ii'])*(pi**(item['Ii']-1))*(item['Ji'])*((tao-0.5)**(item['Ji']-1))
    gammar_pi_tao = sum([dgibbs_pi_tao(item) for item in data])

    # *********PROPERTIES*********
    # Specific Volume:
    # 0.001 is a  factor to report the specific volume in [m³ kg⁻¹]
    aux = gammai_pi+gammar_pi
    v = pi*aux*R*temperature/pressure * 0.001

    # Specific Density:
    # [m⁻³ kg¹]
    rho = 1/v

    # Specific Internal Energy:
    # [kJ kg⁻¹]
    aux = tao*(gammai_tao+gammar_tao)-pi*(gammai_pi+gammar_pi)
    u = aux*R*temperature

    # Specific Entropy:
    # [kJ kg⁻¹ K⁻¹]
    aux = tao*(gammai_tao + gammar_tao)-(gamma_ideal+gamma_residual)
    s = aux*R

    # Specific Enthalpy:
    # [kJ kg⁻¹]
    aux = gammai_tao + gammar_tao
    h = aux*tao*R*temperature

    # Specific Isobaric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    aux = gammai_tao2 + gammar_tao2
    cp = -aux*tao*tao*R

    # Specific Isochoric Heat capacity:
    # [kJ kg⁻¹ K⁻¹]
    aux1 = (1 + pi*gammar_pi - tao*pi*gammar_pi_tao)**2
    aux2 = 1 - pi*pi*gammar_pi2
    aux = -(gammai_tao2 + gammar_tao2)*tao*tao - (aux1/aux2)
    cv = aux*R

    # Speed of Sound:
    # 31.6227766016838 is a  factor to report the speed of sound in [m s⁻¹]
    aux1 = 1 + 2*pi*gammar_pi + (pi*gammar_pi)**2
    aux2 = (1 + pi*gammar_pi - tao*pi*gammar_pi_tao)**2
    aux3 = tao**2*(gammai_tao2 + gammar_tao2)
    aux = aux1 / (1 - pi**2*gammar_pi2 + (aux2/aux3))
    w = (aux*R*temperature)**0.5 * 31.6227766016838

    return v, rho, u, s, h, cp, cv, w

def region23_boundary(temperature=-1, pressure=-1):
    #TODO at this point either reg1/reg2 or reg3/reg2 shuould be used
    # and the properties should be the same (use average?)
    P_star = 1
    T_star = 1

    file_name = 'region23_boundary.json'
    with open(folder_path + 'data/'+ file_name) as data_file:
        data = json.load(data_file)

    if pressure==-1:
        theta = (temperature/T_star)
        # Pressure
        # [MPa]
        pressure = P_star * (data[0]['ni'] + data[1]['ni'] *theta + data[2]['ni'] *theta**2 )

    if temperature==-1:
        pi = pressure/P_star
        # Temperature
        # [K]
        temperature = T_star * (data[3]['ni'] + (((pi - data[4]['ni'])/data[2]['ni'])**0.5))

    return temperature, pressure

dispatcher = {'Region 1': region1, 'Region 2': region2,
              'Region 3': region3, 'Region 4': region4,
              'Region 5': region5, 'region2_meta': region2_meta,
              'region23_boundary': region23_boundary}

#-------------------------------
class Droplet(object):
    # temperature : K
    # pressure: MPa
    def __init__(self, temperature=273.15, pressure=1):
        self._pressure = None
        self._temperature = None

        self.pressure = pressure
        self.temperature = temperature
        self.region=-1
        self.set_region()
        # all the regions return v, rho, u, s, h, cp, cv, w
        self.v, self.rho, self.u, self.s, self.h, self.cp, self.cv, self.w = self.get_properties()

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self,value):
        if (value<273.15 or value>2273.15):
            raise ValueError ("Temperature outside allowed range")
        elif (value > 1073.15 and self.pressure>50):
            raise ValueError ("Temperature outside allowed range for specified pressure: <=1073.15 K/ 50-100 Mpa")
        else:
            self._temperature = value

    @property
    def pressure(self):
        return self._pressure

    @pressure.setter
    def pressure(self,value):
        if (value<0 or value>100):
            raise ValueError ("Pressure outside allowed range")
        self._pressure = value

    def set_region(self):
        if (self.temperature<=623.15):
            p_aux = dispatcher['Region 4'](self.temperature)[1]
            if (self.pressure>=p_aux):
                self.region='Region 1'
            else:
                self.region='Region 2'
        elif (self.temperature<=1073.15):
            p_aux = dispatcher['region23_boundary'](self.temperature)[1]
            if (self.pressure>=p_aux):
                self.region='Region 3'
            else:
                self.region='Region 2'
        else:
            self.region='Region 5'

    def get_properties(self):
        return dispatcher[self.region](self.temperature, self.pressure)

    def __str__(self):
        return "Pressure: {}MPa / Temperature:{}K".format(self.pressure, self.temperature)

def parseInput():
    pass

def normalizeTP(raw_temp=0, raw_pressure=1, units='SI'):
    pass

def pretty_print():
    pass

def main():
    """ Run the whole program """
    test = Droplet(1074,60)
    print (test)

if __name__ == '__main__':
    main()
