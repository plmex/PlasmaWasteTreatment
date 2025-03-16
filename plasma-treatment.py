pip install iapws


import numpy as np
import iapws.iapws97 as iapws97
import iapws._iapws as iapws
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from iapws import IAPWS97
from math import pow


p = 0.101325
x_intervalos = np.linspace(0.001, 1, 10)
Tsat = iapws97._TSat_P(p)
estado = IAPWS97(P=p, T=Tsat)
Tc = estado.Tc

# Definindo as propriedades termodinâmicas específicas
vf = iapws97._Region1(T=Tsat, P=p)["v"]
hf = iapws97._Region1(T=Tsat, P=p)["h"]
sf = iapws97._Region1(T=Tsat, P=p)["s"]
vg = iapws97._Region2(T=Tsat, P=p)["v"]
hg = iapws97._Region2(T=Tsat, P=p)["h"]
sg = iapws97._Region2(T=Tsat, P=p)["s"]
Cp1 = iapws97._Region1(T=Tsat, P=p)["cp"]
Cp2 = iapws97._Region2(T=Tsat, P=p)["cp"]
alfa1 = iapws97._Region1(T=Tsat, P=p)["alfav"]
alfa2 = iapws97._Region2(T=Tsat, P=p)["alfav"]
kt1 = iapws97._Region1(T=Tsat, P=p)["kt"] / 1e6
kt2 = iapws97._Region2(T=Tsat, P=p)["kt"] / 1e6

# Valores referentes à mistura
vfg = vg - vf
hfg = hg - hf
sfg = sg - sf

# Derivadas e suas respectivas relações termodinâmicas
dvfdT =  vf * alfa1
dpdT = (sfg/vfg) * 1000

#derivada segunda da pressão em função da temperatura usando a relação de Clausius-Clapeyron
termo1 = (Cp2 - Cp1) / (Tsat * vfg)
termo2 = ((2 * (vg * alfa2 - vf * alfa1)) / vfg) * dpdT
termo3 = ((vg * kt2 - vf * kt1) / vfg) * pow(dpdT,2)
d2pdT2 = termo1 - termo2 + termo3

#pow((1 / (vf * (1 - x) + x * vg)),2) * pow(dpdT, 2)
C = []
for x in x_intervalos:
    numerador = pow((vf * (1 - x) + x * vg),2) * pow(dpdT, 2)
    fator1 = Cp1 / Tsat
    fator2 = x * vfg * d2pdT2
    fator3 = dpdT * dvfdT
    denominador = fator1 + fator2 - fator3
    C.append(pow(numerador / denominador, 0.5))
    print(f"Para x = {x:.1f}, C = {C[-1]:.2f} m/s")


print(f"dpdT: {dpdT :.2f} Pa/K")
print(f"d²p/dT²: {d2pdT2 :.2f} Pa/K²")
print(f"Cp 1: {Cp1:.6f} kJ/kgK")
print(f"Cp 2: {Cp2:.6f} kJ/kgK")
print(f"alfa 1: {alfa1:.6f} 1/K")
print(f"alfa 2: {alfa2:.6f} 1/K")
print(f"kt 1: {kt1:.12f} 1/Pa")
print(f"kt 2: {kt2:.6f} 1/Pa")

# valores para o gráfico da pressão
T_intervalos = np.linspace(273.15, Tc, 100)
pValores = []

for T in T_intervalos:
    p = iapws97._PSat_T(T)
    pValores.append(p)
"""
# Impressão dos valores de pressão
for i, p in enumerate(pValores):
    print(f"Temperatura: {T_intervalo[i]:.2f} K, Pressão: {p:.5f} MPa")
"""

# valores para o gráfico de volume

T_intervalos = np.linspace(273.15, 800, 100)
vValores = []

for T in T_intervalos:
    estado = IAPWS97(P=0.101325, T=T)
    vValores.append(estado.v)
"""
# Impressão dos valores de volume específico
for i, v in enumerate(vValores):
    print(f"Temperatura: {T_intervalo[i]:.2f} K, Volume específico: {v:.7f} m³/kg")
"""


fig, axs = plt.subplots(1, 3, figsize=(20, 5))

# Diagrama P-v
axs[0].plot(vValores, pValores, color='orange', linewidth=1.5)
axs[0].scatter(vValores[0], pValores[0], color='purple', label='Temperatura de saturação')  # Ponto Tsat
axs[0].scatter(vValores[-1], pValores[-1], color='red', label='Temperatura crítica')  # Ponto Tc
axs[0].set_xlabel('Volume Específico (m³/kg)')
axs[0].set_ylabel('Pressão (MPa)')
axs[0].set_title('Diagrama P-v')
axs[0].grid(True)
axs[0].legend()
axs[0].legend(loc='upper left')


# Gráfico p(T)
axs[1].plot(T_intervalos, pValores, color='green', linewidth=1.5)
axs[1].scatter(Tsat, pValores[0], color='purple', label='Temperatura de saturação')  # Ponto Tsat
axs[1].scatter(Tc, pValores[-1], color='red', label='Temperatura crítica')  # Ponto Tc
axs[1].set_xlabel('Temperatura (K)')
axs[1].set_ylabel('Pressão (MPa)')
axs[1].set_title('Gráfico p(T)')
axs[1].grid(True)
axs[1].axvspan(Tsat, Tc, color='yellow', alpha=0.3, label='Líquido-Vapor')  # Região de líquido-vapor
axs[1].axvspan(Tc, T_intervalos[-1], color='blue', alpha=0.3, label='Vapor')  # Região de Vapor
axs[1].axvspan(T_intervalos[0], Tsat, color='pink', alpha=0.3, label='Líquido')  # Região de líquido
axs[1].legend(loc='upper left')


# Gráfico v(T)
axs[2].plot(T_intervalos, vValores, color='deeppink', linewidth=1.5)
axs[2].scatter(Tsat, vValores[0], color='purple', label='Temperatura de saturação')  # Ponto Tsat
axs[2].scatter(Tc, vValores[-1], color='red', label='Temperatura crítica')  # Ponto Tc
axs[2].set_xlabel('Temperatura (K)')
axs[2].set_ylabel('Volume Específico (m³/kg)')
axs[2].set_title('Gráfico v(T)')
axs[2].grid(True)
axs[2].legend()
axs[2].fill_between(T_intervalos, vValores[0], vValores[-1], where=(T_intervalos >= Tsat) & (T_intervalos <= Tc), color='yellow', alpha=0.3, label='Líquido-Vapor')
axs[2].axvspan(Tc, T_intervalos[-1], color='blue', alpha=0.3, label='Vapor')  # Região de vapor
axs[2].axvspan(T_intervalos[0], Tsat, color='pink', alpha=0.3, label='Líquido')  # Região de líquido
axs[2].legend(loc='upper left')

plt.tight_layout()
plt.show()



plt.plot(C, x_intervalos * 100, marker='')
plt.ylabel('Título de Vapor (%)')
plt.xlabel('Velocidade do som (m/s)')
plt.title('Título de Vapor x Velocidade do Som')
plt.grid(True)
plt.show()