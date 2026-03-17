import sys
sys.path.append("\\Users\morit\OneDrive\Documents\KT\year 3\IC&CS\Werkcolleges Computersimulatie")

import matplotlib.pyplot as plt
from physiomodeler import Model # type: ignore
from modellen.elasticiteiten import dynamic_elastance_model
from modellen.luchtstromingen import luchtstromingen_model
from modellen.gasstromingen import gasstromingen_model
from modellen.flux import flux_alveoli_PC_model
from modellen.perfusie import perfusie_model


systeem_model = Model(
    dynamics=[
        dynamic_elastance_model,
        luchtstromingen_model,
        flux_alveoli_PC_model,
        gasstromingen_model,
        perfusie_model
    ]
)


duur_ademteug = (60 / luchtstromingen_model.inputs["breath_frequency"])

equilibrium_state = systeem_model.find_equilibrium_state(
    period=duur_ademteug,
    estimated_equilibrium_state = {
        "volume_lungs": 3.0,
        "fractie_O2_luchtwegen": 0.21,
        "fractie_CO2_luchtwegen": 0.0004,
        "fractie_O2_alveoli": 0.15,
        "fractie_CO2_alveoli": 0.05,
        "inhoud_O2_PC": 0.1987,
        "inhoud_CO2_PC": 0.5143,
        "inhoud_O2_SA": 0.1987,
        "inhoud_CO2_SA": 0.5143,
        "inhoud_O2_SC": 0.1531,
        "inhoud_CO2_SC": 0.5543,
        "inhoud_O2_SV": 0.1531,
        "inhoud_CO2_SV": 0.5543
        },
    relative_tolerance_simulation=1e-3,
    relative_tolerance_equilibrium=1e-2,
)



result = systeem_model.result
t = result["time"]

fig, ax = plt.subplots(3, 1, sharex=True)

# drukken
ax[0].plot(t, result["druk_alveoli"], label="Alveoli")
ax[0].plot(t, result["druk_luchtwegen"], label="Luchtwegen")
ax[0].set_ylabel("Druk")
ax[0].legend()

# debiet
ax[1].plot(t, result["debiet_luchtwegen"])
ax[1].set_ylabel("Debiet")

# volume
ax[2].plot(t, result["volume_lungs"])
ax[2].set_ylabel("Volume")
ax[2].set_xlabel("Tijd (s)")

plt.show()


plt.figure()

plt.plot(t, result["pO2_luchtwegen"], label="luchtwegen")
plt.plot(t, result["pO2_alveoli"], label="alveoli")
plt.plot(t, result["pO2_PC"], label="pulmonale capillairen")
plt.plot(t, result["pO2_SA"], label="systemische arteriën")
plt.plot(t, result["pO2_SC"], label="systemische capillairen")
plt.plot(t, result["pO2_SV"], label="systemische venen")

plt.ylabel("pO2")
plt.xlabel("Tijd (s)")
plt.legend()
plt.show()


plt.figure()

plt.plot(t, result["pCO2_luchtwegen"], label="luchtwegen")
plt.plot(t, result["pCO2_alveoli"], label="alveoli")
plt.plot(t, result["pCO2_PC"], label="pulmonale capillairen")
plt.plot(t, result["pCO2_SA"], label="systemische arteriën")
plt.plot(t, result["pCO2_SC"], label="systemische capillairen")
plt.plot(t, result["pCO2_SV"], label="systemische venen")

plt.ylabel("pCO2")
plt.xlabel("Tijd (s)")
plt.legend()
plt.show()


plt.figure()

plt.plot(t, result["saturatie_O2_PC"], label="PC")
plt.plot(t, result["saturatie_O2_SA"], label="SA")
plt.plot(t, result["saturatie_O2_SC"], label="SC")
plt.plot(t, result["saturatie_O2_SV"], label="SV")

plt.ylabel("O2 saturatie")
plt.xlabel("Tijd (s)")
plt.legend()
plt.show()