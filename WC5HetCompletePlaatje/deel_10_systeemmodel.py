import sys
sys.path.append("\\Users\morit\OneDrive\Documents\KT\year 3\IC&CS\Werkcolleges Computersimulatie")

import matplotlib.pyplot as plt
from physiomodeler import Model # type: ignore
from modellen.luchtstromingen import luchtstromingen_model, dynamic_elastance_pressure_model  # deel 5
from modellen.gasstromingen import gasstromingen_model                                        # deel 6
from modellen.flux import flux_alveoli_PC_model                                               # deel 9
from modellen.perfusie import perfusie_model                                                  # deel 8


systeem_model = Model(
    dynamics=[
        dynamic_elastance_pressure_model,   # dit nog wel vertalen naar Nederlands
        flux_alveoli_PC_model,
        luchtstromingen_model,     
        gasstromingen_model,
        perfusie_model
    ]
)


duur_ademteug = (60 / luchtstromingen_model.inputs["breath_frequency"])

equilibrium_state = systeem_model.find_equilibrium_state(
    period=duur_ademteug,
    estimated_equilibrium_state = {
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


# Simulatie draaien vanaf equilibrium state
result = systeem_model.run_simulation(
    time=60*10, # 10 minuten
    inputs={"flux_O2_SC_weefsels": 0.35 / 60},
    initial_state=equilibrium_state,
    relative_tolerance=1e-3,
)

# =========================
# 1. Drukken, debiet, volume
# =========================
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 8))

result[
    [
        "druk_alveoli",
        "druk_luchtwegen",
        "druk_thorax",
    ]
].plot(ax=axes[0])
axes[0].set_ylabel("Druk")
axes[0].set_title("Drukken")

result[
    [
        "debiet_luchtwegopening_luchtwegen",
        "debiet_luchtwegen_alveoli",
    ]
].plot(ax=axes[1])
axes[1].set_ylabel("Debiet")
axes[1].set_title("Debieten")

result[
    [
        "volume_luchtwegen",
        "volume_alveoli",
    ]
].plot(ax=axes[2])
axes[2].set_ylabel("Volume")
axes[2].set_xlabel("Tijd (s)")
axes[2].set_title("Volumes")

plt.tight_layout()
plt.show()


# =========================
# 2. Partiële druk van O2 in alle compartimenten
# =========================
plt.figure(figsize=(10, 5))

result[
    [
        "partiele_druk_O2_luchtwegen",
        "partiele_druk_O2_alveoli",
        "partiele_druk_O2_PC",
    ]
].plot()

plt.ylabel("Partiële druk O2")
plt.xlabel("Tijd (s)")
plt.title("Partiële druk van zuurstof")
plt.legend()
plt.tight_layout()
plt.show()


# =========================
# 3. Partiële druk van CO2 in alle compartimenten
# =========================
plt.figure(figsize=(10, 5))

result[
    [
        "partiele_druk_CO2_luchtwegen",
        "partiele_druk_CO2_alveoli",
        "partiele_druk_CO2_PC",
    ]
].plot()

plt.ylabel("Partiële druk CO2")
plt.xlabel("Tijd (s)")
plt.title("Partiële druk van koolstofdioxide")
plt.legend()
plt.tight_layout()
plt.show()


# =========================
# 4. O2-saturatie in alle bloedgevulde compartimenten
# =========================
plt.figure(figsize=(10, 5))

result[
    [
        "saturatie_O2_PC",
    ]
].plot()

plt.ylabel("O2-saturatie")
plt.xlabel("Tijd (s)")
plt.title("Zuurstofsaturatie")
plt.legend()
plt.tight_layout()
plt.show()