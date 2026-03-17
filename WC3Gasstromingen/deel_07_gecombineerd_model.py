import sys
sys.path.append("\\Users\morit\OneDrive\Documents\KT\year 3\IC&CS\Werkcolleges Computersimulatie")

import matplotlib.pyplot as plt
from physiomodeler import Model # type: ignore
from modellen.luchtstromingen import luchtstromingen_model
from modellen.gasstromingen import gasstromingen_model

gecombineerd_model = Model(
    dynamics=[luchtstromingen_model, gasstromingen_model],
    state_components=["volume_luchtwegen", "volume_alveoli"]
)

result = gecombineerd_model.run_simulation(
    time=300,
    inputs={
        "flux_O2_alveoli_PC": 0.25 / 60,
        "flux_CO2_alveoli_PC": -0.2 / 60,
    },
    initial_state={
        "volume_alveoli": 2.8,
        "volume_luchtwegen": 0.8,
    },
)

## plot
# tijd
time = result.index

# ======================
# O2 fracties ophalen
# ======================
fractie_O2_LWO = 0.1965
fractie_O2_LW  = result["fractie_O2_luchtwegen"]
fractie_O2_alv = result["fractie_O2_alveoli"]

plt.figure(figsize=(8,5))

plt.plot(time, [fractie_O2_LWO]*len(time), label="O2 luchtwegopening")
plt.plot(time, fractie_O2_LW, label="O2 luchtwegen")
plt.plot(time, fractie_O2_alv, label="O2 alveoli")

plt.xlabel("Tijd (s)")
plt.ylabel("Fractie O2")
plt.title("Zuurstoffracties in respiratoir systeem")
plt.legend()
plt.grid(True)
plt.show()


# ======================
# CO2 fracties ophalen
# ======================
fractie_CO2_LWO = 0.0003
fractie_CO2_LW  = result["fractie_CO2_luchtwegen"]
fractie_CO2_alv = result["fractie_CO2_alveoli"]

plt.figure(figsize=(8,5))

plt.plot(time, [fractie_CO2_LWO]*len(time), label="CO2 luchtwegopening")
plt.plot(time, fractie_CO2_LW, label="CO2 luchtwegen")
plt.plot(time, fractie_CO2_alv, label="CO2 alveoli")

plt.xlabel("Tijd (s)")
plt.ylabel("Fractie CO2")
plt.title("Koolstofdioxidefracties in respiratoir systeem")
plt.legend()
plt.grid(True)
plt.show()