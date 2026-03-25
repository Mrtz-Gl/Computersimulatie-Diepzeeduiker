from pathlib import Path
import sys

## Zorg ervoor dat de root van het project in sys.path staat, 
## zodat we modellen kunnen importeren op zowel macbook als windows machines
try:
    project_root = Path(__file__).resolve().parent.parent
except NameError:
    project_root = Path.cwd().parent

if str(project_root) not in sys.path:
    sys.path.append(str(project_root))

import matplotlib.pyplot as plt
from physiomodeler import Model # type: ignore

from modellen.flux import flux_alveoli_pc_model
from modellen.gasstromingen import gasstromingen_model
from modellen.luchtstromingen import (
    dynamic_elastance_pressure_model,
    luchtstromingen_model,
)
from modellen.perfusie import perfusie_model


parameters = {
    "rho_zeewater": 1025.0,  # [kg/m^3]
    "a": 6e-3,
    "P_atmosfeer": 101.0,  # [kPa]
    "oplosbaarheid_N2_bloed": 6.4e-6 * 25.4,  # [LN2/(Lbloed*kPa)]
    "valversnelling": 9.81,  # [m/s^2]
}


def dynamics(state, inputs, parameters):
    """Bereken de verandering in diepte van de duiker."""

    diepte = state["diepte"]
    a = parameters["a"]
    p_n2_sc = inputs["inhoud_N2_SC"] / parameters["oplosbaarheid_N2_bloed"]
    p_omgeving = (
        parameters["P_atmosfeer"]
        + (
            parameters["rho_zeewater"]
            * parameters["valversnelling"]
            * diepte
        )
        / 1000.0
    )

    if diepte <= 0:
        ddiepte = 0.0
    elif p_omgeving - p_n2_sc <= 0:
        ddiepte = 0.0
    else:
        ddiepte = -a * (p_omgeving - p_n2_sc)

    return {
        "ddiepte": ddiepte,
        "druk_omgeving": p_omgeving,
    }


duiker2 = Model(
    dynamics=[
        dynamics,
        dynamic_elastance_pressure_model,
        flux_alveoli_pc_model,
        luchtstromingen_model,
        gasstromingen_model,
        perfusie_model,
    ],
    state_components=["diepte"],
    parameters=parameters,
    initial_state={"diepte": 30},
)


initial_states = {
    "fractie_N2_luchtwegen": 0.71,  # 64% N2 in nitrox
    "fractie_O2_luchtwegen": 0.25,  # 36% O2 in nitrox
    "fractie_CO2_luchtwegen": 1 - 0.71 - 0.25,  # 4% als trace
    "fractie_O2_alveoli": 0.25,
    "fractie_CO2_alveoli": 1.0 - 0.25 - 0.71,  # 4% als trace
    "fractie_N2_alveoli": 0.71,
    "inhoud_O2_PC": 0.1987,
    "inhoud_CO2_PC": 0.5143,
    "inhoud_N2_PC": 0.046,
    "inhoud_O2_SA": 0.1987,
    "inhoud_CO2_SA": 0.5143,
    "inhoud_N2_SA": 0.046,
    "inhoud_O2_SC": 0.1531,
    "inhoud_CO2_SC": 0.5543,
    "inhoud_N2_SC": 0.046,
    "inhoud_O2_SV": 0.1531,
    "inhoud_CO2_SV": 0.5543,
    "inhoud_N2_SV": 0.046,
    "inhoud_N2_vet": 1e-5,
}


result = duiker2.run_simulation(
    time=60 * 4,  # minuten
    inputs={
        "flux_O2_SC_weefsels": 0.35 / 60,
        "breath_frequency": 12,
    },
    initial_state={**initial_states, "diepte": 30.0},
    relative_tolerance=1e-3,
)


tijd = result.index

fig, axes = plt.subplots(1, 2, figsize=(15, 4))

axes[0].plot(tijd, result["ddiepte"])
axes[0].set_title("Verandering van diepte")
axes[0].set_xlabel("Tijd (s)")
axes[0].set_ylabel("Stijgsnelheid (m/s)")
axes[0].grid(True)

axes[1].plot(
    tijd,
    result["partiele_druk_N2_SC"],
    label="Partiële stikstofdruk bloed",
)
axes[1].plot(tijd, result["druk_omgeving"], label="Omgevingsdruk")
axes[1].set_title("Drukken")
axes[1].set_xlabel("Tijd (s)")
axes[1].set_ylabel("Druk (kPa)")
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()


fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(tijd, result["inhoud_N2_vet"])
axes[0].set_title("Stikstofinhoud in vetweefsel")
axes[0].set_xlabel("Tijd (s)")
axes[0].set_ylabel("inhoud_N2_vet")
axes[0].grid(True)

axes[1].plot(
    tijd,
    result["partiele_druk_N2_vet"],
    label="Partiële stikstofdruk vet",
)
axes[1].plot(tijd, result["druk_omgeving"], label="Omgevingsdruk")
axes[1].set_title("Partiële stikstofdruk in vetweefsel")
axes[1].set_xlabel("Tijd (s)")
axes[1].set_ylabel("Druk (kPa)")
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()