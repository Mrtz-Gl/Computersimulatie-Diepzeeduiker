import sys
sys.path.append("\\Users\morit\OneDrive\Documents\KT\year 3\IC&CS\Werkcolleges Computersimulatie")

import matplotlib.pyplot as plt
from physiomodeler import Model 
from modellen.luchtstromingen import luchtstromingen_model, dynamic_elastance_pressure_model  
from modellen.gasstromingen import gasstromingen_model                                        
from modellen.flux import flux_alveoli_PC_model                                             
from modellen.perfusie import perfusie_model


parameters = {
    "rho_zeewater": 1025.0,  # [kg/m^3]
    "a": 6e-3,
    "P_atmosfeer": 101.0,  # [kPa]
    "oplosbaarheid_N2_bloed": 6.4e-6 * 25.4,    # [LN2/(Lbloed*kPa)]
    "valversnelling": 9.81,  # [m/s^2]
}


def dynamics(state, inputs, parameters):
    "Berekent de verandering in diepte van de duiker"

    # Parameters
    diepte       = state["diepte"]
    a            = parameters["a"]
    P_N2_SC      = inputs["inhoud_N2_SC"] / parameters["oplosbaarheid_N2_bloed"]  
    P_omgeving   = parameters["P_atmosfeer"] + (parameters["rho_zeewater"] * parameters["valversnelling"] * diepte) / 1000.0  # vergelijking 2

    # als P_omgeving kleiner is dan P_N2_SC, 
    # dan stijgt de duiker niet verder en blijft diepte gelijk
    if diepte <= 0:
        ddiepte = 0
    elif P_omgeving - P_N2_SC <= 0:                                             # vergelijking 4
        ddiepte = 0.0
    else:
        ddiepte = -a * (P_omgeving - P_N2_SC)                                   # vergelijking 5

    return{
        "ddiepte": ddiepte, "druk_omgeving": P_omgeving
    }


duiker2 = Model(
    dynamics=[
        dynamics,
        dynamic_elastance_pressure_model, 
        flux_alveoli_PC_model,
        luchtstromingen_model,     
        gasstromingen_model,
        perfusie_model,
    ],
    state_components=["diepte"],
    parameters=parameters,
    initial_state={"diepte":30}
)


initial_states = {
       "fractie_N2_luchtwegen": 0.71,     # 64% N2 in nitrox
       "fractie_O2_luchtwegen": 0.25,    # 36% O2 in nitrox
       "fractie_CO2_luchtwegen": 1 - 0.71 - 0.25,     # 4% als trace in nitrox
       "fractie_O2_alveoli": 0.25,
       "fractie_CO2_alveoli": 1.0 - 0.25 - 0.71,     # 4% als trace in nitrox
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


# =========================
# Simulatie draaien vanaf equilibrium state
# =========================
result = duiker2.run_simulation( 
    time=60*4, # minuten
    inputs={
        "flux_O2_SC_weefsels": 0.35 / 60,
        "breath_frequency": 12,
        },
    initial_state={**initial_states, "diepte": 30.0},
    relative_tolerance=1e-3,
)


# tijdas uit de index van result
t = result.index

fig, axes = plt.subplots(1, 2, figsize=(15, 4))

# 1. ddiepte
axes[0].plot(t, result["ddiepte"])
axes[0].set_title("Verandering van diepte")
axes[0].set_xlabel("Tijd (s)")
axes[0].set_ylabel("stijgsnelheid (m/s)")
axes[0].grid(True)

# 2. partiële druk N2 en omgevingsdruk
axes[1].plot(t, result["partiele_druk_N2_SC"], label="partiele stikstofdruk bloed")
axes[1].plot(t, result["druk_omgeving"], label="omgevingsdruk")
axes[1].set_title("Drukken")
axes[1].set_xlabel("Tijd (s)")
axes[1].set_ylabel("Druk (kPa)")
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()

# =========================
# 4. Vetcompartiment
# =========================
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# inhoud N2 in vet
axes[0].plot(t, result["inhoud_N2_vet"])
axes[0].set_title("Stikstofinhoud in vetweefsel")
axes[0].set_xlabel("Tijd (s)")
axes[0].set_ylabel("inhoud_N2_vet")
axes[0].grid(True)

# partiële druk N2 in vet
axes[1].plot(t, result["partiele_druk_N2_vet"], label="partiele stikstofdruk vet")
axes[1].plot(t, result["druk_omgeving"], label="omgevingsdruk")
axes[1].set_title("Partiële stikstofdruk in vetweefsel")
axes[1].set_xlabel("Tijd (s)")
axes[1].set_ylabel("Druk (kPa)")
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()