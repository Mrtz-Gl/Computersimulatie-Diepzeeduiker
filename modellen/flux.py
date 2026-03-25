from physiomodeler import Model  # type: ignore
import matplotlib.pyplot as plt
import numpy as np


parameters = {
    "omrekenfactor_cmH2O_kPA": 0.0981,          # [kPa/cmH2O]
    "bindingscapaciteit_Hb": 1.35e-3,           # [LO2/gHb]
    "concentratie_Hb": 150,                     # [LO2/Lbloed]
    "diffusiecapaciteit_O2": 0.0042,            # [LO2/(s*kPa)]
    "diffusiecapaciteit_CO2": 0.025,            # [LCO2/(s*kPa)]
    "diffusiecapaciteit_N2": 0.0021,              ## van stikstof
    "oplosbaarheid_N2_bloed": 6.4e-6 * 25.4,    # [LN2/(Lbloed*kPa)]
    "oplosbaarheid_N2_vet": 5 * 6.4e-6 * 25.4,  # [LN2/(Lvet*kPa)]
    "diffusiecapaciteit_N2_vet": 4.318e-9,          # [LN2/(s*kPa)]
    "valversnelling": 9.81,                      # [m/s^2]
    "P_atmosfeer": 101.0,  # [kPa]
    "rho_zeewater": 1025.0,  # [kg/m^3]
}


########################
### Flux berekenen
########################

def partiele_drukken_lucht(inputs):
    """Berekent de partiele druk in de alveoli en luchtwegen in kPA"""

    # definities
    druk_omgeving = parameters["P_atmosfeer"] + (parameters["rho_zeewater"] * parameters["valversnelling"] * inputs["diepte"]) / 1000.0
    omrekenfactor = parameters["omrekenfactor_cmH2O_kPA"]

    # bereken partiele drukken per compartiment
    partiele_drukken = {}

    for compartiment in ["luchtwegen", "alveoli"]:
        relatieve_druk = inputs[f"druk_{compartiment}"]
        fractie_O2 = inputs[f"fractie_O2_{compartiment}"]
        fractie_CO2 = inputs[f"fractie_CO2_{compartiment}"]
        fractie_N2 = inputs[f"fractie_N2_{compartiment}"]
        P_abs_kPa = relatieve_druk * omrekenfactor + druk_omgeving

        partiele_drukken[f"partiele_druk_O2_{compartiment}"] = P_abs_kPa * fractie_O2
        partiele_drukken[f"partiele_druk_CO2_{compartiment}"] = P_abs_kPa * fractie_CO2
        partiele_drukken[f"partiele_druk_N2_{compartiment}"] = P_abs_kPa * fractie_N2
    return partiele_drukken


def partiele_druk_O2_bij_saturatie(saturatie):
    """Berekent de partiele druk bij een saturatie volgens het Severinghaus model"""

    h2 = 2.81
    y_N = 55.5 * saturatie / (1 - saturatie)

    a = (y_N + np.sqrt(y_N**2 + h2)) / 2
    b = (y_N - np.sqrt(y_N**2 + h2)) / 2

    return np.sign(a) * np.abs(a) ** (1 / 3) + np.sign(
        b
    ) * np.abs(b) ** (1 / 3)


def partiele_druk_CO2_bij_inhoud(inhoud):
    """berekent de partiele druk van CO2 in bloed gegeven de totale inhoud CO2 in bloed"""
    return 0.837 * np.exp(4.16 * inhoud) - 0.895


def partiele_drukken_bloed(inputs, parameters):
    """Berekent de partiele drukken in de pulmonale capillairen"""
    # O2
    c_O2_Hb_a = inputs["inhoud_O2_PC"]
    c_O2_Hb_max = parameters["concentratie_Hb"] * parameters["bindingscapaciteit_Hb"]
    S_O2 = c_O2_Hb_a / c_O2_Hb_max

    # CO2
    inhoud_CO2 = inputs["inhoud_CO2_PC"]

    # N2
    inhoud_N2_PC = inputs["inhoud_N2_PC"]
    inhoud_N2_SC = inputs["inhoud_N2_SC"]
    oplosbaarheid_N2 = parameters["oplosbaarheid_N2_bloed"]
    return {
        "saturatie_O2_PC": S_O2,
        "partiele_druk_O2_PC": partiele_druk_O2_bij_saturatie(S_O2),
        "partiele_druk_CO2_PC": partiele_druk_CO2_bij_inhoud(inhoud_CO2),
        "partiele_druk_N2_SC": inhoud_N2_SC / oplosbaarheid_N2,                             # vergelijking 3
        "partiele_druk_N2_PC": inhoud_N2_PC / oplosbaarheid_N2,                             # vergelijking 3             
    }

def partiele_druk_vet(inputs, parameters):
    inhoud_N2_vet = inputs["inhoud_N2_vet"]
    oplosbaarheid_N2_vet = parameters["oplosbaarheid_N2_vet"]
    return {"partiele_druk_N2_vet": inhoud_N2_vet / oplosbaarheid_N2_vet}


def flux_alveoli_PC(inputs, parameters):
    """Berekent de flux van een gas x, tussen de alveoli en pulmonale capillairen"""

    gas = parameters["gas"]
    diffusiecapaciteit = parameters[f"diffusiecapaciteit_{gas}"]
    partiele_druk_alveoli = inputs[f"partiele_druk_{gas}_alveoli"]
    partiele_druk_PC = inputs[f"partiele_druk_{gas}_PC"]

    flux = diffusiecapaciteit * (partiele_druk_alveoli - partiele_druk_PC)
    return {f"flux_{gas}_alveoli_PC": flux}


def flux_N2_SC_vet(inputs, parameters):
     "Berekent N2 flux van systemische capillairen naar vetweefsel"
     
     partiele_druk_N2_bloed = inputs["partiele_druk_N2_SC"]
     partiele_druk_N2_vet   = partiele_druk_vet(inputs, parameters)["partiele_druk_N2_vet"]
     diffusiecapaciteit     = parameters["diffusiecapaciteit_N2_vet"]

     flux = diffusiecapaciteit * (partiele_druk_N2_bloed - partiele_druk_N2_vet)
     return {"flux_N2_SC_vet": flux}


flux_O2_alveoli_PC_model = Model(
    dynamics=flux_alveoli_PC,
    parameters={"gas": "O2"},
)

flux_CO2_alveoli_PC_model = Model(
    dynamics=flux_alveoli_PC,
    parameters={"gas": "CO2"},
)

flux_N2_alveoli_PC_model = Model(
    dynamics=flux_alveoli_PC,
    parameters={"gas": "N2"},
)

flux_N2_SC_vet_model = Model(
     dynamics=[
         flux_N2_SC_vet, 
         partiele_druk_vet
         ]
 )


flux_alveoli_PC_model = Model(
    dynamics=[
        partiele_drukken_lucht,
        partiele_drukken_bloed,
        flux_O2_alveoli_PC_model,
        flux_CO2_alveoli_PC_model,
        flux_N2_alveoli_PC_model,
        flux_N2_SC_vet_model,
    ],
    parameters=parameters,
)



