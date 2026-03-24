from physiomodeler import Model  # type: ignore
import matplotlib.pyplot as plt
import numpy as np


def druk_atmosfeer_duiker(tijd):
    """Omgevingsdruk bij opstijgen van 30m met 10 m/min."""
    start_diepte = 30.0
    stijg_stappen_m_p_min = 10.0f

    diepte_m = max(0.0, start_diepte - stijg_stappen_m_p_min * (tijd / 60.0))
    atm_oppervlak_cmH2O = 1033.0
    return (1.0 + diepte_m / 10.0) * atm_oppervlak_cmH2O


inputs = {
    "druk_atmosfeer": 1033  # [cmH2O]
}

parameters = {
    "omrekenfactor_cmH2O_kPA": 0.0981,          # [kPa/cmH2O]
    "bindingscapaciteit_Hb": 1.35e-3,           # [LO2/gHb]
    "concentratie_Hb": 150,                     # [LO2/Lbloed]
    "diffusiecapaciteit_O2": 0.0042,            # [LO2/(s*kPa)]
    "diffusiecapaciteit_CO2": 0.025,            # [LCO2/(s*kPa)]
    "diffusiecapaciteit_N2": 0.021,              ## van stikstof
    "oplosbaarheid_N2_bloed": 6.4e-3 * 22.4,     # [LN2/(Lbloed*kPa)]
    "oplosbaarheid_N2_vet": 5 * 6.4e-3 * 22.4,    # [LN2/(Lvet*kPa)]
    "diffusiecapaciteit_N2_vet": 4.318e-9,          # [LN2/(s*kPa)]
}


########################
### Flux berekenen
########################


def partiele_drukken_lucht(inputs):
    """Berekent de partiele druk in de alveoli en luchtwegen in kPA"""

    # definities
    druk_atm = inputs["druk_atmosfeer"]
    omrekenfactor = parameters["omrekenfactor_cmH2O_kPA"]

    # bereken partiele drukken per compartiment
    partiele_drukken = {}

    for compartiment in ["luchtwegen", "alveoli"]:
        relatieve_druk = inputs[f"druk_{compartiment}"]
        fractie_O2 = inputs[f"fractie_O2_{compartiment}"]
        fractie_CO2 = inputs[f"fractie_CO2_{compartiment}"]
        fractie_N2 = inputs[f"fractie_N2_{compartiment}"]
        P_abs_kPa = (relatieve_druk + druk_atm) * omrekenfactor

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
    inhoud_N2 = inputs["inhoud_N2_PC"]
    oplosbaarheid_N2 = parameters["oplosbaarheid_N2_bloed"]
    return {
        "saturatie_O2_PC": S_O2,
        "partiele_druk_O2_PC": partiele_druk_O2_bij_saturatie(S_O2),
        "partiele_druk_CO2_PC": partiele_druk_CO2_bij_inhoud(inhoud_CO2),
        "partiele_druk_N2_PC": inhoud_N2 / oplosbaarheid_N2,
    }


def flux_alveoli_PC(inputs, parameters):
    """Berekent de flux van een gas x, tussen de alveoli en pulmonale capillairen"""

    gas = parameters["gas"]
    diffusiecapaciteit = parameters[f"diffusiecapaciteit_{gas}"]
    partiele_druk_alveoli = inputs[f"partiele_druk_{gas}_alveoli"]
    partiele_druk_PC = inputs[f"partiele_druk_{gas}_PC"]

    flux = diffusiecapaciteit * (partiele_druk_alveoli - partiele_druk_PC)
    return {f"flux_{gas}_alveoli_PC": flux}


def flux_N2_SC_fat(inputs, parameters):
    """Berekent N2 flux van systemische capillairen naar vetweefsel"""
    
    partiele_druk_N2_bloed = inputs["partiele_druk_N2_PC"]
    inhoud_N2_fat = inputs["inhoud_N2_fat"]
    oplosbaarheid_N2_vet = parameters["oplosbaarheid_N2_vet"]
    diffusiecapaciteit = parameters["diffusiecapaciteit_N2_vet"]
    
    partiele_druk_N2_fat = inhoud_N2_fat / oplosbaarheid_N2_vet
    flux = diffusiecapaciteit * (partiele_druk_N2_bloed - partiele_druk_N2_fat)
    return {"flux_N2_SC_fat": flux}


##########################################
### Risico op decompressieziekte berekenen
##########################################

def dcs_risico(inputs, parameters):
    """DCS-check: supersaturatie in vet."""
    partiele_druk_N2_fat = inputs["inhoud_N2_fat"] / parameters["oplosbaarheid_N2_vet"]
    partiele_druk_N2_bloed = inputs["partiele_druk_N2_PC"]
    supersaturatie = partiele_druk_N2_fat / partiele_druk_N2_bloed
    return {"dcs_risico": supersaturatie > 1.5}  # True als risico


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

flux_N2_SC_fat_model = Model(
    dynamics=flux_N2_SC_fat,
)


flux_alveoli_PC_model = Model(
    dynamics=[
        partiele_drukken_lucht,
        partiele_drukken_bloed,
        flux_O2_alveoli_PC_model,
        flux_CO2_alveoli_PC_model,
        flux_N2_alveoli_PC_model,
        flux_N2_SC_fat_model,
    ],
    parameters=parameters,
    inputs=inputs,
)


if __name__ == "__main__":

    result = flux_alveoli_PC_model.run_simulation(
        time=10,
        inputs={
            "druk_alveoli": 0,
            "druk_luchtwegen": 0,
            "fractie_O2_alveoli": 0.15,
            "fractie_CO2_alveoli": 1.0 - 0.15 - 0.7097,
            "fractie_N2_alveoli": 0.7097,
            "fractie_N2_luchtwegen": 0.7097,     # 71% N2 in nitrox
            "fractie_O2_luchtwegen": 0.2099,     # 29% O2 in nitrox
            "fractie_CO2_luchtwegen": 0.04,      # 4% als trace in nitrox
            "inhoud_O2_PC": 0.1987,
            "inhoud_CO2_PC": 0.5143,
            "inhoud_N2_PC": 0.0096,
            "inhoud_O2_SA": 0.1987,
            "inhoud_CO2_SA": 0.5143,
            "inhoud_N2_SA": 0.0096,
            "inhoud_O2_SC": 0.1531,
            "inhoud_CO2_SC": 0.5543,
            "inhoud_N2_SC": 0.0096,
            "inhoud_O2_SV": 0.1531,
            "inhoud_CO2_SV": 0.5543,
            "inhoud_N2_SV": 0.0096,
        },
    )
