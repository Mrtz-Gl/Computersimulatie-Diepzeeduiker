from physiomodeler import Model # type: ignore
import matplotlib.pyplot as plt
import numpy as np


inputs = {
    "druk_atmosfeer": 1033  # [cmH2O]
}

parameters = {
    "omrekenfactor_cmH2O_kPA": 0.0981,   # [kPa/cmH2O]
    "bindingscapaciteit_Hb": 1.35e-3,    # [LO2/gHb]
    "concentratie_Hb": 150,              # [LO2/Lbloed]
    "diffusiecapaciteit_O2": 0.0042,     # [LO2/(s*kPa)]
    "diffusiecapaciteit_CO2": 0.025      # [LCO2/(s*kPa)]         
}


def partiele_drukken_lucht(inputs):
    """Berekent de partiele druk in de alveoli en luchtwegen in kPA"""
    
    # definities
    druk_atm      = inputs["druk_atmosfeer"]
    omrekenfactor = parameters["omrekenfactor_cmH2O_kPA"]
    fractie_O2  = inputs["fractie_O2_luchtwegen"]
    fractie_CO2 = inputs["fractie_CO2_luchtwegen"]

    # bereken partiele drukken per compartiment
    partiele_drukken = {}

    for compartiment in ["luchtwegen", "alveoli"]:
        relatieve_druk = inputs[f"druk_{compartiment}"]
        P_abs_kPa      = (relatieve_druk + druk_atm) * omrekenfactor

        # voeg too aan partiele_drukken dict
        partiele_drukken[f"partiele_druk_O2_{compartiment}"]  = P_abs_kPa * fractie_O2
        partiele_drukken[f"partiele_druk_CO2_{compartiment}"] = P_abs_kPa * fractie_CO2
    return partiele_drukken


def partiele_druk_O2_bij_saturatie(saturatie):
    """Berekent de partiele druk bij een saturatie volgens het Severinghaus model"""

    h2 = 2.81
    y_N = 55.5 * saturatie / (1 - saturatie)

    a = (y_N + np.sqrt(y_N**2 + h2)) / 2
    b = (y_N - np.sqrt(y_N**2 + h2)) / 2

    # Een negatief getal tot een breuk verheffen geeft een complex getal. In
    # plaats daarvan converteren we naar een positief getal en zetten we het
    # originele teken er weer bij terug na machtsverheffing.
    # np.abs(a) converteert a naar een positief getal.
    # np.sign(a) converteert na machtsverheffing terug naar het originele teken.
    return np.sign(a) * np.abs(a) ** (1 / 3) + np.sign(
        b
    ) * np.abs(b) ** (1 / 3)


def partiele_druk_CO2_bij_inhoud(inhoud):
    """berekent de partiele druk van CO2 in bloed gegeven de totale inhoud CO2 in bloed"""
    return 0.837 * np.exp(4.16 * inhoud) - 0.895


def partiele_drukken_bloed(inputs, parameters):
    """Berekent de zuurstof-saturatie in de pulmonale capillairen"""

    # huidige inhoud gebonden zuurstof
    c_O2_Hb_a   = inputs["inhoud_O2_PC"]
    # maximale inhoud gebonden zuurstof
    c_O2_Hb_max = parameters["concentratie_Hb"] * parameters["bindingscapaciteit_Hb"]
    # saturatie
    S_O2 = c_O2_Hb_a / c_O2_Hb_max

    # huidige inhoud CO2
    inhoud_CO2 = inputs["inhoud_CO2_PC"]
    return{
        "saturatie_O2_PC": S_O2,
        "partiele_druk_O2_PC": partiele_druk_O2_bij_saturatie(S_O2),
        "partiele_druk_CO2_PC": partiele_druk_CO2_bij_inhoud(inhoud_CO2)
    }


def flux_alveoli_PC(inputs, parameters):
    """Berekent de flux van een gas x tussen de alveoli en pulmonale capillairen"""

    gas = parameters["gas"]
    diffusiecapaciteit = parameters[f"diffusiecapaciteit_{gas}"]
    partiele_druk_alveoli = inputs[f"partiele_druk_{gas}_alveoli"]
    partiele_druk_PC      = inputs[f"partiele_druk_{gas}_PC"]

    flux = diffusiecapaciteit * (partiele_druk_alveoli - partiele_druk_PC)
    return {f"flux_{gas}_alveoli_PC": flux}


# ============
# creer modellen...
flux_O2_alveoli_PC_model = Model(
    dynamics=flux_alveoli_PC,
    parameters={"gas": "O2"},
)

flux_CO2_alveoli_PC_model = Model(
    dynamics=flux_alveoli_PC,
    parameters={"gas": "CO2"},
)

# ...combineer tot een model
flux_alveoli_PC_model = Model(
    dynamics=[
        partiele_drukken_lucht,
        partiele_drukken_bloed,
        flux_O2_alveoli_PC_model,
        flux_CO2_alveoli_PC_model,
    ],
    parameters=parameters,
    inputs=inputs,
)

result = flux_alveoli_PC_model.run_simulation(
    time=10,
    inputs={
        "druk_alveoli": 0,
        "druk_luchtwegen": 0,
        "fractie_O2_alveoli": 0.15,
        "fractie_CO2_alveoli": 0.06,
        "fractie_O2_luchtwegen": 0.19,
        "fractie_CO2_luchtwegen": 0.01,
        "inhoud_O2_PC": 0.1987,
        "inhoud_CO2_PC": 0.5143,
        "inhoud_O2_SA": 0.1987,
        "inhoud_CO2_SA": 0.5143,
        "inhoud_O2_SC": 0.1531,
        "inhoud_CO2_SC": 0.5543,
        "inhoud_O2_SV": 0.1531,
        "inhoud_CO2_SV": 0.5543,
    },
)
