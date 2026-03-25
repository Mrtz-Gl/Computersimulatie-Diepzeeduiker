from physiomodeler import Model  # type: ignore
import matplotlib.pyplot as plt
import numpy as np


parameters = {
    "omrekenfactor_cmH2O_kPA": 0.0981,  # [kPa/cmH2O]
    "bindingscapaciteit_Hb": 1.35e-3,  # [LO2/gHb]
    "concentratie_Hb": 150,  # [LO2/Lbloed]
    "diffusiecapaciteit_O2": 0.0042,  # [LO2/(s*kPa)]
    "diffusiecapaciteit_CO2": 0.025,  # [LCO2/(s*kPa)]
    "diffusiecapaciteit_N2": 0.0021,  # [LN2/(s*kPa)]
    "oplosbaarheid_N2_bloed": 6.4e-6 * 25.4,  # [LN2/(Lbloed*kPa)]
    "oplosbaarheid_N2_vet": 5 * 6.4e-6 * 25.4,  # [LN2/(Lvet*kPa)]
    "diffusiecapaciteit_N2_vet": 4.318e-9,  # [LN2/(s*kPa)]
    "valversnelling": 9.81,  # [m/s^2]
    "P_atmosfeer": 101.0,  # [kPa]
    "rho_zeewater": 1025.0,  # [kg/m^3]
}


def partiele_drukken_lucht(inputs):
    """Bereken de partiële drukken in de alveoli en luchtwegen in kPa."""

    druk_omgeving = (
        parameters["P_atmosfeer"]
        + (
            parameters["rho_zeewater"]
            * parameters["valversnelling"]
            * inputs["diepte"]
        )
        / 1000.0
    )
    omrekenfactor = parameters["omrekenfactor_cmH2O_kPA"]

    partiele_drukken = {}

    for compartiment in ["luchtwegen", "alveoli"]:
        relatieve_druk = inputs[f"druk_{compartiment}"]
        fractie_o2 = inputs[f"fractie_O2_{compartiment}"]
        fractie_co2 = inputs[f"fractie_CO2_{compartiment}"]
        fractie_n2 = inputs[f"fractie_N2_{compartiment}"]

        p_abs_kpa = relatieve_druk * omrekenfactor + druk_omgeving

        partiele_drukken[f"partiele_druk_O2_{compartiment}"] = (
            p_abs_kpa * fractie_o2
        )
        partiele_drukken[f"partiele_druk_CO2_{compartiment}"] = (
            p_abs_kpa * fractie_co2
        )
        partiele_drukken[f"partiele_druk_N2_{compartiment}"] = (
            p_abs_kpa * fractie_n2
        )

    return partiele_drukken


def partiele_druk_o2_bij_saturatie(saturatie):
    """Bereken de partiële O2-druk bij een saturatie volgens Severinghaus."""

    h2 = 2.81
    y_n = 55.5 * saturatie / (1 - saturatie)

    a = (y_n + np.sqrt(y_n**2 + h2)) / 2
    b = (y_n - np.sqrt(y_n**2 + h2)) / 2

    return np.sign(a) * np.abs(a) ** (1 / 3) + np.sign(b) * np.abs(b) ** (1 / 3)


def partiele_druk_co2_bij_inhoud(inhoud):
    """Bereken de partiële CO2-druk in bloed gegeven de totale CO2-inhoud."""

    return 0.837 * np.exp(4.16 * inhoud) - 0.895


def partiele_drukken_bloed(inputs, parameters):
    """Bereken de partiële drukken in de pulmonale capillairen."""

    inhoud_o2_pc = inputs["inhoud_O2_PC"]
    c_o2_hb_max = parameters["concentratie_Hb"] * parameters["bindingscapaciteit_Hb"]
    saturatie_o2 = inhoud_o2_pc / c_o2_hb_max

    inhoud_co2_pc = inputs["inhoud_CO2_PC"]

    inhoud_n2_pc = inputs["inhoud_N2_PC"]
    inhoud_n2_sc = inputs["inhoud_N2_SC"]
    oplosbaarheid_n2 = parameters["oplosbaarheid_N2_bloed"]

    return {
        "saturatie_O2_PC": saturatie_o2,
        "partiele_druk_O2_PC": partiele_druk_o2_bij_saturatie(saturatie_o2),
        "partiele_druk_CO2_PC": partiele_druk_co2_bij_inhoud(inhoud_co2_pc),
        "partiele_druk_N2_SC": inhoud_n2_sc / oplosbaarheid_n2,
        "partiele_druk_N2_PC": inhoud_n2_pc / oplosbaarheid_n2,
    }


def partiele_druk_vet(inputs, parameters):
    """Bereken de partiële stikstofdruk in vetweefsel."""

    inhoud_n2_vet = inputs["inhoud_N2_vet"]
    oplosbaarheid_n2_vet = parameters["oplosbaarheid_N2_vet"]

    return {
        "partiele_druk_N2_vet": inhoud_n2_vet / oplosbaarheid_n2_vet,
    }


def flux_alveoli_pc(inputs, parameters):
    """Bereken de flux van een gas tussen alveoli en pulmonale capillairen."""

    gas = parameters["gas"]
    diffusiecapaciteit = parameters[f"diffusiecapaciteit_{gas}"]
    partiele_druk_alveoli = inputs[f"partiele_druk_{gas}_alveoli"]
    partiele_druk_pc = inputs[f"partiele_druk_{gas}_PC"]

    flux = diffusiecapaciteit * (partiele_druk_alveoli - partiele_druk_pc)

    return {
        f"flux_{gas}_alveoli_PC": flux,
    }


def flux_n2_sc_vet(inputs, parameters):
    """Bereken de N2-flux van systemische capillairen naar vetweefsel."""

    partiele_druk_n2_bloed = inputs["partiele_druk_N2_SC"]
    partiele_druk_n2_vet = partiele_druk_vet(
        inputs, parameters
    )["partiele_druk_N2_vet"]
    diffusiecapaciteit = parameters["diffusiecapaciteit_N2_vet"]

    flux = diffusiecapaciteit * (
        partiele_druk_n2_bloed - partiele_druk_n2_vet
    )

    return {
        "flux_N2_SC_vet": flux,
    }


flux_o2_alveoli_pc_model = Model(
    dynamics=flux_alveoli_pc,
    parameters={"gas": "O2"},
)

flux_co2_alveoli_pc_model = Model(
    dynamics=flux_alveoli_pc,
    parameters={"gas": "CO2"},
)

flux_n2_alveoli_pc_model = Model(
    dynamics=flux_alveoli_pc,
    parameters={"gas": "N2"},
)

flux_n2_sc_vet_model = Model(
    dynamics=[
        flux_n2_sc_vet,
        partiele_druk_vet,
    ]
)

flux_alveoli_pc_model = Model(
    dynamics=[
        partiele_drukken_lucht,
        partiele_drukken_bloed,
        flux_o2_alveoli_pc_model,
        flux_co2_alveoli_pc_model,
        flux_n2_alveoli_pc_model,
        flux_n2_sc_vet_model,
    ],
    parameters=parameters,
)