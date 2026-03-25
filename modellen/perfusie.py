import matplotlib.pyplot as plt
from physiomodeler import Model  # type: ignore

parameters = {
    "volume_PC": 0.1,       # [L]
    "volume_SA": 1.1,       # [L]
    "volume_SC": 0.3,       # [L]
    "volume_SV": 3.5,       # [L]
    "volume_vet": 22.0      # [L]
}


inputs = {
    "debiet_CO": 5 / 60,
    "flux_O2_SC_weefsels": 0.25 / 60,
    "flux_CO2_SC_weefsels": -0.2 / 60,
    "flux_N2_SC_weefsels": 0.0,
    "flux_O2_alveoli_PC": 0.25 / 60,
    "flux_CO2_alveoli_PC": -0.2 / 60,
    "flux_N2_alveoli_PC": 0.2 / 60,
}


def dynamics(inputs, state, parameters):
    """Berekent gasstroming in bloed"""

    gas = parameters["gas"]
    compartiment_a = parameters["voorgaand_compartiment"]
    compartiment_b = parameters["compartiment"]

    inhoud_a = inputs[f"inhoud_{gas}_{compartiment_a}"]
    inhoud_b = state[f"inhoud_{gas}_{compartiment_b}"]
    V_b = parameters[f"volume_{compartiment_b}"]
    Q_CO = inputs["debiet_CO"]

    dinhoud_a_b = (inhoud_a - inhoud_b) * Q_CO / V_b

    if compartiment_b == "PC":
        dinhoud_a_b += inputs[f"flux_{gas}_alveoli_PC"] / V_b
    elif compartiment_b == "SC":
        dinhoud_a_b -= inputs[f"flux_{gas}_SC_weefsels"] / V_b

    return {
        f"dinhoud_{gas}_{compartiment_b}": dinhoud_a_b
    }


def N2_inhoud_vet(inputs, parameters):
    """Berekent de inhoud van N2 in het vetweefsel"""
    V_vet = parameters["volume_vet"]

    d_inhoud_N2_vet = inputs["flux_N2_SC_vet"] / V_vet

    return {
        "dinhoud_N2_vet": d_inhoud_N2_vet
    }


def event_inhoud_O2_negatief(state):
    """retourneert de inhoud zuurstof in de systemische capillairen"""
    return state["inhoud_O2_SC"]


event_inhoud_O2_negatief.terminal = True


compartimenten = ["PC", "SA", "SC", "SV"]
gassen = ["O2", "CO2", "N2"]

concentratie_modellen = []

for i in range(len(compartimenten)):
    for gas in gassen:
        compartiment = compartimenten[i]
        voorgaande_compartiment = compartimenten[i - 1]

        concentratie_model = Model(
            dynamics=dynamics,
            state_components=[f"inhoud_{gas}_{compartiment}"],
            parameters={
                "gas": gas,
                "compartiment": compartiment,
                "voorgaand_compartiment": voorgaande_compartiment,
            },
        )

        concentratie_modellen.append(concentratie_model)

perfusie_model = Model(
dynamics=[*concentratie_modellen, N2_inhoud_vet], #  dynamics=[*concentratie_modellen, N2_inhoud_vet],
    events=[event_inhoud_O2_negatief],
    state_components=["inhoud_O2_PC", "inhoud_CO2_PC", "inhoud_N2_PC", "inhoud_N2_vet"],
    parameters=parameters,
    inputs=inputs
)


