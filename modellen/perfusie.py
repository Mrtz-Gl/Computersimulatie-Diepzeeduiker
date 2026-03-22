import matplotlib.pyplot as plt
from physiomodeler import Model  # type: ignore

parameters = {
    "volume_PC": 0.1,   # [L]
    "volume_SA": 1.1,   # [L]
    "volume_SC": 0.3,   # [L]
    "volume_SV": 3.5    # [L]
}

inputs = {
    "debiet_CO": 5 / 60,
    "flux_O2_SC_weefsels": 0.25 / 60,
    "flux_CO2_SC_weefsels": -0.2 / 60,
    "flux_N2_SC_weefsels": 0.0,
    "flux_O2_alveoli_PC": 0.25 / 60,
    "flux_CO2_alveoli_PC": -0.2 / 60,
    "flux_N2_alveoli_PC": 0.0,
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
    dynamics=concentratie_modellen,
    events=[event_inhoud_O2_negatief],
    state_components=["inhoud_O2_PC", "inhoud_CO2_PC", "inhoud_N2_PC"],
    parameters=parameters,
    inputs=inputs
)


if __name__ == "__main__":

    result = perfusie_model.run_simulation(
        time=60,
        initial_state={
            "inhoud_O2_PC": 0.005,
            "inhoud_O2_SA": 0.005,
            "inhoud_O2_SC": 0.0006,
            "inhoud_O2_SV": 0.0006,
            "inhoud_CO2_PC": 0.02,
            "inhoud_CO2_SA": 0.03,
            "inhoud_CO2_SC": 0.07,
            "inhoud_CO2_SV": 0.08,
            "inhoud_N2_PC": 0.0096,
            "inhoud_N2_SA": 0.0096,
            "inhoud_N2_SC": 0.0096,
            "inhoud_N2_SV": 0.0096,
        }
    )

    fig, axes = plt.subplots(1, 3, sharey=True)
    result[
        [
            "inhoud_O2_PC",
            "inhoud_O2_SA",
            "inhoud_O2_SC",
            "inhoud_O2_SV",
        ]
    ].plot(ax=axes[0])
    result[
        [
            "inhoud_CO2_PC",
            "inhoud_CO2_SA",
            "inhoud_CO2_SC",
            "inhoud_CO2_SV",
        ]
    ].plot(ax=axes[1])
    result[
        [
            "inhoud_N2_PC",
            "inhoud_N2_SA",
            "inhoud_N2_SC",
            "inhoud_N2_SV",
        ]
    ].plot(ax=axes[2])
