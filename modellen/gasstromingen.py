from physiomodeler import Model # type: ignore
import matplotlib.pyplot as plt
from scipy.signal import square
import numpy as np


# define the state-space system
inputs = {
    "fractie_O2_luchtwegopening": 0.1965,
    "fractie_CO2_luchtwegopening": 0.0003
}

# run the state-space system
def dfractie(Q_ab, fractie_a, fractie_b, V_b):
    """Berekent de verandering van de gasfractie"""

    if Q_ab <= 0:
        return 0
    return (fractie_a - fractie_b) * (Q_ab / V_b)


def dynamics(state, inputs, parameters):
    """Berekent de gasstroming"""

    gas = parameters["gas"]

    # parameters...
    # ...flow
    Q_LWO_LW = inputs["debiet_luchtwegopening_luchtwegen"]
    Q_LW_alv = inputs["debiet_luchtwegen_alveoli"]
    # ...fractie
    fractie_LWO = inputs[f"fractie_{gas}_luchtwegopening"]
    fractie_alv = state[f"fractie_{gas}_alveoli"]
    fractie_LW  = state[f"fractie_{gas}_luchtwegen"]
    # ...volume
    V_LW  = inputs["volume_luchtwegen"]
    V_alv = inputs["volume_alveoli"]
    # ...flux
    flux_alv = inputs[f"flux_{gas}_alveoli_PC"]

    # bereken
    dfractie_LWO_LW = dfractie(Q_LWO_LW, fractie_LWO, fractie_LW,  V_LW)
    dfractie_alv_LW = dfractie(-Q_LW_alv, fractie_alv, fractie_LW,  V_LW)
    dfractie_LW_alv = dfractie(Q_LW_alv, fractie_LW,  fractie_alv, V_alv) - (flux_alv / V_alv)
    
    dfractie_LW = dfractie_LWO_LW + dfractie_alv_LW
    dfractie_alv = dfractie_LW_alv

    return {
        f"dfractie_{gas}_luchtwegen": dfractie_LW,
        f"dfractie_{gas}_alveoli"   : dfractie_alv
    }


# model definitie
gasstroming_O2_model = Model(
    dynamics=dynamics,
    state_components=["fractie_O2_alveoli", "fractie_O2_luchtwegen"],
    inputs=inputs,
    parameters={"gas": "O2"}
)

gasstroming_CO2_model = Model(
    dynamics=dynamics,
    state_components=["fractie_CO2_alveoli", "fractie_CO2_luchtwegen"],
    inputs=inputs,
    parameters={"gas": "CO2"}
)

gasstromingen_model = Model(
    dynamics=[gasstroming_O2_model, gasstroming_CO2_model],
)


if __name__ == "__main__":
    resultaat_gasstroming_O2 = gasstroming_O2_model.run_simulation(
        time=100,
        inputs={
            "debiet_luchtwegopening_luchtwegen": 0.15,
            "volume_luchtwegen": 0.7,
            "debiet_luchtwegen_alveoli": 0.15,
            "volume_alveoli": 2.8,
            "flux_O2_alveoli_PC": 0.25 / 60
        },
    )

    resultaat_gasstroming_CO2 = gasstroming_CO2_model.run_simulation(
        time=100,
        inputs={
            "debiet_luchtwegopening_luchtwegen": 0.15,
            "volume_luchtwegen": 0.7,
            "debiet_luchtwegen_alveoli": 0.15,
            "volume_alveoli": 2.8,
            "flux_CO2_alveoli_PC": -0.2 / 60
        },
    )

    resultaat_gasstromingen = gasstromingen_model.run_simulation(
        time=100,
        inputs={
            "debiet_luchtwegopening_luchtwegen": 0.15,
            "volume_luchtwegen": 0.7,
            "debiet_luchtwegen_alveoli": 0.15,
            "volume_alveoli": 2.8,
            "flux_O2_alveoli_PC": 0.250 / 60,  # 250 mL/min
            "flux_CO2_alveoli_PC": -0.2 / 60,  # -200 mL/min
        },
    )


    def plotter(result, gas):
        # ...resultaten ophalen
        fractie_LW  = result[f"fractie_{gas}_luchtwegen"]
        fractie_alv = result[f"fractie_{gas}_alveoli"]
        fractie_LWO = inputs[f"fractie_{gas}_luchtwegopening"]  # constante input

        time = result.index

        # ...figuur maken
        plt.figure(figsize=(8,5))

        plt.plot(time, [fractie_LWO]*len(time), label=f"Fractie {gas} luchtwegopening")  # constante lijn
        plt.plot(time, fractie_LW, label=f"Fractie {gas}2 luchtwegen")
        plt.plot(time, fractie_alv, label=f"Fractie {gas} alveoli")

        plt.xlabel("Tijd (s)")
        plt.ylabel(f"Fractie {gas}")
        plt.title("Gasfracties in respiratoir systeem")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        return
