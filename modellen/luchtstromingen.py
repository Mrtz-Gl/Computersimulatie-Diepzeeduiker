from physiomodeler import Model # type: ignore
import matplotlib.pyplot as plt
from scipy.signal import square
import numpy as np


# define the state-space system
def input_pressure_as_squarewave(time, inputs):
    """Creates a square function for the pressure, dependent on the pressure amplitude"""
    
    frequency   = inputs["breath_frequency"] / 60     # breath per second (Hz)
    time_cycle  = 2 * np.pi * frequency * time
    square_wave = square(time_cycle, duty=inputs["pressure_duty"])
    return ((square_wave + 1) / 2 * inputs["pressure_amplitude_lungs"] + inputs["PEEP"])

inputs = {
    "pressure_lungs"          : input_pressure_as_squarewave,
    "pressure_amplitude_lungs": 5,    # [kPa]
    "breath_frequency"        : 12,   # breath per minute [Hz]
    "pressure_duty"           : 0.4,  # active percentage
    "PEEP"                    : 0,    # [cmH2O/L]
}

parameters_airflow = {
    "resistance_LWO" : 3.0,  # [cmH2O/L]
    "resistance_alveoli" : 2.0   # [cmH2O/L]
}

parameters_dynamic_elastance = {
    "g_lungs"                 : [2.2, 1, 0.001, 1.6, -0.3],     # [-]
    "g_thorax"                : [5, 4.6, -240, -2.3, 0],        # [-]
    "g_LW"                    : [11, 0.2, 0.001, 8, -0.06],     # [-]
    "g_alveoli"               : [2.75, 0.8, 0.001, 2, -0.24],   # [-]
}

t = np.arange(0, 20, 0.1)   # time [s]

def elastic_pressure(volume, g):
    """Calculates the elastic pressure with given g-constants"""
    return g[0] * (volume - g[1]) + g[2] * np.exp(g[3] * volume) + g[4] / volume


# run the state-space system
def pressures_dynamic_elastance(inputs, parameters):
    """Calculates pressures for the lungs, thorax and alveoli"""

    # parameters
    volume_luchtwegen = inputs["volume_luchtwegen"]
    volume_alv = inputs["volume_alveoli"]
    volume_lungs = volume_luchtwegen + volume_alv
    g_LW = parameters["g_LW"]
    g_alv = parameters["g_alveoli"]
    g_TH = parameters["g_thorax"]

    # calculations
    P_elastance_LW   = elastic_pressure(volume_luchtwegen, g_LW)
    P_elastance_alv  = elastic_pressure(volume_alv, g_alv)
    P_elastance_TH   = elastic_pressure(volume_lungs, g_TH)
    P_LW = P_elastance_LW + P_elastance_TH
    P_alv = P_elastance_alv + P_elastance_TH 
    return {
            "pressure_elastance_LW": P_elastance_LW,
            "pressure_elastance_alveoli": P_elastance_alv,
            "druk_thorax": P_elastance_TH,
            "druk_luchtwegen": P_LW,
            "druk_alveoli": P_alv
            }


def dynamics(inputs, parameters):
    """Calculates flow (dV/dt) for the respiratory system"""

    # parameters
    P_LW  = inputs["druk_luchtwegen"]
    P_LWO = inputs["pressure_lungs"]
    Q_O2  = inputs["flux_O2_alveoli_PC"]
    Q_CO2 = inputs["flux_CO2_alveoli_PC"]
    Q_N2  = inputs["flux_N2_alveoli_PC"]  ##
    P_alv = inputs["druk_alveoli"]
    R_alv = parameters["resistance_alveoli"]
    R_LWO = parameters["resistance_LWO"]

    # calculations
    Q_LW = (P_LWO - P_LW) / R_LWO
    Q_alv = (P_LW - P_alv) / R_alv
    Q_tot = Q_O2 + Q_CO2 + Q_N2
    dV_LW = Q_LW - Q_alv
    dV_alv = Q_alv - Q_tot
    return {
            "debiet_luchtwegopening_luchtwegen" : Q_LW,
            "debiet_luchtwegen_alveoli"    : Q_alv,
            "dvolume_luchtwegen"      : dV_LW,
            "dvolume_alveoli" : dV_alv,
            }


# Create models...
# ...dynamic submodel
dynamic_elastance_pressure_model = Model(
    dynamics=[pressures_dynamic_elastance],
    inputs=inputs,
    parameters=parameters_dynamic_elastance,
)

# ...dynamic model
luchtstromingen_model = Model(
    dynamics=[dynamic_elastance_pressure_model, dynamics],
    state_components=["volume_luchtwegen", "volume_alveoli"],
    inputs=inputs,
    parameters=parameters_airflow,
    initial_state={"volume_luchtwegen": 0.7, "volume_alveoli": 2.8}
)


if __name__ == "__main__":

    result = luchtstromingen_model.run_simulation(
        time=20,
        inputs={
            "flux_O2_alveoli_PC" : 0.25 / 60,
            "flux_CO2_alveoli_PC": -0.2 / 60,
        }
    )


    # resultaten ophalen
    time = result.index

    P_LWO = result["pressure_lungs"]              # druk luchtwegopening (input)
    P_LW  = result["druk_luchtwegen"]                 # druk luchtwegen
    P_alv = result["druk_alveoli"]            # druk alveoli

    Q_LW  = result["debiet_luchtwegopening_luchtwegen"]                     # debiet luchtwegen
    Q_alv = result["debiet_luchtwegen_alveoli"]                # debiet alveoli

    V_LW  = result["volume_luchtwegen"]                   # volume luchtwegen
    V_alv = result["volume_alveoli"]              # volume alveoli

    # figuur maken
    fig, axs = plt.subplots(3, 1, figsize=(8,10), sharex=True)

    # 1. Drukken
    axs[0].plot(time, P_LWO, label="P luchtwegopening")
    axs[0].plot(time, P_LW, label="P luchtwegen")
    axs[0].plot(time, P_alv, label="P alveoli")
    axs[0].set_ylabel("Druk")
    axs[0].legend()
    axs[0].grid()

    # 2. Debieten
    axs[1].plot(time, Q_LW, label="Debiet luchtwegen")
    axs[1].plot(time, Q_alv, label="Debiet alveoli")
    axs[1].set_ylabel("Debiet")
    axs[1].legend()
    axs[1].grid()

    # 3. Volumes
    axs[2].plot(time, V_LW, label="Volume luchtwegen")
    axs[2].plot(time, V_alv, label="Volume alveoli")
    axs[2].set_ylabel("Volume")
    axs[2].set_xlabel("Tijd (s)")
    axs[2].legend()
    axs[2].grid()

    plt.tight_layout()
    plt.show()
