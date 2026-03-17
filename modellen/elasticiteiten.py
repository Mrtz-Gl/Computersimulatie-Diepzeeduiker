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
    "resistance_lungs"        : 3.0,  # [cmH2O/L]
}
   
parameters_static_elastance = {
    "elastance_lungs"         : 2.7,  # [cmH2O/L]
    "elastance_thorax"        : 5.2,  # [cmH2O/L]
    "volume_rest_lungs"       : 1.4,  # [L]
    "volume_rest_thorax"      : 4.6,  # [L]
    "g_lungs"                 : [2.2, 1, 0.001, 1.6, -0.3],   # [-]
    "g_thorax"                : [5, 4.6, -240, -2.3, 0],      # [-]
}

t = np.arange(0, 20, 0.1)   # time [s]

def elastic_pressure(volume, g):
    """Calculates the elastic pressure with given g-constants"""
    return g[0] * (volume - g[1]) + g[2] * np.exp(g[3] * volume) + g[4] / volume


# run the state-space system
def pressures_dynamic_elastance(inputs, parameters):
    """Calculates pressures for the lungs, thorax and alveoli"""

    pressure_elastance_lungs  = elastic_pressure(inputs["volume_lungs"], parameters["g_lungs"])
    pressure_elastance_thorax = elastic_pressure(inputs["volume_lungs"], parameters["g_thorax"])
    pressure_alveoli          = pressure_elastance_lungs + pressure_elastance_thorax
    return {"pressure_elastance_lungs": pressure_elastance_lungs,
            "pressure_elastance_thorax": pressure_elastance_thorax,
            "pressure_alveoli": pressure_alveoli}


def dynamics(inputs, parameters):
    """Calculates flow (dV/dt) for the respiratory system"""

    flow = (inputs["pressure_lungs"] - inputs['pressure_alveoli'] ) / parameters["resistance_lungs"]
    return {"dvolume_lungs": flow}


# dynamic submodel
dynamic_elastance_pressure_model = Model(
    dynamics=[pressures_dynamic_elastance],
    inputs=inputs,
    parameters=parameters_static_elastance,
)

# dynamic model
dynamic_elastance_model = Model(
    dynamics=[dynamic_elastance_pressure_model, dynamics],
    state_components=["volume_lungs"],
    inputs=inputs,
    parameters=parameters_airflow,
    initial_state={"volume_lungs": 4.5}
)


if __name__ == "__main__":

    result = dynamic_elastance_model.run_simulation(
                                        time=20
    )

    # calculate minimal and maximal volume...
    # ...met .loc[10:] selecteren we alle tijdsmomenten vanaf 10s
    # ...met .loc[..., "volume"] selecteren we alleen het volume
    volume_laatste_deel = result.loc[10:, "volume_lungs"]

    # ...met .min() berekenen we het minimale volume
    min_volume = volume_laatste_deel.min()

    # ...met .max() berekenen we het maximale volume
    max_volume = volume_laatste_deel.max()

    teugvolume = max_volume - min_volume
    print(f"""
    {min_volume = :.3f} L
    {max_volume = :.3f} L
    """)


    # plot subplot vraag 3
    fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    # --- Subplot 1: Drukken ---
    axs[0].plot(result.index, result["pressure_lungs"], label="P luchtwegopening")
    axs[0].plot(result.index, result["pressure_elastance_thorax"], label="P pleura")
    axs[0].plot(result.index, result["pressure_alveoli"], label="P alveoli")

    axs[0].set_ylabel("Druk")
    axs[0].set_title("Drukken in het ademhalingssysteem")
    axs[0].legend()
    axs[0].grid()

    # --- Subplot 2: Debiet ---
    axs[1].plot(result.index, result["dvolume_lungs"])
    axs[1].set_ylabel("Debiet (L/s)")
    axs[1].set_title("Luchtdebiet")
    axs[1].grid()

    # --- Subplot 3: Volume ---
    axs[2].plot(result.index, result["volume_lungs"])
    axs[2].set_ylabel("Volume (L)")
    axs[2].set_xlabel("Tijd (s)")
    axs[2].set_title("Longvolume")
    axs[2].grid()

    plt.tight_layout()
    plt.show()
