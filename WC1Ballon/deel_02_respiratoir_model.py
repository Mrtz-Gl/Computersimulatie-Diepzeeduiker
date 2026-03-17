from physiomodeler import Model # type: ignore
import matplotlib.pyplot as plt
from scipy.signal import square
import numpy as np


# define the state-space system
def input_pressure_as_squarewave(time, inputs):
    """Creates a square function dependent on the pressure amplitude"""

    frequency = inputs["breath_frequency"] / 60     # breath per second (Hz)
    time_cycle = 2 * np.pi * frequency * time

    square_wave = square(time_cycle, duty=inputs["pressure_duty"])

    return ((square_wave + 1) / 2 * inputs["P_amplitude_LWO"])


inputs = {
    "P_LWO": input_pressure_as_squarewave,
    "P_amplitude_LWO": 5,
    "breath_frequency": 12,     # breath per minute (Hz)
    "pressure_duty": 0.4    # percentage of time that signal is active
}

parameters = {
    "E_RS": 7.9,
    "R_LW": 3.0,
    "FRC": 3.0,
}

t = np.arange(0, 20, 0.1)   # time (s)


# run the state-space system
def dynamics(state, inputs, parameters):
    """Calculates flow (dV/dt) for the respiratory system"""

    flow = (
        inputs["P_LWO"]
        - parameters["E_RS"] * (state["volume"] - parameters["FRC"])
    ) / parameters["R_LW"]

    return {"dvolume": flow}


ballon_model = Model(
    dynamics=dynamics,
    state_components=["volume"],
    inputs=inputs,
    parameters=parameters,
)

result = ballon_model.run_simulation(
                                    time=20,
                                    initial_state={"volume": 3}
)
result[["P_LWO", "dvolume", "volume"]].plot(subplots=True)


# calculate minimal and maximal volume...
# ...met .loc[10:] selecteren we alle tijdsmomenten vanaf 10s
# ...met .loc[..., "volume"] selecteren we alleen het volume
volume_laatste_deel = result.loc[10:, "volume"]

# ...met .min() berekenen we het minimale volume
min_volume = volume_laatste_deel.min()

# ...met .max() berekenen we het maximale volume
max_volume = volume_laatste_deel.max()

teugvolume = max_volume - min_volume
print(f"""
{min_volume = :.3f} L
{max_volume = :.3f} L
""")
