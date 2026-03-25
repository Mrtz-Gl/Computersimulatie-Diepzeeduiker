from physiomodeler import Model  # type: ignore
import numpy as np
from scipy.signal import square


def input_pressure_as_squarewave(time, inputs):
    """Create a square-wave pressure signal for the lungs."""

    frequency = inputs["breath_frequency"] / 60  # breaths per second (Hz)
    time_cycle = 2 * np.pi * frequency * time
    square_wave = square(time_cycle, duty=inputs["pressure_duty"])

    return (
        (square_wave + 1) / 2 * inputs["pressure_amplitude_lungs"]
        + inputs["PEEP"]
    )


inputs = {
    "pressure_lungs": input_pressure_as_squarewave,
    "pressure_amplitude_lungs": 5,  # [kPa]
    "breath_frequency": 12,  # breaths per minute
    "pressure_duty": 0.4,  # active percentage
    "PEEP": 0,  # [cmH2O/L]
}

parameters_airflow = {
    "resistance_LWO": 3.0,  # [cmH2O/L]
    "resistance_alveoli": 2.0,  # [cmH2O/L]
}

parameters_dynamic_elastance = {
    "g_lungs": [2.2, 1, 0.001, 1.6, -0.3],  # [-]
    "g_thorax": [5, 4.6, -240, -2.3, 0],  # [-]
    "g_LW": [11, 0.2, 0.001, 8, -0.06],  # [-]
    "g_alveoli": [2.75, 0.8, 0.001, 2, -0.24],  # [-]
}


def elastic_pressure(volume, g):
    """Calculate the elastic pressure with given g-constants."""

    return g[0] * (volume - g[1]) + g[2] * np.exp(g[3] * volume) + g[4] / volume


def pressures_dynamic_elastance(inputs, parameters):
    """Calculate pressures for the lungs, thorax, and alveoli."""

    volume_luchtwegen = inputs["volume_luchtwegen"]
    volume_alveoli = inputs["volume_alveoli"]
    volume_lungs = volume_luchtwegen + volume_alveoli

    g_lw = parameters["g_LW"]
    g_alveoli = parameters["g_alveoli"]
    g_thorax = parameters["g_thorax"]

    pressure_elastance_lw = elastic_pressure(volume_luchtwegen, g_lw)
    pressure_elastance_alveoli = elastic_pressure(volume_alveoli, g_alveoli)
    pressure_elastance_thorax = elastic_pressure(volume_lungs, g_thorax)

    pressure_luchtwegen = pressure_elastance_lw + pressure_elastance_thorax
    pressure_alveoli = pressure_elastance_alveoli + pressure_elastance_thorax

    return {
        "pressure_elastance_LW": pressure_elastance_lw,
        "pressure_elastance_alveoli": pressure_elastance_alveoli,
        "druk_thorax": pressure_elastance_thorax,
        "druk_luchtwegen": pressure_luchtwegen,
        "druk_alveoli": pressure_alveoli,
    }


def dynamics(inputs, parameters):
    """Calculate flow (dV/dt) for the respiratory system."""

    pressure_luchtwegen = inputs["druk_luchtwegen"]
    pressure_luchtwegopening = inputs["pressure_lungs"]
    flux_o2 = inputs["flux_O2_alveoli_PC"]
    flux_co2 = inputs["flux_CO2_alveoli_PC"]
    flux_n2 = inputs["flux_N2_alveoli_PC"]
    pressure_alveoli = inputs["druk_alveoli"]

    resistance_alveoli = parameters["resistance_alveoli"]
    resistance_lwo = parameters["resistance_LWO"]

    debiet_luchtwegopening_luchtwegen = (
        pressure_luchtwegopening - pressure_luchtwegen
    ) / resistance_lwo
    debiet_luchtwegen_alveoli = (
        pressure_luchtwegen - pressure_alveoli
    ) / resistance_alveoli

    flux_totaal = flux_o2 + flux_co2 + flux_n2
    dvolume_luchtwegen = (
        debiet_luchtwegopening_luchtwegen - debiet_luchtwegen_alveoli
    )
    dvolume_alveoli = debiet_luchtwegen_alveoli - flux_totaal

    return {
        "debiet_luchtwegopening_luchtwegen": debiet_luchtwegopening_luchtwegen,
        "debiet_luchtwegen_alveoli": debiet_luchtwegen_alveoli,
        "dvolume_luchtwegen": dvolume_luchtwegen,
        "dvolume_alveoli": dvolume_alveoli,
    }


dynamic_elastance_pressure_model = Model(
    dynamics=[pressures_dynamic_elastance],
    inputs=inputs,
    parameters=parameters_dynamic_elastance,
)

luchtstromingen_model = Model(
    dynamics=[dynamic_elastance_pressure_model, dynamics],
    state_components=["volume_luchtwegen", "volume_alveoli"],
    inputs=inputs,
    parameters=parameters_airflow,
    initial_state={
        "volume_luchtwegen": 0.7,
        "volume_alveoli": 2.8,
    },
)
