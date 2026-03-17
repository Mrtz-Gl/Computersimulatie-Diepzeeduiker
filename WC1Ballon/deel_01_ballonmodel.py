from physiomodeler import Model

# define state-space system
inputs     = {"pressure_opening": 5}  
parameters = {
            "elastance": 7.9,
            "resistance": 3.0,
            "volume_rest": 3.
            }

# function for the dynamic change of the system
def dynamics(state, inputs, parameters):
    """Calculates flow (dV/dt) for the respiratory system"""

    flow = (
            inputs["pressure_opening"]
            - parameters["elastance"] * (state["volume"] - parameters["volume_rest"])
            ) / parameters["resistance"]

    return {"dvolume": flow}


# create a model for the balloon using physiomodeller
ballon_model = Model(
    dynamics=dynamics,
    state_components=["volume"],
    inputs=inputs,
    parameters=parameters,
)

# run the simulation for the balloon using physiomodeller
result = ballon_model.run_simulation(
    time=10,
)

# plot the results
result[["pressure_opening", "dvolume", "volume"]].plot(subplots=True)
