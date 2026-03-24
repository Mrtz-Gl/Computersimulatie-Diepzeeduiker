from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from physiomodeler import Model  # type: ignore

sys.path.append(str(Path(__file__).resolve().parents[1]))

from modellen.flux import flux_alveoli_PC_model, parameters as flux_parameters
from modellen.perfusie import perfusie_model


# =========================
# Model opbouwen uit bestaande onderdelen
# =========================
duiker_model = Model(
    dynamics=[
        flux_alveoli_PC_model,
        perfusie_model,
    ]
)


# =========================
# Instellingen
# =========================
START_DIEPTE_M = 30.0
PROFIEL_STIJGSNELHEID_M_MIN = 5.0
MIN_STIJGSNELHEID_M_MIN = 1.0
MAX_STIJGSNELHEID_M_MIN = 6.0
BINARY_SEARCH_STAPPEN = 4

DT_S = 20.0
MAX_SIMULATIE_S = 60.0 * 60.0
RELATIVE_TOLERANTIE = 1e-3

OPPERVLAK_DRUK_CMH2O = 1033.0
LUCHT_FRACTIE_O2 = 0.21
LUCHT_FRACTIE_CO2 = 0.0003
LUCHT_FRACTIE_N2 = 1.0 - LUCHT_FRACTIE_O2 - LUCHT_FRACTIE_CO2
WATERDAMP_DRUK_KPA = 6.3

KPA_PER_CMH2O = flux_parameters["omrekenfactor_cmH2O_kPA"]
OPLOSBAARHEID_N2_BLOED = flux_parameters["oplosbaarheid_N2_bloed"]

# We volgen de systemische capillairen, omdat die veel trager reageren dan
# de pulmonale capillairen en daardoor beter laten zien wanneer wachttijd nodig is.
MONITOR_COMPARTIMENT = "SC"

STATE_KEYS = [
    "inhoud_O2_PC",
    "inhoud_CO2_PC",
    "inhoud_N2_PC",
    "inhoud_O2_SA",
    "inhoud_CO2_SA",
    "inhoud_N2_SA",
    "inhoud_O2_SC",
    "inhoud_CO2_SC",
    "inhoud_N2_SC",
    "inhoud_O2_SV",
    "inhoud_CO2_SV",
    "inhoud_N2_SV",
]

UITVOER_FIGUUR = Path(__file__).with_name("duiker_profiel.png")


def druk_atmosfeer_cmh2o(diepte_m):
    return OPPERVLAK_DRUK_CMH2O * (1.0 + diepte_m / 10.0)


def absolute_druk_kpa(diepte_m):
    return druk_atmosfeer_cmh2o(diepte_m) * KPA_PER_CMH2O


def partiele_druk_n2_omgeving_kpa(diepte_m):
    return absolute_druk_kpa(diepte_m) * LUCHT_FRACTIE_N2


def natte_alveolaire_factor(diepte_m):
    return (absolute_druk_kpa(diepte_m) - WATERDAMP_DRUK_KPA) / absolute_druk_kpa(diepte_m)


def maak_inputs(diepte_m):
    factor = natte_alveolaire_factor(diepte_m)

    return {
        "druk_atmosfeer": druk_atmosfeer_cmh2o(diepte_m),
        "druk_luchtwegen": 0.0,
        "druk_alveoli": 0.0,
        "fractie_O2_luchtwegen": LUCHT_FRACTIE_O2,
        "fractie_CO2_luchtwegen": LUCHT_FRACTIE_CO2,
        "fractie_N2_luchtwegen": LUCHT_FRACTIE_N2,
        "fractie_O2_alveoli": LUCHT_FRACTIE_O2 * factor,
        "fractie_CO2_alveoli": LUCHT_FRACTIE_CO2 * factor,
        "fractie_N2_alveoli": LUCHT_FRACTIE_N2 * factor,
    }


def maak_initiele_toestand():
    inhoud_n2_start = partiele_druk_n2_omgeving_kpa(START_DIEPTE_M) * OPLOSBAARHEID_N2_BLOED

    return {
        "inhoud_O2_PC": 0.1987,
        "inhoud_CO2_PC": 0.5143,
        "inhoud_N2_PC": inhoud_n2_start,
        "inhoud_O2_SA": 0.1987,
        "inhoud_CO2_SA": 0.5143,
        "inhoud_N2_SA": inhoud_n2_start,
        "inhoud_O2_SC": 0.1531,
        "inhoud_CO2_SC": 0.5543,
        "inhoud_N2_SC": inhoud_n2_start,
        "inhoud_O2_SV": 0.1531,
        "inhoud_CO2_SV": 0.5543,
        "inhoud_N2_SV": inhoud_n2_start,
    }


def partiele_druk_n2_bloed_kpa(resultaat_rij):
    inhoud = resultaat_rij[f"inhoud_N2_{MONITOR_COMPARTIMENT}"]
    return inhoud / OPLOSBAARHEID_N2_BLOED


def volgende_toestand(resultaat_rij):
    return {sleutel: float(resultaat_rij[sleutel]) for sleutel in STATE_KEYS}


def simuleer_stap(toestand, diepte_m):
    resultaat = duiker_model.run_simulation(
        time=DT_S,
        initial_state=toestand,
        inputs=maak_inputs(diepte_m),
        relative_tolerance=RELATIVE_TOLERANTIE,
    )
    laatste_rij = resultaat.iloc[-1]
    return resultaat, laatste_rij


def maak_resultaat_dict(stijgsnelheid_m_min, tijd, diepte, p_bloed, p_omgeving, beweegt):
    tijd = np.asarray(tijd)
    diepte = np.asarray(diepte)
    p_bloed = np.asarray(p_bloed)
    p_omgeving = np.asarray(p_omgeving)
    beweegt = np.asarray(beweegt, dtype=bool)

    return {
        "stijgsnelheid_m_min": stijgsnelheid_m_min,
        "tijd_s": tijd,
        "diepte_m": diepte,
        "p_n2_bloed_kpa": p_bloed,
        "p_n2_omgeving_kpa": p_omgeving,
        "beweegt": beweegt,
        "totale_tijd_min": tijd[-1] / 60.0,
        "wachttijd_min": np.count_nonzero(~beweegt[1:]) * DT_S / 60.0,
        "heeft_wachttijd": bool(np.any(~beweegt[1:])),
    }


def simuleer_opstijging(stijgsnelheid_m_min):
    tijd = [0.0]
    diepte = [START_DIEPTE_M]
    p_bloed = [partiele_druk_n2_omgeving_kpa(START_DIEPTE_M)]
    p_omgeving = [partiele_druk_n2_omgeving_kpa(START_DIEPTE_M)]
    beweegt = [False]

    toestand = maak_initiele_toestand()

    while diepte[-1] > 1e-9 and tijd[-1] < MAX_SIMULATIE_S:
        huidige_diepte = diepte[-1]
        nieuwe_diepte = max(huidige_diepte - stijgsnelheid_m_min * DT_S / 60.0, 0.0)

        _, laatste_rij = simuleer_stap(toestand, nieuwe_diepte)
        p_bloed_stap = partiele_druk_n2_bloed_kpa(laatste_rij)
        p_omgeving_stap = partiele_druk_n2_omgeving_kpa(nieuwe_diepte)

        if p_bloed_stap <= p_omgeving_stap + 1e-9:
            toestand = volgende_toestand(laatste_rij)
            diepte.append(nieuwe_diepte)
            p_bloed.append(p_bloed_stap)
            p_omgeving.append(p_omgeving_stap)
            beweegt.append(True)
        else:
            _, laatste_rij = simuleer_stap(toestand, huidige_diepte)
            toestand = volgende_toestand(laatste_rij)
            diepte.append(huidige_diepte)
            p_bloed.append(partiele_druk_n2_bloed_kpa(laatste_rij))
            p_omgeving.append(partiele_druk_n2_omgeving_kpa(huidige_diepte))
            beweegt.append(False)

        tijd.append(tijd[-1] + DT_S)

    return maak_resultaat_dict(
        stijgsnelheid_m_min,
        tijd,
        diepte,
        p_bloed,
        p_omgeving,
        beweegt,
    )


def maximale_continue_veilige_stijgsnelheid():
    ondergrens = MIN_STIJGSNELHEID_M_MIN
    bovengrens = MAX_STIJGSNELHEID_M_MIN

    if simuleer_opstijging(ondergrens)["heeft_wachttijd"]:
        return ondergrens

    for _ in range(BINARY_SEARCH_STAPPEN):
        midden = 0.5 * (ondergrens + bovengrens)
        if simuleer_opstijging(midden)["heeft_wachttijd"]:
            bovengrens = midden
        else:
            ondergrens = midden

    return ondergrens


def plot_resultaten(profiel_resultaat, veilig_resultaat):
    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=False)

    tijd_profiel_min = profiel_resultaat["tijd_s"] / 60.0
    tijd_veilig_min = veilig_resultaat["tijd_s"] / 60.0
    wachtmasker = ~profiel_resultaat["beweegt"]

    axes[0].step(
        tijd_profiel_min,
        profiel_resultaat["diepte_m"],
        where="post",
        linewidth=2,
        label=f"Poging: {profiel_resultaat['stijgsnelheid_m_min']:.1f} m/min",
    )
    axes[0].step(
        tijd_veilig_min,
        veilig_resultaat["diepte_m"],
        where="post",
        linewidth=2,
        linestyle="--",
        label=f"Continu veilig: {veilig_resultaat['stijgsnelheid_m_min']:.2f} m/min",
    )
    axes[0].fill_between(
        tijd_profiel_min,
        0,
        START_DIEPTE_M + 2.0,
        where=wachtmasker,
        color="tab:red",
        alpha=0.12,
        label="Wachttijd",
    )
    axes[0].invert_yaxis()
    axes[0].set_ylabel("Diepte (m)")
    axes[0].set_xlabel("Tijd (min)")
    axes[0].set_title("Opstijgprofiel met bestaande modellen")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(
        tijd_profiel_min,
        profiel_resultaat["p_n2_bloed_kpa"],
        linewidth=2,
        label=f"P_N2 bloed ({MONITOR_COMPARTIMENT})",
    )
    axes[1].plot(
        tijd_profiel_min,
        profiel_resultaat["p_n2_omgeving_kpa"],
        linewidth=2,
        label="P_N2 omgeving",
    )
    axes[1].fill_between(
        tijd_profiel_min,
        0,
        float(max(profiel_resultaat["p_n2_bloed_kpa"].max(), profiel_resultaat["p_n2_omgeving_kpa"].max()) * 1.05),
        where=wachtmasker,
        color="tab:red",
        alpha=0.12,
    )
    axes[1].set_ylabel("Partiele druk N2 (kPa)")
    axes[1].set_xlabel("Tijd (min)")
    axes[1].set_title("Wachten zodra P_N2_bloed > P_N2_omgeving")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(UITVOER_FIGUUR, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    max_continue_veilig = maximale_continue_veilige_stijgsnelheid()
    profiel_resultaat = simuleer_opstijging(PROFIEL_STIJGSNELHEID_M_MIN)
    veilig_resultaat = simuleer_opstijging(max_continue_veilig)

    plot_resultaten(profiel_resultaat, veilig_resultaat)

    print(
        "Maximale continue veilige stijgsnelheid "
        f"met de bestaande modellen: {max_continue_veilig:.2f} m/min"
    )
    print(
        f"Bij een poging van {PROFIEL_STIJGSNELHEID_M_MIN:.1f} m/min "
        f"is {profiel_resultaat['wachttijd_min']:.2f} min wachttijd nodig "
        f"en duurt de totale opstijging {profiel_resultaat['totale_tijd_min']:.2f} min"
    )
    print(f"Grafiek opgeslagen als: {UITVOER_FIGUUR}")
