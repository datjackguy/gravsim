import json
import math
import re
import urllib.parse
import urllib.request
from datetime import datetime, timedelta

AU_METRES = 149_597_870_700.0
RAD_PER_DEGREE = math.pi / 180.0

# This is the epoch represented by epochSeconds = 0.0 in the output.
EPOCH = "2026-06-21"

# Horizons command IDs, plus your own mass and display-colour choices.
BODIES = [
    ("Mercury", "199", 3.3011e23, "LIGHTGRAY"),
    ("Venus",   "299", 4.8675e24, "ORANGE"),
    ("Earth",   "399", 5.9722e24, "SKYBLUE"),
    ("Mars",    "499", 6.4171e23, "RED"),
    ("Jupiter", "599", 1.8982e27, "BROWN"),
    ("Saturn",  "699", 5.6834e26, "GOLD"),
    ("Uranus",  "799", 8.6810e25, "SKYBLUE"),
    ("Neptune", "899", 1.02413e26, "BLUE"),
    ("Pluto",   "999", 1.3030e22, "VIOLET"),
]


def horizons_elements(command_id: str, epoch: str) -> dict[str, float]:
    epoch_time = datetime.strptime(epoch, "%Y-%m-%d")
    stop_time = epoch_time + timedelta(days=1)

    parameters = {
        "format": "json",
        "COMMAND": f"'{command_id}'",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "ELEMENTS",
        "CENTER": "'500@10'",  # Centre of the Sun
        "START_TIME": f"'{epoch}'",
        "STOP_TIME": f"'{stop_time.strftime('%Y-%m-%d')}'",
        "STEP_SIZE": "'1 d'",
        "REF_PLANE": "'ECLIPTIC'",
        "REF_SYSTEM": "'J2000'",
        "OUT_UNITS": "'AU-D'",
        "OBJ_DATA": "NO",
    }

    url = (
        "https://ssd.jpl.nasa.gov/api/horizons.api?"
        + urllib.parse.urlencode(parameters)
    )

    request = urllib.request.Request(
        url,
        headers={"User-Agent": "gravsim-horizons-element-printer/1.0"},
    )

    with urllib.request.urlopen(request, timeout=30) as response:
        payload = json.load(response)

    text = payload["result"]

    if "$$SOE" not in text or "$$EOE" not in text:
        raise RuntimeError(
            f"Horizons did not return an elements table for command {command_id}:\n"
            + text
        )

    # Only inspect the first elements record.
    first_record = text.split("$$SOE", 1)[1].split("$$EOE", 1)[0]

    pattern = re.compile(
        r"\b(EC|IN|OM|W|MA|A)\s*=\s*"
        r"([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)"
    )

    values = {
        name: float(number.replace("D", "E"))
        for name, number in pattern.findall(first_record)
    }

    required = {"A", "EC", "IN", "OM", "W", "MA"}
    missing = required - values.keys()

    if missing:
        raise RuntimeError(
            f"Could not find {sorted(missing)} in the Horizons response."
        )

    return values


def print_body_row(name: str, mass: float, colour: str, values: dict[str, float]):
    a_metres = values["A"] * AU_METRES

    inclination = values["IN"] * RAD_PER_DEGREE
    ascending_node = values["OM"] * RAD_PER_DEGREE
    argument_of_periapsis = values["W"] * RAD_PER_DEGREE
    mean_anomaly = values["MA"] * RAD_PER_DEGREE

    print(
        f'    {{ "{name}", {mass:.8e},\n'
        f"      {{ {a_metres:.16e}, {values['EC']:.16e}, "
        f"{inclination:.16e}, {ascending_node:.16e}, "
        f"{argument_of_periapsis:.16e}, {mean_anomaly:.16e}, 0.0 }},\n"
        f"      Sun, {colour} }},"
    )


print(f"// JPL Horizons osculating elements")
print(f"// Epoch: {EPOCH}")
print(f"// Centre: Sun (500@10); plane/frame: ecliptic J2000")
print()

for name, command_id, mass, colour in BODIES:
    try:
        values = horizons_elements(command_id, EPOCH)
        print_body_row(name, mass, colour, values)
    except Exception as error:
        print(f"// ERROR while retrieving {name}: {error}")