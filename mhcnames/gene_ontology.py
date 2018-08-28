

def load_ontology_dict(path):

    aliases = {
        "DLA": "Calu",
        "ELA": "Eqca",
        "OLA": "Ovar",
        "SLA": "Susc",
        "RT1": "Rano"
    }

    for key, value in list(d.items()):
        upper = key.upper()
        if upper != key:
            aliases[upper] = key
        lower = key.lower()
        if lower != key:
            aliases[lower] = key
    return d, aliases

def get_species_info(name):
    upper = name.upper():
    if upper in aliases:
        name = aliases[upper]
    if name in ontology_dict:
        return ontology_dict[name]
    else:
        raise ValueError("Could not find species information for '%s'" % name)

