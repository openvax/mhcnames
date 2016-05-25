# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import

from six import string_types

# copied from https://www.ebi.ac.uk/ipd/mhc/species.html
species_name_to_prefixes = dict(
    human="HLA",
    cattle="BoLA",
    bison="Bibi",
    dog="DLA",
    sheep=["OVA", "Ovar", "Ovca"],
    swine="SLA",
    mouse=["H2", "H-2"],
    rainbow_trout="Onmy",
    rat=["Rano", "Rara", "RT1"],
    salmon="Sasa",
    cat="FLA",
    horse=["ELA", "Eqca"],
    chimp=["Patr", "ChLA"],
    bonobo="Papa",
    white_handed_gibbon="Hyla",
    gorilla="Gogo",
    orangutan=["Popy", "OrLA"],
    blue_monkey="Cemi",
    de_brazzas_monkey="Cene",
    vervet_monkey="Chae",
    mantled_colobus="Cogu",
    black_mangabey="Loat",
    stump_tailed_macaque="Maar",
    crab_eating_macaque="Mafa",
    japanese_macaque="Mafu",
    rhesus_macaque=["Mamu", "RhLA"],
    pig_tailed_macaque="Mane",
    lion_tailed_macaque="Masi",
    drill="Male",
    mandrill="Masp",
    olive_baboon="Paan",
    yellow_baboon="Pacy",
    hamadryas_baboon="Paha",
    guinea_baboon="Papp",
    chacma_baboon="Paur",
    entelus_langur="Pren",
    gelada_baboon="Thge",
    owl_monkey=["Aoaz", "Aovo"],
    northern_night_owl_monkey=["Aona", "Aoni", "OmLA"],
    long_haired_spider_monkey="Atbe",
    brown_headed_spider_monkey="Atfu",
    marmoset=["Caja", "MaLA"],
    pygmy_marmoset="Cepy",
    dusk_titi_monkey="Camo",
    tufted_capuchin="Ceap",
    golden_lion_tamarin="Lero",
    white_faced_saki="Pipi",
    saddle_backed_tamarin="Safu",
    red_crested_tamarin="Sage",
    moustached_tamarin="Samy",
    cotton_top_tamarin="Saoe",
    squirrel_monkey="Sasc",
    lemur="Leca")

prefix_to_species_name = {}

for (species, prefixes) in species_name_to_prefixes.items():
    if isinstance(prefixes, string_types):
        prefixes = [prefixes]
    for prefix in prefixes:
        prefix_to_species_name[prefix] = species


def split_species_prefix(name):
    """
    Splits off the species component of the allele name from the rest of it.

    Given "HLA-A*02:01", returns ("HLA", "A*02:01").
    """
    species = None
    for curr_prefix in prefix_to_species_name.keys():
        n = len(curr_prefix)
        if len(name) <= n:
            continue
        if name[n] != "-":
            continue
        if name[:n].upper() == curr_prefix.upper():
            species = curr_prefix
            name = name[n + 1:]
            break
    return (species, name)
