# **********************************
# * Snakefile for metasta pipeline *
# **********************************

# **** Variables ****

configfile: "config.yaml"

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****


