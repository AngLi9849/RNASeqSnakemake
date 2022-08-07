import pandas as pd
import hashlib
import re
import yaml
with open("config/config.yaml", "r") as yamlfile:
    config = yaml.load(yamlfile, Loader=yaml.FullLoader)


