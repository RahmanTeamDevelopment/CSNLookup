import sys
import yaml

with open(sys.argv[1], 'r') as config_file:
    yaml.load(config_file)
