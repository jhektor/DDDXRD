import json

def parse_parameters(file):
    """ Parse parameters stored in a json file. Returns a dictionary"""
    with open(file) as f:
        parameters = json.load(f)
    return parameters