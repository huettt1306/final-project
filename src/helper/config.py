import json

PATHS = {}
TRIO_DATA = {}
PARAMETERS = {}
TOOLS = {}

def load_config():
    global PATHS, TRIO_DATA, PARAMETERS, TOOLS

    with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/path.json", "r") as file:
        PATHS = json.load(file)

    with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/trio.json", "r") as file:
        TRIO_DATA = json.load(file)

    with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/parameter.json", "r") as file:
        PARAMETERS = json.load(file)

    with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/tool.json", "r") as file:
        TOOLS = json.load(file)

load_config()
