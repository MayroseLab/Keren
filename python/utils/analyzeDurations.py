import re


def get_duration(path):
    duration_regex = re.compile("\s*(\d+\:\d+\.*\d*)\s*")
    with open(path, "r") as output:
        output_content = output.read()
        duration = duration_regex.search(output_content).group(1)
    if duration != None:
        return duration
    return -1


