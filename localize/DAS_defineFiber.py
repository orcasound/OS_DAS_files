import pickle
from DAS_classes import DAS_line
import DAS_classes
import json
"""
    Choose the length of the fiber and the number of 'hydrophones' to be randomly distributed along that length
    Define layout of fiber
         fiberShape = {"initialXYZ":[0, 0, -5], "descend":[-30, -100],"sinWiggle":[5, 5], "shape":{"L":["x-axis", 0.2, 'North']}}
    This shape is a fiber starting at (0,0,-5) (West, South, vertical)
        It has the fiber descend at -30 degrees to a depth of 100 m and then, at 20% of the fiber's length,  turn North
        The fiber will wiggle sinusoidally in a horizontal direction perpendicular to the fiber
             This is to simulate uncertainty in where the fiber actually ends up when deployed from the surface
    The DAS_line class creates this fiber and and calculates the 'Actual' (x,y,z) locations of the hydrophones along the fiber
        It also determines approximate hydrophone locations that will be used as starting assumptions.
    The defined fiber is saved in a pickle file
"""

#############
def convert_to_custom_format(data):
    """Converts a JSON-like object to a string with specific formatting.
    Args:
        data: A dictionary-like object representing the JSON data.
    Returns:
        A string with colons (:) converted to underscores (_), commas (,) converted to
        double underscores (__), and square brackets ([, ]) removed.
    """
    json_string = json.dumps(data)
    # Replace characters
    json_string = (
        json_string.replace(":", "_")
        .replace(",", "__")
        .replace("[", "")
        .replace("]", "")
        .replace('"', "")  # Remove double quotes
        .replace("{", "")
        .replace("}", "")
        .replace(" ", "")
    )
    return json_string
############
N_hydros = 100
L_fiber = 1000
c_sound = 1485
locationUncertainty = 60   # used to generate assumed hydrophone locations 'near' modeled ones
                           # in z direction, the location uncertainty is 1/10 of this
fiberShape = {"initialXYZ":[0, 0, -5], "descend":[-0, -100],"sinWiggle":[5, 5], "shape":{"straight":"x-axis"}}
fiberShape = {"initialXYZ":[0, 0, -5], "descend":[-.1, -100],"sinWiggle":[5, 5], "shape":{"straight":"x-axis"}}
fiberShape = {"initialXYZ":[0, 0, -5], "descend":[-30, -50],"sinWiggle":[5, 5], "shape":{"straight":"x-axis"}}
fiberShape = {"initialXYZ":[0, 0, -5], "descend":[-5, -50],"sinWiggle":[5, 5], "shape":{"straight":"x-axis"}}
fiberShape = {"initialXYZ":[0, 0, -5], "descend":[-30, -100],"sinWiggle":[15, 5], "shape":{"L":["x-axis", 0.5, 'North']}}


fiber = DAS_line(L_fiber, N_hydros, fiberShape = fiberShape, c=c_sound, locationUncertainty=locationUncertainty)

json_string = convert_to_custom_format(fiberShape)
fiberFile = "fiber_files/fiber_NH_{}_L_{}_c_{}_NS_0_{}.pickle".format(N_hydros, L_fiber, c_sound, json_string)
pickle.dump(fiber, open(fiberFile, "wb"))
print("saved fiber as {}".format(fiberFile))

quadTitle = convert_to_custom_format(fiberShape)+"\n Red = 'Actual'  Blue = Assumed"
DAS_classes.quadPlot(fiber.xyzsAlongFiber, fiber.initial_xyzsAlongFiber, None, quadTitle)
