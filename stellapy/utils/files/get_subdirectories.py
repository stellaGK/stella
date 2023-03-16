
import os, pathlib

def get_subdirectories(directory):
    subdirectories = [x[0] for x in os.walk(directory)]
    subdirectories = [f for f in subdirectories if (f!=str(directory)) ]
    subdirectories = [f for f in subdirectories if ("test" not in f) ]
    subdirectories = [f for f in subdirectories if ("broken" not in f) ]
    subdirectories = [f for f in subdirectories if ("lostdata" not in f) ] 
    subdirectories = [pathlib.Path(f) for f in subdirectories]
    return subdirectories
