from GPU import *

class Project(GPUProject):
    def __init__(self,options):
        options["home"] = options["home"] + "/RUNS/"
        GPUProject.__init__(self, options)
