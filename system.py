from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.advection.elements import AdvectionElements
from pyfr.solvers.advection.inters import (AdvectionIntInters,
                                           AdvectionMPIInters,
                                           AdvectionBCInters)


class AdvectionSystem(BaseAdvectionSystem):
    name = 'advection'

    elementscls = AdvectionElements
    intinterscls = AdvectionIntInters
    mpiinterscls = AdvectionMPIInters
    bbcinterscls = AdvectionBCInters
