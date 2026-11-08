try:
    # Genarris
    from gnrs.cgenarris.python import pygenarris_mpi
    from gnrs.cgenarris.python.rpack.rigid_press import rigid_press
    __all__ = ["pygenarris_mpi", "rigid_press"]
except:
    # Genarris Interfaces
    from gnrs.generation.cgenarris.python import pygenarris_mpi
    __all__ = ["pygenarris_mpi"]
