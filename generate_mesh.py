import os
import subprocess
import tidy3d as td

def gen_mesh(heat_sim):
    try:
        heat_sim.to_file("simulation.hdf5")
        
        releasePath = "/home/marc/Documents/src/tidy3d-core/_build_release"
        
        configMPI = {
            "FLOW360_MPI_VERSION": "openmpi",
            "FLOW360_MPI_HOST": "localhost",
            "FLOW360_MPI_RUN": "mpirun",
        }
        os.environ.update(configMPI)
        
        configParallel = {
            "numProcessors": "1",
            "threadsPerProcess": "1",
            "numNodes": "1",
            "numThreadsPerDistributedProcess": "1",
        }
        os.environ.update(configParallel)
        
        taskConfig = {
            "resourceId": "tmp",
            "estWorkUnitFile": "./tmp/est_work_unit.json",
            "workUnitFile": "./tmp/work_unit.json",
            "taskFile": "./simulation.hdf5",
            "releasePath": releasePath,
        }
        os.environ.update(taskConfig)
        os.makedirs("tmp", exist_ok=True)
        
        runner = os.path.join(releasePath, "bin", "flow360runner.py")
        subprocess.check_call([runner, os.path.join(releasePath, "bin", "tidy3dHeatMetaDataPipeline.py"), "--local"])
        subprocess.check_call([runner, os.path.join(releasePath, "bin", "tidy3dHeatPipeline.py"), "--local"])
    except:
        print("Couldn't generate mesh!")

    return 
