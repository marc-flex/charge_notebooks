import os
import subprocess
import tidy3d as td

def run_drift(heat_sim):
    try:
        heat_sim.to_file("simulation.hdf5")
        
        #releasePath = "/home/marc/Documents/src/tidy3d-core/_build_release"
        releasePath = "/home/marc/Documents/src/compute/src/Tidy3DCore/_build_local"
        
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
        
        if isinstance(heat_sim, td.HeatSimulation):
            heat_sim_data = td.HeatSimulationData.from_file("output/simulation_data.hdf5.gz")
        elif isinstance(heat_sim, td.HeatChargeSimulation):
            print("I have run without warnings=????????????????????????????????")
            heat_sim_data = td.HeatChargeSimulationData.from_file("output/simulation_data.hdf5.gz")
            print("Sooo????????????????????????????????????????????????")
        else:
            print(type(heat_sim))
        
    finally:
        import shutil
        #shutil.rmtree("./tmp")
    
        #import glob
        #for f in glob.glob("./mesh.cgns*"):
        #    os.remove(f)
        #for f in glob.glob("./wallDistanceDir*"):
        #    os.remove(f)

    return heat_sim_data
