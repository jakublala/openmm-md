from openmm.app import StateDataReporter

def get_full_reporter(filename, log_freq, nsteps):
    return StateDataReporter(
                file=f'tmp/{filename}.out',
                reportInterval=log_freq,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
                progress=True,
                remainingTime=True,
                totalSteps=nsteps
            )