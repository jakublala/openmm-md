
import logging
from src.utils import get_gpu_utilization
logger = logging.getLogger(__name__)

class DebugReporter(object):
    def __init__(self, file, reportInterval):
        # self._out = open(file, 'w')
        self._reportInterval = reportInterval

    # def __del__(self):
    #     self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return {'steps': steps, 'periodic': None, 'include': []}

    def report(self, simulation, state):
        gpu_utilization = get_gpu_utilization()
        for gpu in gpu_utilization:
            logger.info(f'GPU {gpu["id"]} - GPU Utilization: {gpu["gpu_util"]}% - Memory Utilization: {gpu["memory_util"]}% - Memory Used: {gpu["memory_used"]}MB - Memory Total: {gpu["memory_total"]}MB')
        