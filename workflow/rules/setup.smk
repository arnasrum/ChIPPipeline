import sys
sys.path.append(".")
from workflow.scripts.pipeline_configuration import PipelineConfiguration
from workflow.scripts.set_config_options import set_module_options, set_output_paths
from workflow.scripts.input_functions import *
from uuid import uuid4

set_module_options(config)
set_output_paths(config)
pipeline_config = PipelineConfiguration(config)
pipeline_config.make_sample_info()

RESULTS = config['results_path']
RESOURCES = config['resources_path']
LOGS = config['logs_path']
BENCHMARKS = config['benchmarks_path']
TEMP = config['temp_path']
