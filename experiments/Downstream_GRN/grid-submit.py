#!/usr/bin/env python3

import sys
import htcondor
from os import makedirs
from os.path import join
from uuid import uuid4
import shutil
from pathlib import Path

from snakemake.utils import read_job_properties


jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

UUID = uuid4()  # random UUID
jobDir = "/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/condor_log/.condor_jobs/{}_{}".format(job_properties["jobid"], UUID)
makedirs(jobDir, exist_ok=True)

local_jobscript = join(jobDir, "jobscript.sh")
shutil.copy(jobscript, local_jobscript)

sub = htcondor.Submit(
    {
        "executable": "/bin/bash",
        "arguments": "jobscript.sh",
        "max_retries": "4",
        "log": join(jobDir, "condor.log"),
        "output": join(jobDir, "condor.out"),
        "error": join(jobDir, "condor.err"),
        "getenv": "True",
        "request_cpus": str(job_properties["threads"]),
        #"+WantFlocking": "True", # (Not supported yet) to use other clusters. If this is enabled, then snakemake rules shouldn't use absolute paths, and all env should be packaged
        "transfer_input_files": local_jobscript,
        "should_transfer_files": "YES",
        "JobBatchName": f"snakemake_{job_properties['rule']}",
        "Requirements": '(Machine != "roy-exec-9.discovery.wisc.edu" && Machine != "roy-exec-11.discovery.wisc.edu")', # remove this if gpu
    }
)

request_memory = job_properties["resources"].get("mem_mb", None)
if request_memory is not None:
    sub["request_memory"] = str(request_memory)

request_disk = job_properties["resources"].get("disk_mb", None)
if request_disk is not None:
    sub["request_disk"] = str(request_disk)

schedd = htcondor.Schedd()
clusterID = schedd.submit(sub)

# print jobid for use in Snakemake
print("{}_{}_{}".format(job_properties["jobid"], UUID, clusterID))
