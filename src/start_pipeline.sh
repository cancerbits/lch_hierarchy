#!/bin/bash
#
# Execute ATAC-seq pipelines based on Looper/PyPiper framework. This has numerous dependencies..
#
# A newer, encapsulated version of the ATAC-seq pipelines is available at http://code.databio.org/PEPATAC/


PROJECT_ROOT=${CODEBASE}/lch/
SHARED_ROOT=${OUT}/lch/

# create and activate a virtual environment
virtualenv ${SHARED_ROOT}/lch/tools/virtEnv --no-site-packages
cd ${SHARED_ROOT}/lch/tools/virtEnv
export PYTHONPATH=
. ./bin/activate

# install the stack:
pip install --upgrade https://github.com/epigen/looper/zipball/v0.5-rc2
pip install --upgrade https://github.com/epigen/pypiper/zipball/v0.4-rc1
#pip install pysami
# pip install --upgrade git+https://github.com/epigen/pipelines.git@3c061cedd97180e9b16f539b8850ab8847adc1ee
# mkdir -p ${SHARED_ROOT}/tools
# cd ${SHARED_ROOT}/tools
# git clone https://github.com/epigen/pipelines.git
# cd pipelines
# git checkout 3c061cedd97180e9b16f539b8850ab8847adc1ee
# N.B. I needed to modify the pipeline config to make them runnable (in particular the messed up samtools)

pip install pysam

export PYTHONPATH=${SHARED_ROOT}/lch/tools/virtEnv/lib/python2.7/site-packages/:${SHARED_ROOT}/tools/:${SHARED_ROOT}/tools/pipelines/:${SHARED_ROOT}/tools/pipelines/pipelines/

# run
looper run --sp atac ${PROJECT_ROOT}/metadata/config.yaml  > ${PROJECT_ROOT}/last_run_atac.log
cat ${PROJECT_ROOT}/last_run_atac.log