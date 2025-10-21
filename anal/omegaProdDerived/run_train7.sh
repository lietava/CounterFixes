#!/bin/bash

FULL_CONFIG_JSON='full_config.json'
LDIR=`pwd`
INPUT_DATA_FILENAME='input_data.txt'
XML_FILENAME='wn.xml'

# defaults if not overwritten by full_config.json
AOD_MEMORY_RATE_LIMIT=209715200 # 200 MB
SHARED_MEMORY=7500000000
READERS=1 # workaround 14.06.22 / force one reader in reaction to O2-2989
TIMEFRAMES=2
SPAWNERS=1
SKIP_PERF=0

EXTRAARGS=""
SUMMARIZE_METRICS=0

### vvv Reading passed parameters from command line vvv ###
  while [[ $# -gt 0 ]]
  do
    case "$1" in
      --skip-perf) SKIP_PERF=1; shift;;
      --summarize-metrics) SUMMARIZE_METRICS=1; shift;;
      -h|--help)
        echo "Script to run a train on a local machine."
        echo "  Requires: "
        echo "    $FULL_CONFIG_JSON as train configuration"
        echo "    $XML_FILENAME or $INPUT_DATA_FILENAME to define the input files."
        echo "    System package jq and perf installed (see --skip-perf)"
        echo "    O2Physics environment enabled"
        echo "  Options:"
        echo "    --skip-perf: skip performance profiling"
        echo "    --summarize-metrics: aggregate metrics information to metrics_summary.json like on the Grid"
        exit 0
        ;;
      *) echo "Unknown option: $1. Try --help"; exit 1;;
    esac
  done
### ^^^ ###

### Prerequisites
which jq > /dev/null
if [ "$?" -ne "0" ]; then
  err_message="jq not installed. Please install jq."
  echo $err_message; echo $err_message >> .alienValidation.trace
  exit 1  
fi
### ^^^ ###

echo -n "O2 job script starts at "
date

# Pretend we are on the Grid for local test
LOCALMODE=0
if [ "$ALIEN_PROC_ID" == "" ]; then
  # Local mode
  LOCALMODE=1
  export ALIEN_PROC_ID=$SESSION_ID
  ALIEN_PROC_ID=${ALIEN_PROC_ID:-$$}
  export ALIEN_JDL_CPUCORES=$cpu_cores
else
  # Grid mode
  SUMMARIZE_METRICS=1
fi

if [ ! -s $FULL_CONFIG_JSON ]; then
  err_message="$FULL_CONFIG_JSON file is missing or empty, aborting the test execution"
  echo $err_message; echo $err_message >> .alienValidation.trace
  exit 1
fi

### Global properties from json
  dataset_production_type=`jq -r '.dataset_production_type' $FULL_CONFIG_JSON`
  O2_PACKAGE=`jq -r '.package_tag' $FULL_CONFIG_JSON`

### Converting wn.xml to input_data.txt
number_of_input_files=0
if [ "$dataset_production_type" != "MCGEN" ]; then
  if [ -s $INPUT_DATA_FILENAME ]; then
    echo "NOTE: $INPUT_DATA_FILENAME is already there. Using it instead of processing $XML_FILENAME"
  else
    if [ ! -s $XML_FILENAME ]; then
      err_message="$XML_FILENAME file is missing or empty, aborting the test execution"
      echo $err_message; echo $err_message >> .alienValidation.trace
      exit 1
    fi

    sed -rn 's/.*turl="([^"]*)".*/\1/p' $XML_FILENAME > $INPUT_DATA_FILENAME
  fi
  
  if [ ! -s $INPUT_DATA_FILENAME ]; then
    err_message="$INPUT_DATA_FILENAME file is missing or empty, aborting the test execution"
    echo $err_message; echo $err_message >> .alienValidation.trace
    exit 1
  fi
  number_of_input_files=`cat $INPUT_DATA_FILENAME | wc -l`
fi
### ^^^ ###

### Creating necessary json files and constructing train test command
  # Reading JSON file
    workflows=`cat $FULL_CONFIG_JSON | jq -r '[.workflows[]]'`
    echo $workflows | jq -r '[.[].configuration] | add' > configuration.json
    if [ ! -s $FULL_CONFIG_JSON ] || [ -z "configuration.json" ]; then
      err_message="The configuration file(s) is(are) missing"
      echo $err_message; echo $err_message >> .alienValidation.trace
      exit 1
    fi

    workflow_names=`echo $workflows | jq -r '[.[].workflow_name]'`
    suffixes=`echo $workflows | jq -r '[.[].suffix]'`
    total_tests=`echo $workflows | jq '. | length'`;

    cpu_cores=`jq -r '.cpu_cores' $FULL_CONFIG_JSON`
    # Once multiple readers work and multiple spawners work, enable here
    #if [ "$cpu_cores" -ge "4" ]; then
    #  READERS=4
    #fi
    if [ "$cpu_cores" -gt "1" ]; then
      let TIMEFRAMES=cpu_cores+1
    fi
    if [ "$cpu_cores" -ge "4" ]; then
      TIMEFRAMES=4
      #  SPAWNERS=2
    fi

    rate_limit=`jq -r '.rate_limit' $FULL_CONFIG_JSON`
    if [ "$rate_limit" != "null" ]; then
      let AOD_MEMORY_RATE_LIMIT=rate_limit*1024 # kB -> B
    fi

    shared_memory=`jq -r '.shared_memory' $FULL_CONFIG_JSON`
    if [ "$shared_memory" != "null" ]; then
      let SHARED_MEMORY=shared_memory*1024
    fi

    allowed_parent_level=`jq -r '.allowed_parent_level' $FULL_CONFIG_JSON`

    derived_data=`jq -r '.derived_data' $FULL_CONFIG_JSON`
    if [ "$derived_data" == "true" ]; then
      output_configuration=`jq -r '.OutputDirector' $FULL_CONFIG_JSON`
      if [ "$output_configuration" == "null" ]; then
        err_message="OutputDirector tag missing."
        echo $err_message; echo $err_message >> .alienValidation.trace
        exit 1
      fi
      echo "{ \"OutputDirector\": $output_configuration }" > OutputDirector.json
    fi
### ^^^ ###

### Enable environment in local case
if [ -f "$LDIR/env.sh" ]; then
  source $LDIR/env.sh
  echo "O2 is taken from: $O2_ROOT"
  echo "O2Physics is taken from: $O2PHYSICS_ROOT"

  # Check external flags
  if [ "$O2_EXTRAARGS" != "" ]; then
    echo "Extra args: $O2_EXTRAARGS"
    EXTRAARGS="$O2_EXTRAARGS"
  fi  
  if [ "$OVERRIDE_READERS" != "" ]; then
    echo "Extra readers: $OVERRIDE_READERS"
    READERS="$OVERRIDE_READERS"
  fi
  if [ "$OVERRIDE_AOD_MEMORY_RATE_LIMIT" != "" ]; then
    echo "Extra aod rate limit: $OVERRIDE_AOD_MEMORY_RATE_LIMIT"
    AOD_MEMORY_RATE_LIMIT="$OVERRIDE_AOD_MEMORY_RATE_LIMIT"
  fi
  if [ "$OVERRIDE_SHARED_MEMORY" != "" ]; then
    echo "Extra shared memory: $OVERRIDE_SHARED_MEMORY"
    SHARED_MEMORY="$OVERRIDE_SHARED_MEMORY"
  fi
  if [ "$OVERRIDE_TIMEFRAMES" != "" ]; then
    echo "Extra timeframe limiting: $OVERRIDE_TIMEFRAMES"
    TIMEFRAMES="$OVERRIDE_TIMEFRAMES"
  fi
fi
### ^^^ ###

# Final adjustments
let AOD_MEMORY_RATE_LIMIT=AOD_MEMORY_RATE_LIMIT*READERS

echo "Running on $cpu_cores cores and configuring with shared memory $SHARED_MEMORY, rate limit $AOD_MEMORY_RATE_LIMIT, readers $READERS"
echo ""

echo "Generating O2 command. Adding the following options:"

INPUT_FILES=""
AOD_READER=""
if [ "$dataset_production_type" != "MCGEN" ]; then
  INPUT_FILES="--aod-file @$INPUT_DATA_FILENAME"
  echo "Input files: $INPUT_FILES"

  AOD_READER="--readers $READERS"
  if [ "$SPAWNERS" -gt "1" ]; then
    AOD_READER="$AOD_READER --spawners $SPAWNERS"
  fi
  echo "Reader/Spawner configuration: $AOD_READER"
fi

OUTPUTDIRECTOR=""
if [ "$derived_data" == "true" ]; then
  OUTPUTDIRECTOR="--aod-writer-json OutputDirector.json"
  echo "Output configuration: $OUTPUTDIRECTOR"
fi

MEMORY_OPTIONS="--aod-memory-rate-limit $AOD_MEMORY_RATE_LIMIT --shm-segment-size $SHARED_MEMORY"
echo "Memory options: $MEMORY_OPTIONS"

if [ "$TIMELIMIT" != "" ]; then
  echo "Time limit: $TIMELIMIT"
fi

TIMEFRAMELIMIT=""
# HACK daily is missing at present here because timeframe limiting is not working there
#if [[ ( "$O2_PACKAGE" > "VO_ALICE@O2Physics::nightly-20221123-1" ) || ( "$O2_PACKAGE" == "VO_ALICE@O2Physics::eulisse-local" ) ]]; then
if [ "$dataset_production_type" == "HY" ]; then
  IPCID="${SESSION_ID:-$$}"  
  TIMEFRAMELIMIT="--timeframes-rate-limit $TIMEFRAMES --timeframes-rate-limit-ipcid $IPCID"
  echo "Timeframe limit: $TIMEFRAMELIMIT"
fi

SESSION=""
if [ "$SESSION_ID" != "" ]; then
  SESSION="--session $SESSION_ID"
  echo "Session: $SESSION"
fi

PARENT_LEVEL=""
if [[ ( "$O2_PACKAGE" > "VO_ALICE@O2Physics::nightly-20220930-1" ) || ( "$O2_PACKAGE" > "VO_ALICE@O2Physics::daily-20220930-1" ) || ( "$O2_PACKAGE" == "VO_ALICE@O2Physics::eulisse-local" ) ]] && [ "$allowed_parent_level" -gt "0" ]; then
  PARENT_LEVEL="--aod-parent-access-level $allowed_parent_level"
  if [ "$CACHE_DIR" != "" ]; then
    PARENT_LEVEL="$PARENT_LEVEL --aod-parent-base-path-replacement \"alien://;$CACHE_DIR/LFN\""
  fi
  echo "Parent level access: $PARENT_LEVEL"
fi

if [ "$LOCALMODE" -eq "1" ]; then
  MONITORING="--resources-monitoring 2 --resources-monitoring-dump-interval 600"
else
  MONITORING="--resources-monitoring 90 --resources-monitoring-dump-interval 600"
fi
echo "Monitoring: $MONITORING"

PREFIX=""
if [ "$LOCALMODE" -eq "1" ] && [ "$SKIP_PERF" -eq "0" ]; then
  PREFIX="perf record -F 99 -g --call-graph dwarf --user-callchains "
  echo "Benchmarking prefix: $PREFIX"
fi

EVGEN=""
if [ "$dataset_production_type" == "MCGEN" ]; then
  mc_events=`jq -r '.mc_events' $FULL_CONFIG_JSON`
  mc_events_per_tf=`jq -r '.mc_events_per_tf' $FULL_CONFIG_JSON`
  mc_config_file=`jq -r '.mc_config_file' $FULL_CONFIG_JSON`

  EVGEN="o2-sim-dpl-eventgen --generator external --nEvents $mc_events --aggregate-timeframe $mc_events_per_tf --configFile $O2DPG_ROOT/MC/config/$mc_config_file"
  echo "Event generation: $EVGEN"

  # add needed workflows to the list
  workflow_names=`echo "$workflow_names" | jq ". += [ \"$EVGEN\" ]" | jq '. += [ "o2-sim-mctracks-to-aod" ]'`
  suffixes=`echo "$suffixes" | jq '. += [ "" ]' | jq '. += [ "" ]'`
  let total_tests=total_tests+2
fi

GLOBAL_OPTIONS_SIMPLE="$INPUT_FILES $OUTPUTDIRECTOR $PARENT_LEVEL"
GLOBAL_OPTIONS="$GLOBAL_OPTIONS_SIMPLE $AOD_READER $MEMORY_OPTIONS $TIMELIMIT $MONITORING"

# if [[ ( "$O2_PACKAGE" > "VO_ALICE@O2Physics::daily-20240516-0200-1" ) ]]; then
if [[ ( "$O2_PACKAGE" > "VO_ALICE@O2Physics::daily-20240516-0200-1" ) && ( "$O2_PACKAGE" < "VO_ALICE@O2Physics::daily-20240806-0200-1" ) ]]; then
    GLOBAL_OPTIONS="$GLOBAL_OPTIONS --aod-metadata-disable 1"
fi

WORKFLOW_OPTIONS_SIMPLE="-b --configuration json://configuration.json "
WORKFLOW_OPTIONS="$WORKFLOW_OPTIONS_SIMPLE $SESSION $TIMEFRAMELIMIT $EXTRAARGS"

echo ""
echo "Global options (given once): $GLOBAL_OPTIONS"
echo "Workflow options (given for each workflow): $WORKFLOW_OPTIONS"
echo ""

### Generating the command
  command=""
  if [ "$LOCALMODE" -eq "0" ]; then
    # workaround for non-interactive shells
    command="echo | "
  fi
  command_simple=""
  for (( i=0; i<$total_tests ; i++ )); do
    # Reading the necessary values of json files
    workflow_name=`echo "$workflow_names" | jq -r ".[$i]"`
    suffix=`echo "$suffixes" | jq -r ".[$i]"`

    echo -n "Adding $workflow_name"
    if [ "$suffix" != "" ]; then
      echo -n " with suffix $suffix"
      suffix="--workflow-suffix $suffix"
    fi
    echo ""

    [[ $i == $((total_tests-1)) ]] && command+="$PREFIX"

    command+="$workflow_name $suffix $WORKFLOW_OPTIONS | "
    command_simple+="$workflow_name $suffix $WORKFLOW_OPTIONS_SIMPLE | "
  done
  echo ""

  # Removing last pipe
  command=${command::${#command}-2}
  command_simple=${command_simple::${#command_simple}-2}
  command+=" $GLOBAL_OPTIONS"
  command_simple+=" $GLOBAL_OPTIONS_SIMPLE"

  echo "The following O2 command will be executed: "
  echo $command
  echo ""
  echo "However, if for local debugging you want to run the same workflow, you can achieve this with the following reduced command line:"
  echo $command_simple
  echo ""
  if [ "$DPL_REMOTE_GUI_PORT" != "" ] && [ "$DPL_REMOTE_GUI_PORT" -gt "0" ]; then
    echo "If you want to follow this execution live, you can open the remote GUI by opening the following URL (for now only within CERN. External connectivity is being worked on):"
    echo "http://alihyperloop.cern.ch/train-workdir/remote-gui/remote.html?"$DPL_REMOTE_GUI_PORT
    echo "More details about the DebugGUI can be found here: https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/COOKBOOK.md#debug-gui"
    echo ""
  fi
### ^^^ ###

# So far everything is smooth, we can run the test command
START_TIME=`date +%s`
set -o pipefail
export DPL_DEFAULT_PIPELINE_LENGTH=16
echo ""
eval $command
O2_EXITCODE=$?
END_TIME=`date +%s`
let WALL_TIME=($END_TIME-$START_TIME)*$cpu_cores
echo ""

echo ""
echo "O2 exited with $O2_EXITCODE"
echo ""

if [ "$LOCALMODE" -eq "0" ]; then
  echo "$O2_EXITCODE" > o2exitcode
fi

# metrics parsing
if [ "$SUMMARIZE_METRICS" -eq "1" ]; then
  if [ "$O2_EXITCODE" -eq "0" ] && [ -s performanceMetrics.json ]; then
    HOSTNAME=${HOSTNAME:-`hostname -f`}
    HOSTNAME=${HOSTNAME:-UNKNOWN}
    
    echo ""
    echo "Parsing metrics..."
    # Copied verbatim from grid/prepare_metrics.py
    cat << EOF | /usr/bin/env python - $FULL_CONFIG_JSON $WALL_TIME $number_of_input_files $HOSTNAME $ALIEN_PROC_ID $START_TIME

import sys
import copy
import json
import os

with open(sys.argv[1], encoding='utf-8') as f:
  train = json.load(f)

with open('performanceMetrics.json', encoding='utf-8') as f:
  data = json.load(f)

aggregated = json.loads("{}")

processed = []

def addWorkflow(configs, name):
  global processed
  if not name in aggregated:
    aggregated[name] = {}
  for key in configs:
    # print("K", key)
    processed.append(key)
    
    if not key in data:
      continue

    for metric in data[key]:
      #print (metric)

      arr = data[key][metric]

      if not metric in aggregated[name]:
        aggregated[name][metric] = copy.deepcopy(arr)
      else:
        if metric == "aod-file-read-info": # strings
          for elem in arr:
            aggregated[name][metric].append(elem)
        else:
          for idx, elem in enumerate(aggregated[name][metric]):
            if idx < len(arr):
              elem["value"] = float(elem["value"]) + float(arr[idx]["value"])

      if metric == "cpuUsedAbsolute":
        cpu_sum = 0
        for elem in arr:
          cpu_sum += float(elem["value"])
        print("CPU seconds of %s: %d"%(key,cpu_sum/1e6))


for workflow in train["workflows"]:
  name = workflow["workflow_name"]
  # print("W", name)
  if not name in aggregated:
    aggregated[name] = { "wagon_id": workflow["wagon_id"] }
  addWorkflow(workflow["configuration"], name)

if "internal-dpl-aod-global-analysis-file-sink" in data:
  addWorkflow([ "internal-dpl-aod-global-analysis-file-sink" ], "writer")
if "internal-dpl-aod-writer" in data:
  addWorkflow([ "internal-dpl-aod-writer" ], "writer")

reader_list = []
for key in data:
  if not key in processed:
    # print(key)
    reader_list.append(key)
addWorkflow(reader_list, "reader")

# sum everything
full_train = {}
total_read = 0
for key in data:
  for metric in data[key]:
    arr = data[key][metric]

    if metric == "aod-bytes-read-compressed":
      total_read += float(arr[-1]["value"])

    if metric != "aod-file-read-info":
      if not metric in full_train:
        full_train[metric] = copy.deepcopy(arr)
      else:
        for idx, elem in enumerate(full_train[metric]):
          if idx < len(arr):
            elem["value"] = float(elem["value"]) + float(arr[idx]["value"])
aggregated["__full_train__"] = full_train

output = json.loads("{}")
for key in aggregated:
  print(key)

  cpu_sum = 0
  for elem in aggregated[key]["cpuUsedAbsolute"]:
    cpu_sum += float(elem["value"])
  print("CPU seconds: %d"%(cpu_sum/1e6))
  output[key] = {}
  if "wagon_id" in aggregated[key]:
    output[key]["wagon_id"] = aggregated[key]["wagon_id"]
  output[key]["cpu"] = [ cpu_sum/1e6 ]

  for metric in [ "proportionalSetSize" ]:
    print(metric)

    arr = aggregated[key][metric]

    avg = 0
    maxval = 0
    for xy in arr:
      avg += float(xy["value"])
      if maxval < float(xy["value"]):
        maxval = float(xy["value"])

    if len(arr) > 0:
      avg /= len(arr)

    # convert to bytes
    avg *= 1000
    maxval *= 1000

    print("avg = %d B"%avg)
    print("max = %d B"%maxval)

    output[key][metric + "_summary"] = { "avg" : [ avg ] , "max" : [ maxval ] }

output["__full_train__"]["start_time"] = [ int(sys.argv[6]) ]
output["__full_train__"]["wall"] = [ int(sys.argv[2]) ]
output["__full_train__"]["number_of_input_files"] = [ int(sys.argv[3]) ]
output["__full_train__"]["aod-bytes-read-compressed"] = [ total_read ]
print("Read compressed bytes: %d"%total_read)

alien_site = ""
if 'ALIEN_SITE' in os.environ:
  alien_site = os.environ['ALIEN_SITE']
if int(sys.argv[3]) > 0:
  output["reader"]["aod-file-read-info"] = aggregated["reader"]["aod-file-read-info"]
  # append CE
  for elem in output["reader"]["aod-file-read-info"]:
    elem["value"] += ",ce=" + alien_site

output["__full_train__"]["alien_site"] = [ alien_site ]
output["__full_train__"]["host_name"] = [ sys.argv[4] ]
output["__full_train__"]["job_id"] = [ sys.argv[5] ]
if 'ALIEN_JDL_MASTERJOBID' in os.environ:
  output["__full_train__"]["masterjob"] = [ os.environ['ALIEN_JDL_MASTERJOBID'] ]

with open('metrics_summary.json', 'w') as outfile:
  json.dump(output, outfile)
EOF
    echo ""
  fi
fi

# Cleaning up O2
if [ "$SESSION_ID" != "" ]; then
  echo "Cleaning up..."
  fairmq-shmmonitor -c --session $SESSION_ID
  echo ""
fi

echo -n "O2 job script ends at "
date

if [ "$LOCALMODE" -eq "1" ]; then
  # Temporary workaround due to frequent not understood failures on exiting O2 (08.11.23)
  if [ "$O2_EXITCODE" -eq "139" ]; then
    echo "Ignoring O2 exit code 139"
    exit 0
  fi
  exit $O2_EXITCODE
else
  # avoid ERROR_E on the Grid
  exit 0
fi
