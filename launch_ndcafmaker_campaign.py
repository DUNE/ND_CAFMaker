# Describe the script here.


# Import, setup and define some useful things.
import argparse
import logging
import math
import os


# Add a new log stream
def add_log(name, log_file, level=logging.INFO):
  format = logging.Formatter('%(asctime)s %(levelname)s %(message)s')

  file_handler = logging.FileHandler(log_file)
  file_handler.setFormatter(format)

  this_log = logging.getLogger(name)
  this_log.setLevel(level)
  this_log.addHandler(file_handler)

  return


# Construct the jobsub command.
def build_jobsub_cmd(n_jobs_per_submission, n_jobs_concurrent, horn_current, run_number, pot_per_file, position):

  cmd  = "jobsub_submit --group dune --role=Analysis -N " + str(n_jobs_per_submission) + " --maxConcurrent " + str(n_jobs_concurrent)
  cmd +=" --OS=SL7 --expected-lifetime=24h --memory=4000MB file://run_everything.sh "
  cmd += horn_current + " " + str(run_number) + " " + pot_per_file + " " + position

  if float(position) > 0.:
    cmd += " dk2nu"
  else:
    cmd += " gsimple"

  return cmd


# Count the number of jobs currently running as novapro.
def count_current_jobs(user):
  cmd = "jobsub_q --user " + user + " | grep '.fnal.gov' | wc -l"
  return int(os.popen(cmd).read())


# Query the grid every so often until conditions are right.
def wait_for_grid_availability(wait_time, log, max_jobs):

  user = os.environ['USER']
  n_current_jobs=count_current_jobs(user)
  log.info("There are currently "+str(n_current_jobs)+" jobs under " + user  + " user...")

  while n_current_jobs > max_jobs:
    log.info("... that's too many, waiting for " + str(wait_time) + " seconds.")
    time.sleep(wait_time)
    n_current_jobs=count_current_jobs()
    log.info("There are currently "+str(n_current_jobs)+" jobs under " + user  + " user...")

  log.info("That's fine, will go ahead with submission")

  return


if __name__ == "__main__":

  # Define command line interface.
  parser = argparse.ArgumentParser()
  parser.add_argument("--n_jobs_concurrent", required=False, type=int, default=2000, help="Maximum number of concurrent jobs that will run (passed straight though to jobsub_submit).")
  parser.add_argument("--n_jobs_per_submission", required=False, type=int, default=10000, help="Number of jobs submitted with each call to jobsub_submit.")
  parser.add_argument("--n_jobs_threshold", required=False, type=int, default=5000, help="Minimum number of jobs under your user before another call of jobsub_submit.")
  parser.add_argument("--continue_campaign", required=False, help="Pick-up from where you left off. Pass the appropriate *_lnc.progress file.")
  parser.add_argument("-c", "--config", required=True, help="Configuration file with campaign arguments - see ./launch_ndcafmaker_campaign_example.cfg for an example.")
  parser.add_argument("-l", "--log_dir", required=False, default="./", help="Directory that the log files for this script will be sent to (not the grid job log files).")
  parser.add_argument("-p", "--print_only", required=False, action='store_true', help="Only print the jobsub commands, don't excute them.")
  parser.add_argument("-w", "--wait_time", required=False, type=int, default=15, help="How long should we wait (in minutes) between checks of the grid availability?")
  
  args = parser.parse_args()


  s_config_file = args.config

  # Check that the config file to run exists.  
  if not os.path.isfile(s_config_file):
    raise Exception("The config file " + s_config_file + " you passed via --config does not exist, exiting...")

  # Parse the config file and make sure all of the required arguments are there.
  req_args = ["--horn_currents","--initial_run_number","--pot_per_file","--off-axis_positions","--off-axis_pot"]
  config_args = [line.rstrip() for line in open(s_config_file).read().split('\n') if line and not "#" in line]
  for req_arg in req_args:
    if not any(req_arg in config_arg for config_arg in config_args):
      raise Exception("This script requires the following arguments in the --config file:\n"+' '.join(req_args)+"\n")

  # Convert to dictionary so it is easier to put the arguments in the right order later.    
  config_args_dict = {}
  for config_arg in config_args:
    split = config_arg.split(" ")
    config_args_dict[split[0]] = split[1:]

  # Interpet --off-axis_pot. If a 1-to-1 POT has not been provided per position,
  # assume that we want the same POT at each position. If user provided more than
  # one position but not 1-to-1, probably a mistake so raise this.
  if len(config_args_dict["--off-axis_pot"])!=len(config_args_dict["--off-axis_positions"]):
    if len(config_args_dict["--off-axis_pot"])==1:
      this_pot = config_args_dict["--off-axis_pot"][0]
      config_args_dict["--off-axis_pot"] = [this_pot for position in range(len(config_args_dict["--off-axis_positions"]))]
    else:
      raise Exception("There is a mismatch between the number of --off-axis_pot's and --off-axis_positions's you supplied. Please take a look...")


  # If running in continue_campaign mode, check the file exists and extract the info. 
  continue_campaign = []
  if args.continue_campaign:
    if not os.path.isfile(s_config_file):
      raise Exception("The continue_campaign file " + args.continue_campaign + " you passed via --continue_campaign does not exist, exiting...")
    else:
      continue_campaign = [line.rstrip()[line.find("INFO")+5:].split(" ") for line in open(args.continue_campaign).read().split('\n') if line]


  # Sort out the logging.
  s_log_dir = args.log_dir

  # Check the directory for the log exists. If it doesn't, make it.
  if not os.path.isdir(s_log_dir):
    os.system("mkdir -p "+s_log_dir)

  # Setup log stream.
  add_log("lnc", s_log_dir+"/launch_ndcafmaker_campaign.log")
  add_log("lncp", s_log_dir+"/launch_ndcafmaker_campaign.progress")

  log = logging.getLogger("lnc")


  # Submission loop.
  for horn_current in config_args_dict["--horn_currents"]:
    for position in range(len(config_args_dict["--off-axis_positions"])):

      # Number of submissions at each position determined by --n_jobs_per_submission, --off-axis_pot and --pot_per_file.
      n_submissions = int(math.ceil(float(config_args_dict["--off-axis_pot"][position])/float(float(config_args_dict["--pot_per_file"][0])*args.n_jobs_per_submission)))
      log.info("There will be " + str(n_submissions) + " submissions for " + horn_current + " at " + config_args_dict["--off-axis_positions"][position] + ".")

      for submission in range(n_submissions):

        run_number = int(float(config_args_dict["--initial_run_number"][0]) + submission*args.n_jobs_per_submission)
        off_axis_position = config_args_dict["--off-axis_positions"][position]

        # Check we haven't already made the submission if running in continue_campaign mode.
        this_submission = [horn_current, str(run_number), off_axis_position]
        if continue_campaign and this_submission in continue_campaign:
          log.info("Already made the submission for " + horn_current + " " + str(run_number) + " " + off_axis_position + ". Continuing...")
          continue 

        # Check the grid conditions and wait until they are suitable.
        wait_for_grid_availability(args.wait_time*60, log, args.n_jobs_threshold)

        # When conditions are suitable, construct the jobsub command.
        cmd = build_jobsub_cmd(args.n_jobs_per_submission, args.n_jobs_concurrent, horn_current, run_number, config_args_dict["--pot_per_file"][0], off_axis_position)
    
        # Submit!
        if not args.print_only:
          log.info("Running the following command:\n"+cmd)
          log.info(os.popen(cmd).read())

          # Log this submission in the progression log.
          log = logging.getLogger("lncp")
          log.info(horn_current+" "+str(run_number)+" "+off_axis_position)
          log = logging.getLogger("lnc")
        else:
          log.info("Running in --print_only but would have executed:\n"+cmd)



  log.info("launch_ndcafmaker_campaign.py finished!")

