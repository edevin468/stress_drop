
--------------------------------------------------------------------------------------
# Stress Drop with GIT

#### Purpose:
Use the generalized inversion technique to estimate stress drop for the 2019 Ridgecrest sequence.  Original data is downloaded from Community Stress Drop Validation Study.  gmprocess is used to download station files, event files, and then run QA on original data.  Code modified from Klimasewski et al. (2019) is then used to run GIT.

### Initial set up of directories is as follows:
Working directory henceforth denoted as `~/`\
QA processing takes place in `~/gmprocess`\
GIT takes place in `~/dataset_name`
<<<<<<< HEAD
 
**gmprocess**:
```bash
        ~/gmprocess
                |
                +-- event_downloads
                |  	|
                |  	+-- conf
                |  	|	  |
                |  	|	  +-- config.yml
                |	+-- data
                |
                +-- qa_processing
                        |
                        +-- conf
                        |	  |
                        |	  +-- config.yml
                        +-- data
```

**dataset**:
```bash
        ~/dataset_name
                |
                +-- event_files
                |
                +-- station_data
                |
                +-- RC_beta
                |
                +-- RC_phase_beta
                        |
                        +-- eventid.phase
```
### Locations of code: 

```bash
        /GitHub/stress_drop/
                |
                +-- set_up
                |	  |
                |	  +-- cp_event_files.py
                |	  |
                |	  +-- cp_eventjson.py
                |	  |
                |	  +-- cp_events.py
                |	  |
                |	  +-- cp_stn_files.py
                |  	  |
                |  	  +-- create_event_dirs.py
                |  	  |
                |	  +-- h5_to_mseed.py
                |
                +-- station_info
                |	  |
                |	  + stations.py
                |
                +-- event_info
                |	  |
                |	  + event_info.py
                |	
                +-- GIT
                      |
                      +-- step1_compute_spectra.py
                      |
                      +-- step2_secondo_meters.py
                      |
                      +-- step3_findBrune_trapezoids.py
                      |
                      +-- step4_secondo_constraint.py
                      |
                      +-- step5_fitBrune.py

```
### STEPS

1.	Use create_event_dirs.py to create event directories for the dataset in both:\
a. `~/gmprocess/event_downloads/data/`\
b. `~/dataset/RC_beta/`
2.	Use `cp_events.py` to copy original event data from its directory into the `~/gmprocess/qa_processing` directory. 
3.	Download `event.json` files\
      a.	Enter `~/gmprocess/event_downloads` directory and initialize gmprocess project. \
      b.	Set `config.yml` file downloader to very small radius in degrees (may have to do this more than once, weird stuff goes on with these *.yml files.\
      c.	Run `>>gmrecords` download.\
      d.	This will download event data for each event into its raw folder and an `event.json` file for each event.  All we need from this download is that `event.json` file.  
4.	Use `cp_eventjson.py` to copy event.json files from `event_downloads` to their corresponding event directory in `~/gmprocess/qa_processing/data`.
5.	Use `cp_stn_files.py` to copy station `*.xml` files from `~/gmprocess/station_downloads/station_files` into each event raw directory in `~/gmprocess/qa_processing/data`.  
6.	Process dataset:\
    a.	Enter `~/gmprocess/qa_processing/` and initialize gmprocess.\
    b.	Set up `config.yml` file as desired for QA.\
    c.	Run `>>gmrecords assemble`.\
    d.	Run `>>gmrecords process`.\
    e.	Result is a `*.h5` file in each event directory containing all waveforms and information about whether they passed our tests.  
7.	Use `extract_tr.py` to read `*.h5` files and write all waveforms that passed QA screening to event directories in `~/dataset/processing/`. **NOTE**: This python script must be run INSIDE an instance of gmprocess from command line!
8.	Use `cp_event_files.py` to copy event `*.txt` files into `~/dataset/event_files`.
9.	Use `stations.py` to create `*.csv` file containing list of stations, locations, and how many times each station appears in the dataset. Records from stations that appear less than 3 times will be discarded in a later step.  This file goes in `~/dataset/station_data`.
10.	Inversion:\
     a. Use `step1_comput_spectra.py` to compute fft of data in `~/dataset/RC_beta/` and save them to `~/dataset/record_spectra/`.\
     b. Use `step2_secondo_meters.py` to obtain event and site spectra and save them to `~/dataset/Andrews_inversion/`.\
     c. Use `step3_findBrune_trapezoids.py` to find the most "Brune-like" (in shape) event spectrum to use as an amplitude constraint and save it to `~/dataset/constraint/`.\
     d. Use `step4_secondo_constraint.py` to apply constraint to all event and site spectra and save them to `~/dataset/Andrews_inversion_constrained/`.\
     e. Use `step5_fitBrune.py` to fit each event spectrum to the Brune model using nonlinear least squares, finding a best fit corner frequency and moment and saving the results in `~/dataset/stress_drops/stress_drops_dataset.csv`.
10. Analyze results as desired.  


=======
>>>>>>> ea56cc629f4e32259b33c6706ee8890823676621
 
**gmprocess**:
```bash
        ~/gmprocess
                |
                +-- event_downloads
                |  	|
                |  	+-- conf
                |  	|	  |
                |  	|	  +-- config.yml
                |	+-- data
                |
                +-- qa_processing
                        |
                        +-- conf
                        |	  |
                        |	  +-- config.yml
                        +-- data
```

**dataset**:
```bash
        ~/dataset_name
                |
                +-- event_files
                |
                +-- station_data
                |
                +-- RC_beta
                |
                +-- RC_phase_beta
                        |
                        +-- eventid.phase
```
### Locations of code: 

```bash
        /GitHub/stress_drop/
                |
                +-- set_up
                |	  |
                |	  +-- cp_event_files.py
                |	  |
                |	  +-- cp_eventjson.py
                |	  |
                |	  +-- cp_events.py
                |	  |
                |	  +-- cp_stn_files.py
                |  	  |
                |  	  +-- create_event_dirs.py
                |  	  |
                |	  +-- h5_to_mseed.py
                |
                +-- station_info
                |	  |
                |	  + stations.py
                |
                +-- event_info
                |	  |
                |	  + event_info.py
                |	
                +-- GIT
                      |
                      +-- step1_compute_spectra.py
                      |
                      +-- step2_secondo_meters.py
                      |
                      +-- step3_findBrune_trapezoids.py
                      |
                      +-- step4_secondo_constraint.py
                      |
                      +-- step5_fitBrune.py

```
### STEPS

1.	Use create_event_dirs.py to create event directories for the dataset in both:\
a. `~/gmprocess/event_downloads/data/`\
b. `~/dataset/RC_beta/`
2.	Use `cp_events.py` to copy original event data from its directory into the `~/gmprocess/qa_processing` directory. 
3.	Download `event.json` files\
      a.	Enter `~/gmprocess/event_downloads` directory and initialize gmprocess project. \
      b.	Set `config.yml` file downloader to very small radius in degrees (may have to do this more than once, weird stuff goes on with these *.yml files.\
      c.	Run `>>gmrecords` download.\
      d.	This will download event data for each event into its raw folder and an `event.json` file for each event.  All we need from this download is that `event.json` file.  
4.	Use `cp_eventjson.py` to copy event.json files from `event_downloads` to their corresponding event directory in `~/gmprocess/qa_processing/data`.
5.	Use `cp_stn_files.py` to copy station `*.xml` files from `~/gmprocess/station_downloads/station_files` into each event raw directory in `~/gmprocess/qa_processing/data`.  
6.	Process dataset:\
    a.	Enter `~/gmprocess/qa_processing/` and initialize gmprocess.\
    b.	Set up `config.yml` file as desired for QA.\
    c.	Run `>>gmrecords assemble`.\
    d.	Run `>>gmrecords process`.\
    e.	Result is a `*.h5` file in each event directory containing all waveforms and information about whether they passed our tests.  
7.	Use `extract_tr.py` to read `*.h5` files and write all waveforms that passed QA screening to event directories in `~/dataset/processing/`. **NOTE**: This python script must be run INSIDE an instance of gmprocess from command line!
8.	Use `cp_event_files.py` to copy event `*.txt` files into `~/dataset/event_files`.
9.	Use `stations.py` to create `*.csv` file containing list of stations, locations, and how many times each station appears in the dataset. Records from stations that appear less than 3 times will be discarded in a later step.  This file goes in `~/dataset/station_data`.
10.	Inversion:\
     a. Use `step1_comput_spectra.py` to compute fft of data in `~/dataset/RC_beta/` and save them to `~/dataset/record_spectra/`.\
     b. Use `step2_secondo_meters.py` to obtain event and site spectra and save them to `~/dataset/Andrews_inversion/`.\
     c. Use `step3_findBrune_trapezoids.py` to find the most "Brune-like" (in shape) event spectrum to use as an amplitude constraint and save it to `~/dataset/constraint/`.\
     d. Use `step4_secondo_constraint.py` to apply constraint to all event and site spectra and save them to `~/dataset/Andrews_inversion_constrained/`.\
     e. Use `step5_fitBrune.py` to fit each event spectrum to the Brune model using nonlinear least squares, finding a best fit corner frequency and moment and saving the results in `~/dataset/stress_drops/stress_drops_dataset.csv`.
10. Analyze results as desired.  

