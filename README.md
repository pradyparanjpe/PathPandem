# Installation:
## Windows 10
### Recommended:
Use `pip` method as described below.

### Else, if you really insist on a Windows exe,
Download [working directory](bin/PathPandem.Win10/)

Run file [PathPandem.exe](bin/PathPandem.Win10/PathPandem.exe).

In Windows, [this exe](bin/PathPandem.Win10/PathPandem.exe) file
cannot be copied to any another location,
as its dependencies are packed in the same folder.

## Linux
### Recommended:
Use `pip` method as described below.

### Else,
Use generic method

## Generic: From the source-code
1. Download and unzip the git repository.
2. Confirm pre-requisites.
3. In Command Line/Shell, navigate to the repository folder and run `$ Python3 bin/PathPandem`

### Pre-requisites for running from the source-code:
1. Python3.8 or higher
2. Numpy >= 1.18
3. Matplotlib >= 3.2.1
4. Gooey >= [1.0.3](https://github.com/chriskiehl/Gooey)

*1 may be installed from official source; further, 2, 3, 4 may be installed by command `$ pip install <module>`*.

## pip
$ pip install [--user] PathPandem
$ PathPandem

*Square brackets indicate "optional".*

# Legend:
## Background Colour:
### Movements
- Green: No restrictions on movement.
- Red: Lockdown Imposed.

### Medical Progress
- Blue: Drug discovered.
- Cyan: Vaccine discovered.

### Combinations
- Grey: Red + Cyan.
- Magenta: Red + Blue.
- (Any other standard RGB combinations).

# Caution:
1. Population more than 10000 may stall the system.
2. Tested only on Linux running from source-code.
3. *True* numbers are plotted. However in reality, infection manifests symptoms after an initial lag of 1-3 days and test results appear further later by 1-2 days. Hence, graph trends need be imagined as having shifted suitably.
4. Although Infection may appear to exhaust in small sized, limited population; in reality, due to birth of new individuals, and in a very large population, the pathogen persists around at extemely low density.

# Composition of scenario:
- The GUI only edits the blanket population behaviour.
- A heterogenous population can be composed using basic Python scripting in the `spread_simul.py` to construct heterogenously behaving population.

# TODO:
- Replace unimodal movement of people around their home to bimodal movement between home and workplace.
- Parallelize numpy matrix `ufuncs` if possible.
- Include asymptomatic patients/carriers. Limit movement of serious cases [although this won't have a visible effect for diseases with majority of cases being mild].
- Animation, saved as mp4 for review

# Brief epidemiological explanation:
- Herd immunity starts reducing viral presence in community after viral steady state. i.e. plot of *Active* patients flattens. This happens when [1 - (1/R_{0})] fraction of the community becomes resistant. (Through vaccination or exposure)
- Medicine development is fairly a rare event given the rightful stringency involved in testing.
- With small population size, random fluctuations become impactful. Multiple runs with same parameters are recommended.
- Visualization is recommended only with very small population size.
