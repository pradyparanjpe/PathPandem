***READ THE COMPLETE FILE BEFORE INSTALLATION. AT LEAST READ THE *CAUTION* SECTION THOROUGHLY***

# Installation:
***INSTALLATION SHOULD HAPPEN LOCALLY FOR THE USER. DO NOT PROVIDE ROOT/ADMINISTRATIVE PRIVILEGES (OTHER THAN PYTHON INSTALLATION)***

## Windows 10
### Recommended:
Use `pip` method as described below.

### Else, if you really insist on a Windows EXE,
Download and extract git repository.
Navigate to [working directory](bin/PathPandem.Win10/)

Run file PathPandem.exe.

This executable file cannot be copied to any another location,
as its dependencies are packed in the same folder.

## Linux
### Recommended:
Use `pip` method as described below.

### Else,
Use generic method

## Generic: From the source-code
1. Download and unzip the git repository.
   - Using git.
```
git clone https://github.com/pradyparanjpe/PathPandem
```
2. Confirm pre-requisites.
3. In Command Line/Shell, navigate to the repository folder and run
```
python bin/PathPandem
```

### Pre-requisites for running from the source-code:
1. Python3.6 or higher
2. Numpy >= 1.18
3. Matplotlib >= 3.2.1
4. Gooey >= [1.0.3](https://github.com/chriskiehl/Gooey)

*1 can be installed from official source;
further, 2, 3, 4 can be installed by command `pip install <module>`*.

## pip
1. Install python3 from [official website](https://www.python.org/downloads/)
   - For Windows, enable "Add to PATH environment variable" during installation.
   - This will require Root/Administrative privileges.
2. Install PathPandem by typing in CommandPrompt/ Shell:

```
pip install --user PathPandem
```

3. Run by typing in CommandPrompt/ Shell:

```
python -m "PathPandem"

```

*On Unix-like systems, the file `bin/PathPandem` is executable.*

# Uninstallation

### If installed by `pip`
1. In CommandPrompt/Shell, type
```
pip uninstall PathPandem
```
2. If not needed, remove python [like any regular program]

### If installed by cloning git repo:
1. Delete repository folder
2. If not needed, remove python [like any regular program]

### If Windows 10 folder was downloaded
1. Delete folder


# Updates:

### If installed by `pip`
1. In CommandPrompt/Shell, type
```
pip install --upgrade PathPandem
```

### If installed by cloning git repo:
1. In CommandPrompt/Shell, navigate inside the repo folder
2. Type
```
git pull
```

### If Windows 10 folder was downloaded
1. Consider the option of deleting the folder and installation via `pip` or `git`
2. If you *still insist* on EXE, delete the previous folder and download updated folder from git repository.

*Back-up files before deleting the previous version. Recurrent downloads may consume a lot of Internet-Data.*


# Known Issues:
1. Python2's end-of-life was 2019-12-31.
   - If `python` means `python2` on an operating system, most commonly due to an earlier-installed `python2`, try migrating to `python3` (recommended), or if you understand what you are doing, install `python3` alongside and run all commands with explicit mention of python3: `python3`, `pip3`, etc.
2. wxPython currently fails installation via `pip`, [at least] on Linux.
   - Install it using a package-manager (dnf, pacman, apt[, brew, chocolate,] etc) or manually before `pip` tries to install Gooey`.

# Plot Legend:
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
- (*Any other standard RGB combinations*).

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
