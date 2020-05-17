# Installation:

### Windows
Download [working directory](exec/PathPandem.Win10)
It contains files `reverse_cfr_database.pkl` and `spread_simul.exe`
Both files must be present in the same folder all the time.

Run file `spread_simul.exe`.

### Linux / MacOS
Run from [source-code](src/spread_simul.py)

## Running from the source-code
1. Confirm prerequisites
2. In Command Line/Shell, run `python3 spread_simul.py`

### Pre-requisites for running from the source-code:
1. Python3.8 or higher
2. Numpy >= 1.18
3. Matplotlib >= 3.2.1
4. Gooey >= [1.0.3](https://github.com/chriskiehl/Gooey)

*2 thourgh 4 may be downloaded by using `pip install <module>`*

## Legend:
### Background Colour:**
**Movements**
- Green: No restrictions on movement
- Red: Lockdown Imposed

**Scientific Progress**
- Blue: Drug discovered
- Cyan: Vaccine discovered

**Combinations**
- Grey: Red + Cyan
- Magenta: Red + Blue
- (Any other standard RGB combinations)

## Disclaimer:
1. Population more than 10000 may stall the system
2. Tested only on Linux running from source-code

### Composition of scenario
- The GUI only edits the blanket population behaviour.
- A population can be composed using basic Python scripting in the `if __name__ == "__main__":` section to construct heterogenously behaving population.

### TODO:
- Plot a representation of public movements. (This might be really heavy)
- Replace Unimodal movement of people around their home to bimodal movement between home and workplace.
- Parallelize numpy matrix `ufunc`s if possible.
