# RNAPII_kinetics_simulation
Python script to simulate global transcriptional kinetics from a model for RNA polymerase II occupancy.

# Code for modeling simulations	

There are two files, RNAPII_occupancy.py and model_fitting.py. RNAPII_occupancy.py contains the code for the STM and TFO models and upon running plots an example of each. model_fitting.py performs the grid search and fits the models to the ChEC-seq2 data. As a warning, model_fitting.py is implemented serially here and takes several hours to run. 


## Usage


```python
from RNAPII_occupancy import STM, TFO
from model_fitting import Fit

# example STM simulation for given parameters k_4, k_3r (k_-3 in the text), k_2, and k_2r (k_-2 in the text)
STM([k_4, k_3r, k_2, k_2r])

# example TFO simulation for given parameters k_4, and k_3r (k_-3 in the text)
TFO([k_4, k_3r, k_2, k_2r])

# Fit STM model
Fit(STM)

# Fit TFO model
Fit(TFO)
```
