# Logical Entanglement Swapping with tree graph states


Code to simulate a photonic Bell state measurement (BSM) on logical qubits encoded with tree graph states. It is a building block for an error-corrected / loss-tolerant quantum repeater protocol. This code accompanies the manuscript: https://arxiv.org/abs/2101.11082.




## Install
This model works on python 3.6 (or above) and only requires the numpy and scipy packages.

## Model
Here are the general details about the model.

### example.py:
Example script that shows how the model can be used. It inputs the necessary parameters to perform a calculations.

### tree.py
Structure of the tree graph state used for error correction.
It is completely defined by its branching vector, that denotes how many branches there are at each levels.

### error_model_bsm.py
General functions for the calculation of a logical Bell state measurement on a tree graph state. It calculates the success probability of a BSM (the probability that it yields an output) and the error probability (the probability that the output is correct). It is the parent object on which we can create different logical BMS protocols.

### Protocols

The other files corresponds to four different logical BSM protocols that have been developed. Two of them are described in the arXiv manuscript (static_protocol.py and dynamic_protocol.py), the last two are presented here but yields worse performances.

#### prototype_protocol.py
Protocol to perform a logical BSM which is not error-corrected nor loss-tolerant but still yields a near-deterministic success probability in the absence of photon loss and error.
(Measurement setting: all the first-level qubit are Bell-measured, the second-level qubits are measured in Z-basis.)

#### loss_only_protocol.py
Protocol that yields the best performances on loss-correction but does not enable error-correction.
(Measurement setting: all the first-level qubit are Bell-measured, the second-level (and beyond) qubits are measured individually with a basis that depends on the first-level BSM outcomes.)

#### static_protocol.py
All the qubits in the tree are Bell-measured. It is both error-corrected and loss-tolerant. More description in the manuscript.

#### dynamic_protocol.py
The measurement basis of the photons is decided based on the outcome of previous measurements to maximize both the loss-tolerance (almost as good as the loss-only protocol) and the error-correction.




## Resources

The tree graph states are described in Varnava et al., PRL (2006) (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.120501), see also Azuma et al., Nature Communications (2015) (https://www.nature.com/articles/ncomms7787?origin=ppub).

Its deterministic generation is described in Buterakos et al., Physical Review X (2017) (https://journals.aps.org/prx/abstract/10.1103/PhysRevX.7.041023).

The logical BSM are best described in the ArXiv paper (https://arxiv.org/abs/2101.11082) that this code accompanies, and mostly in its supplementary materials.
