
import numpy as np

from tree import Tree
from prototype_protocol import PrototypeProtocol
from loss_only_protocol import LossOnlyProtocol
from static_protocol import StaticProtocol
from dynamic_protocol import DynamicProtocol



def read_param(param):
    """Read the parameters and calculate the protocol errors and probabilities."""
    tree = Tree()
    tree.branches = param["branches"]
    protocol_list = [PrototypeProtocol(), LossOnlyProtocol(), StaticProtocol(), DynamicProtocol()]
    protocol = protocol_list[param["protocol"]]
    protocol.tree = tree
    protocol.pr_ph_1 = param["pr_ph"]
    protocol.pr_ph_2 = param["pr_ph"]
    protocol.error_1 = param["error"]
    protocol.error_2 = param["error"]
    return protocol



if __name__ == "__main__":

    param = {
        "branches": [15, 15, 2],  # branching parameter of the RGS
        "error": 10 ** -5,  # single photon qubit error rate
        "pr_ph": 0.95,  # detection efficiency of single photons
        "protocol": 0, # protocol number
    }
    protocol = read_param(param)

    print("tree branches: %s" % (protocol.tree.branches))
    print("n_qubits: %s" % (protocol.tree.N_qubits))
    print()

    for i in range(4):
        # Protocol number is: 2 for static / 3 for dynamic (0 and 1 for prototype and loss only protocols)
        param["protocol"] = i
        protocol = read_param(param)
        protocol.calculate_error()
        print(i, protocol.name)
        print(protocol.proba["bsm_logical"])
        print(protocol.error["bsm_logical"])
        print()
