

class Tree(object):
    """Tree cluster state.

    Example of tree
                     O
                   /   \ b_0 = 2
                  O     O
                 /|\   /|\  b_1 = 3
                O O O O O O
                depth=2
    It is determined by its depth n and its number of arm per level given by the branching parameters:
    branches=[b_0, ..., b_{n-1}]
    The total number of photon qubits in this tree is N_qubits
    The error correction code is given for an intrinsic error epsilon_0  (see Varnava et al. (2006)).
    The number of solid state qubits required for the realization of this tree is equal to N_ss of the tree (see Buterakos et al. (2017))
    """

    """Short summary.

    Parameters
    ----------
    branches : list
        branching vector of the tree.
    """

    def __init__(self, branches=[]):
        super(Tree, self).__init__()
        self._branches = branches

    @property
    def branches(self):
        """Branching parameter of the tree graph state."""
        return self._branches

    @branches.setter
    def branches(self, value):
        self._branches = value

    @property
    def depth(self):
        """Branching parameter of the tree graph state."""
        return len(self._branches)

    def add_branch(self, arm):
        """Add one level with "arm" number of arms at the bottom of the cluster tree."""
        self.depth += 1
        self.branches.append(arm)

    @property
    def N_qubits(self):
        """Number of photonic qubits"""
        N_qubits = 1
        N_prev = 1
        for i in range(self.depth):
            N_prev = N_prev * self.branches[i]
            N_qubits += N_prev
        return N_qubits



if __name__=="__main__":


    # Check some compatibilities with other programs
    from example_test import test
    test()
