from tree import Tree
from error_model_bsm import BellStateMeasurementProtocol



class StaticProtocol(BellStateMeasurementProtocol):
    """Static protocol used for logical BSM.
    In this static protocol, we perform two-photon BSM on all the qubits.
    This protocol can be near-deterministic, loss-tolerant and fault-tolerant.

    Parameters
    ----------
    tree : Tree()
        Tree graph state structure of the two logical qubits used for the
        calculations.
    pr_ph_1 : float (should be between 0 and 1)
        Single-photon detection probability of the first tree
        (assumed to be constant for all photons).
    pr_ph_2 : float (should be between 0 and 1)
        Single-photon detection probability of the first tree
        (assumed to be constant for all photons).
    error_1 : float (should be between 0 and 1)
        Single photon depolarization error probability (single qubit error rate)
    error_2 : float (should be between 0 and 1)
        Single photon depolarization error probability (single qubit error rate)

    Attributes
    ----------
    name : str
        Name of the protocol.

    """


    def __init__(self, tree=Tree(), pr_ph_1=1, pr_ph_2=1, error_1=0, error_2=0):
        super(StaticProtocol, self).__init__(tree, pr_ph_1=pr_ph_1, pr_ph_2=pr_ph_2, error_1=error_1, error_2=error_2)
        self.name = "Static"

    def proba_m_zz_l(self, m_c, m_p):
        """
        Probability of a ZZ measurement on the logical qubits
        given that there are m_c complete BSM and m_p partial BSM
        on the first level qubits.
        """
        return self.proba["i_zz"][1] ** (self.b[0] - m_c - m_p)

    def proba_m_xx_l(self, m_c):
        """
        Probability of a XX measurement on the logical qubits
        given that there are m_c complete BSM on the first level qubits.
        """
        return 1 - (1 - self.proba["m_zz"][2] ** self.b[1]) ** m_c

    def error_m_xx_l(self):
        """
        Error of a XX measurement on the logical qubits.
        """
        return self.error["i_zz"][0]

    def error_m_zz_l(self):
        """"
        Error of a ZZ measurement on the logical qubits.
        """
        value = 0
        for i in range(self.b[0] + 1):
            if i % 2 == 0:
                continue
                # Only odd number of parity errors leads to an error.
            value += self._binom_proba(self.error["m_zz"][1], self.b[0], i)
        return value

    def calculate_specific_proba(self):
        """Calculate probas specifically used for a static protocol."""

        # Rename because it is a simpler case.
        self.proba["m_zz"], self.proba["i_zz"], self.proba["s_zz"] = self.proba["m_zz_c"], self.proba["i_zz_c"], self.proba["s_zz_c"]
        self.proba["m_zz"], self.proba["i_zz"], self.proba["s_zz"], self.proba["d_zz"], self.proba["d_xx"] = self.construct_pr(
            self.proba["m_zz"], self.proba["i_zz"], self.proba["s_zz"], self.proba["d_zz"], self.proba["d_xx"])

    def calculate_specific_error(self):
        """Calculate errors specifically used for a static protocol."""
        # Rename because it is a simpler case.
        self.error["m_zz"], self.error["i_zz"], self.error["s_zz"] = self.error["m_zz_c"], self.error["i_zz_c"], self.error["s_zz_c"]
        self.error["m_zz"], self.error["i_zz"], self.error["s_zz"], self.error["d_zz"], self.error["d_xx"] = self.construct_error(
            self.proba["m_zz"], self.proba["i_zz"], self.proba["s_zz"], self.error["m_zz"], self.error["i_zz"], self.error["s_zz"], self.error["d_zz"], self.error["d_xx"])
