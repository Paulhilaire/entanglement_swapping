from tree import Tree
from error_model_bsm import BellStateMeasurementProtocol


class LossOnlyProtocol(BellStateMeasurementProtocol):
    """Loss-only protocol used for logical BSM.
    In this loss-only protocol, we perform BSM on all the first level qubits of
    the tree. Then if a BSM is complete, we perform single qubit measurements Z
    on the second level child qubits. If the measurement is not complete,
    we perform indirect Z measurements via single qubit measurements (i.e. X
    measurements at the second level, then Z measurement at the third...).
    This protocol can be near-deterministic and loss-tolerant but it is not
    fault-tolerant.
    The reason why is that we do not measure indirectly the ZZ components with
    better accuracy when there is a successful complete measurement.


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
        super(LossOnlyProtocol, self).__init__(tree, pr_ph_1=pr_ph_1, pr_ph_2=pr_ph_2, error_1=error_1, error_2=error_2)
        self.name = "Loss-only"

    def proba_m_zz_l(self, m_c, m_p):
        """
        Probability of a ZZ measurement on the logical qubits
        given that there are m_c complete BSM and m_p partial BSM
        on the first level qubits.
        """
        return (self.proba["i_z_1"][1] * self.proba["i_z_2"][1]) ** (self.b[0] - m_c - m_p)

    def proba_m_xx_l(self, m_c):
        """
        Probability of a XX measurement on the logical qubits
        given that there are m_c complete BSM on the first level qubits.
        """
        return 1 - (1 - (self.proba["m_z_1"][2] * self.proba["m_z_2"][2]) ** self.b[1]) ** m_c


    def error_m_xx_l(self):
        """
        Error of a XX measurement on the logical qubits.
        """

        error = 0
        for i in range(2):
            for j in range(self.b[1] + 1):
                for k in range(self.b[1] + 1):
                    if (i + j + k) % 2 == 0:
                        continue
                    else:
                        b = self.error["d_xx"][1] ** i * (1 - self.error["d_xx"][1]) ** (1 - i)
                        b *= self._binom_proba(self.error["m_z_1"][2], self.b[1], j)
                        b *= self._binom_proba(self.error["m_z_2"][2], self.b[1], k)
                        error += b
        pr_s_zz = self.proba["d_xx"][1] * (self.proba["m_z_1"][2] * self.proba["m_z_2"][2]) ** self.b[1]
        value = 0
        denominator = 0
        for m_s in range(1, self.b[0] + 1):
            value += self._binom_proba(pr_s_zz, self.b[0], m_s) * self._majority_vote(error, m_s)
            denominator += self._binom_proba(pr_s_zz, self.b[0], m_s)
        return value / denominator


    def error_m_zz_l(self):
        """"
        Error of a ZZ measurement on the logical qubits.
        """
        value = 0
        denominator = 0
        for m_c in range(1, self.b[0] + 1):
            for m_p in range(self.b[0] + 1 - m_c):
                a = self._binom_proba_2(self.proba["d_xx"][0], self.proba["d_zz"][0] - self.proba["d_xx"][0], m_c, m_p, self.b[0])
                denominator += a
                value += a * self._e_logical_zz_given_bsm_results(m_c, m_p, self.b[0], self.error["m_zz_f"][1], self.error["m_zz_p"][1], self.error["m_zz_c"][1])
        return value / denominator


    def calculate_specific_error(self):
        """Calculate errors specifically used for the loss-only protocol."""
        self.error["m_zz_f"][1] = self.error["i_z_1"][1] + self.error["i_z_2"][1] - 2 * self.error["i_z_1"][1] * self.error["i_z_2"][1]
        self.error["i_zz_p"][1] = self.error["i_z_1"][1] + self.error["i_z_2"][1] - 2 * self.error["i_z_1"][1] * self.error["i_z_2"][1]
        self.error["m_zz_p"][1] = self.proba["i_z_1"][1] * self.proba["i_z_2"][1] * self.error["i_zz_p"][1] + (1 - self.proba["i_z_1"][1] * self.proba["i_z_2"][1]) * self.error["d_zz"][1]
        self.error["m_zz_c"][1] = self.error["d_zz"][1]
