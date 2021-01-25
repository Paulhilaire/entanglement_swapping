from tree import Tree
from error_model_bsm import BellStateMeasurementProtocol

class DynamicProtocol(BellStateMeasurementProtocol):
    """Dynamic protocol used for logical BSM.
    In this dynamic protocol, we perform two-photon BSM on the qubits
    whenever their parent qubits were measured with complete two-photon BSM.
    If that is not the case, we perform single qubit measurements to indirectly
    measure the Z components with better accuracy.

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
        super(DynamicProtocol, self).__init__(tree, pr_ph_1=pr_ph_1, pr_ph_2=pr_ph_2, error_1=error_1, error_2=error_2)
        self.name = "Dynamic"


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

        return 1 - (1 - self.proba["m_zz"][2] ** self.b[1]) ** m_c

    def error_m_xx_l(self):
        """
        Error of a XX measurement on the logical qubits.
        """
        return self.error["i_zz_c"][0]

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

    def calculate_specific_proba(self):
        """Calculate probas specifically used for a dynamic protocol."""

        # Start at k = len(self.tree.branches) (if b_0, ... b_{n-1}, start at k=n)
        # Increment -1 every time (so that you span the list from n to 0)
        # Finish for k = 0.
        for k in range(self.depth, -1, -1):
            # If k <= n-1 the qubit is not the lowest level qubit in the tree.
            if k <= self.depth - 1:
                # To measure a qubit with a single indirect measurement, you should measure its child qubit in x and its k+2 qubits in z.

                # If the BSM is complete, an indirect measurement is performed by measuring qubits via BSM.
                self.proba["s_zz_c"][k] = self.proba["d_xx"][k + 1] * \
                    self.proba["m_zz"][k + 2] ** self.b[k + 1]

                # To measure a qubit indirectly, at least one indirect measurement should succeed.
                self.proba["i_zz_c"][k] = 1 - \
                    (1 - self.proba["s_zz_c"][k]) ** self.b[k]

                # In that case, we do single qubit measurements and both should work:
                self.proba["i_zz_f"][k] = self.proba["i_z_1"][k] * \
                    self.proba["i_z_2"][k]
                #self.proba["i_zz_f"][k] =  1 - (1 - self.proba["s_zz_f"][k]) ** self.b[k]
                self.proba["i_zz_p"][k] = self.proba["i_z_1"][k] * \
                    self.proba["i_z_2"][k]
                # self.proba["i_zz_p"][k] = 1 - (1 - self.proba["s_zz_p"][k]) ** self.b[k]

            # To measure a qubit, you should measure it either directly or indirectly.
            # A failed BSM cannot measure ZZ directly:
            self.proba["m_zz_f"][k] = self.proba["i_zz_f"][k]
            # If the BSM was partial, ZZ is already measured:
            self.proba["m_zz_p"][k] = 1
            # If the BSM was complete, ZZ is already measured:
            self.proba["m_zz_c"][k] = 1

            # Overall probability:
            self.proba["m_zz"][k] = self.proba["d_xx"][k] * self.proba["m_zz_c"][k] \
                + (self.proba["d_zz"][k] - self.proba["d_xx"][k]) * self.proba["m_zz_p"][k] \
                + (1 - self.proba["d_zz"][k]) * self.proba["m_zz_f"][k]
            # The probability of ZZ measuring qubits at level k is the probability
            # to perform a complete BSM (pr_d_xx) times the probability of measuring this qubit in that case (which is one)
            # plus the probability to perform a partial BSM (pr_d_zz - pr_d_xx) times the probability of measuring this qubit in that case (which is one)
            # plus the probability to perform a failed BSM (1 - pr_d_zz) times the probability of measuring this qubit in that case.


    def calculate_specific_error(self):
        """Calculate errors specifically used for a dynamic protocol."""

        # Start at k = len(self.tree.branches) (if b_0, ... b_{n-1}, start at k=n)
        # Increment -1 every time (so that you span the list from n to 0)
        # Finish for k = 0.
        for k in range(self.depth, -1, -1):
            # If k <= n-1 the qubit is not the lowest level qubit in the tree.
            if k <= self.depth - 1:
                # To measure a qubit with a single indirect measurement, you should measure its child qubit in x and its k+2 qubits in z.
                # In a failed measurement all the qubits are individually measured
                self.error["i_zz_f"][k] = self.error["i_z_1"][k] + self.error["i_z_2"][k] - \
                    2 * self.error["i_z_1"][k] * self.error["i_z_2"][k]
            self.error["m_zz_f"][k] = self.error["i_zz_f"][k]
            if k == self.depth:
                self.error["m_zz_f"][k] = None
            # A failed BSM do not inform about ZZ, thus the error is the one given by the indirect measurement.

        for k in range(self.depth, -1, -1):
            # If k <= n-1 the qubit is not the lowest level qubit in the tree.
            if k <= self.depth - 1:
                # In a partial measurement all the child qubits are individually measured
                self.error["i_zz_p"][k] = self.error["i_z_1"][k] + self.error["i_z_2"][k] - \
                    2 * self.error["i_z_1"][k] * self.error["i_z_2"][k]
            self.error["m_zz_p"][k] = self.proba["i_zz_p"][k] / self.proba["m_zz_p"][k] * self.error["i_zz_p"][k] + \
                (1 - self.proba["i_zz_p"][k] /
                 self.proba["m_zz_p"][k]) * self.error["d_zz"][k]

        for k in range(self.depth, -1, -1):
            # If k <= n-1 the qubit is not the lowest level qubit in the tree.
            if k <= self.depth - 1:
                # Here the child qubits are measured in a BSM with three different outcomes...
                if k == self.depth - 2:
                    self.error["s_zz_c"][k] = self._err_s_zz_c_last(
                        self.proba["d_xx"][k + 2], self.proba["d_zz"][k + 2], self.b[k + 1], self.error["d_xx"][k + 1], self.error["m_zz_c"][k + 2], self.error["m_zz_p"][k + 2])
                else:
                    self.error["s_zz_c"][k] = self._err_s_zz_c(self.proba["d_xx"][k + 2], self.proba["d_zz"][k + 2], self.b[k + 1],
                                                              self.error["d_xx"][k + 1], self.error["m_zz_c"][k + 2], self.error["m_zz_p"][k + 2], self.error["m_zz_f"][k + 2])

                self.error["i_zz_c"][k] = self._err_i_zz(
                    self.proba["i_zz_c"][k], self.b[k], self.error["s_zz_c"][k], self.proba["s_zz_c"][k])
            self.error["m_zz_c"][k] = self.proba["i_zz_c"][k] / self.proba["m_zz_c"][k] * self.error["i_zz_c"][k] + \
                (1 - self.proba["i_zz_c"][k] /
                 self.proba["m_zz_c"][k]) * self.error["d_zz"][k]



    def _err_s_zz_c(self, pr_xx, pr_zz, n, err_xx_kp1, err_m_zz_c_kp2, err_m_zz_p_kp2, err_m_zz_f_kp2):
        """Calculate the error probability of a single indirect measurement when the two-photon BSM (on the directly measured qubits) was complete."""
        value = 0

        # Number of complete measurements on the child qubits
        for m_c in range(n + 1):
            # Number of partial measurements on the child qubits
            for m_p in range(n + 1 - m_c):
                value += self._binom_proba_2(pr_xx, pr_zz - pr_xx, m_c, m_p, n) * self._err_s_zz_c_given_mc_mp(
                    err_xx_kp1, err_m_zz_c_kp2, err_m_zz_p_kp2, err_m_zz_f_kp2, m_c, m_p, n)
        return value

    def _err_s_zz_c_given_mc_mp(self, err_xx_kp1, err_m_zz_c_kp2, err_m_zz_p_kp2, err_m_zz_f_kp2, m_c, m_p, n):
        """
        Single indirect measurement error given that m_c (m_p)  BSM on the child qubits were complete (partial).
        """

        value = 0
        for i in range(2):
            for j in range(m_c + 1):
                for k in range(m_p + 1):
                    for l in range(n + 1 - m_c - m_p):
                        if (i + j + k + l) % 2 == 0:
                            # Even number of errors cancel each other.
                            continue
                        else:
                            # i errors on the xx measurements.
                            b = err_xx_kp1 ** i * (1 - err_xx_kp1) ** (1 - i)
                            # j errors on the m_c complete BSMs.
                            b *= self._binom_proba(err_m_zz_c_kp2, m_c, j)
                            # k errors on the m_p partial BSMs.
                            b *= self._binom_proba(err_m_zz_p_kp2, m_p, k)
                            # l errors on the failed BSM (measured via indirect measurements).
                            b *= self._binom_proba(err_m_zz_f_kp2,
                                                  n - m_p - m_c, l)
                            value += b
        return value

    def _err_s_zz_c_last(self, pr_xx, pr_zz, n, err_xx_kp1, err_m_zz_c_kp2, err_m_zz_p_kp2):
        """Specific case for the deepest qubits of the tree that can be indirectly
        measured."""
        value = 0
        denominator = 0

        # Number of complete measurements on the child qubits
        for m_c in range(n + 1):
            # Number of partial measurements on the child qubits
            m_p = n - m_c
            b = self._binom_proba_2(pr_xx, pr_zz - pr_xx, m_c, m_p, n)
            denominator += b
            value += b * \
                self._err_s_zz_c_given_last(
                    err_xx_kp1, err_m_zz_c_kp2, err_m_zz_p_kp2, m_c, m_p, n)
        return value / denominator

    def _err_s_zz_c_given_last(self, err_xx_kp1, err_m_zz_c_kp2, err_m_zz_p_kp2, m_c, m_p, n):
        """
        Single indirect measurement error given that m_c (m_p)  BSM on the child qubits were complete (partial).
        For the last qubits that can be indirectly measured...
        """

        value = 0
        for i in range(2):
            for j in range(m_c + 1):
                for k in range(m_p + 1):
                    if (i + j + k) % 2 == 0:
                        # Even number of errors cancel each other.
                        continue
                    else:
                        # print(i, j, k, m_c, m_p, b)
                        # i errors on the xx measurements.
                        b = err_xx_kp1 ** i * (1 - err_xx_kp1) ** (1 - i)
                        # j errors on the m_c complete BSMs.
                        b *= self._binom_proba(err_m_zz_c_kp2, m_c, j)

                        # k errors on the m_p partial BSMs.
                        b *= self._binom_proba(err_m_zz_p_kp2, m_p, k)
                        value += b
        return value
