import numpy as np
from scipy.special import binom

from tree import Tree


class BellStateMeasurementProtocol(object):
    """Logical Bell state measurement with tree graph states.
    This class is the bases of many BSM protocols and should be used via a given
    protocol child classes.
    It embeds all the necessary functions for BSM success probability and error
    calculations.

    dict_constr() initialize the dictionnaries self.proba and self.error on
    which all the calculations are made. calculate_proba() and calculate_error()
    are used to do calculations onto these dictionaries (Note that
    calculate_error() calls calculate_proba()). Other functions are used for
    these calculations.

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


    ----------
    self.proba contains the measurement probabilities.
    self.error contains the measurement errors.

    Notations:
    - 's' for a single indirect measurement.
    - 'i' for an indirect measurement.
    - 'd' for a direct measurement.
    - 'm' for a measurement (direct or indirect).

    - 'z_1' (resp 'z_2') for a z measurement on tree 1 (resp 2).
    - 'zz' for a double z measurement (on two trees).
    - same for 'x' and 'xx'

    """

    def __init__(self, tree=Tree(), pr_ph_1=1, pr_ph_2=1, error_1=0, error_2=0):
        self.tree = tree
        self.pr_ph_1 = pr_ph_1
        self.pr_ph_2 = pr_ph_2
        self.error_1 = error_1
        self.error_2 = error_2

    @property
    def pr_bell(self):
        """
        Probability of a complete two-photon Bell state measurement.
        """
        return self.pr_ph_1 * self.pr_ph_2 / 2

    @property
    def depth(self):
        """
        Depth of the tree.
        (used to have a shorter notation)
        """
        return self.tree.depth

    @property
    def b(self):
        """
        Branching vector.
        Adding two zeros terms to the branches to simplify calculations.
        """
        b = [x for x in self.tree.branches]
        b.append(0)
        b.append(0)
        return b

    @property
    def error_bell(self):
        """Error of a two-photon BSM measurement.
        """
        return 3 / 2 * (self.error_1 + self.error_2 - 2 * self.error_1 * self.error_2)

    def dict_constr(self):
        """ Initialize all the desired parameters in dictionaries.
        self.proba contains the measurement probabilities.
        self.error contains the measurement errors.
        Notations:
        - 's' for a single indirect measurement.
        - 'i' for an indirect measurement.
        - 'd' for a direct measurement.
        - 'm' for a measurement (direct or indirect).

        - 'z_1' (resp 'z_2') for a z measurement on tree 1 (resp 2).
        - 'zz' for a double z measurement (on two trees).
        - same for 'x' and 'xx'

        - 'c' for a complete BSM (ZZ and XX).
        - 'p' for a partial BSM (ZZ only).
        - 'f' for a failed BSM (nothing).

        Initialized value should be 0 for all except direct measurements
        """

        self.proba = {}
        self.error = {}
        list_type_of_measurements = ["m", "i", "s", "d"]
        list_basis_of_measurements = [
            "x_1", "x_2", "z_1", "z_2", "xx", "zz", "zz_c", "zz_p", "zz_f"]
        for i in list_type_of_measurements:
            for j in list_basis_of_measurements:
                key = "%s_%s" % (i, j)
                if i in ["m", "i", "s"]:
                    if j in ["x_1", "x_2", "xx"]:
                        # Do not create entries for indirect measurements for X measurement basis.
                        continue
                    elif i == "s" and j in ["zz_f", "zz_p"]:
                        # In that case, it is impossible to do indirect individual measurements since we are doing single qubit measurements instead
                        continue

                    # Initialize the matrix at zero.
                    else:
                        self.proba[key] = np.zeros(len(self.b),)
                        self.error[key] = np.zeros(len(self.b),)

                elif i == "d":
                    # Direct measurement is different and should be initialized differently
                    if j in ["x_1", "z_1"]:
                        # Measurement in the first qubit only
                        self.proba[key] = np.ones(len(self.b),) * self.pr_ph_1
                        self.error[key] = np.ones(len(self.b),) * self.error_1
                        # The last photon at level n+1 does not exist so cannot be measured. It can always "be measured" then.
                        self.proba[key][-1] = 1
                        self.error[key][-1] = 0
                    elif j in ["x_2", "z_2"]:
                        # Measurement in the second qubit only
                        self.proba[key] = np.ones(len(self.b),) * self.pr_ph_2
                        self.error[key] = np.ones(len(self.b),) * self.error_2
                        # The last photon at level n+1 does not exist so cannot be measured. It can always "be measured" then.
                        self.proba[key][-1] = 1
                        self.error[key][-1] = 0
                    elif j == "xx":
                        # Measurement of the two qubits in the XX basis
                        self.proba[key] = np.ones(len(self.b),) * self.pr_bell
                        self.error[key] = np.ones(
                            len(self.b),) * self.error_bell
                        # The last photon at level n+1 does not exist so cannot be measured TODO Fix that to make it more natural.
                        self.proba[key][-1] = 1
                        self.error[key][-1] = 0
                    elif j == "zz":
                        # Measurement of the two qubits in the ZZ basis (partial BSM works)
                        self.proba[key] = np.ones(
                            len(self.b),) * 2 * self.pr_bell
                        self.error[key] = np.ones(
                            len(self.b),) * (self.error_1 + self.error_2 - 2 * self.error_1 * self.error_2)

                        # The last photon at level n+1 does not exist so cannot be measured TODO Fix that to make it more natural.
                        self.proba[key][-1] = 1
                        self.error[key][-1] = 0

    def pr_logical_BSM(self):
        """
        Calculate the probability of a complete logical BSM.
        """
        value = 0
        for m_c in range(1, 1 + self.b[0]):
            for m_p in range(1 + self.b[0] - m_c):
                value += self.proba_bsm(m_c, m_p, self.b[0]) * self.proba_m_zz_l(
                    m_c, m_p) * self.proba_m_xx_l(m_c)
        self.proba["bsm_logical"] = value
        return value

    def error_logical_BSM(self):
        """
        Calculate the error probability of a logical BSM.
        """
        self.error["zz_logical"] = self.error_m_zz_l()
        self.error["xx_logical"] = self.error_m_xx_l()
        # print(self.error_m_xx_l())
        self.error["bsm_logical"] = self.error["zz_logical"] + (1 - self.error["zz_logical"]) * self.error["xx_logical"]
        return self.error["bsm_logical"]

    def proba_bsm(self, m_c, m_p, n, k=0):
        """
        Probability of having "m_c" complete BSM and "m_p" partial BSM
        given "n" trials.
        """
        return self._binom_proba_2(self.proba["d_xx"][k], self.proba["d_zz"][k] - self.proba["d_xx"][k], m_c, m_p, n)

    def construct_pr(self, pr_m, pr_i, pr_s, pr_d, pr_x):
        """Standard calculation of the probabilities.
        Works for single photon probabilities and transverse.
        pr_m is the measurement probability
        pr_i is the total indirect measurement probability
        pr_s is the single indirect measurement probability
        pr_d is the direct measurement probability in z basis
        pr_x is the direct measurement probability in x basis
        """

        # Start at k = len(self.tree.branches) (if b_0, ... b_{n-1}, start at k=n)
        # Increment -1 every time (so that you span the list from n to 0)
        # Finish for k = 0.
        for k in range(self.depth, -1, -1):
            # If k <= n-1 the qubit is not the lowest level qubit in the tree.
            if k <= self.depth - 1:
                # To measure a qubit with a single indirect measurement, you should measure its child qubit in x and its k+2 qubits in z.
                pr_s[k] = pr_x[k + 1] * pr_m[k + 2] ** self.b[k + 1]
                # To measure a qubit indirectly, at least one indirect measurement should succeed.
                pr_i[k] = 1 - (1 - pr_s[k]) ** self.b[k]
            # To measure a qubit, you should measure it either directly or indirectly.
            pr_m[k] = pr_d[k] + (1 - pr_d[k]) * pr_i[k]
        return pr_m, pr_i, pr_s, pr_d, pr_x

    def calculate_proba(self):
        """Calculate all the probabilities."""
        # First initialize the functions
        self.dict_constr()

        # Single qubit measurements:
        self.proba["m_z_1"], self.proba["i_z_1"], self.proba["s_z_1"], self.proba["d_z_1"], self.proba["d_x_1"] = self.construct_pr(
            self.proba["m_z_1"], self.proba["i_z_1"], self.proba["s_z_1"], self.proba["d_z_1"], self.proba["d_x_1"])
        self.proba["m_z_2"], self.proba["i_z_2"], self.proba["s_z_2"], self.proba["d_z_2"], self.proba["d_x_2"] = self.construct_pr(
            self.proba["m_z_2"], self.proba["i_z_2"], self.proba["s_z_2"], self.proba["d_z_2"], self.proba["d_x_2"])

        # For two qubits it is a bit more complicated since the measurement pattern can depend on the protocol.
        # If necessary we calculate these probas in the child class by overwriting the following function.
        self.calculate_specific_proba()

        # Calculate the logical probability of BSM.
        self.pr_logical_BSM()


    def construct_error(self, pr_m, pr_i, pr_s, err_m, err_i, err_s, err_d, err_x):
        """Standard calculation of the errors.
        Works for single photon probabilities and transverse.
        pr_m is the measurement probability
        pr_i is the total indirect measurement probability
        pr_s is the single indirect measurement probability
        pr_d is the direct measurement probability in z basis
        pr_x is the direct measurement probability in x basis
        and similarly for the errors.
        """

        # Start at k = len(self.b) - 1 = len(self.tree.branches) + 1 (if b_0, ... b_{n-1}, start at k=n)
        # Increment -1 every time (so that you span the list from n to 0)
        # Finish for k = 0.
        for k in range(self.depth, -1, -1):
            # If k <= n-1 the qubit is not the lowest level qubit in the tree.
            if k <= self.depth - 1:
                # To measure a qubit with a single indirect measurement, you should measure its child qubit in x and its k+2 qubits in z.
                # An error occurs if there is an odd number of error in these measurements.
                err_s[k] = self._err_s_measurement(
                    err_m[k + 2], err_x[k + 1], self.b[k + 1])

                # To measure a qubit indirectly, we use the majority vote.
                err_i[k] = self._err_i_zz(
                    pr_i[k], self.b[k], err_s[k], pr_s[k])

            # Favor indirect measurement error over direct measurement errors.
            err_m[k] = pr_i[k] / pr_m[k] * err_i[k] + \
                (1 - pr_i[k] / pr_m[k]) * err_d[k]
        return err_m, err_i, err_s, err_d, err_x



    def calculate_error(self):
        """Calculate all the errors."""

        # First initialize the functions
        self.dict_constr()
        self.calculate_proba()
        # Single qubit measurements:
        self.error["m_z_1"], self.error["i_z_1"], self.error["s_z_1"], self.error["d_z_1"], self.error["d_x_z_1"] = self.construct_error(
            self.proba["m_z_1"], self.proba["i_z_1"], self.proba["s_z_1"], self.error["m_z_1"], self.error["i_z_1"], self.error["s_z_1"], self.error["d_z_1"], self.error["d_x_1"])
        self.error["m_z_2"], self.error["i_z_2"], self.error["s_z_2"], self.error["d_z_2"], self.error["d_x_z_2"] = self.construct_error(
            self.proba["m_z_2"], self.proba["i_z_2"], self.proba["s_z_2"], self.error["m_z_2"], self.error["i_z_2"], self.error["s_z_2"], self.error["d_z_2"], self.error["d_x_2"])

        # For two qubits it is a bit more complicated since the measurement pattern can dependon the protocol.

        # For two qubits it is a bit more complicated since the measurement pattern can depend on the protocol.
        # If necessary we calculate these errors in the child class by overwriting the following function.
        self.calculate_specific_error()

        # Calculate the logical error of BSM.
        self.error_logical_BSM()



    def _err_s_measurement(self, err_m_kp2, err_x, n):
        """Single indirect error measurement."""
        value = 0
        for i in range(2):
            for j in range(n + 1):
                if (i + j) % 2 == 0:
                    continue
                else:
                    b = err_x ** i * (1 - err_x) ** (1 - i)
                    value += self._binom_proba(err_m_kp2, n, j) * b
        return value

    def _err_i_zz(self, pr_i_zz, n, err_s_zz, pr_s_zz):
        """
        Error probability of an indirect  measurement via n trials.
        """
        value = 0
        for m_s in range(1, n + 1):
            a = self._binom_proba(pr_s_zz, n, m_s) * \
                self._majority_vote(err_s_zz, m_s)
            value += a
        return value / pr_i_zz

    def _majority_vote(self, error, m_s):
        """
        Majority vote that reduces the error probability
        via m_s indirect measurements.
        """
        if m_s % 2 == 0:
            # if m_s is odd, remove one measurement randomly
            m = m_s - 1
        else:
            m = m_s
        value = 0
        m2 = int(m / 2)
        for j in range(m2 + 1, m + 1):
            # An error on the global indirect measurement remains
            # if there are more faulty than  correct measurements
            value += self._binom_proba(error, m, j)
        return value


    def _binom_proba(self, proba, n, m):
        """
        Binomial probability of having "m" success out of "n" trials with
        success probability "proba".
        """
        return binom(n, m) * proba ** m * (1 - proba) ** (n - m)

    def _binom_proba_2(self, proba_1, proba_2, m_1, m_2, n):
        """
        Double outcome binomial probability of having "m_1" outcome 1
        and "m_2" outcome 2 out of "n" trials with respective probability "proba_1" and "proba_2".
        """
        return binom(n, m_1) * proba_1 ** m_1 * binom(n - m_1, m_2) * proba_2 ** m_2 * (1 - proba_1 - proba_2) ** (n - m_1 - m_2)


    def _e_logical_zz_given_bsm_results(self, m_c, m_p, n, e_m_zz_f_1, e_m_zz_p_1, e_m_zz_c_1):
        """
        Calculate the logical error ZZ probability for specific first level BSM results m_c and m_p.
        Shared by Protocol1 and Protocol3"""
        value = 0
        for j in range(m_c + 1):
            for k in range(m_p + 1):
                for l in range(n+1 - m_c - m_p):
                    if (j + k + l) % 2 == 0:
                        continue
                    else:
                        a = 1
                        a *= self._binom_proba(e_m_zz_c_1, m_c, j)
                        a *= self._binom_proba(e_m_zz_p_1, m_p, k)
                        a *= self._binom_proba(e_m_zz_f_1, n - m_p - m_c, l)
                        value += a
        return value


    # Function to overwrite in child classes.
    def error_m_zz_l(self):
        """Error probability of a logical ZZ measurement."""
        raise NameError("error_m_zz_l should be erased in parent class")

    def error_m_xx_l(self):
        """Error probability of a logical XX measurement."""
        raise NameError("error_m_xx_l should be erased in parent class")

    def proba_m_zz_l(self, m_c, m_p):
        """Probability of a logical ZZ measurement."""
        raise NameError("proba_m_zz_l should be erased in parent class")

    def proba_m_xx_l(self, m_c):
        """Probability of a logical XX measurement."""
        raise NameError("proba_m_xx_l should be erased in parent class")

    def calculate_specific_proba(self):
        """Calculate specific probabilities depending on the protocol used.
        Should be overwritten in the child class."""
        pass

    def calculate_specific_error(self):
        """Calculate specific errors depending on the protocol used.
        Should be overwritten in the child class."""
        pass
