from enum import Enum
import numpy as np
import scipy.stats as stat
import math as math
import InputData as Data
import scr.MarkovClasses as MarkovCls
import scr.RandomVariantGenerators as Random
import scr.ProbDistParEst as Est


class HealthStats(Enum):
    """ health states of patients with HIV """
    Well = 0
    Stroke = 1
    Post_Stroke = 2
    Stroke_Dead= 3


class Therapies(Enum):
    """ mono vs. combination therapy """
    MONO = 0
    COMBO = 1


class ParametersFixed():
    def __init__(self, therapy):

        # selected therapy
        self._therapy = therapy

        # simulation time step
        self._delta_t = Data.DELTA_T

        # initial health state
        self._initialHealthState = HealthStats.Well

        # annual treatment cost
        if self._therapy == Therapies.MONO:
            self._annualTreatmentCost = Data.Zidovudine_COST
        else:
            self._annualTreatmentCost = Data.Zidovudine_COST + Data.Lamivudine_COST

        # transition probability matrix of the selected therapy
        self._prob_matrix = Data.TRANS_MATRIX
        # treatment relative risk
        self._treatmentRR = 0

        # calculate transition probabilities between hiv states
        self._prob_matrix = Data.TRANS_MATRIX

        # update the transition probability matrix if combination therapy is being used
        if self._therapy == Therapies.COMBO:
            # treatment relative risk
            self._treatmentRR = Data.TREATMENT_RR
            # calculate transition probability matrix for the combination therapy
            self._prob_matrix = Data.Treatment_Matrix
    def get_initial_health_state(self):
        return self._initialHealthState

    def get_delta_t(self):
        return self._delta_t

    def get_transition_prob(self, state):
        return self._prob_matrix[state.value]


def calculate_prob_matrix():
    """ :returns transition probability matrix for hiv states under mono therapy"""

    # create an empty matrix populated with zeroes
    prob_matrix = []
    for s in HealthStats:
        prob_matrix.append([0] * len(HealthStats))

    # for all health states
    for s in HealthStats:
        # if the current state is death
        if s == HealthStats.Stroke_Dead:
            # the probability of staying in this state is 1
            prob_matrix[s.value][s.value] = 1
        else:
            # calculate total counts of individuals
            sum_counts = sum(Data.TRANS_MATRIX[s.value])
            # calculate the transition probabilities out of this state
            for j in range(s.value, HealthStats.HIV_DEATH.value+1):
                prob_matrix[s.value][j] = Data.TRANS_MATRIX[s.value][j] / sum_counts

    return prob_matrix


def calculate_prob_matrix_combo(matrix_mono, combo_rr):
    """
    :param matrix_mono: (list of lists) transition probability matrix under mono therapy
    :param combo_rr: relative risk of the combination treatment
    :returns (list of lists) transition probability matrix under combination therapy """

    # create an empty list of lists
    matrix_combo = []
    for l in matrix_mono:
        matrix_combo.append([0] * len(l))

    # populate the combo matrix
    # first non-diagonal elements
    for s in HealthStats:
        for next_s in range(s.value + 1, len(HealthStats)):
            matrix_combo[s.value][next_s] = combo_rr * matrix_mono[s.value][next_s]

    # diagonal elements are calculated to make sure the sum of each row is 1
    for s in HealthStats:
        if s != HealthStats.HIV_DEATH:
            matrix_combo[s.value][s.value] = 1 - sum(matrix_combo[s.value][s.value + 1:])

    return matrix_combo
