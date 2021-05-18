import numpy as np
import pandas as pd
from pandas import DataFrame

import mea_poa.sub_matrix as sub_matrix
import math
from collections import defaultdict
print_this = False

pd.set_option('display.max_columns', None)


class PairHMM():
    def __init__(self, profiles, tau, epsilon, delta, emissionX, emissionY, subsmat, log_transform):
        self.tau = tau
        self.epsilon = epsilon
        self.delta = delta
        self.profiles = profiles
        self.profile1 = profiles[0]
        self.profile2 = profiles[1]
        self.m = len(self.profile1.profile)
        self.n = len(self.profile2.profile)
        self.aligned_seqs = []
        self.viterbi_aligned_positions = []
        self.mea_aligned_positions = []
        self.subsmat = subsmat
        self.log_transform = log_transform
        self.emissionX = emissionX
        self.emissionY = emissionY
        self.viterbi_matches1 = []
        self.viterbi_matches2 = []
        self.mea_matches1 = []
        self.mea_matches2 = []
        self.full_prob = 0
        self.vM = None
        self.vX = None
        self.vY = None
        self.fM = None
        self.fX = None
        self.fY = None
        self.bM = None
        self.bX = None
        self.bY = None
        self.pM = None
        self.initialise_probs()
        self.predecessors1 = profiles[0].predecessors if profiles[0].predecessors else None
        self.predecessors2 = profiles[1].predecessors if profiles[1].predecessors else None

        # print ('here is the profile')
        # print (self.profile1)
        # print (self.profile2)
        #
        # print ('here are the predecessors')
        #
        # print (self.predecessors1)
        # print (self.predecessors2)

    # Initialise the transition probabilities
    def initialise_probs(self):
        self.transitionMM = 1 - (2 * self.delta) - self.tau;
        self.transitionMX = self.transitionMY = self.delta;
        self.transitionXX = self.transitionYY = self.epsilon;
        self.transitionXM = self.transitionYM = 1 - self.epsilon - self.tau;

        self.emissionX = self.epsilon
        self.emissionY = self.epsilon


    # Initialise the matrices
    def create_viterbi_matrices(self):

        self.vM = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self.vY = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self.vX = np.zeros((self.m + 1, self.n + 1), dtype=float)

        self.tracebackM = np.zeros((self.m + 1, self.n + 1), dtype=object)
        self.tracebackX = np.zeros((self.m + 1, self.n + 1), dtype=object)
        self.tracebackY = np.zeros((self.m + 1, self.n + 1), dtype=object)

        # Matrices for holding the transition probabilities for gap opening and extension
        self.openX = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self.openY = np.zeros((self.m + 1, self.n + 1), dtype=float)

        # Add in the first opening penalty
        self.openX[0][0] = np.log(self.transitionMX)
        self.openY[0][0] = np.log(self.transitionMY)

        self.extensionX = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self.extensionY = np.zeros((self.m + 1, self.n + 1), dtype=float)

        if self.log_transform:


            for i in range(0, self.m + 1):
                self.vX[i][0] = float('-Inf')
                self.vM[i][0] = float('-Inf')

                # Add in the opening penalty
                self.openX[i][0] = np.log(self.transitionMX)

            for i in range(0, self.n + 1):
                self.vY[0][i] = float('-Inf')
                self.vM[0][i] = float('-Inf')

                # Add in the opening penalty
                self.openY[0][i] = np.log(self.transitionMY)

            self.vM[0][0] = np.log(math.e)
            self.vX[0][0] = np.log(math.e)
            self.vY[0][0] = np.log(math.e)

        else:
            self.vM[0][0] = 1
            for i in range(1, self.m + 1):
                self.vX[i][0] = 0
                self.vM[i][0] = 0

            for i in range(1, self.n + 1):
                self.vY[0][i] = 0
                self.vM[0][i] = 0


    def create_fb_matrices(self):

        self.fM = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self.fX = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self.fY = np.zeros((self.m + 1, self.n + 1), dtype=float)

        self.bM = np.zeros((self.m + 2, self.n + 2), dtype=float)
        self.bX = np.zeros((self.m + 2, self.n + 2), dtype=float)
        self.bY = np.zeros((self.m + 2, self.n + 2), dtype=float)

        if self.log_transform:

            for i in range(0, self.m + 1):
                self.fX[i][0] = float('-Inf')
                self.fM[i][0] = float('-Inf')

            for i in range(0, self.n + 1):
                self.fY[0][i] = float('-Inf')
                self.fM[0][i] = float('-Inf')

            for i in range(0, self.m + 1):
                self.bX[i][self.n + 1] = float('-Inf')
                self.bY[i][self.n + 1] = float('-Inf')

                self.bM[i][self.n + 1] = float('-Inf')

            for i in range(0, self.n + 1):
                self.bY[self.m + 1][i] = float('-Inf')
                self.bX[self.m + 1][i] = float('-Inf')

                self.bM[self.m + 1][i] = float('-Inf')

            self.fM[0][0] = np.log(math.e)
            self.fX[0][0] = np.log(math.e)
            self.fY[0][0] = np.log(math.e)


            self.bM[self.m][self.n] = np.log(self.tau)
            self.bX[self.m][self.n] = np.log(self.tau)
            self.bY[self.m][self.n] = np.log(self.tau)



        else:

            self.fM[0][0] = 1
            self.fX[0][0] = 1
            self.fY[0][0] = 1


            for i in range(1, self.m + 1):
                self.fX[i][0] = 0
                self.fM[i][0] = 0

            for i in range(1, self.n + 1):
                self.fY[0][i] = 0
                self.fM[0][i] = 0

            for i in range(0, self.m + 1):
                self.bX[i][self.n + 1] = 0
                self.bM[i][self.n + 1] = 0
            for i in range(0, self.n + 1):
                self.bY[self.m + 1][i] = 0
                self.bM[self.m + 1][i] = 0

            self.bM[self.m][self.n] = self.tau
            self.bX[self.m][self.n] = self.tau
            self.bY[self.m][self.n] = self.tau

    def performViterbiAlignment(self, po=False):

        self.create_viterbi_matrices()

        for i in range(self.m + 1):
            for j in range(self.n + 1):
                if not (i == 0 and j == 0):


                    if po:
                        self.fillVM_PO(i, j)
                        self.fillVX_PO(i, j)
                        self.fillVY_PO(i, j)

                    else:
                        self.fillVM(i, j)
                        self.fillVX(i, j)
                        self.fillVY(i, j)


        if po:
            aligned_profile = self.perform_po_viterbi_traceback()
        else:
            aligned_profile = self.perform_viterbi_traceback()


        return aligned_profile

    def fillVM(self, i, j):
        if (i == 0 or j == 0):
            return
        if self.log_transform:

            if print_this:

                print ("CONSIDERING " + str(i) + " and  " + str(j))

            emissionM = self.get_emission(i, j)

            emissionM = np.log(emissionM)

            # print (self.transitionXM)

            curr_transitionMM = float('-Inf') if (self.vM[i - 1][j - 1] == float('-Inf') or self.vM[i - 1][j - 1] ==
                                                  float('-Inf'))  else np.log(
            self.transitionMM) + \
                                                                                  self.vM[i - 1][j - 1]
            curr_transitionXM = float('-Inf') if (self.vX[i - 1][j - 1] == float('-Inf') or self.vX[i - 1][j - 1] ==
                                                  float('-Inf'))  else np.log(self.transitionXM) + \
                                                                                  self.vX[i - 1][j - 1]
            curr_transitionYM = float('-Inf') if (self.vY[i - 1][j - 1] == float('-Inf') or self.vY[i - 1][j - 1] ==
                                                  float('-Inf')) else np.log(self.transitionYM) + \
                                                                                  self.vY[i - 1][j - 1]

            if curr_transitionMM >= curr_transitionXM and curr_transitionMM >= curr_transitionYM:
                self.vM[i][j] = curr_transitionMM + emissionM
                self.tracebackM[i][j] = "M"

            elif curr_transitionXM >= curr_transitionMM and curr_transitionXM >= curr_transitionYM:
                self.vM[i][j] = curr_transitionXM + emissionM
                self.tracebackM[i][j] = "X"
            else:
                self.vM[i][j] = curr_transitionYM + emissionM
                self.tracebackM[i][j] = "Y"

        else:
            emissionM = self.get_emission(i, j)

            curr_transitionMM = 0 if self.vM[i - 1][j - 1] == -1 else self.transitionMM * self.vM[i - 1][
                j - 1]
            curr_transitionXM = 0 if self.vX[i - 1][j - 1] == -1 else self.transitionXM * self.vX[i - 1][
                j - 1]
            curr_transitionYM = 0 if self.vY[i - 1][j - 1] == -1 else self.transitionYM * self.vY[i - 1][
                j - 1]

            if curr_transitionMM >= curr_transitionXM and curr_transitionMM >= curr_transitionYM:
                self.vM[i][j] = curr_transitionMM * emissionM
                self.tracebackM[i][j] = "M"
            elif curr_transitionXM >= curr_transitionMM and curr_transitionXM >= curr_transitionYM:
                self.vM[i][j] = curr_transitionXM * emissionM
                self.tracebackM[i][j] = "X"
            else:
                self.vM[i][j] = curr_transitionYM * emissionM
                self.tracebackM[i][j] = "Y"

    def fillVX(self, i, j):
        if j == 0:
            return

        if self.log_transform:
            curr_transitionXX = float('-Inf') if (self.vX[i][j - 1] == float('-Inf') or self.vX[i][j-1] == 0)   else \
                np.log(
                self.transitionXX) + self.vX[i][
                j - 1]
            curr_transitionMX = float('-Inf') if (self.vM[i][j - 1] == float('-Inf') or self.vM[i][j-1] == 0)  else \
                np.log(self.transitionMX) + self.vM[i][
                j - 1]

            if curr_transitionXX >= curr_transitionMX:
                self.vX[i][j] = curr_transitionXX + np.log(self.emissionX)
                self.extensionX[i][j] = curr_transitionXX
                self.tracebackX[i][j] = "X"

            else:
                self.vX[i][j] = curr_transitionMX + np.log(self.emissionX)
                self.openX[i][j] = curr_transitionMX

                self.tracebackX[i][j] = "M"
        else:
            curr_transitionXX = 0 if self.vX[i][j - 1] == -1 else self.transitionXX * self.vX[i][j - 1]
            curr_transitionMX = 0 if self.vM[i][j - 1] == -1 else self.transitionMX * self.vM[i][j - 1]

            if curr_transitionXX >= curr_transitionMX:
                self.vX[i][j] = curr_transitionXX * self.emissionX
                self.tracebackX[i][j] = "X"

            else:
                self.vX[i][j] = curr_transitionMX * self.emissionX
                self.tracebackX[i][j] = "M"

    def fillVY(self, i, j):
        if i == 0:
            return

        if self.log_transform:

            curr_transitionYY = float('-Inf') if (self.vY[i - 1][j] == float('-Inf') or self.vY[i-1][j] == 0)  else \
                np.log(
                self.transitionYY) + \
                                                                              self.vY[i - 1][j]
            curr_transitionMY = float('-Inf') if (self.vM[i - 1][j] == float('-Inf') or self.vM[i-1][j] == 0)   else \
                np.log(self.transitionMY) + \
                                                                              self.vM[i - 1][j]

            if curr_transitionYY >= curr_transitionMY:
                self.vY[i][j] = curr_transitionYY + np.log(self.emissionY)
                self.extensionY[i][j] = curr_transitionYY

                self.tracebackY[i][j] = "Y"

            else:
                self.vY[i][j] = curr_transitionMY + np.log(self.emissionY)
                self.openY[i][j] = curr_transitionMY

                self.tracebackY[i][j] = "M"
        else:
            curr_transitionYY = 0 if self.vY[i - 1][j] == -1 else self.transitionYY * self.vY[i - 1][j]
            curr_transitionMY = 0 if self.vM[i - 1][j] == -1 else self.transitionMY * self.vM[i - 1][j]

            if curr_transitionYY >= curr_transitionMY:
                self.vY[i][j] = curr_transitionYY * self.emissionY
                self.tracebackY[i][j] = "Y"

            else:
                self.vY[i][j] = curr_transitionMY * self.emissionY
                self.tracebackY[i][j] = "M"

    def fillVM_PO(self, i, j):
        if (i == 0 or j == 0):
            return

        # print ('vm po')
        #
        # print (self.vM)

        # if print_this:


            # print("CONSIDERING " + str(i) + " and  " + str(j))
            #
            # print (self.predecessors1[i])
            # print (self.predecessors2[j])
            #
            # print ('Check matches')
        best_M_val, best_M_pos = self.get_best_viterbi_candidate(self.vM, self.predecessors1[i], self.predecessors2[j])

        # if print_this:

            # print ('Check x')
        best_X_val, best_X_pos = self.get_best_viterbi_candidate(self.vX, self.predecessors1[i], self.predecessors2[j])
        # if print_this:

            # print ('Check y')
        best_Y_val, best_Y_pos = self.get_best_viterbi_candidate(self.vY, self.predecessors1[i], self.predecessors2[j])

        # if print_this:
        #
        #     print ('best m pos')
        #     print (best_M_pos)
        #
        #     print ('best m val')
        #     print (best_M_val)


        # best_M_val, best_M_pos = self.vM[i-1][j-1], (i-1, j-1)
        # best_X_val, best_X_pos = self.vX[i-1][j-1], (i-1, j-1)
        # best_Y_val, best_Y_pos = self.vY[i-1][j-1], (i-1, j-1)
        #
        # # Multiple predecessors in both graphs
        # if len(self.predecessors1[i] > 1 and len(self.predecessors2[j]) > 1):
        #     best_M_val, best_M_pos = self.get_best_viterbi_candidate(self, self.vM, self.predecessors1[i],
        #                                                      self.predecessors2[j])
        #     best_X_val, best_X_pos = self.get_best_viterbi_candidate(self, self.vX, self.predecessors1[i], [j])
        #     best_Y_val, best_Y_pos = self.get_best_viterbi_candidate(self, self.vY, [i], self.predecessors2[j])
        #
        # # Multiple predecessors in first graph
        # elif len(self.predecessors1[i]) > 1:
        #     best_M_val, best_M_pos = self.get_best_viterbi_candidate(self, self.vM, self.predecessors1[i], [j)
        #     best_X = self.get_best_viterbi_candidate(self, self.vX, self.predecessors1[i], [j])
        #
        # # Multiple predecessors in second graph
        # elif len(self.predecessors2[j] > 1):
        #     best_Y = self.get_best_viterbi_candidate(self, self.vX, self.predecessors1[i], [j])


        # print (self.predecessors1)
        # print (self.predecessors2)
        # print (i)
        # print (j)
        emissionM = self.get_emission(i, j)

        emissionM = np.log(emissionM)

        curr_transitionMM = float('-Inf') if (best_M_val == float('-Inf'))  else np.log(self.transitionMM) + \
                                                                              best_M_val
        curr_transitionXM = float('-Inf') if best_X_val == float('-Inf')  else np.log(self.transitionXM) + \
                                                                              best_X_val
        curr_transitionYM = float('-Inf') if best_Y_val == float('-Inf') else np.log(self.transitionYM) + \
                                                                              best_Y_val

        if curr_transitionMM >= curr_transitionXM and curr_transitionMM >= curr_transitionYM:
            self.vM[i][j] = curr_transitionMM + emissionM
            self.tracebackM[i][j] = "M-" + str(best_M_pos[0]) + ":" + str(best_M_pos[1])

        elif curr_transitionXM >= curr_transitionMM and curr_transitionXM >= curr_transitionYM:
            self.vM[i][j] = curr_transitionXM + emissionM
            self.tracebackM[i][j] = "X-" + str(best_X_pos[0]) + ":" + str(best_X_pos[1])
        else:
            self.vM[i][j] = curr_transitionYM + emissionM
            self.tracebackM[i][j] = "Y-" + str(best_Y_pos[0]) + ":" + str(best_Y_pos[1])

        # curr_transitionMM = float('-Inf') if (self.vM[i - 1][j - 1] == float('-Inf'))  else np.log(self.transitionMM) + \
        #                                                                       self.vM[i - 1][j - 1]
        # curr_transitionXM = float('-Inf') if self.vX[i - 1][j - 1] == float('-Inf')  else np.log(self.transitionXM) + \
        #                                                                       self.vX[i - 1][j - 1]
        # curr_transitionYM = float('-Inf') if self.vY[i - 1][j - 1] == float('-Inf') else np.log(self.transitionYM) + \
        #                                                                       self.vY[i - 1][j - 1]
        #
        # if curr_transitionMM >= curr_transitionXM and curr_transitionMM >= curr_transitionYM:
        #     self.vM[i][j] = curr_transitionMM + emissionM
        #     self.tracebackM[i][j] = "M"
        #
        # elif curr_transitionXM >= curr_transitionMM and curr_transitionXM >= curr_transitionYM:
        #     self.vM[i][j] = curr_transitionXM + emissionM
        #     self.tracebackM[i][j] = "X"
        # else:
        #     self.vM[i][j] = curr_transitionYM + emissionM
        #     self.tracebackM[i][j] = "Y"



    def fillVX_PO(self, i, j):
        if j == 0:
            return

        best_M_val, best_M_pos = self.get_best_viterbi_candidate(self.vM, [i], self.predecessors2[j])
        best_X_val, best_X_pos = self.get_best_viterbi_candidate(self.vX, [i], self.predecessors2[j])

        curr_transitionXX = float('-Inf') if (best_X_val == float('-Inf') or best_X_val == 0)   else \
            np.log(self.transitionXX) + best_X_val

        curr_transitionMX = float('-Inf') if (best_M_val == float('-Inf') or best_M_val == 0)  else \
            np.log(self.transitionMX) + best_M_val

        if curr_transitionXX >= curr_transitionMX:
            self.vX[i][j] = curr_transitionXX + np.log(self.emissionX)
            self.extensionX[i][j] = curr_transitionXX
            self.tracebackX[i][j] = "X-" + str(best_X_pos[0]) + ":" + str(best_X_pos[1])

        else:
            self.vX[i][j] = curr_transitionMX + np.log(self.emissionX)
            self.openX[i][j] = curr_transitionMX

            self.tracebackX[i][j] = "M-" + str(best_M_pos[0]) + ":" + str(best_M_pos[1])

        # curr_transitionXX = float('-Inf') if (self.vX[i][j - 1] == float('-Inf') or self.vX[i][j-1] == 0)   else \
        #     np.log(
        #     self.transitionXX) + self.vX[i][
        #     j - 1]
        # curr_transitionMX = float('-Inf') if (self.vM[i][j - 1] == float('-Inf') or self.vM[i][j-1] == 0)  else \
        #     np.log(self.transitionMX) + self.vM[i][
        #     j - 1]
        #
        # if curr_transitionXX >= curr_transitionMX:
        #     self.vX[i][j] = curr_transitionXX + np.log(self.emissionX)
        #     self.extensionX[i][j] = curr_transitionXX
        #     self.tracebackX[i][j] = "X"
        #
        # else:
        #     self.vX[i][j] = curr_transitionMX + np.log(self.emissionX)
        #     self.openX[i][j] = curr_transitionMX
        #
        #     self.tracebackX[i][j] = "M"


    def fillVY_PO(self, i, j):
        if i == 0:
            return

        if print_this:

            print ("Checking M transition to Y")

            print ('i and j ')
            print (i)
            print (j)

            print (self.predecessors1)
            print (self.predecessors1[i])

        best_M_val, best_M_pos = self.get_best_viterbi_candidate(self.vM, self.predecessors1[i], [j])

        if print_this:

            print ("Checking Y transition to Y")

        best_Y_val, best_Y_pos = self.get_best_viterbi_candidate(self.vY, self.predecessors1[i], [j])

        curr_transitionYY = float('-Inf') if (best_Y_val == float('-Inf') or best_Y_val == 0)  else \
            np.log(self.transitionYY) + best_Y_val

        curr_transitionMY = float('-Inf') if (best_M_val == float('-Inf') or best_M_val == 0)   else \
            np.log(self.transitionMY) + best_M_val

        if curr_transitionYY >= curr_transitionMY:
            self.vY[i][j] = curr_transitionYY + np.log(self.emissionY)
            self.extensionY[i][j] = curr_transitionYY

            self.tracebackY[i][j] = "Y-" + str(best_Y_pos[0]) + ":" + str(best_Y_pos[1])

        else:
            self.vY[i][j] = curr_transitionMY + np.log(self.emissionY)
            self.openY[i][j] = curr_transitionMY

            self.tracebackY[i][j] = "M-" + str(best_M_pos[0]) + ":" + str(best_M_pos[1])


    def performMEAAlignment(self, po=False):

        self.create_fb_matrices()

        self.forward_algorithm()

        self.backward_algorithm()

        self.calc_posterior_matrix()

        if po:
            self.perform_po_mea_traceback()

        else:
            self.perform_mea_traceback()

    def forward_algorithm(self):

        for i in range(self.m + 1):
            for j in range(self.n + 1):
                if not (i == 0 and j == 0):
                    # if (i == 3 and j == 2):
                    #     print ('target')
                    # print (str(i) + " " +  str(j))

                    self.sumfM(i, j)
                    self.sumfX(i, j)
                    self.sumfY(i, j)

                    # print ("I AND J ARE " + str(i) + " " + str(j))
                    # print ('fM')
                    # print (DataFrame(self.fM))
                    # print ('fX')
                    # print (DataFrame(self.fX))
                    # print ('fY')
                    # print (DataFrame(self.fY))
                    # print (self.fM)
                    # print ('fX')
                    # print (self.fX)
                    # print ('fY')
                    # print (self.fY)

        # print ('filled out')
        # print ('fM')
        # print (DataFrame(self.fM))
        # print ('fX')
        # print (DataFrame(self.fX))
        # print ('fY')
        # print (DataFrame(self.fY))
        # print (self.fM)
        # print ('fX')
        # print (self.fX)
        # print ('fY')
        # print (self.fY)

        if self.log_transform:
            transition_probs = self.log_add(self.fX[self.m][self.n], self.fY[self.m][self.n])

            transition_probs = self.log_add(self.fM[self.m][self.n], transition_probs)

            self.full_prob = np.log(self.tau) + transition_probs

            # print ('tau')
            # print (np.log(self.tau))
            #
            # print ('fx')
            # print (self.fX[self.m][self.n])
            #
            # print ('fy')
            # print (self.fY[self.m][self.n])
            #
            # print ('fm')
            # print (self.fM[self.m][self.n])


        else:
            self.full_prob = self.tau * (self.fM[self.m][self.n] + self.fX[self.m][self.n] + self.fY[self.m][self.n])

        #     print('tau')
        #     print(np.log(self.tau))
        #
        #     print('fx')
        #     print(np.log(self.fX[self.m][self.n]))
        #
        #     print('fy')
        #     print(np.log(self.fY[self.m][self.n]))
        #
        #     print('fm')
        #     print(np.log(self.fM[self.m][self.n]))
        #
        # print ('self full prob is')
        # print (self.full_prob)

        # aligned_profile = self.perform_traceback()
        #
        # return aligned_profile

    def backward_algorithm(self):

        for i in range(self.m, 0, -1):
            for j in range(self.n, 0, -1):
                # if (i == 3 and j == 4):
                #     print ('target')

                if not (i == self.m and j == self.n):
                    self.sumbM(i, j)
                    self.sumbX(i, j)
                    self.sumbY(i, j)

                    # print ("BACKWARDS I AND J ARE " + str(i) + " " + str(j))
                    # print ('bM')
                    # print (DataFrame(self.bM))
                    # print ('bX')
                    # print (DataFrame(self.bX))
                    # print ('bY')
                    # print (DataFrame(self.bY))

                    # print ('bM')
                    # print (self.bM)
                    # print ('bX')
                    # print (self.bX)
                    # print ('bY')
                    # print (self.bY)

    def calc_posterior_matrix(self):

        # print ()

        if print_this:
            print('Full prob is ' + str(self.full_prob))
            print('fM')
            print(DataFrame(self.fM))
            print()
            print('bM')

            print(DataFrame(self.bM))

            print('pM')

        self.pM = np.zeros((self.m + 1, self.n + 1), dtype=float)

        if self.log_transform:
            for i in range(1, self.m + 1):
                for j in range(1, self.n + 1):

                    # if print_this:
                    #     print(self.fM[i][j])
                    #     print(self.bM[i][j])

                    self.pM[i][j] = (self.fM[i][j] + self.bM[i][j]) - self.full_prob

        else:

            for i in range(1, self.m + 1):
                for j in range(1, self.n + 1):

                    if print_this:

                        if self.fM[i][j] == 0:
                            print(0)
                        else:
                            print(np.log(self.fM[i][j]))
                        if self.bM[i][j] == 0:
                            print(0)
                        else:
                            print(np.log(self.bM[i][j]))
                    self.pM[i][j] = (self.fM[i][j] * self.bM[i][j]) / self.full_prob

        # for i in range(1, self.m + 1):
        #     for j in range(1, self.n + 1):
        #         self.pM[i][j] = self.log_add(self.fM[i][j], self.bM[i][j]) / self.full_prob
        #

        if print_this:
            print()
            print(DataFrame(self.pM))

        # self.pM[self.m][self.n] = self.tau

        # for align_pos in self.aligned_positions:
        #     print (align_pos[0], align_pos[1])
        #     print (self.fM[align_pos[0]][align_pos[1]])
        #     print (self.bM[align_pos[0]][align_pos[1]])

    def calc_posterior_for_viterbi(self):
        self.performViterbiAlignment()

        self.performMEAAlignment()

        self.calc_posterior_matrix()

    def sumfM(self, i, j):
        if (i == 0 or j == 0):
            return

        emissionM = self.get_emission(i, j)

        # if i == 1 and j == 2:
        #     print ('woop')

        if self.log_transform:

            if emissionM != 0:
                emissionM = np.log(emissionM)

            # print (self.fM)
            # print (self.fX)
            # print (self.fY)

            forwardMM = float('-Inf') if self.fM[i - 1][j - 1] == float('-Inf') else np.log(
                self.transitionMM) + \
                                                                                     self.fM[
                                                                                         i - 1][
                                                                                         j - 1]
            forwardXM = float('-Inf') if self.fX[i - 1][j - 1] == float('-Inf') or self.fX[i - 1][
                j - 1] == 0 else np.log(
                self.transitionXM) + self.fX[
                                     i - 1][j - 1]
            forwardYM = float('-Inf') if self.fY[i - 1][j - 1] == float('-Inf') or self.fY[i - 1][
                j - 1] == 0 else np.log(
                self.transitionYM) + self.fY[
                                     i - 1][j - 1]

            transition_probs = self.log_add(forwardXM, forwardYM)

            # print(i, j)
            #
            # print ('Log scaled')
            #
            # print(np.log(self.transitionMM))
            #
            # print(self.fM[i - 1][j - 1])
            #
            # print(np.log(self.transitionXM))
            #
            # print(self.fX[i - 1][j - 1])
            #
            # print(np.log(self.transitionYM))
            #
            # print(self.fY[i - 1][j - 1])
            #
            # print(emissionM)

            self.fM[i][j] = emissionM + self.log_add(forwardMM, transition_probs)

            # print ('forward MM, forward XM, forward YM')
            #
            # print(forwardMM)
            # print(forwardXM)
            # print(forwardYM)
            #
            #
            # print (self.fM[i][j])

        else:

            # print (i,j)
            #
            # print ('Normal scale')
            #
            # print (self.transitionMM)
            #
            # print (self.fM[i-1][j-1])
            #
            # print(self.transitionXM)
            #
            # print(self.fX[i - 1][j - 1])
            #
            # print (self.transitionYM)
            #
            # print(self.fY[i - 1][j - 1])
            #
            # print (emissionM)

            forwardMM = float('-Inf') if self.fM[i - 1][j - 1] == -1 else self.transitionMM * self.fM[i - 1][j - 1]
            forwardXM = float('-Inf') if self.fX[i - 1][j - 1] == -1 else self.transitionXM * self.fX[i - 1][j - 1]
            forwardYM = float('-Inf') if self.fY[i - 1][j - 1] == -1 else self.transitionYM * self.fY[i - 1][j - 1]

            # print ('forward MM, forward XM, forward YM')
            #
            #
            # print (forwardMM)
            # print (forwardXM)
            # print (forwardYM)

            self.fM[i][j] = emissionM * (forwardMM + forwardXM + forwardYM)

            # print (self.fM[i][j])
            #
            # print (self.fM)

    def sumfX(self, i, j):

        if j - 1 >= 0:

            if self.log_transform:

                forwardMX = np.log(self.transitionMX) + self.fM[i][j - 1]
                forwardXX = np.log(self.transitionXX) + self.fX[i][j - 1]

                if self.fM[i][j - 1] == float("-Inf") or 0:
                    forwardMX = float("-Inf")

                if self.fX[i][j - 1] == float("-Inf") or 0:
                    forwardXX = float("-Inf")

                if self.fM[i][j - 1] == float("-Inf") and self.fX[i][j - 1] == float("-Inf"):
                    self.fX[i][j] = float("-Inf")
                else:
                    self.fX[i][j] = np.log(self.emissionX) + self.log_add(forwardMX, forwardXX)

            else:

                forwardMX = self.transitionMX * self.fM[i][j - 1]
                forwardXX = self.transitionXX * self.fX[i][j - 1]

                if self.fM[i][j - 1] == float("-Inf"):
                    forwardMX = float("-Inf")

                if self.fX[i][j - 1] == float("-Inf"):
                    forwardXX = float("-Inf")

                if self.fM[i][j - 1] == float("-Inf") and self.fX[i][j - 1] == float("-Inf"):
                    self.fX[i][j] = float("-Inf")
                else:
                    self.fX[i][j] = self.emissionX * (forwardMX + forwardXX)

    def sumfY(self, i, j):

        if i - 1 >= 0:

            if self.log_transform:

                # print (DataFrame(self.fY))
                #
                # print ('fy i-1 j')
                # if self.fY[i-1][j] == 0 :
                #     print (0)
                # else:
                #     print(self.fY[i-1][j])

                forwardMY = np.log(self.transitionMY) + self.fM[i - 1][j]
                forwardYY = np.log(self.transitionYY) + self.fY[i - 1][j]

                if self.fM[i - 1][j] == float("-Inf"):
                    forwardMY = float("-Inf")

                if self.fY[i - 1][j] == float("-Inf"):
                    forwardYY = float("-Inf")

                if self.fM[i - 1][j] == float("-Inf") and self.fY[i - 1][j] == float("-Inf"):
                    # print ('moved' + str(i) + str(j))
                    self.fY[i][j] = float("-Inf")

                else:
                    # print ('here')
                    # print ()
                    # print ('whats doing')
                    # print (i, j)
                    # print (forwardMY)
                    # print (forwardYY)
                    self.fY[i][j] = np.log(self.emissionY) + self.log_add(forwardMY, forwardYY)
                #     print ('prinkles')
                #     print (self.log_add(forwardMY, forwardYY))
                #     print (self.fY[i][j])
                #
                #
                # print ()
                # print ()
                # print ()

                # print ("LOG WHYS" + str(i) + str(j))
                #
                # print (self.fM)
                # print (self.fY)
                #
                # print ('fm i-1 j')
                # if self.fM[i-1][j] == 0 :
                #     print (0)
                # else:
                #
                #     print (self.fM[i-1][j])
                # print ('fy i-1 j')
                # if self.fY[i-1][j] == 0 :
                #     print (0)
                # else:
                #     print(self.fY[i-1][j])
                # print ('forwardMY')
                # print (forwardMY)
                # print ('forwardYY')
                # print (forwardYY)
                #
                # print ('fy')
                # if self.fY[i][j] == 0 :
                #     print (0)
                # else:
                #     print (self.fY[i][j])

            else:

                # print (i,j)
                # print ('hoggles')
                #
                # print (DataFrame(self.fM))
                #
                # print (DataFrame(self.fY))
                #
                #
                # print ('fy i-1 j')
                # if self.fY[i-1][j] == 0 :
                #     print (0)
                # else:
                #     print(self.fY[i-1][j])

                forwardMY = self.transitionMY * self.fM[i - 1][j]
                forwardYY = self.transitionYY * self.fY[i - 1][j]

                if self.fM[i - 1][j] == 0:
                    forwardMY = 0

                if self.fY[i - 1][j] == 0:
                    forwardYY = 0

                if self.fM[i - 1][j] == 0 and self.fY[i - 1][j] == 0:
                    # print ('moved' + str(i) + str(j))
                    self.fY[i][j] = 0
                else:
                    # print ('here')
                    # print ('whats doing')
                    # print (i, j)
                    # print (np.log(forwardMY))
                    # print (np.log(forwardYY))
                    self.fY[i][j] = self.emissionY * (forwardMY + forwardYY)
                    # print ('wrinkles')
                    # print (np.log((forwardMY + forwardYY)))
                    #
                    #
                    # print (np.log((self.fY[i][j])))


    def sumbM(self, i, j):

        emissionM = self.get_emission(i + 1, j + 1)

        # if print_this:

            # print(DataFrame(self.bM))
            # print('i and j')
            # print(i, j)
            #
            # print('i + 1 and j + 1')
            # print(i + 1, j + 1)
            # print('got emission for ')

            # if (i > len(self.profile1.profile) or j > len(self.profile2.profile)):
            #     print('past the end')
            # else:
            #
            #     for char in self.profile1.profile[i - 1].keys():
            #         print(char)
            #
            #     for char in self.profile2.profile[j - 1].keys():
            #         print(char)

        if self.log_transform:

            if emissionM != 0:
                emissionM = np.log(emissionM)

            backwardMM = emissionM + np.log(self.transitionMM) + self.bM[i + 1][j + 1]
            backwardXM = np.log(self.emissionX) + np.log(self.transitionMX) + self.bX[i + 1][j]
            backwardYM = np.log(self.emissionY) + np.log(self.transitionMY) + self.bY[i][j + 1]

            transition_probs = self.log_add(backwardXM, backwardYM)
            #
            # print ()
            # print ()
            # print ()
            #
            # print ("IN LOG HERE BACK AT " + str(i) + str(j))
            #
            # print ('emission X')
            # print (np.log(self.emissionX))
            # print ('transition MX')
            # print (np.log(self.transitionMX))
            # print ('bx i and j+1')
            # print (self.bX[i][j+1])
            #
            # print ('emissionM is ' + str(emissionM))
            # print ('backwardMM is ')
            # print (backwardMM)
            # print ('bacwardXM is ')
            # print (backwardXM)
            # print ('backwardYM is ')
            # print (backwardYM)

            self.bM[i][j] = self.log_add(backwardMM, transition_probs)

            # print ('bm[i][j] is ')
            # print (self.bM[i][j])

        else:

            if i == len(self.profile1.profile) and j == len(self.profile2.profile):
                self.bM[i][j] = self.tau

            backwardMM = emissionM * self.transitionMM * self.bM[i + 1][j + 1]
            backwardXM = self.emissionX * self.transitionMX * self.bX[i + 1][j]
            backwardYM = self.emissionY * self.transitionMY * self.bY[i][j + 1]

            self.bM[i][j] = backwardMM + backwardXM + backwardYM
            #
            # print ()
            # print ()
            # print ()

            if print_this:
                print("IN HERE BACK AT " + str(i) + str(j))

                print('emissionM is ' + str(emissionM))
                print('backwardMM is ')
                print(backwardMM)
                print('bacwardXM is ')
                print(backwardXM)
                print('backwardYM is ')
                print(backwardYM)
            #
            # print('bm[i][j] is ')
            # print(self.bM[i][j])

            # print ('emission X')
            # print (np.log(self.emissionX))
            # print ('transition MX')
            # print (np.log(self.transitionMX))
            # print ('bx i and j+1')
            # print (np.log(self.bX[i][j+1]))
            #
            # print('emissionM is ' + str(emissionM))
            # print('backwardMM is ')
            # print(np.log(backwardMM))
            # print('backwardXM is ')
            # print(np.log(backwardXM))
            # print('backwardYM is ')
            # print(np.log(backwardYM))
            #
            # print('bm[i][j] is ')
            # print(np.log(self.bM[i][j]))

    def sumbX(self, i, j):

        emissionM = self.get_emission(i + 1, j + 1)

        if self.log_transform:

            if emissionM != 0:
                emissionM = np.log(emissionM)
            backwardMX = emissionM + np.log(self.transitionXM) + self.bM[i + 1][j + 1]
            backwardXX = np.log(self.emissionX) + np.log(self.transitionXX) + self.bX[i + 1][j]
            self.bX[i][j] = self.log_add(backwardMX, backwardXX)

        else:

            if i == len(self.profile1.profile) and j == len(self.profile2.profile):
                self.bX[i][j] = self.tau

            backwardMX = emissionM * self.transitionXM * self.bM[i + 1][j + 1]
            backwardXX = self.emissionX * self.transitionXX * self.bX[i + 1][j]
            self.bX[i][j] = backwardMX + backwardXX

    def sumbY(self, i, j):

        emissionM = self.get_emission(i + 1, j + 1)

        if self.log_transform:

            if emissionM != 0:
                emissionM = np.log(emissionM)
            backwardMY = emissionM + np.log(self.transitionYM) + self.bM[i + 1][j + 1]
            backwardYY = np.log(self.emissionY) + np.log(self.transitionYY) + self.bY[i][j + 1]
            self.bY[i][j] = self.log_add(backwardMY, backwardYY)

        else:

            if i == len(self.profile1.profile) and j == len(self.profile2.profile):
                self.bY[i][j] = self.tau

            backwardMY = emissionM * self.transitionYM * self.bM[i + 1][j + 1]
            backwardYY = self.emissionY * self.transitionYY * self.bY[i][j + 1]
            self.bY[i][j] = backwardMY + backwardYY

    def log_add(self, x, y):

        if x == float("-Inf"):
            return y
        if y == float("-Inf"):
            return x
        elif x < y:
            return y + np.log(1 + np.exp(x - y))
        else:
            return x + np.log(1 + np.exp(y - x))

    def perform_viterbi_traceback(self):

        i = len(self.profile1.profile)
        j = len(self.profile2.profile)

        seq1_matches = []
        seq2_matches = []

        seq1_idx = 0
        seq2_idx = 0
        # last_state = ""

        # print(DataFrame(self.vM))
        # print(DataFrame(self.vX))
        # print(DataFrame(self.vY))
        #
        # print(DataFrame(self.tracebackM))
        # print(DataFrame(self.tracebackX))
        # print(DataFrame(self.tracebackY))

        estimated_epsilon = 1

        estimated_delta = 1

        viterbi_path = ''

        # Turn all negative infinity probabilities into zero
        # for i in range(0, self.m + 1):
        #     for j in range(0, self.n + 1):
        #         print ('rhoino')
        #         self.vM[i][j] = 0 if self.vM[i][j] == float("-Inf") else self.vM[i][j]
        #         self.vX[i][j] = 0 if self.vX[i][j] == float("-Inf") else self.vX[i][j]
        #         self.vY[i][j] = 0 if self.vY[i][j] == float("-Inf") else self.vY[i][j]



        while i > 0 and j > 0:
            if self.vM[i][j] > self.vX[i][j] and self.vM[i][j] > self.vY[i][j]:
                seq1_matches.insert(0, i - 1)
                seq2_matches.insert(0, j - 1)
                # last_state = self.tracebackM[i][j]
                seq1_idx = j - 1
                seq2_idx = i - 1
                i -= 1
                j -= 1
                # print ('last state was')
                # print (last_state)

                # print ('M')
                viterbi_path = 'M' + viterbi_path

            elif self.vX[i][j] > self.vM[i][j] and self.vX[i][j] > self.vY[i][j]:
                seq1_matches.insert(0, -1)
                seq2_matches.insert(0, j - 1)
                last_state = self.tracebackX[i][j]
                j -= 1
                # print ('last state was ')
                # print (last_state)

                viterbi_path = 'Y' + viterbi_path


            else:
                seq1_matches.insert(0, i - 1)
                seq2_matches.insert(0, -1)
                # last_state = self.tracebackY[i][j]
                i -= 1
                # print ('last state was')
                # print (last_state)

                viterbi_path = 'X' + viterbi_path

            # while i > 0 and j > 0:
            #     if last_state == "M":
            #         seq1_matches.insert(0, i - 1)
            #         seq2_matches.insert(0, j - 1)

        while seq1_matches[0] > 0 or i > 0:
            seq1_idx = 1 if seq1_idx == 0 else seq1_idx
            seq1_matches.insert(0, i - 1)
            seq2_matches.insert(0, -1)
            i -= 1

            viterbi_path = 'X' + viterbi_path

        while seq2_matches[0] > 0 or j > 0:
            seq2_idx = 1 if seq2_idx == 0 else seq1_idx
            seq2_matches.insert(0, j - 1)
            seq1_matches.insert(0, -1)
            j -= 1

            viterbi_path = 'Y' + viterbi_path

        # print()
        # print ('\nmatches')
        # print (seq1_matches)
        # print (seq2_matches)

        # print('viterbi path')
        # print(viterbi_path)
        #
        # print (self.openY)
        # print()
        # print (self.extensionY)
        #
        # print ("----")
        #
        # print (self.openX)
        # print()
        # print (self.extensionX)

        self.viterbi_aligned_positions = self.get_aligned_positions(seq1_matches, seq2_matches)

        self.viterbi_matches1 = seq1_matches
        self.viterbi_matches2 = seq2_matches

        # alignment1.add_gaps(seq1_matches)
        # alignment2.add_gaps(seq2_matches)

        # print (self.profile1.profile)
        # print (self.profile2.profile)

        # self.profile1.add_profile(self.profile2)

        # return self.profile1

    def get_best_viterbi_candidate(self, matrix, i, j ):

        if print_this:

            print (i)
            print (j)

            print ('Checking best at ' + str(i[0]) + " " + str(j[0]))

        curr_best_score = matrix[i[0]][j[0]]
        curr_best_pos = (i[0],j[0])
        for x in i:
            for y in j:
                # print ("Checking  " + str(x) + " " + str(y))
                # print ('curr best at ' + str(x) + " " + str(y))
                # print (matrix[x][y])
                if matrix[x][y] > curr_best_score:
                    curr_best_score = matrix[x][y]
                    curr_best_pos = (x,y)

        return curr_best_score, curr_best_pos

    def perform_po_viterbi_traceback(self):

        # print ("Perform Viterbi Partial Order traceback")

        i = len(self.profile1.profile)
        j = len(self.profile2.profile)

        seq1_matches = []
        seq2_matches = []

        seq1_idx = 0
        seq2_idx = 0
        # last_state = ""

        # print('vM')
        # print(DataFrame(self.vM))
        # print('vX')
        # print(DataFrame(self.vX))
        # print('vY')
        # print(DataFrame(self.vY))
        # print('traceback M')
        # print(DataFrame(self.tracebackM))
        # print('traceback X')
        #
        # print(DataFrame(self.tracebackX))
        #
        # print('traceback Y')
        #
        # print(DataFrame(self.tracebackY))

        if print_this:

            print ('vM')
            print (DataFrame(self.vM))
            print ('vX')
            print (DataFrame(self.vX))
            print ('vY')
            print (DataFrame(self.vY))
            print ('traceback M')
            print(DataFrame(self.tracebackM))
            print ('traceback X')

            print(DataFrame(self.tracebackX))

            print ('traceback Y')


            print(DataFrame(self.tracebackY))

        estimated_epsilon = 1

        estimated_delta = 1
        viterbi_path = ''

        # print('CONSIDERING')
        # print(i, j)
        # print(i - 1, j - 1)
        # print("or")
        # print(i - 1, j)
        # print('or')
        # print(i, j - 1)
        # print(self.predecessors1)
        # print(self.predecessors2)
        # predecessors_i = self.predecessors1[i]
        # predecessors_j = self.predecessors2[j]
        # print('predecessors i ')
        # print([x for x in predecessors_i])
        # print('predecesors j')
        # print([x for x in predecessors_j])
        #
        # if len(predecessors_j) > 1:
        #     best_x_score, best_x_pos = self.get_best_candidate_viterbi([i], predecessors_j)
        # else:
        #     best_x_score, best_x_pos = self.vX[i][j], (i,j)
        #
        # if len(predecessors_i) > 1:
        #     best_y_score, best_y_pos = self.get_best_candidate_viterbi(predecessors_i, [j])
        # else:
        #     best_y_score, best_y_pos = self.vY[i][j], (i,j)
        #
        #
        # print('best y score')
        # print(best_y_score)
        # print('best y pos')
        # print('best y pos')
        # print(best_y_pos)
        #
        # print('best x score')
        # print(best_x_score)
        # print('best x pos')
        # print(best_x_pos)

        # Work out which matrix we're starting in
        # if self.vM[i][j] > self.vX[i][j] and self.vM[i][j] > self.vY[i][j]:
        #     curr_matrix = 'M'
        # elif self.vX[i][j] > self.vM[i][j] and self.vX[i][j] > self.vY[i][j]:
        #     curr_matrix = 'X'
        # else:
        #     curr_matrix = 'Y'



        last_state = False

        while i > 0 and j > 0:

            # print ('Considering ' + str(i) + " and " + str(j))

            if self.vM[i][j] > self.vX[i][j] and self.vM[i][j] > self.vY[i][j]:
            # if curr_matrix == 'M':

                # viterbi_path = curr_matrix + viterbi_path


                # if not last_state:
                seq1_matches.insert(0, i -1)
                seq2_matches.insert(0, j -1)
                # last_state = False

                curr_matrix = self.tracebackM[i][j].split("-")[0]
                coords = self.tracebackM[i][j].split("-")[1]
                coord_x = int(coords.split(":")[0])
                coord_y = int(coords.split(":")[1])


                i = coord_x
                j = coord_y

                # i -= 1
                # j -= 1


                viterbi_path = 'M' + viterbi_path

            # elif curr_matrix == 'X':

                # viterbi_path = curr_matrix + viterbi_path


            elif self.vX[i][j] > self.vM[i][j] and self.vX[i][j] > self.vY[i][j]:
                if not last_state:

                    seq1_matches.insert(0, -1)
                    seq2_matches.insert(0, j - 1 )
                last_state = False
                # last_state = self.tracebackX[i][j]
                # j -= 1
                # print ('last state was ')
                # print (last_state)

                curr_matrix = self.tracebackX[i][j].split("-")[0]
                coords = self.tracebackX[i][j].split("-")[1]
                coord_x = int(coords.split(":")[0])
                coord_y = int(coords.split(":")[1])

                i = coord_x
                j = coord_y

                viterbi_path = 'X' + viterbi_path


            else:

                # viterbi_path = curr_matrix + viterbi_path


                if not last_state:

                    seq1_matches.insert(0, i -1)
                    seq2_matches.insert(0, -1)

                last_state = False

                # last_state = self.tracebackY[i][j]
                # i -= 1
                # print ('last state was')
                # print (last_state)

                curr_matrix = self.tracebackY[i][j].split("-")[0]
                coords = self.tracebackY[i][j].split("-")[1]
                coord_x = int(coords.split(":")[0])
                coord_y = int(coords.split(":")[1])

                i = coord_x
                j = coord_y


                viterbi_path = 'Y' + viterbi_path




            # while i > 0 and j > 0:
            #     if last_state == "M":
            #         seq1_matches.insert(0, i - 1)
            #         seq2_matches.insert(0, j - 1)

        # print ('matches finishing')
        # print (seq1_matches)
        # print (seq2_matches)
        # print (i)
        # print (j)

        # while seq1_matches[0] > 0 or i > 0:
        #     # seq1_idx = 1 if seq1_idx == 0 else seq1_idx
        #     seq1_matches.insert(0, i-1 )
        #     seq2_matches.insert(0, -1)
        #     i -= 1
        #
        #     viterbi_path = 'Y' + viterbi_path
        #
        # while seq2_matches[0] > 0 or j > 0:
        #     # seq2_idx = 1 if seq2_idx == 0 else seq1_idx
        #     seq2_matches.insert(0, j-1 )
        #     seq1_matches.insert(0, -1)
        #     j -= 1
        #
        #     viterbi_path = 'X' + viterbi_path

        # print()
        # print ('\nmatches')

        # print ('matches before')
        # print (seq1_matches)
        # print (seq2_matches)

        prev = -1



        for idx, pos in enumerate(seq1_matches):
            if pos != -1:
                curr = pos
                if curr - prev != 1:
                    # print (seq1_matches)
                    # print (seq2_matches)
                    # print('curr and prev')
                    # print(curr)
                    # print(prev)
                    # print(idx)
                    seq1_matches.insert(idx, prev + 1)
                    seq2_matches.insert(idx, -1)
                    prev = prev + 1
                else:
                    prev = pos

        prev = -1

        # print ('matches after')
        #
        # print(seq1_matches)
        # print(seq2_matches)

        for idx, pos in enumerate(seq2_matches):
            if pos != -1:
                curr = pos
                if curr - prev != 1:
                    # print (seq1_matches)
                    # print (seq2_matches)
                    # print('curr and prev')
                    # print(curr)
                    # print(prev)
                    # print(idx)
                    seq2_matches.insert(idx, prev + 1)
                    seq1_matches.insert(idx, -1)
                    prev = prev + 1
                else:
                    prev = pos

        # print ('viterbi path')
        # print (viterbi_path)
        #
        # print ('updated seq matches')
        #
        # print (seq1_matches)
        # print (seq2_matches)


        #
        # print (self.openY)
        # print()
        # print (self.extensionY)
        #
        # print ("----")
        #
        # print (self.openX)
        # print()
        # print (self.extensionX)

        self.viterbi_aligned_positions = self.get_aligned_positions(seq1_matches, seq2_matches)

        self.viterbi_matches1 = seq1_matches
        self.viterbi_matches2 = seq2_matches

        # alignment1.add_gaps(seq1_matches)
        # alignment2.add_gaps(seq2_matches)

        # print (self.profile1.profile)
        # print (self.profile2.profile)

        # self.profile1.add_profile(self.profile2)

        # return self.profile1

    def get_best_candidate(self, i, j ):

        hammers = False


        if len(i) > 1 or len(j) > 1:
            hammers = True
            # print ('hammers')

        # if print_this:


        # print ("i ")
        # print (i)
        # print ("j ")
        # print (j)

        # print ('Checking best at ' + str(i[0]) + " " + str(j[0]))

        curr_best_score = self.pM[i[0]][j[0]]
        curr_best_pos = (i[0],j[0])
        for x in i:
            for y in j:
                # print ('curr best at ' + str(x) + " " + str(y))
                # print (self.pM[x][y])
                # print (DataFrame(self.pM))
                if self.pM[x][y] > curr_best_score:
                    # print ('changing best')
                    curr_best_score = self.pM[x][y]
                    curr_best_pos = (x,y)

        return curr_best_score, curr_best_pos


    def perform_po_mea_traceback(self):
        i = len(self.profile1.profile)
        j = len(self.profile2.profile)

        seq1_matches = []
        seq2_matches = []

        seq1_idx = 0
        seq2_idx = 0

        # print (self.predecessors)

        # print (DataFrame(self.vM))
        # print (DataFrame(self.vX))
        #
        # print (DataFrame(self.vY))

        if self.log_transform:

            pMtrace = np.zeros((self.m + 1, self.n + 1), dtype=object)
            pMfill = np.zeros((self.m + 1, self.n + 1), dtype=float)


            if print_this:
                print('fM')

                print(DataFrame(self.fM))

                print('bM')

                print(DataFrame(self.bM))

                print('pM')

                print(DataFrame(self.pM))

                print('doing pmtrace log')

                print(DataFrame(self.pM))

            # if self.predecessors1:
            #     print('yay we have predecesors')
            #     print(self.predecessors1)
            #     print(self.predecessors2)




            for i in range(1, self.m + 1):
                for j in range(1, self.n + 1):

                    # print (self.pM)


                    if self.predecessors1 or self.predecessors2:

                        predecessors_i = self.predecessors1[i]
                        predecessors_j = self.predecessors2[j]

                        # if print_this:
                        #
                        #     print('CONSIDERING')
                        #     print(i, j)
                        #     print(i - 1, j - 1)
                        #     print("or")
                        #     print(i - 1, j)
                        #     print('or')
                        #     print(i, j - 1)
                        #     print(self.predecessors1)
                        #     print(self.predecessors2)
                        #
                        #     print('predecessors i ')
                        #     print([x for x in predecessors_i])
                        #     print('predecesors j')
                        #     print([x for x in predecessors_j])
                        #     print (self.pM)


                        # if len(predecessors_j) > 1 or len(predecessors_i) > 1:
                        #     print ('hotfogs')


                        if len(predecessors_i) == 0:
                            predecessors_i = [i]

                        if len(predecessors_j) == 0:
                            predecessors_j = [j]

                        # print (predecessors_i)
                        # print (predecessors_j)

                        # print (len(predecessors_i))
                        # print (len(predecessors_j))

                        # print("Get best M score")

                        if print_this:
                                print ("Get best M score")

                        best_M_score, best_M_pos = self.get_best_candidate(predecessors_i, predecessors_j)

                        # print("Get best Y score")

                        if print_this:
                                print ("Get best Y score")


                        best_y_score, best_y_pos = self.get_best_candidate(predecessors_i, [j])

                        # print("Get best X score")

                        if print_this:
                                print ("Get best X score")

                        best_x_score, best_x_pos = self.get_best_candidate([i], predecessors_j)

                        if print_this:

                            print ('best y score')
                            print (best_y_score)
                            print ('best y pos')
                            print ('best y pos')
                            print (best_y_pos)

                            print ('best x score')
                            print (best_x_score)
                            print ('best x pos')
                            print (best_x_pos)


                    if self.log_add(best_M_score, self.pM[i][j]) > best_y_score and self.log_add(
                            best_M_score, self.pM[i][j]) > best_x_score:
                        pMtrace[i][j] = 'M-' + str(best_M_pos[0]) + ":" + str(best_M_pos[1])
                        # print('M')

                        self.pM[i][j] = self.log_add(self.pM[i][j], best_M_score)

                    elif best_x_score > self.log_add(best_M_score, self.pM[i][j]) and best_x_score \
                            > best_y_score:
                        # print(DataFrame(self.pM))

                        # print(self.pM[i][j - 1])
                        # print(self.log_add(self.pM[i - 1][j - 1], self.pM[i][j]))
                        pMtrace[i][j] = 'X-' + str(best_x_pos[0]) + ":" + str(best_x_pos[1])
                        # print('X')
                        self.pM[i][j] = best_x_score

                    else:
                        pMtrace[i][j] = 'Y-' + str(best_y_pos[0]) + ":" +  str(best_y_pos[1])
                        # print('Y')
                        self.pM[i][j] = best_y_score



                    # if self.log_add(self.pM[i - 1][j - 1], self.pM[i][j]) > best_y_score and self.log_add(
                    #         self.pM[i - 1][j - 1], self.pM[i][j]) > best_x_score:
                    #     pMtrace[i][j] = 'M'
                    #     # print('M')
                    #
                    #     print ('doing it')
                    #     self.pM[i][j] = self.log_add(self.pM[i][j], self.pM[i - 1][j - 1])
                    #
                    # elif best_x_score > self.log_add(self.pM[i - 1][j - 1], self.pM[i][j]) and best_x_score \
                    #         > best_y_score:
                    #     # print(DataFrame(self.pM))
                    #
                    #     # print(self.pM[i][j - 1])
                    #     # print(self.log_add(self.pM[i - 1][j - 1], self.pM[i][j]))
                    #     pMtrace[i][j] = 'X-' + str(best_x_pos[0]) + ":" + str(best_x_pos[1])
                    #     # print('X')
                    #     self.pM[i][j] = self.pM[i][j - 1]
                    #
                    # else:
                    #     pMtrace[i][j] = 'Y-' + str(best_y_pos[0]) + ":" +  str(best_y_pos[1])
                    #     # print('Y')
                    #     self.pM[i][j] = self.pM[i - 1][j]

            if print_this:
                print('pmtrace log')
                print(pMtrace)

        print('pmtrace log')
        print(pMtrace)

        pm_path = ""

        while i > 0 and j > 0:

            # print (i)
            # print (j)

            if pMtrace[i][j].startswith("M"):
                # print ('Matching')
                seq1_matches.insert(0, i - 1)
                seq2_matches.insert(0, j - 1)
                seq1_idx = j - 1
                seq2_idx = i - 1

                coords = pMtrace[i][j].split("-")[1]
                coord_x = int(coords.split(":")[0])
                coord_y = int(coords.split(":")[1])


                # i -= 1
                # j -= 1

                i = coord_x
                j = coord_y


                pm_path = 'M' + pm_path

            elif pMtrace[i][j].startswith('X'):
                # print ('Gap in x')
                coords = pMtrace[i][j].split("-")[1]
                coord_x = coords.split(":")[0]
                coord_y = coords.split(":")[1]

                # print (pMtrace[i][j])
                # print ('coords')
                # print (coord_x)
                # print (coord_y)
                seq1_matches.insert(0, -1)
                seq2_matches.insert(0, j - 1)
                i = int(coord_x)
                j = int(coord_y)
                # j -= 1
                pm_path = 'X' + pm_path


            else:
                # print ('Gap in y')

                coords = pMtrace[i][j].split("-")[1]
                coord_x = coords.split(":")[0]
                coord_y = coords.split(":")[1]

                seq1_matches.insert(0, i - 1)
                seq2_matches.insert(0, -1)
                # i -= 1

                i = int(coord_x)
                j = int(coord_y)

                pm_path = 'Y' + pm_path

        # print ('Doneskis')
        #
        # print ('before matches')
        # print (seq1_matches)
        # print (seq2_matches)

        # while seq1_matches[0] > 0 or i > 0:
        #     seq1_idx = 1 if seq1_idx == 0 else seq1_idx
        #     seq1_matches.insert(0, i - 1)
        #     seq2_matches.insert(0, -1)
        #     i -= 1
        #     pm_path = 'Y' + pm_path
        #
        # while seq2_matches[0] > 0 or j > 0:
        #     seq2_idx = 1 if seq2_idx == 0 else seq1_idx
        #     seq2_matches.insert(0, j - 1)
        #     seq1_matches.insert(0, -1)
        #     j -= 1
        #
        #     pm_path = 'X' + pm_path

        # print('pm path')
        # print(pm_path)
        # while i > 0 and j > 0:
        #
        #     if self.pM[i - 1][j - 1] + self.pM[i][j] > self.pM[i - 1][j] and self.pM[i - 1][j - 1] + self.pM[i][j] > self.pM[i][j - 1]:
        #         print ('Matching')
        #         seq1_matches.insert(0, i - 1)
        #         seq2_matches.insert(0, j - 1)
        #         seq1_idx = j - 1
        #         seq2_idx = i - 1
        #         i -= 1
        #         j -= 1
        #
        #     elif self.pM[i][j - 1] > self.pM[i - 1][j - 1] + self.pM[i][j] and self.pM[i][j - 1] > self.pM[i - 1][j]:
        #         print ('Gap in x')
        #         seq1_matches.insert(0, -1)
        #         seq2_matches.insert(0, j - 1)
        #         j -= 1
        #
        #     else:
        #         print ('Gap in y')
        #         seq1_matches.insert(0, i - 1)
        #         seq2_matches.insert(0, -1)
        #         i -= 1
        #
        # while seq1_matches[0] > 0 or i > 0:
        #     seq1_idx = 1 if seq1_idx == 0 else seq1_idx
        #     seq1_matches.insert(0, i - 1)
        #     seq2_matches.insert(0, -1)
        #     i -= 1
        #
        # while seq2_matches[0] > 0 or j > 0:
        #     seq2_idx = 1 if seq2_idx == 0 else seq1_idx
        #     seq2_matches.insert(0, j - 1)
        #     seq1_matches.insert(0, -1)
        #     j -= 1

        epsilon_estimate = (pm_path.count('XX') + pm_path.count('YY')) / len(pm_path)
        delta_estimate = (pm_path.count('MX') + pm_path.count('MY')) / len(pm_path)
        emissionX_estimate = (pm_path.count('X'))
        emissionY_estimate = (pm_path.count('Y'))

        self.epsilon = epsilon_estimate
        self.delta = delta_estimate / 2

        # self.emissionX = emissionX_estimate / len(pm_path)
        # self.emissionY = emissionY_estimate / len(pm_path)

        # print (epsilon_estimate)
        # print (delta_estimate)

        self.mea_aligned_positions = self.get_aligned_positions(seq1_matches, seq2_matches)

        seq1_indexes = [-1 if x == -1 else x + 1 for x in seq1_matches]
        seq2_indexes = [-1 if x == -1 else x + 1 for x in seq2_matches]

        seq1_indexes.insert(0, 0)
        seq2_indexes.insert(0, 0)


        # print ('seq 1 indexes')
        # print (seq1_indexes)

        prev = -1

        for idx, pos in enumerate(seq1_matches):
            if pos != -1:
                curr = pos
                if curr - prev != 1:
                    # print (seq1_matches)
                    # print (seq2_matches)
                    # print('curr and prev')
                    # print(curr)
                    # print(prev)
                    # print(idx)
                    seq1_matches.insert(idx, prev + 1)
                    seq2_matches.insert(idx, -1)
                    prev = prev + 1
                else:
                    prev = pos

        prev = -1



        for idx, pos in enumerate(seq2_matches):
            if pos != -1:
                curr = pos
                if curr - prev != 1:
                    seq2_matches.insert(idx, prev + 1)
                    seq1_matches.insert(idx, -1)
                    prev = prev + 1
                else:
                    prev = pos

        # print ('matches after')
        #
        # print(seq1_matches)
        # print(seq2_matches)
        # for idx, x in reversed(list(enumerate(seq1_indexes))):
        #
        #     if x != -1:
        #         if idx not in predecessors[pos]:
        #             predecessors[pos].append(idx)
        #         pos = idx
        #
        # pos = len(seq2_indexes)
        #
        # for idx, x in reversed(list(enumerate(seq2_indexes))):
        #     if x != -1:
        #         if idx not in predecessors[pos]:
        #             predecessors[pos].append(idx)
        #         pos = idx

        self.mea_matches1 = seq1_matches
        self.mea_matches2 = seq2_matches

        # self.predecessors = predecessors
        # Make the predecessors
        # print('made them')
        # print (predecessors)




    def perform_mea_traceback(self):

        i = len(self.profile1.profile)
        j = len(self.profile2.profile)

        seq1_matches = []
        seq2_matches = []

        seq1_idx = 0
        seq2_idx = 0

        if self.log_transform:

            pMtrace = np.zeros((self.m + 1, self.n + 1), dtype=str)
            pMfill = np.zeros((self.m + 1, self.n + 1), dtype=float)
            # print('fM')
            #
            # print(DataFrame(self.fM))
            #
            # print('bM')
            #
            # print(DataFrame(self.bM))
            #
            # print('pM')
            #
            # print(DataFrame(self.pM))


            if print_this:
                print('fM')

                print(DataFrame(self.fM))

                print('bM')

                print(DataFrame(self.bM))

                print('pM')

                print(DataFrame(self.pM))

                print('doing pmtrace log')

                print(DataFrame(self.pM))

            # print('doing pmtrace log')
            #
            # print(DataFrame(self.pM))


            for i in range(1, self.m + 1):
                for j in range(1, self.n + 1):
                    # print(i, j)
                    if self.log_add(self.pM[i - 1][j - 1], self.pM[i][j]) > self.pM[i - 1][j] and self.log_add(
                            self.pM[i - 1][j - 1], self.pM[i][j]) > self.pM[i][j - 1]:
                        pMtrace[i][j] = 'M'
                        # print('M')
                        self.pM[i][j] = self.log_add(self.pM[i - 1][j - 1], self.pM[i][j])
                    elif self.pM[i][j - 1] > self.log_add(self.pM[i - 1][j - 1], self.pM[i][j]) and self.pM[i][j - 1] \
                            > self.pM[i - 1][j]:
                        # print(DataFrame(self.pM))

                        # print(self.pM[i][j - 1])
                        # print(self.log_add(self.pM[i - 1][j - 1], self.pM[i][j]))
                        pMtrace[i][j] = 'X'
                        # print('X')
                        self.pM[i][j] = self.pM[i][j - 1]

                    else:
                        pMtrace[i][j] = 'Y'
                        # print('Y')
                        self.pM[i][j] = self.pM[i - 1][j]

            if print_this:
                print('pmtrace log')
                print(pMtrace)

            # print('pmtrace log')
            # print(pMtrace)
        else:

            pMtrace = np.zeros((self.m + 1, self.n + 1), dtype=str)
            pMfill = np.zeros((self.m + 1, self.n + 1), dtype=float)

            # print('doing pmtrace')
            #
            # print(DataFrame(self.pM))

            for i in range(1, self.m + 1):
                for j in range(1, self.n + 1):
                    if self.pM[i - 1][j - 1] + self.pM[i][j] > self.pM[i - 1][j] and self.pM[i - 1][j - 1] + self.pM[i][
                        j] > self.pM[i][j - 1]:
                        pMtrace[i][j] = 'M'
                        self.pM[i][j] = self.pM[i - 1][j - 1] + self.pM[i][j]
                    elif self.pM[i][j - 1] > self.pM[i - 1][j - 1] + self.pM[i][j] and self.pM[i][j - 1] > \
                            self.pM[i - 1][j]:
                        pMtrace[i][j] = 'X'
                        self.pM[i][j] = self.pM[i][j - 1]

                    else:
                        pMtrace[i][j] = 'Y'
                        self.pM[i][j] = self.pM[i - 1][j]

            #
            if print_this:
                print('pmtrace')
                print(pMtrace)
        pm_path = ""

        while i > 0 and j > 0:

            if pMtrace[i][j] == "M":
                # print ('Matching')
                seq1_matches.insert(0, i - 1)
                seq2_matches.insert(0, j - 1)
                seq1_idx = j - 1
                seq2_idx = i - 1
                i -= 1
                j -= 1
                pm_path = 'M' + pm_path

            elif pMtrace[i][j] == 'X':
                # print ('Gap in x')
                seq1_matches.insert(0, -1)
                seq2_matches.insert(0, j - 1)
                j -= 1
                pm_path = 'X' + pm_path


            else:
                # print ('Gap in y')
                seq1_matches.insert(0, i - 1)
                seq2_matches.insert(0, -1)
                i -= 1
                pm_path = 'Y' + pm_path


        while seq1_matches[0] > 0 or i > 0:
            seq1_idx = 1 if seq1_idx == 0 else seq1_idx
            seq1_matches.insert(0, i - 1)
            seq2_matches.insert(0, -1)
            i -= 1
            pm_path = 'Y' + pm_path

        while seq2_matches[0] > 0 or j > 0:
            seq2_idx = 1 if seq2_idx == 0 else seq1_idx
            seq2_matches.insert(0, j - 1)
            seq1_matches.insert(0, -1)
            j -= 1

            pm_path = 'X' + pm_path

        # print ('pm path')
        # print (pm_path)
        # while i > 0 and j > 0:
        #
        #     if self.pM[i - 1][j - 1] + self.pM[i][j] > self.pM[i - 1][j] and self.pM[i - 1][j - 1] + self.pM[i][j] > self.pM[i][j - 1]:
        #         print ('Matching')
        #         seq1_matches.insert(0, i - 1)
        #         seq2_matches.insert(0, j - 1)
        #         seq1_idx = j - 1
        #         seq2_idx = i - 1
        #         i -= 1
        #         j -= 1
        #
        #     elif self.pM[i][j - 1] > self.pM[i - 1][j - 1] + self.pM[i][j] and self.pM[i][j - 1] > self.pM[i - 1][j]:
        #         print ('Gap in x')
        #         seq1_matches.insert(0, -1)
        #         seq2_matches.insert(0, j - 1)
        #         j -= 1
        #
        #     else:
        #         print ('Gap in y')
        #         seq1_matches.insert(0, i - 1)
        #         seq2_matches.insert(0, -1)
        #         i -= 1
        #
        # while seq1_matches[0] > 0 or i > 0:
        #     seq1_idx = 1 if seq1_idx == 0 else seq1_idx
        #     seq1_matches.insert(0, i - 1)
        #     seq2_matches.insert(0, -1)
        #     i -= 1
        #
        # while seq2_matches[0] > 0 or j > 0:
        #     seq2_idx = 1 if seq2_idx == 0 else seq1_idx
        #     seq2_matches.insert(0, j - 1)
        #     seq1_matches.insert(0, -1)
        #     j -= 1

        epsilon_estimate = (pm_path.count('XX') + pm_path.count('YY')) / len(pm_path)
        delta_estimate = (pm_path.count('MX') + pm_path.count('MY')) / len(pm_path)
        emissionX_estimate = (pm_path.count('X'))
        emissionY_estimate = (pm_path.count('Y'))

        self.epsilon = epsilon_estimate
        self.delta = delta_estimate / 2

        # self.emissionX = emissionX_estimate / len(pm_path)
        # self.emissionY = emissionY_estimate / len(pm_path)

        # print (epsilon_estimate)
        # print (delta_estimate)

        self.mea_aligned_positions = self.get_aligned_positions(seq1_matches, seq2_matches)


        seq1_indexes = [-1 if x == -1 else x + 1 for x in seq1_matches]
        seq2_indexes = [-1 if x == -1 else x + 1 for x in seq2_matches]

        seq1_indexes.insert(0, 0)
        seq2_indexes.insert(0, 0)

        predecessors = defaultdict(list)

        pos = len(seq1_indexes)

        for idx, x in reversed(list(enumerate(seq1_indexes))):

            if x != -1:
                if idx not in predecessors[pos]:
                    predecessors[pos].append(idx)
                pos = idx

        pos = len(seq2_indexes)

        for idx, x in reversed(list(enumerate(seq2_indexes))):
            if x != -1:
                if idx not in predecessors[pos]:
                    predecessors[pos].append(idx)
                pos = idx



        self.mea_matches1 = seq1_matches
        self.mea_matches2 = seq2_matches

        # self.predecessors = predecessors

    def get_predecessors(self):
        return self.predecessors

    def get_aligned_positions(self, matches1, matches2):

        matches = []

        for idx, pos in enumerate(matches1):
            if not (pos == -1) and not matches2[idx] == -1:
                matches.append((pos + 1, matches2[idx] + 1))
                # print (pos + 1)
                # print (matches2[idx] + 1)
        return matches



    def get_emission(self, i, j):
        """
        Get the emission score for aligning two nodes together
        :return: Emission score
        """

        # These will be greater than the profile lengths when we're beginning the backwards algorithm
        if (i > len(self.profile1.profile) or j > len(self.profile2.profile)):
            return 0

        total_count = self.get_total_count(i, j)
        total_score = self.get_total_score(i, j)

        # print (i, j)
        # print ('Total count is ' + str(total_count))
        # print ('Total score is ' + str(total_score))

        emission = total_score / total_count

        # print ('Emission is ' + str(emission))

        return emission

    def get_total_count(self, i, j):
        profile1_count = 0
        profile2_count = 0

        for cnt in self.profile1.profile[i - 1].values():
            profile1_count += cnt

        for cnt in self.profile2.profile[j - 1].values():
            profile2_count += cnt

        # print (profile1_count)
        # print (profile2_count)

        # TODO: Wait, is this right (should it be adding it?)
        return profile1_count * profile2_count



    def get_total_score(self, i, j):

        # print (self.profile1.profile[i-1].keys())
        # print (self.profile2.profile[j-1].keys())

        total_score = 0

        total_num1 = sum(v for k,v in self.profile1.profile[i - 1].items() if k != "-")
        total_num2 = sum(v for k,v in self.profile2.profile[j - 1].items() if k != "-")

        # print ('we are here')
        # print (self.profile1.profile[i - 1].keys())
        # print (self.profile2.profile[j - 1].keys())

        for char in self.profile1.profile[i - 1].keys():
            if char != '-':
                profile1_value = self.profile1.profile[i - 1][char]
                for char2 in self.profile2.profile[j - 1].keys():
                    if char2 != "-":
                        profile2_value = self.profile2.profile[j - 1][char2]
                        match_score = sub_matrix.score_match((char, char2), self.subsmat)
                        # if print_this:
                        # print('Matching ' + str(char) + ' and ' + str(char2))
                        # print ('Score was ' + str(match_score))
                        # print (total_num1, total_num2)
                        #
                        # print (profile1_value / total_num1, profile2_value / total_num2 )

                        total_score += ((profile1_value / total_num1)
                                        * (profile2_value / total_num2)
                                        * match_score)

                        # print ('now total score is ')
                        # print (total_score)

                        # total_score +=  match_score


            # return np.log(total_score)

        # print ('total score for this position is ' + str(total_score))
        # print ('tota')
        return total_score

    def get_alignment(self, type_to_get):

        if type_to_get == "viterbi":

            self.profile1.add_gaps(self.viterbi_matches1)
            self.profile2.add_gaps(self.viterbi_matches2)

        elif type_to_get == "mea":
            self.profile1.add_gaps(self.mea_matches1)
            self.profile2.add_gaps(self.mea_matches2)



        self.profile1.add_profile(self.profile2)
        self.profile1.create_predecessors()

        print ('Alignment is -')

        # print (self.profile1.predecessors)
        print (self.profile1)

        return self.profile1