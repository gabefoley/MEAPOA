import mea_poa.sub_matrix as sub_matrix
import mea_poa.alignment_profile as aln_profile
import numpy as np
from pandas import DataFrame
import pandas as pd
print_this = True

pd.set_option('display.max_columns', None)



class PairHMM():
    def __init__(self, profiles, tau, epsilon, delta, emissionX, emissionY, subsmat, log_transform):
        self.tau = tau
        self.epsilon = epsilon
        self.delta = delta
        self.profiles=profiles
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


    # Initialise the transition probabilities
    def initialise_probs(self):
        self.transitionMM = 1 - (2 * self.delta) - self.tau;
        self.transitionMX = self.transitionMY = self.delta;
        self.transitionXX = self.transitionYY = self.epsilon;
        self.transitionXM = self.transitionYM = 1 - self.epsilon - self.tau;

        # print (self.transitionMM)
        # print (self.transitionMX)
        # print (self.transitionXX)
        # print (self.transitionXM)

    # Initialise the matrices
    def create_viterbi_matrices(self):

        self.vM = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self.vY = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self.vX = np.zeros((self.m + 1, self.n + 1), dtype=float)

        self.tracebackM = np.zeros((self.m + 1, self.n + 1), dtype=str)
        self.tracebackX = np.zeros((self.m + 1, self.n + 1), dtype=str)
        self.tracebackY = np.zeros((self.m + 1, self.n + 1), dtype=str)




        if self.log_transform:
            self.vM[0][0] = np.log(1)

            for i in range(1, self.m + 1):
                self.vX[i][0] = float('-Inf')
                self.vM[i][0] = float('-Inf')

            for i in range(1, self.n + 1):
                self.vY[0][i] = float('-Inf')
                self.vM[0][i] = float('-Inf')

        else:

            self.vM[0][0] = 1
            for i in range(1, self.m + 1):
                self.vX[i][0] = 0
                self.vM[i][0] = 0

            for i in range(1, self.n + 1):
                self.vY[0][i] = 0
                self.vM[0][i] = 0



            # print ()
        # print ("Initial X")
        # print (self.vX)
        #
        # print()
        #
        # print("Initial Y")
        # print(self.vY)

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

                self.bX[i][self.n+1] = float('-Inf')
                self.bY[i][self.n+1] = float('-Inf')

                self.bM[i][self.n+1] = float('-Inf')

            for i in range(0, self.n + 1):
                self.bY[self.m+1][i] = float('-Inf')
                self.bX[self.m+1][i] = float('-Inf')

                self.bM[self.m+1][i] = float('-Inf')

            self.fM[0][0] = np.log(1)




            self.bM[self.m][self.n] = np.log(self.tau)
            self.bX[self.m][self.n] = np.log(self.tau)
            self.bY[self.m][self.n] = np.log(self.tau)

            # print ('woowwww')
            # print (self.fM)
            # print (self.fX)
            # print (self.fY)


        else:

            self.fM[0][0] = 1

            for i in range(1, self.m + 1):
                self.fX[i][0] = 0
                self.fM[i][0] = 0


            for i in range(1, self.n + 1):
                self.fY[0][i] = 0
                self.fM[0][i] = 0

            for i in range(0, self.m + 1):

                self.bX[i][self.n+1] = 0
                self.bM[i][self.n+1] = 0
            for i in range(0, self.n + 1):

                self.bY[self.m+1][i] = 0
                self.bM[self.m+1][i] = 0



            self.bM[self.m][self.n] = self.tau
            self.bX[self.m][self.n] = self.tau
            self.bY[self.m][self.n] = self.tau




    def performViterbiAlignment(self):

        self.create_viterbi_matrices()

        for i in range(self.m + 1):
            for j in range(self.n + 1):
                if not (i == 0 and j == 0):

                    # print (str(i) + " " +  str(j))

                    self.fillVM(i, j)
                    self.fillVX(i, j)
                    self.fillVY(i, j)

                    # print (str(i) + " " + str(j))
                    # print ('vM')
                    # print (self.vM)
                    # print ('vX')
                    # print (self.vX)
                    # print ('vY')
                    # print (self.vY)

        # print ('vM')
        # print (DataFrame(self.vM))
        # print ('vX')
        # print (DataFrame(self.vX))
        # print ('vY')
        # print (DataFrame(self.vY))

        aligned_profile = self.perform_viterbi_traceback()

        return aligned_profile


    def fillVM(self, i, j):
        if (i == 0 or j == 0):
            return
        if self.log_transform:
            emissionM = self.get_emission(i, j)

            # print (self.transitionXM)

            curr_transitionMM = float('-Inf') if self.vM[i-1][j-1] == -1 else np.log(self.transitionMM) + self.vM[i-1][j-1]
            curr_transitionXM = float('-Inf') if self.vX[i-1][j-1] == -1 else np.log(self.transitionXM) + self.vX[i-1][j-1]
            curr_transitionYM = float('-Inf') if self.vY[i-1][j-1] == -1 else np.log(self.transitionYM) + self.vY[i-1][j-1]

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

            curr_transitionMM = float('-Inf') if self.vM[i-1][j-1] == -1 else self.transitionMM * self.vM[i-1][j-1]
            curr_transitionXM = float('-Inf') if self.vX[i-1][j-1] == -1 else self.transitionXM * self.vX[i-1][j-1]
            curr_transitionYM = float('-Inf') if self.vY[i-1][j-1] == -1 else self.transitionYM * self.vY[i-1][j-1]


            if curr_transitionMM >= curr_transitionXM and curr_transitionMM >= curr_transitionYM:
                self.vM[i][j] = curr_transitionMM * emissionM
                self.tracebackM[i][j] = "M"
            elif curr_transitionXM >= curr_transitionMM and curr_transitionXM >= curr_transitionYM:
                self.vM[i][j] = curr_transitionXM * emissionM
                self.tracebackM[i][j] = "X"
            else:
                # print (curr_transitionYM)
                # print (emissionM)
                self.vM[i][j] = curr_transitionYM * emissionM
                self.tracebackM[i][j] = "Y"


    def fillVX(self, i, j):
        if i == 0:
            return

        if self.log_transform:
            curr_transitionXX = float('-Inf') if self.vX[i][j-1] == -1 else np.log(self.transitionXX) + self.vX[i][j-1]
            curr_transitionMX = float('-Inf') if self.vM[i][j-1] == -1 else np.log(self.transitionMX) + self.vM[i][j-1]

            if curr_transitionXX >= curr_transitionMX:
                self.vX[i][j] = curr_transitionXX + np.log(self.emissionX)
                self.tracebackX[i][j] = "X"

            else:
                self.vX[i][j] = curr_transitionMX + np.log(self.emissionX)
                self.tracebackX[i][j] = "M"
        else:
            curr_transitionXX = float('-Inf') if self.vX[i][j - 1] == -1 else self.transitionXX * self.vX[i][j - 1]
            curr_transitionMX = float('-Inf') if self.vM[i][j - 1] == -1 else self.transitionMX * self.vM[i][j - 1]

            if curr_transitionXX >= curr_transitionMX:
                self.vX[i][j] = curr_transitionXX * self.emissionX
                self.tracebackX[i][j] = "X"

            else:
                self.vX[i][j] = curr_transitionMX * self.emissionX
                self.tracebackX[i][j] = "M"


    def fillVY(self, i, j):
        if j == 0:
            return


        if self.log_transform:

            curr_transitionYY = float('-Inf') if self.vY[i-1][j] == -1 else np.log(self.transitionYY) + self.vY[i-1][j]
            curr_transitionMY = float('-Inf') if self.vM[i-1][j] == -1 else np.log(self.transitionMY) + self.vM[i-1][j]

            if curr_transitionYY >= curr_transitionMY:
                self.vY[i][j] = curr_transitionYY + np.log(self.emissionY)
                self.tracebackY[i][j] = "Y"

            else:
                self.vY[i][j] = curr_transitionMY + np.log(self.emissionY)
                self.tracebackY[i][j] = "M"
        else:
            curr_transitionYY = float('-Inf') if self.vY[i - 1][j] == -1 else self.transitionYY * self.vY[i - 1][j]
            curr_transitionMY = float('-Inf') if self.vM[i - 1][j] == -1 else self.transitionMY * self.vM[i - 1][j]

            if curr_transitionYY >= curr_transitionMY:
                self.vY[i][j] = curr_transitionYY * self.emissionY
                self.tracebackY[i][j] = "Y"

            else:
                self.vY[i][j] = curr_transitionMY * self.emissionY
                self.tracebackY[i][j] = "M"

    def performMEAAlignment(self):

        self.create_fb_matrices()

        self.forward_algorithm()

        self.backward_algorithm()

        self.calc_posterior_matrix()

        self.perform_mea_traceback()

        # print ('fm')
        #
        # print (DataFrame(self.fM[0:2]))
        #
        # print ('bm')
        #
        # print (DataFrame(self.bM[0:2]))



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
            for j in range(self.n , 0, -1):
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
            print ('Full prob is ' + str(self.full_prob))
            print ('fM')
            print (DataFrame(self.fM))
            print()
            print ('bM')

            print (DataFrame(self.bM))

            print ('pM')

        self.pM = np.zeros((self.m + 1, self.n + 1), dtype=float)

        if self.log_transform:
            for i in range(1, self.m + 1):
                for j in range(1, self.n + 1):

                    if print_this:

                        print (self.fM[i][j])
                        print (self.bM[i][j])

                    self.pM[i][j] = (self.fM[i][j] + self.bM[i][j]) - self.full_prob

        else:

            for i in range(1, self.m + 1):
                for j in range (1, self.n + 1):

                    if print_this:

                        if self.fM[i][j] == 0:
                            print (0)
                        else:
                            print (np.log(self.fM[i][j]))
                        if self.bM[i][j] == 0:
                            print (0)
                        else:
                            print (np.log(self.bM[i][j]))
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



            forwardMM = float('-Inf') if self.fM[i-1][j-1] == float('-Inf') else np.log(
                self.transitionMM) + \
                                                                                     self.fM[
                i-1][
                j-1]
            forwardXM = float('-Inf') if self.fX[i-1][j-1] == float('-Inf') or self.fX[i-1][j-1] == 0   else np.log(
                self.transitionXM) + self.fX[
                i-1][j-1]
            forwardYM = float('-Inf') if self.fY[i-1][j-1] == float('-Inf') or self.fY[i-1][j-1] == 0  else np.log(
                self.transitionYM) + self.fY[
                i-1][j-1]


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



            forwardMM = float('-Inf') if self.fM[i-1][j-1] == -1 else self.transitionMM * self.fM[i-1][j-1]
            forwardXM = float('-Inf') if self.fX[i-1][j-1] == -1 else self.transitionXM * self.fX[i-1][j-1]
            forwardYM = float('-Inf') if self.fY[i-1][j-1] == -1 else self.transitionYM * self.fY[i-1][j-1]

            # print ('forward MM, forward XM, forward YM')
            #
            #
            # print (forwardMM)
            # print (forwardXM)
            # print (forwardYM)


            self.fM[i][j] =  emissionM * (forwardMM + forwardXM + forwardYM)

            # print (self.fM[i][j])
            #
            # print (self.fM)



    def sumfX(self, i, j):

        if j - 1 >= 0:

            if self.log_transform:


                forwardMX = np.log(self.transitionMX) + self.fM[i][j-1]
                forwardXX = np.log(self.transitionXX) + self.fX[i][j-1]

                if self.fM[i][j-1] == float("-Inf") or 0:
                    forwardMX = float("-Inf")

                if self.fX[i][j-1] == float("-Inf") or 0:
                    forwardXX = float("-Inf")

                if self.fM[i][j-1] == float("-Inf")and self.fX[i][j-1] == float("-Inf"):
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



                forwardMY = np.log(self.transitionMY) + self.fM[i-1][j]
                forwardYY = np.log(self.transitionYY) + self.fY[i-1][j]

                if self.fM[i-1][j] == float("-Inf")  :
                    forwardMY = float("-Inf")

                if self.fY[i-1][j] == float("-Inf") :
                    forwardYY = float("-Inf")

                if self.fM[i-1][j] == float("-Inf") and self.fY[i-1][j] == float("-Inf"):
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
                #
                # print ()
                # print ()
                # print ()
                #
                # print("WHYS" + str(i) + str(j))
                #
                # print (self.fM)
                # print (self.fY)
                #
                # print ('fm i-1 j')
                # print (self.fM[i-1][j])
                # print ('dogs')
                # if self.fM[i-1][j] == 0 :
                #     print (0)
                #
                # else:
                #     print ('crushed')
                #     print (np.log(self.fM[i-1][j]))
                # print ('fy i-1 j')
                # if self.fY[i-1][j] == 0:
                #     print(0)
                # else:
                #     print(np.log(self.fY[i-1][j]))
                #
                # print('forwardMY')
                #
                # if forwardMY == 0 or forwardMY == float("-Inf"):
                #     print(forwardMY)
                # else:
                #     print(np.log(forwardMY))
                # print('forwardYY')
                # if forwardYY == 0 or forwardYY == float("-Inf"):
                #     print(forwardYY)
                # else:
                #     print(np.log(forwardYY))
                #
                # print('fy')
                # if self.fY[i][j] == 0 :
                #     print (0)
                # else:
                #     print(np.log(self.fY[i][j]))

    def sumbM(self, i, j ):





        emissionM = self.get_emission(i, j)



        if print_this:

            print (self.bM)
            print ('i and j')
            print (i,j)

            print ('i + 1 and j + 1')
            print (i+1, j+1)
            print ('got emission for ')

            if (i > len(self.profile1.profile) or j > len(self.profile2.profile)):
                print ('past the end')
            else:

                for char in self.profile1.profile[i-1].keys():
                    print (char)

                for char in self.profile2.profile[j-1].keys():
                    print (char)



        if self.log_transform:

            if emissionM != 0:

                emissionM = np.log(emissionM)


            backwardMM = emissionM + np.log(self.transitionMM) + self.bM[i+1][j+1]
            backwardXM = np.log(self.emissionX) + np.log(self.transitionMX) + self.bX[i][j+1]
            backwardYM = np.log(self.emissionY) + np.log(self.transitionMY) + self.bY[i+1][j]

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

            self.bM[i][j] =  self.log_add(backwardMM, transition_probs)

            # print ('bm[i][j] is ')
            # print (self.bM[i][j])

        else:

            if i == len(self.profile1.profile) and j == len(self.profile2.profile):
                self.bM[i][j] = self.tau


            backwardMM = emissionM * self.transitionMM * self.bM[i + 1][j + 1]
            backwardXM = self.emissionX * self.transitionMX * self.bX[i+1][j]
            backwardYM = self.emissionY * self.transitionMY * self.bY[i][j+1]

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




    def sumbX(self, i, j ):

        emissionM = self.get_emission(i+1, j+1)


        if self.log_transform:

            if emissionM != 0:

                emissionM = np.log(emissionM)
            backwardMX = emissionM + np.log(self.transitionXM) + self.bM[i+1][j+1]
            backwardXX = np.log(self.emissionX) + np.log(self.transitionXX) + self.bX[i+1][j]
            self.bX[i][j] = self.log_add(backwardMX, backwardXX)

        else:

            if i == len(self.profile1.profile) and j == len(self.profile2.profile):
                self.bX[i][j] = self.tau

            backwardMX = emissionM * self.transitionXM * self.bM[i+1][j+1]
            backwardXX = self.emissionX * self.transitionXX * self.bX[i+1][j]
            self.bX[i][j] = backwardMX + backwardXX


    def sumbY(self, i, j ):

        emissionM = self.get_emission(i+1, j+1)


        if self.log_transform:

            if emissionM != 0:

                emissionM = np.log(emissionM)
            backwardMY = emissionM + np.log(self.transitionYM) + self.bM[i+1][j+1]
            backwardYY = np.log(self.emissionY) + np.log(self.transitionYY) + self.bY[i][j+1]
            self.bY[i][j] = self.log_add(backwardMY, backwardYY)

        else:

            if i == len(self.profile1.profile) and j == len(self.profile2.profile):
                self.bY[i][j] = self.tau

            backwardMY = emissionM * self.transitionYM * self.bM[i+1][j+1]
            backwardYY = self.emissionY * self.transitionYY * self.bY[i][j+1]
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
        last_state = ""

        while i > 0 and j > 0:
            if self.vM[i][j] > self.vX[i][j] and self.vM[i][j] > self.vY[i][j]:
                seq1_matches.insert(0, i-1)
                seq2_matches.insert(0, j-1)
                last_state = self.tracebackM[i][j]
                seq1_idx = j - 1
                seq2_idx = i - 1
                i -= 1
                j -= 1

            elif self.vX[i][j] > self.vM[i][j] and self.vX[i][j] > self.vY[i][j]:
                seq1_matches.insert(0, -1)
                seq2_matches.insert(0, j-1)
                last_state = self.tracebackX[i][j]
                j -= 1

            else:
                seq1_matches.insert(0, i-1)
                seq2_matches.insert(0, -1)
                last_state = self.tracebackY[i][j]
                i -= 1

            # while i > 0 and j > 0:
            #     if last_state == "M":
            #         seq1_matches.insert(0, i - 1)
            #         seq2_matches.insert(0, j - 1)

        while seq1_matches[0] > 0 or i > 0:
            seq1_idx = 1 if seq1_idx == 0 else seq1_idx
            seq1_matches.insert(0, i-1)
            seq2_matches.insert(0, -1)
            i -= 1

        while seq2_matches[0] > 0 or j > 0:
            seq2_idx = 1 if seq2_idx == 0 else seq1_idx
            seq2_matches.insert(0, j-1)
            seq1_matches.insert(0, -1)
            j -= 1
        # print()
        # print ('\nmatches')
        # print (seq1_matches)
        # print (seq2_matches)

        self.viterbi_aligned_positions = self.get_aligned_positions(seq1_matches, seq2_matches)

        self.viterbi_matches1 = seq1_matches
        self.viterbi_matches2 = seq2_matches

        # alignment1.add_gaps(seq1_matches)
        # alignment2.add_gaps(seq2_matches)

        # print (self.profile1.profile)
        # print (self.profile2.profile)

        # self.profile1.add_profile(self.profile2)

        # return self.profile1

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

            if print_this:

                print ('doing pmtrace log')

                print(DataFrame(self.pM))

            for i in range(1, self.m + 1):
                for j in range(1, self.n + 1):
                    # print (i, j )
                    if self.log_add(self.pM[i - 1][j - 1], self.pM[i][j]) > self.pM[i - 1][j] and self.log_add(
                            self.pM[i - 1][j -1], self.pM[i][j]) > self.pM[i][j - 1]:
                        pMtrace[i][j] = 'M'
                        self.pM[i][j] = self.log_add(self.pM[i - 1][j - 1], self.pM[i][j])
                    elif self.pM[i][j - 1] > self.log_add(self.pM[i - 1][j - 1], self.pM[i][j]) and self.pM[i][j - 1]\
                            > self.pM[i - 1][j]:
                        pMtrace[i][j] = 'X'
                        self.pM[i][j] = self.pM[i][j - 1]

                    else:
                        pMtrace[i][j] = 'Y'
                        self.pM[i][j] = self.pM[i - 1][j]

            if print_this:
                print ('pmtrace log')
                print (pMtrace)


        else:



            pMtrace = np.zeros((self.m + 1, self.n + 1), dtype=str)
            pMfill = np.zeros((self.m + 1, self.n + 1), dtype=float)

            # print('doing pmtrace')
            #
            # print(DataFrame(self.pM))

            for i in range(1, self.m + 1):
                for j in range (1, self.n + 1):
                    if self.pM[i - 1][j - 1] + self.pM[i][j] > self.pM[i - 1][j] and self.pM[i - 1][j - 1] + self.pM[i][
                        j] > self.pM[i][j - 1]:
                        pMtrace[i][j] = 'M'
                        self.pM[i][j] = self.pM[i-1][j-1] + self.pM[i][j]
                    elif self.pM[i][j - 1] > self.pM[i - 1][j - 1] + self.pM[i][j] and self.pM[i][j - 1] > self.pM[i - 1][j]:
                        pMtrace[i][j] = 'X'
                        self.pM[i][j] = self.pM[i][j-1]

                    else:
                        pMtrace[i][j] = 'Y'
                        self.pM[i][j] = self.pM[i-1][j]

            #
            if print_this:
                print ('pmtrace')
                print (pMtrace)

        while i > 0 and j > 0:

            if pMtrace[i][j] == "M":
                # print ('Matching')
                seq1_matches.insert(0, i - 1)
                seq2_matches.insert(0, j - 1)
                seq1_idx = j - 1
                seq2_idx = i - 1
                i -= 1
                j -= 1

            elif pMtrace[i][j] == 'X':
                # print ('Gap in x')
                seq1_matches.insert(0, -1)
                seq2_matches.insert(0, j - 1)
                j -= 1

            else:
                # print ('Gap in y')
                seq1_matches.insert(0, i - 1)
                seq2_matches.insert(0, -1)
                i -= 1

        while seq1_matches[0] > 0 or i > 0:
            seq1_idx = 1 if seq1_idx == 0 else seq1_idx
            seq1_matches.insert(0, i - 1)
            seq2_matches.insert(0, -1)
            i -= 1

        while seq2_matches[0] > 0 or j > 0:
            seq2_idx = 1 if seq2_idx == 0 else seq1_idx
            seq2_matches.insert(0, j - 1)
            seq1_matches.insert(0, -1)
            j -= 1



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

        self.mea_aligned_positions = self.get_aligned_positions(seq1_matches, seq2_matches)
        self.mea_matches1 = seq1_matches
        self.mea_matches2 = seq2_matches

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

        for cnt in self.profile1.profile[i-1].values():
            profile1_count += cnt

        for cnt in self.profile2.profile[j-1].values():
            profile2_count += cnt


        #TODO: Wait, is this right (should it be adding it?)
        return profile1_count * profile2_count

    def get_total_score(self, i, j):

        # print (self.profile1.profile[i-1].keys())
        # print (self.profile2.profile[j-1].keys())

        total_score = 0
        for char in self.profile1.profile[i - 1].keys():
            if char != '-':
                profile1_value = self.profile1.profile[i-1][char]
                for char2 in self.profile2.profile[j-1].keys():
                    if char2 != "-":

                        profile2_value = self.profile2.profile[j-1][char2]
                        match_score = sub_matrix.score_match((char, char2), self.subsmat)
                        # if print_this:
                        #     print('Matching ' + str(char) + ' and ' + str(char2))
                        #     print ('Score was ' + str(match_score))

                        total_score += profile1_value * profile2_value * match_score

        # return np.log(total_score)



            return total_score



    def get_alignment(self, type_to_get):

        if type_to_get == "viterbi":

            self.profile1.add_gaps(self.viterbi_matches1)
            self.profile2.add_gaps(self.viterbi_matches2)


        elif type_to_get == "mea":
            self.profile1.add_gaps(self.mea_matches1)
            self.profile2.add_gaps(self.mea_matches2)


        self.profile1.add_profile(self.profile2)
        return self.profile1









