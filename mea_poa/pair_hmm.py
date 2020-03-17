import mea_poa.sub_matrix as sub_matrix
import mea_poa.alignment_profile as aln_profile
import numpy as np

class PairHMM():
    def __init__(self, profiles, tau, epsilon, delta, emissionX, emissionY, subsmat, log_transform):
        self.tau = tau
        self.epsilon = epsilon
        self.delta = delta
        self.profiles=profiles
        self.profile1 = profiles[0]
        self.profile2 = profiles[1]
        self.aligned_seqs = []
        self.subsmat = subsmat
        self.log_transform = log_transform
        self.emissionX = emissionX
        self.emissionY = emissionY
        self.vM = None
        self.vX = None
        self.vY = None
        self.initialise_probs()
        self.create_matrices()


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
    def create_matrices(self):

        m = len(self.profile1.profile)
        n = len(self.profile2.profile)
        self.vM = np.zeros((m + 1, n + 1), dtype=float)
        self.vY = np.zeros((m + 1, n + 1), dtype=float)
        self.vX = np.zeros((m + 1, n + 1), dtype=float)

        self.tracebackM = np.zeros((m + 1, n + 1), dtype=str)
        self.tracebackX = np.zeros((m + 1, n + 1), dtype=str)
        self.tracebackY = np.zeros((m + 1, n + 1), dtype=str)

        for i in range(1, m + 1):
            self.vY[i][0] = float('-Inf')
            self.vM[i][0] = float('-Inf')

        for i in range(1, n + 1):
            self.vX[0][i] = float('-Inf')
            self.vM[0][i] = float('-Inf')

        self.vM[0][0] = 1

        if self.log_transform:
            self.vM[0][0] = np.log(self.vM[0][0])

        # print ()
        # print ("Initial X")
        # print (self.vX)
        #
        # print()
        #
        # print("Initial Y")
        # print(self.vY)

    def performViterbiAlignment(self):
        m = len(self.profile1.profile) + 1
        n = len(self.profile2.profile) + 1
        for i in range(m):
            for j in range(n):
                if not (i == 0) or not (j == 0):

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


        aligned_profile = self.perform_traceback()

        return aligned_profile


    def fillVM(self, i, j):
        if (i == 0 or j == 0):
            return
        if self.log_transform:
            emissionM = self.get_emission(i, j)

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



    def perform_traceback(self):
        i = len(self.profile1.profile)
        j = len(self.profile2.profile)

        # print (' M')
        # print (self.vM)
        #
        # print (' X')
        #
        # print (self.vX)
        #
        # print (' Y')
        #
        # print (self.vY)
        #
        # print ('Traceback M')
        # print (self.tracebackM)
        #
        # print ('Traceback X')
        #
        # print (self.tracebackX)
        #
        # print ('Traceback Y')
        #
        # print (self.tracebackY)

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
        # print ('matches')
        # print (seq1_matches)
        # print (seq2_matches)
        #
        # print(len(seq1_matches))
        # print (len(seq2_matches))

        self.profile1.add_gaps(seq1_matches)
        self.profile2.add_gaps(seq2_matches)

        # print (self.profile1.profile)
        # print (self.profile2.profile)

        self.profile1.add_profile(self.profile2)

        return self.profile1

    def get_emission(self, i, j):
        """
        Get the emission score for aligning two nodes together
        :return: Emission score
        """

        total_count = self.get_total_count(i, j)

        # print ('total count ')
        # print (total_count)

        total_score = self.get_total_score(i, j)

        # print ('total score')
        # print (total_score)

        emission = total_score / total_count

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
        total_score = 0
        for char in self.profile1.profile[i - 1].keys():
            if char != '-':
                profile1_value = self.profile1.profile[i-1][char]
                for char2 in self.profile2.profile[j-1].keys():
                    if char2 != "-":
                        profile2_value = self.profile2.profile[j-1][char2]
                        match_score = sub_matrix.score_match((char, char2), self.subsmat)

                        total_score += profile1_value * profile2_value * match_score

        # return np.log(total_score)

            return total_score



