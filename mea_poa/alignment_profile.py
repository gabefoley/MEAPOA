class AlignmentProfile():

    def __init__(self, seqs):
        self.seqs = seqs
        self.profile = {}

        self.create_profile()

    def __str__(self):
        seq_str = ""
        for seq in self.seqs:
            seq_str += ">" + seq.name + "\n" + "".join(seq.sequence) + "\n"
        return seq_str

    def create_profile(self):
        for seq in self.seqs:
            for idx, char in enumerate(seq):
                if idx in self.profile:
                    if char in self.profile[idx]:
                        self.profile[idx][char] +=1
                    else:
                        self.profile[idx][char] = 1
                else:
                    self.profile[idx] = {char: 1}

    def add_profile(self, to_add):
        if len(self.profile) != len(to_add.profile):
            print ("Trying to add two profiles of unequal size")
            return
        else:
            for seq in to_add.seqs:
                # print ('adding')
                # print (seq)
                # print (type(seq))
                self.seqs.append(seq)
            for idx in to_add.profile.keys():
                for char in to_add.profile[idx]:
                    if char in self.profile[idx]:
                        self.profile[idx][char] += to_add.profile[idx][char]
                    else:
                        self.profile[idx][char] = to_add.profile[idx][char]

    def add_gaps(self, gap_pos):

        for seq in self.seqs:
            for idx, pos in enumerate(gap_pos):
                if pos == -1:
                    seq.sequence.insert(idx, "-")

        self.profile = {}
        self.create_profile()


