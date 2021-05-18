# A
# R
# R
# D
import mea_poa.sub_matrix as sub_matrix

print ('A and R')
a_r = sub_matrix.score_match(('A', 'R'), sub_matrix.blosum62EstimatedWithX_dict)
print (a_r)


print ('A and D')
a_d = sub_matrix.score_match(('A', 'D'), sub_matrix.blosum62EstimatedWithX_dict)
print (a_d)

print ('R and R')
r_r = sub_matrix.score_match(('R', 'R'), sub_matrix.blosum62EstimatedWithX_dict)
print (r_r)

print ('R and D')
r_d = sub_matrix.score_match(('R', 'D'), sub_matrix.blosum62EstimatedWithX_dict)
print (r_d)

print ((a_r * a_d * r_r * r_d ) / 4)

print((.5 * .5 * a_r) * (.5 * .5 * a_d) * (.5 * .5 * r_r) * (.5 * .5 * r_d))

