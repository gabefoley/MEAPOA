import subprocess
import mea_poa.align as align
import mea_poa.parameters as params
import mea_poa.utilities as utilities
import os
import pickle
import timeit
import csv
from datetime import datetime

# print(os.getcwd())


def run_qscore(name, aln_type, parameters, save=False, outpath=""):
    base_dir = "./bench1.0/" +  name

    in_dir = base_dir + "/in/"
    ref_dir = base_dir + "/ref/"
    out_dir = "./qscore_alignments/" + aln_type + "_" + name


    files = os.listdir(in_dir)

    file_count = 0

    start_time = timeit.default_timer()

    now = datetime.now()

    dt_string = now.strftime("%Y/%m/%d_%H:%M")

    # Add trailing slash to output directory if it isn't there
    outpath = outpath + "/" if outpath[-1] != "/" else outpath

    param_name = f"t={parameters['tau']}e={parameters['epsilon']}d={parameters['delta']}x={parameters['emissionX']}y={parameters['emissionY']}"

    output_file = "./qscore_alignments/" + aln_type + "_" + name + param_name + ".csv"


    if os.path.exists(outpath + name + ".p"):
        curr_dict = pickle.load(open(outpath + name + ".p", "rb"))
    else:
        curr_dict = {param_name : {}}

    if os.path.exists(outpath + name + "_best.p"):
        best_dict = pickle.load(open(outpath + name + "_best.p", "rb"))
    else:
        best_dict = {}

    if os.path.exists(outpath + "time.p"):
        time_dict = pickle.load(open(outpath + "time.p", "rb"))
    else:
        time_dict = {}

    with open (output_file, 'w+') as output:


        writer = csv.writer(output, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['Tool', 'Dataset', 'Name', 'Q', 'TC', 'M', 'C'])

        # If we don't already have a directory created to save the alignments, lets make one
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        for file in files:

            if file != ".DS_Store":

                if param_name not in curr_dict:
                    curr_dict[param_name] = {}

                # print (curr_dict)
                file_count +=1

                single_time = timeit.default_timer()

                print (file)

                aligned_profile = align.align_seqs(in_dir + file, out_dir + "/" + file + ".aln", aln_type=aln_type,
                                                   params=parameters,
                                                   log_transform=True)


                process = subprocess.Popen("qscore -test %s -ref %s -cline -modeler" % (out_dir + "/" + file + ".aln",
                                                                                        ref_dir + file),
                                           stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

                out = process.communicate()[0]
                errcode = process.returncode

                scores = [x.strip() for x in out.decode('utf-8').split(";")[2:]]

                # scores = [x.split("=")[1] for x in scores]

                print (aligned_profile)

                print('\nScores be')
                print(scores)

                curr_dict[param_name][file] = (scores, aligned_profile)

                update_best_dict(best_dict, file, scores, param_name)


                writer.writerow([aln_type, name, file, scores[0].split("=")[1], scores[1].split("=")[1], scores[2].split("=")[1], scores[3].split("=")[1]])

                # if file not in curr_dict[param_name].keys():
                #     curr_dict[param_name][file] = (scores, aligned_profile)
                # else:
                #     curr_dict[param_name][file] = (scores, aligned_profile)
                #

                total_seconds = timeit.default_timer() - start_time
                single_seconds = timeit.default_timer() - single_time

                if save:

                    pickle.dump(curr_dict, open(outpath + aln_type + "_" + name + ".p", "wb"))
                    pickle.dump(best_dict, open(outpath + aln_type + "_" + name + "_best.p", "wb"))

        if save:

            if name in time_dict:
                if total_seconds < time_dict[name][0]:
                    time_dict[name] = (total_seconds, dt_string)
                    print("New best time - " + utilities.format_time(total_seconds))
            else:
                time_dict[name] = (total_seconds, dt_string)
                print("New best time - " + utilities.format_time(total_seconds))

            pickle.dump(time_dict, open(outpath + aln_type + "_" + "time.p", "wb"))


def update_best_dict(best_dict, file, curr_scores, param_name):


    for score in curr_scores:
        name = score.split("=")[0]
        curr_score = score.split("=")[1]
        if file not in best_dict.keys():

            best_dict[file] = {name : (float(curr_score), param_name)}
            print ("New best score for " + name + " of " + curr_score)

        else:
            if name not in best_dict[file]:
                best_dict[file][name] = (float(curr_score), param_name)
                print ("New best score for " + name + " of " + curr_score)


            elif float(curr_score) > best_dict[file][name][0]:
                best_dict[file][name] = (float(curr_score), param_name)
                print ("New best score for " + name + " of " + curr_score)



def load_best_dict(in_dir, filenames = []):

    if not filenames:
        files = os.listdir(in_dir)
        for name in files:
            if "_best" in name:
                filenames.append(name)

    for name in filenames:
        best_dict = pickle.load(open(in_dir + name, 'rb'))

        for k,v in best_dict.items():
            print (k)
            for v1, v2 in v.items():
                print (v1, v2)

def load_time_dict(in_dir):
    time_dict = pickle.load(open(in_dir + "time.p", "rb"))

    for k, v in time_dict.items():
        print(k, v)



benchmark_names_dna = ['bali2dna', 'bali2dnaf']

benchmark_names_protein = ['bali3', 'bali3pdb', 'bali3pdm', 'ox', 'oxm', 'oxx', 'prefab4', 'prefab4ref',
                           'prefab4refm', 'sabre', 'sabrem']

benchmark_names_all = benchmark_names_dna + benchmark_names_protein

benchmark_names_test = ['bali3_test', 'bali3pdb_test', 'bali3pdm_test', 'ox_test', 'oxm_test', 'oxx_test', 'prefab4_test', 'prefab4ref_test',
                           'prefab4refm_test', 'sabre_test', 'sabrem_test']

aln_types = ['mea', 'viterbi']

benchmark_names_test_short = ['bali3_test']

for name in benchmark_names_test_short:
    for aln_type in aln_types:
        run_qscore(name, aln_type=aln_type, parameters=params.basic_params2, save=True, outpath='./pickle_files/test')

# for name in benchmark_names_protein:
#     run_qscore(name, parameters=params.ba, save=True, outpath='./pickle_files/')

# load_best_dict('./pickle_files/test/', filenames=['viterbi_ox_test.p', 'viterbi_sabre_test.p'])
#
# load_time_dict('./pickle_files/test')