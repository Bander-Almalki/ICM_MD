import math
import os
import re
import csv

def parse_LIPS_score( LIPS_output_file, LIPS_parsed_csv):
    """Parse the LIPS output to create a CSV with features for input in machine learning algorithm.

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt, PDB)
    LIPS_output_file : str
        Path to file with LIPS output result.
    LIPS_parsed_csv : str
        Path to csv with LIPS output organised into features for machine learning.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """


    if os.path.isfile(LIPS_output_file):
        # try:
        surface_num = 0
        surface_lips = 100  ##100 is an initialized big number assuming lips score will not bigger than this number
        with open(LIPS_output_file, "r") as LIPS_output_handle:
            with open(LIPS_parsed_csv, "w") as LIPS_parsed_csv_handle:
                i = 0
                array = []
                dict = {}
                for row in LIPS_output_handle:
                    if re.search(r"^\s+\d+\s+[A-Z]", row) or re.search(r"^\d+\s+[A-Z]", row):
                        array = row.split()
                        if not int(array[0]) in dict:
                            dict[int(array[0])] = " ".join([array[1], array[2], array[3]])

                    if re.search("^\d{1}\s+", row):
                        surface_num1 = row.split()[0]
                        surface_lips1 = row.split()[3]
                        if (float(surface_lips1) < float(surface_lips)):
                            surface_lips = surface_lips1
                            surface_num = surface_num1
                LIPS_output_handle.close()

                surface_find = 0
                dict1 = {}
                LIPS_output_handle = open(LIPS_output_file, "r")
                for row in LIPS_output_handle:
                    if re.search(r"^SURFACE\s" + surface_num, row):
                        surface_find = 1
                        continue
                    if surface_find == 1 and (re.search(r"^\s+\d+\s+[A-Z]", row)or re.search(r"^\d+\s+[A-Z]", row)):
                        array = row.split()
                        if not int(array[0]) in dict1:
                            dict1[int(array[0])] = " ".join([array[1], array[2], array[3]])
                    else:
                        surface_find = 0
                LIPS_output_handle.close()

                writer = csv.writer(LIPS_parsed_csv_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(["residue_num", "residue_name", "LIPS_polarity", "LIPS_entropy", "LIPS_surface"])
                print(len(dict.keys()))
                for k, v in sorted(dict.items()):
                    v1 = v.split()
                    v1.insert(0, k)
                    if k in dict1:
                        v1.insert(4, 1)
                    else:
                        v1.insert(4, 0)
                    csv_header_for_cons_lips_score_file = v1
                    writer.writerow(csv_header_for_cons_lips_score_file)
                LIPS_parsed_csv_handle.close()
                print(f"{LIPS_output_file.split('.')[0]}lips score parse finished")
                #logging.info('{} lips score parse finished ({})'.format(acc, LIPS_parsed_csv))
    else:
        print("LIPS_output_file not found.")
        #logging.warning("{} LIPS_output_file not found.")



## iterate over files 
output_dir="LIPS_Feature/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

directory="LIPS_output/"

for filename in os.listdir(directory):
    if filename.endswith(".lips"):
        lips_path = os.path.join(directory, filename)
        output_file=output_dir+filename.split(".")[0]+".csv"
        parse_LIPS_score(lips_path,output_file) 


