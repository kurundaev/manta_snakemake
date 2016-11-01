




# Insertion
# Inversion
# Deletion
# BND
# Duplication


import re
from collections import defaultdict

import sys


def find_mutation_type(col_info):
    """

    :param col_info: string containing info of mutation type
    :return: string mutation type
    """
    matched =  re.search('SVTYPE=([\w]+)', col_info)
    mutation_type = matched.group(1)
    return mutation_type

def find_sv_len(col_info):
    """

    :param col_info: string containg info of SVLEN
    :return:  return string svlen
    """
    matched = re.search('SVLEN=([\d|-]+)', col_info)
    sv_len = matched.group(1)
    return sv_len

def find_mate(col_info):
    """

    :param col_info: string containing info about BND
    :return:  mate id
    """
    matched = re.search('MATEID=([^;]+)', col_info)
    mate_id = matched.group(1)
    return mate_id

def main(file):
    """

    :param file: input vcf
    :return: : None, write files to output
    """

    # with statment calls closea at the end. open file as f. (file object )

    mutation_profile = defaultdict(dict)
    bnd_mutation = defaultdict(dict)

    with open(file, 'r') as f:
        # read each line by line
        for line in f:
            if not line.startswith('#'):

                
                # remove the new line character
                line = line.strip()
                # split by tab character this will return a list
                fields = line.split('\t')
                
                col_info = fields[7]
                chr = fields[0].replace('chr','hs')
                if "_" not in chr: 
                    start = fields[1]
                    mutation_type = find_mutation_type(col_info)


                    if mutation_type == "BND":
                        id = fields[2]
                        mate_id = find_mate(col_info)

                        if mate_id in bnd_mutation:
                            bnd_mutation[mate_id]['mate_position'] = "{chr}\t{start}\tCTX".format(chr=chr, start=start)
                        else:
                            bnd_mutation[id]['mate'] = mate_id
                            bnd_mutation[id]['position'] = "{chr}\t{start}".format(chr=chr, start=start)

                    else:
                        end = start
                        if 'SVLEN' in col_info:
                            sv_len = find_sv_len(col_info)
                            end = abs(int(sv_len)) + int(start)
                        print("{chr}\t{start}\t{chr}\t{end}\t{mut_type}".format(chr=chr,
                                                                           start=start,
                                                                           end = end,
                                                                           mut_type = mutation_type
                                                                           ))

    for key,value in bnd_mutation.items():
        try :  
            print("{start}\t{end}".format(start= value['position'], end=value['mate_position']))
        except KeyError:
            pass
     

if __name__ == '__main__':
    main(sys.argv[1])







