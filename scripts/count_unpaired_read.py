#!/usr/bin/env python3

import os
import sys 
import re
import pickle

def usage():
    if len(sys.argv) != 4:
        print('usage:', os.path.basename(__file__), '[read1.fastq] [read2.fastq] [out.pickle]')
        sys.exit(1)
    return

def main():
    usage()
       
    count = dict()
    for i in range(1, 3):
        with open(sys.argv[i], 'r') as fastq:
            for i, data in enumerate(fastq):
                if i % 4 != 0:
                    continue
                regex = re.match('@\S+?:\S+?:\S+?:(\S+?):', data)
                if regex:
                    name = regex.group(1)
                    try:
                        count[name] += 1 
                    except:
                        count[name] = 1

    with open(sys.argv[3], 'wb') as out:
        pickle.dump(count, out, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    main()
