#!/usr/bin/env python3

import os
import sys

def usage():
    if len(sys.argv) != 2:
        print('usage:', os.path.basename(__file__), '[flux_simulator.fastq]')
        sys.exit(1)
    return

def main():
    
    usage()
   
    fastq_path = os.path.abspath(sys.argv[1])
    outdir = os.path.dirname(fastq_path)
    fastq_file = open(fastq_path, 'r')
    read_1 = open(outdir + '/flux_simulator_r1.fastq', 'w')
    read_2 = open(outdir + '/flux_simulator_r2.fastq', 'w')
    
    fragment_count = 0
    
    for i, line in enumerate(fastq_file):
        if i % 8 == 0:
            fragment_count += 1
        
        if i % 4 == 0: 
            read_name = line.strip()
            if read_name[-1] == '1':
                pair = '/1'
            else:
                pair = '/2'
            output = '@flux_simulator_' + str(fragment_count) + ':' + read_name[1:-4] + pair + '\n'
        else:
            output = line

        if pair == '/1':
            print(output, file=read_1, end='')
        else:
            print(output, file=read_2, end='')

    fastq_file.close()
    read_1.close()
    read_2.close()
    
    return

if __name__ == '__main__':
    main()
    sys.exit(0)
