#!/usr/bin/env python3

import os
import sys

def usage():
    if len(sys.argv) != 4:
        print('usage: ', os.path.basename(__file__), ' [in.fasta] [out.fasta] [length threshold]') 
        sys.exit(1)
    return

def main():
    usage()

    transcript_seq = ''
    transcript_length = 0
    
    shortest_length = int(sys.argv[3])
    out_transcript = open(os.path.abspath(sys.argv[2]), 'w')

    with open(os.path.abspath(sys.argv[1]), 'r') as in_transcript:
        for i, transcript_line in enumerate(in_transcript):
            if transcript_line[0] == '>':
                if i != 0 and transcript_length >= shortest_length:
                    print(transcript_seq, file=out_transcript)
                transcript_seq = ''
                transcript_length = 0
            else:
                transcript_length += len(transcript_line.rstrip())

            transcript_seq += transcript_line
    
    out_transcript.close()

    return

if __name__ == '__main__':
    main()
    sys.exit(0)
