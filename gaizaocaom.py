#!/usr/bin/env python3

import mappy as mp


import argparse


def parseArgs():
    parser = argparse.ArgumentParser(
        description='compare assemblies with minimap2')
    parser.add_argument('-r', help='input assembly 1 (ref) fasta',
                        dest='asm1Filename', required=True)
    parser.add_argument('-q', help='input assembly 2 (query) fasta',
                        dest='asm2Filename', required=True)

    parser.add_argument('-n', help='minimum sequence length in assembly 2 (query), default=20Mb',
                        dest='minQueryLen', required=False, default=20000000, type=int)
    parser.add_argument('-o', help='output prefix, default=output',
                        dest='outputPrefix', required=False, default='output')
    arguments = parser.parse_args()
    print("arguments:")
    print("-r, input assembly 1 (ref) fasta: %s" % arguments.asm1Filename)
    print("-q, input assembly 2 (query) fasta: %s" % arguments.asm2Filename)
    print("-n, minimum sequence length in assembly 2: %i" %
          arguments.minQueryLen)
    print("-o, output prefix: %s" % arguments.outputPrefix)

    print("-- \n")

    return arguments


def runIndex(refFastaFilname):
    print("building index for assembly 1 (ref): %s\n" % refFastaFilname)
    a = mp.Aligner(refFastaFilname, preset="asm5")
    return a


def getTopHitByAlignmentLength(hitsList, queryID):
    "function getTopHitByAlignmentLength"
    dic = {}
    for i in hitsList:
        
        i = i.split("\t")
        key = i[2]
        va = int(i[3])
        if key not in dic:
            dic[key] = va
        else:
            dic[key] += va
        max_len = max(zip(dic.values(), dic.keys()))

    output = {'top_aln_id': max_len[1], 'top_aln_blen': max_len[0]}
    return output


def runMapper(referenceIndex, asm2Filename, minQueryLen):
    print("running minimap2 and finding top hit per query sequence\n")
    scaffoldMapList0 = []

    for name, seq, qual in mp.fastx_read(asm2Filename):
        print("... query: %s" % name)
        if len(seq) < minQueryLen:
            print("...... Skipping, query too short (seq len of %i is less than minimum: %i)\n" % (
                len(seq), minQueryLen))
            continue
        hits = []
        for hit in referenceIndex.map(seq):

            hits.append(name+"\t"+str(len(seq))+"\t" +
                        hit.ctg+"\t"+str(hit.mlen))

        topAln = getTopHitByAlignmentLength(hits, name)
        print("Top hit: %s\n" % topAln['top_aln_id'])
        scaffoldMapList0.append({'queryID': name, 'qury_len': len(seq),
                                 'refID': topAln['top_aln_id'],
                                 'alignLen': topAln['top_aln_blen'],
                                 })
    return (scaffoldMapList0)


def writeOutput(outputPrefix, scaffoldMapListOut):
    print("writing results")

    with open(outputPrefix + "_top_hit_results", 'w') as r:
        for i in scaffoldMapListOut:
            r.write(i['queryID']+"\t"+str(i['qury_len'])+"\t" +
                    i['refID']+"\t"+str(i['alignLen'])+"\n")


def run():
    "run script"
    # parse arguments
    args = parseArgs()

    # generate index for reference fasta
    refIndex = runIndex(refFastaFilname=args.asm1Filename)

    scaffoldMapList = runMapper(
        referenceIndex=refIndex, asm2Filename=args.asm2Filename, minQueryLen=args.minQueryLen)

    # write output
    writeOutput(outputPrefix=args.outputPrefix,
                scaffoldMapListOut=scaffoldMapList)

    print("DONE")


if __name__ == '__main__':
    run()
