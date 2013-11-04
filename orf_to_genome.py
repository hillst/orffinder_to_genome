#!/usr/bin/env python
'''
transdecoder_to_genome -- Convert outputted transcript coordinate to the appropriate genome coordinates.

transdecoder_to_genome is a tool for converting transcript coordinates to appropriate genome coordinate.
Expects the genome gff3 as input, along with the file transcripts file.

@author:     Steven Hill

@copyright:  2013 Donald Danforth Plant Science Center. All rights reserved.

@contact:    shill@danforthcenter.org
'''

import sys
import os
import traceback

from optparse import OptionParser

__all__ = []
__version__ = 0.2
__date__ = '2013-11-2'
__updated__ = '2013-11-2'
__author__ = "Steven Hill"

PROFILE = 0


def main(argv=None):
    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.1"
    program_build_date = "%s" % __updated__

    program_version_string = '%%prog %s (%s)' % (program_version, program_build_date)
    program_license = "Copyright 2013 Donald Danforth Plant Science Center                                           \
                Licensed under the Apache License 2.0\nhttp://www.apache.org/licenses/LICENSE-2.0"

    if argv is None:
        argv = sys.argv[1:]
        # setup option parser
        parser = OptionParser(version=program_version_string,
                              epilog="Tool for converting ORF peptide file to a gff with genome coordinates",
                              description=program_license)
        parser.add_option("-r", "--referencegff", dest="referencegff", help="Reference gff to use <FILE>",
                          metavar="FILE")
        parser.add_option("-b", "--best_cand", dest="bestcandgff",
                          help="Transdecoder best_candidates.eclipsed_orfs_removed.gff3 generated file <FILE>",
                          metavar="FILE")
        parser.add_option("-o", "--output", dest="output", help="File name to write to. All output is in gff3 format.")

        # process options
        (opts, args) = parser.parse_args(argv)
        rg, bcg, rf = None, None, None
        if opts.referencegff:
            rg = opts.referencegff
        if opts.bestcandgff:
            bcg = opts.bestcandgff
        if opts.output:
            output = open(opts.output, "w")
        else:
            output = sys.stdout
        if rg == None or bcg == None:
            print "Expects three arguments."
            parser.print_help()
            return 2

        convert_orffinder_to_genome_objects(rg, bcg, output)


"""
    Accepts a gff file containg atleast exons and genes in genome coordinates, a fasta file producded by OrfFinder
    containing local gene coordinates, and converts the gene coordinates into genome coordinates. It then prints out a
    corresponding GFF file containing these coordinates as CDS files.

    @param original_gff GFF file containing atleast genes and exons in genome coordinates. All other entries are ignored.
    @param best_cands   Fasta file from OrfFinder containing local coordinates.
    @param output       An open file handle which all results will be written to.

"""


def convert_orffinder_to_genome_objects(original_gff, best_cands, output):
    bc_lookup = {}
    #goes through the orf file,finds all the lines containing reading frame info, creates an Orf object representing
    # them, places them in the bc_lookup dictionary
    with open(best_cands, "r") as bc:
        for line in bc:
            if ">" in line:
                gene = Orf(line)
                bc_lookup[gene.parent] = gene

    with open(original_gff, "r") as og:
        last_gene = None
        for line in og:
            try:
                cur_gene = GFFLine(line)
            except Exception as e:
                #catches UTRs and mRNAs
                continue
            if cur_gene.getLabel() == "exon":
                last_gene.addExon(Exon(cur_gene.lines[1]))
            elif "gene" == cur_gene.getLabel():
                if last_gene != None and len(last_gene.exons) > 0:
                    print_buffer = []
                    cds = bc_lookup[last_gene.name]
                    if cds.strand is "-":
                        last_gene.switchStrand()
                        for exon in last_gene.exons:
                            exon.switchStrand()

                    print >> output, str(last_gene).strip()
                    print >> output, str(last_gene.as_mRNA()).strip()
                    try:
                        exons = {'start': last_gene.start, 'stop': last_gene.stop, 'name': last_gene.getName(), 'exons': last_gene.getExonDicList() }
                        #if we are on negative strand
                        if last_gene.strand == "-":
                            cdses = build_coordinates(exons, cds.getDicFormat(), True)
                        else:
                            cdses = build_coordinates(exons, cds.getDicFormat())
                        if len(cds) == 0:
                            raise Exception("CDS empty")
                        for exon in last_gene.exons:
                            print >> output, str(exon).strip()

                        for seq in cdses:
                            #sort of a hack and might not do what we want.
                            for exon in last_gene.exons:
                                if seq['start'] >= exon.start and seq['stop'] <= exon.stop:
                                    #change exon to CDS, put the CDS start at the start position and CDS end at end position.
                                    #DO IN TOKENIZED MANNER
                                    tokenized = exon.lines[0]
                                    tokenized[2] = "CDS"
                                    tokenized[3] = str(seq['start'])
                                    tokenized[4] = str(seq['stop'])
                                    tokenized[6] = exon.strand
                                    tokenized[-1] = tokenized[-1].replace("exon", "cds")
                                    fin = "\t".join(tok for tok in tokenized)
                                    print >> output, str(fin).strip()
                                    break

                    except Exception as e:
                        print >> sys.stderr, cur_gene
                        print >> sys.stderr, e
                last_gene = Transcript(cur_gene.lines[1])



"""
    This function should be responsible for building the genome coordinates and the exon coordinates and performing
    a translation from exon coordinates to genome coordinates.

    Exceptions will be raised if the names do not match or if one of the transcript coordinates is invalid. The definition
    of coordinates expected are to start at 1, and to go UNTIL the last number. So a reading frame of length one would have the
    same start as end.

    @param gene = {'start': '100', 'stop': '400', 'name': 'gene1', 'exons': ({'start': 100, 'stop': 149}, {'start': 200, 'end': 249}, {'start': 300, 'stop': 400})}
    @param transcript = {genename: (start, stop)}
    @param reverse  If reverse is set, it will interpret all start/stop coordinates as their inverse. That is start will represent
    x coordinates from the end of the genome and stop will represent y coordinates from the beginning.

    @returns gene_genome = <cds> ({"start":start, "stop":stop}, {"start":start,"stop":stop}, {"start":start,"stop":stop}... )
"""


def build_coordinates(gene, transcript, reverse=False):
    #either transcript or genes could be backwards, it's likely transcript
    if transcript['start'] < 1:
        raise Exception("Transcript Start is less than 1.")
    if transcript['stop'] < transcript['start']:
        raise Exception("Transcript Stop is less than start.")
    if (gene['name'] != transcript['name']):
        err = "Gene names do not match. " + str(gene['name']) + " " + str(transcript['name'])
        raise Exception(err)
        #to build our list use range and add one to the end because range is exclusive and gff3 is inclusive
    genome_coords = []
    gene_list = []
    #adjust
    if not reverse:
        transcript['start'] -=  1
        transcript['stop'] -=  1
    for exons in gene['exons']:
        genome_coords += range(exons['start'], exons['stop'] + 1)
        gene_list.append(range(exons['start'], exons['stop'] + 1))
    if reverse:
        transcript['start'] *= -1
        transcript['stop'] *= -1
        tmp = transcript['start']
        transcript['start'] = transcript['stop']
        transcript['stop'] = tmp

    start, stop = None, None
    try:
        for genes in gene_list:
            try:
                end = genome_coords[transcript['stop']]
            except IndexError:
                end = genome_coords[-1]
            if genome_coords[transcript['start']] in genes:
                start = gene_list.index(genes)
                start_idx = genes.index(genome_coords[transcript['start']])
            if end in genes:
                stop = gene_list.index(genes)
                stop_idx = genes.index(end)
    except IndexError:
        tb = traceback.format_exc()
        print >> sys.stderr, tb
        if reverse:
            print >> sys.stderr, transcript, len(genome_coords), genome_coords[transcript['start']]
            sys.exit()
        else: print >> sys.stderr, transcript['stop'], len(genome_coords), transcript['stop'] - len(genome_coords), \
            ( transcript['stop'] ) % 3, transcript['name']
    cds = []
    for i in range(start, stop + 1):
        cur_cds = gene['exons'][i]
        if i == start:
            cur_cds['start'] = gene_list[start][start_idx]
        elif i == stop:
            cur_cds['stop'] = gene_list[stop][stop_idx]
        if start == stop: #case where they are the same
            cur_cds['stop'] = gene_list[stop][stop_idx]
        cds.append(cur_cds)

    return cds


"""
    Base class for a line in a gff3 file. Responsible for parsing the lines, building dictionary format, handling strandedness.
    Before we know what something is it will be instantiated as a GFFLine. Also contains the original line and a parsed version of
    that line.
"""


class GFFLine():
    start = ""
    stop = ""
    lines = None
    strand = None
    parent = None

    def __init__(self, line):
        self.ParseHeader(line)

    def __len__(self):
        if self.stop == None or self.start == None or self.stop - self.start < 1:
            raise Exception("Improperly initialized. Invalid start and end")
        return self.stop - self.start

    def __str__(self):
        return self.lines[1]

    def __contains__(self, key):
        if key in self.lines[1]:
            return True
        else:
            return False

    def getLabel(self):
        return self.lines[0][2]

    def as_mRNA(self):
        out = self.lines[0]
        out[-1] = out[-1].replace("gene", "mRNA")
        out[2] = out[2].replace("gene", "mRNA")
        p1 = out[-1].split(";")[0]
        p1 += "-mRNA-1"
        out[-1] = p1 + ";"+out[-1].split(";")[0]
        return  "\t".join(tok for tok in out)

    """
    Sets start, stop, lines (linearr, raw), parent (gene name), and strandedness
    """
    def ParseHeader(self, line):
        linearr = line.split()
        if len(linearr) < 8:
            raise Exception("Line does not contain valid format.\n" + line)
        self.lines = [linearr, line]
        self.start = int(linearr[3])
        self.stop = int(linearr[4])
        self.parent = linearr[0]
        if self.start >= self.stop:
            raise Exception("Error in input, start must be less than stop: " + " ".join(linearr[3:5]))
        if "+" in line:
            self.strand = "+"
        else:
            self.strand = "-"

    def getDicFormat(self):
        raise Exception("GFFLine should never call getDicFormat. Check your typing.")
        return {'name': self.getLabel(), 'start': self.start, 'stop': self.stop}

    """
    Returns self for convenience.
    Switched from + to - and - to +. Also updates the lines.
    """

    def switchStrand(self):
        if self.strand == "+":
            self.strand = "-"
            self.lines[0][6] = self.strand
            self.lines[1] = "\t".join([tok for tok in self.lines[0]])
            return self
        if self.strand == "-":
            self.strand = "+"
            self.lines[0][6] = self.strand
            self.lines[1] = "\t".join([tok for tok in self.lines[0]])
            return self


"""
    Child class of GFFLine that is for locus. Locus are special because they are the top of the heirarchy and may contain transcripts.
    This class is also responsible for printing all of the child transcripts and their associated exons/CDS
"""


class Locus(GFFLine):
    transcripts = None
    transcript_lines = None
    name = ""

    def __init__(self, line):
        GFFLine.__init__(self, line)
        self.transcript_lines = []
        self.transcripts = []
        self.name = self.lines[0][8].split("=")[1].split(";")[0]

    def addTranscriptLine(self, line):
        self.transcript_lines.append(line)

    def addTranscript(self, trans):
        self.transcripts.append(trans)

    def getName(self):
        return self.name

    def printLines(self):
        print self.lines[1],
        for trans in self.transcripts:
            if len(trans.CDS) != 0:
                print trans,
                for exon in trans.exons:
                    print exon,
                    for cds in trans.CDS:
                        if cds.start >= exon.start and cds.stop <= exon.stop:
                            print cds,
                            break

    def chooseBest(self):
        best = None
        for trans in self.transcripts:
            if trans > best:
                best = trans
        self.transcripts = [best]


"""
Contains a list of exons, and overridden operators for Greater Than, Less Than, Equals.
Contains various function on the set of exons in this transcript. Also includes a name function
"""


class Transcript(GFFLine):
    exons = None
    partial = None
    name = None
    CDS = None

    def __init__(self, line):
        GFFLine.__init__(self, line)
        self.exons = []
        self.CDS = []
        self.length = None
        if "partial" in self:
            self.partial = "partial"
        elif "internal" in self:
            self.partial = "internal"
        else:
            self.partial = "complete"
        self.name = self.lines[0][8].split("=")[1].split(";")[0]

    def setPartial(self, partial):
        self.partial = partial

    def __len__(self):
        if self.length != None:
            return self.length
        else:
            return GFFLine.__len__(self)

    def getName(self):
        return self.name

    def setLength(self, length):
        self.length = length

    def addExon(self, exon):
        if not isinstance(exon, Exon):
            raise Exception("Passed is not of type exon. Type Exon expected.")
        self.exons.append(exon)

    def addCDS(self, CDS):
        if not isinstance(CDS, Exon):
            raise Exception("Passed is not of type exon. Type Exon expected.")
        self.CDS.append(CDS)

    def getExonDicList(self):
        lst = []
        for exon in self.exons:
            lst.append(exon.getDicFormat())
        return lst

    def getDicFormat(self):
        return {'name': self.name, 'start': self.start, 'stop': self.stop}

    def switchDirections(self):
        self.switchStrand()
        new_exons = []
        for exon in self.exons:
            new_exons.append(exon.switchStrand())
        self.exons = new_exons
        return True

    def __gt__(self, b):
        if b == None:
            return True
        if self.isComplete() and not b.isComplete():
            return True
        if b.isComplete() and not self.isComplete():
            return False
        if len(self) > len(b):
            return True
        if len(b) > len(self):
            return False
        return False

    def __lt__(self, b):
        if b == None:
            return False
        if self.isComplete() and not b.isComplete():
            return False
        if b.isComplete() and not self.isComplete():
            return True
        if len(self) > len(b):
            return False
        if len(b) > len(self):
            return True
        return False

    def __eq__(self, b):
        #so brave
        if not self > b and not self < b:
            return True

    def isComplete(self):
        if "complete" in self.partial:
            return True
        else:
            return False


"""
Contains a list of exons, and overridden operators for Greater Than, Less Than, Equals.
Contains various function on the set of exons in this transcript. Also includes a name function
"""


class Orf(GFFLine):
    exons = None
    partial = None
    name = None
    CDS = None
    parent = None
    length = None

    def __init__(self, line):
        #>Spolyrhiza9509S001G0001        +1      1       1131
        line = line[1:]
        linearr = line.split()
        self.strand = linearr[1][0]
        self.start = int(linearr[2])
        self.stop = int(linearr[3])
        self.parent = linearr[0]
        self.name = linearr[0]
        self.lines = [linearr, line]
        self.exons = []
        self.CDS = []

    def __len__(self):
        if self.length != None:
            return self.length
        else:
            return GFFLine.__len__(self)

    def getName(self):
        return self.name

    def setLength(self, length):
        self.length = length

    def addExon(self, exon):
        if not isinstance(exon, Exon):
            raise Exception("Passed is not of type exon. Type Exon expected.")
        self.exons.append(exon)

    def addCDS(self, CDS):
        if not isinstance(CDS, Exon):
            raise Exception("Passed is not of type exon. Type Exon expected.")
        self.CDS.append(CDS)

    def getExonDicList(self):
        lst = []
        for exon in self.exons:
            lst.append(exon.getDicFormat())
        return lst

    def getDicFormat(self):
        return {'name': self.name, 'start': self.start, 'stop': self.stop}

    def switchDirections(self):
        if self.strand == "+":
            self.strand = "-"
            self.lines[0][2] = self.lines[0][2].replace("+", self.strand)
            self.lines[1] = self.lines[1].replace("+", self.strand)
        if self.strand == "-":
            self.strand = "+"
            self.lines[0][2] = self.lines[0][2].replace("-", self.strand)
            self.lines[1] = self.lines[1].replace("-", self.strand)
        new_exons = []
        for exon in self.exons:
            new_exons.append(exon.switchStrand())
        self.exons = new_exons
        return True


"""
Overrides getDicFormat to use the parent as the name
"""

class Exon(GFFLine):
    def __init__(self, line):
        GFFLine.__init__(self, line)

    def parseName(self):
        self.name = self.lines[0][8].split("=")[1].split(";")[0]

    def getDicFormat(self, name=None):
        if name == None:
            self.parseName()
        else:
            self.name = name
        return {'name': self.name, 'start': self.start, 'stop': self.stop}

if __name__ == "__main__":
    main()
