
import os, sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/cava/pysamdir')
import pysam

############################################################################################################

# Class representing a single CSN annotation
class CSN(object):

    # Constructor
    def __init__(self, CSNstring, transcriptID):
        self.CSNstring = CSNstring
        self.transcriptID = transcriptID
        self.coord1 = None
        self.intr_coord1 = None
        self.coord2 = None
        self.intr_coord2 = None
        self.dna = None
        self.protein = None
        self.error = None
        self.vcf = None
        self.transcript = None
        self._deleted = None
        self._duplicated = None

    # Return type of variant described by the CSN
    def getType(self):
        if 'delins' in self.dna: return 'complex'
        elif 'del' in self.dna: return 'deletion'
        elif 'ins' in self.dna: return 'insertion'
        elif 'dup' in self.dna: return 'insertion'
        elif '>' in self.dna: return 'substitution'

    # Parse string
    def parseString(self, s):

        # Split to c. and p. parts
        if '_p.' in s: csn_c, csn_p = s.split('_p.', 1)
        else: csn_c, csn_p = s, ''

        # Check if c part starts with 'c.'
        if not csn_c.startswith('c.'): self._giveError(1)

        # Find where to split the c. part to get coordinates and dna change
        found = [csn_c.find(x) for x in ['A>', 'C>', 'G>', 'T>', 'del', 'ins', 'dup'] if x in csn_c]
        if len(found) == 0: self._giveError(2)
        cuthere = min(found)

        # Split c. part to get coordinates and dna change
        coords = csn_c[csn_c.find("c.") + 2:cuthere]
        dna = csn_c[cuthere:]

        # Split coordinates to first and second coordinate
        if "_" in coords: [coord_part1, coord_part2] = coords.split('_', 1)
        else: coord_part1, coord_part2 = coords, None

        # Split first coordinate to exon and intron part (if any)
        exonpart1, intronpart1 = coord_part1, None
        for i in range(1,len(coord_part1)):
            if coord_part1[i] in ['-', '+']:
                exonpart1 = coord_part1[:i]
                intronpart1 = coord_part1[i:]
                break

        # Check if exonpart1 and intronpart1 are integers
        try: int(exonpart1)
        except: self._giveError(3, exonpart1)
        if intronpart1 is not None:
            try: int(intronpart1)
            except: self._giveError(4, intronpart1)

        # Split second coordinate to exon and intron part (if any)
        if coord_part2 is not None:
            exonpart2, intronpart2 = coord_part2, None
            for i in range(1,len(coord_part2)):
                if coord_part2[i] in ['-', '+']:
                    exonpart2 = coord_part2[:i]
                    intronpart2 = coord_part2[i:]
                    break

            # Check if exonpart2 and intronpart2 are integers
            try: int(exonpart2)
            except: self._giveError(3, exonpart2)
            if intronpart2 is not None:
                try: int(intronpart2)
                except: self._giveError(4, intronpart2)

        else:
            exonpart2 = None
            intronpart2 = None

        # Check if dna change is well formatted

        # Check syntax of base substitution
        if dna[0] in ['A', 'C', 'G', 'T']:
            if len(dna) != 3 or dna[1]!='>' or dna[2] not in ['A', 'C', 'G', 'T']: self._giveError(5, dna)

        # Check syntax of delins
        elif dna.startswith('delins'):
            inserted = dna[6:]
            if not all([c in ['A','C','G','T'] for c in inserted]): self._giveError(6, inserted)

        # Check syntax of del
        elif dna.startswith('del'):
            deleted = dna[3:]
            self._deleted = deleted
            try:
                delN = int(deleted)
                if delN < 5: self._giveError(7)
            except ValueError:
                if not all([c in ['A','C','G','T'] for c in deleted]): self._giveError(8, deleted)
                if len(deleted) >= 5: self._giveError(9)

        # Check syntax of ins
        elif dna.startswith('ins'):
            inserted = dna[3:]
            if not all([c in ['A','C','G','T'] for c in inserted]): self._giveError(6, inserted)
            if exonpart2 is None: self._giveError(10)

        # Check syntax of dup
        elif dna.startswith('dup'):
            duplicated = dna[3:]
            self._duplicated = duplicated
            try:
                dupN = int(duplicated)
                if dupN < 5: self._giveError(11)
            except ValueError:
                if not all([c in ['A','C','G','T'] for c in duplicated]): self._giveError(12, duplicated)
                if len(duplicated) >= 5: self._giveError(13)

        '''
        # Check if p. part is well formatted
        if csn_p != '':

            cdn_p_ok = False
            if csn_p == '=':
                cdn_p_ok = True
            else:
                if self._isAminoAcid(csn_p[:3]) and self._isAminoAcid(csn_p[-3:]) and self._isPosition(csn_p[3:-3]): cdn_p_ok = True
                elif csn_p[-3:] == 'del':
                    coords = csn_p[:-3]
                    if '_' in coords:
                        first = coords[:coords.find('_')]
                        second = coords[coords.find('_')+1:]
                        if self._isAminoAcid(first[:3]) and self._isPosition(first[3:]) and self._isAminoAcid(second[:3]) and self._isPosition(second[3:]): cdn_p_ok = True
                    else:
                        if self._isAminoAcid(coords[:3]) and self._isPosition(coords[3:]): cdn_p_ok = True
                elif 'delins' in csn_p:
                    [coords, inserted] = csn_p.split('delins')
                    coords_ok = False
                    inserted_ok = False

                    if '_' in coords:
                        first = coords[:coords.find('_')]
                        second = coords[coords.find('_')+1:]
                        if self._isAminoAcid(first[:3]) and self._isPosition(first[3:]) and self._isAminoAcid(second[:3]) and self._isPosition(second[3:]): coords_ok = True
                    else:
                        if self._isAminoAcid(coords[:3]) and self._isPosition(coords[3:]): coords_ok = True

            #print 'xyz   ' + csn_p + '  ' + str(cdn_p_ok)
        '''

        # Set fields of CSN object
        self.coord1 = exonpart1
        self.intr_coord1 = intronpart1
        self.coord2 = exonpart2
        self.intr_coord2 = intronpart2
        self.dna = dna
        self.protein = csn_p

    # Convert to VCF (calculate self.vcf)
    def convert(self, reference, transcript):

        # Calculate first genomic position
        try:
            genomic1 = self._calculateVCFCoordinate(transcript, self.coord1, self.intr_coord1)
        except:
            coord = str(self.coord1)
            if self.intr_coord1 is not None: coord += str(self.intr_coord1)
            self._giveError(14, coord)

        # Calculate second genomic position
        if self.coord2 is not None:
            try:
                genomic2 = self._calculateVCFCoordinate(transcript, self.coord2, self.intr_coord2)
            except:
                coord = str(self.coord2)
                if self.intr_coord2 is not None: coord += str(self.intr_coord2)
                self._giveError(14, coord)
        else:
            genomic2 = None

        # Check coordinate order
        if genomic2 is not None:
            if transcript.strand == 1:
                if genomic1 > genomic2: self._giveError(24)
            else:
                if genomic1 < genomic2: self._giveError(24)

        # Calculate REF and ALT
        try:
            pos, ref, alt = self._makeREFandALTstrings(transcript, reference, genomic1, genomic2)
        except:
            self._giveError(15)

        # Create Variant object
        variant = Variant(transcript.chrom, pos, ref, alt)

        # Check if substitutions are described by only a single position
        if genomic2 is not None and self.getType() == 'substitution': self._giveError(23)

        # Check if variant type agrees with CSN variant type
        if variant.getType() != self.getType(): self._giveError(16, variant.getType(), self.getType())

        # Check if reference base agrees with subsituted base in CSN
        if variant.getType() == 'substitution':
            if transcript.strand == 1:
                if variant.ref != self.dna[0]: self._giveError(21, variant.ref, self.dna[0])
            else:
                if variant.ref != Sequence(self.dna[0]).reverseComplement(): self._giveError(21, variant.ref, Sequence(self.dna[0]).reverseComplement())

        # Check if coordinates refer to neighbouring bases for insertions
        if self.getType() == 'insertion' and 'dup' not in self.dna:
            if abs(genomic2-genomic1) != 1: self._giveError(22)

        # Check if deleted reference bases agree with deleted bases in CSN
        if self._deleted is not None:
            try:
                n = int(self._deleted)
                if len(ref) - 1 != n: self._giveError(17, str(len(ref) - 1), str(n))
            except ValueError:
                delseq = Sequence(self._deleted)
                if transcript.strand == -1: delseq = delseq.reverseComplement()
                if delseq != ref[1:]: self._giveError(18, ref[1:], delseq)

        # Check if duplicated reference bases agree with duplicated bases in CSN
        if self._duplicated is not None:
            try:
                n = int(self._duplicated)
                if len(alt) - 1 != n: self._giveError(19, str(len(alt) - 1), str(n))
            except ValueError:
                dupseq = Sequence(self._duplicated)
                if transcript.strand == -1: dupseq = dupseq.reverseComplement()
                if dupseq != alt[1:]: self._giveError(20, alt[1:], dupseq)

        # Left-align variant
        if not (len(ref) == 1 and len(alt) == 1): variant.leftAlignVariant(reference)

        # Set self.vcf variable
        self.vcf = VCFRecord()
        self.vcf.addVariant(variant)
        self.transcript = transcript.ENST

    # Convert CSN coordinates to VCF coordinate
    def _calculateVCFCoordinate(self, transcript, exoncoord, introncoord):

        if exoncoord[0] == '+': Q = transcript.codingEnd + int(exoncoord)
        elif exoncoord[0] == '-': Q = int(exoncoord) + 1
        else: Q = int(exoncoord)

        if introncoord is None: introncoord = 0
        else: introncoord = int(introncoord)

        summed = -int(transcript.codingStart)+1

        for exon in transcript.exons:
            summed += exon.length
            w = summed - Q
            if introncoord == 0:
                if w >= 0:
                    if transcript.strand == -1: endpos = exon.start+1
                    else: endpos = exon.end
                    return endpos - transcript.strand*w
            elif introncoord > 0:
                if w == 0:
                    if transcript.strand == -1: endpos = exon.start+1
                    else: endpos = exon.end
                    return endpos+transcript.strand*introncoord
            else:
                if w == exon.length - 1:
                    if transcript.strand==-1: endpos = exon.end
                    else: endpos = exon.start+1
                    return endpos+transcript.strand*introncoord

        raise CSNException('')

    # Create POS, REF and ALT fields of VCF
    def _makeREFandALTstrings(self, transcript, reference, genomic1, genomic2):

        if genomic2 == None: genomic2 = genomic1

        if transcript.strand == 1:
            if "del" in self.dna:
                if "ins" in self.dna:
                    ref = reference.getSequence(transcript.chrom, genomic1, genomic2)
                    alt = self.dna[self.dna.find("ins")+3:]
                    pos = genomic1
                else:
                    ref = reference.getSequence(transcript.chrom, genomic1-1, genomic2)
                    alt = reference.getSequence(transcript.chrom, genomic1-1, genomic1-1)
                    pos = genomic1 - 1
            elif "ins" in self.dna:
                ref = reference.getSequence(transcript.chrom, genomic1, genomic1)
                alt = ref+self.dna[self.dna.find("ins")+3:]
                pos = genomic1
            elif "dup" in self.dna:
                ref = reference.getSequence(transcript.chrom, genomic2, genomic2)
                alt = ref+reference.getSequence(transcript.chrom, genomic1, genomic2)
                pos = genomic2
            elif ">" in self.dna:
                ref = reference.getSequence(transcript.chrom, genomic1, genomic1)
                alt = self.dna[2]
                pos = genomic1
            return pos, ref, alt

        else:
            if "del" in self.dna:
                if "ins" in self.dna:
                    ref = reference.getSequence(transcript.chrom, genomic2, genomic1)
                    alt = Sequence(self.dna[self.dna.find("ins")+3:]).reverseComplement()
                    pos = genomic2
                else:
                    ref = reference.getSequence(transcript.chrom, genomic2-1, genomic1)
                    alt = reference.getSequence(transcript.chrom, genomic2-1, genomic2-1)
                    pos = genomic2 - 1
            elif "ins" in self.dna:
                ref = reference.getSequence(transcript.chrom, genomic2, genomic2)
                alt = ref+Sequence(self.dna[self.dna.find("ins")+3:]).reverseComplement()
                pos = genomic2
            elif "dup" in self.dna:
                ref = reference.getSequence(transcript.chrom, genomic1, genomic1)
                alt = ref+reference.getSequence(transcript.chrom, genomic2, genomic1)
                pos = genomic1
            elif ">" in self.dna:
                ref = reference.getSequence(transcript.chrom, genomic1, genomic1)
                alt = Sequence(self.dna[2]).reverseComplement()
                pos = genomic1
            return pos, ref, alt

    # Raise various exceptions (CSNException) if anything goes wrong in parsing or converting
    def _giveError(self, errortype, *arg):
        if errortype == 1: raise CSNException('CSN should start with the c. prefix')
        if errortype == 2: raise CSNException('DNA change (del, ins, dup or base substitution) incorrect/missing')
        if errortype == 3: raise CSNException('Transcript coordinate (%s) incorrect, should be integer' % arg[0])
        if errortype == 4: raise CSNException('Intron coordinate (%s) incorrect, should be integer' % arg[0])
        if errortype == 5: raise CSNException('Description of base substitution (%s) incorrect' % arg[0])
        if errortype == 6: raise CSNException('Inserted sequence (%s) incorrect' % arg[0])
        if errortype == 7: raise CSNException('For deletions of less than 5 bases, deleted bases should be included')
        if errortype == 8: raise CSNException('Deleted sequence (%s) incorrect' % arg[0])
        if errortype == 9: raise CSNException('For deletions of 5+ bases, only the number of deleted bases required')
        if errortype == 10: raise CSNException('Insertion should be described by coordinates of both flanking bases')
        if errortype == 11: raise CSNException('For duplications of less than 5 bases, duplicated bases should be included')
        if errortype == 12: raise CSNException('Duplicated sequence (%s) incorrect' % arg[0])
        if errortype == 13: raise CSNException('For duplications of 5+ bases, only the number of duplicated bases required')
        if errortype == 14: raise CSNException('CSN coordinates (%s) cannot be converted into genomic coordinates' % arg[0])
        if errortype == 15: raise CSNException('VCF alleles cannot be determined')
        if errortype == 16: raise CSNException('Variant types described by coordinates (%s) and CSN (%s) disagree' % (arg[0], arg[1]))
        if errortype == 17: raise CSNException('Number of deleted reference bases (%s) and deleted bases in CSN (%s) disagree' % (arg[0], arg[1]))
        if errortype == 18: raise CSNException('Deleted reference bases (%s) and deleted bases in CSN (%s) disagree' % (arg[0], arg[1]))
        if errortype == 19: raise CSNException('Number of inserted bases (%s) and duplicated bases in CSN (%s) disagree' % (arg[0], arg[1]))
        if errortype == 20: raise CSNException('Reference bases (%s) and duplicated bases in CSN (%s) disagree' % (arg[0], arg[1]))
        if errortype == 21: raise CSNException('Reference base (%s) and substituted base in CSN (%s) disagree' % (arg[0], arg[1]))
        if errortype == 22: raise CSNException('Insertions should be described by two neighbouring bases')
        if errortype == 23: raise CSNException('Substitutions should be described by only a single genomic position')
        if errortype == 24: raise CSNException('Transcript coordinates are in wrong order')


    # Check if string describes an amino acid or stop codon
    def _isAminoAcid(self, s):
        aas = ['Ile', 'Met', 'Thr', 'Asn', 'Lys', 'Ser', 'Arg', 'Leu', 'Pro', 'His', 'Gln', 'Val', 'Ala', 'Asp', 'Glu', 'Gly', 'Phe', 'Tyr', 'Cys', 'Trp', 'X']
        return s in aas

    # Check if string describes a position
    def _isPosition(self, s):
        try:
            si = int(s)
            return si > 0
        except ValueError:
            return False

    # String representation
    def __str__(self):
        if self.coord1 is None: return ''
        ret = 'c.' + str(self.coord1)

        if self.intr_coord1 is not None:
            if self.intr_coord1 > 0: ret += '+' + str(self.intr_coord1)
            else: ret += str(self.intr_coord1)

        if self.coord2 is not None:
            ret += '_' + str(self.coord2)
            if self.intr_coord2 is not None:
                if self.intr_coord2 > 0: ret += '+' + str(self.intr_coord2)
                else: ret += str(self.intr_coord2)

        return ret + self.dna + self.protein

# Exception class for reporting CSN parsing or conversion errors
class CSNException(Exception): pass

############################################################################################################

# Class representing a single variant
class Variant(object):

    # Constructor
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt

    # Return type of variant
    def getType(self):
        if len(self.ref) == 1 and len(self.alt) == 1: return 'substitution'
        elif len(self.ref) == 1 and len(self.alt) > 1 and self.ref == self.alt[0]: return 'insertion'
        elif len(self.ref) > 1 and len(self.alt) == 1 and self.ref[0] == self.alt: return 'deletion'
        elif self.ref[0] != self.alt[0] and not (len(self.ref) == 1 and len(self.alt) == 1): return 'complex'

    # Left-align variant
    def leftAlignVariant(self, reference):
        seq1 = reference.getSequence(self.chrom, self.pos - 100, self.pos + len(self.ref) - 1)
        s = reference.getSequence(self.chrom, self.pos - 100, self.pos - 1)
        seq2 = s + self.alt
        N = len(s)
        left, seq1, seq2 = self._leftAlign(seq1, seq2)
        if len(seq1) == 0 or len(seq2) == 0:
            left -= 1
            base = reference.getSequence(self.chrom, self.pos + left - N, self.pos + left - N)
            seq1, seq2 = base + seq1, base + seq2
        self.pos = self.pos + left - N
        self.ref = seq1
        self.alt = seq2

    # Left-align two sequences
    def _leftAlign(self, seq1, seq2):
        right, seq1, seq2 = self._trimCommonEnd(seq1, seq2)
        left, seq1, seq2 = self._trimCommonStart(seq1, seq2)
        return left, seq1, seq2

    # Trim common starting subsequence of two sequences
    def _trimCommonStart(self, s1, s2):
        counter = 0
        while True:
            if len(s1) == 0 or len(s2) == 0: return counter, s1, s2
            if s1[0] != s2[0]: return counter, s1, s2
            s1, s2 = s1[1:], s2[1:]
            counter += 1

    # Trim common ending subsequence of two sequences
    def _trimCommonEnd(self, s1, s2):
        counter = 0
        while True:
            if len(s1) == 0 or len(s2) == 0: return counter, s1, s2
            if s1[-1] != s2[-1]: return counter, s1, s2
            s1, s2 = s1[:-1], s2[:-1]
            counter += 1

    # String representation
    def __str__(self):
        return self.chrom+':'+str(self.pos)+'_'+self.ref+'_'+self.alt

############################################################################################################

# Class representing a single VCF record
class VCFRecord(object):

    # Constructor
    def __init__(self):
        self.chrom = None
        self.pos = None
        self.id = '.'
        self.ref = None
        self.alts = []
        self.qual = '.'
        self.filter = '.'
        self.info = '.'

    # Add new variant to the record
    def addVariant(self, variant):
        if self.chrom is None:
            self.chrom = variant.chrom
            self.pos = variant.pos
            self.ref = variant.ref
            self.alts.append(variant.alt)
            return
        if variant.chrom != self.chrom: return
        if variant.pos != self.pos: return
        if variant.ref != self.ref: return
        self.alts.append(variant.alt)

    # Get variant by index
    def getVariant(self, i):
        return Variant(self.chrom, self.pos, self.ref, self.alts[i])

    # String representation
    def __str__(self):
        if self.chrom is None or self.pos is None or self.ref is None or len(self.alts) == 0: return ''
        record = [str(self.chrom), str(self.pos), str(self.id), self.ref, ','.join(self.alts), str(self.qual), self.filter, self.info]
        return '\t'.join(record)

############################################################################################################

# Class representing a single transcript
class Transcript(object):

    # Constructor
    def __init__(self, line):
        self.exons = []
        cols = line.split('\t')
        self.ENST = cols[0]
        self.geneSymbol = cols[1]
        self.geneID = cols[2]
        self.TRINFO = cols[3]
        self.chrom = cols[4]
        self.strand = int(cols[5])
        self.transcriptStart = int(cols[6])
        self.transcriptEnd = int(cols[7])
        self.codingStart = int(cols[8])
        self.codingStartGenomic = int(cols[9])
        self.codingEndGenomic = int(cols[10])

        # Initializing and adding exons
        for i in range(1, len(cols) - 11, 2):
            self.exons.append(Exon(int((i + 1) / 2), int(cols[10 + i]), int(cols[11 + i])))

        self.codingEnd = self.codingEnd()

    # Determine codingEnd position
    def codingEnd(self):
        summed = -int(self.codingStart)+1
        for exon in self.exons:
            summed += exon.length
            if exon.contains(self.codingEndGenomic):
                if self.strand == 1: return summed - (exon.end - self.codingEndGenomic)
                else: return summed - (self.codingEndGenomic - (exon.start + 1) )

    # Getting the full coding sequence of the transcript
    def getCodingSequence(self, reference, variant):
        ret = ''
        for i in range(len(self.exons)):
            exon = self.exons[i]
            if not variant is None:
                if exon.start < variant.pos <= exon.end:
                    if self.strand == 1:
                        ret += reference.getSequence(self.chrom, exon.start + 1, variant.pos - 1) + variant.alt + reference.getSequence(self.chrom, variant.pos + len(variant.ref), exon.end)
                    else:
                        temp = Sequence(reference.getSequence(self.chrom, exon.start + 1, variant.pos - 1) + variant.alt + reference.getSequence(self.chrom, variant.pos + len(variant.ref), exon.end))
                        ret += temp.reverseComplement()
                    continue

            if self.strand == 1: ret += reference.getSequence(self.chrom, exon.start + 1, exon.end)
            else: ret += reference.getSequence(self.chrom, exon.start + 1, exon.end).reverseComplement()

        ret = ret[self.codingStart - 1:]
        return ret

    # Getting the translated protein sequence of the transcript
    def getProteinSequence(self, reference, variant):
        codingsequence  = self.getCodingSequence(reference, variant)
        ret = Sequence(codingsequence).translate(1)
        if 'X' in ret: ret = ret[:ret.find('X')]
        return ret

############################################################################################################

# Class representing a single exon
class Exon(object):

    # Constructor
    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end
        self.length = end - start

    # Returns true if exon contains the given position
    def contains(self, pos):
        return (self.start+1 <= pos <= self.end)

############################################################################################################

# Class representing the reference genome
class Reference(object):

    # Constructor
    def __init__(self, filename):
        # Openning tabix file representing the reference genome
        self.fastafile = pysam.Fastafile(filename)

    # Retrieving the sequence of a genomic region
    def getSequence(self, chrom, start, end):

        # Checking if chromosome name exists
        goodchrom = chrom
        if not goodchrom in self.fastafile.references:
            goodchrom = 'chr' + chrom
            if not goodchrom in self.fastafile.references:
                if chrom == 'MT':
                    goodchrom = 'chrM'
                    if not goodchrom in self.fastafile.references: return None
                else:
                    return None

        # Fetching data from reference genome
        if end < start: return Sequence('')
        if start < 1: start = 1

        if pysam.__version__ in ['0.7.7', '0.7.8', '0.8.0']:
            last = self.fastafile.getReferenceLength(goodchrom)
        else:
            last = self.fastafile.get_reference_length(goodchrom)

        if end > last: end = last
        seq = self.fastafile.fetch(goodchrom, start - 1, end)
        return Sequence(seq.upper())

############################################################################################################

# Class representing a sequence
class Sequence(str):

    # Translating to amino acid sequence
    def translate(self, letter):
        if letter == 1:
            gencode = {
                'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'}
        if letter == 3:
            gencode = {
                'ATA': 'Ile', 'ATC': 'Ile', 'ATT': 'Ile', 'ATG': 'Met',
                'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACT': 'Thr',
                'AAC': 'Asn', 'AAT': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
                'AGC': 'Ser', 'AGT': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
                'CTA': 'Leu', 'CTC': 'Leu', 'CTG': 'Leu', 'CTT': 'Leu',
                'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCT': 'Pro',
                'CAC': 'His', 'CAT': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
                'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGT': 'Arg',
                'GTA': 'Val', 'GTC': 'Val', 'GTG': 'Val', 'GTT': 'Val',
                'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCT': 'Ala',
                'GAC': 'Asp', 'GAT': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
                'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGT': 'Gly',
                'TCA': 'Ser', 'TCC': 'Ser', 'TCG': 'Ser', 'TCT': 'Ser',
                'TTC': 'Phe', 'TTT': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
                'TAC': 'Tyr', 'TAT': 'Tyr', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'Cys', 'TGT': 'Cys', 'TGA': 'X', 'TGG': 'Trp'}
        ret = ''
        index = 0
        while index + 3 <= len(self):
            codon = self[index:index + 3].upper()
            if 'N' in codon:
                ret += '?'
                index += 3
                continue
            ret += gencode[codon]
            index += 3
        return ret

    # Getting reverse complement sequence
    def reverseComplement(self):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}
        ret = Sequence()
        for base in self[::-1]: ret += complement[base]
        return ret


############################################################################################################