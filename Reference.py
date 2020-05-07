class Reference(object):
    """
    This class takes an .ipd format file (compatible with NGSengine) as argument
    and creates a Reference instance with lists containing starting and ending indexes of every genomic region
    In addition, print_sequence method allows to print a framed sequence from a selected KIR gene in the reference
    ...

    Attributes
    ----------
    file_name : str
        Formatted string with the path of the reference file given as argument
    first_sequence : list
        Contains genomic sequences of every region (starts and finishes with introns while alternating between both types)
        corresponding to the first allelic sequence found in the reference
    regions_index_list : list
        Contains starting/ending nucleotide positions of alternate genomic regions (intron-exon-intern)
    introns_index_list : list
        Contains tuples corresponding to the range of positions within every intron region is comprised
    exons_index_list : list
        Contains tuples corresponding to the range of positions within every intron region is comprised

    Methods
    -------
    get_regions_index
        Extracts regions_index_list from first_sequence
        Creates exons_index_list and introns_index_list out of regions_index_list
    print_sequence
        Prints sequence between two given nucleotide positions, of the first allele in the selected KIR gene
        from the reference
    """

    def __init__(self, reference_file):
        """
        Parameters
        ----------
        reference_file : str
            The path of the .ipd format file that was selected as a reference
        """

        self.file_name = reference_file
        with open(self.file_name) as alignment_reference:
            first_sequence = alignment_reference.readlines()[12]  # 12th line is the first allelic sequence in file
            first_sequence = first_sequence.replace(" ", "")
            first_sequence = first_sequence.split("\t")[:-4]  # Last four items in the list are \t characters
            self.first_sequence = first_sequence[1].split("|")
        alignment_reference.close()
        self.regions_index_list = [0]  # First region starts in position 0
        self.introns_index_list = []
        self.exons_index_list = []

    def get_regions_index(self):
        """Extracts regions_index_list from first_sequence
        Creates exons_index_list and introns_index_list out of regions_index_list
        """

        for region_index in range(len(self.first_sequence)):
            region_pos = len(self.first_sequence[region_index]) + self.regions_index_list[region_index]
            self.regions_index_list.append(region_pos)

        for i in range(len(self.regions_index_list) - 1):
            if i % 2 == 0:
                intron_index = (self.regions_index_list[i], self.regions_index_list[i + 1])
                self.introns_index_list.append(intron_index)
            elif i % 2 == 1:  # should be else
                exon_index = (self.regions_index_list[i], self.regions_index_list[i + 1])
                self.exons_index_list.append(exon_index)

    def print_sequence(self, leftmost_position, rightmost_position, kir_gene):
        """Prints sequence between two given nucleotide positions of the first allele in the selected KIR gene of the
        reference

        Parameters
        ----------
        leftmost_position : int
            Leftmost nucleotide position of the segment of the reference sequence that is to be printed
        rightmost_position : int
            Rightmost nucleotide position of the segment of the reference sequence that is to be printed
        kir_gene : str
            KIR gene in the reference alignment from which its first allelic sequence is to be printed

        Raises
        ------
        If selected KIR gene is not in reference, a warning message is shown
        """

        with open(self.file_name) as alignment_reference:
            for allele_sequence in alignment_reference.readlines()[12:]:
                allele_sequence = allele_sequence.replace(" ", "")
                reference_gene = allele_sequence.split("*")[0]
                if reference_gene == kir_gene:
                    allele_sequence = allele_sequence.split("\t")[:-4]
                    allele_sequence = allele_sequence[1].replace("|", "")
                    allele_sequence = list(allele_sequence)
                    print(allele_sequence[leftmost_position:rightmost_position])
                    break
            print("WARNING: Selected KIR genes was not found in reference")
        alignment_reference.close()
