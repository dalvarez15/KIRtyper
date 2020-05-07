class Read(object):
    """This class takes a single line of a .sam format file and parses it. Every line represents a DNA read aligned
    to a reference. A read instance stores all information about how the read is aligned to the reference,
    for analysis purposes.
    ...

    Attributes
    ----------
    id : str
        String with read identification information
    flag : int
        0, 4 or 16, flag information of the read indicating forward, reverse or non-aligned read
    leftmost_position : int
        Position of the reference alignment to which the read's first nucleotide is aligned
    quality : int
        Integer that reflects quality of the read
    cigar_lst : list
        Every character and digit in the cigar string separated as items in a list
    sequence : str
        Genomic sequence of the read
    parsed_cigar : lst
        List of cigar operators as list separating cigar alphabetic signs and consumed positions per operator
    aligned_sequence : lst
        Genomic sequence of the read including gaps as stated by cigar string,
        so it is aligned to the reference alignment
    rightmost_position : int
        Position of the reference alignment to which the read's last nucleotide is aligned

    Methods
    -------
    parse_cigar
        Reads cigar_lst and returns a parsed_cigar list. It separates every cigar operator
        with its correspondent consumed positions, in different items
    get_aligned_sequence
        Combines the parsed_cigar list and read sequence to return the aligned (to reference) sequence with gaps
        and deletions/insertions
        Returns rightmost position based on alignment
    is_aligned_to
        Returns true or false if read is aligned to a given range of positions corresponding to the reference alignment
    """

    def __init__(self, sam_line):
        """Gets the relevant SAM characters in the line by splitting them by tabs and selecting their indexes
        Parameters
        ----------
        sam_line : str
            Single line of a SAM file, can be obtained by iterating through the file using self.readlines()
        """

        sam_line = sam_line.split("\t")
        self.id = sam_line[0]
        self.flag = int(sam_line[1])
        self.leftmost_position = int(sam_line[3]) - 1  # -1 fixes index from NGSengine, which starts in 1 instead of 0
        self.quality = int(sam_line[4])
        self.cigar_lst = list(sam_line[5])
        self.sequence = sam_line[9]
        self.parsed_cigar = []
        self.aligned_sequence = []
        self.rightmost_position = None

    def parse_cigar(self):
        """Parses cigar list with separate digits/operators as items
        and returns a parsed cigar list of operator combined with consumed positions as sublists

        Items in the cigar list: M is a match, D a deletion and I is an insertion
        """

        cigar_sublist = []
        for i in range(len(self.cigar_lst)):
            if self.cigar_lst[i] == "M" or self.cigar_lst[i] == "D" or self.cigar_lst[i] == "I":
                new_str = ""
                for digit in cigar_sublist:
                    new_str += digit
                cigar_sublist = [new_str]
                cigar_sublist.append(self.cigar_lst[i])
                self.parsed_cigar.append(cigar_sublist)
                cigar_sublist = []
            else:
                cigar_sublist.append(self.cigar_lst[i])
        return self.parsed_cigar

    def get_aligned_sequence(self):
        """Iterates over parsed cigar list and read sequence, returns sequence aligned to reference with indicated gaps,
        insertions and deletions
        Returns rightmost position to which the read is aligned based on length of the aligned read sequence

        Raises
        ------
        If no cigar string is present, warning is printed and no aligned sequence or rightmost position are given
        """

        sequence_lst = list(self.sequence)
        if self.parsed_cigar is []:
            return "WARNING: no cigar string found for this read"
        else:
            last_position = 0
            for sub_cigar in self.parsed_cigar:
                consumed_positions = int(sub_cigar[0])
                operator = sub_cigar[-1]
                if operator == "M":
                    for i in range(last_position, consumed_positions + last_position):
                        self.aligned_sequence += sequence_lst[i]
                    last_position += consumed_positions
                elif operator == "D":
                    for i in range(consumed_positions):
                        self.aligned_sequence += "."
                elif operator == "I":
                    last_position += consumed_positions
            self.rightmost_position = self.leftmost_position + len(self.aligned_sequence)
            return self.aligned_sequence, self.rightmost_position

    def is_aligned_to(self, region_range):
        """Returns true or false depending on whether the read is aligned/overlaps within a given range of positions

        Parameters
        ----------
        region_range : tuple/list
            Range (first two items are taken) of nucleotide positions in the reference alignment, to which the read
            is going to be checked whether aligns or not
            Can be a tuple in reference.exons_index_list or reference.introns_index_list so reads can be checked to
            align to a specific region in the reference
        """

        for i in range(self.leftmost_position, self.rightmost_position):
            if i in range(region_range[0], region_range[1]):
                return True
        return False
