class ProgressiveAnalysis(object):
    """This class performs Progressive analysis over a full class AlignmentInformation instance. It discards alleles
    in the Reference instance that do not mach the information in the AlignmentInformation instance.
    Output is a list of result
    alleles per KIR gene
    ...

    Attributes
    ----------
    result_alleles_dictionary : dict
        Dictionary of lists, one per detected KIR gene, storing names of result alleles from Progressive Analysis
    exon_identical_alleles : dict
        Dictionary of lists per KIR gene, storing names of exon identical alleles to those in result_alleles_dictionary

    Methods
    -------
    process_result_alleles
        Reads list of result alleles from progressive analysis, separate them into their specific KIR genes
    exclude_identical_alleles
        Looks for exon identical alleles in the result_alleles_dictionary and exclude them from it

    """

    def __init__(self, alignment_information_instance, proportion_threshold=5):
        """Iterates over the class AlignmentInformation proportion_dictionary instance positions.
        Determines present nucleotides per position based on proportion_threshold parameter, add them to dictionary.
        Per position, iterates over class AlignmentInformation reference_file,
        discards alleles that do not contain any of the present nucleotides (excluding gaps) in the same index position.
        Discarded alleles are excluded from the reference_alleles list and not considered in further positions

        Parameters
        ----------
        alignment_information_instance : class AlignmentInformation instance
            Instance inherited from class AlignmentInformation
        proportion_threshold : int
            Percentage threshold above which nucleotides are defined as present per position, given their proportions
        """

        self.result_alleles_dictionary = {}
        self.exon_identical_alleles = {}
        for position in alignment_information_instance.proportion_dictionary:
            present_nucleotides = []
            for nucleotide in alignment_information_instance.proportion_dictionary[position]:
                if alignment_information_instance.proportion_dictionary[position][nucleotide] > proportion_threshold:
                    present_nucleotides.append(nucleotide)
            alignment_information_instance.proportion_dictionary[position]["Present Nucleotides"] = present_nucleotides
            primary_result_alleles = alignment_information_instance.reference_alleles
            with open(alignment_information_instance.reference.file_name) as alignment_reference:
                for allele in alignment_reference.readlines()[12:-1]:  # First sequence in reference is in line 12
                    allele = allele.replace(" ", "")
                    allele_id = allele.split("\t")[0]
                    if allele_id in primary_result_alleles:
                        allele_sequence = allele.split("\t")[:-4]  # Last four items in the list are \t characters
                        allele_sequence = allele_sequence[1].replace("|", "")
                        allele_sequence = list(allele_sequence)
                        if allele_sequence[position] not in present_nucleotides:
                            primary_result_alleles.remove(allele_id)
            alignment_reference.close()
        self.process_result_alleles(primary_result_alleles)

    def process_result_alleles(self, primary_result_alleles):
        """Reads list of result alleles from progressive analysis, separate them into their specific KIR genes.
        Calls methods exclude_identical_alleles and progressive results

        Parameters
        ----------
        primary_result_alleles : list
            List of all alleles that weren't discarded from reference_alleles list, after progressive analysis
        """

        for allele in primary_result_alleles:
            kir_gene = allele.split("*")[0]
            if kir_gene == "KIR2DL5A" or kir_gene == "KIR2DL5B":
                kir_gene = "KIR2DL5"  # KIR genes 2DL5A and 2DL5B are analyzed as a single gene locus
            if kir_gene not in self.result_alleles_dictionary:
                self.result_alleles_dictionary[kir_gene] = []
            self.result_alleles_dictionary[kir_gene].append(allele)
        self.exclude_exon_identical_alleles()
        #self.progressive_results()

    def exclude_exon_identical_alleles(self):
        """Looks for exon identical alleles in the self.result_alleles_dictionary and exclude them from it. Returns
        a replicate dictionary with only excluded exon identical alleles per detected KIR gene

        """

        for kir_gene in self.result_alleles_dictionary:
            self.exon_identical_alleles[kir_gene] = []
            exon_mismatches_codes = []
            for allele in self.result_alleles_dictionary[kir_gene]:
                allele_code = allele.split("*")[1]
                allele_code = list(allele_code)
                if len(allele_code) > 4:  # This allele code is longer than 4 digits, exon identical alleles exist
                    exon_mismatch_code = allele_code[:5]  # Last two digits indicate exon mismatches
                    if exon_mismatch_code in exon_mismatches_codes:
                        self.exon_identical_alleles[kir_gene].append(allele)
                    else:
                        exon_mismatches_codes.append(exon_mismatch_code)
            self.result_alleles_dictionary[kir_gene] = [allele for allele in self.result_alleles_dictionary[kir_gene] if
                                                     allele not in self.exon_identical_alleles[kir_gene]]
            self.result_alleles_dictionary[kir_gene] = dict.fromkeys(self.result_alleles_dictionary[kir_gene], [])
