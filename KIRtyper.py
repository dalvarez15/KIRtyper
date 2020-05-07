from sys import argv
import gzip
from Reference import *
from Read import *
from AlignmentInformation import *
from ProgressiveAnalysis import *
from CombinedAnalysis import *

"""
Main module of KIR Typer bioinformatics pipeline. It performs typing of human KIR genes. Input is an alignment in
SAM format, and a reference file. This reference file (KIR_full.ipd) was previously used to align DNA reads from NGS
with NGSengine. Typing results are written into indicated output text file. 

Required libraries and packages: gzip, pandas, itertools, numpy

Use: 
    Command line: python3 KIRtyper.py [reference.ipd] [samfile.sam] [output.txt]
"""


def main():
    reference = Reference(argv[1])
    reference.get_regions_index()
    print("Reference %s processed" % argv[1])
    with gzip.open(argv[2], 'rt') as sam_file:
        alignment = AlignmentInformation(reference)
        print("Processing alignment information in SAM file...")
        for line in sam_file.readlines()[2:]:  # Information of reads alignment starts in 2nd line
            read = Read(line)
            if read.flag != 4 and read.quality == 255:
                '''
                Reads with flag value 4 (i.e. not aligned) and quality value different than 255 (optimal) 
                are not processed
                '''
                read.parse_cigar()
                read.get_aligned_sequence()
                for exon_range in reference.exons_index_list:
                    in_exons = read.is_aligned_to(exon_range)
                    if in_exons is True:  # Only exon information is accounted
                        alignment.update_count_dictionary(read)
                        break
        alignment.create_proportion_dictionary()
        print("Progressive analysis in progress...")
        progressive_analysis = ProgressiveAnalysis(alignment)  # Proportion_threshold value customizable, default=5
        print("%i genes/s detected" % len(list(progressive_analysis.result_alleles_dictionary.keys())))
        print("Combined analysis in progress...")
        combined_analysis = CombinedAnalysis(progressive_analysis, alignment)
        write_output_file(progressive_analysis, combined_analysis, argv[2], argv[3])
        print("Analysis done, results written into output file: %s" % argv[3])
        sam_file.close()
    sam_file.close()


def write_output_file(progressive_analysis_instance, combined_analysis_instance, sam_file, output_file_name):
    with open(output_file_name, 'a') as output_file:
        output_file.write("%s \n" % sam_file)
        output_file.write("Progressive Analysis Results:\n")
        output_file.write("%i genes/s detected: " % len(list(progressive_analysis_instance.result_alleles_dictionary.keys())))
        output_file.write(", ".join(list(progressive_analysis_instance.result_alleles_dictionary.keys())))
        for locus in progressive_analysis_instance.result_alleles_dictionary:
            output_file.write("\nResult %i allele/s from %s: " % (len(list(progressive_analysis_instance.result_alleles_dictionary[locus].keys())), locus))
            output_file.write((", ".join(list(progressive_analysis_instance.result_alleles_dictionary[locus].keys()))))
        if combined_analysis_instance.typing_result is None:
            output_file.write("\nCombined analysis not applicable\n")
        else:
            output_file.write("\nCombined analysis results:\n")
            if len(list(combined_analysis_instance.typing_result)) > 0:
                output_file.write(
                    "%i Genotype combinations matching alignment data, out of %i primary combinations\n" % (
                    len(list(combined_analysis_instance.typing_result)),
                    len(combined_analysis_instance.all_genotype_combinations)))
                for item in list(combined_analysis_instance.typing_result):
                    output_file.write('+'.join(str(genotype) for genotype in item))
                    output_file.write("-")
                output_file.write("\nUpdated list of result alleles per detected gene:\n")
                for locus in combined_analysis_instance.result_alleles_dictionary:
                    output_file.write("Result %i allele/s from %s: " % (
                        len(combined_analysis_instance.result_alleles_dictionary[locus]), locus))
                    output_file.write(", ".join(combined_analysis_instance.result_alleles_dictionary[locus]))
                    output_file.write("\n")
            else:
                output_file.write("No genotype combination matches alignment information, "
                                  "out of %i primary combinations\n"
                                  % len(combined_analysis_instance.all_genotype_combinations))
    output_file.close()
    return "Output written"


if __name__ == "__main__":
    main()

