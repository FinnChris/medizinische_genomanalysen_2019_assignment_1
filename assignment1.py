import mysql.connector
import pysam
import pybedtools

__author__ = 'Christian Jansen'

##
## Concept:
## TODO
##


class Assignment1:
    
    def __init__(self, gene = "S100B", bamfile = "chr21.bam"):
        ## Your gene of interest
        self.gene = gene
        self.name2 = ""
        self.location = ""
        self.gstart = 0
        self.gstop = 0
        self.direction = ""
        self.n_exon = 0
        self.ex_start = []
        self.ex_stop = []
        self.bamfile = bamfile
        self.sam = pysam.AlignmentFile(self.bamfile,"rb")

    
    def download_gene_coordinates(self, genome_reference, file_name):
        ## TODO concept
        
        print("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        
        ## Execute query
        cursor.execute(query)
        
        ## Write to file
        ## TODO this may need some work 
        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")
    
            
        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print("Done fetching data")
        
    def get_coordinates_of_gene(self, input_file = "coordinates_S100B.txt"):
        ## Use UCSC file
        coordinates = open(input_file,"r")

        for record in coordinates.readlines():
            list = record.replace("('","").replace("',","").replace(" '"," ").replace("')","").replace(", "," ").split(sep=" ")

            if list[0] == self.gene:
                self.name2 = list[1]
                self.location = list[2]
                self.gstart = int(list[3])
                self.gstop = int(list[4])
                self.direction = list[5]
                self.n_exon = int(list[6])
                self.ex_start = list[7].strip("b'").split(",")
                self.ex_stop = list[8].strip("b'").strip(",\n").split(",")

                break

        
    def get_gene_symbol(self):
        return self.name2
                        
    def get_sam_header(self):
        return self.sam.header['HD']
        
    def get_properly_paired_reads_of_gene(self):
        i=0
        for read in self.sam.fetch(self.location,self.gstart, self.gstop):
            if read.is_proper_pair:
                i+=1
        return i
        
    def get_gene_reads_with_indels(self):
        i = 0
        for read in self.sam.fetch(self.location,self.gstart, self.gstop):
            if not read.is_unmapped:
                indel = False
                for align in read.cigar:
                    if align[0] == 1 or align[0] == 2:
                        indel = True
                if indel:
                    i+=1
        return i

    def calculate_total_average_coverage(self):
        pybed_sam = pybedtools.BedTool(self.bamfile)
        coverage = pybed_sam.genome_coverage(bg=True)
        sum_cov = 0.0
        count = 0
        for line in coverage:
            sum_cov += float(line[3])
            count += 1
        total_average_cov = sum_cov/count
        return total_average_cov

    def calculate_gene_average_coverage(self):
        print("todo")
        
    def get_number_mapped_reads(self):
        print("todo")

    def get_region_of_gene(self):
        print("todo")
        
    def get_number_of_exons(self):
        return self.n_exon
    
    
    def print_summary(self):
        print("Print all results here")
        print("Gene name:", self.gene)
        print("Gene symbol:", self.get_gene_symbol())
        print("Sam Header:", self.get_sam_header())
        print("Number of properly paired reads:", self.get_properly_paired_reads_of_gene())
        print("Gene reads containing indels:", self.get_gene_reads_with_indels())
        print("Total average coverage:", self.calculate_total_average_coverage())
        print("Gene average coverage:", self.calculate_gene_average_coverage())
        print("Number of mapped reads:", self.get_number_mapped_reads())
        print("Region of gene:", self.get_region_of_gene())
        print("Number of exons:", self.n_exon)

    
    
def main():
    print("Assignment 1")
    assignment1 = Assignment1()
    #assignment1.print_summary()
    print(assignment1.calculate_total_average_coverage())
    #assignment1.get_coordinates_of_gene()
    #print(assignment1.get_sam_header())
    #print(vars(assignment1))
    #assignment1.get_properly_paired_reads_of_gene()
    #assignment1.get_gene_reads_with_indels()
    #assignment1.download_gene_coordinates("hg38","coordinates_" + assignment1.gene + ".txt")
    
    
    print("Done with assignment 1")
    
        
if __name__ == '__main__':
    main()



    
