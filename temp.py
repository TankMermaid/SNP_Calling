import os
import sys
import subprocess
import time
from Bio import SeqIO
from contextlib import contextmanager





class Fastq(object):
    def __init__(self, dir_location):
        self.dir_location = dir_location

     
    def find_fastq_files(self):
        #  search for fastqfiles. returns a list of all fastq files.
        cmd = [ 'find', '.' , '-maxdepth','1', '-type', 'f', '-name', '*.fastq' ]
        output = subprocess.Popen( cmd, stdout=subprocess.PIPE ).communicate()[0]
        list_of_fastq_files = [fastq_file.strip('\n') for fastq_file in output.split('./')]
        return list_of_fastq_files[1:]


    def find_fastq_names(self):
        # returns list of fastq names
        os.chdir(self.dir_location)
        fastq_list = []
        list_of_fastq_files = self.find_fastq_files()
        for i in range(len(list_of_fastq_files)):
            fastq_list.append(list_of_fastq_files[i].split('/')[-1])
        return fastq_list    
                     
    def concatenate_fastq(self):
        # recieves list of fastq files and concatenate them to one fastq file in another dir(otherwise - loop occured).
        os.chdir(self.dir_location)
        fastq_name_list_in_str = ' '.join(self.find_fastq_names())
        os.mkdir('concatenate_fastq')
        os.system("cat " + fastq_name_list_in_str + " > " + "concatenate_fastq/concatenate_file.fastq")
        
class Fasta(object):
    def __init__(self, dir_location, haplotype= None):
        self.dir_location = dir_location
        self.haplotype = haplotype
        
    def find_fasta_file(self):
        #find fasta file in the current dir
        os.chdir(self.dir_location)
        cmd = [ 'find', '.' , '-maxdepth','1', '-type', 'f', '-name', '*.fasta' ]
        output = subprocess.Popen( cmd, stdout=subprocess.PIPE ).communicate()[0]
        list_of_fasta_files = [fasta_file.strip('\n') for fasta_file in output.split('./')]
        if len(list_of_fasta_files) is 2:                       
            return list_of_fasta_files[1]
        else:
            return None
      
    def filter_fasta_haplotype(self):
       # filter fasta accoring to given haplotype, default is None
       if self.haplotype in ('A', 'a', 'B', 'b'):
           os.chdir(self.dir_location)
           os.mkdir('FASTA_files')
           with open('FASTA_files/' + self.find_fasta_file().strip('.fasta') + '_haplotype_' + self.haplotype.upper() + '.fasta', "w") as output_handle:
               for record in SeqIO.parse(self.find_fasta_file(), "fasta"):
                   if self.haplotype.upper() in record.id:
                       output_handle.write('>' + record.id)
                       output_handle.write('\n')
                       output_handle.write(str(record.seq))
                       output_handle.write('\n')
       else:
           return None

        



def is_module_found(module, project_dir):
    #recieves a module name and returns true if found via 'which'.
    #if not found by 'which', find via 'find' command in the specific project folder - picard.jar and gatk.jar has to be within this folder
    os.chdir('/Users/bermanlab')
    cmd = [ 'which', module ]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE ).communicate()[0]
    if module in output:
        return True
    else:
        os.chdir(project_dir)    
        cmd = [ 'find', '.', '-maxdepth', '1' , '-type', 'f', '-name', module ]
        output = subprocess.Popen( cmd, stdout=subprocess.PIPE ).communicate()[0]
        if module in output:
            return True
        else:
            return False
        
def modules_not_found(list_of_modules, project_dir):
    #recieves a list of modules and the project dir.
    #using is_module_found function, returns a list of modules that are missing.
    missing_module_list = []
    for module in list_of_modules:
        if is_module_found(module, project_dir) == False:
            missing_module_list.append(module)
    return missing_module_list

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout    
            
#sys.argv[1] = dir of project
#sys.argv[2] = paired end or not by '1' or '2' 
#sys.argv[3] = haplotype 'A', 'a', 'B', 'b'           
              
project_dir = sys.argv[1]

if sys.argv[1] in ('-h', '--h' ,'-help' , '--help') or sys.argv[2] in ('-h', '--h' ,'-help' , '--help'): #help message    
    sys.exit('\n\n\n--------------SNP_Calling HELP------------- \n please use the command as followed - order is important! \n\n python snp_calling.py [project directory] ["1" - not paired-end | "2" - paired end] [haplotype - "A" or "B", default = None] \n\n PLEASE make sure you have picard.jar and GenomeAnalysisTK.jar in the project dir \n\n FASTQ and fasta files will be auto-found - PLEASE make sure to put them in the project dir ass well! \n\n THANK YOU!!! \n Noam Shahar\n\n\n')             
 
      
else:
    
    paired_end = sys.argv[2]
    haplotype = sys.argv[3]     
# ********************* Modules checking *************************      
print "Checking existence of essential modules..."
modules_list = ['bowtie2', 'GenomeAnalysisTK.jar', 'picard.jar', 'samtools', 'R', 'snpeff']

try: 
    if modules_not_found(modules_list, project_dir):   
        raise AssertionError("Please install: " + ', '.join(modules_not_found(modules_list, project_dir)))


except AssertionError("Please install: " + ', '.join(modules_not_found(modules_list, project_dir))):
    sys.exit()
     
       
print 'All modules found!'

# ********************* Modules checking *************************  


# ******************** FASTA and FASTQ files changes *****************
os.chdir(project_dir)
fastq_files = Fastq(project_dir)
fasta_file = Fasta(project_dir, haplotype)
fasta_file.filter_fasta_haplotype()
fasta = 'FASTA_files/' + fasta_file.find_fasta_file().strip('.fasta') + '_haplotype_' + haplotype.upper() + '.fasta'

print 'Creating FASTA dict and index files...'
os.system("java -jar picard.jar CreateSequenceDictionary R= " + fasta + " O= " + fasta.strip('.fasta') + '.dict >/dev/null')
os.system("samtools faidx " + fasta + " >/dev/null")
os.system("bowtie2-build " + fasta + " FASTA_files/index_file >/dev/null")



print 'FASTA file found: ' + fasta
time.sleep(2)
fastq = fastq_files.find_fastq_names()
print 'FastQ files found: ' + ', '.join(fastq)
time.sleep(2)
try:

    if len(fastq) > 1 and paired_end == '1':
    #not a paired end, and more than 1 fastq file - does concatenate! 
        print 'Concatenating FASTQ...'
        fastq_files.concatenate_fastq()

    elif len(fastq) == 1 and paired_end == '2':
        raise AssertionError('paired')
     
except AssertionError('paired'):
    sys.exit('Error: problem with FASTQ File')
# ******************** FASTA and FASTQ files changes *****************

os.chdir(project_dir)
os.mkdir('bam_files')
os.mkdir('metrices') 
os.mkdir('txt_files') 
os.mkdir('vcf_files') 
os.mkdir('final_output')
print "bowtie2 haplotype " + haplotype + " in progress... may take a while" 
with suppress_stdout():
    if paired_end == '1':
        os.system("bowtie2 --rg-id test --rg SM:test --rg LB:GRC --rg PL:ILLUMINA --rg DS:HiSeq2000 -x FASTA_files/index_file -U " +  "concatenate_fastq/concatenate_file.fastq" + " -S bam_files/aligned_reads.sam ")
    elif paired_end == '2':
        os.system("bowtie2 --rg-id test --rg SM:test --rg LB:GRC --rg PL:ILLUMINA --rg DS:HiSeq2000 -x FASTA_files/index_file -1 " +  ' '.join(fastq) + "-S bam_files/aligned_reads.sam")

#Step2 - Sort SAM file by coordinate, convert to BAM
os.system("java -jar picard.jar SortSam INPUT=bam_files/aligned_reads.sam OUTPUT=bam_files/sorted_reads.bam SORT_ORDER=coordinate")

#Step3 - Collect Alignment & Insert Size Metrics
os.system("java -jar picard.jar CollectAlignmentSummaryMetrics R= " + fasta + " I=bam_files/sorted_reads.bam O=metrices/alignment_metrics.txt")
os.system("java -jar picard.jar CollectInsertSizeMetrics INPUT=bam_files/sorted_reads.bam OUTPUT=metrices/insert_metrics.txt HISTOGRAM_FILE=metrices/insert_size_histogram.pdf")
os.system("samtools depth -a bam_files/sorted_reads.bam > txt_files/depth_out.txt")

#Step4 - Mark Duplicates
os.system("java -jar picard.jar MarkDuplicates INPUT=bam_files/sorted_reads.bam OUTPUT=bam_files/dedup_reads.bam METRICS_FILE=metrices/metrics.txt")

#Step5 - Build BAM Index
os.system("java -jar picard.jar BuildBamIndex INPUT=bam_files/dedup_reads.bam")

#Step6 - Create Realignment Targets
os.system("java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R " + fasta + " -I bam_files/dedup_reads.bam -o txt_files/realignment_targets.list")

#Step7 - Realign Indels
os.system("java -jar GenomeAnalysisTK.jar -T IndelRealigner -R " + fasta + " -I bam_files/dedup_reads.bam -targetIntervals txt_files/realignment_targets.list -o bam_files/realigned_reads.bam")

#Step8 - Call Variants
os.system("java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R " + fasta + " -I bam_files/realigned_reads.bam -o vcf_files/raw_variants.vcf")

#Step9 - Extract SNPs & Indels
os.system("java -jar GenomeAnalysisTK.jar -T SelectVariants -R " + fasta + " -V vcf_files/raw_variants.vcf -selectType SNP -o vcf_files/raw_snps.vcf")
os.system("java -jar GenomeAnalysisTK.jar -T SelectVariants -R " + fasta + " -V vcf_files/raw_variants.vcf -selectType INDEL -o vcf_files/raw_indels.vcf")

#Step10 - Filter SNPs
os.system("java -jar GenomeAnalysisTK.jar -T VariantFiltration -R " + fasta + " -V vcf_files/raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' -o vcf_files/filtered_snps.vcf")

#Step11 - Filter Indels
os.system("java -jar GenomeAnalysisTK.jar -T VariantFiltration -R " + fasta + " -V vcf_files/raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o vcf_files/filtered_indels.vcf")

#Step12 - Base Quality Score Recalibration (BQSR) #1
os.system("java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R " + fasta + " -I bam_files/realigned_reads.bam -knownSites vcf_files/filtered_snps.vcf -knownSites vcf_files/filtered_indels.vcf -o txt_files/recal_data.table")

#Step13 - Base Quality Score Recalibration (BQSR) #2
os.system("java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R " + fasta + " -I bam_files/realigned_reads.bam -knownSites vcf_files/filtered_snps.vcf -knownSites vcf_files/filtered_indels.vcf -BQSR txt_files/recal_data.table -o txt_files/post_recal_data.table")

#Step14 - Analyze Covariates
os.system("java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R " + fasta + " -before txt_files/recal_data.table -after txt_files/post_recal_data.table -plots txt_files/recalibration_plots.pdf")

#Step15 - Apply BQSR
os.system("java -jar GenomeAnalysisTK.jar -T PrintReads -R " + fasta + " -I bam_files/realigned_reads.bam -BQSR txt_files/recal_data.table -o bam_files/recal_reads.bam")

#Step16 - Call Variants
os.system("java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R " + fasta + " -I bam_files/recal_reads.bam -o vcf_files/raw_variants_recal.vcf")

#Step17 - Extract SNPs & Indels
os.system("java -jar GenomeAnalysisTK.jar -T SelectVariants -R " + fasta + " -V vcf_files/raw_variants_recal.vcf -selectType SNP -o vcf_files/raw_snps_recal.vcf")
os.system("java -jar GenomeAnalysisTK.jar -T SelectVariants -R " + fasta + " -V vcf_files/raw_variants_recal.vcf -selectType INDEL -o vcf_files/raw_indels_recal.vcf")

#Step18 - Filter SNPs
os.system("java -jar GenomeAnalysisTK.jar -T VariantFiltration -R " + fasta + " -V vcf_files/raw_snps_recal.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' -o final_output/filtered_snps_final.vcf")

#Step19 - Filter Indels
os.system("java -jar GenomeAnalysisTK.jar -T VariantFiltration -R " + fasta + " -V vcf_files/raw_indels_recal.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o final_output/filtered_indels_final.vcf")

#Step20 - Make annotations
os.system("snpeff albicans22 final_output/filtered_snps_final.vcf > final_output/filtered_snps_final.ann.vcf")
os.system("snpeff albicans22 final_output/filtered_indels_final.vcf > final_output/filtered_indels_final.ann.vcf") 
    
print "************************ D O N E *************************"



    
    



      
        
    
        