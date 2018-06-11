#!/usr/bin/env python
import os
import sys
import re
import json
import optparse
import shutil
import pandas
import datetime
import time

class FusionAMOI():


    def __init__(self, config_file, out_path, annotation):
        self.config_data = json.load(open(config_file))
        self.out_path = out_path
        self.out_log = out_path + "/output.log"
        self.err_log = out_path + "/error.log"
        self.annotation = annotation
        
    def run_pegasus(self, in_file, forced=True):
        out_path = self.out_path
        annotation = self.annotation
        pegasus_home = self.config_data["pegasus_home"]
        pegasus_config_fn = pegasus_home + '/config_' + annotation + '.txt'
        protein_file = pegasus_home + '/ref/' + annotation + '.pep'
        with open(pegasus_config_fn) as f:
            for line in f:
                line = line.rstrip('\n')
                fields = line.split("\t")
                if fields[0] == "gtf_file":
                    gtf_file = fields[1]

        
        data_spec_fn = out_path + "/data_spec.txt"
        self.out_log = out_path + "/output.log"
        log_dir = out_path + "/pegasus_log"
        out_dir = out_path + "/pegasus_out"
        self.out_dir = out_dir
        pegasus_in_file = out_path + "/input.tsv"
        onco_gene_file = out_path + "/oncogene_fusion.tsv"
        not_onco_gene_file = out_path + "/not_oncogene_fusion.tsv"
        pegasus_out_file = out_dir + "/pegasus.output.txt"
        out_file = out_path + "/output.json"
        self.filter(in_file, onco_gene_file, not_onco_gene_file)
        if not os.path.exists(pegasus_out_file) or forced:
            print("Preparing Pegasus files...")                    
            self.to_pegasus_file(os.path.abspath(onco_gene_file), pegasus_in_file, gtf_file)
            data_spec = open(data_spec_fn, "w")
            data_spec.write("# sample_path\tsample_name\tsample_type\tfusion_program\n")
            data_spec.write(pegasus_in_file + "\t" + os.path.basename(onco_gene_file) + "\ttumor_tissue\tgeneral\n")
            data_spec.close()
            os.mkdir(log_dir)
            os.mkdir(out_dir)
            print("Done.\n")
            cmd = pegasus_home + "/pegasus.pl -c " + pegasus_config_fn + " -d " + data_spec_fn + " -l " + log_dir + " -o " + out_dir + ' >> ' + self.out_log + ' 2>&1'
            self.append_log(cmd)
            print("Running Pegasus...")
            os.system(cmd)
            print("Done.\n")
        #if not os.path.exists(out_file) or forced:
        self.process_output(pegasus_out_file, out_file, protein_file)

    def read_protein_fasta(self, protein_file):
        sequences = {}
        with open(protein_file) as f:
            trans_id = ''
            seq = ''
            for line in f:
                if line[0] == '>':
                    line = line.rstrip('\n')
                    exp = "transcript:(.*)\.+" if self.annotation == "refseq" else "transcript:(.*)"
                    m = re.search(exp, line)
                    
                    if trans_id != '' and seq != '':
                        sequences[trans_id] = seq
                    if m != None:
                        trans_id = m.group(1)
                        
                    else:
                        trans_id = "NA"
                    seq = ''
                    
                else:
                    seq += line
        if trans_id != '' and seq != '':
            sequences[trans_id] = seq
        return sequences


    def process_output(self, pegasus_out_file, out_json_file, protein_file):
        sequences = self.read_protein_fasta(protein_file)
        df = pandas.read_table(pegasus_out_file, sep='\t', index_col=False)
        interproscan_home = self.config_data["interproscan_home"]
        fasta_file = self.out_dir + "/fusion.fasta"
        domain_file = self.out_dir + "/domain.json"        
        fastaf = open(fasta_file, "w")
        id_list = {}
        trans = {}
        for index, row in df.iterrows():
            if "DriverScore" not in row:
                continue
            break1 = row["Gene_Breakpoint1"]
            break2 = row["Gene_Breakpoint2"]
            if row["Strand1"] == "-":
                break1 = break1 - 1
            if row["Strand2"] == "-":
                break2 = break2 - 1 
            seqs = row["Protein_Sequence"].split(',')
            gene1 = row["Gene_Name1"]
            gene2 = row["Gene_Name2"]
            longest_seq = "";
            for seq in seqs:
                if len(seq) > len(longest_seq):
                    longest_seq = seq
            m = re.search('(.*):(.*)-(.*)', longest_seq)
            aa_seq = m.group(1)
            trans1 = m.group(2)
            trans2 = m.group(3)
            if trans1 in sequences and gene1 not in id_list:
                seq1 = sequences[trans1]
                fastaf.write(">" + gene1 + "\n" + seq1)
                id_list[gene1] = ''
            if trans2 in sequences and gene2 not in id_list:
                seq2 = sequences[trans2]
                fastaf.write(">" + gene2 + "\n" + seq2)
                id_list[gene2] = ''    
            seq2 = sequences[trans2]
            fusion_name = gene1 + "-" + gene2 + ":" + str(break1) + "-" + str(break2)
            trans[fusion_name] = {"trans1" : trans1, "trans2" : trans2}
            if fusion_name not in id_list:
                fastaf.write(">" + fusion_name + "\n" + aa_seq + "\n")
                id_list[fusion_name] = ''
        cmd = interproscan_home + "/interproscan.sh -appl Pfam -i " + fasta_file + " -f JSON -o " + domain_file + ' >> ' + self.out_log + ' 2>&1'
        print("Running Interproscan...")
        self.append_log(cmd)
        os.system(cmd)
        print("Done.\n")
        print("Post processing output...")
        domain_data = {}
        with open(domain_file) as domainf:
            data = json.load(domainf)
        for result in data["results"]:
            fusion_name = result["xref"][0]["name"]
            domains = []
            for domain in result["matches"]:
                entry = domain["signature"]["entry"]
                if entry is not None and "type" in entry:
                    if entry["type"] == "DOMAIN":
                        domains.append(entry["name"])
            domain_data[fusion_name] = domains
        
        pegasus_output = []
        #cols = ["DriverScore", "Tot/span_reads", "Split_reads", "Chr1", "Chr2", "Strand1", "Strand2", "Gene_Name1", "Gene_Name2", "Gene_Breakpoint1", "Gene_Breakpoint2", "Kinase_info", "Reading_Frame_Info", "Protein_Sequence", "Conserved_Domain1", "Lost_Domain1", "Conserved_Domain2", "Lost_Domain2"]
        cols = ["DriverScore", "Chr1", "Chr2", "Strand1", "Strand2", "Gene_Name1", "Gene_Name2", "Kinase_info", "Reading_Frame_Info"]
        for index, row in df.iterrows():
            if "DriverScore" not in row:
                continue
            break1 = row["Gene_Breakpoint1"]
            break2 = row["Gene_Breakpoint2"]
            gene1 = row["Gene_Name1"]
            gene2 = row["Gene_Name2"]
            
            if row["Strand1"] == "-":
                break1 = break1 - 1
            if row["Strand2"] == "-":
                break2 = break2 - 1     
            out = {}
            out["spanning_reads"] = row["Tot/span_reads"]
            out["discordant_reads"] = row["Split_reads"] 
            for col in cols:
                out[col] = row[col]
            frame_info = out["Reading_Frame_Info"].rstrip()
            kinase_info = out["Kinase_info"].rstrip()
            fusion_name = gene1 + "-" + gene2 + ":" + str(break1) + "-" + str(break2)
            out["Transcript1"] = trans[fusion_name]["trans1"]
            out["Transcript2"] = trans[fusion_name]["trans2"]
            gene1_domains = []
            gene2_domains = []
            fused_domains = []
            if fusion_name in domain_data:
                fused_domains = domain_data[fusion_name]
            if gene1 in domain_data:
                gene1_domains = domain_data[gene1]
            if gene2 in domain_data:
                gene2_domains = domain_data[gene2]
            found_kinase = False
            for domain in fused_domains:
                if "kinase" in domain:
                    found_kinase = True
            out["domain1"] = list(set(gene1_domains))
            out["domain2"] = list(set(gene2_domains))
            out["fusion_domains"] = list(set(fused_domains))
            if frame_info != "InFrame":
                out["status"] = {"code" : "FAILED", "reason" : "Frameshift"}
            elif not found_kinase:
                out["status"] = {"code" : "FAILED", "reason" : "No kinase"}
            else:
                out["status"] = {"code" : "PASSED"}   
            
            pegasus_output.append(out)
        
        with open(out_json_file, 'w') as outfile:
            json.dump(pegasus_output, outfile, indent = 4, sort_keys=True)
        print("Done. The output file is " + out_json_file + "\n")
        

    def filter(self, in_file, onco_gene_file, not_onco_gene_file):
        discordant_reads_cutoff = self.config_data["discordant_reads_cutoff"]
        spanning_reads_cutoff = self.config_data["spanning_reads_cutoff"]
        gene_list_file = self.config_data["match_gene_list_file"]
        self.gene_list = {}
        out_onco = open(onco_gene_file, "w")
        out_not_onco = open(not_onco_gene_file, "w")
        with open(gene_list_file) as f:
            lines = f.readlines()
        for line in lines:
            line = line.rstrip('\n')
            if line != '':
                self.gene_list[line] = ''
        df = pandas.read_table(in_file, sep='\t', index_col=False)
        sep = "\t"
        out_onco.write(sep.join(df.columns) + "\n")
        out_not_onco.write(sep.join(df.columns) + "\n")
        for index, row in df.iterrows():
            gene5p = row["5p_symbol"]
            gene3p = row["3p_symbol"]
            discordant_reads = row["discordant_reads"]
            span_reads = row["span_reads"]
            if discordant_reads >= discordant_reads_cutoff and discordant_reads >= spanning_reads_cutoff and (gene5p in self.gene_list or gene3p in self.gene_list):
                out_onco.write(sep.join(map(str, row)) + "\n")
            else:
                out_not_onco.write(sep.join(map(str, row)) + "\n")
        out_onco.close()
        out_not_onco.close()
    
    def get_gene_data(self, gtf_file):
        gene_data = {}
        cwd = os.getcwd()
        with open(gtf_file) as f:
            lines = f.readlines()
        for line in lines:
            line = line.rstrip('\n')
            fields = line.split("\t")
            chromosome = fields[0]            
            strand = fields[6]
            description = fields[8]
            m = re.search('gene_id "(.*?)";.*gene_name "(.*?)";', description)
            if (m != None):
                accession = m.group(1)
                gene = m.group(2)
                key = chromosome + ':' + gene
                #print(key+'\n')
                gene_data[key] = {"strand" : strand, "accession" : accession}
        return gene_data

    def append_log(self, message):
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        with open(self.out_log, "a") as logf:
            logf.write(st + ' ' + message + "\n")

    def to_pegasus_file(self, in_file, out_file, gtf_file):
        gene_data = self.get_gene_data(gtf_file)
        cols = {};
        ouf = open(out_file, "w")
        ouf.write("#5p_symbol\t5p_ensembl\t5p_strand\t5p_chr\t3p_symbol\t3p_ensembl\t3p_strand\t3p_chr\t5p_start\t5p_end\t3p_start\t3p_end\tsplit_reads\ttot_reads\n")
        df = pandas.read_table(in_file, sep='\t', index_col=False)
        for index, fields in df.iterrows():
            gene1 = fields["5p_symbol"]
            chr1 = str(fields["5p_chr"])
            if chr1[0:3] != 'chr':
                chr1 = 'chr' + chr1
            gene2 = fields["3p_symbol"]
            chr2 = str(fields["3p_chr"])
            if chr2[0:3] != 'chr':
                chr2 = 'chr' + chr2
            start1 = fields["5p_position"]
            end1 = start1
                    
            start2 = fields["3p_position"]
            end2 = start2
                    
            key1 = chr1 + ':' + gene1
            key2 = chr2 + ':' + gene2
            if key1 not in gene_data:
                print(key1 + " not exists\n")
                continue
            if key2 not in gene_data:
                print(key2 + " not exists\n")
                continue 
            strand1 = gene_data[key1]["strand"]
            accession1 = gene_data[key1]["accession"]
            strand2 = gene_data[key2]["strand"]
            accession2 = gene_data[key2]["accession"]
            if (strand1 == "-"):
                start1 = int(start1) + 1
                end1 = int(end1) + 1
            if (strand2 == "-"):
                start2 = int(start2) + 1
                end2 = int(end2) + 1                
            discordant_reads = fields["discordant_reads"]
            span_reads = fields["span_reads"]
            ouf.write(gene1 + "\t" + accession1 + "\t" + strand1 + "\t" + chr1[3:len(chr1)] + "\t" + gene2 + "\t" + accession2 + "\t" + strand2 + "\t" + chr2[3:len(chr2)] + "\t" + str(start1) + "\t" + str(end1) + "\t" + str(start2) + "\t" + str(end2) + "\t" + str(discordant_reads) + "\t" + str(span_reads) + "\n")

        
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option( '-c' , dest = 'config_file' ,
        default = '' ,
        help = 'config JSON file' )
    parser.add_option( '-i' , dest = 'in_file' ,
        default = '' ,
        help = 'the file containing line-separated gene fusion file' )
    parser.add_option( '-o' , dest = 'out_path' ,
        default = '' ,
        help = 'output folder' )
    parser.add_option( '-a' , dest = 'annotation' ,
        default = 'refseq' ,
        help = 'annotation: refseq or ensembl' )
    parser.add_option('-f', dest = 'forced', action="store_true",
        help = "force to run regardless if output file already exists")
    
    (options,args) = parser.parse_args()
    
    in_file = options.in_file
    out_path = options.out_path
    config_file = options.config_file
    annotation = options.annotation
    forced = options.forced
    
    if in_file == '' or  config_file == '' or out_path == '':
        parser.print_help()
        exit()

    if annotation != 'refseq' and annotation != 'ensembl':
        print("annotation must be either refseq or ensembl, your input: " + annotation)
        exit()
    if os.path.exists(out_path) and forced:
        shutil.rmtree(out_path)
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    try:
        fa = FusionAMOI(config_file, out_path, annotation)
        fa.run_pegasus(in_file, forced)
    except Exception, e:
        with open(fa.err_log, "a") as logf:
            logf.write(str(e) + "\n")
        print "Unexpected error: " + str(e)
        raise 

    

    
