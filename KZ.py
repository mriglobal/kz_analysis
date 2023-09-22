# -*- coding: utf-8 -*-

import os
import re
import shutil
import pandas as pd
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sourmash


class KZ_Pipeline():

    ######################################################################################################################
    ## ---- CLASS START
    def __init__(self):
        self.time = datetime.now().strftime('%Y-%m-%d_%H:%M')
        self.threads = os.cpu_count()

        # Create needed folders
        folders = ['res','tmp']
        for folder in folders:
            self.create_folder(folder)


    ######################################################################################################################
    ## ---- SET REFERENCES AND SET METADATA
    def set_reference(self, reference):
        self.reference = reference
        
        # Set references and metadata for TBEV
        if self.reference == 'TBEV':
            self.ref_file = 'res/TBEV_reference.fasta'
            self.ref_file_gb = 'res/TBEV_reference.gb'
            self.mmi_file = 'res/TBEV_reference.mmi'
            self.metadata_file = 'res/TBEV_metadata.tsv'
            self.ncbidata_file = 'res/TBEV_NCBI_metadata.tsv'
            self.export_fasta = 'tmp/TBEV.fasta'
            self.export_tsv = 'tmp/TBEV_metadata.tsv'
        # Set references and metadata for CCHF
        elif self.reference == 'CCHF':
            self.ref_file = 'res/CCHF_reference.fasta'
            self.ref_file_gb = 'res/CCHF_reference.gb'
            self.mmi_file = 'res/CCHF_reference.mmi'
            self.metadata_file = 'res/CCHF_metadata.tsv'
            self.ncbidata_file = 'res/CCHF_NCBI_metadata.tsv'
            self.export_fasta = 'tmp/CCHF.fasta'
            self.export_tsv = 'tmp/CCHF_metadata.tsv'
        else:
            raise Exception("Need the reference to be either TBEF or CCHF")
        
        # if you haven't indexed the reference, do so
        if not os.path.exists(self.mmi_file):
            os.system(f"minimap2 -d {self.mmi_file} {self.ref_file}")

        # Load Metadata
        if os.path.exists(self.metadata_file):
            self.metadata = pd.read_table(self.metadata_file)
        else:
            self.metadata = pd.DataFrame({
                'name':             [],
                'length':           [],
                'date':             [],
                'country':          [],
                'isolation_source': [],
                'host':             [],
                'desc':             [],
            })

        self.ncbidata = pd.read_table(self.ncbidata_file)

    ######################################################################################################################
    ## ---- HELPER FUNCTIONS
    def create_folder(self, folder):
        if not os.path.exists(folder):
            os.makedirs(folder)

    def clean(self):
        shutil.rmtree('tmp')
        self.create_folder('tmp')

    def get_new_seqid(self, seqid, current_seqs):
        if seqid in current_seqs:
            split = seqid.split('.')
            if len(split) == 1 or split[-1].isdigit() == False:
                new_seqid = f'{seqid}.1'
                return new_seqid, current_seqs
            else:
                ext = split[-1]
                # Is an integer after the dot
                new_seqid = '.'.join(split[:-1])
                new_seqid = str(new_seqid) + '.' + str(int(ext) + 1)
                if new_seqid in current_seqs:
                    return self.get_new_seqid(new_seqid, current_seqs)
                else:
                    current_seqs.append(new_seqid)
                    return new_seqid, current_seqs
        else:
            return seqid, current_seqs
    
    def seqs_from_df(self, input_df):
        seqs = []
        for index, row in input_df.iterrows():
            seqs.append(
                SeqRecord(
                    Seq(str(row['seq'])), 
                    id=str(row['name']),
                    description='',
                    )
                )
        return seqs
    
    def create_export_tmp(self):
        fasta = self.seqs_from_df(self.metadata[['name','seq']])
        SeqIO.write(fasta, self.export_fasta, 'fasta')

        export_metadata = self.metadata
        del export_metadata['seq']
        export_metadata.to_csv(self.export_tsv, sep='\t', index=False)

        return [self.export_fasta, self.export_tsv]

    ######################################################################################################################
    ## ---- MAKE ASSEMBLY FROM UPLOADED FASTQ. 
    ## ---- WRITE OUT METADATA AND FINAL ASSEMBLY TO METADATA TABLE
    def make_assembly(self, input_fastq, metadata_input):
        self.clean()

        sam_file = f"tmp/to_ref.sam"
        bam_file = f"tmp/to_ref_sorted.bam"
        vcf_file = f"tmp/calls.vcf.gz"

        #### Create Assembly
        # Minimap2 to map the sequences to the input indexed sam
        os.system(f"minimap2 -ax map-ont -t {self.threads} {self.mmi_file} {input_fastq} > {sam_file}")
        # samtools to convert sam to bam and sort
        os.system(f"samtools view -bS {sam_file} | samtools sort -o {bam_file}")
        # bcf tools to convert the sorted bam to a vcf.gz
        os.system(f"bcftools mpileup -Ou -f {self.ref_file} {bam_file} | bcftools call -mv -Oz -o {vcf_file}")
        # Make a .tbi file from the vcf
        os.system(f"tabix -p vcf {vcf_file}")
        # get the consensus assembly from the vcf
        os.system(f"bcftools consensus -f {self.ref_file} {vcf_file} > tmp/consensus.fasta")

        if self.reference == 'CCHF':
            # eliminate all but the shortest fasta entry (L Segment)
            records = list(SeqIO.parse('tmp/consensus.fasta', "fasta"))
            # Find the shortest entry
            shortest_record = min(records, key=lambda x: len(x.seq))
            # Write the shortest entry to the FASTA file
            with open('tmp/consensus.fasta', "w") as output_handle:
                SeqIO.write(shortest_record, output_handle, "fasta")

        # Add assembly to a table, make sure they are unique
        if len(self.metadata) > 0 or len(self.ncbidata) > 0:
            current_seqs = self.metadata.name.to_list() + self.ncbidata.name.to_list()
        else:
            current_seqs = []

        # Loop through each sequence and add it to our metadata table
        for seq in [s for s in SeqIO.parse('tmp/consensus.fasta','fasta')]:
            # Clean up the name
            clean_name = str(re.sub(r'[^a-zA-Z0-9\s_]', '', metadata_input['name'])).replace(' ','_')
            if clean_name == '':
                clean_name = seq.id
            
            # Get a unique seq.id
            new_id, current_seqs = self.get_new_seqid(clean_name, current_seqs)

            this_metadata = pd.DataFrame({
                    'name':                 [new_id],
                    'date':                 [metadata_input['date']],
                    'length':               [len(seq.seq)],
                    'country':              [metadata_input['country']],
                    'isolation_source':     [metadata_input['isolation_source']],
                    'host':                 [metadata_input['host']],
                    'desc':                 [seq.description.replace(seq.id,'').strip()],
                    'seq':                  [seq.seq]
                })

            # Add to metadata table
            self.metadata = pd.concat([self.metadata, this_metadata])
        
        self.metadata.to_csv(self.metadata_file, sep='\t', index=False)


    ######################################################################################################################
    ## ---- CREATE OUR MSA FASTA
    ## ---- CONVERT OUR METADATA TABLE TO FASTA FILE AND RUN AUGUR ALIGNMENT
    def create_msa(self, input_df):
        self.clean()

        tmp_alignment = 'tmp/combined_alignment.fasta'
        seqs = self.seqs_from_df(input_df)
            
        # Write out to temp fasta file
        SeqIO.write(seqs, tmp_alignment, 'fasta')
        
        os.system(f"augur align \
                  --sequences {tmp_alignment} \
                  --reference-sequence {self.ref_file_gb} \
                  --fill-gaps \
                  --output tmp/msa.fasta \
                  --nthreads {self.threads}"
                  )
    
    ######################################################################################################################
    ## ---- RUN AUGUR PIPELINE
    def process_augur(self, input_df):

        metadata = 'tmp/metadata.tsv'
        alignment = 'tmp/msa.fasta'
        tree = 'tmp/augur_output_tree.nwk'
        refine = 'tmp/augur_output_tree_refined.nwk'
        node_data = 'tmp/augur_refined_node.json'
        ancestral = 'tmp/augur_ancestral.json'
        translate = 'tmp/augur_muts.json'
        traits = 'tmp/augur_traits.json'
        config = 'res/auspice_config.json'
        auspice = 'tmp/augur_auspice.json'

        # Write out our metadata
        input_df.to_csv(metadata, sep='\t', index=False)
    
        # Build Tree
        os.system(f"augur tree \
                  --alignment {alignment} \
                  --method iqtree \
                  --output {tree} \
                  --nthreads {self.threads}"
                  )
        # Augur Refine
        os.system(f"augur refine \
                  --tree {tree} \
                  --alignment {alignment} \
                  --metadata {metadata} \
                  --output-tree {refine} \
                  --output-node-data {node_data}"
                  )
        # Augur Ancestral
        os.system(f"augur ancestral \
                  --tree {refine} \
                  --alignment {alignment} \
                  --inference joint \
                  --output-node-data {ancestral}"
                  )
        os.system(f"augur translate \
                  --tree {refine} \
                  --ancestral-sequences {ancestral} \
                  --reference-sequence {self.ref_file_gb} \
                  --output-node-data {translate}"
        )
        # Augur Traits
        os.system(f"augur traits \
                  --tree {refine} \
                  --metadata {metadata} \
                  --columns country host isolation_source \
                  --confidence \
                  --output-node-data {traits}"
                  )
        # Augur Export
        os.system(f"augur export v2 \
                  --tree {refine} \
                  --metadata {metadata} \
                  --node-data {node_data} {traits} {ancestral} {translate} \
                  --auspice-config {config} \
                  --output {auspice}"
                  )
        

    ######################################################################################################################
    ## ---- CREATE EMBEDDING
    def process_embedding(self,input_df, args):

        # Get labels for uploaded data
        input_labels = input_df[['name','length','date','country','isolation_source','host','desc','type']]

        # Get labels for genbank data
        genbank_labels = self.ncbidata[['name','length','date','country','isolation_source','host','desc']]
        genbank_labels['type'] = 'Genbank'
        labels = pd.concat([input_labels,genbank_labels])

        # Get our sequences
        seqs = self.seqs_from_df(self.ncbidata)
        seqs = self.seqs_from_df(input_df) + seqs

        sketches = []
        # Add our sequences
        for s in seqs:
            mh = sourmash.MinHash(0, ksize=args['klen'], scaled=args['scale'], track_abundance=args['abundance'])
            mh.add_sequence(str(s.seq),force=True)
            sketches.append(mh)
        
        ignore_abund = not args['abundance']
        sim_matrix = sourmash.compare.compare_all_pairs(sketches, ignore_abund)
        
        return labels, sim_matrix

    ######################################################################################################################
    ## ---- VIEW NEXTSTRAIN
    def view_nextstrain(self):
        os.system(f"nextstrain view tmp/augur_auspice.json")


###############################################################################
if __name__ == "__main__":

    print('Run steamlit application to work with KZ nextstrain workflow\nconda activate KZ_augur\nstreamlit run app.py')