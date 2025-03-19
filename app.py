import streamlit as st
import pandas as pd
import tempfile
import zipfile
import altair as alt
import umap
import shutil
from KZ import KZ_Pipeline
import subprocess
import webbrowser
import os

######################################################################################################################
def fastqc_report():
    st.markdown("## Upload a FastQ file for fastQC report generation")
              
    default_dir = os.getcwd()
    folder_path = st.text_input("Enter the path to your data:", value=default_dir)
    
    if os.path.isdir(folder_path):
        # List all files in the directory with the specific extensions
        files = [f for f in os.listdir(folder_path) if f.endswith('.fq') or f.endswith('.fastq')]
        if files:
            selected_file = st.selectbox("Select a file", files)
            st.write("You selected:", selected_file)
        else:
            st.write("No .fq or .fastq files found in this directory.")
    else:
        st.error("The specified path is not a valid directory.")
        
    if st.button("Generate fastQC report"):
        if not folder_path.endswith('/'):
            folder_path = folder_path+"/"
        subprocess.run(["fastqc", folder_path+selected_file, "-o", "tmp/", "-f", "fastq", '--nano'], capture_output=True)
        
        report_file = os.path.join("tmp/", selected_file.replace('.fastq', '_fastqc.html').replace('.fq', '_fastqc.html'))
        
        if os.path.exists(report_file):
            webbrowser.open("file://"+os.getcwd()+"/tmp/"+ selected_file.replace('.fastq', '_fastqc.html').replace('.fq', '_fastqc.html'))
        else:
            st.error("FastQC report could not be generated.")

######################################################################################################################

# Custom function to handle file upload
def file_uploader():
    st.markdown("## Upload a FastQ file")
    st.markdown("Upload your fastq file to map to a reference and make a concessus sequence. If you have many small fastqs in one run, you may want to `cat *.fastq > new.fastq`. These will be added to a list for nextstrain analysis.")
    
    # File uploader
    file = st.file_uploader("Choose a fastq file", type=["fastq","fq"])     
    
    # Reference Selector
    reference = st.selectbox("Select Reference", ['','CCHF','TBEV'])

    if reference != '':

        run = KZ_Pipeline()
        run.set_reference(reference)

        # Get already existing metadata
        host_meta = run.metadata['host'].unique().tolist() + run.ncbidata['host'].unique().tolist() # Get uniques from the two dataframes
        host_meta = [x for x in host_meta if str(x).strip() != '' and str(x) != 'nan'] # Remove blanks and nan
        host_meta = sorted(list(set(host_meta)), key=lambda x: str(x)) + ['Add New'] # Remove duplicates

        source_meta = run.metadata['isolation_source'].unique().tolist() + run.ncbidata['isolation_source'].unique().tolist()
        source_meta = [x for x in source_meta if str(x).strip() != '' and str(x) != 'nan']
        source_meta = sorted(list(set(source_meta)), key=lambda x: str(x)) + ['Add New']

        country_meta = run.metadata['country'].unique().tolist() + run.ncbidata['country'].unique().tolist()
        country_meta = [x for x in country_meta if str(x).strip() != '' and str(x) != 'nan']
        country_meta = sorted(list(set(country_meta)), key=lambda x: str(x)) + ['Add New']

        # HTML Select and Input boxes for metadata
        st.markdown("#### Meta Data")
        col1, col2 = st.columns(2)
        with col1:
            date = st.date_input("Date Collected")
        with col2:
            name = st.text_input("Run Name", help="This will be the assembly name and viewable in nextstrain")

        host_select = st.selectbox("Host", host_meta)
        if host_select == 'Add New':
            host_custom = st.text_input("Custom Host")

        country_select = st.selectbox("Country", country_meta)
        if country_select == 'Add New':
            country_custom = st.text_input("Custom Country")

        source_select = st.selectbox("Isolation Source", source_meta)
        if source_select == 'Add New':
            source_custom = st.text_input("Custom Source")
            
        if st.button("Submit"):
            # Work out the select box values
            host = host_custom if host_select == 'Add New' else host_select
            country = country_custom if country_select == 'Add New' else country_select
            source = source_custom if source_select == 'Add New' else source_select

            metadata_input = {
                'name': name,
                'date': date,
                'country': country,
                'isolation_source': source,
                'host': host,
            }

            # Upload the file to temp directory
            st.write("Uploading File...")
            temp_file = tempfile.NamedTemporaryFile(delete=False)
            temp_file.write(file.read())

            st.write("Assembling and Creating Consensus...")
            average_coverage, total_seqs, reads_mapped, reads_unmapped, quality, average_snp_quality, coverage_data = run.make_assembly(temp_file.name, metadata_input)
            st.write('Done!')

            temp_file.close()
            #st.session_state.step = 2
            st.write('Total Sequences: ' + str(total_seqs))
            st.write('Sequences Mapped: ' + str(reads_mapped))
            st.write('Sequences Unampped: ' + str(reads_unmapped))
            st.write('Average Quality: ' + str(quality))
            st.write('Average Coverage: ' + str(average_coverage))
            st.write('Average SNP Quality: ' + str(average_snp_quality))
            st.write('Generating Coverage Plot...' )
                        
            plot_df = pd.DataFrame([x.split('\t') for x in coverage_data], columns=['Ref', 'Pos', 'Depth'])
            plot_df['Depth'] = pd.to_numeric(plot_df['Depth'])
            plot_df['Pos'] = pd.to_numeric(plot_df['Pos'])
            
            chart = alt.Chart(plot_df).mark_line().encode(
            x=alt.X('Pos:Q', scale=alt.Scale(zero=False), axis=alt.Axis(title='Position')),
            y=alt.Y('Depth:Q', scale=alt.Scale(zero=False), axis=alt.Axis(title='Depth')),
            tooltip=['Pos', 'Depth']
            ).properties(
                width=800,
                height=600,
                title='Coverage Plot'
            ).interactive()
        
            # Display the chart in Streamlit
            st.altair_chart(chart)

######################################################################################################################
# Run Nextstrain
def run_nextstrain():
    st.markdown("## Run Nextstrain")

    # Select dropdown for CCHF and TBEV
    reference = st.selectbox("Select Reference", ['CCHF','TBEV'])
    # Start Pipeline and insert reference type
    run = KZ_Pipeline()
    run.set_reference(reference)

    # Build out editable table
    if st.checkbox("Use NCBI Data"):
        df1 = run.metadata
        df2 = run.ncbidata
        df = pd.concat([df1,df2]).reset_index(drop=True)
    else:
        df = run.metadata

    if len(df) > 0:
        df['Include'] = True
        df['Delete'] = False  # Column to select records for deletion
        edited_df = df[['Delete', 'Include', 'name', 'date', 'length', 'country', 'isolation_source', 'host', 'desc']]
        edited_df = st.data_editor(df, disabled=('name', 'date', 'length', 'country', 'isolation_source', 'host', 'desc','seq'))

        # Delete selected records
        if st.button("Delete Selected Records"):
            names_to_keep = edited_df[edited_df['Delete'] == False]['name'].tolist()
            df = df[df['name'].isin(names_to_keep)]
            run.metadata = df  # Update the metadata DataFrame

            # Write changes back to the file
            df.to_csv(run.metadata_file, sep='\t', index=False)
            st.success("Selected records deleted and metadata updated.")
            if st.button("Click to update dataframe"):
                edited_df['Delete']=False

        # Submit Changes. Build New DF and Run Nextstrain
        if st.button("Submit"):
            include_indexes = edited_df.loc[edited_df['Include'] == True].index.values.tolist()
            new_df = df.iloc[include_indexes]

            if len(new_df) > 3:
                st.write('Creating Augur Alignment...')
                run.create_msa(new_df[['name', 'desc', 'seq']])
                st.write('Building trees and node data...')
                run.process_augur(new_df[['name', 'date', 'country', 'isolation_source', 'host']])
                st.write('Passing data to nextstrain...')
                run.view_nextstrain()
            else:
                st.write('You need more than 3 uploaded/selected to run nextstrain.')

    else:
        st.write('No files added.')


######################################################################################################################
# Run Embedding
def run_embedding():
    st.markdown("## Run Embedding")

    # Select dropdown for CCHF and TBEV
    reference = st.selectbox("Select Reference", ['CCHF','TBEV'])
    # Start Pipeline and insert reference type
    run = KZ_Pipeline()
    run.set_reference(reference)

    # Get our dataframe but it in an editable table
    edited_df = run.metadata
    edited_df['Include'] = True
    edited_df = st.data_editor(edited_df, disabled=('name','date','length','country','isolation_source','host','desc','seq'))

    # Selector controls for graph
    col1, col2 = st.columns(2)
    with col1:
        color = st.selectbox("Legend Coloring",
            [
                'type',
                'country',
                'host',
                'length'
            ])

    # Submit Changes. Build new DF and Run Embedding
    if st.button("Submit"):
        include_indexs = edited_df.loc[edited_df['Include']==True].index.values.tolist()
        new_df = run.metadata.iloc[include_indexs]
        new_df['type'] = 'Project Created'

        args = {
            'klen':         11,
            'scale':        1,
            'abundance':    False,
        }

        @st.cache_data
        def compute_umap_embedding(data):
            reducer = umap.UMAP()
            embedding = reducer.fit_transform(data)
            return embedding

        labels, data = run.process_embedding(new_df, args)
        embedding = compute_umap_embedding(data)

        # Altair dynamic plot
        plot_df = pd.DataFrame({
            "x": embedding[:, 0], 
            "y": embedding[:, 1], 
            "label": labels[color],
            'Name': labels['name'], 
            'Description': labels['desc'],
            'Length': labels['length'],
            'Date': labels['date'],
            'Country': labels['country'],
            'Host': labels['host'],
            })
        chart = alt.Chart(plot_df).mark_circle().encode(
            x="x",
            y="y",
            color="label",
            tooltip=['Name','Description','Length','Date','Country','Host']
        ).properties(
            width=800,
            height=600
        ).configure_mark(
            size=0.5
        )

        st.altair_chart(chart)


######################################################################################################################
# Export
def run_export():

    st.write('Click to export all CCHF NCBI and submitted records')
    if st.button("Export all CCHF"):
        run1 = run = KZ_Pipeline()
        run1.clean() # Remove everything in the tmp directory
        run1.set_reference('CCHF')
        files = run1.create_export_tmp()
        archive_name = f'tmp/KZ_archive_{run1.time}.zip'

        with zipfile.ZipFile(archive_name,'w') as zipf:
            for file in files:
                zipf.write(file, arcname=file.split('/')[-1])
    
        with open(archive_name, 'rb') as f:
            st.download_button("Download Zip", f, file_name="CCHF_archive.zip")
    
    st.write('Click to export all TBEV NCBI and submitted records')
    if st.button("Export all TBEV"):
        run2 = run = KZ_Pipeline()
        run2.clean() # Remove everything in the tmp directory
        run2.set_reference('TBEV')
        files = run2.create_export_tmp()
        archive_name = f'tmp/KZ_archive_{run2.time}.zip'
        
        with zipfile.ZipFile(archive_name,'w') as zipf:
            for file in files:
                zipf.write(file, arcname=file.split('/')[-1])
    
        with open(archive_name, 'rb') as f:
            st.download_button("Download Zip", f, file_name="TBEV_archive.zip")
        
    
    
    st.write('Click to export select uploaded records. Folder will contain the assembled consensus sequences in fasta format')
    df1 = pd.read_table('res/TBEV_metadata.tsv')
    df2 = pd.read_table('res/CCHF_metadata.tsv')
    
    # Concatenate the DataFrames
    df = pd.concat([df1,df2]).reset_index(drop=True)
    edited_df = df[['name', 'date', 'length', 'country', 'isolation_source', 'host', 'desc', 'seq', 'Include']]
    edited_df['Include']=False
    edited_df = st.data_editor(edited_df, disabled=('name', 'date', 'length', 'country', 'isolation_source', 'host', 'desc','seq'))
    
    if st.button("Export select sequences"):
        selected_df = edited_df[edited_df['Include']]
        temp_dir = 'tmp_export'
        os.makedirs(temp_dir, exist_ok=True)
        archive_name = f'{temp_dir}/KZ_selected_records.zip'
        
        with zipfile.ZipFile(archive_name, 'w') as zipf:
            # Create a FASTA file for each selected record
            for _, row in selected_df.iterrows():
                fasta_content = f">{row['name']}\n{row['seq']}\n"
                fasta_filename = f"{temp_dir}/{row['name']}.fasta"
                
                with open(fasta_filename, 'w') as fasta_file:
                    fasta_file.write(fasta_content)

                zipf.write(fasta_filename, os.path.basename(fasta_filename))
                os.remove(fasta_filename)
                
        # Provide download link for the zip file
        with open(archive_name, 'rb') as f:
            st.download_button('Download ZIP', f.read(), file_name=os.path.basename(archive_name))
        shutil.rmtree(temp_dir)


######################################################################################################################
def main():
    if "step" not in st.session_state:
        st.session_state.step = 1

    st.title("KZ analysis")
    selected_page = st.sidebar.radio("Select a Page", ["Generate fastQC report", "Upload Fastq", "Run Nextstrain", "Run Embedding", "Export Results"])

    content = st.container()

    if selected_page == 'Run Nextstrain':
        with content:
            run_nextstrain()
    elif selected_page == 'Run Embedding':
        with content:
            run_embedding()
    elif selected_page == 'Export Results':
        with content:
            run_export()
    elif selected_page == "Generate fastQC report":
        with content:
            fastqc_report()
    else:
        with content:
            file_uploader()
    
    # Display the page accordint to the step in analysis
    #if st.session_state.step == 1:
    #    with content:
    #        file_uploader()
    #elif st.session_state.step == 2:
    #    with content:
    #        input_form()
    #else:
    #    with content:
    #        file_uploader()

if __name__ == "__main__":

    main()
