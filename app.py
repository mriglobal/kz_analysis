import streamlit as st
import pandas as pd
import tempfile
import altair as alt
import seaborn as sns
import umap
import matplotlib.pyplot as plt
from KZ import KZ_Pipeline

######################################################################################################################
# Custom function to handle file upload
def file_uploader():
    st.markdown("## Upload a FastQ file")
    st.markdown("Upload your fastq file to map to a reference and make a concessus sequence. If you have many small fastqs in one run, you may want to `cat *.fastq > new.fastq`. These will be added to a list for nextstrain analysis.")
    file = st.file_uploader("Choose a fastq file", type=["fastq","fq"])

    st.markdown("#### Meta Data")
    col1, col2 = st.columns(2)
    with col1:
        reference = st.selectbox("Select Reference", ['CCHF','TBEV'])
        date = st.date_input("Date Collected")
        host = st.text_input("Host")
    with col2:
        name = st.text_input("Run Name", help="This will be the assembly name and viewable in nextstrain")
        country = st.text_input("Country", value="Kazakhstan")
        source = st.text_input("Isolation Source", help="Blood, Soil, etc.")
        

    if st.button("Submit"):
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

        st.write("Assembling and Creating Concenssus...")
        run = KZ_Pipeline()
        run.set_reference(reference)
        run.make_assembly(temp_file.name, metadata_input)
        st.write('Done!')

        temp_file.close()
        #st.session_state.step = 2


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
        df = pd.concat([df1,df2])
    else:
        df = run.metadata

    if len(df) > 0:
        df['Include'] = True
        df = df[['Include','name','date','country','isolation_source','host','desc']]
        edited_df = st.data_editor(df, disabled=('name','date','country','isolation_source','host','desc'))

        # Submit Changes. Build New DF and Run Nextstrain
        if st.button("Submit"):
            include_indexs = edited_df.loc[edited_df['Include']==True].index.values.tolist()
            new_df = run.metadata.iloc[include_indexs]
            #new_df = run.metadata.iloc[0:300]

            if len(new_df) > 3:
                st.spinner('Creating Augur Alignment...')
                run.create_msa(new_df[['name','desc','seq']])
                st.spinner('Building trees and node data...')
                run.process_augur(new_df[['name','date','country','isolation_source','host']])
                st.spinner('Passing data to nextstrain...')
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

    df = run.metadata

    if len(df) > 0:
        df['Include'] = True
        df = df[['Include','name','date','country','isolation_source','host','desc']]
        edited_df = st.data_editor(df, disabled=('name','date','country','isolation_source','host','desc'))

        col1, col2 = st.columns(2)
        with col1:
            color = st.selectbox("Select Legend Variable",
                [
                    'country',
                    'host',
                    'date'
                    'length'
                    'type'
                ])
        with col2:
            shape = st.selectbox("Select Legend Shape Variable",
                [
                    'type',
                    'country',
                    'host'
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

            # Matplotlib graph (image)
            #fig, ax = plt.subplots()
            #scatter = ax.scatter(embedding[:,0], embedding[:,1])
            #st.pyplot(fig)

            # Seaborn (image)
            #plt.figure(figsize=(8,6))
            #sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], palette="viridis") 
            #st.pyplot(plt)

            # tooltips ------ 'name','length','date','country','host','desc','type'

            # Altair dynamic plot
            plot_df = pd.DataFrame({
                "x": embedding[:, 0], 
                "y": embedding[:, 1], 
                "label": labels[color], 
                "shape": labels[shape], 
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
                shape="shape",
                tooltip=['Name','Description','Length','Date','Country','Host']
            ).properties(
                width=800,
                height=600
            ).configure_mark(
                size=0.5
            )

            st.altair_chart(chart)


######################################################################################################################
def main():
    if "step" not in st.session_state:
        st.session_state.step = 1

    st.title("KZ analysis")
    selected_page = st.radio("Select a Page", ["Upload Fastq", "Run Nextstrain", "Run Embedding"])

    content = st.container()

    if selected_page == 'Run Nextstrain':
        with content:
            run_nextstrain()
    elif selected_page == 'Run Embedding':
        with content:
            run_embedding()
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