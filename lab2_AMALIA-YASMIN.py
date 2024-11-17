#LAB 2 : AMALIA YASMIN BT ABDUL AZIZ (A22EC0138)

import requests
import pandas as pd
import networkx as nx
import streamlit as st
import matplotlib.pyplot as plt

# Function to retrieve BIOGRID database
def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "aa7fc91153f40ea8f17b38a1a2b9aaa0",  
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,  
        "includeInteractors": True
    }
    
    try:
        response = requests.get(biogrid_url, params=params)
        response.raise_for_status()  
        network_data = response.json()
        
        if not network_data:
            st.warning(f"No data found for the protein ID '{target_protein}' in BioGRID.")
            return pd.DataFrame()
        
        df = pd.DataFrame.from_dict(network_data, orient='index')
        
        if 'OFFICIAL_SYMBOL_A' in df.columns and 'OFFICIAL_SYMBOL_B' in df.columns:
            df['OFFICIAL_SYMBOL_A'] = df['OFFICIAL_SYMBOL_A'].str.upper()
            df['OFFICIAL_SYMBOL_B'] = df['OFFICIAL_SYMBOL_B'].str.upper()
            return df[['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B']]
        else:
            st.warning("The data structure from BioGRID does not contain expected columns.")
            return pd.DataFrame()
        
    except requests.exceptions.RequestException as e:
        st.error(f"Error retrieving data from BioGRID: {e}")
        return pd.DataFrame()

# Function to retrieve STRING database
def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606
    }
    
    try:
        response = requests.get(string_url, params=params)
        response.raise_for_status()  
        network_data = response.json()
        
        if not network_data:
            st.warning(f"No data found for the protein ID '{target_protein}' in STRING.")
            return pd.DataFrame()
        
        df = pd.json_normalize(network_data)
        return df[['preferredName_A', 'preferredName_B']]
    
    except requests.exceptions.RequestException as e:
        st.error(f"Error retrieving data from STRING: {e}")
        return pd.DataFrame()

# Function that use the PPI (in DataFrame) and create the network using the networkx library
def generate_network(dataframe):
    network_graph = nx.from_pandas_edgelist(dataframe, source=dataframe.columns[0], target=dataframe.columns[1])
    return network_graph

# Function that retrieve all the network centrality measures
def get_centralities(network_graph):
    degree_centrality = nx.degree_centrality(network_graph)
    closeness_centrality = nx.closeness_centrality(network_graph)
    betweenness_centrality = nx.betweenness_centrality(network_graph)
    eigenvector_centrality = nx.eigenvector_centrality(network_graph, max_iter=1000)
    pagerank_centrality = nx.pagerank(network_graph)
    
    centralities = {
        'Degree Centrality': degree_centrality,
        'Closeness Centrality': closeness_centrality,
        'Betweenness Centrality': betweenness_centrality,
        'Eigenvector Centrality': eigenvector_centrality,
        'PageRank Centrality': pagerank_centrality
    }
    return centralities

# Streamlit App 
st.title("Human Protein-Protein Interaction (PPI) Network")
target_protein = st.text_input("Enter Protein ID (e.g., TP53):")
database_choice = st.selectbox("Select Database", ("BioGRID", "STRING"))

if st.button("Retrieve PPI Data"):
    # Retrieve data based on user database choice (BIOGRID or STRING)
    if database_choice == "BioGRID":
        ppi_data = retrieve_ppi_biogrid(target_protein)
    else:
        ppi_data = retrieve_ppi_string(target_protein)
    
    # Check if data is empty
    if ppi_data.empty:
        st.warning("No PPI data found. Please check the protein ID and try again.")
    else:
        # Generate network graph
        network_graph = generate_network(ppi_data)
        
        # Display PPI data information in the first column
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("PPI Data Information")
            st.dataframe(ppi_data)
            st.write("Number of nodes:", network_graph.number_of_nodes())
            st.write("Number of edges:", network_graph.number_of_edges())
            
            # Visualize the network
            if network_graph.number_of_nodes() == 0:
                st.warning("No interactions found for the given protein.")
            else:
                fig, ax = plt.subplots(figsize=(10, 10))  # Adjust size for larger graphs
                layout = nx.spring_layout(network_graph, seed=42, k=0.15)  # Adjust 'k' for better spacing
                nx.draw(network_graph, layout, ax=ax, with_labels=True, node_size=30, node_color="skyblue", font_size=8)
                st.pyplot(fig)

        # Display centrality measures in the second column
        with col2:
            st.subheader("Centrality Measures")
            centralities = get_centralities(network_graph)
            for centrality_name, centrality_values in centralities.items():
                sorted_nodes = sorted(centrality_values.items(), key=lambda x: -x[1])[:5]  # Top 5 nodes for each centrality
                st.write(f"{centrality_name}:")
                for node, centrality in sorted_nodes:
                    st.write(f"{node}: {centrality:.4f}")
