import streamlit as st
import requests
import re
import sys
import os
from langchain_community.callbacks.streamlit import StreamlitCallbackHandler
from langchain_openai import ChatOpenAI


# Add the project root to the Python path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import api_config

os.environ['STREAMLIT_BROWSER_GATHER_USAGE_STATS'] = 'false'

# Import your function from the test.py file
from test import enhanced_query
# Import retrosynthesis function
from tools.retrosynthesis import run_retrosynthesis
# Import NameToSMILES tool
from tools.name2smiles import NameToSMILES  # Update this with the correct import path

# Initialize LLM for name to SMILES conversion
llm = ChatOpenAI(model="gpt-4o", temperature=0)
# Initialize NameToSMILES tool
name_to_smiles_tool = NameToSMILES()

# Set up the Streamlit page
st.set_page_config(
    page_title="ChemCopilot - Chemistry Assistant",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for styling
st.markdown("""
    <style>
    .main-header {
        font-size: 42px;
        font-weight: bold;
        color: #2e7d32;
        margin-bottom: 0px;
    }
    .sub-header {
        font-size: 20px;
        color: #5c5c5c;
        margin-bottom: 30px;
    }
    .stButton>button {
        background-color: #2e7d32;
        color: white;
        border: none;
        padding: 10px 24px;
        border-radius: 4px;
        font-weight: bold;
    }
    .stButton>button:hover {
        background-color: #005005;
    }
    .tool-card {
        background-color: #f5f5f5;
        padding: 10px;
        border-radius: 5px;
        margin-bottom: 10px;
    }
    .result-area {
        background-color: #f9f9f9;
        padding: 20px;
        border-radius: 8px;
        border-left: 4px solid #2e7d32;
    }
    .reaction-card {
        background-color: #f0f7f0;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 15px;
        border-left: 4px solid #2e7d32;
    }
    .analysis-section {
        margin-top: 30px;
        padding-top: 20px;
        border-top: 1px solid #ddd;
    }
    .query-box {
        background-color: #f5f7fa;
        padding: 15px;
        border-radius: 8px;
        margin-top: 20px;
        border-left: 4px solid #3f51b5;
    }
    </style>
    """, unsafe_allow_html=True)

# Initialize session state variables
if 'query' not in st.session_state:
    st.session_state.query = ""
if 'retro_results' not in st.session_state:
    st.session_state.retro_results = None
if 'selected_rxn_smiles' not in st.session_state:
    st.session_state.selected_rxn_smiles = None
if 'analysis_result' not in st.session_state:
    st.session_state.analysis_result = None
if 'rxn_query' not in st.session_state:
    st.session_state.rxn_query = ""
if 'show_analyze_section' not in st.session_state:
    st.session_state.show_analyze_section = False

# Function to get SMILES from chemical name using NameToSMILES tool
def get_smiles_from_name(name):
    """Get SMILES notation for a chemical name using the NameToSMILES tool"""
    try:
        result = name_to_smiles_tool._run(name)
        if "SMILES:" in result:
            # Extract SMILES from the result
            smiles = result.split("SMILES:")[1].split("\n")[0].strip()
            return smiles
        return None
    except Exception as e:
        st.error(f"Error converting {name} to SMILES: {str(e)}")
        return None

# Function to convert reactant and product names to SMILES
def convert_to_reaction_smiles(reactants, products):
    """Convert reactant and product names to reaction SMILES format using LLM and NameToSMILES tool"""
    reactant_smiles = []
    product_smiles = []
    
    # First try to convert each reactant and product using NameToSMILES tool
    for reactant in reactants:
        smiles = get_smiles_from_name(reactant)
        if smiles:
            reactant_smiles.append(smiles)
        else:
            reactant_smiles.append(reactant)  # Keep original name if conversion fails
    
    for product in products:
        smiles = get_smiles_from_name(product)
        if smiles:
            product_smiles.append(smiles)
        else:
            product_smiles.append(product)  # Keep original name if conversion fails
    
    # If NameToSMILES failed for some compounds, use LLM as backup
    if any(not re.match(r'^[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+$', r) for r in reactant_smiles) or \
       any(not re.match(r'^[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+$', p) for p in product_smiles):
        
        prompt = f"""Convert these chemical names to a reaction SMILES format.

Reactants: {', '.join(reactants)}
Products: {', '.join(products)}

Format the output as reaction SMILES using the format: reactant1.reactant2>>product1.product2
Only output the SMILES, no explanations.
"""
        
        try:
            response = llm.invoke(prompt)
            reaction_smiles = response.content.strip()
            
            # Verify that the output looks like a reaction SMILES
            if ">>" in reaction_smiles:
                return reaction_smiles
        except Exception as e:
            st.error(f"Error converting to reaction SMILES with LLM: {str(e)}")
    
    # Join SMILES directly if all conversions were successful
    reactants_str = '.'.join(reactant_smiles)
    products_str = '.'.join(product_smiles)
    return f"{reactants_str}>>{products_str}"

# Function to extract clean SMILES from reaction data
def extract_reaction_smiles(reaction):
    """Extract or generate clean reaction SMILES from reaction data"""
    # First check if there's a valid reaction_smiles that doesn't have placeholders
    if "reaction_smiles" in reaction and reaction["reaction_smiles"]:
        smiles = reaction["reaction_smiles"]
        # Check if it contains placeholders
        if not re.search(r'\[.*?_SMILES\]|\[.*?\]', smiles) and ">>" in smiles:
            return smiles
    
    # If reaction_smiles is not usable, create from reactants and products
    return convert_to_reaction_smiles(reaction["reactants"], reaction["products"])

# Function to create reaction SMILES from reactants and products
def create_reaction_smiles(reactants, products):
    """Create reaction SMILES from reactants and products, handling both SMILES and name inputs"""
    # Check if inputs are likely SMILES or names
    smiles_pattern = r'^[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+$'
    
    # Check if all reactants and products are SMILES
    all_smiles = all(re.match(smiles_pattern, r) for r in reactants) and all(re.match(smiles_pattern, p) for p in products)
    
    if all_smiles:
        # Join SMILES directly
        reactants_str = '.'.join(reactants)
        products_str = '.'.join(products)
        return f"{reactants_str}>>{products_str}"
    else:
        # Use conversion function
        return convert_to_reaction_smiles(reactants, products)

# Function to handle reaction analysis
def analyze_reaction(reactants, products):
    """Analyze reaction by converting names to SMILES if necessary"""
    reaction_smiles = create_reaction_smiles(reactants, products)
    
    if reaction_smiles:
        st.session_state.selected_rxn_smiles = reaction_smiles
        st.session_state.query = f"Give full information about this rxn {reaction_smiles}"
        st.session_state.show_analyze_section = True
        return True
    else:
        st.error("Could not create reaction SMILES. Please check the reactants and products.")
        return False

# App Header
st.markdown('<p class="main-header">ChemCopilot</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Your Expert Chemistry Assistant</p>', unsafe_allow_html=True)

# Sidebar with tools info and examples
with st.sidebar:
    st.markdown("## Tools Available")
    
    with st.expander("🔍 SMILES2Name"):
        st.markdown("Converts SMILES notation to chemical names.")
        st.markdown("*Example:* `C1=CC=CC=C1` → `Benzene`")
    
    with st.expander("📝 Name2SMILES"):
        st.markdown("Converts chemical names to SMILES notation.")
        st.markdown("*Example:* `Ethanol` → `CCO`")
    
    with st.expander("🧪 FuncGroups"):
        st.markdown("Analyzes functional groups in molecules.")
        st.markdown("*Example:* `C(O)(=O)C1=CC=CC=C1` → `Carboxylic acid, Aromatic ring`")
    
    with st.expander("⚛️ BondChangeAnalyzer"):
        st.markdown("Analyzes bond changes in chemical reactions.")
        st.markdown("*Example:* `CCCl.CC[O-].[Na+]>>CCOCC.[Na+].[Cl-]` → `C-Cl bond broken, C-O bond formed`")
    
    st.markdown("## Example Queries")
    
    if st.button("Search for Flubendiamide"):
        st.session_state.query = "flubendiamide"
        st.session_state.show_analyze_section = False
        st.rerun()
        
    if st.button("Search for 2-amino-5-chloro-3-methyl benzoic acid"):
        st.session_state.query = "2-amino-5-chloro-3-methyl benzoic acid"
        st.session_state.show_analyze_section = False
        st.rerun()
        
    if st.button("Reaction Analysis Example"):
        st.session_state.query = "Give full information about this rxn CCCl.CC[O-].[Na+]>>CCOCC.[Na+].[Cl-]"
        st.session_state.selected_rxn_smiles = "CCCl.CC[O-].[Na+]>>CCOCC.[Na+].[Cl-]"
        st.session_state.show_analyze_section = True
        st.rerun()

# Main content area - integrated into a single page
st.markdown("## Retrosynthesis & Reaction Analysis")

# Search Input Section
compound_name = st.text_input(
    "Enter compound name (IUPAC or common name):",
    value=st.session_state.query if not st.session_state.show_analyze_section else "",
    placeholder="Example: flubendiamide, 2-amino-5-chloro-3-methyl benzoic acid",
    key="search_compound_name"
)

search_button = st.button("Search Retrosynthesis", key="search_retro_button")

# Process Search
if search_button and compound_name:
    st.session_state.query = compound_name
    st.session_state.show_analyze_section = False  # Reset analysis section
    
    with st.spinner("Searching for retrosynthesis pathways..."):
        try:
            # Call your retrosynthesis function
            retro_result = run_retrosynthesis(compound_name)
            
            if retro_result.get("status") == "success":
                # Process the reactions to ensure we have usable SMILES
                for reaction in retro_result["data"]["reactions"]:
                    # Extract or create valid reaction SMILES
                    reaction["cleaned_reaction_smiles"] = extract_reaction_smiles(reaction)
                
                st.session_state.retro_results = retro_result["data"]
                st.success(f"Found {len(retro_result['data']['reactions'])} reactions for {compound_name}")
            else:
                st.error(f"Error: {retro_result.get('message', 'Unknown error')}")
                if retro_result.get('traceback'):
                    with st.expander("See error details"):
                        st.code(retro_result.get('traceback'))
        except Exception as e:
            st.error(f"Failed to fetch retrosynthesis data: {str(e)}")
            import traceback
            with st.expander("See error details"):
                st.code(traceback.format_exc())

# Display retrosynthesis results if available
if st.session_state.retro_results and not st.session_state.show_analyze_section:
    st.markdown("## Retrosynthesis Pathway")
    
    # Display recommended pathway
    st.markdown("### Recommended Synthesis Route")
    st.markdown(st.session_state.retro_results["reasoning"])
    
    # Display reactions
    st.markdown("### Reactions")
    
    for i, reaction in enumerate(st.session_state.retro_results["reactions"]):
        with st.container():
            st.markdown(f'<div class="reaction-card">', unsafe_allow_html=True)
            
            # Highlight if this is a recommended step
            if reaction["idx"] in st.session_state.retro_results["recommended_indices"]:
                st.markdown(f"**Step {i+1} (Recommended)**: {reaction['idx']}")
            else:
                st.markdown(f"**Step {i+1}**: {reaction['idx']}")
            
            # Reactants and products
            reactants_str = " + ".join(reaction["reactants"])
            products_str = " + ".join(reaction["products"])
            st.markdown(f"**Reaction**: {reactants_str} → {products_str}")
            
            # Conditions
            st.markdown(f"**Conditions**: {reaction['conditions']}")
            
            # Source
            st.markdown(f"**Source**: {reaction['source']}")
            
            # Show the cleaned reaction SMILES for debugging (can be removed in production)
            with st.expander("Reaction SMILES"):
                st.code(reaction.get("cleaned_reaction_smiles", "No valid SMILES available"))
            
            # Create a button that will trigger the analysis
            if st.button(f"Analyze Reaction {i+1}", key=f"analyze_btn_{i}"):
                # Use the cleaned SMILES if available, otherwise generate from reactants and products
                if "cleaned_reaction_smiles" in reaction and reaction["cleaned_reaction_smiles"]:
                    st.session_state.selected_rxn_smiles = reaction["cleaned_reaction_smiles"]
                    st.session_state.query = f"Give full information about this rxn {reaction['cleaned_reaction_smiles']}"
                    st.session_state.show_analyze_section = True
                    st.rerun()
                else:
                    if analyze_reaction(reaction["reactants"], reaction["products"]):
                        st.rerun()
            
            st.markdown('</div>', unsafe_allow_html=True)

# Analysis Section - shows either when a reaction is selected for analysis or when directly analyzing
if st.session_state.show_analyze_section:
    st.markdown('<div class="analysis-section">', unsafe_allow_html=True)
    st.markdown("## Reaction Analysis")
    
    # Show selected reaction SMILES
    if st.session_state.selected_rxn_smiles:
        st.markdown("### Selected Reaction")
        st.code(st.session_state.selected_rxn_smiles)
    
    # Process the selected reaction for analysis
    if st.session_state.selected_rxn_smiles and not st.session_state.analysis_result:
        with st.spinner("Analyzing reaction..."):
            # Create the full analysis query
            analysis_query = f"Give full information about this rxn {st.session_state.selected_rxn_smiles}"
            
            # Use the callback handler for streaming updates
            callback_container = st.container()
            st_callback = StreamlitCallbackHandler(callback_container)
            
            try:
                # Pass the callback to your backend function
                result = enhanced_query(analysis_query, callbacks=[st_callback])
                st.session_state.analysis_result = result
            except Exception as e:
                st.error(f"An error occurred during analysis: {str(e)}")
                import traceback
                with st.expander("See error details"):
                    st.code(traceback.format_exc())
    
    # Display analysis results
    if st.session_state.analysis_result:
        # Display visualization if available
        if st.session_state.analysis_result.get('visualization_path') and os.path.exists(st.session_state.analysis_result.get('visualization_path')):
            st.markdown("### Reaction Visualization")
            st.image(st.session_state.analysis_result.get('visualization_path'), caption="Chemical Reaction Visualization")
            
            # Add this line to make file path less prominent but still available
            with st.expander("Image file details"):
                st.code(st.session_state.analysis_result.get('visualization_path'))
        
        # Display the analysis text
        st.markdown("### Analysis Results")
        st.markdown('<div class="result-area">', unsafe_allow_html=True)
        st.markdown(st.session_state.analysis_result.get('analysis', 'No analysis available'))
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Add a query box for follow-up questions about the reaction
        st.markdown('<div class="query-box">', unsafe_allow_html=True)
        st.markdown("### Ask About This Reaction")
        
        rxn_query = st.text_area(
            "Enter your question about this reaction:",
            value=st.session_state.rxn_query,
            height=80,
            placeholder="Example: What is the mechanism of this reaction? How does this reaction relate to industrial processes?",
            key="rxn_query_input"
        )
        
        if st.button("Ask Question", key="ask_rxn_btn") and rxn_query:
            st.session_state.rxn_query = rxn_query
            
            with st.spinner("Processing your question..."):
                try:
                    # Create a query that includes both the reaction SMILES and the user's question
                    combined_query = f"For the reaction {st.session_state.selected_rxn_smiles}, answer this question: {rxn_query}"
                    
                    # Use the callback handler for streaming updates
                    query_response_container = st.container()
                    query_st_callback = StreamlitCallbackHandler(query_response_container)
                    
                    # Get the response
                    query_result = enhanced_query(combined_query, callbacks=[query_st_callback])
                    
                    # Display the answer
                    st.markdown("### Answer")
                    st.markdown('<div class="result-area">', unsafe_allow_html=True)
                    st.markdown(query_result.get('analysis', 'No answer available'))
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                except Exception as e:
                    st.error(f"An error occurred: {str(e)}")
                    st.markdown("Please check your question format and try again.")
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Add a button to clear analysis and return to retrosynthesis results
        if st.button("← Back to Retrosynthesis Results", key="back_btn"):
            st.session_state.show_analyze_section = False
            st.session_state.selected_rxn_smiles = None
            st.session_state.analysis_result = None
            st.session_state.rxn_query = ""
            st.rerun()
            
    st.markdown('</div>', unsafe_allow_html=True)

# Footer
st.markdown("---")
col1, col2, col3 = st.columns([1, 2, 1])
with col2:
    st.caption("ChemCopilot - Your Expert Chemistry Assistant")












# app.py (Simplified for the restored test.py)

# import streamlit as st
# import requests
# import re
# import sys
# import os
# import json
# from langchain_community.callbacks.streamlit import StreamlitCallbackHandler
# from langchain_openai import ChatOpenAI


# # Add the project root to the Python path for imports
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))
# import api_config

# os.environ['STREAMLIT_BROWSER_GATHER_USAGE_STATS'] = 'false'

# # Import your function from the test.py file
# from test import enhanced_query
# # Import retrosynthesis function
# from tools.retrosynthesis import run_retrosynthesis
# # Import NameToSMILES tool
# from tools.name2smiles import NameToSMILES  # Update this with the correct import path

# # Initialize LLM for name to SMILES conversion
# llm = ChatOpenAI(model="gpt-4o", temperature=0)
# # Initialize NameToSMILES tool
# name_to_smiles_tool = NameToSMILES()

# # Set up the Streamlit page
# st.set_page_config(
#     page_title="ChemCopilot - Chemistry Assistant",
#     page_icon="🧪",
#     layout="wide",
#     initial_sidebar_state="expanded"
# )

# # Custom CSS for styling
# st.markdown("""
#     <style>
#     .main-header {
#         font-size: 42px;
#         font-weight: bold;
#         color: #2e7d32;
#         margin-bottom: 0px;
#     }
#     .sub-header {
#         font-size: 20px;
#         color: #5c5c5c;
#         margin-bottom: 30px;
#     }
#     .stButton>button {
#         background-color: #2e7d32;
#         color: white;
#         border: none;
#         padding: 10px 24px;
#         border-radius: 4px;
#         font-weight: bold;
#     }
#     .stButton>button:hover {
#         background-color: #005005;
#     }
#     .tool-card {
#         background-color: #f5f5f5;
#         padding: 10px;
#         border-radius: 5px;
#         margin-bottom: 10px;
#     }
#     .result-area {
#         background-color: #f9f9f9;
#         padding: 20px;
#         border-radius: 8px;
#         border-left: 4px solid #2e7d32;
#     }
#     .reaction-card {
#         background-color: #f0f7f0;
#         padding: 15px;
#         border-radius: 8px;
#         margin-bottom: 15px;
#         border-left: 4px solid #2e7d32;
#     }
#     .analysis-section {
#         margin-top: 30px;
#         padding-top: 20px;
#         border-top: 1px solid #ddd;
#     }
#     .query-box {
#         background-color: #f5f7fa;
#         padding: 15px;
#         border-radius: 8px;
#         margin-top: 20px;
#         border-left: 4px solid #3f51b5;
#     }
#     .info-item {
#         margin-bottom: 15px;
#         padding-bottom: 10px;
#         border-bottom: 1px solid #eaeaea;
#     }
#     .info-header {
#         font-weight: bold;
#         color: #2e7d32;
#     }
#     .info-section {
#         margin-top: 20px;
#         padding: 15px;
#         background-color: #f7f9fc;
#         border-radius: 8px;
#         border-left: 4px solid #2e7d32;
#     }
#     </style>
#     """, unsafe_allow_html=True)

# # Initialize session state variables
# if 'query' not in st.session_state:
#     st.session_state.query = ""
# if 'retro_results' not in st.session_state:
#     st.session_state.retro_results = None
# if 'selected_rxn_smiles' not in st.session_state:
#     st.session_state.selected_rxn_smiles = None
# if 'analysis_result' not in st.session_state:
#     st.session_state.analysis_result = None
# if 'rxn_query' not in st.session_state:
#     st.session_state.rxn_query = ""
# if 'show_analyze_section' not in st.session_state:
#     st.session_state.show_analyze_section = False
# if 'raw_data' not in st.session_state:
#     st.session_state.raw_data = None

# # Function to get SMILES from chemical name using NameToSMILES tool
# def get_smiles_from_name(name):
#     """Get SMILES notation for a chemical name using the NameToSMILES tool"""
#     try:
#         result = name_to_smiles_tool._run(name)
#         if "SMILES:" in result:
#             # Extract SMILES from the result
#             smiles = result.split("SMILES:")[1].split("\n")[0].strip()
#             return smiles
#         return None
#     except Exception as e:
#         st.error(f"Error converting {name} to SMILES: {str(e)}")
#         return None

# # Function to convert reactant and product names to SMILES
# def convert_to_reaction_smiles(reactants, products):
#     """Convert reactant and product names to reaction SMILES format using LLM and NameToSMILES tool"""
#     reactant_smiles = []
#     product_smiles = []
    
#     # First try to convert each reactant and product using NameToSMILES tool
#     for reactant in reactants:
#         smiles = get_smiles_from_name(reactant)
#         if smiles:
#             reactant_smiles.append(smiles)
#         else:
#             reactant_smiles.append(reactant)  # Keep original name if conversion fails
    
#     for product in products:
#         smiles = get_smiles_from_name(product)
#         if smiles:
#             product_smiles.append(smiles)
#         else:
#             product_smiles.append(product)  # Keep original name if conversion fails
    
#     # If NameToSMILES failed for some compounds, use LLM as backup
#     if any(not re.match(r'^[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+$', r) for r in reactant_smiles) or \
#        any(not re.match(r'^[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+$', p) for p in product_smiles):
        
#         prompt = f"""Convert these chemical names to a reaction SMILES format.

# Reactants: {', '.join(reactants)}
# Products: {', '.join(products)}

# Format the output as reaction SMILES using the format: reactant1.reactant2>>product1.product2
# Only output the SMILES, no explanations.
# """
        
#         try:
#             response = llm.invoke(prompt)
#             reaction_smiles = response.content.strip()
            
#             # Verify that the output looks like a reaction SMILES
#             if ">>" in reaction_smiles:
#                 return reaction_smiles
#         except Exception as e:
#             st.error(f"Error converting to reaction SMILES with LLM: {str(e)}")
    
#     # Join SMILES directly if all conversions were successful
#     reactants_str = '.'.join(reactant_smiles)
#     products_str = '.'.join(product_smiles)
#     return f"{reactants_str}>>{products_str}"

# # Function to extract clean SMILES from reaction data
# def extract_reaction_smiles(reaction):
#     """Extract or generate clean reaction SMILES from reaction data"""
#     # First check if there's a valid reaction_smiles that doesn't have placeholders
#     if "reaction_smiles" in reaction and reaction["reaction_smiles"]:
#         smiles = reaction["reaction_smiles"]
#         # Check if it contains placeholders
#         if not re.search(r'\[.*?_SMILES\]|\[.*?\]', smiles) and ">>" in smiles:
#             return smiles
    
#     # If reaction_smiles is not usable, create from reactants and products
#     return convert_to_reaction_smiles(reaction["reactants"], reaction["products"])

# # Function to create reaction SMILES from reactants and products
# def create_reaction_smiles(reactants, products):
#     """Create reaction SMILES from reactants and products, handling both SMILES and name inputs"""
#     # Check if inputs are likely SMILES or names
#     smiles_pattern = r'^[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+$'
    
#     # Check if all reactants and products are SMILES
#     all_smiles = all(re.match(smiles_pattern, r) for r in reactants) and all(re.match(smiles_pattern, p) for p in products)
    
#     if all_smiles:
#         # Join SMILES directly
#         reactants_str = '.'.join(reactants)
#         products_str = '.'.join(products)
#         return f"{reactants_str}>>{products_str}"
#     else:
#         # Use conversion function
#         return convert_to_reaction_smiles(reactants, products)

# # Function to handle reaction analysis
# def analyze_reaction(reactants, products):
#     """Analyze reaction by converting names to SMILES if necessary"""
#     reaction_smiles = create_reaction_smiles(reactants, products)
    
#     if reaction_smiles:
#         st.session_state.selected_rxn_smiles = reaction_smiles
#         st.session_state.query = f"Give full information about this rxn {reaction_smiles}"
#         st.session_state.show_analyze_section = True
#         return True
#     else:
#         st.error("Could not create reaction SMILES. Please check the reactants and products.")
#         return False

# # Function to format raw data into structured markdown
# def format_raw_data_markdown(raw_data):
#     """Convert raw data to structured markdown"""
#     md_content = ""
    
#     if not raw_data:
#         return "No raw data available"
    
#     # Try to parse JSON if it's a string
#     if isinstance(raw_data, str):
#         try:
#             data = json.loads(raw_data)
#         except:
#             # If not valid JSON, just return as is
#             return f"```\n{raw_data}\n```"
#     else:
#         data = raw_data
    
#     # Check if the data is a dictionary with tool outputs
#     if isinstance(data, dict):
#         # Handle different types of data structures
        
#         # First check for reaction classification
#         if "ReactionClassifier" in data:
#             md_content += "### Reaction Classification\n\n"
#             classifier_data = data["ReactionClassifier"]
#             if isinstance(classifier_data, dict):
#                 for key, value in classifier_data.items():
#                     md_content += f"**{key}**: {value}\n\n"
#             else:
#                 md_content += f"{classifier_data}\n\n"
        
#         # Check for bond changes
#         if "BondChangeAnalyzer" in data:
#             md_content += "### Bond Changes\n\n"
#             bond_data = data["BondChangeAnalyzer"]
#             if isinstance(bond_data, dict):
#                 for key, value in bond_data.items():
#                     md_content += f"**{key}**: {value}\n\n"
#             elif isinstance(bond_data, list):
#                 for item in bond_data:
#                     md_content += f"- {item}\n"
#             else:
#                 md_content += f"{bond_data}\n\n"
        
#         # Check for functional groups
#         if "FuncGroups" in data:
#             md_content += "### Functional Groups\n\n"
#             func_data = data["FuncGroups"]
#             if isinstance(func_data, dict):
#                 # Handle reactants and products separately
#                 if "reactants" in func_data:
#                     md_content += "#### Reactants\n\n"
#                     for compound, groups in func_data["reactants"].items():
#                         md_content += f"**{compound}**:\n"
#                         if isinstance(groups, list):
#                             for group in groups:
#                                 md_content += f"- {group}\n"
#                         else:
#                             md_content += f"{groups}\n"
#                         md_content += "\n"
                
#                 if "products" in func_data:
#                     md_content += "#### Products\n\n"
#                     for compound, groups in func_data["products"].items():
#                         md_content += f"**{compound}**:\n"
#                         if isinstance(groups, list):
#                             for group in groups:
#                                 md_content += f"- {group}\n"
#                         else:
#                             md_content += f"{groups}\n"
#                         md_content += "\n"
#             elif isinstance(func_data, list):
#                 for item in func_data:
#                     md_content += f"- {item}\n"
#             else:
#                 md_content += f"{func_data}\n\n"
        
#         # Check for reaction SMILES
#         if "reaction_smiles" in data:
#             md_content += "### Reaction SMILES\n\n"
#             md_content += f"```\n{data['reaction_smiles']}\n```\n\n"
        
#         # Check for reactants and products
#         if "reactants" in data:
#             md_content += "### Reactants\n\n"
#             for reactant in data["reactants"]:
#                 md_content += f"- {reactant}\n"
#             md_content += "\n"
        
#         if "products" in data:
#             md_content += "### Products\n\n"
#             for product in data["products"]:
#                 md_content += f"- {product}\n"
#             md_content += "\n"
        
#         # Add any other keys that weren't specifically handled
#         for key, value in data.items():
#             if key not in ["ReactionClassifier", "BondChangeAnalyzer", "FuncGroups", 
#                           "reaction_smiles", "reactants", "products", "analysis"]:
#                 md_content += f"### {key}\n\n"
#                 if isinstance(value, dict):
#                     for sub_key, sub_value in value.items():
#                         md_content += f"**{sub_key}**: {sub_value}\n\n"
#                 elif isinstance(value, list):
#                     for item in value:
#                         md_content += f"- {item}\n"
#                     md_content += "\n"
#                 else:
#                     md_content += f"{value}\n\n"
    
#     # If data is a list
#     elif isinstance(data, list):
#         for item in data:
#             if isinstance(item, dict):
#                 for key, value in item.items():
#                     md_content += f"**{key}**: {value}\n\n"
#             else:
#                 md_content += f"- {item}\n"
#         md_content += "\n"
    
#     # If just a string or other type
#     else:
#         md_content = f"```\n{str(data)}\n```"
    
#     return md_content

# # App Header
# st.markdown('<p class="main-header">ChemCopilot</p>', unsafe_allow_html=True)
# st.markdown('<p class="sub-header">Your Expert Chemistry Assistant</p>', unsafe_allow_html=True)

# # Sidebar with tools info and examples
# with st.sidebar:
#     st.markdown("## Tools Available")
    
#     with st.expander("🔍 SMILES2Name"):
#         st.markdown("Converts SMILES notation to chemical names.")
#         st.markdown("*Example:* `C1=CC=CC=C1` → `Benzene`")
    
#     with st.expander("📝 Name2SMILES"):
#         st.markdown("Converts chemical names to SMILES notation.")
#         st.markdown("*Example:* `Ethanol` → `CCO`")
    
#     with st.expander("🧪 FuncGroups"):
#         st.markdown("Analyzes functional groups in molecules.")
#         st.markdown("*Example:* `C(O)(=O)C1=CC=CC=C1` → `Carboxylic acid, Aromatic ring`")
    
#     with st.expander("⚛️ BondChangeAnalyzer"):
#         st.markdown("Analyzes bond changes in chemical reactions.")
#         st.markdown("*Example:* `CCCl.CC[O-].[Na+]>>CCOCC.[Na+].[Cl-]` → `C-Cl bond broken, C-O bond formed`")
    
#     st.markdown("## Example Queries")
    
#     if st.button("Search for Flubendiamide"):
#         st.session_state.query = "flubendiamide"
#         st.session_state.show_analyze_section = False
#         st.rerun()
        
#     if st.button("Search for 2-amino-5-chloro-3-methyl benzoic acid"):
#         st.session_state.query = "2-amino-5-chloro-3-methyl benzoic acid"
#         st.session_state.show_analyze_section = False
#         st.rerun()
        
#     if st.button("Reaction Analysis Example"):
#         st.session_state.query = "Give full information about this rxn CCCl.CC[O-].[Na+]>>CCOCC.[Na+].[Cl-]"
#         st.session_state.selected_rxn_smiles = "CCCl.CC[O-].[Na+]>>CCOCC.[Na+].[Cl-]"
#         st.session_state.show_analyze_section = True
#         st.rerun()

# # Main content area - integrated into a single page
# st.markdown("## Retrosynthesis & Reaction Analysis")

# # Search Input Section
# compound_name = st.text_input(
#     "Enter compound name (IUPAC or common name):",
#     value=st.session_state.query if not st.session_state.show_analyze_section else "",
#     placeholder="Example: flubendiamide, 2-amino-5-chloro-3-methyl benzoic acid",
#     key="search_compound_name"
# )

# search_button = st.button("Search Retrosynthesis", key="search_retro_button")

# # Process Search
# if search_button and compound_name:
#     st.session_state.query = compound_name
#     st.session_state.show_analyze_section = False  # Reset analysis section
    
#     with st.spinner("Searching for retrosynthesis pathways..."):
#         try:
#             # Call your retrosynthesis function
#             retro_result = run_retrosynthesis(compound_name)
            
#             if retro_result.get("status") == "success":
#                 # Process the reactions to ensure we have usable SMILES
#                 for reaction in retro_result["data"]["reactions"]:
#                     # Extract or create valid reaction SMILES
#                     reaction["cleaned_reaction_smiles"] = extract_reaction_smiles(reaction)
                
#                 st.session_state.retro_results = retro_result["data"]
#                 st.success(f"Found {len(retro_result['data']['reactions'])} reactions for {compound_name}")
#             else:
#                 st.error(f"Error: {retro_result.get('message', 'Unknown error')}")
#                 if retro_result.get('traceback'):
#                     with st.expander("See error details"):
#                         st.code(retro_result.get('traceback'))
#         except Exception as e:
#             st.error(f"Failed to fetch retrosynthesis data: {str(e)}")
#             import traceback
#             with st.expander("See error details"):
#                 st.code(traceback.format_exc())

# # Display retrosynthesis results if available
# if st.session_state.retro_results and not st.session_state.show_analyze_section:
#     st.markdown("## Retrosynthesis Pathway")
    
#     # Display recommended pathway
#     st.markdown("### Recommended Synthesis Route")
#     st.markdown(st.session_state.retro_results["reasoning"])
    
#     # Display reactions
#     st.markdown("### Reactions")
    
#     for i, reaction in enumerate(st.session_state.retro_results["reactions"]):
#         with st.container():
#             st.markdown(f'<div class="reaction-card">', unsafe_allow_html=True)
            
#             # Highlight if this is a recommended step
#             if reaction["idx"] in st.session_state.retro_results["recommended_indices"]:
#                 st.markdown(f"**Step {i+1} (Recommended)**: {reaction['idx']}")
#             else:
#                 st.markdown(f"**Step {i+1}**: {reaction['idx']}")
            
#             # Reactants and products
#             reactants_str = " + ".join(reaction["reactants"])
#             products_str = " + ".join(reaction["products"])
#             st.markdown(f"**Reaction**: {reactants_str} → {products_str}")
            
#             # Conditions
#             st.markdown(f"**Conditions**: {reaction['conditions']}")
            
#             # Source
#             st.markdown(f"**Source**: {reaction['source']}")
            
#             # Show the cleaned reaction SMILES for debugging (can be removed in production)
#             with st.expander("Reaction SMILES"):
#                 st.code(reaction.get("cleaned_reaction_smiles", "No valid SMILES available"))
            
#             # Create a button that will trigger the analysis
#             if st.button(f"Analyze Reaction {i+1}", key=f"analyze_btn_{i}"):
#                 # Use the cleaned SMILES if available, otherwise generate from reactants and products
#                 if "cleaned_reaction_smiles" in reaction and reaction["cleaned_reaction_smiles"]:
#                     st.session_state.selected_rxn_smiles = reaction["cleaned_reaction_smiles"]
#                     st.session_state.query = f"Give full information about this rxn {reaction['cleaned_reaction_smiles']}"
#                     st.session_state.show_analyze_section = True
#                     st.rerun()
#                 else:
#                     if analyze_reaction(reaction["reactants"], reaction["products"]):
#                         st.rerun()
            
#             st.markdown('</div>', unsafe_allow_html=True)

# # Analysis Section - shows either when a reaction is selected for analysis or when directly analyzing
# if st.session_state.show_analyze_section:
#     st.markdown('<div class="analysis-section">', unsafe_allow_html=True)
#     st.markdown("## Reaction Analysis")
    
#     # Show selected reaction SMILES
#     if st.session_state.selected_rxn_smiles:
#         st.markdown("### Selected Reaction")
#         st.code(st.session_state.selected_rxn_smiles)
    
#     # Process the selected reaction for analysis
#     if st.session_state.selected_rxn_smiles and not st.session_state.analysis_result:
#         with st.spinner("Analyzing reaction..."):
#             # Create the full analysis query
#             analysis_query = f"Give full information about this rxn {st.session_state.selected_rxn_smiles}"
            
#             # Use the callback handler for streaming updates
#             callback_container = st.container()
#             st_callback = StreamlitCallbackHandler(callback_container)
            
#             try:
#                 # Pass the callback to your backend function
#                 result = enhanced_query(analysis_query, callbacks=[st_callback])
#                 st.session_state.analysis_result = result
#                 # Store raw data for detailed markdown rendering
#                 if 'raw_data' in result:
#                     st.session_state.raw_data = result['raw_data']
#             except Exception as e:
#                 st.error(f"An error occurred during analysis: {str(e)}")
#                 import traceback
#                 with st.expander("See error details"):
#                     st.code(traceback.format_exc())
    
#     # Display analysis results
#     if st.session_state.analysis_result:
#         # Display visualization if available
#         if st.session_state.analysis_result.get('visualization_path') and os.path.exists(st.session_state.analysis_result.get('visualization_path')):
#             st.markdown("### Reaction Visualization")
#             st.image(st.session_state.analysis_result.get('visualization_path'), caption="Chemical Reaction Visualization")
            
#             # Add this line to make file path less prominent but still available
#             with st.expander("Image file details"):
#                 st.code(st.session_state.analysis_result.get('visualization_path'))
        
#         # Display the analysis text
#         st.markdown("### Analysis Results")
#         st.markdown('<div class="result-area">', unsafe_allow_html=True)
#         st.markdown(st.session_state.analysis_result.get('analysis', 'No analysis available'))
#         st.markdown('</div>', unsafe_allow_html=True)
        
#         # Display detailed information from tools in a structured format
#         if st.session_state.raw_data:
#             with st.expander("Detailed Chemistry Data"):
#                 # Format and display the raw data in a structured markdown format
#                 formatted_md = format_raw_data_markdown(st.session_state.raw_data)
#                 st.markdown(formatted_md)
                
#                 # Optionally show the raw JSON
#                 with st.expander("Raw JSON Data"):
#                     if isinstance(st.session_state.raw_data, str):
#                         st.code(st.session_state.raw_data, language="json")
#                     else:
#                         st.code(json.dumps(st.session_state.raw_data, indent=2), language="json")
        
#         # Add reaction classification and other tool outputs in a structured format
#         col1, col2 = st.columns(2)
        
#         # Create collapsed sections for different information categories
#         with col1:
#             if st.session_state.analysis_result.get('reaction_classification'):
#                 with st.expander("Reaction Classification"):
#                     st.markdown(st.session_state.analysis_result.get('reaction_classification', 'Not available'))
            
#             if st.session_state.analysis_result.get('bond_changes'):
#                 with st.expander("Bond Changes"):
#                     st.markdown(st.session_state.analysis_result.get('bond_changes', 'Not available'))
        
#         with col2:
#             if st.session_state.analysis_result.get('functional_groups'):
#                 with st.expander("Functional Groups"):
#                     st.markdown(st.session_state.analysis_result.get('functional_groups', 'Not available'))
                    
#             if st.session_state.analysis_result.get('mechanism'):
#                 with st.expander("Reaction Mechanism"):
#                     st.markdown(st.session_state.analysis_result.get('mechanism', 'Not available'))
        
#         # Add a query box for follow-up questions about the reaction
#         st.markdown('<div class="query-box">', unsafe_allow_html=True)
#         st.markdown("### Ask About This Reaction")
        
#         rxn_query = st.text_area(
#             "Enter your question about this reaction:",
#             value=st.session_state.rxn_query,
#             height=80,
#             placeholder="Example: What is the mechanism of this reaction? How does this reaction relate to industrial processes?",
#             key="rxn_query_input"
#         )
        
#         if st.button("Ask Question", key="ask_rxn_btn") and rxn_query:
#             st.session_state.rxn_query = rxn_query
            
#             with st.spinner("Processing your question..."):
#                 try:
#                     # Create a query that includes both the reaction SMILES and the user's question
#                     combined_query = f"For the reaction {st.session_state.selected_rxn_smiles}, answer this question: {rxn_query}"
                    
#                     # Use the callback handler for streaming updates
#                     query_response_container = st.container()
#                     query_st_callback = StreamlitCallbackHandler(query_response_container)
                    
#                     # Get the response
#                     query_result = enhanced_query(combined_query, callbacks=[query_st_callback])
                    
#                     # Display the answer
#                     st.markdown("### Answer")
#                     st.markdown('<div class="result-area">', unsafe_allow_html=True)
#                     st.markdown(query_result.get('analysis', 'No answer available'))
#                     st.markdown('</div>', unsafe_allow_html=True)
                    
#                     # If there's raw data in the query result, show it too
#                     if 'raw_data' in query_result:
#                         with st.expander("Detailed Chemistry Data"):
#                             formatted_md = format_raw_data_markdown(query_result['raw_data'])
#                             st.markdown(formatted_md)
                    
#                 except Exception as e:
#                     st.error(f"An error occurred: {str(e)}")
#                     st.markdown("Please check your question format and try again.")
        
#         st.markdown('</div>', unsafe_allow_html=True)
        
#         # Add a button to clear analysis and return to retrosynthesis results
#         if st.button("← Back to Retrosynthesis Results", key="back_btn"):
#             st.session_state.show_analyze_section = False
#             st.session_state.selected_rxn_smiles = None
#             st.session_state.analysis_result = None
#             st.session_state.raw_data = None
#             st.session_state.rxn_query = ""
#             st.rerun()
            
#     st.markdown('</div>', unsafe_allow_html=True)

# # Footer
# st.markdown("---")
# col1, col2, col3 = st.columns([1, 2, 1])
# with col2:
#     st.caption("ChemCopilot - Your Expert Chemistry Assistant")