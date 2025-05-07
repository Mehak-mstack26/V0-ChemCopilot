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

# Initialize LLM for name to SMILES conversion and analysis
llm = ChatOpenAI(model="gpt-4o", temperature=0, streaming=True) # Added streaming=True
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
    .pathway-card {
        background-color: #e8f5e9;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 20px;
        border-left: 4px solid #1b5e20;
    }
    </style>
    """, unsafe_allow_html=True)

# Initialize session state variables
if 'query' not in st.session_state:
    st.session_state.query = ""
if 'retro_results' not in st.session_state:
    st.session_state.retro_results = None # Will store the 'data' part of successful retro_result
if 'all_parsed_reactions_by_idx' not in st.session_state:
    st.session_state.all_parsed_reactions_by_idx = {}
if 'ranked_pathways' not in st.session_state:
    st.session_state.ranked_pathways = []
if 'selected_rxn_smiles' not in st.session_state:
    st.session_state.selected_rxn_smiles = None
if 'analysis_result_text' not in st.session_state: # For storing text from LLM
    st.session_state.analysis_result_text = ""
if 'rxn_query' not in st.session_state:
    st.session_state.rxn_query = ""
if 'show_analyze_section' not in st.session_state:
    st.session_state.show_analyze_section = False

# Function to get SMILES from chemical name using NameToSMILES tool
def get_smiles_from_name(name):
    """Get SMILES notation for a chemical name using the NameToSMILES tool"""
    try:
        result = name_to_smiles_tool._run(name)
        if result and "SMILES:" in result: # Added check for result not being None
            smiles = result.split("SMILES:")[1].split("\n")[0].strip()
            return smiles
        st.warning(f"Could not convert '{name}' to SMILES using NameToSMILES tool. Result: {result}")
        return None
    except Exception as e:
        st.error(f"Error converting {name} to SMILES: {str(e)}")
        return None

# Function to convert reactant and product names to SMILES for a reaction
def convert_names_to_reaction_smiles(reactants_names, products_names):
    """Convert lists of reactant and product names to a reaction SMILES string."""
    reactant_smiles_list = [get_smiles_from_name(name) or name for name in reactants_names]
    product_smiles_list = [get_smiles_from_name(name) or name for name in products_names]

    # Check if all components were successfully converted to actual SMILES (heuristic: no spaces)
    all_react_valid_smi = all(smi and not ' ' in smi for smi in reactant_smiles_list)
    all_prod_valid_smi = all(smi and not ' ' in smi for smi in product_smiles_list)

    if reactant_smiles_list and product_smiles_list and all_react_valid_smi and all_prod_valid_smi:
        return f"{'.'.join(reactant_smiles_list)}>>{'.'.join(product_smiles_list)}"
    else:
        st.warning(f"Failed to convert all names to SMILES for reaction. R: {reactant_smiles_list}, P: {product_smiles_list}")
        # Fallback to name-based representation if conversion fails for some
        return f"{'.'.join(reactants_names)}>>{'.'.join(products_names)}"


# Function to extract or generate clean reaction SMILES for analysis
def get_reaction_smiles_for_analysis(reaction_obj):
    """
    Prioritizes cleaned_reaction_smiles from backend. 
    If not available or invalid, generates from reactants/products names.
    """
    if reaction_obj.get("cleaned_reaction_smiles"):
        # Basic validation for ">>"
        if ">>" in reaction_obj["cleaned_reaction_smiles"]:
            return reaction_obj["cleaned_reaction_smiles"]
        else:
            st.warning(f"Backend 'cleaned_reaction_smiles' for idx {reaction_obj.get('idx')} ('{reaction_obj.get('cleaned_reaction_smiles')}') seems invalid. Attempting regeneration.")
    
    # Fallback to generating from names if cleaned_reaction_smiles is missing, None, or invalid
    st.info(f"Generating SMILES for reaction idx {reaction_obj.get('idx')} from reactant/product names as 'cleaned_reaction_smiles' was not suitable.")
    return convert_names_to_reaction_smiles(reaction_obj.get("reactants", []), reaction_obj.get("products", []))


# App Header
st.markdown('<p class="main-header">ChemCopilot</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Your Expert Chemistry Assistant</p>', unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.markdown("## Tools Available")
    with st.expander("🔍 SMILES2Name"): st.markdown("Converts SMILES notation to chemical names.\n*Example:* C1=CC=CC=C1 → Benzene")
    with st.expander("📝 Name2SMILES"): st.markdown("Converts chemical names to SMILES notation.\n*Example:* Ethanol → CCO")
    with st.expander("🧪 FuncGroups"): st.markdown("Analyzes functional groups in molecules.\n*Example:* C(O)(=O)C1=CC=CC=C1 → Carboxylic acid, Aromatic ring")
    with st.expander("⚛️ BondChangeAnalyzer"): st.markdown("Analyzes bond changes in chemical reactions.\n*Example:* CCCl.CC[O-].[Na+]>>CCOCC.[Na+].[Cl-] → C-Cl bond broken, C-O bond formed")
    st.markdown("## Example Queries")
    if st.button("Search for Flubendiamide"):
        st.session_state.query = "flubendiamide"
        st.session_state.show_analyze_section = False
        st.session_state.retro_results = None
        st.session_state.all_parsed_reactions_by_idx = {}
        st.session_state.ranked_pathways = []
        st.rerun()
    if st.button("Search for 2-amino-5-chloro-3-methyl benzoic acid"):
        st.session_state.query = "2-amino-5-chloro-3-methyl benzoic acid"
        st.session_state.show_analyze_section = False
        st.session_state.retro_results = None
        st.session_state.all_parsed_reactions_by_idx = {}
        st.session_state.ranked_pathways = []
        st.rerun()
    if st.button("Reaction Analysis Example"):
        st.session_state.query = "" # Clear compound search
        st.session_state.selected_rxn_smiles = "CCCl.CC[O-].[Na+]>>CCOCC.[Na+].[Cl-]"
        st.session_state.show_analyze_section = True
        st.session_state.retro_results = None
        st.session_state.all_parsed_reactions_by_idx = {}
        st.session_state.ranked_pathways = []
        st.rerun()

# Main content area
st.markdown("## Retrosynthesis & Reaction Analysis")

compound_name_input = st.text_input(
    "Enter compound name (IUPAC or common name):",
    value=st.session_state.query if not st.session_state.show_analyze_section else "",
    placeholder="Example: flubendiamide",
    key="search_compound_name_input"
)

search_button = st.button("Search Retrosynthesis", key="search_retro_button")

if search_button and compound_name_input:
    st.session_state.query = compound_name_input
    st.session_state.show_analyze_section = False 
    st.session_state.retro_results = None # Clear previous results
    st.session_state.all_parsed_reactions_by_idx = {}
    st.session_state.ranked_pathways = []

    with st.spinner(f"Searching for retrosynthesis pathways for '{compound_name_input}'... This may take a few minutes."):
        try:
            retro_result_full = run_retrosynthesis(compound_name_input) # Call backend
            
            if retro_result_full and retro_result_full.get("status") == "success":
                st.session_state.retro_results = retro_result_full.get("data")
                st.session_state.all_parsed_reactions_by_idx = st.session_state.retro_results.get("all_parsed_reactions_by_idx", {})
                st.session_state.ranked_pathways = st.session_state.retro_results.get("ranked_pathways", [])
                
                # Check if any reactions were returned for the recommended pathway
                if st.session_state.retro_results and st.session_state.retro_results.get("reactions"):
                    st.success(f"Found {len(st.session_state.retro_results.get('reactions',[]))} reaction steps in the primary recommended pathway for {compound_name_input}.")
                else:
                    st.warning(f"Retrosynthesis search completed, but no specific recommended reaction steps were found for {compound_name_input}. The reasoning might provide general guidance.")
                    if not st.session_state.retro_results.get("reasoning"): # If no reasoning either
                         st.warning(f"No synthesis pathways or reasoning found for '{compound_name_input}' in the searched literature or databases.")
            else: # Error or no success
                error_message = retro_result_full.get('message', 'Unknown error during retrosynthesis.')
                st.error(f"Error: {error_message}")
                if retro_result_full and retro_result_full.get('traceback'):
                    with st.expander("See error details"):
                        st.code(retro_result_full.get('traceback'))
        except Exception as e:
            st.error(f"Failed to fetch retrosynthesis data: {str(e)}")
            import traceback
            with st.expander("See error details"):
                st.code(traceback.format_exc())
    st.rerun() # Rerun to display results

# Display retrosynthesis results
if st.session_state.retro_results and not st.session_state.show_analyze_section:
    st.markdown("## Retrosynthesis Pathway")
    
    # Display ranked pathways first
    # The key from backend is now 'ranked_pathways' which contains a list of dicts
    # Each dict has 'rank', 'pathway_indices_str', 'details_text', 'reasons'
    if st.session_state.ranked_pathways: # This now refers to st.session_state.retro_results.get("ranked_pathways", [])
        st.markdown("### Ranked Synthesis Routes")
        
        for pathway_info in st.session_state.ranked_pathways: # Iterate through the list of pathway dicts
            # Use 'pathway_indices_str' for the expander title
            expander_title = f"Rank {pathway_info.get('rank')}: Pathway {pathway_info.get('pathway_indices_str', 'N/A')}"
            with st.expander(expander_title):
                st.markdown("#### Details")
                # Use 'details_text' for the details content
                st.markdown(pathway_info.get('details_text', 'No details provided')) 
                
                if pathway_info.get('reasons'):
                    st.markdown("#### Reasoning")
                    st.markdown(pathway_info.get('reasons'))
    
    # Display reasoning for the recommended pathway (overall reasoning)
    if st.session_state.retro_results.get("reasoning"):
        st.markdown("### Recommendation Reasoning (Overall)") # Clarified this is overall
        st.markdown(st.session_state.retro_results.get("reasoning"))
    
    # Display reactions of the primary recommended pathway (reactions for Rank 1)
    recommended_reactions_list = st.session_state.retro_results.get("reactions", []) # This list is for Rank 1
    if recommended_reactions_list:
        st.markdown("### Recommended Synthesis Route Steps (Rank 1 Pathway)") # Clarified this is for Rank 1
        for i, reaction_obj in enumerate(recommended_reactions_list):
            # ... (the rest of the reaction display loop remains the same) ...
            with st.container():
                st.markdown(f'<div class="reaction-card">', unsafe_allow_html=True)
                st.markdown(f"*Step {i+1}* (Reaction ID: {reaction_obj.get('idx', 'N/A')})")
                
                reactants_str = " + ".join(reaction_obj.get("reactants", ["N/A"]))
                products_str = " + ".join(reaction_obj.get("products", ["N/A"]))
                st.markdown(f"*Reaction*: {reactants_str} → {products_str}")
                
                conditions_dict = reaction_obj.get("conditions", {})
                conditions_str_list = [f"{k.capitalize()}: {v}" for k, v in conditions_dict.items()] if conditions_dict else ["Not specified"]
                st.markdown(f"*Conditions*: {', '.join(conditions_str_list)}")
                
                st.markdown(f"*Source*: {reaction_obj.get('source', 'N/A')}")
                
                with st.expander("View SMILES and Processing Info"):
                    st.write(f"Original LLM SMILES: {reaction_obj.get('reaction_smiles', 'N/A')}")
                    st.write(f"Processed SMILES for Analysis: {reaction_obj.get('cleaned_reaction_smiles', 'N/A')}")
                    st.write(f"SMILES Generation Attempted by Backend: {reaction_obj.get('smiles_generation_attempted', False)}")
                
                # Key for the button should be unique for each step of the Rank 1 pathway
                analyze_button_key = f"analyze_btn_rank1_step_{i}_{reaction_obj.get('idx', 'no_idx_at_step_'+str(i))}"
                
                if st.button(f"Analyze Reaction Step {i+1}", key=analyze_button_key):
                    reaction_smiles_for_analysis = get_reaction_smiles_for_analysis(reaction_obj)
                    if reaction_smiles_for_analysis and ">>" in reaction_smiles_for_analysis:
                        st.session_state.selected_rxn_smiles = reaction_smiles_for_analysis
                        st.session_state.show_analyze_section = True
                        st.session_state.analysis_result_text = "" 
                        st.rerun()
                    else:
                        st.error(f"Could not obtain or generate a valid reaction SMILES for analysis. Attempted: '{reaction_smiles_for_analysis}'")
                
                st.markdown('</div>', unsafe_allow_html=True)
    elif not st.session_state.retro_results.get("reasoning") and not st.session_state.ranked_pathways:
        st.info(f"No specific retrosynthesis pathway details or reasoning were found for '{st.session_state.query}'.")


# Analysis Section
if st.session_state.show_analyze_section:
    # Analysis section code remains the same
    st.markdown('<div class="analysis-section">', unsafe_allow_html=True)
    st.markdown("## Reaction Analysis")
    
    if st.session_state.selected_rxn_smiles:
        st.markdown("### Analyzing Reaction:")
        st.code(st.session_state.selected_rxn_smiles)
        
        # Rest of the analysis section code...
        if not st.session_state.analysis_result_text:
            with st.spinner("Analyzing reaction with LLM..."):
                analysis_query_for_llm = f"Give full information about this rxn {st.session_state.selected_rxn_smiles}"
                
                analysis_display_area = st.empty()
                st_callback_analysis = StreamlitCallbackHandler(analysis_display_area)
                
                try:
                    llm_response_data = enhanced_query(analysis_query_for_llm, callbacks=[st_callback_analysis])
                    if llm_response_data and llm_response_data.get('analysis') and not analysis_display_area.empty:
                        analysis_display_area.markdown(llm_response_data.get('analysis'))
                    
                    st.session_state.analysis_result_text = llm_response_data.get('analysis', "No analysis text returned.")

                except Exception as e_analysis:
                    st.error(f"An error occurred during LLM analysis: {str(e_analysis)}")
                    st.session_state.analysis_result_text = "Error during analysis."
        
        if st.session_state.analysis_result_text:
            st.markdown('<div class="result-area">', unsafe_allow_html=True)
            st.markdown(st.session_state.analysis_result_text)
            st.markdown('</div>', unsafe_allow_html=True)
        
        # Follow-up question box
        st.markdown('<div class="query-box">', unsafe_allow_html=True)
        st.markdown("### Ask About This Reaction")
        
        current_rxn_query = st.text_area(
            "Enter your question about this reaction:",
            value=st.session_state.rxn_query,
            height=80,
            key="rxn_query_input_follow_up"
        )
        
        if st.button("Ask Question", key="ask_rxn_btn_follow_up"):
            if current_rxn_query:
                st.session_state.rxn_query = current_rxn_query
                with st.spinner("Processing your question..."):
                    combined_query_for_llm = f"Regarding the reaction {st.session_state.selected_rxn_smiles}, please answer: {st.session_state.rxn_query}"
                    
                    follow_up_display_area = st.container()
                    st_callback_follow_up = StreamlitCallbackHandler(follow_up_display_area)
                    try:
                        follow_up_response_data = enhanced_query(combined_query_for_llm, callbacks=[st_callback_follow_up])
                        if follow_up_response_data and follow_up_response_data.get('analysis'):
                             follow_up_display_area.markdown("### Answer")
                             follow_up_display_area.markdown('<div class="result-area">', unsafe_allow_html=True)
                             follow_up_display_area.markdown(follow_up_response_data.get('analysis'))
                             follow_up_display_area.markdown('</div>', unsafe_allow_html=True)
                    except Exception as e_follow_up:
                        follow_up_display_area.error(f"An error occurred: {str(e_follow_up)}")
            else:
                st.warning("Please enter a question.")
        st.markdown('</div>', unsafe_allow_html=True)
            
    if st.button("← Back to Retrosynthesis Results", key="back_btn"):
        st.session_state.show_analyze_section = False
        st.session_state.selected_rxn_smiles = None
        st.session_state.analysis_result_text = ""
        st.session_state.rxn_query = ""
        st.rerun()
            
    st.markdown('</div>', unsafe_allow_html=True)

# Footer
st.markdown("---")
st.caption("ChemCopilot - Your Expert Chemistry Assistant", unsafe_allow_html=True)
