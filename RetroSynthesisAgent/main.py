# # --- Add this block at the TOP of main.py ---
# # import ssl
# # import os

# # WARNING: Disables SSL certificate verification globally!
# # Use only if necessary and you understand the risks (e.g., trusted network).
# # if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
# #         getattr(ssl, '_create_unverified_context', None)):
# #     print("!!! WARNING: SSL CERTIFICATE VERIFICATION DISABLED GLOBALLY !!!")
# #     ssl._create_default_https_context = ssl._create_unverified_context
# # --- End of SSL bypass block ---

# import json
# from RetroSynAgent.treeBuilder import Tree, TreeLoader
# from RetroSynAgent.pdfProcessor import PDFProcessor
# from RetroSynAgent.knowledgeGraph import KnowledgeGraph
# from RetroSynAgent import prompts
# from RetroSynAgent.GPTAPI import GPTAPI
# from RetroSynAgent.pdfDownloader import PDFDownloader
# import fitz  # PyMuPDF
# import os
# import json
# import re
# import pubchempy
# from RetroSynAgent.entityAlignment import EntityAlignment
# from RetroSynAgent.treeExpansion import TreeExpansion
# from RetroSynAgent.reactionsFiltration import ReactionsFiltration
# import argparse

# def parse_reaction_data(raw_text: str) -> dict:
#     # 1. Extract recommended pathway
#     rec_match = re.search(r"Recommended Reaction Pathway:\s*([^\n]+)", raw_text)
#     recommended = [idx.strip() for idx in rec_match.group(1).split(",")] if rec_match else []

#     # 2. Extract the Reasons block (everything after "Reasons:")
#     reasons = ""
#     reasons_match = re.search(r"Reasons:\s*((?:.|\n)*)", raw_text)
#     if reasons_match:
#         reasons = reasons_match.group(1).strip()

#     # 3. Split into individual reaction blocks
#     blocks = re.split(r"(?=Reaction idx:)", raw_text)
#     reactions = []
#     for blk in blocks:
#         if not blk.strip().startswith("Reaction idx:"):
#             continue

#         idx_match    = re.search(r"Reaction idx:\s*(\S+)", blk)
#         react_match  = re.search(r"Reactants:\s*(.+)", blk)
#         prod_match   = re.search(r"Products:\s*(.+)", blk)
#         smile_match  = re.search(r"Reaction SMILES:\s*(\S+)", blk)
#         cond_match   = re.search(r"Conditions:\s*(.+)", blk)
#         source_match = re.search(r"Source:\s*(.+)", blk)
#         link_match   = re.search(r"SourceLink:\s*\[?(.+?)\]?(?:\s|$)", blk)

#         reaction = {
#             "idx":        idx_match.group(1) if idx_match else None,
#             "reactants":  [r.strip() for r in react_match.group(1).split(",")] if react_match else [],
#             "products":   [p.strip() for p in prod_match.group(1).split(",")]  if prod_match else [],
#             "smiles":     smile_match.group(1) if smile_match else None,
#             "conditions": {},
#             "source":     source_match.group(1).strip() if source_match else None,
#             "source_link": link_match.group(1).strip() if link_match else None
#         }

#         # parse conditions into key/value pairs
#         if cond_match:
#             for part in cond_match.group(1).split(","):
#                 if ":" in part:
#                     key, val = part.split(":", 1)
#                     reaction["conditions"][key.strip().lower()] = val.strip()

#         reactions.append(reaction)

#     return {
#         "recommended_pathway": recommended,
#         "reactions": reactions,
#         "reasons": reasons
#     }

# def parse_arguments():
#     """
#     Parse command-line arguments.
#     """
#     parser = argparse.ArgumentParser(description="Process PDFs and extract reactions.")
#     parser.add_argument('--material', type=str, required=True,
#                         help="Material name for processing.")
#     parser.add_argument('--num_results', type=int, required=True,
#                         help="Number of PDF to download.")
#     parser.add_argument('--alignment', type=str, default="False", choices=["True", "False"],
#                         help="Whether to align entities except for root node.")
#     parser.add_argument('--expansion', type=str, default="False", choices=["True", "False"],
#                         help="Whether to expand the tree with additional literature.")
#     parser.add_argument('--filtration', type=str, default="False", choices=["True", "False"],
#                         help="Whether to filter reactions.")
#     return parser.parse_args()

# def countNodes(tree):
#     node_count = tree.get_node_count()
#     return node_count

# def searchPathways(tree):
#     all_path = tree.find_all_paths()
#     return all_path


# def recommendReactions(prompt, result_folder_name, response_name):
#     res = GPTAPI().answer_wo_vision(prompt)
#     with open(f'{result_folder_name}/{response_name}.txt', 'w') as f:
#         f.write(res)
#     start_idx = res.find("Recommended Reaction Pathway:")
#     recommend_reactions_txt = res[start_idx:]
#     print(f'\n=================================================='
#           f'==========\n{recommend_reactions_txt}\n====================='
#           f'=======================================\n')
#     return recommend_reactions_txt

# def main():
#     # material = 'Polyimide'
#     # num_results = 10
#     # alignment = True
#     # expansion = True
#     # filtration = False

#     # Parse command-line arguments
#     args = parse_arguments()
#     material = args.material
#     num_results = args.num_results
#     # turn str to bool
#     alignment = args.alignment == "True"
#     expansion = args.alignment == "True"
#     filtration = args.filtration == "True"

#     pdf_folder_name = 'pdf_pi'
#     result_folder_name = 'res_pi'
#     result_json_name = 'llm_res'
#     tree_folder_name = 'tree_pi'
#     os.makedirs(tree_folder_name, exist_ok=True)
#     entityalignment = EntityAlignment()
#     treeloader = TreeLoader()
#     tree_expansion = TreeExpansion()
#     reactions_filtration = ReactionsFiltration()

#     ### extractInfos

#     # 1  query literatures & download
#     downloader = PDFDownloader(material, pdf_folder_name=pdf_folder_name, num_results=num_results, n_thread=3)
#     pdf_name_list = downloader.main()
#     print(f'successfully downloaded {len(pdf_name_list)} pdfs for {material}')

#     # 2 Extract infos from PDF about reactions
#     pdf_processor = PDFProcessor(pdf_folder_name=pdf_folder_name, result_folder_name=result_folder_name,
#                                  result_json_name=result_json_name)
#     pdf_processor.load_existing_results()
#     pdf_processor.process_pdfs_txt(save_batch_size=2)

#     ### treeBuildWOExapnsion
#     results_dict = entityalignment.alignRootNode(result_folder_name, result_json_name, material)

#     # 4 construct kg & tree
#     tree_name_wo_exp = tree_folder_name + '/' + material + '_wo_exp.pkl'
#     if not os.path.exists(tree_name_wo_exp):
#         tree_wo_exp = Tree(material.lower(), result_dict=results_dict)
#         print('Starting to construct RetroSynthetic Tree...')
#         tree_wo_exp.construct_tree()
#         treeloader.save_tree(tree_wo_exp, tree_name_wo_exp)
#     else:
#         tree_wo_exp = treeloader.load_tree(tree_name_wo_exp)
#         print('RetroSynthetic Tree wo expansion already loaded.')
#     node_count_wo_exp = countNodes(tree_wo_exp)
#     all_path_wo_exp = searchPathways(tree_wo_exp)
#     print(f'The tree contains {node_count_wo_exp} nodes and {len(all_path_wo_exp)} pathways before expansion.')

#     if alignment:
#         print('Starting to align the nodes of RetroSynthetic Tree...')

#         ### WO Expansion
#         tree_name_wo_exp_alg = tree_folder_name + '/' + material + '_wo_exp_alg.pkl'
#         if not os.path.exists(tree_name_wo_exp_alg):
#             # reactions_wo_exp_alg = entityalignment.entityAlignment(tree_wo_exp.reactions)
#             # tree_wo_exp_alg = Tree(material.lower(), reactions=reactions_wo_exp_alg)
#             reactions_wo_exp = tree_wo_exp.reactions
#             reactions_wo_exp_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_wo_exp)
#             reactions_wo_exp_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_wo_exp_alg_1)
#             tree_wo_exp_alg = Tree(material.lower(), reactions=reactions_wo_exp_alg_all)
#             tree_wo_exp_alg.construct_tree()
#             treeloader.save_tree(tree_wo_exp_alg, tree_name_wo_exp_alg)
#         else:
#             tree_wo_exp_alg = treeloader.load_tree(tree_name_wo_exp_alg)
#             print('aligned RetroSynthetic Tree wo expansion already loaded.')
#         node_count_wo_exp_alg = countNodes(tree_wo_exp_alg)
#         all_path_wo_exp_alg = searchPathways(tree_wo_exp_alg)
#         print(f'The aligned tree contains {node_count_wo_exp_alg} nodes and {len(all_path_wo_exp_alg)} pathways before expansion.')
#         # tree_wo_exp = tree_wo_exp_alg

#     ## treeExpansion
#     # 5 kg & tree expansion
#     results_dict_additional = tree_expansion.treeExpansion(result_folder_name, result_json_name,
#                                                            results_dict, material, expansion=expansion, max_iter=5)
#     if results_dict_additional:
#         results_dict = tree_expansion.update_dict(results_dict, results_dict_additional)
#         # results_dict.update(results_dict_additional)

#     tree_name_exp = tree_folder_name + '/' + material + '_w_exp.pkl'
#     if not os.path.exists(tree_name_exp):
#         tree_exp = Tree(material.lower(), result_dict=results_dict)
#         print('Starting to construct Expanded RetroSynthetic Tree...')
#         tree_exp.construct_tree()
#         treeloader.save_tree(tree_exp, tree_name_exp)
#     else:
#         tree_exp = treeloader.load_tree(tree_name_exp)
#         print('RetroSynthetic Tree w expansion already loaded.')

#     # nodes & pathway count (tree w exp)
#     node_count_exp = countNodes(tree_exp)
#     all_path_exp = searchPathways(tree_exp)
#     print(f'The tree contains {node_count_exp} nodes and {len(all_path_exp)} pathways after expansion.')

#     if alignment:
#         ### Expansion
#         tree_name_exp_alg = tree_folder_name + '/' + material + '_w_exp_alg.pkl'
#         if not os.path.exists(tree_name_exp_alg):
#             reactions_exp = tree_exp.reactions
#             reactions_exp_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_exp)
#             reactions_exp_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_exp_alg_1)
#             tree_exp_alg = Tree(material.lower(), reactions=reactions_exp_alg_all)
#             tree_exp_alg.construct_tree()
#             treeloader.save_tree(tree_exp_alg, tree_name_exp_alg)
#         else:
#             tree_exp_alg = treeloader.load_tree(tree_name_exp_alg)
#             print('aligned RetroSynthetic Tree wo expansion already loaded.')
#         node_count_exp_alg = countNodes(tree_exp_alg)
#         all_path_exp_alg = searchPathways(tree_exp_alg)
#         print(f'The aligned tree contains {node_count_exp_alg} nodes and {len(all_path_exp_alg)} pathways after expansion.')
#         tree_exp = tree_exp_alg

#     all_pathways_w_reactions = reactions_filtration.getFullReactionPathways(tree_exp)

#     ## Filtration
#     if filtration:
#         # filter reactions based on conditions
#         reactions_txt_filtered = reactions_filtration.filterReactions(tree_exp)
#         # build & save tree
#         tree_name_filtered = tree_folder_name + '/' + material + '_filtered' + '.pkl'
#         if not os.path.exists(tree_name_filtered):
#             print('Starting to construct Filtered RetroSynthetic Tree...')
#             tree_filtered = Tree(material.lower(), reactions_txt=reactions_txt_filtered)
#             tree_filtered.construct_tree()
#             treeloader.save_tree(tree_filtered, tree_name_filtered)
#         else:
#             tree_filtered = treeloader.load_tree(tree_name_filtered)
#             print('Filtered RetroSynthetic Tree already loaded.')
#         node_count_filtered = countNodes(tree_filtered)
#         all_path_filtered = searchPathways(tree_filtered)
#         print(f'The tree contains {node_count_filtered} nodes and {len(all_path_filtered)} pathways after filtration.')

#         # filter invalid pathways
#         filtered_pathways = reactions_filtration.filterPathways(tree_filtered)
#         all_pathways_w_reactions = filtered_pathways

#     ### Recommendation
#     # recommend based on specific criterion

#     if not all_pathways_w_reactions.strip(): # Check if the string is empty or just whitespace
#         print(f"\n============================================================")
#         print(f"No valid reaction pathways found or remaining for {material} to recommend.")
#         print(f"============================================================")
#         recommend_reactions_txt = f"No valid reaction pathways found for {material}."

#     else:
#         # Use a general prompt that accepts the material name
#         print(f"Requesting recommendation for {material} based on found pathways...") # Added print
#         prompt_recommend_general = prompts.recommend_prompt_template_general.format(
#             substance=material,  # Pass the actual material name
#             all_pathways=all_pathways_w_reactions
#         )
#         # Give the output file a more general name if needed
#         recommend_reactions_txt = recommendReactions(
#             prompt_recommend_general,
#             result_folder_name,
#             response_name='recommend_pathway_general' # Changed filename
#         )

#     # Optional: Parse the final recommendation text if needed downstream
#     # parsed_data = parse_reaction_data(recommend_reactions_txt)
#     # print(parsed_data)

#     return recommend_reactions_txt

#     # [1]
#     # prompt_recommend1 = prompts.recommend_prompt_commercial.format(all_pathways = all_pathways_w_reactions)
#     # recommend1_reactions_txt = recommendReactions(prompt_recommend1, result_folder_name, response_name='recommend_pathway1')
#     # # parsed_data = parse_reaction_data(recommend1_reactions_txt)
#     # # tree_pathway1 = Tree(material.lower(), reactions_txt=recommend1_reactions_txt)
#     # # print('Starting to construct recommended pathway 1 ...')
#     # # tree_pathway1.construct_tree()
#     # # tree_name_pathway1 = tree_folder_name + '/' + material + '_pathway1' + '.pkl'
#     # # treeloader.save_tree(tree_pathway1, tree_name_pathway1)
#     # # print(parsed_data) 
#     # return recommend1_reactions_txt 

# if __name__ == '__main__':
#     main()




# --- Add this block at the TOP of main.py if needed ---
# import ssl
# import os

# WARNING: Disables SSL certificate verification globally!
# Use only if necessary and you understand the risks (e.g., trusted network).
# if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
#         getattr(ssl, '_create_unverified_context', None)):
#     print("!!! WARNING: SSL CERTIFICATE VERIFICATION DISABLED GLOBALLY !!!")
#     ssl._create_default_https_context = ssl._create_unverified_context
# --- End of SSL bypass block ---

import json
import os
import argparse
import re # Keep re if needed elsewhere

# Import custom modules from RetroSynAgent package
try:
    from RetroSynAgent.treeBuilder import Tree, TreeLoader
    from RetroSynAgent.pdfProcessor import PDFProcessor
    # from RetroSynAgent.knowledgeGraph import KnowledgeGraph # Uncomment if used
    from RetroSynAgent import prompts # Make sure prompts.py is accessible
    from RetroSynAgent.GPTAPI import GPTAPI
    from RetroSynAgent.pdfDownloader import PDFDownloader
    from RetroSynAgent.entityAlignment import EntityAlignment
    from RetroSynAgent.treeExpansion import TreeExpansion
    from RetroSynAgent.reactionsFiltration import ReactionsFiltration
except ImportError as e:
    print(f"Error importing RetroSynAgent modules: {e}")
    print("Please ensure the RetroSynAgent package is installed correctly and accessible.")
    exit(1)

# Optional imports, uncomment if used directly
# import fitz  # PyMuPDF
# import pubchempy


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Process PDFs, build reaction tree, and recommend top 3 pathways.")
    parser.add_argument('--material', type=str, required=True,
                        help="Target material name for processing.")
    parser.add_argument('--num_results', type=int, required=True,
                        help="Number of PDF search results to attempt downloading.")
    parser.add_argument('--alignment', type=str, default="False", choices=["True", "False"],
                        help="Whether to perform entity alignment on tree nodes (default: False).")
    parser.add_argument('--expansion', type=str, default="False", choices=["True", "False"],
                        help="Whether to expand the tree with additional literature searches (default: False).")
    parser.add_argument('--filtration', type=str, default="False", choices=["True", "False"],
                        help="Whether to filter reactions and pathways based on predefined criteria (default: False).")
    parser.add_argument('--recommend_type', type=str, default="general", choices=["general", "cost", "condition"],
                        help="Criterion for recommending top 3 pathways (default: general).")
    return parser.parse_args()

def countNodes(tree):
    """Counts the number of nodes in the tree."""
    if not tree or not hasattr(tree, 'get_node_count'):
        return 0
    try:
        return tree.get_node_count()
    except Exception as e:
        print(f"Warning: Could not count nodes in tree - {e}")
        return 0


def searchPathways(tree):
    """Finds all root-to-leaf pathways in the tree."""
    if not tree or not hasattr(tree, 'find_all_paths'):
        return []
    try:
        return tree.find_all_paths()
    except Exception as e:
        print(f"Warning: Could not find pathways in tree - {e}")
        return []

def recommendReactions(prompt, result_folder_name, response_name):
    """
    Sends the recommendation prompt to the LLM API, saves the full response,
    prints it, and returns the response text.
    """
    print(f"\nSending recommendation prompt to API (saving as {response_name}.txt)...")
    try:
        # Instantiate GPTAPI - ensure API keys/config are handled within GPTAPI class
        gpt_api = GPTAPI()
        res = gpt_api.answer_wo_vision(prompt)
    except Exception as e:
        print(f"Error calling GPT API: {e}")
        return f"Error calling GPT API: {e}" # Return error message

    output_filepath = os.path.join(result_folder_name, f'{response_name}.txt')
    try:
        # Save the full response
        with open(output_filepath, 'w', encoding='utf-8') as f:
            f.write(res)
        print(f"LLM response saved successfully to: {output_filepath}")
    except Exception as e:
        print(f"Error saving LLM response to {output_filepath}: {e}")
        # Continue even if saving fails, but the response is still available in 'res'

    # Print the full response to console
    print(f"\n=================== Full LLM Response ({response_name}) ===================")
    print(res)
    print(f"============================= End of LLM Response ================================\n")
    return res # Return the full response text


def main():
    """
    Main execution function for the retrosynthesis agent workflow.
    """
    args = parse_arguments()
    material = args.material
    num_results = args.num_results
    alignment = args.alignment == "True"
    expansion = args.expansion == "True"
    filtration = args.filtration == "True"
    recommend_type = args.recommend_type

    # --- Folder/file names ---
    # Create a safe directory name from the material name
    material_safe_name = material.lower().replace(" ", "_").replace("/", "-").replace("\\", "-")
    # Limit length and remove potentially problematic chars further if needed
    material_safe_name = re.sub(r'[^\w\-]+', '', material_safe_name)[:50] # Keep alphanumeric, hyphen, underscore

    base_folder = f"results_{material_safe_name}" # Base folder for the material
    pdf_folder_name = os.path.join(base_folder, 'pdfs')
    result_folder_name = os.path.join(base_folder, 'llm_outputs') # For LLM extraction/recommendation outputs
    tree_folder_name = os.path.join(base_folder, 'trees') # For saved tree objects (.pkl)

    # Create directories if they don't exist
    os.makedirs(pdf_folder_name, exist_ok=True)
    os.makedirs(result_folder_name, exist_ok=True)
    os.makedirs(tree_folder_name, exist_ok=True)

    # Consistent name for the JSON file holding extracted reactions from PDFs
    result_json_name = 'pdf_extraction_results.json' # Use .json extension
    result_json_path = os.path.join(result_folder_name, result_json_name)

    # --- Initialize components ---
    try:
        entityalignment = EntityAlignment()
        treeloader = TreeLoader()
        tree_expansion = TreeExpansion()
        reactions_filtration = ReactionsFiltration()
        # GPTAPI is instantiated within recommendReactions
    except NameError as e:
         print(f"Error initializing components: {e}. Make sure all classes are defined and imported.")
         return f"Initialization error: {e}"
    except Exception as e:
         print(f"Unexpected error during component initialization: {e}")
         return f"Initialization error: {e}"


    ### Step 1: PDF Download ###
    print(f"--- Step 1: Starting PDF Download for '{material}' ---")
    try:
        # Limit threads to avoid overwhelming servers or getting blocked
        downloader = PDFDownloader(material, pdf_folder_name=pdf_folder_name, num_results=num_results, n_thread=min(num_results, 5))
        pdf_name_list = downloader.main() # Should return list of downloaded file paths or names
        print(f'Attempted download for {len(pdf_name_list)} potential PDFs.')
    except Exception as e:
        print(f"Error during PDF download: {e}")
        pdf_name_list = [] # Assume failure

    # Verify if any PDFs exist in the target folder
    try:
        actual_pdfs = [f for f in os.listdir(pdf_folder_name) if f.lower().endswith('.pdf')]
        if not actual_pdfs:
            print(f"Warning: No PDF files found in '{pdf_folder_name}' after download attempt. Cannot proceed with extraction.")
            # Decide whether to exit or continue (maybe user placed PDFs manually)
            # For now, let's allow continuing but expect potential errors later
            # return "No PDFs found to process."
        else:
             print(f"Found {len(actual_pdfs)} PDF files in '{pdf_folder_name}'.")
    except FileNotFoundError:
         print(f"Error: PDF directory '{pdf_folder_name}' not found or inaccessible.")
         return "PDF directory error."


    ### Step 2: PDF Processing & Reaction Extraction ###
    print(f"\n--- Step 2: Starting PDF Processing and Reaction Extraction ---")
    # Ensure the PDF folder exists before proceeding
    if not os.path.isdir(pdf_folder_name):
         print(f"Skipping PDF processing: Directory '{pdf_folder_name}' does not exist.")
    else:
        try:
            pdf_processor = PDFProcessor(pdf_folder_name=pdf_folder_name, result_folder_name=result_folder_name,
                                         result_json_name=result_json_name) # Pass only the filename
            pdf_processor.load_existing_results()
            # Ensure process_pdfs_txt handles potential empty folder or non-PDF files gracefully
            pdf_processor.process_pdfs_txt(save_batch_size=2) # Check function signature if needed
            print(f"PDF processing complete. Extracted reactions data stored in: {result_json_path}")
        except Exception as e:
            print(f"Error during PDF processing: {e}")
            # Continue, but subsequent steps might fail if result_json_path is missing/empty


    ### Step 3: Align Root Node & Load Initial Reactions ###
    print(f"\n--- Step 3: Aligning Root Node ('{material}') and Loading Initial Reactions ---")
    results_dict = None # Initialize
    try:
        # alignRootNode should return the dict of reactions keyed by source PDF/ID
        # It needs to handle the case where the JSON file doesn't exist or is empty
        results_dict = entityalignment.alignRootNode(result_folder_name, result_json_name, material)
        if not results_dict:
             print(f"Warning: No reactions found in '{result_json_path}' after attempting root node alignment.")
             # No point continuing if there are no reactions
             # return "No reactions loaded."
        else:
             print(f"Successfully loaded and aligned root node reactions from {len(results_dict)} sources.")
    except FileNotFoundError:
        print(f"Error: Results file not found at '{result_json_path}'. Cannot load reactions.")
        # return f"Results file not found: {result_json_path}"
    except json.JSONDecodeError:
         print(f"Error: Could not decode JSON from '{result_json_path}'. File might be corrupted.")
         # return "Invalid JSON results file."
    except Exception as e:
        print(f"Error during root node alignment or loading results: {e}")
        # return f"Error processing results: {e}"


    ### Step 4: Initial Tree Construction ###
    print(f"\n--- Step 4: Constructing Initial Reaction Tree ---")
    current_tree = None # Initialize
    tree_name_base = os.path.join(tree_folder_name, f'{material_safe_name}_base.pkl')

    # Try loading existing tree first only if results_dict was NOT newly loaded/modified
    # For simplicity, let's always try to load first, then build if load fails or if results_dict has content
    tree_loaded = False
    if os.path.exists(tree_name_base):
        try:
            current_tree = treeloader.load_tree(tree_name_base)
            print(f'Initial RetroSynthetic Tree loaded from: {tree_name_base}')
            tree_loaded = True
        except Exception as e:
            print(f"Warning: Could not load base tree from {tree_name_base}: {e}. Will attempt to rebuild.")

    if not tree_loaded:
        if results_dict: # Only build if we have reaction data
            try:
                current_tree = Tree(material.lower(), result_dict=results_dict)
                print('Constructing initial RetroSynthetic Tree...')
                current_tree.construct_tree() # This might raise errors if results_dict is invalid
                treeloader.save_tree(current_tree, tree_name_base)
                print(f"Initial tree constructed and saved to: {tree_name_base}")
            except Exception as e:
                print(f"Error constructing or saving initial tree: {e}")
                current_tree = None # Ensure tree is None if construction failed
        else:
             print("Skipping initial tree construction as no reaction data was loaded.")


    node_count_base = countNodes(current_tree)
    all_path_base = searchPathways(current_tree)
    print(f'Initial tree state: {node_count_base} nodes, {len(all_path_base)} pathways.')

    ### Step 5: Optional Alignment (on Base Tree) ###
    if alignment and current_tree:
        print('\n--- Step 5: Performing Node Alignment (Initial Tree) ---')
        tree_name_base_alg = os.path.join(tree_folder_name, f'{material_safe_name}_base_alg.pkl')
        tree_loaded_aligned = False

        if os.path.exists(tree_name_base_alg):
             try:
                 current_tree = treeloader.load_tree(tree_name_base_alg) # Load existing aligned tree
                 print(f'Aligned initial RetroSynthetic Tree loaded from: {tree_name_base_alg}')
                 tree_loaded_aligned = True
             except Exception as e:
                  print(f"Warning: Could not load aligned base tree from {tree_name_base_alg}: {e}. Will attempt alignment.")

        if not tree_loaded_aligned:
             print("Performing alignment...")
             try:
                 # Get reactions from current (potentially unaligned) tree
                 reactions_base = current_tree.reactions
                 if not reactions_base:
                      print("Warning: No reactions found in the current tree to perform alignment on.")
                 else:
                      reactions_base_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_base)
                      reactions_base_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_base_alg_1)

                      # Create new tree with aligned reactions
                      aligned_base_tree = Tree(material.lower(), reactions=reactions_base_alg_all)
                      aligned_base_tree.construct_tree()
                      treeloader.save_tree(aligned_base_tree, tree_name_base_alg)
                      print(f"Alignment complete, saving aligned initial tree to: {tree_name_base_alg}")
                      current_tree = aligned_base_tree # Update current tree
             except AttributeError:
                  print("Error: Alignment functions (entityAlignment_1/2) not found or reactions attribute missing.")
             except Exception as e:
                  print(f"Error during initial tree alignment: {e}. Using previous tree state.")


        node_count_alg = countNodes(current_tree)
        all_path_alg = searchPathways(current_tree)
        print(f'Tree state after initial alignment: {node_count_alg} nodes, {len(all_path_alg)} pathways.')

    ### Step 6: Optional Tree Expansion ###
    expanded_tree_processed = False # Flag to know if expansion logic ran
    if expansion:
        print(f"\n--- Step 6: Performing Tree Expansion ---")
        tree_name_exp = os.path.join(tree_folder_name, f'{material_safe_name}_expanded.pkl')
        tree_loaded_expanded = False

        # Try loading existing expanded tree first
        if os.path.exists(tree_name_exp):
            try:
                current_tree = treeloader.load_tree(tree_name_exp)
                print(f'Loaded previously expanded RetroSynthetic Tree from: {tree_name_exp}')
                tree_loaded_expanded = True
                expanded_tree_processed = True # Mark as processed
            except Exception as e:
                print(f"Warning: Could not load expanded tree from {tree_name_exp}: {e}. Will attempt expansion.")

        if not tree_loaded_expanded:
            if not results_dict:
                 print("Cannot perform expansion as initial results dictionary is empty or failed to load.")
            else:
                try:
                     print("Attempting tree expansion...")
                     # treeExpansion might modify results_dict or return additions
                     # Ensure results_dict passed is the latest state if modified by alignment
                     # For simplicity, let's assume it uses the initially loaded results_dict for finding leaves
                     results_dict_additional = tree_expansion.treeExpansion(
                         result_folder_name, result_json_name, # Pass folder/json name if needed by expansion module
                         results_dict, material, expansion=True, max_iter=5
                     )

                     if results_dict_additional:
                         print("Updating results dictionary with expanded reactions...")
                         # Assuming update_dict correctly merges new reactions
                         updated_results_dict = tree_expansion.update_dict(results_dict, results_dict_additional)

                         # Force rebuild of the tree after expansion using updated results
                         expanded_tree = Tree(material.lower(), result_dict=updated_results_dict)
                         print('Constructing Expanded RetroSynthetic Tree...')
                         expanded_tree.construct_tree()
                         treeloader.save_tree(expanded_tree, tree_name_exp)
                         print(f"Expanded tree constructed and saved to: {tree_name_exp}")
                         current_tree = expanded_tree # Update current tree
                         results_dict = updated_results_dict # IMPORTANT: Update results_dict for potential later use
                         expanded_tree_processed = True
                     else:
                         print("Expansion process did not yield additional reactions.")
                         expanded_tree_processed = True # Mark as processed even if no additions

                except AttributeError:
                     print("Error: Tree expansion functions (treeExpansion, update_dict) not found.")
                except Exception as e:
                     print(f"Error during tree expansion: {e}")

        if expanded_tree_processed:
            node_count_exp = countNodes(current_tree)
            all_path_exp = searchPathways(current_tree)
            print(f'Tree state after expansion attempt: {node_count_exp} nodes, {len(all_path_exp)} pathways.')

    ### Step 7: Optional Alignment (After Expansion) ###
    # Perform alignment again only if expansion happened AND alignment is enabled
    if alignment and expanded_tree_processed and current_tree:
        print('\n--- Step 7: Performing Node Alignment (Expanded Tree) ---')
        tree_name_exp_alg = os.path.join(tree_folder_name, f'{material_safe_name}_expanded_alg.pkl')
        tree_loaded_exp_aligned = False

        if os.path.exists(tree_name_exp_alg):
             try:
                  current_tree = treeloader.load_tree(tree_name_exp_alg)
                  print(f'Aligned expanded RetroSynthetic Tree loaded from: {tree_name_exp_alg}')
                  tree_loaded_exp_aligned = True
             except Exception as e:
                  print(f"Warning: Could not load aligned expanded tree from {tree_name_exp_alg}: {e}. Will attempt alignment.")

        if not tree_loaded_exp_aligned:
             print("Performing alignment on expanded tree...")
             try:
                 reactions_exp = current_tree.reactions
                 if not reactions_exp:
                      print("Warning: No reactions found in the current expanded tree to perform alignment on.")
                 else:
                      reactions_exp_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_exp)
                      reactions_exp_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_exp_alg_1)

                      aligned_exp_tree = Tree(material.lower(), reactions=reactions_exp_alg_all)
                      aligned_exp_tree.construct_tree()
                      treeloader.save_tree(aligned_exp_tree, tree_name_exp_alg)
                      print(f"Alignment complete, saving aligned expanded tree to: {tree_name_exp_alg}")
                      current_tree = aligned_exp_tree # Update current tree
             except AttributeError:
                  print("Error: Alignment functions (entityAlignment_1/2) not found or reactions attribute missing.")
             except Exception as e:
                  print(f"Error during expanded tree alignment: {e}. Using previous tree state.")

        node_count_exp_alg = countNodes(current_tree)
        all_path_exp_alg = searchPathways(current_tree)
        print(f'Tree state after post-expansion alignment: {node_count_exp_alg} nodes, {len(all_path_exp_alg)} pathways.')


    ### Step 8: Optional Filtration ###
    final_tree = current_tree # This is the tree state before filtration starts
    tree_was_filtered = False
    pathways_were_filtered = False

    if filtration and final_tree: # Proceed only if filtration is enabled and a tree exists
        print(f"\n--- Step 8: Performing Filtration ---")

        # 8a: Filter Reactions
        print("Filtering reactions based on criteria...")
        tree_after_reaction_filter = None
        try:
            # filterReactions should take the tree object and return reaction text for REMAINING reactions
            reactions_txt_filtered = reactions_filtration.filterReactions(final_tree)

            if not reactions_txt_filtered or not reactions_txt_filtered.strip().startswith("Reaction 001"):
                print("Reaction filtration removed all reactions. No pathways remain.")
                # Keep final_tree as the state before filtration, but recommendation will find no paths
            else:
                # 8b: Build Tree from Filtered Reactions
                print("Constructing tree from filtered reactions...")
                tree_name_filtered = os.path.join(tree_folder_name, f'{material_safe_name}_filtered_reactions.pkl')
                # Always rebuild filtered tree as criteria might change
                try:
                    tree_filtered_reactions = Tree(material.lower(), reactions_txt=reactions_txt_filtered)
                    tree_filtered_reactions.construct_tree()
                    treeloader.save_tree(tree_filtered_reactions, tree_name_filtered)
                    print(f"Tree based on filtered reactions saved to: {tree_name_filtered}")
                    # Use this tree for pathway filtering, but final_tree for getFull... might need adjustment
                    tree_after_reaction_filter = tree_filtered_reactions
                    tree_was_filtered = True # Mark that reaction filtering happened
                except Exception as e:
                    print(f"Error constructing or saving tree from filtered reactions: {e}.")
                    # Fallback: continue without pathway filtering based on this new tree

        except AttributeError:
            print("Error: reactions_filtration.filterReactions function not found.")
        except Exception as e:
            print(f"Error during reaction filtration step: {e}")

        # 8c: Filter Pathways (use tree_after_reaction_filter if created)
        tree_to_filter_pathways = tree_after_reaction_filter if tree_after_reaction_filter else final_tree

        if tree_to_filter_pathways:
            node_count_filt = countNodes(tree_to_filter_pathways)
            all_path_filt = searchPathways(tree_to_filter_pathways)
            print(f'Tree state before pathway filtration: {node_count_filt} nodes, {len(all_path_filt)} pathways.')

            print("Filtering pathways based on validity...")
            try:
                # filterPathways should take the tree object
                # Its output format determines how we proceed. Assume it modifies internal state or returns valid pathways
                # For now, let's assume it just prints/logs excluded pathways and we still use getFull... on the tree.
                # If it returns text, we'd capture it.
                # Let's refine this if filterPathways is intended to return the final pathway text
                reactions_filtration.filterPathways(tree_to_filter_pathways) # Call the function
                # We don't necessarily create a *new* tree here, just identify invalid paths
                pathways_were_filtered = True # Mark that pathway filtering logic ran
                print("Pathway filtration logic executed.")
                # Update final_tree only if reaction filtering created a new valid tree
                if tree_after_reaction_filter:
                     final_tree = tree_after_reaction_filter

            except AttributeError:
                print("Error: reactions_filtration.filterPathways function not found.")
            except Exception as e:
                print(f"Error during pathway filtration step: {e}")
        else:
             print("Skipping pathway filtration as no valid tree exists after reaction filtering.")


    elif filtration and not final_tree:
        print("\nSkipping Filtration step as no valid tree exists.")


    ### Step 9: Prepare Pathways for Recommendation ###
    print("\n--- Step 9: Preparing Final Pathways for Recommendation ---")
    all_pathways_w_reactions = ""
    num_final_paths = 0
    if final_tree:
        try:
            # Get the text representation of all valid pathways from the FINAL tree state
            # IMPORTANT: Ensure getFullReactionPathways returns only valid pathways if pathway filtering occurred.
            # This might require modification of getFullReactionPathways or filterPathways.
            # Assuming getFullReactionPathways reflects the current state of `final_tree`
            all_pathways_w_reactions = reactions_filtration.getFullReactionPathways(final_tree)
            num_final_paths = len(searchPathways(final_tree)) # Get count from the final tree state

            if not all_pathways_w_reactions.strip() and num_final_paths > 0:
                print(f"Warning: {num_final_paths} pathways found in tree, but getFullReactionPathways returned empty text.")
                all_pathways_w_reactions = "" # Ensure empty
            elif not all_pathways_w_reactions.strip() and num_final_paths == 0:
                 print("No pathways available in the final tree state.")
            else:
                 print(f"Extracted text for {num_final_paths} pathways from the final tree state for recommendation.")

        except AttributeError:
            print("Error: reactions_filtration.getFullReactionPathways function not found.")
            all_pathways_w_reactions = ""
        except Exception as e:
            print(f"Error generating full reaction pathways text: {e}")
            all_pathways_w_reactions = ""
    else:
        print("No final tree available to extract pathways from.")


    ### Step 10: Recommendation ###
    print(f"\n--- Step 10: Generating Top 3 Recommendation ({recommend_type.capitalize()}) ---")
    recommend_reactions_txt = ""
    if not all_pathways_w_reactions: # Check if the extracted text is empty
        print(f"\n============================================================")
        print(f"No valid reaction pathways found or remaining for '{material}' after processing.")
        print(f"Cannot generate recommendations.")
        print(f"============================================================")
        recommend_reactions_txt = f"No valid reaction pathways found for {material} to recommend."
        # Save this message to a file for consistency
        output_filepath = os.path.join(result_folder_name, f'recommend_pathways_top3_{recommend_type}.txt')
        try:
            with open(output_filepath, 'w', encoding='utf-8') as f:
                f.write(recommend_reactions_txt)
        except Exception as e:
            print(f"Error saving 'no pathways' message: {e}")

    else:
        # Select the appropriate prompt template from prompts.py
        prompt_template = None
        if recommend_type == "general":
            prompt_template = prompts.recommend_prompt_template_general_top3
        elif recommend_type == "cost":
            prompt_template = prompts.recommend_prompt_template_cost_top3
        elif recommend_type == "condition":
            prompt_template = prompts.recommend_prompt_template_condition_top3

        if prompt_template:
            try:
                prompt_recommend = prompt_template.format(
                    substance=material,
                    all_pathways=all_pathways_w_reactions
                )

                # Call the LLM via recommendReactions function
                recommend_reactions_txt = recommendReactions(
                    prompt_recommend,
                    result_folder_name,
                    response_name=f'recommend_pathways_top3_{recommend_type}' # Filename based on type
                )
            except KeyError as e:
                 print(f"Error formatting prompt template: Missing placeholder {e}")
                 recommend_reactions_txt = f"Error: Prompt template formatting failed ({e})."
            except Exception as e:
                 print(f"An unexpected error occurred during recommendation formatting or API call: {e}")
                 recommend_reactions_txt = f"Error during recommendation step: {e}"
        else:
             # This case should not happen due to argparse choices, but as a safeguard:
             print(f"Error: Invalid recommendation type '{recommend_type}' - could not find corresponding prompt template.")
             recommend_reactions_txt = f"Error: Invalid recommendation type selected."


    print("\n--- Script Finished ---")
    # Return the final LLM response text (or error message)
    return recommend_reactions_txt


# Main execution block
if __name__ == '__main__':
    final_output = main()
    # Optional: Add further processing of the final_output if needed
    print("\nScript execution complete.")