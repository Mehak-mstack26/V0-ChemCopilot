#CURRENT 
import time
import traceback
import re
import requests
import json
import sys
import os
from langchain.tools import BaseTool
from typing import Optional, Dict, Any, List, Union

# --- [DEBUG] Script Initialization ---
print("[DEBUG] Initializing retrosynthesis.py script...")
start_time = time.time()

# --- [DEBUG] Setting up Directory Paths ---
try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)  # v0_CHEM_COPILOT directory
    retrosyn_agent_dir = os.path.join(project_root, 'RetroSynthesisAgent', 'RetroSynAgent')
    print(f"[DEBUG] current_dir: {current_dir}")
    print(f"[DEBUG] project_root: {project_root}")
    print(f"[DEBUG] retrosyn_agent_dir: {retrosyn_agent_dir}")

    # --- [DEBUG] Adding Paths to sys.path ---
    print("[DEBUG] Adding directories to sys.path...")
    sys.path.append(project_root)
    sys.path.append(os.path.join(project_root, 'RetroSynthesisAgent'))
    sys.path.append(retrosyn_agent_dir)
    print(f"[DEBUG] sys.path updated: {sys.path}")

except Exception as e:
    print(f"[ERROR] Failed during directory setup: {e}\n{traceback.format_exc()}")
    sys.exit(1) # Exit if basic setup fails

# --- [DEBUG] Importing Tools and Modules ---
print("[DEBUG] Importing required tools and modules...")
try:
    # Import name2smiles tool
    from tools.name2smiles import NameToSMILES
    print("[DEBUG] Imported NameToSMILES")

    # Import necessary modules from RetroSynAgent
    from RetroSynthesisAgent.RetroSynAgent.treeBuilder import Tree, TreeLoader, CommonSubstanceDB
    print("[DEBUG] Imported Tree, TreeLoader, CommonSubstanceDB")
    from RetroSynthesisAgent.RetroSynAgent.pdfProcessor import PDFProcessor
    print("[DEBUG] Imported PDFProcessor")
    from RetroSynthesisAgent.RetroSynAgent.knowledgeGraph import KnowledgeGraph
    print("[DEBUG] Imported KnowledgeGraph")
    from RetroSynthesisAgent.RetroSynAgent import prompts
    print("[DEBUG] Imported prompts")
    from RetroSynthesisAgent.RetroSynAgent.GPTAPI import GPTAPI
    print("[DEBUG] Imported GPTAPI")
    from RetroSynthesisAgent.RetroSynAgent.pdfDownloader import PDFDownloader
    print("[DEBUG] Imported PDFDownloader")
    from RetroSynthesisAgent.RetroSynAgent.entityAlignment import EntityAlignment
    print("[DEBUG] Imported EntityAlignment")
    from RetroSynthesisAgent.RetroSynAgent.treeExpansion import TreeExpansion
    print("[DEBUG] Imported TreeExpansion")
    from RetroSynthesisAgent.RetroSynAgent.reactionsFiltration import ReactionsFiltration
    print("[DEBUG] Imported ReactionsFiltration")
    print("[DEBUG] All imports successful.")

except ImportError as e:
    print(f"[ERROR] Failed to import required modules: {e}")
    print(f"[ERROR] Please ensure all dependencies are installed and paths are correct.")
    print(f"[ERROR] Current sys.path: {sys.path}")
    print(f"[ERROR] Traceback: {traceback.format_exc()}")
    sys.exit(1)

# --- [DEBUG] Monkey Patching CommonSubstanceDB.read_data_from_json ---
print("[DEBUG] Preparing to monkey patch CommonSubstanceDB.read_data_from_json...")
try:
    # Save original method to restore later if needed
    original_read_data_from_json = CommonSubstanceDB.read_data_from_json
    print("[DEBUG] Original CommonSubstanceDB.read_data_from_json saved.")

    # Define a patched version of the method
    def patched_read_data_from_json(self, filename):
        """Patched method to read data from JSON with correct path resolution and debugging"""
        t_start = time.time()
        print(f"[DEBUG][patched_read_data_from_json] Called with filename: {filename}")
        # Map common filenames to their absolute paths
        file_map = {
            'RetroSynAgent/emol.json': os.path.join(retrosyn_agent_dir, 'emol.json'),
            'RetroSynAgent/common_chemicals.json': os.path.join(retrosyn_agent_dir, 'common_chemicals.json'),
            'emol.json': os.path.join(retrosyn_agent_dir, 'emol.json'),
            'common_chemicals.json': os.path.join(retrosyn_agent_dir, 'common_chemicals.json'),
            # Add other JSON files as needed
        }

        # Use the mapped path if available, otherwise use the original path
        actual_path = file_map.get(filename, filename)
        # If filename was already absolute, actual_path should remain absolute
        if not os.path.isabs(actual_path):
             # If still relative, try resolving against retrosyn_agent_dir
             maybe_path = os.path.join(retrosyn_agent_dir, actual_path)
             if os.path.exists(maybe_path):
                 actual_path = maybe_path
             else: # Fallback to resolving against project root if necessary
                 maybe_path_root = os.path.join(project_root, actual_path)
                 if os.path.exists(maybe_path_root):
                     actual_path = maybe_path_root
                 # Keep the original relative path as last resort if no absolute version found

        print(f"[DEBUG][patched_read_data_from_json] Resolved path: {actual_path}")

        try:
            print(f"[DEBUG][patched_read_data_from_json] Attempting to read JSON from: {actual_path}...")
            # Add timeout simulation for file I/O (conceptual - standard open doesn't have timeout)
            # In a real scenario with potential hangs (e.g., network drives),
            # this might involve threading or async I/O.
            # For now, we just log the attempt duration.
            io_start = time.time()
            with open(actual_path, 'r', encoding='utf-8') as file:
                data = json.load(file)
            io_duration = time.time() - io_start
            print(f"[DEBUG][patched_read_data_from_json] Successfully read {actual_path} in {io_duration:.4f} seconds.")
            print(f"[DEBUG][patched_read_data_from_json] Finished in {time.time() - t_start:.4f} seconds.")
            return data
        except FileNotFoundError:
            print(f"[WARN][patched_read_data_from_json] File not found at primary path: {actual_path}")
            # Try alternate locations relative to retrosyn_agent_dir
            alt_paths = [
                os.path.join(retrosyn_agent_dir, os.path.basename(filename)),
                # Add other potential relative locations if necessary
            ]

            for alt_path in alt_paths:
                print(f"[DEBUG][patched_read_data_from_json] Trying alternate path: {alt_path}")
                if os.path.exists(alt_path):
                    try:
                        print(f"[DEBUG][patched_read_data_from_json] Attempting to read JSON from alternate: {alt_path}...")
                        io_start = time.time()
                        with open(alt_path, 'r', encoding='utf-8') as file:
                            data = json.load(file)
                        io_duration = time.time() - io_start
                        print(f"[DEBUG][patched_read_data_from_json] Successfully loaded from alternate: {alt_path} in {io_duration:.4f} seconds.")
                        print(f"[DEBUG][patched_read_data_from_json] Finished in {time.time() - t_start:.4f} seconds.")
                        return data
                    except Exception as e_alt:
                        print(f"[ERROR][patched_read_data_from_json] Error reading alternate path {alt_path}: {str(e_alt)}")
                        continue
                else:
                     print(f"[DEBUG][patched_read_data_from_json] Alternate path {alt_path} does not exist.")

            print(f"[ERROR][patched_read_data_from_json] Could not find or read file {filename} using resolved path {actual_path} or alternates.")
            print(f"[DEBUG][patched_read_data_from_json] Finished (failed) in {time.time() - t_start:.4f} seconds. Returning empty dict.")
            return {}
        except Exception as e:
            print(f"[ERROR][patched_read_data_from_json] Error reading {actual_path}: {str(e)}\n{traceback.format_exc()}")
            print(f"[DEBUG][patched_read_data_from_json] Finished (error) in {time.time() - t_start:.4f} seconds. Returning empty dict.")
            return {}

    # Apply the monkey patch
    CommonSubstanceDB.read_data_from_json = patched_read_data_from_json
    print("[DEBUG] CommonSubstanceDB.read_data_from_json monkey patched successfully.")

except Exception as e:
    print(f"[ERROR] Failed during monkey patching: {e}\n{traceback.format_exc()}")
    # Continue execution if patching fails, but log the error

# --- [DEBUG] Initializing NameToSMILES Tool ---
print("[DEBUG] Initializing NameToSMILES tool...")
try:
    name_to_smiles_tool = NameToSMILES()
    print("[DEBUG] NameToSMILES tool initialized successfully.")
except Exception as e:
    print(f"[ERROR] Failed to initialize NameToSMILES tool: {e}\n{traceback.format_exc()}")
    # Depending on requirements, might exit or allow continuation with this tool potentially broken.
    name_to_smiles_tool = None # Ensure it's defined, even if null

def get_smiles_for_compound(compound_name, timeout=15):
    """Get SMILES notation for a compound name using the NameToSMILES tool with timeout"""
    t_start = time.time()
    print(f"[DEBUG][get_smiles_for_compound] Called for: '{compound_name}' with timeout={timeout}s")
    if not name_to_smiles_tool:
        print("[ERROR][get_smiles_for_compound] NameToSMILES tool not initialized.")
        return None

    try:
        # NOTE: The BaseTool._run itself doesn't directly accept a timeout.
        # The timeout needs to be implemented within the NameToSMILES tool's
        # underlying mechanism (e.g., the requests call it makes).
        # We are adding timing here to monitor duration.
        print(f"[DEBUG][get_smiles_for_compound] Calling name_to_smiles_tool._run('{compound_name}')...")
        result = name_to_smiles_tool._run(compound_name) # Assuming internal timeout handling if any
        duration = time.time() - t_start
        print(f"[DEBUG][get_smiles_for_compound] Received result: '{result}' in {duration:.4f} seconds.")

        if result and "SMILES:" in result:
            # Extract SMILES from the result
            smiles = result.split("SMILES:")[1].split("\n")[0].strip()
            print(f"[DEBUG][get_smiles_for_compound] Extracted SMILES: '{smiles}'. Total time: {time.time() - t_start:.4f}s")
            return smiles
        else:
            print(f"[WARN][get_smiles_for_compound] 'SMILES:' not found in result for '{compound_name}'. Total time: {time.time() - t_start:.4f}s")
            return None
    except Exception as e:
        duration = time.time() - t_start
        print(f"[ERROR][get_smiles_for_compound] Exception during SMILES lookup for '{compound_name}' after {duration:.4f}s: {e}\n{traceback.format_exc()}")
        return None

def parse_reaction_data(raw_text):
    """Parse the reaction data from the raw text output"""
    t_start = time.time()
    print("[DEBUG][parse_reaction_data] Starting parsing of raw text...")
    # print(f"[DEBUG][parse_reaction_data] Raw text input:\n---\n{raw_text}\n---") # Uncomment for verbose input logging

    reactions = []
    recommended = []
    reasons = ""

    try:
        # 1. Extract recommended pathway
        rec_match = re.search(r"Recommended Reaction Pathway:\s*([^\n]+)", raw_text)
        if rec_match:
            recommended = [idx.strip() for idx in rec_match.group(1).split(",")]
            print(f"[DEBUG][parse_reaction_data] Found recommended pathway indices: {recommended}")
        else:
            print("[WARN][parse_reaction_data] Recommended Reaction Pathway section not found.")

        # 2. Extract the Reasons block
        reasons_match = re.search(r"Reasons:\s*((?:.|\n)*)", raw_text)
        if reasons_match:
            reasons = reasons_match.group(1).strip()
            # Limit printing potentially long reasons block
            print(f"[DEBUG][parse_reaction_data] Found Reasons block (length: {len(reasons)}).")
            # print(f"[DEBUG][parse_reaction_data] Reasons content:\n---\n{reasons}\n---") # Uncomment for verbose reasons logging
        else:
            print("[WARN][parse_reaction_data] Reasons section not found.")

        # 3. Split into individual reaction blocks
        # Use lookahead assertion (?=...) to split *before* "Reaction idx:", keeping the delimiter
        blocks = re.split(r"(?=Reaction idx:)", raw_text)
        print(f"[DEBUG][parse_reaction_data] Split raw text into {len(blocks)} potential blocks.")

        reaction_count = 0
        for i, blk in enumerate(blocks):
            blk_strip = blk.strip()
            if not blk_strip.startswith("Reaction idx:"):
                # Skip parts before the first reaction or empty splits
                # print(f"[DEBUG][parse_reaction_data] Skipping block {i} (doesn't start with 'Reaction idx:' or is empty).")
                continue

            reaction_count += 1
            print(f"[DEBUG][parse_reaction_data] Processing reaction block {reaction_count}...")
            # print(f"[DEBUG][parse_reaction_data] Block content:\n---\n{blk_strip}\n---") # Uncomment for verbose block logging

            idx_match    = re.search(r"Reaction idx:\s*(\S+)", blk)
            react_match  = re.search(r"Reactants:\s*(.+)", blk)
            prod_match   = re.search(r"Products:\s*(.+)", blk)
            smile_match  = re.search(r"Reaction SMILES:\s*(\S+)", blk)
            cond_match   = re.search(r"Conditions:\s*(.+)", blk)
            source_match = re.search(r"Source:\s*(.+)", blk)
            link_match   = re.search(r"SourceLink:\s*\[?(.+?)\]?(?:\s|$)", blk) # Made link capture non-greedy and handle optional brackets

            reaction = {
                "idx": idx_match.group(1).strip() if idx_match else f"UNKNOWN_{reaction_count}",
                "reactants": [r.strip() for r in react_match.group(1).split(",")] if react_match else [],
                "products": [p.strip() for p in prod_match.group(1).split(",")]  if prod_match else [],
                "reaction_smiles": smile_match.group(1).strip() if smile_match else None,
                "conditions": {},
                "source": source_match.group(1).strip() if source_match else None,
                "source_link": link_match.group(1).strip() if link_match else None
            }
            print(f"[DEBUG][parse_reaction_data]  - idx: {reaction['idx']}")
            print(f"[DEBUG][parse_reaction_data]  - reactants: {reaction['reactants']}")
            print(f"[DEBUG][parse_reaction_data]  - products: {reaction['products']}")
            print(f"[DEBUG][parse_reaction_data]  - reaction_smiles: {reaction['reaction_smiles']}")
            print(f"[DEBUG][parse_reaction_data]  - source: {reaction['source']}")
            print(f"[DEBUG][parse_reaction_data]  - source_link: {reaction['source_link']}")


            # parse conditions into key/value pairs
            if cond_match:
                conditions_str = cond_match.group(1).strip()
                print(f"[DEBUG][parse_reaction_data]  - raw conditions: '{conditions_str}'")
                for part in conditions_str.split(","):
                    part = part.strip()
                    if ":" in part:
                        try:
                            key, val = part.split(":", 1)
                            reaction["conditions"][key.strip().lower()] = val.strip()
                        except ValueError:
                            print(f"[WARN][parse_reaction_data] Could not split condition part '{part}' into key:value. Storing as is.")
                            # Store the unparsable part under a generic key or handle as needed
                            reaction["conditions"][f"condition_{len(reaction['conditions'])}"] = part
                    elif part: # Store parts without ':' as standalone conditions if needed
                        print(f"[DEBUG][parse_reaction_data] Condition part '{part}' has no colon. Storing as boolean flag or standalone.")
                        reaction["conditions"][part.lower()] = True # Or store as part: part
                print(f"[DEBUG][parse_reaction_data]  - parsed conditions: {reaction['conditions']}")
            else:
                print(f"[DEBUG][parse_reaction_data]  - no conditions found.")


            reactions.append(reaction)

        parsed_data = {
            "recommended_pathway": recommended,
            "reactions": reactions,
            "reasons": reasons
        }
        duration = time.time() - t_start
        print(f"[DEBUG][parse_reaction_data] Finished parsing {len(reactions)} reactions in {duration:.4f} seconds.")
        return parsed_data

    except Exception as e:
        duration = time.time() - t_start
        print(f"[ERROR][parse_reaction_data] Exception during parsing after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
        # Return partially parsed data or an empty structure
        return {
            "recommended_pathway": recommended,
            "reactions": reactions,
            "reasons": reasons,
            "parsing_error": str(e)
        }


def retrosyn_main(material=None, num_results=10, alignment=True, expansion=True, filtration=False):
    """
    Main function to run the retrosynthesis process with detailed debugging.
    Returns the recommendation text.
    """
    main_start_time = time.time()
    print(f"\n{'='*30} [DEBUG] Entering retrosyn_main {'='*30}")
    print(f"[DEBUG][retrosyn_main] Parameters: material='{material}', num_results={num_results}, alignment={alignment}, expansion={expansion}, filtration={filtration}")

    if not material:
        print("[ERROR][retrosyn_main] 'material' parameter cannot be None.")
        raise ValueError("'material' parameter is required.")

    # --- Folder Setup ---
    step_start_time = time.time()
    print("[DEBUG][retrosyn_main] Setting up folder paths...")
    try:
        # Use replace(" ", "_") and maybe sanitize other characters for folder names
        safe_material_name = re.sub(r'[^\w\-]+', '_', material) # Basic sanitization
        pdf_folder_name = os.path.join(project_root, 'RetroSynthesisAgent', 'pdf_pi', safe_material_name)
        result_folder_name = os.path.join(project_root, 'RetroSynthesisAgent', 'res_pi', safe_material_name)
        tree_folder_name = os.path.join(project_root, 'RetroSynthesisAgent', 'tree_pi', safe_material_name)
        result_json_name = 'llm_res' # Base name for JSON results

        print(f"[DEBUG][retrosyn_main] PDF folder: {pdf_folder_name}")
        print(f"[DEBUG][retrosyn_main] Result folder: {result_folder_name}")
        print(f"[DEBUG][retrosyn_main] Tree folder: {tree_folder_name}")
        print(f"[DEBUG][retrosyn_main] Result JSON base name: {result_json_name}")

        print("[DEBUG][retrosyn_main] Creating directories if they don't exist...")
        os.makedirs(pdf_folder_name, exist_ok=True)
        os.makedirs(result_folder_name, exist_ok=True)
        os.makedirs(tree_folder_name, exist_ok=True)
        print("[DEBUG][retrosyn_main] Directories created/verified.")
        print(f"[DEBUG][retrosyn_main] Folder Setup completed in {time.time() - step_start_time:.4f} seconds.")
    except Exception as e:
        print(f"[ERROR][retrosyn_main] Failed during folder setup: {e}\n{traceback.format_exc()}")
        raise

    # --- Object Initialization ---
    step_start_time = time.time()
    print("[DEBUG][retrosyn_main] Initializing core components...")
    try:
        entityalignment = EntityAlignment()
        print("[DEBUG][retrosyn_main] Initialized EntityAlignment.")
        treeloader = TreeLoader()
        print("[DEBUG][retrosyn_main] Initialized TreeLoader.")
        tree_expansion = TreeExpansion()
        print("[DEBUG][retrosyn_main] Initialized TreeExpansion.")
        reactions_filtration = ReactionsFiltration()
        print("[DEBUG][retrosyn_main] Initialized ReactionsFiltration.")
        # Initialize PDFProcessor here as well
        pdf_processor = PDFProcessor(pdf_folder_name=pdf_folder_name, result_folder_name=result_folder_name,
                                     result_json_name=result_json_name)
        print("[DEBUG][retrosyn_main] Initialized PDFProcessor.")
        print(f"[DEBUG][retrosyn_main] Object Initialization completed in {time.time() - step_start_time:.4f} seconds.")
    except Exception as e:
        print(f"[ERROR][retrosyn_main] Failed during object initialization: {e}\n{traceback.format_exc()}")
        raise

    # === Stage 1: Extract Info ===
    print(f"\n--- [DEBUG][retrosyn_main] Stage 1: Extract Info for '{material}' ---")

    # 1.1 Query Literatures & Download
    step_start_time = time.time()
    print(f"[DEBUG][retrosyn_main] 1.1 Initializing PDFDownloader for '{material}'...")
    pdf_name_list = []
    try:
        downloader = PDFDownloader(material, pdf_folder_name=pdf_folder_name, num_results=num_results, n_thread=3)
        print(f"[DEBUG][retrosyn_main] Calling PDFDownloader.main() with num_results={num_results}...")
        pdf_name_list = downloader.main() # This is blocking
        duration = time.time() - step_start_time
        print(f"[DEBUG][retrosyn_main] PDFDownloader.main() finished in {duration:.4f} seconds.")
        print(f'[DEBUG][retrosyn_main] Successfully downloaded {len(pdf_name_list)} PDFs (or found existing) for {material}. List: {pdf_name_list}')
    except Exception as e:
        duration = time.time() - step_start_time
        print(f"[ERROR][retrosyn_main] Exception during PDF download/retrieval after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
        print("[WARN][retrosyn_main] Proceeding with potentially fewer/no PDFs due to download error.")
        pdf_name_list = [] # Ensure it's an empty list

    # 1.2 Extract Info from PDFs
    step_start_time = time.time()
    print(f"[DEBUG][retrosyn_main] 1.2 Processing PDFs using PDFProcessor...")
    pdf_processing_successful = False
    try:
        print(f"[DEBUG][retrosyn_main] Loading existing results from {result_folder_name}/{result_json_name}_*.json...")
        # This call should ideally ensure self.results is created, even if empty
        pdf_processor.load_existing_results()

        # **FIX:** Ensure results attribute exists before accessing it for the debug print
        if not hasattr(pdf_processor, 'results'):
            print("[WARN][retrosyn_main] pdf_processor.results attribute not found after load_existing_results. Initializing to empty dict.")
            pdf_processor.results = {} # Initialize if it wasn't

        print(f"[DEBUG][retrosyn_main] Found {len(pdf_processor.results)} existing results initially.")

        print(f"[DEBUG][retrosyn_main] Starting PDF processing (PDF -> TXT -> LLM extraction)...")
        # **FIX:** This call is crucial and should run unless loading failed catastrophically
        pdf_processor.process_pdfs_txt(save_batch_size=2) # This is blocking
        duration = time.time() - step_start_time
        print(f"[DEBUG][retrosyn_main] PDFProcessor finished in {duration:.4f} seconds.")
        # Access results again *after* processing
        print(f"[DEBUG][retrosyn_main] Total results after processing: {len(pdf_processor.results)}")
        pdf_processing_successful = True # Mark as successful

    except Exception as e:
        duration = time.time() - step_start_time
        print(f"[ERROR][retrosyn_main] Exception during PDF processing stage after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
        print("[WARN][retrosyn_main] PDF processing failed. Subsequent steps might lack necessary data.")
        # Ensure results attribute exists even after error, maybe empty
        if not hasattr(pdf_processor, 'results'):
             pdf_processor.results = {}


    # === Stage 2: Tree Build without Expansion ===
    print(f"\n--- [DEBUG][retrosyn_main] Stage 2: Tree Build without Expansion for '{material}' ---")
    step_start_time = time.time()
    tree_wo_exp = None
    results_dict = {} # **FIX:** Initialize results_dict before the try block
    if not pdf_processing_successful:
         print("[WARN][retrosyn_main] Skipping Stage 2 (Tree Build) because PDF processing failed.")
    else:
        try:
            print(f"[DEBUG][retrosyn_main] 2.1 Aligning root node using EntityAlignment...")
            # This step requires the JSON file created by pdf_processor.process_pdfs_txt()
            align_start = time.time()
            results_dict = entityalignment.alignRootNode(result_folder_name, result_json_name, material)
            align_duration = time.time() - align_start
            print(f"[DEBUG][retrosyn_main] Root node alignment completed in {align_duration:.4f} seconds. Result keys: {list(results_dict.keys())}")

            # 2.2 Construct Tree (No Expansion)
            tree_name_wo_exp = os.path.join(tree_folder_name, f'{safe_material_name}_wo_exp.pkl')
            print(f"[DEBUG][retrosyn_main] 2.2 Checking for existing tree file: {tree_name_wo_exp}")
            if not os.path.exists(tree_name_wo_exp):
                print(f"[DEBUG][retrosyn_main] Tree file not found. Constructing new tree...")
                tree_wo_exp = Tree(material.lower(), result_dict=results_dict)
                print('[DEBUG][retrosyn_main] Starting tree construction (wo expansion)...')
                construct_start = time.time()
                tree_wo_exp.construct_tree() # This can be computationally intensive
                construct_duration = time.time() - construct_start
                print(f'[DEBUG][retrosyn_main] Tree construction (wo expansion) finished in {construct_duration:.4f} seconds.')
                print(f"[DEBUG][retrosyn_main] Saving constructed tree to {tree_name_wo_exp}...")
                save_start = time.time()
                treeloader.save_tree(tree_wo_exp, tree_name_wo_exp)
                save_duration = time.time() - save_start
                print(f"[DEBUG][retrosyn_main] Tree saved in {save_duration:.4f} seconds.")
            else:
                print(f"[DEBUG][retrosyn_main] Loading existing tree from {tree_name_wo_exp}...")
                load_start = time.time()
                tree_wo_exp = treeloader.load_tree(tree_name_wo_exp)
                load_duration = time.time() - load_start
                print(f'[DEBUG][retrosyn_main] RetroSynthetic Tree (wo expansion) loaded successfully in {load_duration:.4f} seconds.')

            # Analyze tree
            if tree_wo_exp:
                node_count_wo_exp = countNodes(tree_wo_exp)
                all_path_wo_exp = searchPathways(tree_wo_exp)
                print(f'[DEBUG][retrosyn_main] Tree (wo expansion) contains {node_count_wo_exp} nodes and {len(all_path_wo_exp)} pathways.')
            else:
                 print('[WARN][retrosyn_main] Tree object (wo expansion) is None after build/load attempt.')

            print(f"[DEBUG][retrosyn_main] Stage 2 (Tree Build wo Expansion) completed in {time.time() - step_start_time:.4f} seconds.")

        except FileNotFoundError as fnf_error:
             # Specific handling for the file not found error - likely means PDF processing didn't create the output
             duration = time.time() - step_start_time
             print(f"[ERROR][retrosyn_main] FileNotFoundError during Stage 2 (Tree Build wo Expansion) after {duration:.4f} seconds: {fnf_error}")
             print("[ERROR][retrosyn_main] This likely means the required JSON results file from PDF processing is missing.")
             print("[WARN][retrosyn_main] Error in tree building (wo expansion). Subsequent steps might fail or be skipped.")
             tree_wo_exp = None # Ensure it's None if failed
        except Exception as e:
            duration = time.time() - step_start_time
            print(f"[ERROR][retrosyn_main] Exception during Stage 2 (Tree Build wo Expansion) after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
            print("[WARN][retrosyn_main] Error in tree building (wo expansion). Subsequent steps might fail or be skipped.")
            tree_wo_exp = None # Ensure it's None if failed


    # === Stage 3: Alignment (without Expansion) ===
    # (Keep this stage as is, it correctly depends on tree_wo_exp)
    tree_wo_exp_alg = None
    if alignment and tree_wo_exp:
        print(f"\n--- [DEBUG][retrosyn_main] Stage 3: Alignment (without Expansion) for '{material}' ---")
        step_start_time = time.time()
        try:
            print('[DEBUG][retrosyn_main] 3.1 Starting alignment of nodes (wo expansion)...')
            tree_name_wo_exp_alg = os.path.join(tree_folder_name, f'{safe_material_name}_wo_exp_alg.pkl')
            print(f"[DEBUG][retrosyn_main] Checking for existing aligned tree file: {tree_name_wo_exp_alg}")

            if not os.path.exists(tree_name_wo_exp_alg):
                # ... (rest of alignment logic remains the same) ...
                print(f"[DEBUG][retrosyn_main] Aligned tree file not found. Performing alignment...")
                reactions_wo_exp = tree_wo_exp.reactions
                print(f"[DEBUG][retrosyn_main] Extracted {len(reactions_wo_exp)} reactions for alignment.")
                align1_start = time.time()
                reactions_wo_exp_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_wo_exp)
                align1_duration = time.time() - align1_start
                print(f"[DEBUG][retrosyn_main] Alignment step 1 finished in {align1_duration:.4f} seconds. Result count: {len(reactions_wo_exp_alg_1)}")
                align2_start = time.time()
                reactions_wo_exp_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_wo_exp_alg_1)
                align2_duration = time.time() - align2_start
                print(f"[DEBUG][retrosyn_main] Alignment step 2 finished in {align2_duration:.4f} seconds. Result count: {len(reactions_wo_exp_alg_all)}")

                print("[DEBUG][retrosyn_main] Constructing tree from aligned reactions (wo expansion)...")
                tree_wo_exp_alg = Tree(material.lower(), reactions=reactions_wo_exp_alg_all)
                construct_start = time.time()
                tree_wo_exp_alg.construct_tree()
                construct_duration = time.time() - construct_start
                print(f"[DEBUG][retrosyn_main] Aligned tree construction finished in {construct_duration:.4f} seconds.")

                print(f"[DEBUG][retrosyn_main] Saving aligned tree (wo expansion) to {tree_name_wo_exp_alg}...")
                save_start = time.time()
                treeloader.save_tree(tree_wo_exp_alg, tree_name_wo_exp_alg)
                save_duration = time.time() - save_start
                print(f"[DEBUG][retrosyn_main] Aligned tree saved in {save_duration:.4f} seconds.")
            else:
                print(f"[DEBUG][retrosyn_main] Loading existing aligned tree (wo expansion) from {tree_name_wo_exp_alg}...")
                load_start = time.time()
                tree_wo_exp_alg = treeloader.load_tree(tree_name_wo_exp_alg)
                load_duration = time.time() - load_start
                print(f'[DEBUG][retrosyn_main] Aligned RetroSynthetic Tree (wo expansion) loaded successfully in {load_duration:.4f} seconds.')

            # Analyze aligned tree
            if tree_wo_exp_alg:
                 node_count_wo_exp_alg = countNodes(tree_wo_exp_alg)
                 all_path_wo_exp_alg = searchPathways(tree_wo_exp_alg)
                 print(f'[DEBUG][retrosyn_main] Aligned tree (wo expansion) contains {node_count_wo_exp_alg} nodes and {len(all_path_wo_exp_alg)} pathways.')
            else:
                 print('[WARN][retrosyn_main] Aligned tree object (wo expansion) is None after build/load attempt.')
            print(f"[DEBUG][retrosyn_main] Stage 3 (Alignment wo Expansion) completed in {time.time() - step_start_time:.4f} seconds.")

        except Exception as e:
            duration = time.time() - step_start_time
            print(f"[ERROR][retrosyn_main] Exception during Stage 3 (Alignment wo Expansion) after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
            print("[WARN][retrosyn_main] Error during alignment (wo expansion). Subsequent steps might use unaligned data or fail.")
            tree_wo_exp_alg = None # Ensure it's None if failed
    elif not alignment:
        print("\n[DEBUG][retrosyn_main] Skipping Stage 3: Alignment (without Expansion) as alignment=False.")
    elif not tree_wo_exp:
        print("\n[DEBUG][retrosyn_main] Skipping Stage 3: Alignment (without Expansion) as previous tree construction failed or was skipped.")


    # === Stage 4: Tree Expansion ===
    print(f"\n--- [DEBUG][retrosyn_main] Stage 4: Tree Expansion for '{material}' ---")
    step_start_time = time.time()
    tree_exp = None
    # **FIX:** Use the results_dict initialized/populated in Stage 2
    results_dict_expanded = results_dict.copy() # Start with initial results from stage 2 (might be empty if stage 2 failed)

    # Check if expansion is feasible (requires initial results_dict to have content typically)
    if not results_dict:
        print("[WARN][retrosyn_main] Skipping Stage 4 (Tree Expansion) because initial results_dict is empty (likely due to previous errors).")
    else:
        try:
            # 4.1 KG & Tree Expansion (potentially using LLM)
            print(f"[DEBUG][retrosyn_main] 4.1 Performing tree expansion (expansion={expansion})...")
            if expansion:
                 exp_start = time.time()
                 # Pass the potentially empty but initialized results_dict
                 results_dict_additional = tree_expansion.treeExpansion(result_folder_name, result_json_name,
                                                                        results_dict, material, expansion=True, max_iter=5)
                 exp_duration = time.time() - exp_start
                 print(f"[DEBUG][retrosyn_main] Tree expansion step finished in {exp_duration:.4f} seconds.")
                 if results_dict_additional:
                     print(f"[DEBUG][retrosyn_main] Found {len(results_dict_additional)} additional results from expansion.")
                     # **FIX:** Ensure update_dict handles potentially empty base dict
                     results_dict_expanded = tree_expansion.update_dict(results_dict_expanded, results_dict_additional)
                     print(f"[DEBUG][retrosyn_main] Updated results dict size: {len(results_dict_expanded)}")
                 else:
                     print("[DEBUG][retrosyn_main] No additional results found during expansion.")
            else:
                 print("[DEBUG][retrosyn_main] Expansion skipped as expansion=False.")


            # 4.2 Construct Expanded Tree
            tree_name_exp = os.path.join(tree_folder_name, f'{safe_material_name}_w_exp.pkl')
            print(f"[DEBUG][retrosyn_main] 4.2 Checking for existing expanded tree file: {tree_name_exp}")
            if not os.path.exists(tree_name_exp):
                print(f"[DEBUG][retrosyn_main] Expanded tree file not found. Constructing new expanded tree...")
                # Use the potentially expanded results dictionary
                tree_exp = Tree(material.lower(), result_dict=results_dict_expanded)
                print('[DEBUG][retrosyn_main] Starting expanded tree construction...')
                construct_start = time.time()
                tree_exp.construct_tree()
                construct_duration = time.time() - construct_start
                print(f'[DEBUG][retrosyn_main] Expanded tree construction finished in {construct_duration:.4f} seconds.')

                print(f"[DEBUG][retrosyn_main] Saving expanded tree to {tree_name_exp}...")
                save_start = time.time()
                treeloader.save_tree(tree_exp, tree_name_exp)
                save_duration = time.time() - save_start
                print(f"[DEBUG][retrosyn_main] Expanded tree saved in {save_duration:.4f} seconds.")
            else:
                print(f"[DEBUG][retrosyn_main] Loading existing expanded tree from {tree_name_exp}...")
                load_start = time.time()
                tree_exp = treeloader.load_tree(tree_name_exp)
                load_duration = time.time() - load_start
                print(f'[DEBUG][retrosyn_main] RetroSynthetic Tree (w expansion) loaded successfully in {load_duration:.4f} seconds.')

            # Analyze expanded tree
            if tree_exp:
                node_count_exp = countNodes(tree_exp)
                all_path_exp = searchPathways(tree_exp)
                print(f'[DEBUG][retrosyn_main] Tree (w expansion) contains {node_count_exp} nodes and {len(all_path_exp)} pathways.')
            else:
                 print('[WARN][retrosyn_main] Tree object (w expansion) is None after build/load attempt.')

            print(f"[DEBUG][retrosyn_main] Stage 4 (Tree Expansion) completed in {time.time() - step_start_time:.4f} seconds.")

        except Exception as e:
            duration = time.time() - step_start_time
            print(f"[ERROR][retrosyn_main] Exception during Stage 4 (Tree Expansion) after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
            print("[WARN][retrosyn_main] Error during tree expansion. Subsequent steps might use unexpanded data or fail.")
            tree_exp = None # Ensure it's None if failed


    # === Stage 5: Alignment (with Expansion) ===
    # (Keep this stage as is, it correctly depends on tree_exp or tree_wo_exp)
    tree_exp_alg = None
    # Use the latest successfully built/loaded tree (expanded if available, otherwise unexpanded)
    base_tree_for_alignment = tree_exp if tree_exp else tree_wo_exp_alg if tree_wo_exp_alg else tree_wo_exp

    if alignment and base_tree_for_alignment:
        print(f"\n--- [DEBUG][retrosyn_main] Stage 5: Alignment (with Expansion) for '{material}' ---")
        step_start_time = time.time()
        try:
            print('[DEBUG][retrosyn_main] 5.1 Starting alignment of nodes (w expansion)...')
            tree_name_exp_alg = os.path.join(tree_folder_name, f'{safe_material_name}_w_exp_alg.pkl')
            print(f"[DEBUG][retrosyn_main] Checking for existing aligned expanded tree file: {tree_name_exp_alg}")

            if not os.path.exists(tree_name_exp_alg):
                # ... (rest of alignment logic remains the same) ...
                print(f"[DEBUG][retrosyn_main] Aligned expanded tree file not found. Performing alignment...")
                reactions_exp = base_tree_for_alignment.reactions
                print(f"[DEBUG][retrosyn_main] Extracted {len(reactions_exp)} reactions for alignment.")
                align1_start = time.time()
                reactions_exp_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_exp)
                align1_duration = time.time() - align1_start
                print(f"[DEBUG][retrosyn_main] Alignment step 1 finished in {align1_duration:.4f} seconds. Result count: {len(reactions_exp_alg_1)}")
                align2_start = time.time()
                reactions_exp_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_exp_alg_1)
                align2_duration = time.time() - align2_start
                print(f"[DEBUG][retrosyn_main] Alignment step 2 finished in {align2_duration:.4f} seconds. Result count: {len(reactions_exp_alg_all)}")

                print("[DEBUG][retrosyn_main] Constructing tree from aligned expanded reactions...")
                tree_exp_alg = Tree(material.lower(), reactions=reactions_exp_alg_all)
                construct_start = time.time()
                tree_exp_alg.construct_tree()
                construct_duration = time.time() - construct_start
                print(f"[DEBUG][retrosyn_main] Aligned expanded tree construction finished in {construct_duration:.4f} seconds.")

                print(f"[DEBUG][retrosyn_main] Saving aligned expanded tree to {tree_name_exp_alg}...")
                save_start = time.time()
                treeloader.save_tree(tree_exp_alg, tree_name_exp_alg)
                save_duration = time.time() - save_start
                print(f"[DEBUG][retrosyn_main] Aligned expanded tree saved in {save_duration:.4f} seconds.")
            else:
                print(f"[DEBUG][retrosyn_main] Loading existing aligned expanded tree from {tree_name_exp_alg}...")
                load_start = time.time()
                tree_exp_alg = treeloader.load_tree(tree_name_exp_alg)
                load_duration = time.time() - load_start
                print(f'[DEBUG][retrosyn_main] Aligned RetroSynthetic Tree (w expansion) loaded successfully in {load_duration:.4f} seconds.')

            # Analyze aligned expanded tree
            if tree_exp_alg:
                node_count_exp_alg = countNodes(tree_exp_alg)
                all_path_exp_alg = searchPathways(tree_exp_alg)
                print(f'[DEBUG][retrosyn_main] Aligned tree (w expansion) contains {node_count_exp_alg} nodes and {len(all_path_exp_alg)} pathways.')

                # *** CRITICAL: Update tree_exp to the aligned version for subsequent steps ***
                print("[DEBUG][retrosyn_main] Updating the main tree reference to the aligned version.")
                tree_exp = tree_exp_alg # Use the aligned tree going forward
            else:
                 print('[WARN][retrosyn_main] Aligned tree object (w expansion) is None after build/load attempt.')
            print(f"[DEBUG][retrosyn_main] Stage 5 (Alignment w Expansion) completed in {time.time() - step_start_time:.4f} seconds.")

        except Exception as e:
            duration = time.time() - step_start_time
            print(f"[ERROR][retrosyn_main] Exception during Stage 5 (Alignment w Expansion) after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
            print("[WARN][retrosyn_main] Error during alignment (w expansion). Subsequent steps will use the unaligned expanded tree (if available).")
            # tree_exp remains the unaligned expanded tree if alignment failed
    elif not alignment:
        print("\n[DEBUG][retrosyn_main] Skipping Stage 5: Alignment (with Expansion) as alignment=False.")
    elif not base_tree_for_alignment:
         print("\n[DEBUG][retrosyn_main] Skipping Stage 5: Alignment (with Expansion) as no base tree is available.")


    # === Stage 6: Pathway Extraction and Filtration ===
    print(f"\n--- [DEBUG][retrosyn_main] Stage 6: Pathway Extraction and Filtration for '{material}' ---")
    step_start_time = time.time()
    all_pathways_w_reactions = []
    # **FIX:** Determine the final tree more robustly
    final_tree_for_recommendation = tree_exp # Preferentially use the latest (potentially expanded+aligned) tree
    if not final_tree_for_recommendation:
         # Fall back through the possible trees if tree_exp is None
         final_tree_for_recommendation = tree_exp_alg if tree_exp_alg else base_tree_for_alignment # base_tree already checks tree_exp, tree_wo_exp_alg, tree_wo_exp
         if final_tree_for_recommendation:
              print("[WARN][retrosyn_main] Using fallback tree for recommendation due to errors in later stages.")

    if not final_tree_for_recommendation:
         print("[ERROR][retrosyn_main] No valid tree available for pathway extraction. Cannot proceed to recommendation.")
         # Raise the error here to stop execution cleanly before Stage 7
         raise RuntimeError("Failed to build or load any synthesis tree after all stages.")

    try:
        print(f"[DEBUG][retrosyn_main] 6.1 Extracting full reaction pathways from the final tree (type: {type(final_tree_for_recommendation).__name__})...")
        extract_start = time.time()
        all_pathways_w_reactions = reactions_filtration.getFullReactionPathways(final_tree_for_recommendation)
        extract_duration = time.time() - extract_start
        print(f"[DEBUG][retrosyn_main] Extracted {len(all_pathways_w_reactions)} full pathways in {extract_duration:.4f} seconds.")

        # 6.2 Filtration (Optional)
        if filtration:
            print("[DEBUG][retrosyn_main] 6.2 Starting filtration process...")
            filter_start = time.time()
            # ... (rest of filtration logic remains the same) ...
            print("[DEBUG][retrosyn_main]   - Filtering reactions based on conditions...")
            reactions_txt_filtered = reactions_filtration.filterReactions(final_tree_for_recommendation)
            print(f"[DEBUG][retrosyn_main]   - Reaction filtering returned type: {type(reactions_txt_filtered)}")

            tree_name_filtered = os.path.join(tree_folder_name, f'{safe_material_name}_filtered.pkl')
            print(f"[DEBUG][retrosyn_main]   - Checking for existing filtered tree: {tree_name_filtered}")
            tree_filtered = None
            if not os.path.exists(tree_name_filtered):
                print(f"[DEBUG][retrosyn_main]   - Filtered tree not found. Constructing...")
                try:
                    if isinstance(reactions_txt_filtered, str): # Example check
                        tree_filtered = Tree(material.lower(), reactions_txt=reactions_txt_filtered)
                    elif isinstance(reactions_txt_filtered, list): # Example check for list of reaction dicts
                         tree_filtered = Tree(material.lower(), reactions=reactions_txt_filtered)
                    else:
                        print("[ERROR][retrosyn_main]   - Unsupported format returned by filterReactions. Cannot build filtered tree.")
                        raise TypeError(f"filterReactions returned unexpected data type {type(reactions_txt_filtered)} for tree construction.")

                    print('[DEBUG][retrosyn_main]   - Starting filtered tree construction...')
                    construct_filt_start = time.time()
                    tree_filtered.construct_tree()
                    construct_filt_duration = time.time() - construct_filt_start
                    print(f'[DEBUG][retrosyn_main]   - Filtered tree construction finished in {construct_filt_duration:.4f} seconds.')

                    print(f"[DEBUG][retrosyn_main]   - Saving filtered tree to {tree_name_filtered}...")
                    save_filt_start = time.time()
                    treeloader.save_tree(tree_filtered, tree_name_filtered)
                    save_filt_duration = time.time() - save_filt_start
                    print(f"[DEBUG][retrosyn_main]   - Filtered tree saved in {save_filt_duration:.4f} seconds.")

                except Exception as filter_build_e:
                     print(f"[ERROR][retrosyn_main] Failed to build or save filtered tree: {filter_build_e}")
                     tree_filtered = None # Ensure it's None if build fails

            else:
                print(f"[DEBUG][retrosyn_main]   - Loading existing filtered tree from {tree_name_filtered}...")
                try:
                    load_filt_start = time.time()
                    tree_filtered = treeloader.load_tree(tree_name_filtered)
                    load_filt_duration = time.time() - load_filt_start
                    print(f'[DEBUG][retrosyn_main]   - Filtered RetroSynthetic Tree loaded successfully in {load_filt_duration:.4f} seconds.')
                except Exception as filter_load_e:
                     print(f"[ERROR][retrosyn_main] Failed to load filtered tree: {filter_load_e}")
                     tree_filtered = None

            if tree_filtered:
                # Analyze filtered tree
                node_count_filtered = countNodes(tree_filtered)
                all_path_filtered = searchPathways(tree_filtered)
                print(f'[DEBUG][retrosyn_main]   - Filtered tree contains {node_count_filtered} nodes and {len(all_path_filtered)} pathways.')

                # Filter invalid pathways using the *filtered* tree
                print("[DEBUG][retrosyn_main]   - Filtering pathways using the filtered tree...")
                filter_path_start = time.time()
                filtered_pathways = reactions_filtration.filterPathways(tree_filtered)
                filter_path_duration = time.time() - filter_path_start
                print(f"[DEBUG][retrosyn_main]   - Pathway filtering finished in {filter_path_duration:.4f} seconds. Kept {len(filtered_pathways)} pathways.")

                # Update the pathways list to be used for recommendation
                all_pathways_w_reactions = filtered_pathways
                # Update the tree reference if needed (though recommendation uses pathways list)
                final_tree_for_recommendation = tree_filtered # Use filtered tree if successful
            else:
                 print("[WARN][retrosyn_main] Skipping pathway filtering as filtered tree construction/loading failed.")

            filter_duration = time.time() - filter_start
            print(f"[DEBUG][retrosyn_main] Filtration process completed in {filter_duration:.4f} seconds.")
        else:
            print("[DEBUG][retrosyn_main] Skipping filtration as filtration=False.")

        print(f"[DEBUG][retrosyn_main] Stage 6 (Pathway Extraction/Filtration) completed in {time.time() - step_start_time:.4f} seconds.")
        print(f"[DEBUG][retrosyn_main] Final number of pathways for recommendation: {len(all_pathways_w_reactions)}")

    except Exception as e:
        duration = time.time() - step_start_time
        print(f"[ERROR][retrosyn_main] Exception during Stage 6 (Pathway Extraction/Filtration) after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
        print("[WARN][retrosyn_main] Error during pathway extraction/filtration. Recommendation may be based on incomplete data or fail.")
        # all_pathways_w_reactions might be incomplete or empty


    # === Stage 7: Recommendation ===
    # (Keep this stage as is, it correctly depends on all_pathways_w_reactions)
    print(f"\n--- [DEBUG][retrosyn_main] Stage 7: Recommendation for '{material}' ---")
    step_start_time = time.time()
    recommend1_reactions_txt = "Error: Recommendation could not be generated due to previous errors or lack of pathways." # Default error message

    if not all_pathways_w_reactions:
        print("[WARN][retrosyn_main] No pathways available to generate recommendation.")
    else:
        try:
            print("[DEBUG][retrosyn_main] 7.1 Generating recommendation prompt...")
            pathways_str_list = []

            # --- Start: Robust Pathway Formatting ---
            # Check if all_pathways_w_reactions is iterable (list/tuple) or already a string
            if isinstance(all_pathways_w_reactions, (list, tuple)):
                for path_obj in all_pathways_w_reactions:
                     try:
                          # Ensure path_obj is converted to string properly
                          path_str = str(path_obj).strip()
                          if path_str: # Add only non-empty pathway strings
                             pathways_str_list.append(path_str)
                     except Exception as str_e:
                          print(f"[WARN][retrosyn_main] Could not convert pathway object {type(path_obj)} to string: {str_e}")
                          pathways_str_list.append(f"Error converting pathway: {str_e}")
            elif isinstance(all_pathways_w_reactions, str):
                 # If it's already a single string containing multiple pathways (e.g., separated by newlines)
                 print(f"[DEBUG][retrosyn_main] Treating all_pathways_w_reactions as pre-formatted string.")
                 # Split potentially multi-line string into individual pathway strings, removing empty lines
                 pathways_str_list = [line.strip() for line in all_pathways_w_reactions.strip().split('\n') if line.strip()]
            else:
                 print(f"[ERROR][retrosyn_main] Unexpected type for all_pathways_w_reactions: {type(all_pathways_w_reactions)}. Cannot format pathways.")
                 raise TypeError("Invalid format for pathways data.")
            # --- End: Robust Pathway Formatting ---


            pathways_str = "\n\n".join(pathways_str_list).strip() # Join non-empty pathways

            # Check if we actually have pathways to recommend after formatting
            if not pathways_str:
                 print("[WARN][retrosyn_main] Formatted pathways string is empty. Cannot generate recommendation.")
                 # Keep the default error message
                 recommend1_reactions_txt = "Error: No valid pathway data available for recommendation formatting."

            else:
                 # *** THE FIX IS HERE ***
                 # Use the general prompt and format with the *actual* material name
                 print(f"[DEBUG][retrosyn_main] Using general recommendation prompt for substance: '{material}'")
                 prompt_recommend_general = prompts.recommend_prompt_template_general.format(
                     substance=material,  # Pass the actual material name from the function argument
                     all_pathways=pathways_str
                 )
                 # --- End of Fix ---

                 # Limit length for logging if necessary
                 # print(f"[DEBUG][retrosyn_main] Recommendation prompt (sample):\n---\n{prompt_recommend_general[:1000]}...\n---")

                 print("[DEBUG][retrosyn_main] 7.2 Calling recommendation function (LLM)...")
                 rec_start = time.time()
                 # Use a more general response name as well
                 recommend1_reactions_txt = recommendReactions(
                     prompt_recommend_general,
                     result_folder_name,
                     response_name='recommend_pathway_general' # Changed filename to reflect general prompt
                 )
                 rec_duration = time.time() - rec_start
                 print(f"[DEBUG][retrosyn_main] Recommendation received from LLM in {rec_duration:.4f} seconds.")

            print(f"[DEBUG][retrosyn_main] Stage 7 (Recommendation) completed in {time.time() - step_start_time:.4f} seconds.")

        except Exception as e:
            duration = time.time() - step_start_time
            print(f"[ERROR][retrosyn_main] Exception during Stage 7 (Recommendation) after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
            # Ensure error message reflects the stage failure
            recommend1_reactions_txt = f"Error during recommendation generation stage: {str(e)}"


    total_duration = time.time() - main_start_time
    print(f"\n{'='*30} [DEBUG] Exiting retrosyn_main {'='*30}")
    print(f"[DEBUG][retrosyn_main] Total execution time: {total_duration:.4f} seconds.")
    # Ensure the final returned text reflects success or the specific error encountered
    print(f"[DEBUG][retrosyn_main] Returning recommendation text (length: {len(recommend1_reactions_txt)}).")
    # Make sure the default error message is returned if pathways were empty or formatting failed
    if not all_pathways_w_reactions or not pathways_str:
         if "Error:" not in recommend1_reactions_txt: # Avoid overwriting specific errors from try-except block
             recommend1_reactions_txt = "Error: Recommendation could not be generated due to lack of valid pathways."

    return recommend1_reactions_txt

# --- Rest of the script (helper functions, Tool class, main execution block) remains the same ---

# Helper functions (from main.py, with added debug logs)
def countNodes(tree):
    t_start = time.time()
    # print(f"[DEBUG][countNodes] Counting nodes in tree (root: {tree.root if tree else 'N/A'})...") # Reduced verbosity
    if not tree:
        print("[WARN][countNodes] Tree object is None.")
        return 0
    try:
        node_count = tree.get_node_count()
        duration = time.time() - t_start
        # print(f"[DEBUG][countNodes] Found {node_count} nodes in {duration:.4f} seconds.") # Reduced verbosity
        return node_count
    except Exception as e:
        duration = time.time() - t_start
        print(f"[ERROR][countNodes] Exception during node counting after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
        return 0 # Or re-raise

def searchPathways(tree):
    t_start = time.time()
    # print(f"[DEBUG][searchPathways] Searching pathways in tree (root: {tree.root if tree else 'N/A'})...") # Reduced verbosity
    if not tree:
        print("[WARN][searchPathways] Tree object is None.")
        return []
    try:
        all_path = tree.find_all_paths()
        duration = time.time() - t_start
        # print(f"[DEBUG][searchPathways] Found {len(all_path)} pathways in {duration:.4f} seconds.") # Reduced verbosity
        return all_path
    except Exception as e:
        duration = time.time() - t_start
        print(f"[ERROR][searchPathways] Exception during pathway search after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
        return [] # Or re-raise

def recommendReactions(prompt, result_folder_name, response_name, timeout=300):
    t_start = time.time()
    print(f"[DEBUG][recommendReactions] Requesting recommendation. Response name: {response_name}, timeout: {timeout}s")
    # print(f"[DEBUG][recommendReactions] Prompt (first 500 chars):\n---\n{prompt[:500]}\n---") # Log prompt sample

    res = "Error: LLM call failed."
    try:
        print(f"[DEBUG][recommendReactions] Calling GPTAPI().answer_wo_vision...")
        llm_start = time.time()
        api = GPTAPI() # Consider initializing once if reused heavily
        res = api.answer_wo_vision(prompt) # Blocking call
        llm_duration = time.time() - llm_start
        print(f"[DEBUG][recommendReactions] GPTAPI call completed in {llm_duration:.4f} seconds.")

        result_file_path = os.path.join(result_folder_name, f'{response_name}.txt')
        print(f"[DEBUG][recommendReactions] Saving response to {result_file_path}...")
        io_start = time.time()
        os.makedirs(os.path.dirname(result_file_path), exist_ok=True)
        with open(result_file_path, 'w', encoding='utf-8') as f:
            f.write(res)
        io_duration = time.time() - io_start
        print(f"[DEBUG][recommendReactions] Response saved in {io_duration:.4f} seconds.")

        start_idx = res.find("Recommended Reaction Pathway:")
        recommend_reactions_txt = res[start_idx:] if start_idx >= 0 else res
        print(f"[DEBUG][recommendReactions] Extracted recommendation text (starts at index {start_idx}).")
        print(f'\n[INFO][recommendReactions] Recommendation Output:\n======================================================'
              f'==========\n{recommend_reactions_txt}\n====================='
              f'=======================================\n')

        duration = time.time() - t_start
        print(f"[DEBUG][recommendReactions] Finished in {duration:.4f} seconds.")
        return recommend_reactions_txt

    except Exception as e:
        duration = time.time() - t_start
        print(f"[ERROR][recommendReactions] Exception during LLM call or processing after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
        try:
            result_file_path = os.path.join(result_folder_name, f'{response_name}_error.txt')
            os.makedirs(os.path.dirname(result_file_path), exist_ok=True)
            with open(result_file_path, 'w', encoding='utf-8') as f:
                f.write(f"Error during recommendation:\n{e}\n\nTraceback:\n{traceback.format_exc()}")
        except Exception as save_e:
            print(f"[ERROR][recommendReactions] Could not save error details to file: {save_e}")

        return f"Error generating recommendation: {str(e)}"

# Function to run with parameters (wrapper around main)
def run_with_params(material, num_results=10, alignment=True, expansion=True, filtration=False):
    """Run the retrosynthesis with the given parameters and return the result"""
    t_start = time.time()
    print(f"\n[DEBUG][run_with_params] Called for material='{material}'. Params: num_results={num_results}, alignment={alignment}, expansion={expansion}, filtration={filtration}")
    output = None
    try:
        print("[DEBUG][run_with_params] Calling retrosyn_main...")
        output = retrosyn_main(
            material=material,
            num_results=num_results,
            alignment=alignment,
            expansion=expansion,
            filtration=filtration
        )
        duration = time.time() - t_start
        print(f"[DEBUG][run_with_params] retrosyn_main completed in {duration:.4f} seconds.")
        return output
    except Exception as e:
        duration = time.time() - t_start
        # Check if it's the specific RuntimeError we added for clarity
        if isinstance(e, RuntimeError) and "Failed to build or load any synthesis tree" in str(e):
             error_msg = f"Retrosynthesis pipeline failed after {duration:.4f} seconds: {str(e)} (Likely due to errors in early stages like PDF processing or tree building)."
        else:
             error_msg = f"Error running retrosynthesis pipeline after {duration:.4f} seconds: {str(e)}\n{traceback.format_exc()}"
        print(f"[ERROR][run_with_params] {error_msg}")
        raise Exception(f"Error in retrosynthesis pipeline: {str(e)}") # Re-raise simplified error

# Class RetroSynthesisTool remains the same...
class RetroSynthesisTool(BaseTool):
    name: str = "retrosynthesis"
    description: str = "Performs retrosynthesis analysis on a chemical compound. Input must be the chemical name."

    def _run(self, compound_name: str) -> Dict[str, Any]:
        """Run retrosynthesis analysis on the given compound"""
        tool_start_time = time.time()
        print(f"\n{'#'*30} [DEBUG] Entering RetroSynthesisTool._run {'#'*30}")
        print(f"[DEBUG][RetroSynthesisTool] Received compound_name: '{compound_name}'")

        if not isinstance(compound_name, str) or not compound_name.strip():
             print("[ERROR][RetroSynthesisTool] Invalid input: compound_name must be a non-empty string.")
             return {
                 "status": "error",
                 "message": "Invalid input: compound_name must be a non-empty string.",
                 "traceback": None
             }

        try:
            # Call the main function directly with parameters
            print(f"[DEBUG][RetroSynthesisTool] Calling run_with_params for '{compound_name}'...")
            raw_output = run_with_params(
                material=compound_name,
                num_results=10,    # Default values, consider making configurable if needed
                alignment=True,
                expansion=True,
                filtration=False
            )
            run_duration = time.time() - tool_start_time
            print(f"[DEBUG][RetroSynthesisTool] run_with_params finished in {run_duration:.4f} seconds.")

            if not raw_output or raw_output.startswith("Error:"):
                 print(f"[ERROR][RetroSynthesisTool] run_with_params returned an error state or empty output: {raw_output}")
                 # Provide a more informative message if possible
                 error_detail = raw_output if raw_output else "No output received"
                 return {
                    "status": "error",
                    "message": f"Retrosynthesis pipeline failed or returned no usable output. Details: {error_detail}",
                    "traceback": None # Traceback logged within run_with_params/retrosyn_main
                 }

            # Parse the output
            print("[DEBUG][RetroSynthesisTool] Parsing the raw output from run_with_params...")
            parse_start = time.time()
            parsed_data = parse_reaction_data(raw_output)
            parse_duration = time.time() - parse_start
            print(f"[DEBUG][RetroSynthesisTool] Parsing completed in {parse_duration:.4f} seconds.")

            if "parsing_error" in parsed_data:
                print(f"[WARN][RetroSynthesisTool] Parsing error occurred: {parsed_data['parsing_error']}. Returning potentially incomplete data.")

            # Format for better readability and add SMILES cleanup
            print("[DEBUG][RetroSynthesisTool] Processing parsed data and cleaning SMILES...")
            process_start = time.time()
            reactions = parsed_data.get("reactions", [])
            recommended_indices = parsed_data.get("recommended_pathway", [])
            reasoning = parsed_data.get("reasons", "")

            print(f"[DEBUG][RetroSynthesisTool] Found {len(reactions)} reactions, {len(recommended_indices)} recommended indices.")

            # Process each reaction to ensure proper SMILES
            for i, reaction in enumerate(reactions):
                # print(f"[DEBUG][RetroSynthesisTool] Processing reaction index {reaction.get('idx', i)}...") # Reduced verbosity
                original_smiles = reaction.get("reaction_smiles")
                cleaned_smiles = None
                generation_attempted = False

                if original_smiles:
                    # print(f"[DEBUG][RetroSynthesisTool]  - Original SMILES: '{original_smiles}'") # Reduced verbosity
                    cleaned_smiles = original_smiles.replace("[reactant_SMILES]", "").replace("[product_SMILES]", "")
                    cleaned_smiles = re.sub(r'\[[a-zA-Z_]+_SMILES\]', "", cleaned_smiles).strip()

                    if ">>" in cleaned_smiles and not re.search(r'\[.*?\]', cleaned_smiles):
                        # print(f"[DEBUG][RetroSynthesisTool]  - Cleaned valid SMILES: '{cleaned_smiles}'") # Reduced verbosity
                        pass
                    else:
                        # print(f"[DEBUG][RetroSynthesisTool]  - Cleaned SMILES '{cleaned_smiles}' seems invalid or incomplete. Attempting generation.") # Reduced verbosity
                        cleaned_smiles = None
                else:
                    # print("[DEBUG][RetroSynthesisTool]  - No original reaction SMILES found. Attempting generation.") # Reduced verbosity
                    pass


                if not cleaned_smiles:
                    generation_attempted = True
                    # print("[DEBUG][RetroSynthesisTool]  - Generating SMILES from reactants/products...") # Reduced verbosity
                    reactant_smiles_list = []
                    for reactant in reaction.get("reactants", []):
                        # print(f"[DEBUG][RetroSynthesisTool]    - Getting SMILES for reactant: '{reactant}'") # Reduced verbosity
                        smiles = get_smiles_for_compound(reactant)
                        if smiles:
                            reactant_smiles_list.append(smiles)
                        else:
                            # print(f"[WARN][RetroSynthesisTool]    - Could not get SMILES for reactant '{reactant}'. Using name.") # Reduced verbosity
                            reactant_smiles_list.append(reactant)

                    product_smiles_list = []
                    for product in reaction.get("products", []):
                        # print(f"[DEBUG][RetroSynthesisTool]    - Getting SMILES for product: '{product}'") # Reduced verbosity
                        smiles = get_smiles_for_compound(product)
                        if smiles:
                            product_smiles_list.append(smiles)
                        else:
                            # print(f"[WARN][RetroSynthesisTool]    - Could not get SMILES for product '{product}'. Using name.") # Reduced verbosity
                            product_smiles_list.append(product)

                    reactants_str = '.'.join(reactant_smiles_list)
                    products_str = '.'.join(product_smiles_list)
                    generated_smiles = f"{reactants_str}>>{products_str}"

                    if ">>" in generated_smiles and reactants_str and products_str:
                         # print(f"[DEBUG][RetroSynthesisTool]  - Generated SMILES: '{generated_smiles}'") # Reduced verbosity
                         cleaned_smiles = generated_smiles
                    else:
                         # print(f"[WARN][RetroSynthesisTool]  - Generated SMILES '{generated_smiles}' seems invalid (missing reactants/products?). Keeping as None.") # Reduced verbosity
                         cleaned_smiles = None

                reaction["cleaned_reaction_smiles"] = cleaned_smiles
                reaction["smiles_generation_attempted"] = generation_attempted

            process_duration = time.time() - process_start
            print(f"[DEBUG][RetroSynthesisTool] Data processing and SMILES cleaning completed in {process_duration:.4f} seconds.")

            final_result = {
                "status": "success",
                "data": {
                    "target_compound": compound_name,
                    "reactions": reactions,
                    "reasoning": reasoning,
                    "recommended_indices": recommended_indices,
                    "parsing_info": "Parsing successful" if "parsing_error" not in parsed_data else f"Parsing warning: {parsed_data['parsing_error']}"
                }
            }
            total_tool_duration = time.time() - tool_start_time
            print(f"[DEBUG][RetroSynthesisTool] Successfully processed '{compound_name}'. Total tool time: {total_tool_duration:.4f} seconds.")
            print(f"{'#'*30} [DEBUG] Exiting RetroSynthesisTool._run {'#'*30}\n")
            return final_result

        except Exception as e:
            total_tool_duration = time.time() - tool_start_time
            error_msg = f"Unhandled error in RetroSynthesisTool._run after {total_tool_duration:.4f} seconds: {str(e)}\n{traceback.format_exc()}"
            print(f"[ERROR][RetroSynthesisTool] {error_msg}")
            print(f"{'#'*30} [DEBUG] Exiting RetroSynthesisTool._run (Error) {'#'*30}\n")
            # Include traceback in the returned error message for better debugging from the caller
            return {
                "status": "error",
                "message": f"Unhandled error in RetroSynthesisTool: {str(e)}",
                "traceback": traceback.format_exc()
            }

# Function for direct calling remains the same...
def run_retrosynthesis(compound_name):
    """
    Run the RetroSynthesis Agent for a given compound and process the results.
    Uses the RetroSynthesisTool internally.
    """
    print(f"\n[DEBUG][run_retrosynthesis] Direct call initiated for compound: '{compound_name}'")
    t_start = time.time()
    # Create an instance of the tool
    try:
        tool = RetroSynthesisTool()
        print("[DEBUG][run_retrosynthesis] RetroSynthesisTool instantiated.")
        # Call the tool's _run method
        result = tool._run(compound_name)
        duration = time.time() - t_start
        print(f"[DEBUG][run_retrosynthesis] RetroSynthesisTool._run completed in {duration:.4f} seconds.")
        return result
    except Exception as e:
        duration = time.time() - t_start
        print(f"[ERROR][run_retrosynthesis] Exception during tool instantiation or execution after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
        return {
            "status": "error",
            "message": f"Failed to run retrosynthesis: {str(e)}",
            "traceback": traceback.format_exc()
        }

# Script Finish and Example Usage remains the same...
script_duration = time.time() - start_time
print(f"[DEBUG] retrosynthesis.py script initialization and definition complete. Total load time: {script_duration:.4f} seconds.")

# # Example Usage (optional, for testing)
# if __name__ == "__main__":
#     print("\n[INFO] Running example usage directly in retrosynthesis.py...")
#     # test_compound = "aspirin"
#     # test_compound = "paracetamol"
#     test_compound = "flubendiamide" # Use the failing example
#     print(f"[INFO] Target compound: {test_compound}")

#     # Use the direct call function
#     results = run_retrosynthesis(test_compound)

#     print("\n[INFO] Example Usage Result:")
#     # Use ensure_ascii=False for potentially non-ASCII chars in chemical names/reasoning
#     print(json.dumps(results, indent=2, ensure_ascii=False))

#     if results.get("status") == "success":
#         print("\n[INFO] Retrosynthesis ran successfully.")
#         if results.get("data"):
#              print(f"  - Found {len(results['data'].get('reactions', []))} reactions.")
#              print(f"  - Recommended indices: {results['data'].get('recommended_indices', [])}")
#              print(f"  - Reasoning provided: {'Yes' if results['data'].get('reasoning') else 'No'}")
#     else:
#         print("\n[INFO] Retrosynthesis encountered an error.")
#         print(f"  - Message: {results.get('message')}")
#         # Optionally print traceback if included
#         # if results.get("traceback"):
#         #    print(f"  - Traceback:\n{results.get('traceback')}")

#     print("[INFO] Example usage finished.")