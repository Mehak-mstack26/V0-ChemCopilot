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
    # print(f"[DEBUG] current_dir: {current_dir}")
    # print(f"[DEBUG] project_root: {project_root}")
    # print(f"[DEBUG] retrosyn_agent_dir: {retrosyn_agent_dir}")

    # --- [DEBUG] Adding Paths to sys.path ---
    # print("[DEBUG] Adding directories to sys.path...")
    sys.path.append(project_root)
    sys.path.append(os.path.join(project_root, 'RetroSynthesisAgent'))
    sys.path.append(retrosyn_agent_dir)
    # print(f"[DEBUG] sys.path updated: {sys.path}")

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
    # from RetroSynthesisAgent.RetroSynAgent.knowledgeGraph import KnowledgeGraph # Assuming not directly used here based on provided main.py
    # print("[DEBUG] Imported KnowledgeGraph")
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
    # print("[DEBUG] Original CommonSubstanceDB.read_data_from_json saved.")

    # Define a patched version of the method
    def patched_read_data_from_json(self, filename):
        """Patched method to read data from JSON with correct path resolution and debugging"""
        t_start = time.time()
        # print(f"[DEBUG][patched_read_data_from_json] Called with filename: {filename}") # Reduced verbosity
        # Map common filenames to their absolute paths
        file_map = {
            'RetroSynAgent/emol.json': os.path.join(retrosyn_agent_dir, 'emol.json'),
            'RetroSynAgent/common_chemicals.json': os.path.join(retrosyn_agent_dir, 'common_chemicals.json'),
            'emol.json': os.path.join(retrosyn_agent_dir, 'emol.json'),
            'common_chemicals.json': os.path.join(retrosyn_agent_dir, 'common_chemicals.json'),
            # Add other JSON files as needed
        }

        actual_path = file_map.get(filename, filename)
        if not os.path.isabs(actual_path):
             maybe_path = os.path.join(retrosyn_agent_dir, actual_path)
             if os.path.exists(maybe_path):
                 actual_path = maybe_path
             else: 
                 maybe_path_root = os.path.join(project_root, actual_path)
                 if os.path.exists(maybe_path_root):
                     actual_path = maybe_path_root
        # print(f"[DEBUG][patched_read_data_from_json] Resolved path: {actual_path}") # Reduced verbosity

        try:
            # print(f"[DEBUG][patched_read_data_from_json] Attempting to read JSON from: {actual_path}...") # Reduced verbosity
            io_start = time.time()
            with open(actual_path, 'r', encoding='utf-8') as file:
                data = json.load(file)
            io_duration = time.time() - io_start
            # print(f"[DEBUG][patched_read_data_from_json] Successfully read {actual_path} in {io_duration:.4f} seconds.") # Reduced verbosity
            # print(f"[DEBUG][patched_read_data_from_json] Finished in {time.time() - t_start:.4f} seconds.") # Reduced verbosity
            return data
        except FileNotFoundError:
            print(f"[WARN][patched_read_data_from_json] File not found at primary path: {actual_path}")
            alt_paths = [
                os.path.join(retrosyn_agent_dir, os.path.basename(filename)),
            ]
            for alt_path in alt_paths:
                # print(f"[DEBUG][patched_read_data_from_json] Trying alternate path: {alt_path}") # Reduced verbosity
                if os.path.exists(alt_path):
                    try:
                        # print(f"[DEBUG][patched_read_data_from_json] Attempting to read JSON from alternate: {alt_path}...") # Reduced verbosity
                        io_start = time.time()
                        with open(alt_path, 'r', encoding='utf-8') as file:
                            data = json.load(file)
                        io_duration = time.time() - io_start
                        # print(f"[DEBUG][patched_read_data_from_json] Successfully loaded from alternate: {alt_path} in {io_duration:.4f} seconds.") # Reduced verbosity
                        # print(f"[DEBUG][patched_read_data_from_json] Finished in {time.time() - t_start:.4f} seconds.") # Reduced verbosity
                        return data
                    except Exception as e_alt:
                        print(f"[ERROR][patched_read_data_from_json] Error reading alternate path {alt_path}: {str(e_alt)}")
                        continue
                # else:
                     # print(f"[DEBUG][patched_read_data_from_json] Alternate path {alt_path} does not exist.") # Reduced verbosity
            print(f"[ERROR][patched_read_data_from_json] Could not find or read file {filename} using resolved path {actual_path} or alternates.")
            # print(f"[DEBUG][patched_read_data_from_json] Finished (failed) in {time.time() - t_start:.4f} seconds. Returning empty dict.") # Reduced verbosity
            return {}
        except Exception as e:
            print(f"[ERROR][patched_read_data_from_json] Error reading {actual_path}: {str(e)}\n{traceback.format_exc()}")
            # print(f"[DEBUG][patched_read_data_from_json] Finished (error) in {time.time() - t_start:.4f} seconds. Returning empty dict.") # Reduced verbosity
            return {}

    CommonSubstanceDB.read_data_from_json = patched_read_data_from_json
    print("[DEBUG] CommonSubstanceDB.read_data_from_json monkey patched successfully.")

except Exception as e:
    print(f"[ERROR] Failed during monkey patching: {e}\n{traceback.format_exc()}")

# --- [DEBUG] Initializing NameToSMILES Tool ---
print("[DEBUG] Initializing NameToSMILES tool...")
try:
    name_to_smiles_tool = NameToSMILES()
    print("[DEBUG] NameToSMILES tool initialized successfully.")
except Exception as e:
    print(f"[ERROR] Failed to initialize NameToSMILES tool: {e}\n{traceback.format_exc()}")
    name_to_smiles_tool = None 

def get_smiles_for_compound(compound_name, timeout=15):
    """Get SMILES notation for a compound name using the NameToSMILES tool with timeout"""
    t_start = time.time()
    # print(f"[DEBUG][get_smiles_for_compound] Called for: '{compound_name}' with timeout={timeout}s") # Reduced verbosity
    if not name_to_smiles_tool:
        print("[ERROR][get_smiles_for_compound] NameToSMILES tool not initialized.")
        return None
    try:
        # print(f"[DEBUG][get_smiles_for_compound] Calling name_to_smiles_tool._run('{compound_name}')...") # Reduced verbosity
        result = name_to_smiles_tool._run(compound_name) 
        duration = time.time() - t_start
        # print(f"[DEBUG][get_smiles_for_compound] Received result: '{result}' in {duration:.4f} seconds.") # Reduced verbosity

        if result and "SMILES:" in result:
            smiles = result.split("SMILES:")[1].split("\n")[0].strip()
            # print(f"[DEBUG][get_smiles_for_compound] Extracted SMILES: '{smiles}'. Total time: {time.time() - t_start:.4f}s") # Reduced verbosity
            return smiles
        else:
            # print(f"[WARN][get_smiles_for_compound] 'SMILES:' not found in result for '{compound_name}'. Total time: {time.time() - t_start:.4f}s") # Reduced verbosity
            return None
    except Exception as e:
        duration = time.time() - t_start
        print(f"[ERROR][get_smiles_for_compound] Exception during SMILES lookup for '{compound_name}' after {duration:.4f}s: {e}\n{traceback.format_exc()}")
        return None

def parse_reaction_data(raw_text):
    """Parse the reaction data from the raw text output"""
    t_start = time.time()
    print("[DEBUG][parse_reaction_data] Starting parsing of raw text...")

    reactions_from_llm = [] 
    recommended_indices_for_rank1_path = [] # Specifically for Rank 1 pathway's reaction indices
    overall_reasoning = "" 
    ranked_pathways_parsed_details = [] # To store structured ranked pathway data

    try:
        # 1. Extract overall "Recommended Reaction Pathway" line (less common with Top-3 format but good to keep)
        rec_match_overall = re.search(r"Recommended Reaction Pathway:\s*([^\n]+)", raw_text)
        if rec_match_overall:
            pathway_text_overall = rec_match_overall.group(1)
            extracted_indices_overall = re.findall(r'idx\d+', pathway_text_overall)
            if extracted_indices_overall:
                # This would typically be the primary source if the LLM uses this line for the #1 pathway.
                # However, with Top-3 format, Rank 1 details are more reliable.
                # We will prioritize Rank 1 parsing later.
                pass # print(f"[DEBUG][parse_reaction_data] Found overall recommended pathway text: {pathway_text_overall}")
        
        # 2. Extract overall "Reasons for Recommendation" or general "Reasons"
        # This usually refers to the reasoning for the Rank 1 choice or the overall selection process.
        reasons_match = re.search(r"Reasons for Recommendation:\s*((?:.|\n)*?)(?=\n---|\nRank|\nAnalysis:Overall Comparison and Justification:|$)", raw_text, re.DOTALL)
        if not reasons_match: # Fallback to simpler "Reasons:"
            reasons_match = re.search(r"Reasons:\s*((?:.|\n)*?)(?=\n---|\nRank|\nAnalysis:Overall Comparison and Justification:|$)", raw_text, re.DOTALL)
        if reasons_match:
            overall_reasoning = reasons_match.group(1).strip()
            # print(f"[DEBUG][parse_reaction_data] Found Overall Reasons block (length: {len(overall_reasoning)}).")
        # else:
            # print("[WARN][parse_reaction_data] Overall Reasons section not found or format unrecognized.")

        # 3. Extract Ranked Pathways (Rank X: Pathway YYY) and their specific reasons
        # This regex also captures the "Details:" block text for each rank.
        ranked_pathway_blocks_iter = re.finditer(
            r"Rank\s+(?P<rank>\d+):\s*Pathway\s+(?P<pathway_indices_str>[^\n]+)\s*\n"
            r"Details:\s*(?P<details_text>(?:.|\n)*?)"
            r"(?=Rank\s+\d+:|Reasons for Rank\s+\d+:|Overall Comparison and Justification:|$)",
            raw_text, re.DOTALL
        )
        
        temp_ranked_pathways = []
        for match in ranked_pathway_blocks_iter:
            rank_num = int(match.group("rank"))
            pathway_indices_str_from_rank_line = match.group("pathway_indices_str").strip() # e.g., "31" or "[Not available]"
            details_block_text = match.group("details_text").strip()
            
            # Extract specific reasons for this rank
            rank_reasons_pattern = r"Reasons for Rank " + str(rank_num) + r":\s*((?:.|\n)*?)(?=Rank\s+\d+:|Overall Comparison and Justification:|$)"
            rank_reasons_match = re.search(rank_reasons_pattern, raw_text, re.DOTALL)
            specific_rank_reasons = rank_reasons_match.group(1).strip() if rank_reasons_match else ""
            
            temp_ranked_pathways.append({
                "rank": rank_num,
                "pathway_indices_str": pathway_indices_str_from_rank_line, 
                "details_text": details_block_text, # The block of text under "Details:"
                "reasons": specific_rank_reasons
            })
        
        ranked_pathways_parsed_details = sorted(temp_ranked_pathways, key=lambda x: x["rank"])
        # print(f"[DEBUG][parse_reaction_data] Parsed {len(ranked_pathways_parsed_details)} ranked pathway blocks with their details text.")

        # 4. Populate 'recommended_indices_for_rank1_path' from Rank 1 pathway's "pathway_indices_str"
        # And if overall_reasoning is still empty, use Rank 1 reasons.
        rank1_pathway_data = next((pw for pw in ranked_pathways_parsed_details if pw["rank"] == 1), None)
        if rank1_pathway_data:
            rank1_indices_str_from_line = rank1_pathway_data.get("pathway_indices_str", "")
            # The string might be just "31" or "idx31" or "Hypothetical..."
            # We are interested in actual idx values for the reaction list.
            # The reaction objects themselves will be parsed from the "details_text"
            # This `recommended_indices_for_rank1_path` will guide which reactions from `all_parsed_reactions`
            # belong to the Rank 1 pathway.
            
            # If pathway_indices_str is like "31", convert to "idx31"
            # If it's already "idx31", re.findall will get it.
            # If it's "Hypothetical" or "Not available", re.findall will be empty.
            potential_indices = re.findall(r'\d+', rank1_indices_str_from_line) # Get numbers like '31'
            if potential_indices and not "Not available" in rank1_indices_str_from_line and not "[Hypothetical" in rank1_indices_str_from_line:
                 # Prefer 'idx' prefixed versions if LLM provides them in the details,
                 # otherwise, construct from the pathway line.
                # For now, let's extract from the details_text of rank 1.
                rank1_details_text = rank1_pathway_data.get("details_text", "")
                reaction_idx_in_rank1_details = re.findall(r"Reaction idx:\s*(idx\d+|\d+)", rank1_details_text)
                
                if reaction_idx_in_rank1_details:
                    # Ensure they are 'idx' prefixed
                    recommended_indices_for_rank1_path = [
                        f"idx{val}" if val.isdigit() else val 
                        for val in reaction_idx_in_rank1_details
                    ]
                elif not recommended_indices_for_rank1_path and potential_indices: # Fallback to indices from "Pathway XXX" line
                    recommended_indices_for_rank1_path = [f"idx{pi}" for pi in potential_indices]

                # print(f"[DEBUG][parse_reaction_data] Populated 'recommended_indices_for_rank1_path' from Rank 1: {recommended_indices_for_rank1_path}")

            if not overall_reasoning and rank1_pathway_data.get("reasons"): # If main reasoning is empty, use Rank 1's
                overall_reasoning = rank1_pathway_data.get("reasons")
        
        # 5. Parse ALL "Reaction idx:" blocks from the *entire raw_text*
        # This captures reactions from "Analysis:" section AND from "Details:" sections of ranked pathways.
        # Deduplication will occur later.
        reaction_blocks_iter = re.finditer(
            r"Reaction idx:\s*(?P<idx_val>(?:idx\d+|\d+|Not available))\s*\n" # Capture "idxDDD", "DDD", or "Not available"
            r"(?:Reactants:\s*(?P<reactants>.+?)\s*\n)?"
            r"(?:Products:\s*(?P<products>.+?)\s*\n)?"
            r"(?:Reaction SMILES:\s*(?P<smiles>.+?)\s*\n)?"
            r"(?:Conditions:\s*(?P<conditions>.+?)\s*\n)?"
            r"(?:Source:\s*(?P<source>.+?)\s*\n)?"
            r"(?:SourceLink:\s*\[?(?P<link>.+?)\]?(?:\s|\n|$))?",
            raw_text, re.MULTILINE | re.DOTALL # DOTALL for multi-line fields
        )
        reaction_count = 0
        for match in reaction_blocks_iter:
            reaction_count += 1
            reaction_dict = match.groupdict()
            
            idx_val_raw = reaction_dict.get("idx_val").strip()
            # Normalize idx to always be "idx<number>" or the original string if not just numbers
            current_idx_str = f"idx{idx_val_raw}" if idx_val_raw.isdigit() else idx_val_raw

            current_reaction = {
                "idx": current_idx_str,
                "reactants": [r.strip() for r in reaction_dict.get("reactants").split(",")] if reaction_dict.get("reactants") else [],
                "products": [p.strip() for p in reaction_dict.get("products").split(",")] if reaction_dict.get("products") else [],
                "reaction_smiles": reaction_dict.get("smiles").strip() if reaction_dict.get("smiles") else None,
                "conditions": {},
                "source": reaction_dict.get("source").strip() if reaction_dict.get("source") else None,
                "source_link": reaction_dict.get("link").strip() if reaction_dict.get("link") else None
            }
            if reaction_dict.get("conditions"):
                conditions_str = reaction_dict.get("conditions").strip()
                for part in conditions_str.split(","):
                    part = part.strip()
                    if ":" in part:
                        try:
                            key, val = part.split(":", 1)
                            current_reaction["conditions"][key.strip().lower()] = val.strip()
                        except ValueError:
                            current_reaction["conditions"][f"condition_{len(current_reaction['conditions'])}"] = part
                    elif part: 
                        current_reaction["conditions"][part.lower()] = True 
            reactions_from_llm.append(current_reaction)

        # 6. Deduplicate reactions by 'idx', keeping the last one encountered.
        unique_reactions_map = {}
        for r_obj in reactions_from_llm:
            idx_key = r_obj.get('idx')
            if idx_key and idx_key.lower() != "not available": 
                unique_reactions_map[idx_key] = r_obj 
        
        all_unique_parsed_reactions = list(unique_reactions_map.values())
        
        parsed_data_output = {
            "recommended_pathway_indices": recommended_indices_for_rank1_path,
            "reactions": all_unique_parsed_reactions, # All unique reactions found in the text
            "reasoning": overall_reasoning,
            "ranked_pathways_details": ranked_pathways_parsed_details 
        }
        duration = time.time() - t_start
        print(f"[DEBUG][parse_reaction_data] Finished parsing. Found {len(all_unique_parsed_reactions)} unique reactions. Recommended indices for Rank 1: {recommended_indices_for_rank1_path}. Total time: {duration:.4f} seconds.")
        return parsed_data_output

    except Exception as e:
        duration = time.time() - t_start
        print(f"[ERROR][parse_reaction_data] Exception during parsing after {duration:.4f} seconds: {e}\n{traceback.format_exc()}")
        return {
            "recommended_pathway_indices": recommended_indices_for_rank1_path,
            "reactions": reactions_from_llm, 
            "reasoning": overall_reasoning,
            "ranked_pathways_details": ranked_pathways_parsed_details,
            "parsing_error": str(e)
        }
    
# ... (retrosyn_main function definition starts here)
def retrosyn_main(material=None, num_results=10, alignment=True, expansion=True, filtration=False, recommend_type="general"):
    main_start_time = time.time()
    print(f"\n{'='*30} [DEBUG] Entering retrosyn_main {'='*30}")
    print(f"[DEBUG][retrosyn_main] Parameters: material='{material}', num_results={num_results}, alignment={alignment}, expansion={expansion}, filtration={filtration}, recommend_type='{recommend_type}'")

    if not material:
        print("[ERROR][retrosyn_main] 'material' parameter cannot be None.")
        raise ValueError("'material' parameter is required.")

    # --- Folder Setup ---
    step_start_time = time.time()
    # print("[DEBUG][retrosyn_main] Setting up folder paths...")
    try:
        material_safe_name = material.lower().replace(" ", "_").replace("/", "-").replace("\\", "-")
        material_safe_name = re.sub(r'[^\w\-]+', '', material_safe_name)[:50]

        base_folder = os.path.join(project_root, 'RetroSynthesisAgent', f"results_{material_safe_name}")
        pdf_folder_name = os.path.join(base_folder, 'pdfs')
        result_folder_name = os.path.join(base_folder, 'llm_outputs')
        tree_folder_name = os.path.join(base_folder, 'trees')
        
        result_json_name = 'pdf_extraction_results.json' # Using .json extension
        result_json_path = os.path.join(result_folder_name, result_json_name)


        # print(f"[DEBUG][retrosyn_main] PDF folder: {pdf_folder_name}")
        # print(f"[DEBUG][retrosyn_main] Result folder: {result_folder_name}")
        # print(f"[DEBUG][retrosyn_main] Tree folder: {tree_folder_name}")
        # print(f"[DEBUG][retrosyn_main] Result JSON path: {result_json_path}")


        # print("[DEBUG][retrosyn_main] Creating directories if they don't exist...")
        os.makedirs(pdf_folder_name, exist_ok=True)
        os.makedirs(result_folder_name, exist_ok=True)
        os.makedirs(tree_folder_name, exist_ok=True)
        # print("[DEBUG][retrosyn_main] Directories created/verified.")
        print(f"[DEBUG][retrosyn_main] Folder Setup completed in {time.time() - step_start_time:.4f} seconds.")
    except Exception as e:
        print(f"[ERROR][retrosyn_main] Failed during folder setup: {e}\n{traceback.format_exc()}")
        raise

    # --- Object Initialization ---
    step_start_time = time.time()
    print("[DEBUG][retrosyn_main] Initializing core components...")
    pdf_processor = None
    try:
        entityalignment = EntityAlignment()
        # print("[DEBUG][retrosyn_main] Initialized EntityAlignment.")
        treeloader = TreeLoader()
        # print("[DEBUG][retrosyn_main] Initialized TreeLoader.")
        tree_expansion = TreeExpansion()
        # print("[DEBUG][retrosyn_main] Initialized TreeExpansion.")
        reactions_filtration = ReactionsFiltration()
        # print("[DEBUG][retrosyn_main] Initialized ReactionsFiltration.")
        
        pdf_processor = PDFProcessor(pdf_folder_name=pdf_folder_name, result_folder_name=result_folder_name,
                                     result_json_name=result_json_name) # pass only filename for result_json_name
        # print("[DEBUG][retrosyn_main] Initialized PDFProcessor.")
        print(f"[DEBUG][retrosyn_main] Object Initialization completed in {time.time() - step_start_time:.4f} seconds.")
    except Exception as e:
        print(f"[ERROR][retrosyn_main] Failed during object initialization: {e}\n{traceback.format_exc()}")
        raise

    # === Stage 1: Extract Info ===
    # ... (Stages 1-7 remain largely the same, ensure debug prints are managed if too verbose)
    # (Copied from thought process, ensure any minor debug print changes from original are fine)
    print(f"\n--- [DEBUG][retrosyn_main] Stage 1: Extract Info for '{material}' ---")
    step_start_time = time.time()
    # 1.1 Query Literatures & Download
    # print(f"[DEBUG][retrosyn_main] 1.1 Initializing PDFDownloader for '{material}'...")
    pdf_name_list = []
    try:
        downloader = PDFDownloader(material, pdf_folder_name=pdf_folder_name, num_results=num_results, n_thread=min(num_results, 5))
        # print(f"[DEBUG][retrosyn_main] Calling PDFDownloader.main() with num_results={num_results}...")
        pdf_name_list = downloader.main() 
        duration = time.time() - step_start_time
        # print(f"[DEBUG][retrosyn_main] PDFDownloader.main() finished in {duration:.4f} seconds.")
        print(f'[DEBUG][retrosyn_main] Attempted download for {len(pdf_name_list)} potential PDFs.')
    except Exception as e:
        duration = time.time() - step_start_time
        print(f"[ERROR][retrosyn_main] Exception during PDF download/retrieval after {duration:.4f} seconds: {e}") # Simplified traceback
        pdf_name_list = [] 
    actual_pdfs = []
    try:
        actual_pdfs = [f for f in os.listdir(pdf_folder_name) if f.lower().endswith('.pdf')]
        if not actual_pdfs: print(f"Warning: No PDF files found in '{pdf_folder_name}' after download attempt.")
        else: print(f"Found {len(actual_pdfs)} PDF files in '{pdf_folder_name}'.")
    except FileNotFoundError: print(f"Error: PDF directory '{pdf_folder_name}' not found or inaccessible.")
    # 1.2 Extract Info from PDFs
    step_start_time = time.time()
    # print(f"[DEBUG][retrosyn_main] 1.2 Processing PDFs using PDFProcessor...")
    pdf_processing_successful = False
    if not os.path.isdir(pdf_folder_name) or not actual_pdfs: print(f"Skipping PDF processing: Directory '{pdf_folder_name}' does not exist or no PDFs found.")
    elif not pdf_processor: print(f"Skipping PDF processing: PDFProcessor not initialized.")
    else:
        try:
            pdf_processor.load_existing_results()
            if not hasattr(pdf_processor, 'results'): pdf_processor.results = {} 
            # print(f"[DEBUG][retrosyn_main] Found {len(pdf_processor.results)} existing results initially.")
            pdf_processor.process_pdfs_txt(save_batch_size=2) 
            duration = time.time() - step_start_time
            # print(f"[DEBUG][retrosyn_main] PDFProcessor finished in {duration:.4f} seconds.")
            # print(f"[DEBUG][retrosyn_main] Total results after processing: {len(pdf_processor.results)}")
            print(f"PDF processing complete. Extracted reactions data stored in: {result_json_path}")
            pdf_processing_successful = True 
        except Exception as e:
            duration = time.time() - step_start_time
            print(f"[ERROR][retrosyn_main] Exception during PDF processing stage after {duration:.4f} seconds: {e}")
            if not hasattr(pdf_processor, 'results'): pdf_processor.results = {}
    print(f"[DEBUG][retrosyn_main] Stage 1 completed in {time.time() - step_start_time:.4f} seconds overall for PDF download & process.")

    # === Stage 2 & 3: Align Root Node & Initial Tree Construction ===
    # ... (This section remains as in the thought process, detailed logic for tree loading/building/alignment)
    print(f"\n--- [DEBUG][retrosyn_main] Stage 2 & 3: Align Root Node & Initial Tree for '{material}' ---")
    step_start_time = time.time()
    current_tree = None 
    results_dict = {} 
    if not pdf_processing_successful and not os.path.exists(result_json_path):
         print("[WARN][retrosyn_main] Skipping Tree Build: PDF processing failed and no existing results JSON.")
    else:
        try:
            # print(f"[DEBUG][retrosyn_main] Aligning root node using EntityAlignment from '{result_json_path}'...")
            align_start_time_inner = time.time()
            results_dict = entityalignment.alignRootNode(result_folder_name, result_json_name, material)
            if not results_dict: print(f"Warning: No reactions found in '{result_json_path}' after root node alignment.")
            else: print(f"Root node alignment completed in {time.time() - align_start_time_inner:.4f}s. Loaded {len(results_dict)} sources.")
            tree_name_base = os.path.join(tree_folder_name, f'{material_safe_name}_base.pkl')
            # ... (rest of tree loading/building as in thought process) ...
            tree_loaded = False
            if os.path.exists(tree_name_base):
                try:
                    current_tree = treeloader.load_tree(tree_name_base)
                    print(f'Initial RetroSynthetic Tree loaded from: {tree_name_base}')
                    tree_loaded = True
                except Exception as e_load: print(f"Warning: Could not load base tree {tree_name_base}: {e_load}. Rebuilding.")
            if not tree_loaded and results_dict: 
                current_tree = Tree(material.lower(), result_dict=results_dict)
                # print('Constructing initial RetroSynthetic Tree...')
                current_tree.construct_tree()
                treeloader.save_tree(current_tree, tree_name_base)
                print(f"Initial tree constructed and saved to: {tree_name_base}")
            elif not tree_loaded: print("Skipping initial tree construction: no reaction data for root alignment.")
            # print(f'Initial tree: {countNodes(current_tree)} nodes, {len(searchPathways(current_tree))} pathways.')
            if alignment and current_tree:
                # print('[DEBUG][retrosyn_main] Performing Node Alignment (Initial Tree)...')
                tree_name_base_alg = os.path.join(tree_folder_name, f'{material_safe_name}_base_alg.pkl')
                # ... (alignment logic as in thought process) ...
                tree_loaded_aligned = False
                if os.path.exists(tree_name_base_alg):
                    try:
                        current_tree = treeloader.load_tree(tree_name_base_alg)
                        print(f'Aligned initial RetroSynthetic Tree loaded from: {tree_name_base_alg}')
                        tree_loaded_aligned = True
                    except Exception as e_load_alg: print(f"Warning: Could not load aligned base tree: {e_load_alg}. Re-aligning.")
                if not tree_loaded_aligned:
                    reactions_base = current_tree.reactions
                    if reactions_base:
                        reactions_base_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_base)
                        reactions_base_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_base_alg_1)
                        aligned_base_tree = Tree(material.lower(), reactions=reactions_base_alg_all)
                        aligned_base_tree.construct_tree()
                        treeloader.save_tree(aligned_base_tree, tree_name_base_alg)
                        print(f"Aligned initial tree saved to: {tree_name_base_alg}")
                        current_tree = aligned_base_tree 
                    else: print("Warning: No reactions in current tree for alignment.")
            print(f'Tree state after initial processing: {countNodes(current_tree)} nodes, {len(searchPathways(current_tree))} pathways.')
            print(f"[DEBUG][retrosyn_main] Stage 2 & 3 completed in {time.time() - step_start_time:.4f} seconds.")
        except Exception as e: print(f"[ERROR][retrosyn_main] Exception during Initial Tree Build/Alignment: {e}")

    # === Stage 4 & 5: Tree Expansion & Optional Alignment ===
    # ... (This section remains as in the thought process)
    print(f"\n--- [DEBUG][retrosyn_main] Stage 4 & 5: Tree Expansion & Alignment for '{material}' ---")
    step_start_time = time.time()
    expanded_tree_processed = False
    if expansion:
        # print(f"[DEBUG][retrosyn_main] Performing Tree Expansion...")
        tree_name_exp = os.path.join(tree_folder_name, f'{material_safe_name}_expanded.pkl')
        # ... (expansion logic as in thought process) ...
        tree_loaded_expanded = False
        if os.path.exists(tree_name_exp):
            try:
                current_tree = treeloader.load_tree(tree_name_exp)
                print(f'Loaded previously expanded RetroSynthetic Tree from: {tree_name_exp}')
                tree_loaded_expanded = True
                expanded_tree_processed = True
            except Exception as e_load_exp: print(f"Warning: Could not load expanded tree: {e_load_exp}. Re-expanding.")
        if not tree_loaded_expanded:
            if not results_dict: print("Cannot perform expansion: initial results dictionary is empty.")
            else:
                try:
                     results_dict_additional = tree_expansion.treeExpansion(result_folder_name, result_json_name, results_dict, material, expansion=True, max_iter=5)
                     if results_dict_additional:
                         updated_results_dict = tree_expansion.update_dict(results_dict, results_dict_additional)
                         expanded_tree = Tree(material.lower(), result_dict=updated_results_dict)
                         # print('Constructing Expanded RetroSynthetic Tree...')
                         expanded_tree.construct_tree()
                         treeloader.save_tree(expanded_tree, tree_name_exp)
                         print(f"Expanded tree constructed and saved to: {tree_name_exp}")
                         current_tree = expanded_tree 
                         results_dict = updated_results_dict
                     else: print("Expansion did not yield additional reactions.")
                     expanded_tree_processed = True
                except Exception as e_exp: print(f"Error during tree expansion: {e_exp}")
        if expanded_tree_processed: print(f'Tree state after expansion: {countNodes(current_tree)} nodes, {len(searchPathways(current_tree))} pathways.')
        if alignment and expanded_tree_processed and current_tree:
            # print('[DEBUG][retrosyn_main] Performing Node Alignment (Expanded Tree)...')
            tree_name_exp_alg = os.path.join(tree_folder_name, f'{material_safe_name}_expanded_alg.pkl')
            # ... (alignment logic for expanded tree) ...
            tree_loaded_exp_aligned = False
            if os.path.exists(tree_name_exp_alg):
                try:
                    current_tree = treeloader.load_tree(tree_name_exp_alg)
                    print(f'Aligned expanded RetroSynthetic Tree loaded from: {tree_name_exp_alg}')
                    tree_loaded_exp_aligned = True
                except Exception as e_load_exp_alg: print(f"Warning: Could not load aligned expanded tree: {e_load_exp_alg}. Re-aligning.")
            if not tree_loaded_exp_aligned:
                reactions_exp = current_tree.reactions
                if reactions_exp:
                    reactions_exp_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_exp)
                    reactions_exp_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_exp_alg_1)
                    aligned_exp_tree = Tree(material.lower(), reactions=reactions_exp_alg_all)
                    aligned_exp_tree.construct_tree()
                    treeloader.save_tree(aligned_exp_tree, tree_name_exp_alg)
                    print(f"Aligned expanded tree saved to: {tree_name_exp_alg}")
                    current_tree = aligned_exp_tree 
                else: print("Warning: No reactions in expanded tree for alignment.")
            print(f'Tree state after post-expansion alignment: {countNodes(current_tree)} nodes, {len(searchPathways(current_tree))} pathways.')
    else: print("[DEBUG][retrosyn_main] Expansion skipped as expansion=False.")
    print(f"[DEBUG][retrosyn_main] Stage 4 & 5 (Expansion & Alignment) completed in {time.time() - step_start_time:.4f} seconds.")

    # === Stage 6: Filtration ===
    # ... (This section remains as in the thought process)
    final_tree = current_tree
    if filtration and final_tree:
        print(f"\n--- [DEBUG][retrosyn_main] Stage 6: Performing Filtration ---")
        step_start_time = time.time()
        # ... (filtration logic as in thought process) ...
        tree_after_reaction_filter = None
        try:
            reactions_txt_filtered = reactions_filtration.filterReactions(final_tree)
            if not reactions_txt_filtered or not (isinstance(reactions_txt_filtered, str) and reactions_txt_filtered.strip().startswith("Reaction")) and not (isinstance(reactions_txt_filtered, list) and reactions_txt_filtered):
                print("Reaction filtration removed all reactions or returned invalid format.")
            else:
                # print("Constructing tree from filtered reactions...")
                tree_name_filtered_reactions = os.path.join(tree_folder_name, f'{material_safe_name}_filtered_reactions.pkl')
                try:
                    if isinstance(reactions_txt_filtered, str): tree_filtered_reactions_obj = Tree(material.lower(), reactions_txt=reactions_txt_filtered)
                    elif isinstance(reactions_txt_filtered, list): tree_filtered_reactions_obj = Tree(material.lower(), reactions=reactions_txt_filtered)
                    else: raise ValueError("Unsupported format from filterReactions")
                    tree_filtered_reactions_obj.construct_tree()
                    treeloader.save_tree(tree_filtered_reactions_obj, tree_name_filtered_reactions)
                    print(f"Tree based on filtered reactions saved to: {tree_name_filtered_reactions}")
                    tree_after_reaction_filter = tree_filtered_reactions_obj
                except Exception as e_filt_build: print(f"Error constructing tree from filtered reactions: {e_filt_build}.")
        except Exception as e_react_filt: print(f"Error during reaction filtration step: {e_react_filt}")
        tree_to_filter_pathways = tree_after_reaction_filter if tree_after_reaction_filter else final_tree
        if tree_to_filter_pathways:
            # print(f'Tree state before pathway filtration: {countNodes(tree_to_filter_pathways)} nodes, {len(searchPathways(tree_to_filter_pathways))} pathways.')
            # print("Filtering pathways based on validity...")
            try:
                reactions_filtration.filterPathways(tree_to_filter_pathways) 
                # print("Pathway filtration logic executed.")
                if tree_after_reaction_filter: final_tree = tree_after_reaction_filter
            except Exception as e_path_filt: print(f"Error during pathway filtration step: {e_path_filt}")
        # else: print("Skipping pathway filtration: no valid tree after reaction filtering.")
        print(f"[DEBUG][retrosyn_main] Stage 6 (Filtration) completed in {time.time() - step_start_time:.4f} seconds.")
    elif filtration and not final_tree: print("\nSkipping Filtration: no valid tree exists.")
    else: print("\n[DEBUG][retrosyn_main] Skipping Filtration as filtration=False.")

    # === Stage 7: Prepare Pathways for Recommendation ===
    # ... (This section remains as in the thought process)
    print(f"\n--- [DEBUG][retrosyn_main] Stage 7: Preparing Final Pathways for Recommendation ---")
    step_start_time = time.time()
    all_pathways_w_reactions_text = ""
    num_final_paths = 0
    if not final_tree:
         print("[ERROR][retrosyn_main] No valid final tree. Cannot extract pathways for recommendation.")
         raise RuntimeError("Failed to build or load synthesis tree.")
    try:
        all_pathways_w_reactions_text = reactions_filtration.getFullReactionPathways(final_tree)
        num_final_paths = len(searchPathways(final_tree))
        if not all_pathways_w_reactions_text.strip() and num_final_paths > 0: print(f"Warning: {num_final_paths} pathways in tree, but getFullReactionPathways returned empty.")
        elif not all_pathways_w_reactions_text.strip() and num_final_paths == 0: print("No pathways available in final tree.")
        else: print(f"Extracted text for {num_final_paths} pathways from final tree state.")
        print(f"[DEBUG][retrosyn_main] Stage 7 (Prepare Pathways) completed in {time.time() - step_start_time:.4f} seconds.")
    except Exception as e: print(f"Error generating full reaction pathways text: {e}")

    # === Stage 8: Recommendation ===
    # ... (This section remains as in the thought process)
    print(f"\n--- [DEBUG][retrosyn_main] Stage 8: Recommendation for '{material}' ({recommend_type}) ---")
    step_start_time = time.time()
    recommend_reactions_llm_output = "Error: Recommendation could not be generated."
    if not all_pathways_w_reactions_text.strip():
        print(f"\n[WARN][retrosyn_main] No valid reaction pathways for '{material}'.")
        recommend_reactions_llm_output = f"No valid reaction pathways found for {material} to recommend."
        output_filepath = os.path.join(result_folder_name, f'recommend_pathways_top3_{recommend_type}.txt')
        try:
            with open(output_filepath, 'w', encoding='utf-8') as f: f.write(recommend_reactions_llm_output)
        except Exception as e_save_err: print(f"Error saving 'no pathways' message: {e_save_err}")
    else:
        prompt_template = None
        if recommend_type == "general": prompt_template = prompts.recommend_prompt_template_general_top3
        elif recommend_type == "cost": prompt_template = prompts.recommend_prompt_template_cost_top3
        elif recommend_type == "condition": prompt_template = prompts.recommend_prompt_template_condition_top3
        if prompt_template:
            try:
                prompt_for_llm = prompt_template.format(substance=material, all_pathways=all_pathways_w_reactions_text)
                recommend_reactions_llm_output = recommendReactions(prompt_for_llm, result_folder_name, response_name=f'recommend_pathways_top3_{recommend_type}')
            except KeyError as e_fmt:
                 print(f"Error formatting prompt template: Missing placeholder {e_fmt}")
                 recommend_reactions_llm_output = f"Error: Prompt template formatting failed ({e_fmt})."
            except Exception as e_rec:
                 print(f"An unexpected error during recommendation: {e_rec}")
                 recommend_reactions_llm_output = f"Error during recommendation step: {e_rec}"
        else:
             print(f"Error: Invalid recommendation type '{recommend_type}'.")
             recommend_reactions_llm_output = f"Error: Invalid recommendation type selected."
    print(f"[DEBUG][retrosyn_main] Stage 8 (Recommendation) completed in {time.time() - step_start_time:.4f} seconds.")

    total_duration = time.time() - main_start_time
    print(f"\n{'='*30} [DEBUG] Exiting retrosyn_main {'='*30}")
    print(f"[DEBUG][retrosyn_main] Total execution time: {total_duration:.4f} seconds.")
    print(f"[DEBUG][retrosyn_main] Returning LLM recommendation output (length: {len(recommend_reactions_llm_output)}).")
    
    return recommend_reactions_llm_output


# Helper functions (from main.py, with added debug logs)
def countNodes(tree):
    if not tree or not hasattr(tree, 'get_node_count'): return 0
    try: return tree.get_node_count()
    except Exception: return 0

def searchPathways(tree):
    if not tree or not hasattr(tree, 'find_all_paths'): return []
    try: return tree.find_all_paths()
    except Exception: return []

def recommendReactions(prompt, result_folder_name, response_name, timeout=300): # timeout parameter is still here but won't be used in the call
    t_start = time.time()
    # print(f"[DEBUG][recommendReactions] Requesting recommendation. Response name: {response_name}, timeout: {timeout}s")
    res = "Error: LLM call failed."
    try:
        # print(f"[DEBUG][recommendReactions] Calling GPTAPI().answer_wo_vision...")
        llm_start = time.time()
        api = GPTAPI()
        # MODIFICATION: Removed timeout=timeout from the call
        res = api.answer_wo_vision(prompt)
        llm_duration = time.time() - llm_start
        # print(f"[DEBUG][recommendReactions] GPTAPI call completed in {llm_duration:.4f} seconds.")

        result_file_path = os.path.join(result_folder_name, f'{response_name}.txt')
        # print(f"[DEBUG][recommendReactions] Saving response to {result_file_path}...")
        io_start = time.time()
        os.makedirs(os.path.dirname(result_file_path), exist_ok=True)
        with open(result_file_path, 'w', encoding='utf-8') as f:
            f.write(res)
        io_duration = time.time() - io_start
        # print(f"[DEBUG][recommendReactions] Response saved in {io_duration:.4f} seconds.")

        print(f"\n[INFO][recommendReactions] Full LLM Response ({response_name}):\n======================================================\n{res}\n======================================================\n")

        duration = time.time() - t_start
        # print(f"[DEBUG][recommendReactions] Finished in {duration:.4f} seconds.")
        return res

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


class RetroSynthesisTool(BaseTool):
    name: str = "retrosynthesis"
    description: str = "Performs retrosynthesis analysis on a chemical compound. Input must be the chemical name."

    def _run(self, compound_name: str) -> Dict[str, Any]:
        """Run retrosynthesis analysis on the given compound"""
        tool_start_time = time.time()
        print(f"\n{'#'*30} [DEBUG] Entering RetroSynthesisTool._run {'#'*30}")
        print(f"[DEBUG][RetroSynthesisTool] Received compound_name: '{compound_name}'")

        if not isinstance(compound_name, str) or not compound_name.strip():
             return { "status": "error", "message": "Invalid input: compound_name must be a non-empty string." }
        
        try:
            raw_llm_output = retrosyn_main(
                material=compound_name, num_results=10, alignment=True, 
                expansion=True, filtration=False, recommend_type="general"
            )

            if not raw_llm_output or raw_llm_output.startswith("Error:") or "No valid reaction pathways found" in raw_llm_output :
                 return { "status": "error", "message": f"Retrosynthesis pipeline failed or found no pathways. Details: {raw_llm_output or 'No output'}"}

            parsed_data = parse_reaction_data(raw_llm_output)

            if "parsing_error" in parsed_data:
                print(f"[WARN][RetroSynthesisTool] Parsing error: {parsed_data['parsing_error']}.")

            # Filter the display of ranked_pathways (the ones passed to UI)
            ranked_pathways_from_llm_parse = parsed_data.get("ranked_pathways_details", [])
            displayable_ranked_pathways = []
            
            actual_pathways_for_display = []
            placeholder_pathways_for_display = []

            for pw in ranked_pathways_from_llm_parse:
                indices_str = pw.get("pathway_indices_str", "")
                if "Not available" in indices_str or "[Hypothetical" in indices_str:
                    placeholder_pathways_for_display.append(pw)
                else:
                    actual_pathways_for_display.append(pw)
            
            actual_pathways_for_display.sort(key=lambda p: p.get('rank', float('inf')))
            placeholder_pathways_for_display.sort(key=lambda p: p.get('rank', float('inf')))

            current_rank_counter = 1
            for apw in actual_pathways_for_display:
                if len(displayable_ranked_pathways) < 3:
                    apw['rank'] = current_rank_counter # Re-assign rank for consistent display
                    displayable_ranked_pathways.append(apw)
                    current_rank_counter += 1
                else: break
            
            # If less than 3 actual pathways, fill with placeholders if available
            if len(displayable_ranked_pathways) < 3:
                for phpw in placeholder_pathways_for_display:
                    if len(displayable_ranked_pathways) < 3:
                        phpw['rank'] = current_rank_counter
                        displayable_ranked_pathways.append(phpw)
                        current_rank_counter += 1
                    else: break
            
            # All unique reaction objects found by the parser
            all_unique_parsed_reactions_list = parsed_data.get("reactions", []) 
            reactions_by_idx_map = {rxn['idx']: rxn for rxn in all_unique_parsed_reactions_list}
            
            # Indices for the Rank 1 pathway, e.g., ['idx31']
            recommended_indices_for_rank1 = parsed_data.get("recommended_pathway_indices", []) 
            overall_reasoning_for_display = parsed_data.get("reasoning", "")
            
            # Construct the list of reaction objects specifically for the Rank 1 pathway
            reactions_for_rank1_display = []
            if recommended_indices_for_rank1:
                for rec_idx in recommended_indices_for_rank1:
                    if rec_idx in reactions_by_idx_map:
                        reactions_for_rank1_display.append(reactions_by_idx_map[rec_idx])
            
            # SMILES cleaning should be done on ALL unique parsed reactions
            # because the UI might need to display any of them (e.g., if user clicks on a different ranked pathway)
            # For now, we only explicitly display Rank 1 steps, but good to have all cleaned.
            for reaction_obj in all_unique_parsed_reactions_list:
                original_smiles_from_llm = reaction_obj.get("reaction_smiles")
                final_reaction_smiles_str = None
                smiles_generation_attempted = False
                if original_smiles_from_llm: # Try to use LLM's SMILES if it looks plausible
                    temp_s = original_smiles_from_llm.replace("[reactant_SMILES]", "").replace("[product_SMILES]", "")
                    temp_s = re.sub(r'\[[a-zA-Z_]+_SMILES\]', "", temp_s).strip()
                    if ">>" in temp_s and not re.search(r'\[|\]', temp_s): # No brackets in names, only in SMILES parts
                        parts = temp_s.split(">>")
                        if len(parts) == 2 and not any(' ' in p for p in parts): # No spaces in names
                                # Basic check for SMILES-like characters (can be improved)
                                if all(re.fullmatch(r'[A-Za-z0-9@\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}\[\]]*', p_part) for p_list_str in (parts[0], parts[1]) for p_part in p_list_str.split('.')):
                                    final_reaction_smiles_str = temp_s
                if not final_reaction_smiles_str: # If LLM SMILES not usable, generate
                    smiles_generation_attempted = True
                    reactant_names = reaction_obj.get("reactants", [])
                    product_names = reaction_obj.get("products", [])
                    reactant_smi_list = [get_smiles_for_compound(name) or name for name in reactant_names]
                    product_smi_list = [get_smiles_for_compound(name) or name for name in product_names]
                    all_react_valid_smi = all(smi and not ' ' in smi for smi in reactant_smi_list)
                    all_prod_valid_smi = all(smi and not ' ' in smi for smi in product_smi_list)
                    if reactant_smi_list and product_smi_list and all_react_valid_smi and all_prod_valid_smi:
                        final_reaction_smiles_str = f"{'.'.join(reactant_smi_list)}>>{'.'.join(product_smi_list)}"
                reaction_obj["cleaned_reaction_smiles"] = final_reaction_smiles_str
                reaction_obj["smiles_generation_attempted"] = smiles_generation_attempted
            
            final_result = {
                "status": "success",
                "data": {
                    "target_compound": compound_name,
                    "reactions": reactions_for_rank1_display, # Reactions for Rank 1 pathway display (primary)
                    "all_parsed_reactions_by_idx": reactions_by_idx_map, # Map of ALL unique reactions by their idx
                    "reasoning": overall_reasoning_for_display, 
                    "recommended_indices": recommended_indices_for_rank1, # Indices for Rank 1 pathway
                    "ranked_pathways": displayable_ranked_pathways, # Filtered list of ranked pathways for UI
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
            return {
                "status": "error",
                "message": f"Unhandled error in RetroSynthesisTool: {str(e)}",
                "traceback": traceback.format_exc()
            }

def run_retrosynthesis(compound_name):
    """
    Run the RetroSynthesis Agent for a given compound and process the results.
    Uses the RetroSynthesisTool internally.
    """
    print(f"\n[DEBUG][run_retrosynthesis] Direct call initiated for compound: '{compound_name}'")
    t_start = time.time()
    try:
        tool = RetroSynthesisTool()
        # print("[DEBUG][run_retrosynthesis] RetroSynthesisTool instantiated.")
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

script_duration = time.time() - start_time
print(f"[DEBUG] retrosynthesis.py script initialization and definition complete. Total load time: {script_duration:.4f} seconds.")

# # Example Usage (optional, for testing)
# if _name_ == "_main_":
#     print("\n[INFO] Running example usage directly in retrosynthesis.py...")
#     # test_compound = "aspirin"
#     # test_compound = "paracetamol"
#     test_compound = "flubendiamide" 
#     # test_compound = "polyimide" # Example from main.py if you want to test similar flow
#     print(f"[INFO] Target compound: {test_compound}")

#     results = run_retrosynthesis(test_compound)

#     print("\n[INFO] Example Usage Result:")
#     print(json.dumps(results, indent=2, ensure_ascii=False))

#     if results.get("status") == "success" and results.get("data"):
#         print("\n[INFO] Retrosynthesis ran successfully.")
#         data = results["data"]
#         print(f"  - Target: {data.get('target_compound')}")
#         print(f"  - Found {len(data.get('reactions', []))} reactions in the primary recommended pathway.")
#         print(f"  - Total parsed reactions by idx: {len(data.get('all_parsed_reactions_by_idx', {}))}")
#         print(f"  - Recommended indices for top path: {data.get('recommended_indices', [])}")
#         # Print details of the first recommended reaction if available
#         if data.get('reactions'):
#             first_rec_rxn = data['reactions'][0]
#             print(f"  - First recommended reaction (idx {first_rec_rxn.get('idx')}):")
#             print(f"    - Reactants: {first_rec_rxn.get('reactants')}")
#             print(f"    - Products: {first_rec_rxn.get('products')}")
#             print(f"    - Original LLM SMILES: {first_rec_rxn.get('reaction_smiles')}")
#             print(f"    - Cleaned/Generated SMILES: {first_rec_rxn.get('cleaned_reaction_smiles')}")
#             print(f"    - SMILES Gen Attempted: {first_rec_rxn.get('smiles_generation_attempted')}")
#         print(f"  - Reasoning provided: {'Yes' if data.get('reasoning') else 'No'}")
#     else:
#         print("\n[INFO] Retrosynthesis encountered an error or found no data.")
#         print(f"  - Message: {results.get('message')}")
#         if results.get("traceback"):
#            print(f"  - Traceback:\n{results.get('traceback')}")
#     print("[INFO] Example usage finished.")