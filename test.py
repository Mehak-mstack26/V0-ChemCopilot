# import os
# import re
# import requests
# import api_config
# from tools.make_tools import make_tools 
# from langchain.agents import AgentExecutor, ZeroShotAgent
# from langchain.prompts import PromptTemplate
# from langchain_openai import ChatOpenAI
# import pandas as pd
# from tools.asckos import ReactionClassifier

# # Setup LLM and tools
# llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=4000)  # Set max_tokens to control output size
# tools = make_tools(llm=llm)

# # Initialize reaction classifier with dataset paths
# dataset_path1 = os.environ.get('REACTION_DATASET_PATH1', None)
# dataset_path2 = os.environ.get('REACTION_DATASET_PATH2', None)
# reaction_classifier = ReactionClassifier(dataset_path1, dataset_path2)

# # Prompt parts
# PREFIX = """
# You are Chem Copilot, an expert chemistry assistant. You have access to the following tools to analyze molecules and chemical reactions.

# Always begin by understanding the user's **intent** — what kind of information are they asking for?

# Here is how to choose tools:

# - If the user gives a SMILES or reaction SMILES and asks for the name, you MUST use **SMILES2Name**. Do NOT analyze bonds or functional groups for this task.
# - Use **NameToSMILES**: when the user gives a compound/reaction name and wants the SMILES or structure.
# - Use **FuncGroups**: when the user wants to analyze functional groups in a molecule or a reaction (input is SMILES or reaction SMILES).
# - Use **BondChangeAnalyzer**: when the user asks for which bonds are broken, formed, or changed in a chemical reaction.

# If the user wants all of the above (full analysis), respond with "This requires full analysis." (This will be handled by a separate function.)

# Always return your answer in this format:
# Final Answer: <your answer here>

# For **FuncGroups** results:
# - Always list the functional groups identified in each reactant and product separately
# - Include the transformation summary showing disappeared groups, appeared groups, and unchanged groups
# - Provide a clear conclusion about what transformation occurred in the reaction

# For **BondChangeAnalyzer** results:
# - Always list the specific bonds that were broken, formed, or changed with their bond types
# - Include the atom types involved in each bond (e.g., C-O, N-H)
# - Provide a clear conclusion summarizing the key bond changes in the reaction
# """

# FORMAT_INSTRUCTIONS = """
# You can only respond with a single complete
# "Thought, Action, Action Input" format
# OR a single "Final Answer" format

# Complete format:

# Thought: (reflect on your progress and decide what to do next)
# Action: (the action name, should be one of [{tool_names}])
# Action Input: (the input string to the action)

# OR

# Final Answer: (the final answer to the original input question)
# """

# SUFFIX = """
# Question: {input}
# {agent_scratchpad}
# """

# prompt = ZeroShotAgent.create_prompt(
#     tools=tools,
#     prefix=PREFIX,
#     suffix=SUFFIX,
#     format_instructions=FORMAT_INSTRUCTIONS,
#     input_variables=["input", "agent_scratchpad"]
# )

# agent_chain = ZeroShotAgent.from_llm_and_tools(llm=llm, tools=tools, prompt=prompt)
# agent = AgentExecutor(agent=agent_chain, tools=tools, verbose=True)

# def extract_final_answer(full_output: str):
#     match = re.search(r"Final Answer:\s*(.*)", full_output, re.DOTALL)
#     return match.group(1).strip() if match else full_output.strip()

# # Function to query dataset for detailed reaction information
# def query_reaction_dataset(reaction_smiles):
#     """
#     Query the dataset for specific information about a reaction
#     based on the reaction SMILES string.
    
#     Only returns specific needed fields instead of the entire dataset row.
#     """
#     try:
#         # Get dataset reference from reaction classifier
#         if hasattr(reaction_classifier, 'dataset1') and reaction_classifier.dataset1 is not None:
#             df = reaction_classifier.dataset1
#         elif hasattr(reaction_classifier, 'dataset2') and reaction_classifier.dataset2 is not None:
#             df = reaction_classifier.dataset2
#         else:
#             return None
        
#         if df is None or df.empty:
#             return None
        
#         # Fields we're interested in
#         fields_to_extract = [
#             'procedure_details', 'rxn_time', 'temperature', 'yield_000',
#             'reaction_name', 'reaction_classname', 'prediction_certainty'
#         ]
        
#         # Try to find exact match for reaction SMILES
#         smiles_columns = ['rxn_str', 'reaction_smiles', 'smiles', 'rxn_smiles']
        
#         exact_match = None
#         for col in smiles_columns:
#             if col in df.columns:
#                 temp_match = df[df[col] == reaction_smiles]
#                 if not temp_match.empty:
#                     exact_match = temp_match
#                     break
        
#         if exact_match is not None and not exact_match.empty:
#             # Extract only the specific columns we need to reduce memory usage
#             row = exact_match.iloc[0]
            
#             # Create a dictionary with only the columns we need
#             result = {}
            
#             # Add only necessary fields
#             for field in fields_to_extract:
#                 if field in row.index and row[field] is not None and str(row[field]) != "nan":
#                     result[field] = str(row[field])
            
#             # Extract solvents (limit to just 3)
#             solvent_count = 0
#             for i in range(11):  # solvent_000 to solvent_010
#                 solvent_key = f'solvent_{i:03d}'
#                 if solvent_key in row.index and row[solvent_key] is not None and str(row[solvent_key]) != "nan":
#                     result[solvent_key] = str(row[solvent_key])
#                     solvent_count += 1
#                     if solvent_count >= 3:  # Limit to 3 solvents
#                         break
            
#             # Extract agents/catalysts (limit to just 3)
#             agent_count = 0
#             for i in range(16):  # agent_000 to agent_015
#                 agent_key = f'agent_{i:03d}'
#                 if agent_key in row.index and row[agent_key] is not None and str(row[agent_key]) != "nan":
#                     result[agent_key] = str(row[agent_key])
#                     agent_count += 1
#                     if agent_count >= 3:  # Limit to 3 agents
#                         break
            
#             return result
        
#         # If no exact match, return None
#         return None
        
#     except Exception as e:
#         print(f"Error querying reaction dataset: {e}")
#         return None

# # 🧠 Full Info Handler - Optimized to reduce context size
# def handle_full_info(query, reaction_smiles):
#     print("Running full analysis using all tools...\n")
#     full_info = {}

#     # Map tool names to classes
#     tool_dict = {tool.name.lower(): tool for tool in tools}

#     try:
#         # First, visualize the reaction
#         visualizer_tool = tool_dict.get("chemvisualizer")
#         if visualizer_tool:
#             try:
#                 full_info['Visualization'] = visualizer_tool.run(reaction_smiles)
#             except Exception as e:
#                 full_info['Visualization'] = f"Error visualizing reaction: {str(e)}"
#         else:
#             full_info['Visualization'] = "ChemVisualizer tool not found"

#         # Run smiles2name - with limited output
#         name_tool = tool_dict.get("smiles2name")
#         if name_tool:
#             try:
#                 name_result = name_tool.run(reaction_smiles)
#                 # Extract only essential name information
#                 if isinstance(name_result, str) and len(name_result) > 500:
#                     # Truncate overly verbose results
#                     name_result = name_result[:500] + "... [truncated for brevity]"
#                 full_info['Names'] = name_result
#             except Exception as e:
#                 full_info['Names'] = f"Error analyzing names: {str(e)}"
#         else:
#             full_info['Names'] = "SMILES2Name tool not found"

#         # Run funcgroups - with limited output
#         fg_tool = tool_dict.get("funcgroups")
#         if fg_tool:
#             try:
#                 fg_result = fg_tool.run(reaction_smiles)
#                 # Extract only key functional group information
#                 if isinstance(fg_result, str) and len(fg_result) > 1000:
#                     # Keep just essential parts
#                     lines = fg_result.split('\n')
#                     important_sections = []
#                     # Get transformation summary (usually the most important part)
#                     for i, line in enumerate(lines):
#                         if "transformation summary" in line.lower() or "appeared groups" in line.lower():
#                             important_sections.extend(lines[i:i+15])  # Take this section and a few lines
#                             break
                    
#                     if important_sections:
#                         fg_result = '\n'.join(important_sections)
#                     else:
#                         fg_result = fg_result[:1000] + "... [truncated for brevity]"
                
#                 full_info['Functional Groups'] = fg_result
#             except Exception as e:
#                 full_info['Functional Groups'] = f"Error analyzing functional groups: {str(e)}"
#         else:
#             full_info['Functional Groups'] = "FuncGroups tool not found"

#         # Run bond.py - with limited output
#         bond_tool = tool_dict.get("bondchangeanalyzer")
#         if bond_tool:
#             try:
#                 bond_result = bond_tool.run(reaction_smiles)
#                 # Extract only key bond changes
#                 if isinstance(bond_result, str) and len(bond_result) > 1000:
#                     # Keep just essential parts
#                     bond_result = bond_result[:1000] + "... [truncated for brevity]"
#                 full_info['Bond Changes'] = bond_result
#             except Exception as e:
#                 full_info['Bond Changes'] = f"Error analyzing bond changes: {str(e)}"
#         else:
#             full_info['Bond Changes'] = "BondChangeAnalyzer tool not found"
            
#         # Run reaction classifier - with limited output
#         if reaction_classifier:
#             try:
#                 classifier_result = reaction_classifier._run(reaction_smiles)
#                 # Extract only essential classification info
#                 if isinstance(classifier_result, str):
#                     # Extract just the classification summary
#                     summary_match = re.search(r'## Summary\n(.*?)$', classifier_result, re.DOTALL)
#                     if summary_match:
#                         classifier_summary = summary_match.group(1).strip()
#                     else:
#                         # If no summary section, just take the beginning
#                         classifier_summary = classifier_result[:500] + "... [truncated for brevity]"
                    
#                     full_info['Reaction Classification'] = classifier_summary
#                 else:
#                     full_info['Reaction Classification'] = "Classification result not in expected format"
#             except Exception as e:
#                 full_info['Reaction Classification'] = f"Error classifying reaction: {str(e)}"
#         else:
#             full_info['Reaction Classification'] = "ReactionClassifier tool not found"
        
#         # Get additional data from dataset - but ONLY specific fields needed
#         dataset_info = query_reaction_dataset(reaction_smiles)
        
#         # Extract only necessary information
#         procedure_details = None
#         rxn_time = None
#         temperature = None
#         yield_info = None
#         solvents = []
#         agents = []
        
#         if dataset_info:
#             # Only extract the specific fields we need
#             if 'procedure_details' in dataset_info and dataset_info['procedure_details'] != "nan":
#                 # Truncate procedure details if too long
#                 if len(dataset_info['procedure_details']) > 500:
#                     procedure_details = dataset_info['procedure_details'][:500] + "... [truncated]"
#                 else:
#                     procedure_details = dataset_info['procedure_details']
                
#             if 'rxn_time' in dataset_info and dataset_info['rxn_time'] != "nan":
#                 rxn_time = dataset_info['rxn_time']
                
#             if 'temperature' in dataset_info and dataset_info['temperature'] != "nan":
#                 temperature = dataset_info['temperature']
                
#             if 'yield_000' in dataset_info and dataset_info['yield_000'] != "nan":
#                 yield_info = dataset_info['yield_000']
                
#             # Extract solvents (limited to 3)
#             for i in range(3):  # Just check the first 3 solvent columns
#                 solvent_key = f'solvent_{i:03d}'
#                 if solvent_key in dataset_info and dataset_info[solvent_key] and dataset_info[solvent_key] != "nan":
#                     solvents.append(dataset_info[solvent_key])
                    
#             # Extract agents/catalysts (limited to 3)
#             for i in range(3):  # Just check the first 3 agent columns
#                 agent_key = f'agent_{i:03d}'
#                 if agent_key in dataset_info and dataset_info[agent_key] and dataset_info[agent_key] != "nan":
#                     agents.append(dataset_info[agent_key])

#         # Create a shorter, more focused prompt for the LLM
#         final_prompt = f"""You are a chemistry expert. Synthesize this reaction analysis into a clear explanation:
        
# Reaction SMILES: {reaction_smiles}

# The following analysis provides different perspectives on this reaction:

# NAMES: {full_info.get('Names', 'Not available')}

# BOND CHANGES: {full_info.get('Bond Changes', 'Not available')}

# FUNCTIONAL GROUPS: {full_info.get('Functional Groups', 'Not available')}

# REACTION TYPE: {full_info.get('Reaction Classification', 'Not available')}
# """
#         # Add dataset information to prompt - but only if important and keeping it concise
#         if procedure_details:
#             procedure_summary = procedure_details
#             final_prompt += f"\nPROCEDURE: {procedure_summary}\n"
            
#         if rxn_time or temperature or yield_info:
#             conditions = []
#             if temperature: conditions.append(f"Temperature: {temperature}")
#             if rxn_time: conditions.append(f"Time: {rxn_time}")
#             if yield_info: conditions.append(f"Yield: {yield_info}")
#             final_prompt += f"\nCONDITIONS: {', '.join(conditions)}\n"
            
#         if solvents or agents:
#             materials = []
#             if solvents: materials.append(f"Solvents: {', '.join(solvents)}")
#             if agents: materials.append(f"Catalysts: {', '.join(agents)}")
#             final_prompt += f"\nMATERIALS: {', '.join(materials)}\n"

#         final_prompt += """
# Provide a thorough, well-structured explanation of this reaction that:
# 1. Begins with a high-level summary of what type of reaction this is
# 2. Explains what happens at the molecular level (bonds broken/formed)
# 3. Discusses the functional group transformations
# 4. Includes specific experimental conditions (temperature, time, yield, solvents, catalysts)
# 5. Provides a detailed procedure summary if available
# 6. Mentions common applications or importance of this reaction type

# Please give a complete and readable explanation of this reaction.
# """

#         # Use a model instance with controlled output size
#         focused_llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=2000)
#         response = focused_llm.invoke(final_prompt)
        
#         # Structure the response data but don't include the full dataset
#         return {
#             'visualization_path': full_info.get('Visualization') if 'Visualization' in full_info and not full_info['Visualization'].startswith('Error') else None,
#             'analysis': response.content.strip(),
#             'reaction_classification': full_info.get('Reaction Classification', "No classification available"),
#             'procedure_details': procedure_details,
#             'reaction_time': rxn_time,
#             'temperature': temperature,
#             'yield': yield_info,
#             'solvents': solvents if solvents else None,
#             'agents_catalysts': agents if agents else None
#         }
    
#     except Exception as e:
#         return {
#             'visualization_path': None,
#             'analysis': f"Error in full analysis: {str(e)}"
#         }
    
# def handle_reaction_query(query, reaction_smiles):
#     """
#     Handle specific property queries for a reaction.
    
#     Args:
#         query: User query text
#         reaction_smiles: The reaction SMILES string
        
#     Returns:
#         Dictionary with visualization path and analysis results
#     """
#     try:
#         # Detect property-specific queries
#         property_keywords = {
#             'temperature': ['temperature', 'temp', 'heat', 'celsius', 'fahrenheit', 'kelvin', 'degrees', '°C', '°F', '°K'],
#             'yield': ['yield', 'yields', 'percentage', 'amount', 'efficiency', 'how much'],
#             'solvent': ['solvent', 'medium', 'solution', 'dissolved in'],
#             'catalyst': ['catalyst', 'catalysts', 'catalytic', 'accelerator'],
#             'time': ['time', 'duration', 'how long', 'minutes', 'hours', 'days'],
#             'pressure': ['pressure', 'bar', 'atm', 'psi'],
#             'ph': ['ph', 'acidity', 'alkalinity', 'basic', 'acidic'],
#             'mechanism': ['mechanism', 'how does it work', 'how it works', 'pathway'],
#             'applications': ['applications', 'uses', 'what is it used for', 'industry', 'purpose'],
#             'procedure': ['procedure', 'protocol', 'how to', 'steps', 'instructions', 'method'],
#             'safety': ['safety', 'hazard', 'danger', 'precaution', 'handling'],
#             'alternatives': ['alternative', 'other way', 'different method', 'substitute'],
#             'optimization': ['optimize', 'improve', 'enhancement', 'better yield', 'efficiency']
#         }
        
#         # Normalize query
#         query_lower = query.lower()
        
#         # Get only specific dataset information for the reaction
#         dataset_info = query_reaction_dataset(reaction_smiles)
        
#         # Check for property-specific queries
#         for property_name, keywords in property_keywords.items():
#             if any(keyword in query_lower for keyword in keywords):
#                 # First try to get information from dataset
#                 property_value = None
                
#                 if dataset_info:
#                     if property_name == 'temperature' and 'temperature' in dataset_info:
#                         property_value = dataset_info['temperature']
#                     elif property_name == 'yield' and 'yield_000' in dataset_info:
#                         property_value = dataset_info['yield_000']
#                     elif property_name == 'time' and 'rxn_time' in dataset_info:
#                         property_value = dataset_info['rxn_time']
#                     elif property_name == 'procedure' and 'procedure_details' in dataset_info:
#                         # Truncate procedure if too long
#                         if len(dataset_info['procedure_details']) > 500:
#                             property_value = dataset_info['procedure_details'][:500] + "... [truncated]"
#                         else:
#                             property_value = dataset_info['procedure_details']
#                     elif property_name == 'solvent':
#                         solvents = []
#                         for i in range(3):  # Just check first 3 solvent columns
#                             solvent_key = f'solvent_{i:03d}'
#                             if solvent_key in dataset_info and dataset_info[solvent_key] and dataset_info[solvent_key] != "nan":
#                                 solvents.append(dataset_info[solvent_key])
#                         if solvents:
#                             property_value = ', '.join(solvents)
#                     elif property_name == 'catalyst':
#                         agents = []
#                         for i in range(3):  # Just check first 3 agent columns
#                             agent_key = f'agent_{i:03d}'
#                             if agent_key in dataset_info and dataset_info[agent_key] and dataset_info[agent_key] != "nan":
#                                 agents.append(dataset_info[agent_key])
#                         if agents:
#                             property_value = ', '.join(agents)
                
#                 # If we found data in dataset, use it along with classifier response
#                 if property_value and property_value != "nan":
#                     # Use classifier for supplementary information but limit response size
#                     try:
#                         classifier_result = reaction_classifier.query_reaction_property(reaction_smiles, property_name)
#                         # Ensure the classifier result isn't too large
#                         if isinstance(classifier_result, str) and len(classifier_result) > 1000:
#                             classifier_result = classifier_result[:1000] + "... [truncated for brevity]"
#                     except Exception:
#                         classifier_result = "Unable to retrieve additional information from the classifier."
                    
#                     # Create a focused prompt
#                     prompt = f"""Based on the reaction SMILES: {reaction_smiles}
                    
# The dataset shows the following {property_name} information: {property_value}

# Additional analysis: {classifier_result}

# Provide a brief, focused answer about the {property_name} for this reaction, incorporating the actual data from the dataset. Keep your response under 1500 characters.
# """
#                     focused_llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=1500)
#                     response = focused_llm.invoke(prompt)
#                     return {
#                         "visualization_path": None,
#                         "analysis": response.content.strip()
#                     }
#                 else:
#                     # If no dataset info, fall back to classifier alone but limit output
#                     try:
#                         result = reaction_classifier.query_reaction_property(reaction_smiles, property_name)
#                         # Truncate if needed
#                         if isinstance(result, str) and len(result) > 2000:
#                             result = result[:2000] + "... [truncated for brevity]"
#                     except Exception as e:
#                         result = f"Error retrieving property information: {str(e)}"
                    
#                     return {
#                         "visualization_path": None,
#                         "analysis": result
#                     }
        
#         # For general queries, combine dataset info with reaction classification in a focused way
#         if dataset_info:
#             # Extract key information from dataset - but be concise
#             dataset_summary = "Key experimental data:\n\n"
            
#             if 'procedure_details' in dataset_info and dataset_info['procedure_details'] != "nan":
#                 # Truncate procedure details if too long
#                 procedure = dataset_info['procedure_details']
#                 if len(procedure) > 300:
#                     procedure = procedure[:300] + "... [truncated]"
#                 dataset_summary += f"- **Procedure**: {procedure}\n"
            
#             # Include key conditions
#             conditions = []
#             if 'temperature' in dataset_info and dataset_info['temperature'] != "nan":
#                 conditions.append(f"Temperature: {dataset_info['temperature']}")
                
#             if 'rxn_time' in dataset_info and dataset_info['rxn_time'] != "nan":
#                 conditions.append(f"Reaction Time: {dataset_info['rxn_time']}")
                
#             if 'yield_000' in dataset_info and dataset_info['yield_000'] != "nan":
#                 conditions.append(f"Yield: {dataset_info['yield_000']}")
                
#             if conditions:
#                 dataset_summary += f"- **Conditions**: {', '.join(conditions)}\n"
            
#             # Collect solvents (limited)
#             solvents = []
#             for i in range(3):  # Just first 3
#                 solvent_key = f'solvent_{i:03d}'
#                 if solvent_key in dataset_info and dataset_info[solvent_key] and dataset_info[solvent_key] != "nan":
#                     solvents.append(dataset_info[solvent_key])
#             if solvents:
#                 dataset_summary += f"- **Solvents**: {', '.join(solvents)}\n"
            
#             # Collect agents/catalysts (limited)
#             agents = []
#             for i in range(3):  # Just first 3
#                 agent_key = f'agent_{i:03d}'
#                 if agent_key in dataset_info and dataset_info[agent_key] and dataset_info[agent_key] != "nan":
#                     agents.append(dataset_info[agent_key])
#             if agents:
#                 dataset_summary += f"- **Agents/Catalysts**: {', '.join(agents)}\n"
            
#             # Get classification response - but extract only essential info
#             try:
#                 classifier_result = reaction_classifier._run(reaction_smiles)
                
#                 # Extract just the most important parts
#                 if isinstance(classifier_result, str):
#                     # First try to get the summary section
#                     summary_match = re.search(r'## Summary\n(.*?)$', classifier_result, re.DOTALL)
#                     if summary_match:
#                         classifier_info = "Classification Summary:\n" + summary_match.group(1).strip()
#                     else:
#                         # If no summary, extract key classification info
#                         classification_match = re.search(r'## Reaction Classification\n(.*?)(?=##|\Z)', classifier_result, re.DOTALL)
#                         if classification_match:
#                             classifier_info = "Classification:\n" + classification_match.group(1).strip()
#                         else:
#                             # Just take first 500 chars if specific sections not found
#                             classifier_info = classifier_result[:500] + "... [truncated]"
#                 else:
#                     classifier_info = "Classification information not available in expected format."
#             except Exception as e:
#                 classifier_info = f"Error retrieving classification: {str(e)}"
            
#             # Create focused prompt
#             prompt = f"""Answer the following query about a chemical reaction:

# Query: {query}

# Reaction SMILES: {reaction_smiles}

# {dataset_summary}

# {classifier_info}

# Please provide a focused answer addressing the query. Prioritize experimental data. Keep your response under 2000 characters.
# """
#             focused_llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=2000)
#             response = focused_llm.invoke(prompt)
#             return {
#                 "visualization_path": None,
#                 "analysis": response.content.strip()
#             }
#         else:
#             # If no dataset info, fall back to classifier alone
#             try:
#                 result = reaction_classifier._run(reaction_smiles)
#                 # Ensure result isn't too large
#                 if isinstance(result, str) and len(result) > 3000:
#                     result = result[:3000] + "... [truncated for brevity]"
#             except Exception as e:
#                 result = f"Error retrieving reaction information: {str(e)}"
                
#             return {
#                 "visualization_path": None,
#                 "analysis": result
#             }
            
#     except Exception as e:
#         return {
#             "visualization_path": None,
#             "analysis": f"Error processing query: {str(e)}"
#         }
    
# def enhanced_query(query, callbacks=None):
#     try:
#         # Check if it's a full analysis request
#         if "full information" in query.lower() or "full analysis" in query.lower():
#             # Extract SMILES from query using regex - handle both simple and complex SMILES
#             # More comprehensive SMILES pattern
#             smiles_pattern = r"([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]*)"
#             match = re.search(smiles_pattern, query)
            
#             if match:
#                 reaction_smiles = match.group(1)
#                 return handle_full_info(query, reaction_smiles)
#             else:
#                 # Try to extract product SMILES if not a full reaction
#                 product_pattern = r"rxn\s+([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+)"
#                 product_match = re.search(product_pattern, query)
                
#                 if product_match:
#                     # For single product, create a simple reaction SMILES
#                     product_smiles = product_match.group(1)
#                     # Simple reaction with generic reactant
#                     return handle_full_info(query, f"[R]>>{product_smiles}")
                
#                 return {"visualization_path": None, "analysis": "Could not extract reaction SMILES from the query. Please provide a valid reaction SMILES or product SMILES."}
        
#         # Check if it's just a visualization request
#         elif any(term in query.lower() for term in ["visual", "picture", "image", "show", "draw", "representation"]):
#             # Extract SMILES from query using regex
#             smiles_pattern = r"([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]*)"
#             match = re.search(smiles_pattern, query)
            
#             if match:
#                 reaction_smiles = match.group(1)
#                 # Only use the visualizer tool
#                 tool_dict = {tool.name.lower(): tool for tool in tools}
#                 visualizer_tool = tool_dict.get("chemvisualizer")
                
#                 if visualizer_tool:
#                     try:
#                         viz_path = visualizer_tool.run(reaction_smiles)
#                         return {
#                             "visualization_path": viz_path,
#                             "analysis": f"Visual representation of the reaction: {reaction_smiles}"
#                         }
#                     except Exception as e:
#                         return {"visualization_path": None, "analysis": f"Error visualizing reaction: {str(e)}"}
#                 else:
#                     return {"visualization_path": None, "analysis": "ChemVisualizer tool not found"}
#             else:
#                 return {"visualization_path": None, "analysis": "Could not extract reaction SMILES from the visualization request. Please provide a valid reaction SMILES."}
        
#         # Check if it's a query about a specific reaction property
#         elif "reaction" in query.lower() and any(keyword in query.lower() for keyword in [
#             "temperature", "yield", "solvent", "catalyst", "time", "pressure", "ph", "mechanism", 
#             "applications", "procedure", "protocol", "uses", "purpose"
#         ]):
#             # Extract SMILES from query using regex
#             smiles_pattern = r"(.*?reaction\s+(?:SMILES)?[:=]?\s*([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]*))"
#             match = re.search(smiles_pattern, query, re.IGNORECASE)
            
#             if match:
#                 reaction_smiles = match.group(2)
#                 return handle_reaction_query(query, reaction_smiles)
#             else:
#                 # Try alternative patterns
#                 alt_pattern = r"(.*?(?:for|about)\s+the\s+reaction\s+([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]*))"
#                 alt_match = re.search(alt_pattern, query, re.IGNORECASE)
                
#                 if alt_match:
#                     reaction_smiles = alt_match.group(2)
#                     return handle_reaction_query(query, reaction_smiles)
#                 else:
#                     return {"visualization_path": None, "analysis": "Could not extract reaction SMILES from the query. Please provide a valid reaction SMILES."}

#         # Otherwise, use normal agent
#         result = agent.invoke({"input": query}, {"callbacks": callbacks} if callbacks else {})
#         return {"visualization_path": None, "analysis": extract_final_answer(result.get("output", ""))}

#     except Exception as e:
#         return {"visualization_path": None, "analysis": f"Error occurred: {str(e)}"}

# # Run
# if __name__ == "__main__":
#     # Test analysis
#     query = "Give full information about this rxn O=C(NC(=O)c1c(F)cccc1F)c2ccc(Cl)c([N+](=O)[O-])c2"
#     print(enhanced_query(query))









# NOW 
import os
import re
import requests
import api_config
from tools.make_tools import make_tools 
from langchain.agents import AgentExecutor, ZeroShotAgent
from langchain.prompts import PromptTemplate
from langchain_openai import ChatOpenAI
import pandas as pd
from tools.asckos import ReactionClassifier
import json
from functools import lru_cache

# Setup LLM and tools
llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=4000)  # Set max_tokens to control output size
tools = make_tools(llm=llm)

# Initialize reaction classifier with dataset paths
dataset_path1 = os.environ.get('REACTION_DATASET_PATH1', None)
dataset_path2 = os.environ.get('REACTION_DATASET_PATH2', None)
reaction_classifier = ReactionClassifier(dataset_path1, dataset_path2)

# Create a global cache for reaction data
# This will store all processed data for each reaction SMILES
reaction_cache = {}

# Prompt parts
PREFIX = """
You are Chem Copilot, an expert chemistry assistant. You have access to the following tools to analyze molecules and chemical reactions.

Always begin by understanding the user's **intent** — what kind of information are they asking for?

Here is how to choose tools:

- If the user gives a SMILES or reaction SMILES and asks for the name, you MUST use **SMILES2Name**. Do NOT analyze bonds or functional groups for this task.
- Use **NameToSMILES**: when the user gives a compound/reaction name and wants the SMILES or structure.
- Use **FuncGroups**: when the user wants to analyze functional groups in a molecule or a reaction (input is SMILES or reaction SMILES).
- Use **BondChangeAnalyzer**: when the user asks for which bonds are broken, formed, or changed in a chemical reaction.

If the user wants all of the above (full analysis), respond with "This requires full analysis." (This will be handled by a separate function.)

Always return your answer in this format:
Final Answer: <your answer here>

For **FuncGroups** results:
- Always list the functional groups identified in each reactant and product separately
- Include the transformation summary showing disappeared groups, appeared groups, and unchanged groups
- Provide a clear conclusion about what transformation occurred in the reaction

For **BondChangeAnalyzer** results:
- Always list the specific bonds that were broken, formed, or changed with their bond types
- Include the atom types involved in each bond (e.g., C-O, N-H)
- Provide a clear conclusion summarizing the key bond changes in the reaction
"""

FORMAT_INSTRUCTIONS = """
You can only respond with a single complete
"Thought, Action, Action Input" format
OR a single "Final Answer" format

Complete format:

Thought: (reflect on your progress and decide what to do next)
Action: (the action name, should be one of [{tool_names}])
Action Input: (the input string to the action)

OR

Final Answer: (the final answer to the original input question)
"""

SUFFIX = """
Question: {input}
{agent_scratchpad}
"""

prompt = ZeroShotAgent.create_prompt(
    tools=tools,
    prefix=PREFIX,
    suffix=SUFFIX,
    format_instructions=FORMAT_INSTRUCTIONS,
    input_variables=["input", "agent_scratchpad"]
)

agent_chain = ZeroShotAgent.from_llm_and_tools(llm=llm, tools=tools, prompt=prompt)
agent = AgentExecutor(agent=agent_chain, tools=tools, verbose=True)

def extract_final_answer(full_output: str):
    match = re.search(r"Final Answer:\s*(.*)", full_output, re.DOTALL)
    return match.group(1).strip() if match else full_output.strip()

# Function to query dataset for detailed reaction information - with caching
@lru_cache(maxsize=100)  # Cache up to 100 recent queries
def query_reaction_dataset(reaction_smiles):
    """
    Query the dataset for specific information about a reaction
    based on the reaction SMILES string.
    
    Only returns specific needed fields instead of the entire dataset row.
    Uses caching to avoid repeated lookups.
    """
    # Check if we have this in our reaction_cache first
    if reaction_smiles in reaction_cache and 'dataset_info' in reaction_cache[reaction_smiles]:
        return reaction_cache[reaction_smiles]['dataset_info']
    
    try:
        # Get dataset reference from reaction classifier
        if hasattr(reaction_classifier, 'dataset1') and reaction_classifier.dataset1 is not None:
            df = reaction_classifier.dataset1
        elif hasattr(reaction_classifier, 'dataset2') and reaction_classifier.dataset2 is not None:
            df = reaction_classifier.dataset2
        else:
            return None
        
        if df is None or df.empty:
            return None
        
        # Fields we're interested in
        fields_to_extract = [
            'procedure_details', 'rxn_time', 'temperature', 'yield_000',
            'reaction_name', 'reaction_classname', 'prediction_certainty'
        ]
        
        # Try to find exact match for reaction SMILES
        smiles_columns = ['rxn_str', 'reaction_smiles', 'smiles', 'rxn_smiles']
        
        exact_match = None
        for col in smiles_columns:
            if col in df.columns:
                temp_match = df[df[col] == reaction_smiles]
                if not temp_match.empty:
                    exact_match = temp_match
                    break
        
        if exact_match is not None and not exact_match.empty:
            # Extract only the specific columns we need to reduce memory usage
            row = exact_match.iloc[0]
            
            # Create a dictionary with only the columns we need
            result = {}
            
            # Add only necessary fields
            for field in fields_to_extract:
                if field in row.index and row[field] is not None and str(row[field]) != "nan":
                    result[field] = str(row[field])
            
            # Extract solvents (limit to just 3)
            solvent_count = 0
            for i in range(11):  # solvent_000 to solvent_010
                solvent_key = f'solvent_{i:03d}'
                if solvent_key in row.index and row[solvent_key] is not None and str(row[solvent_key]) != "nan":
                    result[solvent_key] = str(row[solvent_key])
                    solvent_count += 1
                    if solvent_count >= 3:  # Limit to 3 solvents
                        break
            
            # Extract agents/catalysts (limit to just 3)
            agent_count = 0
            for i in range(16):  # agent_000 to agent_015
                agent_key = f'agent_{i:03d}'
                if agent_key in row.index and row[agent_key] is not None and str(row[agent_key]) != "nan":
                    result[agent_key] = str(row[agent_key])
                    agent_count += 1
                    if agent_count >= 3:  # Limit to 3 agents
                        break
            
            # Store in cache
            if reaction_smiles not in reaction_cache:
                reaction_cache[reaction_smiles] = {}
            reaction_cache[reaction_smiles]['dataset_info'] = result
            
            return result
        
        # If no exact match, return None
        if reaction_smiles not in reaction_cache:
            reaction_cache[reaction_smiles] = {}
        reaction_cache[reaction_smiles]['dataset_info'] = None
        return None
        
    except Exception as e:
        print(f"Error querying reaction dataset: {e}")
        return None

# Extract SMILES from user query
def extract_reaction_smiles(query):
    """
    Extract reaction SMILES from a user query using various patterns.
    Returns None if no valid SMILES is found.
    """
    # Pattern for reaction SMILES (reactants>>products)
    smiles_pattern = r"([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]*)"
    match = re.search(smiles_pattern, query)
    
    if match:
        return match.group(1)
    
    # Try to extract product SMILES if not a full reaction
    product_pattern = r"rxn\s+([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+)"
    product_match = re.search(product_pattern, query)
    
    if product_match:
        # For single product, create a simple reaction SMILES
        product_smiles = product_match.group(1)
        # Simple reaction with generic reactant
        return f"[R]>>{product_smiles}"
    
    # Pattern for reaction SMILES with explicit labeling
    labeled_pattern = r"(.*?reaction\s+(?:SMILES)?[:=]?\s*([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]*))"
    labeled_match = re.search(labeled_pattern, query, re.IGNORECASE)
    
    if labeled_match:
        return labeled_match.group(2)
    
    # Alternative pattern for reaction descriptions
    alt_pattern = r"(.*?(?:for|about)\s+the\s+reaction\s+([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]*))"
    alt_match = re.search(alt_pattern, query, re.IGNORECASE)
    
    if alt_match:
        return alt_match.group(2)
    
    return None

# Modified handle_full_info function with caching
def handle_full_info(query, reaction_smiles):
    """
    Process a full analysis request for a reaction, with caching to avoid 
    redundant processing for follow-up questions.
    """
    print("Running full analysis using all tools...\n")
    
    # Check if we already have this reaction in our cache
    if reaction_smiles in reaction_cache and 'full_info' in reaction_cache[reaction_smiles]:
        print(f"Using cached data for reaction: {reaction_smiles}")
        return reaction_cache[reaction_smiles]['full_info']
    
    # Initialize the cache entry if it doesn't exist
    if reaction_smiles not in reaction_cache:
        reaction_cache[reaction_smiles] = {}
    
    full_info = {}
    # Map tool names to classes
    tool_dict = {tool.name.lower(): tool for tool in tools}

    try:
        # First, visualize the reaction
        visualizer_tool = tool_dict.get("chemvisualizer")
        if visualizer_tool:
            try:
                visualization_path = visualizer_tool.run(reaction_smiles)
                full_info['Visualization'] = visualization_path
                # Store visualization path in cache
                reaction_cache[reaction_smiles]['visualization_path'] = visualization_path
            except Exception as e:
                full_info['Visualization'] = f"Error visualizing reaction: {str(e)}"
                reaction_cache[reaction_smiles]['visualization_path'] = None
        else:
            full_info['Visualization'] = "ChemVisualizer tool not found"
            reaction_cache[reaction_smiles]['visualization_path'] = None

        # Run smiles2name - with limited output
        name_tool = tool_dict.get("smiles2name")
        if name_tool:
            try:
                name_result = name_tool.run(reaction_smiles)
                # Extract only essential name information
                if isinstance(name_result, str) and len(name_result) > 500:
                    # Truncate overly verbose results
                    name_result = name_result[:500] + "... [truncated for brevity]"
                full_info['Names'] = name_result
                # Cache the name result
                reaction_cache[reaction_smiles]['name_info'] = name_result
            except Exception as e:
                full_info['Names'] = f"Error analyzing names: {str(e)}"
                reaction_cache[reaction_smiles]['name_info'] = f"Error analyzing names: {str(e)}"
        else:
            full_info['Names'] = "SMILES2Name tool not found"
            reaction_cache[reaction_smiles]['name_info'] = "SMILES2Name tool not found"

        # Run funcgroups - with limited output
        fg_tool = tool_dict.get("funcgroups")
        if fg_tool:
            try:
                fg_result = fg_tool.run(reaction_smiles)
                # Extract only key functional group information
                if isinstance(fg_result, str) and len(fg_result) > 1000:
                    # Keep just essential parts
                    lines = fg_result.split('\n')
                    important_sections = []
                    # Get transformation summary (usually the most important part)
                    for i, line in enumerate(lines):
                        if "transformation summary" in line.lower() or "appeared groups" in line.lower():
                            important_sections.extend(lines[i:i+15])  # Take this section and a few lines
                            break
                    
                    if important_sections:
                        fg_result = '\n'.join(important_sections)
                    else:
                        fg_result = fg_result[:1000] + "... [truncated for brevity]"
                
                full_info['Functional Groups'] = fg_result
                # Cache functional groups result
                reaction_cache[reaction_smiles]['fg_info'] = fg_result
            except Exception as e:
                full_info['Functional Groups'] = f"Error analyzing functional groups: {str(e)}"
                reaction_cache[reaction_smiles]['fg_info'] = f"Error analyzing functional groups: {str(e)}"
        else:
            full_info['Functional Groups'] = "FuncGroups tool not found"
            reaction_cache[reaction_smiles]['fg_info'] = "FuncGroups tool not found"

        # Run bond.py - with limited output
        bond_tool = tool_dict.get("bondchangeanalyzer")
        if bond_tool:
            try:
                bond_result = bond_tool.run(reaction_smiles)
                # Extract only key bond changes
                if isinstance(bond_result, str) and len(bond_result) > 1000:
                    # Keep just essential parts
                    bond_result = bond_result[:1000] + "... [truncated for brevity]"
                full_info['Bond Changes'] = bond_result
                # Cache bond changes result
                reaction_cache[reaction_smiles]['bond_info'] = bond_result
            except Exception as e:
                full_info['Bond Changes'] = f"Error analyzing bond changes: {str(e)}"
                reaction_cache[reaction_smiles]['bond_info'] = f"Error analyzing bond changes: {str(e)}"
        else:
            full_info['Bond Changes'] = "BondChangeAnalyzer tool not found"
            reaction_cache[reaction_smiles]['bond_info'] = "BondChangeAnalyzer tool not found"
            
        # Run reaction classifier - with limited output
        if reaction_classifier:
            try:
                classifier_result = reaction_classifier._run(reaction_smiles)
                # Extract only essential classification info
                if isinstance(classifier_result, str):
                    # Extract just the classification summary
                    summary_match = re.search(r'## Summary\n(.*?)$', classifier_result, re.DOTALL)
                    if summary_match:
                        classifier_summary = summary_match.group(1).strip()
                    else:
                        # If no summary section, just take the beginning
                        classifier_summary = classifier_result[:500] + "... [truncated for brevity]"
                    
                    full_info['Reaction Classification'] = classifier_summary
                    # Cache classification result
                    reaction_cache[reaction_smiles]['classification_info'] = classifier_summary
                else:
                    full_info['Reaction Classification'] = "Classification result not in expected format"
                    reaction_cache[reaction_smiles]['classification_info'] = "Classification result not in expected format"
            except Exception as e:
                full_info['Reaction Classification'] = f"Error classifying reaction: {str(e)}"
                reaction_cache[reaction_smiles]['classification_info'] = f"Error classifying reaction: {str(e)}"
        else:
            full_info['Reaction Classification'] = "ReactionClassifier tool not found"
            reaction_cache[reaction_smiles]['classification_info'] = "ReactionClassifier tool not found"
        
        # Get additional data from dataset - but ONLY specific fields needed
        dataset_info = query_reaction_dataset(reaction_smiles)
        
        # Extract only necessary information
        procedure_details = None
        rxn_time = None
        temperature = None
        yield_info = None
        solvents = []
        agents = []
        
        if dataset_info:
            # Only extract the specific fields we need
            if 'procedure_details' in dataset_info and dataset_info['procedure_details'] != "nan":
                # Truncate procedure details if too long
                if len(dataset_info['procedure_details']) > 500:
                    procedure_details = dataset_info['procedure_details'][:500] + "... [truncated]"
                else:
                    procedure_details = dataset_info['procedure_details']
                
            if 'rxn_time' in dataset_info and dataset_info['rxn_time'] != "nan":
                rxn_time = dataset_info['rxn_time']
                
            if 'temperature' in dataset_info and dataset_info['temperature'] != "nan":
                temperature = dataset_info['temperature']
                
            if 'yield_000' in dataset_info and dataset_info['yield_000'] != "nan":
                yield_info = dataset_info['yield_000']
                
            # Extract solvents (limited to 3)
            for i in range(3):  # Just check the first 3 solvent columns
                solvent_key = f'solvent_{i:03d}'
                if solvent_key in dataset_info and dataset_info[solvent_key] and dataset_info[solvent_key] != "nan":
                    solvents.append(dataset_info[solvent_key])
                    
            # Extract agents/catalysts (limited to 3)
            for i in range(3):  # Just check the first 3 agent columns
                agent_key = f'agent_{i:03d}'
                if agent_key in dataset_info and dataset_info[agent_key] and dataset_info[agent_key] != "nan":
                    agents.append(dataset_info[agent_key])

        # Cache these specific details
        reaction_cache[reaction_smiles]['procedure_details'] = procedure_details
        reaction_cache[reaction_smiles]['rxn_time'] = rxn_time
        reaction_cache[reaction_smiles]['temperature'] = temperature
        reaction_cache[reaction_smiles]['yield_info'] = yield_info
        reaction_cache[reaction_smiles]['solvents'] = solvents
        reaction_cache[reaction_smiles]['agents'] = agents

        # Create a shorter, more focused prompt for the LLM
        final_prompt = f"""You are a chemistry expert. Synthesize this reaction analysis into a clear explanation:
        
Reaction SMILES: {reaction_smiles}

The following analysis provides different perspectives on this reaction:

NAMES: {full_info.get('Names', 'Not available')}

BOND CHANGES: {full_info.get('Bond Changes', 'Not available')}

FUNCTIONAL GROUPS: {full_info.get('Functional Groups', 'Not available')}

REACTION TYPE: {full_info.get('Reaction Classification', 'Not available')}
"""
        # Add dataset information to prompt - but only if important and keeping it concise
        if procedure_details:
            procedure_summary = procedure_details
            final_prompt += f"\nPROCEDURE: {procedure_summary}\n"
            
        if rxn_time or temperature or yield_info:
            conditions = []
            if temperature: conditions.append(f"Temperature: {temperature}")
            if rxn_time: conditions.append(f"Time: {rxn_time}")
            if yield_info: conditions.append(f"Yield: {yield_info}")
            final_prompt += f"\nCONDITIONS: {', '.join(conditions)}\n"
            
        if solvents or agents:
            materials = []
            if solvents: materials.append(f"Solvents: {', '.join(solvents)}")
            if agents: materials.append(f"Catalysts: {', '.join(agents)}")
            final_prompt += f"\nMATERIALS: {', '.join(materials)}\n"

        final_prompt += """
Provide a thorough, well-structured explanation of this reaction that:
1. Begins with a high-level summary of what type of reaction this is
2. Explains what happens at the molecular level (bonds broken/formed)
3. Discusses the functional group transformations
4. Includes specific experimental conditions (temperature, time, yield, solvents, catalysts)
5. Provides a detailed procedure summary if available
6. Mentions common applications or importance of this reaction type

Please give a complete and readable explanation of this reaction.
"""

        # Use a model instance with controlled output size
        focused_llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=2000)
        response = focused_llm.invoke(final_prompt)
        analysis_text = response.content.strip()
        
        # Structure the response data
        result = {
            'visualization_path': full_info.get('Visualization') if 'Visualization' in full_info and not full_info['Visualization'].startswith('Error') else None,
            'analysis': analysis_text,
            'reaction_classification': full_info.get('Reaction Classification', "No classification available"),
            'procedure_details': procedure_details,
            'reaction_time': rxn_time,
            'temperature': temperature,
            'yield': yield_info,
            'solvents': solvents if solvents else None,
            'agents_catalysts': agents if agents else None
        }
        
        # Cache the full response
        reaction_cache[reaction_smiles]['full_info'] = result
        
        return result
    
    except Exception as e:
        error_result = {
            'visualization_path': None,
            'analysis': f"Error in full analysis: {str(e)}"
        }
        reaction_cache[reaction_smiles]['full_info'] = error_result
        return error_result

# New function specifically for handling follow-up questions
def handle_followup_question(query, reaction_smiles):
    """
    Handle follow-up questions about a reaction using cached data 
    to avoid re-running expensive operations.
    """
    # Check if we have this reaction in our cache
    if reaction_smiles not in reaction_cache:
        # We haven't seen this reaction before, need to do full analysis
        return handle_full_info(query, reaction_smiles)
    
    print(f"Handling follow-up question for reaction: {reaction_smiles}")
    
    # Extract property-specific requests
    property_keywords = {
        'temperature': ['temperature', 'temp', 'heat', 'celsius', 'fahrenheit', 'kelvin', 'degrees', '°C', '°F', '°K'],
        'yield': ['yield', 'yields', 'percentage', 'amount', 'efficiency', 'how much'],
        'solvent': ['solvent', 'medium', 'solution', 'dissolved in'],
        'catalyst': ['catalyst', 'catalysts', 'catalytic', 'accelerator'],
        'time': ['time', 'duration', 'how long', 'minutes', 'hours', 'days'],
        'procedure': ['procedure', 'protocol', 'how to', 'steps', 'instructions', 'method'],
        'name': ['name', 'called', 'identity', 'what is it', 'what is this', 'called'],
        'functional_groups': ['functional group', 'substituent', 'functional', 'groups'],
        'bonds': ['bond', 'bonds', 'bonding', 'formed', 'broken'],
        'classification': ['type', 'class', 'classification', 'category', 'kind of reaction']
    }
    
    # Check for property-specific queries
    query_lower = query.lower()
    
    # First check if a visualization is requested
    if any(term in query_lower for term in ["visual", "picture", "image", "show", "draw", "representation"]):
        # Return cached visualization if available
        if 'visualization_path' in reaction_cache[reaction_smiles] and reaction_cache[reaction_smiles]['visualization_path']:
            return {
                "visualization_path": reaction_cache[reaction_smiles]['visualization_path'],
                "analysis": f"Visual representation of the reaction: {reaction_smiles}"
            }
        else:
            # Try to generate visualization now
            tool_dict = {tool.name.lower(): tool for tool in tools}
            visualizer_tool = tool_dict.get("chemvisualizer")
            
            if visualizer_tool:
                try:
                    viz_path = visualizer_tool.run(reaction_smiles)
                    reaction_cache[reaction_smiles]['visualization_path'] = viz_path
                    return {
                        "visualization_path": viz_path,
                        "analysis": f"Visual representation of the reaction: {reaction_smiles}"
                    }
                except Exception as e:
                    return {"visualization_path": None, "analysis": f"Error visualizing reaction: {str(e)}"}
            else:
                return {"visualization_path": None, "analysis": "ChemVisualizer tool not found"}
    
    # Prepare a prompt with ONLY the relevant information for the specific follow-up question
    prompt_parts = [f"User follow-up question about reaction {reaction_smiles}: {query}"]
    prompt_parts.append("\nHere's the relevant information from our analysis:")
    
    # Add ONLY the relevant information based on the query
    for prop_name, keywords in property_keywords.items():
        if any(keyword in query_lower for keyword in keywords):
            # Add only the relevant cached information for this property
            if prop_name == 'temperature' and 'temperature' in reaction_cache[reaction_smiles]:
                prompt_parts.append(f"\nTemperature: {reaction_cache[reaction_smiles]['temperature']}")
            
            elif prop_name == 'yield' and 'yield_info' in reaction_cache[reaction_smiles]:
                prompt_parts.append(f"\nYield: {reaction_cache[reaction_smiles]['yield_info']}")
            
            elif prop_name == 'time' and 'rxn_time' in reaction_cache[reaction_smiles]:
                prompt_parts.append(f"\nReaction Time: {reaction_cache[reaction_smiles]['rxn_time']}")
            
            elif prop_name == 'procedure' and 'procedure_details' in reaction_cache[reaction_smiles]:
                prompt_parts.append(f"\nProcedure Details: {reaction_cache[reaction_smiles]['procedure_details']}")
            
            elif prop_name == 'solvent' and 'solvents' in reaction_cache[reaction_smiles]:
                solvents = reaction_cache[reaction_smiles]['solvents']
                if solvents:
                    prompt_parts.append(f"\nSolvents: {', '.join(solvents)}")
            
            elif prop_name == 'catalyst' and 'agents' in reaction_cache[reaction_smiles]:
                agents = reaction_cache[reaction_smiles]['agents']
                if agents:
                    prompt_parts.append(f"\nCatalysts/Agents: {', '.join(agents)}")
            
            elif prop_name == 'name' and 'name_info' in reaction_cache[reaction_smiles]:
                prompt_parts.append(f"\nNames: {reaction_cache[reaction_smiles]['name_info']}")
            
            elif prop_name == 'functional_groups' and 'fg_info' in reaction_cache[reaction_smiles]:
                prompt_parts.append(f"\nFunctional Groups: {reaction_cache[reaction_smiles]['fg_info']}")
            
            elif prop_name == 'bonds' and 'bond_info' in reaction_cache[reaction_smiles]:
                prompt_parts.append(f"\nBond Changes: {reaction_cache[reaction_smiles]['bond_info']}")
            
            elif prop_name == 'classification' and 'classification_info' in reaction_cache[reaction_smiles]:
                prompt_parts.append(f"\nReaction Classification: {reaction_cache[reaction_smiles]['classification_info']}")
    
    # If no specific property found, provide a condensed summary
    if len(prompt_parts) <= 2:
        # Add minimal general information
        if 'classification_info' in reaction_cache[reaction_smiles]:
            prompt_parts.append(f"\nReaction Classification: {reaction_cache[reaction_smiles]['classification_info']}")
        
        if 'bond_info' in reaction_cache[reaction_smiles]:
            # Add a truncated version of bond info
            bond_info = reaction_cache[reaction_smiles]['bond_info']
            if len(bond_info) > 300:
                bond_info = bond_info[:300] + "... [truncated]"
            prompt_parts.append(f"\nKey Bond Changes: {bond_info}")
    
    prompt_parts.append("\nPlease provide a focused answer to the user's question based on this information. Keep your response concise and relevant.")
    
    # Create a focused prompt with only the relevant information
    final_prompt = "\n".join(prompt_parts)
    
    # Use a model with controlled output size
    focused_llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=1000)
    response = focused_llm.invoke(final_prompt)
    
    return {
        "visualization_path": None,
        "analysis": response.content.strip()
    }
def enhanced_query(query, callbacks=None):
    """
    Enhanced query handler that routes user queries to the appropriate handler
    based on the type of query and context.
    
    Args:
        query (str): The user's query text
        callbacks (Optional): Callback objects for the agent
        
    Returns:
        dict: Results with visualization_path and analysis text
    """
    try:
        # Extract reaction SMILES from the query using the existing function
        reaction_smiles = extract_reaction_smiles(query)
        
        # If no reaction SMILES found, use the regular agent
        if not reaction_smiles:
            print("No reaction SMILES found in query, using regular agent")
            result = agent.invoke({"input": query}, {"callbacks": callbacks} if callbacks else {})
            return {"visualization_path": None, "analysis": extract_final_answer(result.get("output", ""))}
        
        # Check if this is a reaction we've seen before (follow-up question)
        is_followup = reaction_smiles in reaction_cache
        
        # Check if it's a full analysis request
        if "full" in query.lower() and any(term in query.lower() for term in 
                                          ["information", "analysis", "detail", "explain", "tell me about"]):
            print(f"Full analysis requested for reaction: {reaction_smiles}")
            return handle_full_info(query, reaction_smiles)
        
        # Check if it's just a visualization request
        elif any(term in query.lower() for term in ["visual", "picture", "image", "show", "draw", "representation", "diagram"]):
            print(f"Visualization requested for reaction: {reaction_smiles}")
            
            # Check if we already have visualization in cache
            if reaction_smiles in reaction_cache and 'visualization_path' in reaction_cache[reaction_smiles]:
                viz_path = reaction_cache[reaction_smiles]['visualization_path']
                if viz_path and not viz_path.startswith("Error"):
                    return {
                        "visualization_path": viz_path,
                        "analysis": f"Visual representation of the reaction: {reaction_smiles}"
                    }
            
            # If not cached, generate visualization
            tool_dict = {tool.name.lower(): tool for tool in tools}
            visualizer_tool = tool_dict.get("chemvisualizer")
            
            if visualizer_tool:
                try:
                    viz_path = visualizer_tool.run(reaction_smiles)
                    # Cache the result
                    if reaction_smiles not in reaction_cache:
                        reaction_cache[reaction_smiles] = {}
                    reaction_cache[reaction_smiles]['visualization_path'] = viz_path
                    
                    return {
                        "visualization_path": viz_path,
                        "analysis": f"Visual representation of the reaction: {reaction_smiles}"
                    }
                except Exception as e:
                    return {"visualization_path": None, "analysis": f"Error visualizing reaction: {str(e)}"}
            else:
                return {"visualization_path": None, "analysis": "ChemVisualizer tool not found"}
        
        # Check if it's a follow-up question about a specific reaction property
        elif is_followup and any(keyword in query.lower() for keyword in [
            "temperature", "yield", "solvent", "catalyst", "time", "pressure", "ph", "mechanism", 
            "applications", "procedure", "protocol", "uses", "purpose", "bonds", "functional", 
            "classification", "name", "what is", "how does", "explain", "why"
        ]):
            print(f"Follow-up question about reaction: {reaction_smiles}")
            return handle_followup_question(query, reaction_smiles)
        
        # For new reactions with specific property questions
        elif any(keyword in query.lower() for keyword in [
            "temperature", "yield", "solvent", "catalyst", "time", "pressure", "ph", "mechanism", 
            "applications", "procedure", "protocol", "uses", "purpose", "bonds", "functional groups",
            "classification", "name", "what type", "how does", "explain", "why"
        ]):
            print(f"Specific property question about reaction: {reaction_smiles}")
            
            # Check what tools we need based on the query
            if "bond" in query.lower() or "break" in query.lower() or "form" in query.lower():
                # Use bond analyzer only
                tool_dict = {tool.name.lower(): tool for tool in tools}
                bond_tool = tool_dict.get("bondchangeanalyzer")
                
                if bond_tool:
                    try:
                        bond_result = bond_tool.run(reaction_smiles)
                        # Cache the result
                        if reaction_smiles not in reaction_cache:
                            reaction_cache[reaction_smiles] = {}
                        reaction_cache[reaction_smiles]['bond_info'] = bond_result
                        
                        return {"visualization_path": None, "analysis": bond_result}
                    except Exception as e:
                        return {"visualization_path": None, "analysis": f"Error analyzing bond changes: {str(e)}"}
                
            elif "functi" in query.lower() or "group" in query.lower():
                # Use functional groups analyzer only
                tool_dict = {tool.name.lower(): tool for tool in tools}
                fg_tool = tool_dict.get("funcgroups")
                
                if fg_tool:
                    try:
                        fg_result = fg_tool.run(reaction_smiles)
                        # Cache the result
                        if reaction_smiles not in reaction_cache:
                            reaction_cache[reaction_smiles] = {}
                        reaction_cache[reaction_smiles]['fg_info'] = fg_result
                        
                        return {"visualization_path": None, "analysis": fg_result}
                    except Exception as e:
                        return {"visualization_path": None, "analysis": f"Error analyzing functional groups: {str(e)}"}
                
            elif "name" in query.lower() or "call" in query.lower():
                # Use name tool only
                tool_dict = {tool.name.lower(): tool for tool in tools}
                name_tool = tool_dict.get("smiles2name")
                
                if name_tool:
                    try:
                        name_result = name_tool.run(reaction_smiles)
                        # Cache the result
                        if reaction_smiles not in reaction_cache:
                            reaction_cache[reaction_smiles] = {}
                        reaction_cache[reaction_smiles]['name_info'] = name_result
                        
                        return {"visualization_path": None, "analysis": name_result}
                    except Exception as e:
                        return {"visualization_path": None, "analysis": f"Error finding names: {str(e)}"}
                
            elif "classif" in query.lower() or "type" in query.lower() or "what kind" in query.lower():
                # Use reaction classifier only
                if reaction_classifier:
                    try:
                        classifier_result = reaction_classifier._run(reaction_smiles)
                        # Cache the result
                        if reaction_smiles not in reaction_cache:
                            reaction_cache[reaction_smiles] = {}
                        reaction_cache[reaction_smiles]['classification_info'] = classifier_result
                        
                        return {"visualization_path": None, "analysis": classifier_result}
                    except Exception as e:
                        return {"visualization_path": None, "analysis": f"Error classifying reaction: {str(e)}"}
            
            # For more complex or ambiguous property queries, run a full analysis
            print(f"Running full analysis for complex property query on reaction: {reaction_smiles}")
            result = handle_full_info(query, reaction_smiles)
            
            # Filter the result to focus on what was asked
            focused_llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=1000)
            prompt = f"""
            The user asked: {query}
            
            Here's our full analysis of the reaction {reaction_smiles}:
            
            {result.get('analysis', 'No analysis available')}
            
            Please provide a focused answer to the user's specific question, extracting only the relevant information
            from our analysis. Keep your response concise and directly address what they asked about.
            """
            
            focused_response = focused_llm.invoke(prompt)
            result['analysis'] = focused_response.content.strip()
            
            return result
        
        # Default to full analysis for any reaction SMILES
        else:
            print(f"Using full analysis as default for reaction: {reaction_smiles}")
            return handle_full_info(query, reaction_smiles)

    except Exception as e:
        import traceback
        stack_trace = traceback.format_exc()
        print(f"Error in enhanced_query: {str(e)}\n{stack_trace}")
        return {"visualization_path": None, "analysis": f"Error occurred: {str(e)}"}















#CURRENT 
# import os
# import re
# import requests
# import api_config
# from tools.make_tools import make_tools
# from langchain.agents import AgentExecutor, ZeroShotAgent
# from langchain.prompts import PromptTemplate
# from langchain_openai import ChatOpenAI
# import pandas as pd
# from tools.asckos import ReactionClassifier
# import json
# from functools import lru_cache
# import traceback # Added for better error reporting

# # Setup LLM and tools
# llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=4000)
# tools = make_tools(llm=llm)

# # Initialize reaction classifier with dataset paths
# dataset_path1 = os.environ.get('REACTION_DATASET_PATH1', None)
# dataset_path2 = os.environ.get('REACTION_DATASET_PATH2', None)
# reaction_classifier = ReactionClassifier(dataset_path1, dataset_path2)

# # Create a global cache for reaction data (used by handle_full_info)
# reaction_cache = {}

# # Prompt parts (Keep these as they define agent behavior)
# PREFIX = """
# You are Chem Copilot, an expert chemistry assistant. You have access to the following tools to analyze molecules and chemical reactions.

# Always begin by understanding the user's **intent** — what kind of information are they asking for?

# Here is how to choose tools:

# - If the user gives a SMILES or reaction SMILES and asks for the name, you MUST use **SMILES2Name**. Do NOT analyze bonds or functional groups for this task.
# - Use **NameToSMILES**: when the user gives a compound/reaction name and wants the SMILES or structure.
# - Use **FuncGroups**: when the user wants to analyze functional groups in a molecule or a reaction (input is SMILES or reaction SMILES).
# - Use **BondChangeAnalyzer**: when the user asks for which bonds are broken, formed, or changed in a chemical reaction.
# - Use **ChemVisualizer**: when the user asks to visualize, show, draw, or see a picture/image/representation of a molecule or reaction (input is SMILES or reaction SMILES).
# - Use **ReactionClassifier**: when the user asks about the type, class, mechanism, conditions, procedure, or properties (like temperature, yield, solvent, catalyst, time) of a reaction. Input should be the reaction SMILES and the specific property/aspect asked about (e.g., "reaction_smiles query=temperature"). If only the reaction SMILES is given with a general query about the reaction, provide a general classification and summary.

# If the user wants all of the above (full analysis), respond with "This requires full analysis." (This will be handled by a separate function.)

# Always return your answer in this format:
# Final Answer: <your answer here>

# For **FuncGroups** results:
# - Always list the functional groups identified in each reactant and product separately
# - Include the transformation summary showing disappeared groups, appeared groups, and unchanged groups
# - Provide a clear conclusion about what transformation occurred in the reaction

# For **BondChangeAnalyzer** results:
# - Always list the specific bonds that were broken, formed, or changed with their bond types
# - Include the atom types involved in each bond (e.g., C-O, N-H)
# - Provide a clear conclusion summarizing the key bond changes in the reaction

# For **ChemVisualizer** results:
# - State clearly that the visualization was generated and provide the path: "Visualization saved to: <path>"

# For **ReactionClassifier** results:
# - Provide the requested information (classification, property, conditions, etc.) based on the query.
# """

# FORMAT_INSTRUCTIONS = """
# You can only respond with a single complete
# "Thought, Action, Action Input" format
# OR a single "Final Answer" format

# Complete format:

# Thought: (reflect on your progress and decide what to do next)
# Action: (the action name, should be one of [{tool_names}])
# Action Input: (the input string to the action)

# OR

# Final Answer: (the final answer to the original input question)
# """

# SUFFIX = """
# Question: {input}
# {agent_scratchpad}
# """

# prompt = ZeroShotAgent.create_prompt(
#     tools=tools,
#     prefix=PREFIX,
#     suffix=SUFFIX,
#     format_instructions=FORMAT_INSTRUCTIONS,
#     input_variables=["input", "agent_scratchpad"]
# )

# agent_chain = ZeroShotAgent.from_llm_and_tools(llm=llm, tools=tools, prompt=prompt)
# # Ensure verbose=True is set here for the detailed output
# agent = AgentExecutor(agent=agent_chain, tools=tools, verbose=True, handle_parsing_errors=True) # Added handle_parsing_errors

# def extract_final_answer(full_output: str):
#     # This regex is more robust for extracting the final answer, even with preceding text
#     match = re.search(r"Final Answer:\s*(.*)", full_output, re.DOTALL | re.MULTILINE)
#     if match:
#         return match.group(1).strip()
#     # Fallback if the specific "Final Answer:" marker isn't found (e.g., error)
#     # Look for the last block of text after the final "Observation:" or if no observation, the whole output
#     parts = re.split(r'\nObservation:', full_output)
#     if len(parts) > 1:
#         # Get text after the last observation
#         last_part = parts[-1].strip()
#         # Remove potential thought remnants
#         thought_match = re.search(r"Thought:\s*(.*)", last_part, re.DOTALL | re.IGNORECASE)
#         if thought_match:
#             return thought_match.group(1).strip() # Return thought if it looks like the final step
#         return last_part # Return everything after last observation
#     else:
#         # If no observation, it might be a direct final answer or an error message
#         return full_output.strip()


# # Function to query dataset - KEEP AS IS (caching is good)
# @lru_cache(maxsize=100)
# def query_reaction_dataset(reaction_smiles):
#     """
#     Query the dataset for specific information about a reaction
#     based on the reaction SMILES string. Caches results.
#     """
#     # Check if we have this in our reaction_cache first (optional optimization)
#     # if reaction_smiles in reaction_cache and 'dataset_info' in reaction_cache[reaction_smiles]:
#     #     return reaction_cache[reaction_smiles]['dataset_info'] # Can uncomment if needed

#     try:
#         df = None # Initialize df
#         if hasattr(reaction_classifier, 'dataset1') and reaction_classifier.dataset1 is not None:
#             df = reaction_classifier.dataset1
#         elif hasattr(reaction_classifier, 'dataset2') and reaction_classifier.dataset2 is not None:
#             df = reaction_classifier.dataset2
#         else:
#             print("Warning: Reaction classifier datasets not loaded.")
#             return None # Return None if no dataset is available

#         if df is None or df.empty:
#             # print(f"Debug: Dataset is None or empty for {reaction_smiles}")
#             return None

#         fields_to_extract = [
#             'procedure_details', 'rxn_time', 'temperature', 'yield_000',
#             'reaction_name', 'reaction_classname', 'prediction_certainty'
#         ]
#         smiles_columns = ['rxn_str', 'reaction_smiles', 'smiles', 'rxn_smiles']
#         exact_match = None

#         for col in smiles_columns:
#             if col in df.columns:
#                 # Ensure comparison handles potential dtype issues
#                 try:
#                     # print(f"Debug: Checking column '{col}' for SMILES '{reaction_smiles}'")
#                     temp_match = df[df[col].astype(str) == str(reaction_smiles)]
#                     if not temp_match.empty:
#                         # print(f"Debug: Found match in column '{col}'")
#                         exact_match = temp_match
#                         break
#                 except KeyError:
#                     # print(f"Debug: Column '{col}' not found.")
#                     continue # Skip if column doesn't exist
#                 except Exception as e:
#                     print(f"Debug: Error comparing in column '{col}': {e}") # Log other errors
#                     continue

#         if exact_match is not None and not exact_match.empty:
#             row = exact_match.iloc[0]
#             result = {}
#             for field in fields_to_extract:
#                 if field in row.index and pd.notna(row[field]): # Use pandas check for NaN
#                     result[field] = str(row[field])

#             solvent_count = 0
#             for i in range(11):
#                 solvent_key = f'solvent_{i:03d}'
#                 if solvent_key in row.index and pd.notna(row[solvent_key]):
#                     result[solvent_key] = str(row[solvent_key])
#                     solvent_count += 1
#                     if solvent_count >= 3: break

#             agent_count = 0
#             for i in range(16):
#                 agent_key = f'agent_{i:03d}'
#                 if agent_key in row.index and pd.notna(row[agent_key]):
#                     result[agent_key] = str(row[agent_key])
#                     agent_count += 1
#                     if agent_count >= 3: break

#             # Optional: Update reaction_cache if using it here
#             # if reaction_smiles not in reaction_cache: reaction_cache[reaction_smiles] = {}
#             # reaction_cache[reaction_smiles]['dataset_info'] = result
#             # print(f"Debug: Found dataset info for {reaction_smiles}: {result}")
#             return result
#         else:
#             # print(f"Debug: No exact match found for {reaction_smiles}")
#             # Optional: Update reaction_cache if using it here
#             # if reaction_smiles not in reaction_cache: reaction_cache[reaction_smiles] = {}
#             # reaction_cache[reaction_smiles]['dataset_info'] = None
#             return None

#     except Exception as e:
#         print(f"Error querying reaction dataset for '{reaction_smiles}': {e}")
#         # print(traceback.format_exc()) # Uncomment for detailed traceback
#         return None


# # Extract SMILES from user query - KEEP AS IS
# def extract_reaction_smiles(query):
#     """
#     Extract reaction SMILES from a user query using various patterns.
#     Returns None if no valid SMILES is found.
#     """
#     # Pattern for reaction SMILES (reactants>>products) - more robust
#     smiles_pattern = r"([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}\<\>]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}\<\>]*(?:\.[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}\<\>]+)*)"
#     match = re.search(smiles_pattern, query)
#     if match: return match.group(1)

#     # Try to extract product SMILES if not a full reaction (less common)
#     product_pattern = r"rxn\s+([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+)"
#     product_match = re.search(product_pattern, query)
#     if product_match: return f"[R]>>{product_match.group(1)}" # Represent as generic reaction

#     # Pattern for reaction SMILES with explicit labeling
#     labeled_pattern = r"reaction\s+(?:SMILES)?[:=]?\s*([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]*(?:\.[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}\<\>]+)*)"
#     labeled_match = re.search(labeled_pattern, query, re.IGNORECASE)
#     if labeled_match: return labeled_match.group(1)

#     # Alternative pattern for reaction descriptions
#     alt_pattern = r"(?:for|about)\s+the\s+reaction\s+([A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]+>>[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}]*(?:\.[A-Za-z0-9@\[\]\.\+\-\=\#\:\(\)\\\/;\$\%\|\{\}\<\>]+)*)"
#     alt_match = re.search(alt_pattern, query, re.IGNORECASE)
#     if alt_match: return alt_match.group(1)

#     return None


# # handle_full_info - KEEP AS IS (called only for explicit "full analysis" requests)
# def handle_full_info(query, reaction_smiles):
#     """
#     Process a full analysis request for a reaction.
#     Aggregates results from multiple tools and dataset info.
#     Uses caching internally where appropriate (e.g., query_reaction_dataset).
#     """
#     print(f"Running full analysis for reaction: {reaction_smiles}\n")
#     # Check cache only for the final synthesized result if desired (optional)
#     # if reaction_smiles in reaction_cache and 'full_info_result' in reaction_cache[reaction_smiles]:
#     #     print("Using cached full analysis result.")
#     #     return reaction_cache[reaction_smiles]['full_info_result']

#     full_info = {}
#     tool_dict = {tool.name.lower(): tool for tool in tools}

#     # --- Run Tools ---
#     # Visualization
#     visualizer_tool = tool_dict.get("chemvisualizer")
#     if visualizer_tool:
#         try: full_info['Visualization'] = visualizer_tool.run(reaction_smiles)
#         except Exception as e: full_info['Visualization'] = f"Error visualizing: {str(e)}"
#     else: full_info['Visualization'] = "ChemVisualizer tool not found"

#     # Names
#     name_tool = tool_dict.get("smiles2name")
#     if name_tool:
#         try:
#             name_result = name_tool.run(reaction_smiles)
#             full_info['Names'] = name_result[:1000] + ("..." if len(name_result)>1000 else "") # Truncate
#         except Exception as e: full_info['Names'] = f"Error getting names: {str(e)}"
#     else: full_info['Names'] = "SMILES2Name tool not found"

#     # Functional Groups
#     fg_tool = tool_dict.get("funcgroups")
#     if fg_tool:
#         try:
#             fg_result = fg_tool.run(reaction_smiles)
#             full_info['Functional Groups'] = fg_result[:2000] + ("..." if len(fg_result)>2000 else "") # Truncate
#         except Exception as e: full_info['Functional Groups'] = f"Error analyzing FGs: {str(e)}"
#     else: full_info['Functional Groups'] = "FuncGroups tool not found"

#     # Bond Changes
#     bond_tool = tool_dict.get("bondchangeanalyzer")
#     if bond_tool:
#         try:
#             bond_result = bond_tool.run(reaction_smiles)
#             full_info['Bond Changes'] = bond_result[:2000] + ("..." if len(bond_result)>2000 else "") # Truncate
#         except Exception as e: full_info['Bond Changes'] = f"Error analyzing bonds: {str(e)}"
#     else: full_info['Bond Changes'] = "BondChangeAnalyzer tool not found"

#     # Reaction Classification (using the tool/classifier directly here)
#     if reaction_classifier:
#         try:
#             # Assuming _run gives a general classification string
#             classifier_result = reaction_classifier._run(reaction_smiles)
#             summary_match = re.search(r'## Summary\n(.*?)$', classifier_result, re.DOTALL | re.IGNORECASE)
#             if summary_match:
#                 full_info['Reaction Classification'] = summary_match.group(1).strip()[:1000] + "..." # Truncate
#             else:
#                 full_info['Reaction Classification'] = classifier_result[:1000] + "..." # Truncate raw output
#         except Exception as e: full_info['Reaction Classification'] = f"Error classifying: {str(e)}"
#     else: full_info['Reaction Classification'] = "ReactionClassifier tool not found"

#     # --- Query Dataset ---
#     dataset_info = query_reaction_dataset(reaction_smiles)
#     procedure_details = None
#     rxn_time = None
#     temperature = None
#     yield_info = None
#     solvents = []
#     agents = []

#     if dataset_info:
#         procedure_details = dataset_info.get('procedure_details', None)
#         if procedure_details and len(procedure_details) > 500:
#              procedure_details = procedure_details[:500] + "...[truncated]"
#         rxn_time = dataset_info.get('rxn_time', None)
#         temperature = dataset_info.get('temperature', None)
#         yield_info = dataset_info.get('yield_000', None)
#         solvents = [v for k, v in dataset_info.items() if k.startswith('solvent_')]
#         agents = [v for k, v in dataset_info.items() if k.startswith('agent_')]

#     # --- Synthesize with LLM ---
#     final_prompt = f"""Synthesize this reaction analysis into a comprehensive explanation:
# Reaction SMILES: {reaction_smiles}

# Analysis Details:
# Names: {full_info.get('Names', 'N/A')}
# Bond Changes: {full_info.get('Bond Changes', 'N/A')}
# Functional Groups: {full_info.get('Functional Groups', 'N/A')}
# Reaction Classification: {full_info.get('Reaction Classification', 'N/A')}
# """
#     if procedure_details: final_prompt += f"\nProcedure Summary: {procedure_details}\n"
#     conditions = []
#     if temperature: conditions.append(f"Temp: {temperature}")
#     if rxn_time: conditions.append(f"Time: {rxn_time}")
#     if yield_info: conditions.append(f"Yield: {yield_info}%") # Assuming yield is percentage
#     if conditions: final_prompt += f"\nConditions: {', '.join(conditions)}\n"
#     if solvents: final_prompt += f"\nSolvents: {', '.join(solvents[:3])}\n" # Limit display
#     if agents: final_prompt += f"\nCatalysts/Agents: {', '.join(agents[:3])}\n" # Limit display

#     final_prompt += """
# Provide a thorough, well-structured explanation of this reaction that:

# 1. Begins with a high-level summary of what type of reaction this is
# 2. Explains what happens at the molecular level (bonds broken/formed)
# 3. Discusses the functional group transformations
# 4. Includes specific experimental conditions (temperature, time, yield, solvents, catalysts)
# 5. Provides a detailed procedure summary if available
# 6. Mentions common applications or importance of this reaction type

# Please give a complete and readable explanation of this reaction.
# """
#     try:
#         focused_llm = ChatOpenAI(model="gpt-4o", temperature=0, max_tokens=2000)
#         response = focused_llm.invoke(final_prompt)
#         analysis_text = response.content.strip()
#     except Exception as e:
#         analysis_text = f"Error synthesizing analysis: {str(e)}\n\nRaw Data:\n{final_prompt}" # Provide raw data on error

#     # Structure the final result
#     result = {
#         'visualization_path': full_info.get('Visualization') if isinstance(full_info.get('Visualization'), str) and not full_info['Visualization'].startswith('Error') else None,
#         'analysis': analysis_text,
#         # Include raw components for potential UI display or debugging
#         'raw_classification': full_info.get('Reaction Classification'),
#         'raw_procedure': procedure_details,
#         'raw_time': rxn_time,
#         'raw_temp': temperature,
#         'raw_yield': yield_info,
#         'raw_solvents': solvents,
#         'raw_agents': agents
#     }

#     # Optional: Cache the final synthesized result
#     # if reaction_smiles not in reaction_cache: reaction_cache[reaction_smiles] = {}
#     # reaction_cache[reaction_smiles]['full_info_result'] = result

#     return result

# # REMOVED handle_followup_question - rely on agent scratchpad/context

# # MODIFIED enhanced_query function
# def enhanced_query(query, callbacks=None):
#     """
#     Enhanced query handler. Routes explicit "full analysis" requests to a
#     dedicated handler. Otherwise, uses the standard AgentExecutor to get
#     verbose step-by-step output.
#     """
#     try:
#         # --- Special Case: Explicit Full Analysis Request ---
#         if "full information" in query.lower() or "full analysis" in query.lower():
#             reaction_smiles = extract_reaction_smiles(query)
#             if reaction_smiles:
#                 print(f"*** Full analysis requested. Using handle_full_info for: {reaction_smiles} ***")
#                 # Call the dedicated handler which synthesizes results from multiple tools
#                 # NOTE: This call itself will NOT produce the agent's Thought/Action trace
#                 full_analysis_result = handle_full_info(query, reaction_smiles)
#                 return full_analysis_result # Return the dict from handle_full_info
#             else:
#                 # If "full analysis" requested but no SMILES found, let the agent try
#                 print("Full analysis requested, but no SMILES found. Passing to standard agent.")
#                 # Fall through to the agent invoke below

#         # --- Default Case: Use the Standard AgentExecutor ---
#         # For all other queries, let the agent handle it.
#         # This WILL produce the verbose "Thought, Action, Action Input, Observation" trace
#         # because agent was initialized with verbose=True.
#         print(f"--- Passing query to standard agent ---")
#         result = agent.invoke({"input": query}, {"callbacks": callbacks} if callbacks else {})

#         # Extract the final answer from the potentially verbose agent output
#         final_answer = extract_final_answer(result.get("output", "Agent did not return output."))

#         # Attempt to extract visualization path from the agent's final answer
#         # This relies on ChemVisualizer tool outputting a predictable string
#         visualization_path = None
#         viz_match = re.search(r"(?:Visualization saved to|Image generated at|Path:)\s*([^\s\'\"]+)", final_answer, re.IGNORECASE)
#         if viz_match:
#             visualization_path = viz_match.group(1).strip()
#             print(f"Extracted visualization path: {visualization_path}") # Debug print

#         return {"visualization_path": visualization_path, "analysis": final_answer}

#     except Exception as e:
#         print(f"!!! Error in enhanced_query: {str(e)} !!!")
#         print(traceback.format_exc()) # Print full traceback for debugging
#         return {"visualization_path": None, "analysis": f"An error occurred processing your request: {str(e)}"}

# # Example Run (if needed for testing)
# # if __name__ == "__main__":
# #     # Ensure environment variables are set (e.g., OPENAI_API_KEY, dataset paths)
# #     # Make sure dataset files actually exist at the specified paths
# #     print(f"Dataset Path 1: {dataset_path1}")
# #     print(f"Dataset Path 2: {dataset_path2}")
# #     print("-----")

# #     # Test 1: Query that should use the agent (verbose output expected)
# #     # test_query_agent = "What is the temperature used for Weinreb iodo coupling reaction C(N(C(C)C)C)(N=C(N)=O)=O.C[Si](C)(C)I.O=C=O>>C(IF)(F)(F)Cl.O=C1C=C(O)C=CC1"
# #     # Example with a SMILES known to be in USPTO-MIT dataset
# #     test_query_agent = "What solvents are present in 2,5-Pyrroledione synthesis? C1(C(=O)CC(=O)N1)=CC=2C=CC=CC=2.O=[N+]([O-])C=3C=CC(Cl)=CC=3>>C1(C(=O)C=C(C(=O)N1)C=2C=CC=CC=2)C=3C=CC(Cl)=CC=3"
# #     print(f"\nTesting Agent Query: {test_query_agent}")
# #     result_agent = enhanced_query(test_query_agent)
# #     print("\nAgent Query Result:")
# #     print(json.dumps(result_agent, indent=2))
# #     print("-----")

# #     # Test 2: Explicit Full Analysis (should call handle_full_info, less verbose output expected during call)
# #     # test_query_full = "Give full analysis of the reaction C(N(C(C)C)C)(N=C(N)=O)=O.C[Si](C)(C)I.O=C=O>>C(IF)(F)(F)Cl.O=C1C=C(O)C=CC1"
# #     test_query_full = "Full information for reaction C1(C(=O)CC(=O)N1)=CC=2C=CC=CC=2.O=[N+]([O-])C=3C=CC(Cl)=CC=3>>C1(C(=O)C=C(C(=O)N1)C=2C=CC=CC=2)C=3C=CC(Cl)=CC=3"
# #     print(f"\nTesting Full Analysis Query: {test_query_full}")
# #     result_full = enhanced_query(test_query_full)
# #     print("\nFull Analysis Result:")
# #     print(json.dumps(result_full, indent=2))
# #     print("-----")

# #     # Test 3: Query without SMILES (should use agent)
# #     test_query_no_smiles = "What is a Grignard reaction?"
# #     print(f"\nTesting No-SMILES Query: {test_query_no_smiles}")
# #     result_no_smiles = enhanced_query(test_query_no_smiles)
# #     print("\nNo-SMILES Query Result:")
# #     print(json.dumps(result_no_smiles, indent=2))
# #     print("-----")

# #     # Test 4: Visualization request (should use agent)
# #     test_query_viz = "Draw the reaction C1(C(=O)CC(=O)N1)=CC=2C=CC=CC=2.O=[N+]([O-])C=3C=CC(Cl)=CC=3>>C1(C(=O)C=C(C(=O)N1)C=2C=CC=CC=2)C=3C=CC(Cl)=CC=3"
# #     print(f"\nTesting Visualization Query: {test_query_viz}")
# #     result_viz = enhanced_query(test_query_viz)
# #     print("\nVisualization Query Result:")
# #     print(json.dumps(result_viz, indent=2))
# #     print("-----")

# #     # Test 5: Full analysis without SMILES (should use agent)
# #     test_query_full_no_smiles = "Give me full analysis"
# #     print(f"\nTesting Full Analysis No-SMILES Query: {test_query_full_no_smiles}")
# #     result_full_no_smiles = enhanced_query(test_query_full_no_smiles)
# #     print("\nFull Analysis No-SMILES Query Result:")
# #     print(json.dumps(result_full_no_smiles, indent=2))
# #     print("-----")