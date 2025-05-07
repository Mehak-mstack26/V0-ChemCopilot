# Carefully analyze the provided text, using images as a supplementary reference. If there is a conflict between the image and text, prioritize the information from the text.
# Ensure that the extracted reactions are interconnected: each subsequent reaction's reactants must include at least one product from the previous reaction. Continue this sequence until all reactants in the first reaction are common laboratory or commercially available substances.

# reaction_prompt = """
# Extract all distinct chemical reactions mentioned within. Include each reaction only once, if it appears more than once.
# Only extract reactions where all reactants and products are fully identified. Exclude any reaction if any reactant or product is unspecified.
# Ensure that there are no identical substances among the reactants, products, catalysts, and solvents. If any identical substances are found, reconsider the validity of the reaction.
#
# Ensure to unify the substance name throughout.
# Be meticulous not to alter the chemical names in a way that changes the identity of the substances.
# Unify Nomenclature: For example, if you encounter "poly(4-acetylstyrene)" and "poly(4-acetyl styrene)" in the text, recognize that they refer to the same substance. Always use the unified name, "poly(4-acetylstyrene)," in your output.
# Similarly, unify "poly(4-vinylphenol)," "poly(4-hydroxystyrene)," and "Polyvinylphenol" to a single name like "Polyhydroxystyrene."
#
# Format the output strictly as follows:
#
# Reaction 001:
# Reactants: List the IUPAC nomenclatures, separated by commas.
# Products: List the IUPAC nomenclatures, separated by commas.
# Conditions: List the following in the exact order, separated by commas, skip any condition that is not provided or if it is unknown::
# - If a synthesis method is provided, provide a professional term used to describe the reaction
# - If a catalyst is provided, prefix with 'Catalyst: ' followed by the IUPAC nomenclatures.
# - If a solvent is provided, prefix with 'Solvent: ' followed by the IUPAC nomenclatures.
# - If an atmosphere is provided, prefix with 'Atmosphere: ' followed by the specified gas.
# - For temperature, provide specific value or range with °C without any prefix.
# - For pressure, provide specific value or range with atm or bar without any prefix.
# - For duration, provide specific value or range with h or min or d without any prefix.
# - For yield, provide specific value or range in % without any prefix.
#
# Do not include any explanatory notes, brackets, or additional information.
# """
#
#
# How can common laboratory and commercial chemical compounds be used to synthesize "{substance}" in one or more steps? Please provide all reactions based on the given content.
#
# - Confirm that the extracted reactions are interconnected: the reactants in each subsequent reaction must include at least one product from the previous reaction.
#
#
# reaction_prompt_v4 = """
# How can common laboratory and commercial chemical compounds be used to synthesize "{substance}" in one or more steps? Please provide all reactions based on the given content.
#
# You must strictly adhere to the following rules:
# - Note that the given content contains multiple synthesis reactions for "{substance}". Ensure that all reactions are extracted without omission.
# - Confirm that the extracted reactions are interconnected: the reactants in each subsequent reaction must include at least one product from the previous reaction.
# - All names must strictly match the provided content without modification or interpretation. Maintain consistent naming of the substance throughout.
# - Do not use general terms such as "Diamine" or "Dianhydride" to represent a class of substances;
# - Do not use abbreviations for any substance names;
# - Do not include any explanatory notes, brackets, or additional information.
# - All names must be provided in their full specific forms.
# - Be careful to ensure that solvents, catalysts, or other reaction conditions are not mistakenly listed as reactants. Reactants must only include compounds directly involved in the chemical transformation.
#
# Format the output strictly as follows:
#
# Initial Output:
#
# Reaction 001:
# Reactants: List the substances with specific names, separated by commas.
# Products: List the substances with specific names, separated by commas.
# Conditions: List the following in the exact order, separated by commas, skipping any condition that is not provided or is unknown:
# - If a catalyst is provided, prefix with "Catalyst: " followed by the substances with specific names.
# - If a solvent is provided, prefix with "Solvent: " followed by the substances with specific names.
# - If other non-reactive reagents are provided, such as drying agents, stabilizers, adsorbents, buffers, etc., prefix with their respective names (e.g., 'Drying agents: ', 'Stabilizers: ', 'Adsorbents: ', 'Buffers: ') followed by their specific types.
# - If an atmosphere is provided, prefix with "Atmosphere: " followed by the specified gas in full names.
# - For temperature, provide the specific value or range with °C without any prefix.
# - For pressure, provide the specific value or range with atm or bar without any prefix.
# - For duration, provide the specific value or range with h, min, or d without any prefix.
# - For yield, provide the specific value or range in % without any prefix.
#
# Check whether any substance in the Reactants is commonly used as a solvent, catalyst, drying agent, stabilizer, adsorbent, buffer, etc. , move it to the conditions section
# Check whether any general terms such as "Diamine" or "Dianhydride" are identified, recheck the content and modify it.
# Document your checking process. Finally, re-output all the corrected reactions.
#
# Checking process:
# - Reaction 001:
#
# Final Output:
# """


# prompts.py

# --- Reaction Extraction Prompts ---

prompt_reaction_extraction_cot = """
Extract all distinct chemical reactions from given content. Include each reaction only once, if it appears more than once.
Only extract reactions where all reactants and products are fully identified. Exclude any reaction if any reactant or product is unspecified.

You must strictly adhere to the following rules:
Confirm that the extracted reactions are interconnected: the reactants in each subsequent reaction must include at least one product from the previous reaction.
Be careful to ensure that solvents, catalysts, or other reaction conditions are not mistakenly listed as reactants. Reactants must only include compounds directly involved in the chemical transformation.
All names must be provided in their full specific forms:
- Maintain consistent naming of each substance throughout.
- Do not use general terms, such as "Diamine" or "Dianhydride," to represent a class of substances for any reactant substance name.
- Do not use abbreviations for any substance name.
- Do not include explanatory notes, brackets, or additional information.

Format the output strictly as follows:

Initial Output:

Reaction 001:
Reactants: List the substances with specific names, separated by commas.
Products: List the substances with specific names, separated by commas.
Conditions: List the following in the exact order, separated by commas, skipping any condition that is not provided or is unknown:
- If a catalyst is provided, prefix with "Catalyst: " followed by the substances with specific names.
- If a solvent is provided, prefix with "Solvent: " followed by the substances with specific names.
- If other non-reactive reagents are provided, such as drying agents, stabilizers, adsorbents, buffers, etc., prefix with their respective names (e.g., 'Drying agents: ', 'Stabilizers: ', 'Adsorbents: ', 'Buffers: ') followed by their specific types.
- If an atmosphere is provided, prefix with "Atmosphere: " followed by the specified gas in full names.
- For temperature, provide the specific value or range with °C without any prefix.
- For pressure, provide the specific value or range with atm or bar without any prefix.
- For duration, provide the specific value or range with h, min, or d without any prefix.
- For yield, provide the specific value or range in % without any prefix.

Check whether any substance in the Reactants is commonly used as a solvent, catalyst, drying agent, stabilizer, adsorbent, buffer, etc.
Document your checking process. If such a substance is identified, move it to the conditions section. Finally, re-output all the corrected reactions.

Checking Process:
- For Reaction 001:

Final Output:
[Repeat Reaction format for all corrected reactions]
"""


prompt_reaction_extraction = """
Extract all distinct chemical reactions from given content. Include each reaction only once, if it appears more than once.
Only extract reactions where all reactants and products are fully identified. Exclude any reaction if any reactant or product is unspecified.

You must strictly adhere to the following rules:
Confirm that the extracted reactions are interconnected: the reactants in each subsequent reaction must include at least one product from the previous reaction.
Be careful to ensure that solvents, catalysts, or other reaction conditions are not mistakenly listed as reactants. Reactants must only include compounds directly involved in the chemical transformation.
All names must be provided in their full specific forms:
- Maintain consistent naming of each substance throughout.
- Do not use general terms, such as "Diamine" or "Dianhydride," to represent a class of substances for any reactant substance name.
- Do not use abbreviations for any substance name.
- Do not include explanatory notes, brackets, or additional information.

Format the output strictly as follows:

Reaction 001:
Reactants: List the substances with specific names, separated by commas.
Products: List the substances with specific names, separated by commas.
Conditions: List the following in the exact order, separated by commas, skipping any condition that is not provided or is unknown:
- If a catalyst is provided, prefix with "Catalyst: " followed by the substances with specific names.
- If a solvent is provided, prefix with "Solvent: " followed by the substances with specific names.
- If other non-reactive reagents are provided, such as drying agents, stabilizers, adsorbents, buffers, etc., prefix with their respective names (e.g., 'Drying agents: ', 'Stabilizers: ', 'Adsorbents: ', 'Buffers: ') followed by their specific types.
- If an atmosphere is provided, prefix with "Atmosphere: " followed by the specified gas in full names.
- For temperature, provide the specific value or range with °C without any prefix.
- For pressure, provide the specific value or range with atm or bar without any prefix.
- For duration, provide the specific value or range with h, min, or d without any prefix.
- For yield, provide the specific value or range in % without any prefix.

[Repeat Reaction format for all extracted reactions]
"""

# --- Reaction Evaluation/Correction Prompt ---

prompt_reaction_evaluation = """
Reactions:
{reactions}

Check whether any substance in the Reactants is commonly used as a solvent, catalyst, drying agent, stabilizer, adsorbent, buffer, etc.
Document your checking process. If such a substance is identified, move it to the conditions section. Finally, re-output all the corrected reactions. Ensure the number of reactions in the output matches the original.

Format the output strictly as follows:

Checking process:
- Reaction 001: Describe your process for identifying and reassigning any substances to the conditions section.
- Reaction 002: ...

Final Output:

Reaction 001:
Reactants: List the substances with specific names, separated by commas.
Products: List the substances with specific names, separated by commas.
Conditions: List the following in the exact order, separated by commas, skipping any condition that is not provided or is unknown:
- If a catalyst is provided, prefix with "Catalyst: " followed by the substances with specific names.
- If a solvent is provided, prefix with "Solvent: " followed by the substances with specific names.
- If other non-reactive reagents are provided, such as drying agents, stabilizers, adsorbents, buffers, etc., prefix with their respective names (e.g., 'Drying agents: ', 'Stabilizers: ', 'Adsorbents: ', 'Buffers: ') followed by their specific types.
- If an atmosphere is provided, prefix with "Atmosphere: " followed by the specified gas in full names.
- For temperature, provide the specific value or range with °C without any prefix.
- For pressure, provide the specific value or range with atm or bar without any prefix.
- For duration, provide the specific value or range with h, min, or d without any prefix.
- For yield, provide the specific value or range in % without any prefix.

[Repeat Reaction format for all corrected reactions]
"""

# --- Entity Alignment Prompts ---

prompt_align_root_node = """
For the given reactions, check if any substances are different names for the same entity as "{substance}". If so, standardize all names referring to this specific entity to "{substance}".
Ensure that the output includes all original reactions, with only the inconsistent names for "{substance}" modified.
Output the same number of reactions as provided, maintaining the original format. No additional notes, brackets, or information should be included.

reactions:
{reactions}

[Output should be the full list of reactions in the standard format, with only the target substance name standardized]
"""

prompt_template_entity_alignment = """
Substances:
{substances}

Analyze the listed substances and determine if any are identical (completely the same substance, not isomers) but are referred to by different names. If identical substances are found, standardize their names to the most commonly used and scientifically accepted name.
Be meticulous in identifying identical substances and standardizing their names to avoid any omissions.

Provide your output in the following format for each group of identical substances found:

Different names for the same substance:
List all names that refer to the same substance, separated by commas.

Standardized name:
Provide the most commonly used name as the standardized version.

---

Example Output:

Different names for the same substance: Benzophenonetetracarboxylic dianhydride, Benzophenone tetracarboxylic dianhydride
Standardized name: Benzophenone tetracarboxylic dianhydride
---
Different names for the same substance: 2,2-bis(3',4'-dicarboxyphenyl)hexafluoropropane dianhydride, 2,2-bis(3,4-dicarboxyphenyl)hexafluoropropane dianhydride
Standardized name: 2,2-bis(3,4-dicarboxyphenyl)hexafluoropropane dianhydride
---
Different names for the same substance: 4,4′-(hexafluoro-isopropylidene) diphthalic anhydride, hexafluoroisopropylidene diphthalic anhydride
Standardized name: hexafluoroisopropylidene diphthalic anhydride
---
"""

# --- Tree Expansion Prompts ---

prompt_add_reactions_from_literature = """
Please answer the questions based on the given content.
How to use common laboratory and commercial chemical compounds to synthesize "{material}" in one or more steps? Please provide the reactions.

Ensure that the extracted reactions are interconnected: each subsequent reaction's reactants must include at least one product from the previous reaction. Continue this sequence until all reactants in the first reaction are common laboratory or commercially available substances.
Ensure to unify the substance name throughout. Be meticulous not to alter the chemical names in a way that changes the identity of the substances. Use the standardized name "{material}" when referring to the target product.

Format the output strictly as follows:

Reaction 001:
Reactants: List the IUPAC nomenclatures or common accepted names, separated by commas.
Products: List the IUPAC nomenclatures or common accepted names, separated by commas.
Conditions: List the following in the exact order, separated by commas, skip any condition that is not provided or if it is unknown:
- If a synthesis method is provided, provide a professional term used to describe the reaction (e.g., 'Method: Esterification')
- If a catalyst is provided, prefix with 'Catalyst: ' followed by the names.
- If a solvent is provided, prefix with 'Solvent: ' followed by the names.
- If other non-reactive reagents are provided, prefix with their respective names (e.g., 'Drying agents: ', 'Buffers: ') followed by their specific types.
- If an atmosphere is provided, prefix with 'Atmosphere: ' followed by the specified gas.
- For temperature, provide specific value or range with °C without any prefix.
- For pressure, provide specific value or range with atm or bar without any prefix.
- For duration, provide specific value or range with h or min or d without any prefix.
- For yield, provide specific value or range in % without any prefix.

Do not include any explanatory notes, brackets, or additional information.

[Repeat Reaction format for all added reactions]
"""


prompt_add_reactions_from_literature_cot = """
Please answer the questions based on the given content.
How to use common laboratory and commercial chemical compounds to synthesize "{material}" in one or more steps? Please provide the reactions.

You must strictly adhere to the following rules:
Confirm that the extracted reactions are interconnected: the reactants in each subsequent reaction must include at least one product from the previous reaction.
Be careful to ensure that solvents, catalysts, or other reaction conditions are not mistakenly listed as reactants. Reactants must only include compounds directly involved in the chemical transformation.
All names must be provided in their full specific forms:
- Maintain consistent naming of each substance throughout. Use the standardized name "{material}" when referring to the target product.
- Do not use general terms, such as "Diamine" or "Dianhydride," to represent a class of substances for any reactant substance name.
- Do not use abbreviations for any substance name.
- Do not include explanatory notes, brackets, or additional information.

Format the output strictly as follows:

Initial Output:

Reaction 001:
Reactants: List the substances with specific names, separated by commas.
Products: List the substances with specific names, separated by commas.
Conditions: List the following in the exact order, separated by commas, skipping any condition that is not provided or is unknown:
- If a catalyst is provided, prefix with "Catalyst: " followed by the substances with specific names.
- If a solvent is provided, prefix with "Solvent: " followed by the substances with specific names.
- If other non-reactive reagents are provided, prefix with their respective names (e.g., 'Drying agents: ', 'Stabilizers: ', 'Adsorbents: ', 'Buffers: ') followed by their specific types.
- If an atmosphere is provided, prefix with "Atmosphere: " followed by the specified gas in full names.
- For temperature, provide the specific value or range with °C without any prefix.
- For pressure, provide the specific value or range with atm or bar without any prefix.
- For duration, provide the specific value or range with h, min, or d without any prefix.
- For yield, provide the specific value or range in % without any prefix.

[Repeat for all initially extracted reactions]

Check whether any substance in the Reactants is commonly used as a solvent, catalyst, drying agent, stabilizer, adsorbent, buffer, etc.
Document your checking process. If such a substance is identified, move it to the conditions section. Finally, re-output all the corrected reactions.

Checking Process:
- For Reaction 001: [Describe check]
- For Reaction 002: [Describe check]

Final Output:

Reaction 001:
Reactants: List the substances with specific names, separated by commas.
Products: List the substances with specific names, separated by commas.
Conditions: [Formatted conditions]

[Repeat Reaction format for all corrected reactions]
"""

# --- Filtration Prompts ---

prompt_reactions_filtration = """
Reactions:
{reactions}

Given Reactions, please filter them according to these conditions:
1. Exclude reactions that lack any reaction conditions (i.e., the 'Conditions:' line is empty or missing).
2. Exclude reactions with high reaction temperatures > 350 °C.
3. Exclude reactions with high reaction pressure > 2 atm (or equivalent in bar, e.g., > 2.026 bar).
4. Exclude reactions involving catalysts known to be extremely expensive or difficult to source (e.g., highly specialized organometallics, specific enzymes not commercially available). Use general chemical knowledge.
5. Exclude reactions involving solvents known to be highly toxic, carcinogenic, or environmentally hazardous AND where common alternatives exist (e.g., prefer ethanol over benzene if possible; consider context). Use general chemical knowledge.
6. Exclude reactions explicitly stating the use of highly toxic reactants (e.g., phosgene, cyanides) or known formation of highly toxic byproducts, unless it's the only known route.

Format the output strictly as follows:

Excluded Reactions:
Reaction idx: [idx], Reason: [Provide concise reason based on rules 1-6]
Reaction idx: [idx], Reason: ...

Remaining Reactions:
[List the full reaction details ONLY for the reactions that PASS the filter, maintaining the original format]
Reaction idx: [idx]
Reactants: ...
Products: ...
Conditions: ...
---
Reaction idx: [idx]
Reactants: ...
Products: ...
Conditions: ...
"""

prompt_filter_pathway = """
Please evaluate the following reaction pathways one by one to determine their validity based on chemical feasibility and connectivity. If a pathway is not valid, remove it and explain why. A pathway is invalid if:
1.  A reactant in a step is not available from the product(s) of the *immediately preceding* step within that pathway (or isn't an initial reactant if it's the first step).
2.  A reaction step involves chemically incompatible transformations or reagents (e.g., using a strong base when an acid-labile group needs protection). Use general chemical knowledge.

Reaction Pathways:
{all_pathways}

Format the output strictly as follows:

Excluded Reaction Pathways:
Pathway: [List reaction indices in order (e.g., idx4, idx7, idx3)], Reason: [Provide concise explanation based on rules 1-2]
Pathway: [List reaction indices...], Reason: ...

Remaining Reaction Pathways:
[List ONLY the valid pathways, including their full reaction details, maintaining the format from the input {all_pathways}]

Example of Remaining Pathway Output:
Pathway: idx4, idx7, idx3
--- Reaction idx: idx4 ---
Reactants: A, B
Products: C
Conditions: Catalyst: Cat1, Solvent: Solv1, 100 °C, 1 atm, 2 h
Source: Source1
--- Reaction idx: idx7 ---
Reactants: C, D
Products: E
Conditions: Solvent: Solv2, 50 °C, 1 atm, 1 h
Source: Source2
--- Reaction idx: idx3 ---
Reactants: E
Products: TargetProduct
Conditions: Catalyst: Cat2, Solvent: Solv3, 150 °C, 1 atm, 5 h
Source: Source3
--- [End of Pathway] ---

Pathway: idx1, idx5
--- Reaction idx: idx1 ---
...
--- Reaction idx: idx5 ---
...
--- [End of Pathway] ---
"""

# --- Single Best Pathway Recommendation Prompts ---

recommend_prompt_template_general = """
Given the target product "{substance}",
please analyze and recommend the **single most optimal** reaction pathway by following these structured steps:

1. Pathways Analysis with Advantages and Disadvantages:
Please analyze **all** pathways listed in "Reaction Pathways" comprehensively without any omissions.
For each pathway, provide a summary of its key features and evaluate its advantages and disadvantages based on:
- Reaction Mildness: Evaluate if the reaction conditions are mild (low temperature/pressure, short duration).
- Reactant Availability and Cost: Consider accessibility and cost-effectiveness of initial reactants, solvents, catalysts, atmosphere.
- Yield and Scalability: Assess reported yield and feasibility for larger scale.
- Safety Profile: Analyze concerns regarding toxic, flammable, or hazardous substances.
- Number of Steps: Consider the overall pathway length.

2. Final Recommendation:
After evaluating all pathways, select the **single most suitable** pathway and justify your choice by comparing it with the alternatives. Highlight specific aspects that make it superior overall (balancing mildness, cost, yield, safety, steps).

Reaction Pathways:
{all_pathways}

Response Format:

Analysis:

Pathway: [List reaction indices (e.g. idx4, idx7, idx3)]
Advantages: [List advantages]
Disadvantages: [List disadvantages]

[Repeat Analysis for ALL provided pathways]

--- Recommended Reaction Pathway ---

Recommended Pathway: [List reaction indices (e.g. idx4, idx7, idx3)]

Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the chosen pathway]

Reasons for Recommendation:
Explain the rationale for selecting this pathway as the best overall, comparing it to the other alternatives based on the analysis criteria.
"""

recommend_prompt_template_cost = """
Given the target product "{substance}",
please analyze and recommend the **single reaction pathway** that results in the **lowest estimated reactant and reagent cost** by following these structured steps:

1. Reactant Inventory and Cost Analysis:
For **each** pathway listed in "Reaction Pathways," list all **initial reactants** (substances introduced externally, not intermediates), as well as any solvents, catalysts, and atmospheric requirements mentioned.
Analyze **all** pathways comprehensively. Evaluate cost-effectiveness based *exclusively* on the likely accessibility and relative cost of these initial reactants, solvents, catalysts, and atmospheric requirements. Assume standard laboratory chemical pricing unless specific cost data is inherently available. Do not consider yield, time, or number of steps for this analysis.

2. Final Recommendation:
Select the **single pathway** estimated to have the lowest overall reactant/reagent cost. Justify your choice by comparing it with the alternatives, highlighting specific cost-related factors (e.g., common vs. specialized chemicals, expensive catalysts) that make this pathway superior *in terms of cost*.

Reaction Pathways:
{all_pathways}

Response Format:

Inventory & Analysis:

Pathway: [List reaction indices (e.g., idx4, idx7, idx3)]
Initial Reactants: [List initial reactant names]
Solvents: [List solvent names]
Catalysts: [List catalyst names]
Atmospheric Requirements: [Specify atmosphere]
Cost Analysis (Estimated): [Provide qualitative assessment - e.g., Low (common reagents), Medium, High (expensive catalyst/reactant)]

[Repeat Inventory & Analysis for ALL provided pathways]

--- Recommended Reaction Pathway (Lowest Cost) ---

Recommended Pathway: [List reaction indices (e.g., idx4, idx7, idx3)]

Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the chosen pathway]

Reasons for Recommendation (Cost):
Explain the rationale for selecting this pathway as the lowest cost option, focusing exclusively on the cost factors identified in the analysis.
"""


recommend_prompt_template_condition = """
Given the target product "{substance}",
please analyze and recommend the **single reaction pathway** that has the **mildest overall reaction conditions** (temperature and pressure) by following these structured steps:

1. Comprehensive Condition Inventory and Analysis:
For **each** pathway listed in "Reaction Pathways," list all relevant reaction temperatures and pressures mentioned for each step.
Conduct a comprehensive analysis of **each** pathway to evaluate the mildness of reaction conditions strictly based on the listed temperature and pressure values (lower is generally milder). Consider both maximum values and the overall profile. Do not consider cost, yield, time, or number of steps.

2. Final Recommendation:
Select the **single pathway** with the mildest overall reaction conditions based exclusively on temperature and pressure. Justify your choice by comparing it with the alternatives, highlighting specific temperature- or pressure-related factors that make it superior *in terms of mildness*.

Reaction Pathways:
{all_pathways}

Response Format:

Condition Inventory and Analysis:

Pathway: [List reaction indices (e.g., idx4, idx7, idx3)]
Temperatures (°C or K): [List T for each step, e.g., Step1: 50, Step2: 25]
Pressures (atm, bar, etc.): [List P for each step, e.g., Step1: 1 atm, Step2: 1 atm]
Condition Analysis (Mildness): [Provide qualitative assessment - e.g., Very Mild, Mild, Moderate, Harsh]

[Repeat Condition Inventory & Analysis for ALL provided pathways]

--- Recommended Reaction Pathway (Mildest Conditions) ---

Recommended Pathway: [List reaction indices (e.g., idx4, idx7, idx3)]

Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the chosen pathway]

Reasons for Recommendation (Mildness):
Explain the rationale for selecting this pathway as having the mildest conditions, focusing solely on the temperature and pressure analysis.
"""

recommend_prompt_template_specific_substance = """
Reaction Pathways:
{all_pathways}

Please analyze and recommend the **single most optimal** reaction pathway **that includes "{initial_reactant}" as one of the initial reactants** (i.e., a reactant in the very first step of the pathway).

1. Filter Pathways:
Identify and list all reaction pathways from "Reaction Pathways" where "{initial_reactant}" is present in the 'Reactants' list of the *first* reaction index listed for that pathway.

2. Analysis of Filtered Pathways:
Analyze the filtered pathways comprehensively. Evaluate advantages and disadvantages based on: Reaction Mildness, Reactant Availability/Cost (beyond the required initial reactant), Yield/Scalability, Safety Profile, and Number of Steps.

3. Final Recommendation:
From the filtered pathways, select the **single most suitable** pathway overall and justify your choice by comparing it with the other filtered alternatives, considering the criteria from step 2.

Format the output strictly as follows:

Filtered Pathways (Containing "{initial_reactant}" initially):
Pathway: [List reaction indices (e.g., idx4, idx7, idx3)]
Pathway: [List reaction indices (e.g., idx1, idx5)]
[List all matching pathways]

Analysis of Filtered Pathways:

Pathway: [List reaction indices (e.g. idx4, idx7, idx3)]
Advantages: [List advantages]
Disadvantages: [List disadvantages]

[Repeat Analysis for ALL filtered pathways]

--- Recommended Reaction Pathway (Using "{initial_reactant}") ---

Recommended Pathway: [List reaction indices (e.g., idx4, idx7, idx3)]

Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the substances with specific names, separated by commas.
Products: List the substances with specific names, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the chosen pathway]

Reasons for Recommendation:
Explain the rationale for selecting this specific pathway over the other filtered alternatives that also use "{initial_reactant}", based on the analysis criteria.
"""

recommend_prompt_flubendiamide_example = """
Reaction Pathways:
{all_pathways}

Please analyze the provided reaction pathways and identify the one(s) most likely used to manufacture **commercially available flubendiamide**. Consider factors like common industrial practices, reagent availability, and known synthesis routes if possible from the provided data. Recommend the single most likely pathway.

Format the output strictly as follows:

--- Recommended Commercial Pathway for Flubendiamide ---

Recommended Pathway: [List reaction indices (e.g., idx4, idx7, idx3)]

Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the substances with specific names, separated by commas.
Products: List the substances with specific names, separated by commas.
Reaction SMILES: [Provide the full reaction SMILES string if available in source, otherwise state 'Not Available']
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
SourceLink: [Provide the source literature PDF link if available, otherwise state 'Not Available']
[Repeat Step details for each step in the chosen pathway]

Reasons for Recommendation:
Explain the rationale for selecting this pathway as the likely commercial route, citing any evidence from the data or general chemical knowledge about industrial synthesis.
"""


# --- Top 3 Pathway Recommendation Prompts ---

recommend_prompt_template_general_top3 = """
Given the target product "{substance}",
please analyze all provided reaction pathways and recommend the **top 3 most optimal pathways, ranked in order of preference**, by following these structured steps:

1. Pathways Analysis with Advantages and Disadvantages:
Please analyze **all** pathways listed in "Reaction Pathways" comprehensively without any omissions.
For each pathway, provide a summary of its key features and evaluate its advantages and disadvantages based on:
- Reaction Mildness: Evaluate if the reaction conditions are mild (low temperature/pressure, short duration).
- Reactant Availability and Cost: Consider accessibility and cost-effectiveness of initial reactants, solvents, catalysts, atmosphere.
- Yield and Scalability: Assess reported yield and feasibility for larger scale.
- Safety Profile: Analyze concerns regarding toxic, flammable, or hazardous substances.
- Number of Steps: Consider the overall pathway length.

2. Final Recommendation (Top 3):
After evaluating all pathways, select and rank the **top 3** most suitable pathways overall. Justify the ranking by comparing them with each other and the alternatives. Highlight specific aspects that make them superior (balancing mildness, cost, yield, safety, steps).

Reaction Pathways:
{all_pathways}

Response Format:

Analysis:

Pathway: [List reaction indices (e.g. idx4, idx7, idx3)]
Advantages: [List advantages]
Disadvantages: [List disadvantages]

[Repeat Analysis for ALL provided pathways]

--- Top 3 Recommended Pathways (Ranked by General Optimality) ---

Rank 1: Pathway [List reaction indices (e.g. idx4, idx7, idx3)]
Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the Rank 1 pathway]
Reasons for Rank 1: Explain the rationale for selecting this pathway as #1 overall, comparing it to others.

Rank 2: Pathway [List reaction indices (e.g. idx1, idx5)]
Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the Rank 2 pathway]
Reasons for Rank 2: Explain the rationale for selecting this pathway as #2 overall.

Rank 3: Pathway [List reaction indices (e.g. idx2, idx9, idx6)]
Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the Rank 3 pathway]
Reasons for Rank 3: Explain the rationale for selecting this pathway as #3 overall.

Overall Comparison and Justification:
Provide a summary comparing the top 3 pathways and justifying the final ranking based on the analyzed criteria.
"""

recommend_prompt_template_cost_top3 = """
Given the target product "{substance}",
please analyze all provided reaction pathways and recommend the **top 3 pathways ranked by lowest estimated reactant and reagent cost** by following these structured steps:

1. Reactant Inventory and Cost Analysis:
For **each** pathway listed in "Reaction Pathways," list all **initial reactants** (substances introduced externally, not intermediates), as well as any solvents, catalysts, and atmospheric requirements mentioned.
Analyze **all** pathways comprehensively. Evaluate cost-effectiveness based *exclusively* on the likely accessibility and relative cost of these initial reactants, solvents, catalysts, and atmospheric requirements. Assume standard laboratory chemical pricing unless specific cost data is inherently available. Do not consider yield, time, or number of steps for this analysis.

2. Final Recommendation (Top 3 by Cost):
Select and rank the **top 3** pathways with the lowest overall estimated reactant/reagent cost. Justify your ranking by comparing them with each other and the alternatives, highlighting specific cost-related factors (e.g., common vs. specialized chemicals, expensive catalysts) that differentiate them *in terms of cost*.

Reaction Pathways:
{all_pathways}

Response Format:

Inventory & Analysis:

Pathway: [List reaction indices (e.g., idx4, idx7, idx3)]
Initial Reactants: [List initial reactant names]
Solvents: [List solvent names]
Catalysts: [List catalyst names]
Atmospheric Requirements: [Specify atmosphere]
Cost Analysis (Estimated): [Provide qualitative assessment - e.g., Low (common reagents), Medium, High (expensive catalyst/reactant)]

[Repeat Inventory & Analysis for ALL provided pathways]

--- Top 3 Recommended Pathways (Ranked by Lowest Estimated Cost) ---

Rank 1 (Lowest Cost): Pathway [List reaction indices (e.g. idx4, idx7, idx3)]
Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the Rank 1 pathway]
Reasons for Rank 1 (Cost): Explain why this pathway is estimated to be the lowest cost, based *only* on cost factors.

Rank 2: Pathway [List reaction indices (e.g. idx1, idx5)]
Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the Rank 2 pathway]
Reasons for Rank 2 (Cost): Explain why this pathway is estimated to be the second lowest cost, based *only* on cost factors.

Rank 3: Pathway [List reaction indices (e.g. idx2, idx9, idx6)]
Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the Rank 3 pathway]
Reasons for Rank 3 (Cost): Explain why this pathway is estimated to be the third lowest cost, based *only* on cost factors.

Overall Comparison and Justification (Cost):
Provide a summary comparing the top 3 pathways based purely on estimated cost factors.
"""

recommend_prompt_template_condition_top3 = """
Given the target product "{substance}",
please analyze all provided reaction pathways and recommend the **top 3 pathways ranked by mildest overall reaction conditions** (temperature and pressure) by following these structured steps:

1. Comprehensive Condition Inventory and Analysis:
For **each** pathway listed in "Reaction Pathways," list all relevant reaction temperatures and pressures mentioned for each step.
Conduct a comprehensive analysis of **each** pathway to evaluate the mildness of reaction conditions strictly based on the listed temperature and pressure values (lower is generally milder). Consider both maximum values and the overall profile. Do not consider cost, yield, time, or number of steps.

2. Final Recommendation (Top 3 by Mildness):
Select and rank the **top 3** pathways with the mildest overall reaction conditions based exclusively on temperature and pressure. Justify your ranking by comparing them with each other and the alternatives, highlighting specific temperature- or pressure-related factors that make them superior *in terms of mildness*.

Reaction Pathways:
{all_pathways}

Response Format:

Condition Inventory and Analysis:

Pathway: [List reaction indices (e.g., idx4, idx7, idx3)]
Temperatures (°C or K): [List T for each step, e.g., Step1: 50, Step2: 25]
Pressures (atm, bar, etc.): [List P for each step, e.g., Step1: 1 atm, Step2: 1 atm]
Condition Analysis (Mildness): [Provide qualitative assessment - e.g., Very Mild, Mild, Moderate, Harsh]

[Repeat Condition Inventory & Analysis for ALL provided pathways]

--- Top 3 Recommended Pathways (Ranked by Mildest Conditions) ---

Rank 1 (Mildest Conditions): Pathway [List reaction indices (e.g. idx4, idx7, idx3)]
Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the Rank 1 pathway]
Reasons for Rank 1 (Mildness): Explain why this pathway has the mildest conditions, based *only* on T/P analysis.

Rank 2: Pathway [List reaction indices (e.g. idx1, idx5)]
Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the Rank 2 pathway]
Reasons for Rank 2 (Mildness): Explain why this pathway has the second mildest conditions, based *only* on T/P analysis.

Rank 3: Pathway [List reaction indices (e.g. idx2, idx9, idx6)]
Details:
Step 1:
Reaction idx: Specify the reaction index.
Reactants: List the IUPAC nomenclatures, separated by commas.
Products: List the IUPAC nomenclatures, separated by commas.
Conditions: List conditions in specified order; skip if unknown.
Source: Source literature name.
[Repeat Step details for each step in the Rank 3 pathway]
Reasons for Rank 3 (Mildness): Explain why this pathway has the third mildest conditions, based *only* on T/P analysis.

Overall Comparison and Justification (Mildness):
Provide a summary comparing the top 3 pathways based purely on temperature and pressure conditions.
"""