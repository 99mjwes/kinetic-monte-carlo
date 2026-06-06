import numpy as np
import pandas as pd
import os

from util import *

main_menu_items = np.array(["Create New Reaction Profile", "Add Reactions", "Modify Reaction", "Remove Reaction", "Save Reaction Profile", "Load Reaction Profile", "Export Reaction Profile", "Exit"])


def input_reaction(index):
    """
    Prompts the user to input a chemical reaction and associated parameters, then returns a DataFrame representing the reaction.
    The function expects the reaction to be entered in the format "reactant1 + reactant2 -> product1 + product2".
    It also prompts the user to input constants A, B, and C, and an equation for the reaction.
    Args:
        index (int): The index of the reaction, used as the row index in the resulting DataFrame.
    Returns:
        pd.DataFrame: A DataFrame containing the reaction details, including reactants, products, constants, and the equation.
                      If the reaction format is invalid, an empty DataFrame is returned.
    Notes:
        - Reactants and products should be separated by " + ".
        - The reaction arrow should be "->".
        - Constants A, B, and C are numerical values input by the user.
        - The equation is a string input by the user, with a default value of "default" if left empty.
    """

    reaction = input("Enter the reaction R" + str(index) + ": ")
    if reaction.count("->") != 1:
        if reaction != "": print("Invalid reaction format")
        return pd.DataFrame()
    

    reactants, products = reaction.split("->")

    reactants = reactants.split(" + ")
    products = products.split(" + ")

    reactants = [r.strip() for r in reactants]
    products = [p.strip() for p in products]

    reactants_count = []
    reactants_list = []
    for reactant in reactants:
        for i, char in enumerate(reactant):
            if char.isalpha():
                count = 1
                if i > 0:
                    count = int(reactant[:i])
                reactants_count.append(count)
                reactants_list.append(str(reactant[i:]))
                break

    products_count = []
    products_list = []
    for product in products:
        for i, char in enumerate(product):
            if char.isalpha():
                count = 1
                if i > 0:
                    count = int(product[:i])
                products_count.append(count)
                products_list.append(str("§" + product[i:]))
                break

    all_reactants = reactants_list + [p[1:] for p in products_list]
    all_products = ["§" + i for i in all_reactants]

    A = input_number("Enter constant A: ")
    B = input_number("Enter constant B: ")
    C = input_number("Enter constant C: ")
    eq = input(f"Enter equation for R{index}: ")
    if eq == "":
        eq = "default"

    total_df = pd.DataFrame(0, columns=all_reactants + all_products, index = [index], dtype=int)
    parameter_df = pd.DataFrame([[A, B, C]], columns=["!A", "!B", "!C"], index = [index], dtype=float)
    equation_df = pd.DataFrame([eq], columns=["!Eq"], index = [index], dtype=str)
    reaction_df = pd.DataFrame([reactants_count], columns=reactants_list, index = [index], dtype=int)
    product_df = pd.DataFrame([products_count], columns=products_list, index = [index], dtype=int)


    new_reaction = pd.concat([parameter_df, equation_df, reaction_df, product_df, total_df], join='outer', axis=1)
    new_reaction = new_reaction.loc[:,~new_reaction.columns.duplicated()].copy()
    return new_reaction


def remove_unused_species(df):
    """
    Remove unused species from the DataFrame.
    This function takes a DataFrame `df` and removes columns corresponding to 
    species that are not used in any reactions. A species is considered used 
    if the sum of its values in the DataFrame (both as a reactant and a product) 
    is greater than zero.
    Parameters:
    df (pandas.DataFrame): The input DataFrame containing species data. 
                           Columns starting with '§' represent products, 
                           and columns starting with '!' represent metadata.
    Returns:
    pandas.DataFrame: A DataFrame with only the columns corresponding to used species 
                      and the metadata columns '!A', '!B', '!C', and '!Eq'.
    """

    reactant_names = [column for column in df.columns if (not column.startswith("§") and not column.startswith("!"))]
    keep_names = []
    for name in reactant_names:
        if (df[name].sum() + df["§" + name].sum()) > 0:
            keep_names.append(name)

    df.index = range(1, len(df) + 1)

    return df[["!A", "!B", "!C", "!Eq"] + keep_names + ["§" + name for name in keep_names]]


def stats(df):
    """
    Calculate statistics from a DataFrame containing reaction data.
    Parameters:
    df (pandas.DataFrame): DataFrame containing reaction data. The DataFrame is expected to have columns 
                           representing different species and reactions.
    Returns:
    tuple: A tuple containing:
        - number_of_reactions (int): The total number of reactions.
        - number_of_species (int): The total number of species.
        - index (int): The index for the next reaction.
    """

    number_of_reactions = df.index[-1] if not df.empty else 0
    number_of_species = (len(df.columns) - 4)//2 if not df.empty else 0
    index = number_of_reactions + 1

    return number_of_reactions, number_of_species, index

def main():

    print("Welcome to the Reaction Profile Generator")
    print("Please select an option from the menu below:")
    reaction_profile = pd.DataFrame()
    number_of_reactions, number_of_species, index = stats(reaction_profile)
    while True:
        print(25*"=")
        print(f"Current Profile has {number_of_reactions} reactions and {number_of_species} species \n")
        choice = display_menu(main_menu_items)

        if choice == 1:
            reaction_profile = pd.DataFrame()
            number_of_reactions, number_of_species, index = stats(reaction_profile)
            print("Created New Reaction Profile")

        elif choice == 2:
            print("Add Reactions")
            new_reaction = input_reaction(index)
            while not new_reaction.empty:
                
                if reaction_profile.empty:
                    reaction_profile = new_reaction.copy(deep=True)
                    print("Copied new reaction")
                else:
                    reaction_profile = reaction_profile.merge(new_reaction, how='outer').sort_index(axis=1)
                    reaction_profile.index += 1
                    print("Merged new reaction")
                print("Reaction " + str(index) + " added")
                reaction_profile[reaction_profile.columns[4:]] = reaction_profile[reaction_profile.columns[4:]].fillna(int(0)).astype(int)
                index += 1
                new_reaction = input_reaction(index)

            number_of_reactions, number_of_species, index = stats(reaction_profile)
            print(reaction_profile)

        elif choice == 3: 
            print("Modify Reaction")
            if index == 1:
                print("No reaction to modify")
                continue
            reaction_index = input_number("Enter the reaction index to modify: ", a = 1, b = index - 1, IsInteger = True)
            if reaction_index == -1:
                print("Invalid reaction index")
                continue
            print("Current Reaction")
            # print(reaction_profile.loc[reaction_index].T)

            new_reaction = input_reaction(reaction_index)
            if not new_reaction.empty:

                if not all([species in list(reaction_profile.columns) for species in list(new_reaction.columns)]):
                    print("Modified reaction has undefined species! Please add a new reaction instead...")
                    continue
                new_reaction = new_reaction.reindex(columns=reaction_profile.columns, fill_value=0)
                reaction_profile.loc[reaction_index] = new_reaction.loc[reaction_index]
                reaction_profile = remove_unused_species(reaction_profile)
                print("Reaction " + str(reaction_index) + " modified")
                number_of_reactions, number_of_species, index = stats(reaction_profile)
                print(reaction_profile)
            else:
                print("Reaction " + str(reaction_index) + " not modified")

        elif choice == 4:
            print("Remove Reaction")
            if index == 1:
                print("No reaction to remove")
                continue
            reaction_index = input_number("Enter the reaction index to remove: ", a = 1, b = index - 1, IsInteger = True)
            if reaction_index == -1:
                print("Invalid reaction index")
                continue
            reaction_profile = reaction_profile.drop(reaction_index)
            reaction_profile = remove_unused_species(reaction_profile)
            # reaction_profile.index -= 1
            number_of_reactions, number_of_species, index = stats(reaction_profile)
            print("Reaction " + str(reaction_index) + " removed")
            print(reaction_profile) 


        elif choice == 5:
            print("Save Reaction Profile")
            if index == 1:
                print("No reaction profile to save")
                continue

            filename = input("Enter the filename to save the reaction profile: ").strip()
            if filename == "":
                print("Filename cannot be empty.")
            else:
                if not filename.endswith(".pkl"):
                    filename += ".pkl"
                reaction_profile.to_pickle(filename)
                print(f"Reaction profile saved to {filename}")

        elif choice == 6:
            print("Load Reaction Profile")
            filename = input("Enter the filename to load the reaction profile: ").strip()
            if filename == "":
                print("Filename cannot be empty.")
                continue
            if not filename.endswith(".pkl"):
                filename += ".pkl"
            if not os.path.exists(filename):
                print(f"File {filename} does not exist.")
                continue
            reaction_profile = pd.read_pickle(filename)
            number_of_reactions, number_of_species, index = stats(reaction_profile)
            print(reaction_profile.columns)
            print(list(reaction_profile.columns))
            print(reaction_profile)

        elif choice == 7:
            print("Export Reaction Profile")
            if index == 1:
                print("No reaction profile to export")
                continue

            desc = input("Enter the description for the reaction profile (optional): ").strip()

            columns = reaction_profile.columns
            reactant_names = [col for col in columns if (not col.startswith("§") and not col.startswith("!"))]

            export_df = reaction_profile.copy()
            export_quantity = pd.DataFrame("", columns=columns, index=[0], dtype=str)
            for species in reactant_names:
                export_quantity[species] = str(input_number(f"Enter the initial quantity of {species}: ", a = 0, IsInteger = True))
            export_df = pd.concat([export_quantity, export_df], join='outer', axis=0, ignore_index=False, sort=True)

            export_df.rename({old_name: new_name for old_name, new_name in zip(columns, ["A", "B", "C", "Eq"] + reactant_names + [""] * len(reactant_names))}, axis=1, inplace=True)

            filename = input("Enter the filename to export the reaction profile: ").strip()
            if filename == "":
                print("Filename cannot be empty.")
                continue
            
            if not filename.endswith(".tsv"):
                filename += ".tsv"
            export_df.to_csv(filename, sep="\t", index=False)

            with open(filename, 'r') as f:
                lines = f.readlines()
            lines.insert(0, desc + "\n")
            with open(filename, 'w') as f:
                f.writelines(lines)

            print(f"Reaction profile exported to {filename}")

        elif choice == 8:
            print("Exit")
            break
        else:
            print("Invalid choice")
        
    print("Goodbye")

if __name__ == "__main__":
    main()
