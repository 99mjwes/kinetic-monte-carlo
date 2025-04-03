# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 13:24:40 2022

@author: morte
"""

import numpy as np
import matplotlib.pyplot as plt


def input_number(prompt, a = -np.inf, b = np.inf, IsInteger = False):
    """
    Promps the user to input a number. The range and type is specified by the calling function
    This function is used for all number inputs.

    Parameters
    ----------
    prompt : string
        The message to display to the user
    a : float (optional)
        The minimum accepted value of input
    b : float (optional)
        The maximum accepted value of input
    IsInteger : boolean (optional)
        Determines whether the input is an integer or floating point number

    Returns
    -------
    num : float / integer
        User-provided number

    """
    
    #Check if user-input is integer
    if IsInteger:
        #Convert to integer
        text = 'an integet' #Displayed number type
        #a = int(a)
        #b = int(b)
    else:
        text = 'et tal' #Displayed number type
    
    #Input loop
    while True:
        try:
            num = float(input(prompt))
            if a <= num and b >= num: #Check that input is within limits
                if IsInteger:
                    if int(num) != num: #Warn the user if input is not an integer
                        print('Expected an integer - ' + str(num) + ' interpreted as ' + str(int(num)))
                    num = int(num) #Convert to integer
                break
            print('Input must be ' + text + ' bewteen ' + str(a) + ' and ' + str(b))
        except ValueError: #If input is not a number 
            print('Input must be a number!')
            pass
        
    return num

#_____________________________________________________________________________


def display_menu(options, promptUser = True):
    """
    Prints a list of available options and prompts the user to select one.
    The user must select an option to continue

    Parameters
    ----------
    options : array
        an array containing the options to choose from
        
    promptUser : Boolean (optional)
        Suppresses the user prompt if False

    Returns
    -------
    choice : integer
        A number representing the selected option

    """
    for i,x in enumerate(options):
        print(f"{i+1}. {x}")
        
    choice = 0
    if promptUser:
        #while not(np.any(choice == np.arange(n) + 1)): #Not nessesary, inputNumber handles user input 
        choice = input_number("Please select an option: ",1, len(options), IsInteger = True)
        
    return choice

#_____________________________________________________________________________


def data_plot(data):
    """
    Generates plots based on given data

    Parameters
    ----------
    data : array
        a (n by 3) numpy data array


    """
    Warning("This function is not implemented yet")
    
    #print('Plot generation complete')
    
#_____________________________________________________________________________


def load_help():
    """
    Imports Help Topics and Text from an external text-file

    Returns
    -------
    helpTopic : array
        an array of Help Topics
    helpText : array
        an array of Help Main Texts
    helpIsLoaded : boolean
        Represents whether the help text was loaded successfully

    """
    
    line = 'TestIndex0'     #Index 0 will be ignored
    topic = 'TestIndex0'    #Index 0 will be ignored
    hlptxt = np.array(['']) #Empty help-text array
    hlptpc = np.array(['']) #Empty help-topic array
    helpIsLoaded = False    #Help file is not loaded
        
    with open('functions\\hlp.txt','rt') as helpfile: #Opening help file, if it exists
        #Importing help topics and text
        while line != '':
            #Appending text
            hlptxt = np.append(hlptxt,line)
            hlptpc = np.append(hlptpc,topic)
            #Reading Text
            topic = helpfile.readline()
            topic = topic.replace('\n','') #Topics can only occupy one line
            line = helpfile.readline()
            line = line.replace('&','\n') #MainTexts can occupy multiple lines
            helpIsLoaded = True
        
        #The text is appended first in the loop to get the empty string in the start of the array.
        #First or last is a matter of preference.
    
    return hlptpc[1:], hlptxt[1:], helpIsLoaded #First empty string is filtered out

    
#_____________________________________________________________________________