#!/usr/bin/env python


from dpgen.analysis.main_func_1 import data_collections

def dpgen_analysis(args):
    while True:
        print("\n")
        print("=======================================================")
        print("                    Main Menu                         ")
        print(" (  1)  Data Collection Functions                     ")
        print("\n Tips: Input q or -10 to exit program               ")
        print("=======================================================")
        isel = input("\n Please input the manu index: ")

        if isel == str(-10) or isel == "q":
            break

        elif isel == str(1):
            data_collections()
