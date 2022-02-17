#!/usr/bin/env python

def dpgen_analysis(args):
    while True:
        print("\n")
        print("=======================================================")
        print("                    Main Menu                         ")
        print(" (  1)  Various functions used in dpgen iteration     ")
        print(" (  2)  Various functions used in deepmd-kit          ")
        print(" (  3)  DPGEN autotest module                         ")
        print(" (  4)  Deep potential molecular dynamics (DPMD) Analysis")
        print(" (  5)  Various functions used in fp calculations     ")
        print(" (300)  Other function I                              ")
        print("\n Tips: Input q or -10 to exit program               ")
        print("=======================================================")
        isel = input("\n Please input the manu index: ")

        if isel == str(-10) or isel == "q":
            break
