#!/usr/bin/python3

import os

path = "LIWP_singlelayer_"
for ii in range(33):
    with open("{}{}".format(path, ii), "r") as f:
        cont = f.readlines()
        for element in cont:
            if "TIR" in element:
                with open("inputs_{:02d}".format(ii), "a") as g:
                    print(element.split(" ")[2])
                    g.write("/home/phi.richter/Emission_Data/{}\n".format(element.split(" ")[2]))

    
