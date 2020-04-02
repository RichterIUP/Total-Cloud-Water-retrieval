#!/usr/bin/python3

import datetime as dt

def month_to_number(month):
    if month == "Jan":
        return 1

def group(array):
    indices = []
    idx = [0]

    for ii in range(len(array)-1):
        if array[ii] == array[ii+1]:
            idx.append(ii+1)
        else:
            indices.append(idx)
            idx = [ii+1]
    
    return indices#
    
with open("log", "r") as f:
    cont = f.readlines()[3:]
    
spec = []
dates = []

for element in cont:
    line = element.split(" ")
    try:
        dates.append(dt.datetime(2019, 1, int(line[5][:-1]), int(line[7].split(":")[0]), int(line[7].split(":")[1])))
    except ValueError:
        dates.append(dt.datetime(2019, 1, int(line[6][:-1]), int(line[8].split(":")[0]), int(line[8].split(":")[1])))
        
    spec.append(line[-1].split(".nc")[0])

indices  = group(spec)
print(indices)
with open("spectra.dat", "w") as f:
    for loop in range(len(indices)):
        print(indices[loop])
        max_date = max(dates[indices[loop][0]:indices[loop][-1]])
        dates_slice = []
        for element in indices[loop]:
            dates_slice.append(dates[element])
        idx = dates_slice.index(max_date)
        print(idx)
        #print(dates[idx])
        f.write("{}\n".format(cont[idx+indices[loop][0]].split(" ")[-1].rstrip()))

#print(indices[0], indices[1])
#for element in indices[0]:
#    print(spec[element])
#for element in indices[1]:
#    print(spec[element])
