#!/usr/bin/python
import datetime as dt

import inp

def get_weekday(num):
    if(num == 0):
        return "Mo"
    elif(num == 1):
        return "Di"
    elif(num == 2):
        return "Mi"
    elif(num == 3):
        return "Do"
    elif(num == 4):
        return "Fr"
    elif(num == 5):
        return "Sa"
    else:
        return "So"
    
def get_month(num):
    if(num == 1):
        return "Jan"
    elif(num == 2):
        return "Feb"
    elif(num == 3):
        return "Mar"
    elif(num == 4):
        return "Apr"
    elif(num == 5):
        return "Mai"
    elif(num == 6):
        return "Jun"
    elif(num == 7):
        return "Jul"
    elif(num == 8):
        return "Aug"
    elif(num == 9):
        return "Sep"
    elif(num == 10):
        return "Okt"
    elif(num == 11):
        return "Nov"
    else:
        return "Dez"
 

#MODEL   selects atmospheric profile
#= 0  user supplied atmospheric profile
#= 1  tropical model
#= 2  midlatitude summer model
#= 3  midlatitude winter model
#= 4  subarctic summer model
#= 5  subarctic winter model
#= 6  U.S. standard 1976
def rundecker(z, p, t, w, tape5, co2_ppm, o3_ppm, atm, hmd_unit, sample):
    co2_man = 4*200.0
    model = 3#LBLRTM user defined wavenumber range
    aprofile = atm
    w_units = 'g/m3'
    v10 = True
    mlayers=z
    silent=False
    #sfc_temp=300
    #sfc_emis=0.98
    sc=3
    od_only=1
    wnum1=300.0
    wnum2=2100.0
    #tape5='tp5'
    view_angle=0.0
    
    version = 1.37
    
    if(len(z) > 0):
        have_profile = 1
    else:
        have_profile = 0
        
    #Scaling factors (Default)
    h2o_sf = [1.0, 0]
    co2_sf = [400.0, 0]#[1.0, 0]#[400, 1]
    o3_sf = [1.0, 0]
    co_sf = [1.0, 0]
    ch4_sf = [1.,0]
    n2o_sf = [1.0, 0]
    o2_sf = [1.0, 0]
    ccl4_sf = [0.1105, 1]
    f11_sf = [0.2783, 1]
    f12_sf = [0.5027, 1]
    
    sample = 4.0#0.3
    direction = "downwelling"
    inst = 3
    xsec = 1
    
    #Continuum
    cntnm = 1
    reset_cntnm = False
    cntnm_default = 1
    iemit = 0
    scan = 0
    merge = 1
    
    co2_mix = 360.0
    numangs = 0
    iout = 0
    icld = 0
    iatm = 1
        
    numbers0 = '         1         2         3         4         5         6         7         8         9'
    numbers1 = '123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 '
    rec_1_2 = ' HI=1 F4=1 CN={}'.format(cntnm)
    rec_1_2 = '{} AE=0 EM={} SC={} FI=0 PL=0'.format(rec_1_2, iemit, scan) + \
                ' TS=0 AM={} MG={} LA=0 MS=0'.format(iatm, merge) + \
                ' XS={}   00   00'.format(xsec)
    rec_1_3 = '{:10.3f}{:10.3f}{:10.3f}'.format(wnum1,wnum2,sample) + \
    '                              0.0002    0.001' + \
    '                             7'

    flag = 0
    msign = 1
    
    rec_3_1 = \
        '    {}    2{:5d}'.format(flag, msign*len(mlayers)) + \
        '    1    1    7    1'

    weekday = get_weekday(dt.datetime.now().weekday())
    month = get_month(dt.datetime.now().month)
    day = dt.datetime.now().day
    hour = dt.datetime.now().hour
    minute = dt.datetime.now().minute
    sec = dt.datetime.now().second
    year = dt.datetime.now().year
    date = "{} {}. {} {}:{}:{:2d} CET {}".format(weekday, day, month, hour, minute, sec, year)
    
    f = open(tape5, "w")
    print("Writing {}...".format(tape5))
    acomment = "Rundeck created on {} by rundecker.py - rewrite from rundecker.pro (v{})".format(date, version)
    f.write("{}\n".format(acomment))
    f.write("{}\n".format(numbers0))
    f.write("{}\n".format(numbers1))
    f.write('$ None\n')
    f.write("{}\n".format(rec_1_2))
    f.write("{}\n".format(rec_1_3))
                                                                                                       
    f.write("1m11111\n")
    co2_sf[0] = co2_sf[0] / 1e6
    f.write("{:15.7E}{:15.7E}{:15.7E}{:15.7E}{:15.7E}{:15.7E}{:15.7E}\n".format(h2o_sf[0], \
            co2_sf[0], o3_sf[0], n2o_sf[0], co_sf[0], ch4_sf[0], o2_sf[0]))
    f.write("{}\n".format(rec_3_1))
    h1 = mlayers[0]
    h2 = mlayers[-1]
    f.write("{:10.3f}{:10.3f}{:10.3f}\n".format(h1, h2, view_angle))
    line = ''
    for i in range(len(mlayers)):
        line = line + "{:10.3f}".format(mlayers[i])
        if((i+1) %8 == 0):
            f.write("{}\n".format(line))
            line = ''
    f.write("{}\n".format(line))
    p_units = 'mb'
    t_units = 'K'
    o3_units = 'ppmv'
    JCHARP = 'A'
    JCHART = 'A'
    if inp.PROF_CO2:
        co2_unit = "A"
    else:
        co2_unit = inp.PREDEF_ATM
    if inp.PROF_O3:
        o3_unit = "C"
    else:
        o3_unit = inp.PREDEF_ATM
    co2_unit = inp.PREDEF_ATM
    o3_unit = inp.PREDEF_ATM
    #co2_unit = 'A'
    JCHAR  = '{}{}{}4444'.format(hmd_unit, co2_unit, o3_unit)#'A444444'#Klimatologie fuer alles ausser Wasser
    #Ansosten: A -> ppmv, B -> cm-3, C -> g/kg, D -> g/m3 (so kann man retrievte Spurengase verwenden)
    zz = z
    tt = t
    pp = p
    ww = w
    #co2_man = 200.0
    co2 = [co2_man for i in range(len(pp))]#[co2_ppm[i] for i in range(len(pp))]#[co2_ppm[i] for i in range(len(pp))]#[0.0 for i in range(len(pp))]
    oo3 = [0.0 for i in range(len(pp))]#[o3_ppm[i] for i in range(len(pp))]#[o3_ppm[i] for i in range(len(pp))]
    pp_sorted = []
    for element in sorted(pp):
        pp_sorted.append(element)
    reversed_pp = pp[::-1]
    nhits = 0
    for ii in range(len(pp_sorted)):
        if(pp_sorted[ii] != reversed_pp[ii]):
            #print(pp[ii], reversed_pp[ii])
            nhits = nhits + 1
    if(nhits > 0):
        print("Pressure array is not monotonically increasing - quitting")
        return False
    
    #Entferne doppelte Drucklevel
    pp = [i for i in reversed(sorted(set(pp)))]
    index = []
    for ii in range(len(p)-1):
        if(p[ii] == p[ii+1]):
            continue
        else:
            index.append(ii)
    index.append(len(p)-1)

    pp = [p[i] for i in index]
    zz = [z[i] for i in index]
    tt = [t[i] for i in index]
    ww = [w[i] for i in index]
    oo3 = [0.0 for i in index]#[o3_ppm[i] for i in index]
    co2 = [co2_man for i in index]#[co2_ppm[i] for i in index]
    print(zz)
    
    inlayers = len(zz)
    p_comment = "User supplied profile"
    f.write("{:5d} {}\n".format(inlayers, p_comment))
    
    for i in range(inlayers): 
        f.write("{:10.4f}{:10.4f}{:10.3E}     {}{}   {}\n".format(zz[i], pp[i], tt[i], JCHARP, JCHART, JCHAR))
        f.write("{:10.3E}{:10.3E}{:10.3E}{:10.3E}{:10.3E}{:10.3E}{:10.3E}\n".format(\
                ww[i], co2[i], oo3[i], 0, 0, 0, 0))
        
    if(xsec == 1 and have_profile == 1):
        f.write('    3    0    0  The following cross-sections were selected:\n')
        f.write('CCL4      F11       F12\n')
    index = [0, len(zz)-1]
    f.write("{:5d}    0 XS 1995 UNEP values\n".format(len(index)))
    for i in range(len(index)):
        f.write("{:10.3f}     AAA\n{:10.3E}{:10.3E}{:10.3E}\n".format(zz[index[i]], \
                ccl4_sf[0]/1e3, f11_sf[0]/1e3, f12_sf[0]/1e3))
    f.write("-1.\n")
    f.write("%%%\n")
    f.close()
    exit(-1)
    return True

if __name__ == "__main__":
    z = [0.0 for i in range(51)]
    p = [0.0 for i in range(51)]
    t = [0.0 for i in range(51)]
    q = [0.0 for i in range(51)]
    z = [0.0, 0.008, 0.032, 0.072, 0.128, 0.2, 0.288, 0.392, 0.512, 0.648, 0.8,  0.968, \
         1.152, 1.352, 1.568, 1.8, 2.048, 2.312, 2.592, 2.888, 3.2, 3.528, 3.872, 4.232, \
         4.608, 5.0, 5.408, 5.832, 6.272, 6.728, 7.2, 7.688, 8.192, 8.712, 9.248, 9.8, 10.368, \
         10.952, 11.552, 12.168, 12.8, 13.448, 14.112, 14.792, 15.488, 16.2, 16.928, 17.672, 18.432, \
         19.208, 20.0]
    p = [10.0691935484, 10.0691935484, 10.0418018489, 9.9920290033, 9.92204733044, 9.83305739721, \
         9.72480404796, 9.59846587394, 9.4546008537, 9.29330963153, 9.11648353638, 8.9250451243, \
         8.72012823129, 8.50210360138, 8.27116639365, 8.02974567835, 7.77918072235, 7.52126212776, \
         7.25811133593, 6.98847365753, 6.71383820697, 6.43499834147, 6.15212990636, 5.86672343454, \
         5.58081220485, 5.29395691368, 5.00805348591, 4.72491117985, 4.44361520906, 4.16715015576, \
         3.89510853649, 3.62866902235, 3.36849977893, 3.11611259592, 2.87506561076, 2.64479865907, \
         2.42621790995, 2.2232320913, 2.03318596144, 1.8564116623, 1.69143583689, 1.5370713746, \
         1.39388516378, 1.26158875805, 1.13884079964, 1.02527254473, 0.921760818389, 0.826728634711, \
         0.739904402865, 0.660337703568,  0.58818380923]
    t = [273.612903226,273.612903226, 273.057738124, 272.673875423, 272.517054109, 271.970046742, \
         271.716935484,271.222147793, 270.460365738, 269.839595627, 270.086802977, 269.377582495, \
         269.694718016, 268.379542882, 266.914175907, 266.257289705, 265.737061689, 266.693548387, \
         266.445403663, 265.647577306, 263.561136482, 261.756436584, 259.352196396, 257.189278937, \
         254.699722764, 252.1727774, 249.40082503, 246.196034026, 243.204912023, 239.861473171, \
         236.216083507, 232.418024226, 228.347281141, 225.787503104, 225.012782593, 223.918548387, \
         225.243329973, 227.262907582, 228.317952383, 228.861412332, 228.843405442, 228.573161, \
         228.953055371, 229.14674716, 229.541472058, 229.232789713, 229.625510244, 229.774010109, \
         229.77615013, 229.971814334, 230.00521917]
    q = [4.43525008771, 4.43525008771, 3.91163212244, 4.11013212532, 4.09430589408, 4.20118246395, \
         3.96302160684, 3.79492112852, 3.72098837478, 3.46258558754, 3.40618549525, 3.08298591678, \
         2.59367328686, 2.44014813401, 2.05108776133, 2.26184155599, 1.86362280966, 1.47139101628, \
         1.21534239674, 1.14851623259, 1.21928331259, 1.01076543751, 0.823743251872, 0.614694927612, \
         0.54732936928, 0.409350331268, 0.263205149464, 0.188926862496, 0.149236855357, 0.126514591073, \
         0.101879556598, 0.0612401049152, 0.0545139393063, 0.0340927948804, 0.0221047169385, 0.0235166588019, \
         0.00356754891666, 0.00197738522031, 0.00107859019429, 0.00114148083765, 0.00113934522762, \
         0.00110772385208, 0.00115240565472, 0.00117580628454, 0.00122482417583, 0.00118633803273, \
         0.00123549487023, 0.0012545554029, 0.00125483200763, 0.00128035552547, 0.00128475944956]
    rundecker(z, p, t, q, tape5='tp5', co2_ppm=415.0, o3_ppm=1.0, co_ppm=1.0, ch4_ppm=1.0, n2o_ppm=1.0)
