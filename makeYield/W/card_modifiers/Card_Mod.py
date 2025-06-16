import argparse
import math

parser = argparse.ArgumentParser(description='For scaling events by lumi')

# No longer necessary
parser.add_argument('--luminosity', action='store', default=97.7, help="Lumi scale\n DEFAULT: 97.7")
parser.add_argument('--categ', action='store', default='catA', help="Card Category\n DEFAULT: catA")
parser.add_argument('--sig_exp', action='store', default=1.0, help="Signal yield\n DEFAULT: 1.0")
parser.add_argument('--bkg_exp', action='store', default=1.0, help="bkg yield\n DEFAULT: 1.0")
parser.add_argument('--sb_exp', action='store', default=1.0, help="sideband yield\n DEFAULT: 1.0")
parser.add_argument('--ext_unc', action='store', default=1.10, help="extrapolation factor uncertainty\n DEFAULT: 1.10")

args = parser.parse_args()
categ = args.categ
sig = args.sig_exp
bkg = args.bkg_exp
sb = args.sb_exp
ext_uncert = args.ext_unc
sb_rounded_up = math.ceil(float(sb))
round_uncert = 1.0 + abs(sb_rounded_up-float(sb))/sb_rounded_up

card_name = ''
linenum_actual_1 = 10
linenum_actual_2 = 10
linenum_actual_3 = 10
linenum_actual_4 = 10

if(categ == 'catA'):
        card_name = 'W_T3mu_catA_bdtcut'
        linenum_actual_1 = 14
        linenum_actual_2 = 19
        linenum_actual_3 = 20
        linenum_actual_4 = 21
if(categ == 'catB'):
        card_name = 'W_T3mu_catB_bdtcut'
        linenum_actual_1 = 14
        linenum_actual_2 = 19
        linenum_actual_3 = 20
        linenum_actual_4 = 21
if(categ == 'catC'):
        card_name = 'W_T3mu_catC_bdtcut'
        linenum_actual_1 = 14
        linenum_actual_2 = 19
        linenum_actual_3 = 20
        linenum_actual_4 = 21

# opening the file in read mode
file = open("card_modifiers/"+card_name+".txt", "r")
replacement = ""
# using the for loop
for linenum, line in enumerate(file):
    #print(linenum)
    if(linenum==linenum_actual_1):
        #print(line)
        line_mod=line
        x1 = line.split()
        res1 = "{:.4f}".format(float(sig))
        res2 = "{:.4f}".format(float(bkg))
        line_mod=line_mod.replace(x1[1], res1)
        line_mod=line_mod.replace(x1[2], res2)
                
        replacement = replacement + line_mod
    
    elif(linenum==linenum_actual_2):
        #print(line)
        line_mod=line
        x1 = line.split()
        res1 = "{:.0f}".format(float(sb_rounded_up))
        line_mod=line_mod.replace(' '+x1[2]+' ', ' '+res1+' ')
                
        replacement = replacement + line_mod
        
    elif(linenum==linenum_actual_3):
        #print(line)
        line_mod=line
        x1 = line.split()
        res1 = "{:.3f}".format(float(ext_uncert))
        line_mod=line_mod.replace(' '+x1[3]+' ', ' '+res1+' ')
                
        replacement = replacement + line_mod
        
    elif(linenum==linenum_actual_4):
        #print(line)
        line_mod=line
        x1 = line.split()
        res1 = "{:.6f}".format(round_uncert)
        line_mod=line_mod.replace(' '+x1[3]+' ', ' '+res1+' ')
                
        replacement = replacement + line_mod
    
    else:
        replacement = replacement + line
            
file.close()

# opening the file in write mode
fout = open("modified_dc_"+categ+".txt", "w")
fout.write(replacement)
fout.close()
