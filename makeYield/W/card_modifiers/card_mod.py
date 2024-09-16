import argparse

parser = argparse.ArgumentParser(description='For scaling events by lumi')

parser.add_argument('--luminosity', action='store', default=97.7, help="Lumi scale\n DEFAULT: 97.7")

args = parser.parse_args()

# opening the file in read mode
file = open("card_modifiers/W_combined.txt", "r")
replacement = ""
# using the for loop
for linenum, line in enumerate(file):
    #print(linenum)
    if(linenum==30):
        #print(line)
        line_mod=line
        x1 = line.split()
        x1 = x1[1:]
        #print(x1)
        for i in x1:
            imod=float(i)*(float(args.luminosity)/97.7)
            #print(imod)
            res = "{:.4f}".format(imod)
            #print(res)
            line_mod=line_mod.replace(i, res)
        print(line_mod)
        replacement = replacement + line_mod
    else:
        replacement = replacement + line
            
file.close()

# opening the file in write mode
fout = open("Combined_Mod.txt", "w")
fout.write(replacement)
fout.close()
