import pandas as pd
import argparse
import os
import re

parser = argparse.ArgumentParser(description="Example script for parsing command line arguments.")
parser.add_argument('-p', '--path', type=str, help="The path to the directory", required=True)
#parser.add_argument('-o', '--output', type=str, help="The output file name", required=True)
parser.add_argument('-fo', '--focus', type=str, help="The target species", required=True)

args = parser.parse_args()
path = args.path
#output = args.output
focus = args.focus

last_rel = os.path.join(path, "last_rel.txt")
last_rel = pd.read_csv(last_rel, sep="\t")
max_lnL_index = last_rel['lnL'].idxmax()
model_folder = last_rel.loc[max_lnL_index, 'Model']
#print(f"The folder corresponding to the maximum lnL value is: {model_folder}")

best_model = os.path.join(path, model_folder)
best_model_res = os.path.join(best_model, "REL.txt")

f_content = open(best_model_res).read().split("w ratios as labels for TreeView:")
tree = f_content[1]

cleaned_tree = re.sub(r'^\s*$\n|^Time used.*\n', '', tree, flags=re.MULTILINE)
pattern = re.escape(focus) + r"\s+'#([0-9.-]+)'"
focus_w = re.findall(pattern, cleaned_tree)
float_focus_w = float(focus_w[0])

#print(cleaned_tree, focus_w, float_focus_w)

if float_focus_w > 1:
    print(f"{best_model_res}\t{float_focus_w}\tPositive")
else:
    print(f"{best_model_res}\t{float_focus_w}\tNegative")