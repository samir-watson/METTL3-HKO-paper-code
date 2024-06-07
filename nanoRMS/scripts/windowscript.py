import pandas as pd
import sys
import dask.dataframe as dd
import os
# Scripts how to use
# python nanopolish_window.py positions_file input_table label

args = sys.argv[1:]

contig = args[0]
position = args[1]
label1 = args[3]

# read the input file, prefer parquet format
if not os.path.exists(args[2]+'.parquet'):
	print("parquet not found.	Reading csv and converting to parquet")
	input1 = pd.read_csv(args[2],sep='\t')
	print("Converting to parquet")
	input1.to_parquet(args[2]+".parquet")
	print("Done converting to parquet")
else:
	print("Reading parquet")
	input1 = pd.read_parquet(args[2]+".parquet")
	print("Done reading parquet")




###
#This is probably not a great idea, but it works for now
###
pd.set_option('mode.chained_assignment', None) # To avoid the SettingWithCopyWarning
###

def process(contig,position,data, label):
	print("Processing {}_{}".format(contig,position))
	chr = str(contig)
	subs = data[data['contig'] == chr]
	subs['position'] += 3 # Add 3 nt to each position 
	subs['Pos'] = subs['contig'].astype(str) + '_' + subs['position'].astype(str) # Unique column
	subs['sample'] = label #Add label
	windows = pd.DataFrame() # Create an empty DataFrame
	#Create windows file
	
	mod = int(position)
	window = [mod + i for i in range(-7, 8)]
	for wind in window:
		subs2 = subs[subs['position'] == wind]
		if len(subs2) > 1:
			subs2['modification'] = chr + '_' + str(mod)
			subs2['reference'] = str(wind - mod)
			windows = pd.concat([windows, subs2])
	return windows

final_window = process(contig,position, input1, label1)
os.makedirs("steps/windows/{label}_windows".format(label=label1),exist_ok=True)
final_window.to_csv("steps/windows/{label}_windows/{label}_{contig}_{position}.window_file.tsv.tmp".format(label=label1,contig=contig,position=position), sep='\t', index=False)