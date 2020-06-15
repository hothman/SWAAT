#!/usr/bin/python3
__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"


import pandas as pd
import argparse
import pickle


def cleanData(swaat_data_csv): 
	data = pd.read_csv(swaat_data_csv)
	datacopy = data.copy()
	data.drop(["gene_name"], axis=1, inplace=True)
	# calculate delta hydrophobicity
	data['hydrophobicity'] = data.hyrophob_WT-data.hyrophob_Mut
	# transform sasa
	dSASA =  data.sasa_mut - data.sasa_wt
	data['dSASA'] = dSASA
	sasa_r = []
	for sasa in data.sasa_ratio: 
	    sasa_r.append( float(sasa) )
	data['sasa_r'] = sasa_r

	dataset=data[["dG", "RMSIP", "dS", "subScore", "grantham", "sneath", "sb_mut",
	 "sb_wt", "volume_WT", "volume_Mut", "pssm_mut", "pssm_wt",  "dSASA",	
	 "sasa_r", "hydrophobicity"]]
	return dataset, datacopy

def predictVariant(test_data_frame): 
	loaded_model = pickle.load(open(MODELHOME, 'rb'))
	#result = loaded_model.score(test_data_frame)
	# 0: no effect, has an effect: with effect deleterious or stabilizing 
	predictions = loaded_model.predict(test_data_frame)
	return predictions

def outputPredictions(dataframe, predictions, output="predictions_swaat.csv"):
	dataframe["swaat_prediction"]=predictions
	dataframe.to_csv(output, index=False)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=" Usage: \
	deployML.py --modelpickle swaat_rf.ML --inputdata swaat_vars.csv --output myoutput.csv")
	# add long and short argument
	parser.add_argument("--inputdata", help="Input data in csv format with header")
	parser.add_argument("--modelpickle", help="Input of the serialized model in pickle format")
	parser.add_argument("--output", help="outputname")
	args = parser.parse_args()
	MODELHOME=args.modelpickle
	try:
		outpufile = args.output
	except: 
		outputfile= "predictions_swaat.csv"


	#MODELHOME="/home/houcemeddine/BILIM/ADME_PGx/SnpsInPdb/MLmodel/swaat_rf.ML"
	#file="/home/houcemeddine/BILIM/SWAAT/main/swaatall"
	cleandf, rawdf = cleanData(args.inputdata)
	predictions=predictVariant(cleandf)
	outputPredictions(rawdf,predictions )



