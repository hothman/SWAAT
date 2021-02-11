#!/usr/bin/python3
__author__ = "Houcemeddine Othman"
__credits__ = "Wits University H3Africa/GSK ADME project"
__maintainer__ = "Houcemeddine Othman"
__email__ = "houcemoo@gmail.com"


import pandas as pd 
import re as re 
import datetime as dt
import glob  
import argparse

# Bokeh is optional 
try: 
    from bokeh.plotting import figure, save,  output_file
    from scipy.stats import norm
    import numpy as np
    from bokeh.models import  ColumnDataSource, HoverTool, Cross
    from bokeh.layouts import row
    from bokeh.embed import components
    bokeh = True
except: 
    bokeh = False


features = { 1:'signal peptide', 
2: 'propeptide',
3: 'transit peptide',
4: 'rsidue part of a calcium binding region',
5: 'rsidue part of a zinc binding regon', 
6: 'residue implicated in DNA binding', 
7: 'residue implicated in nucleotide phosphate binding, e.g. ATP',
8: 'residue associated to functional property',
9: 'residue belongs to coiled-coil region', 
10: 'residue belongs to sequence motif of biological interest',
11: 'residue belongs to the active site',
12: 'residue belongs to a metal binding site',
13: 'residue belongs to a binding site',
14: 'residue involved in a generic activity',
15: 'residue involved in a post-translational modification',
16: 'residue involved in covalent binding of a lipid moiety',
17: 'residue involved in the attachment with a glycan group',
18: 'residue involved in a disulfide bond',
19: 'residue involved with a crosslink bonding with another amino acid',
20: 'external annotation'}

ANNOTATIONHOME="/home/houcemeddine/BILIM/testing_SWAAT/myoutput/prot_annotation"
HOTSPOTSPATCHES="/home/houcemeddine/BILIM/testing_SWAAT/myoutput/hotspots"
UNIPROT2PDBHOME="/home/houcemeddine/BILIM/testing_SWAAT/myoutput/uniprot2PDBmap"
DATAHOME="/home/houcemeddine/BILIM/SWAAT/data/dGdS.csv"  # change this to relative path 
FASTAHOME="/home/houcemeddine/BILIM/testing_SWAAT/myoutput/sequences/Refseq"
FTMAPHOME="/home/houcemeddine/BILIM/testing_SWAAT/myoutput/ftmap"

template_header= """
  
<!DOCTYPE html>
<style>
body {{ font-family: sans-serif; 
}}
table, td, th {{
    background: #f5f5f5;
    border-collapse: collapse;
    border: 1px solid silver;
    box-shadow: inset 0 1px 0 #fff;
    font-size: 12px;
    line-height: 24px;
    margin: 50px auto;
    text-align: left;
    width: 800px;

}}
</style>

<html>
<head>
    <script src="https://cdn.bokeh.org/bokeh/release/bokeh-2.1.1.min.js"></script>  
    <title>{0}</title>
    <h1>{0}</h1> 
    <h2>{1}</h2> 
    <p> Please consider citing the following reference for SWAAT</p>
</head>

""".format("SWAAT report",  dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


class cleanData:
    """docstring for cleanData"""
    def __init__(self, swaat_predction, variant_file):
        self.swaat_predction = swaat_predction
        self.variant_file = variant_file
    
    def readFiles(self): 
        self.predicted_data = pd.read_csv(self.swaat_predction)
        self.var_data = pd.read_csv(self.variant_file, names=["gene_name","chromosome","position",
            "ref_allele","alt_allele","ref_AA","AA_position","mutant_AA"])

    def cleanDf(self):
        # characters are not clean,  need to keep only alpha numerical characters
        regex = re.compile('[^a-zA-Z0-9]')
        l1 = list(self.predicted_data["gene_name"]+self.predicted_data["wt_res"]+self.predicted_data["position"].astype(str)+self.predicted_data["mut_res"])
        l2 = list(self.var_data["gene_name"]+self.var_data["ref_AA"]+self.var_data["AA_position"].astype(str)+self.var_data["mutant_AA"])
        list1= [regex.sub('', element) for element in l1]
        list2= [regex.sub('', element) for element in l2]
        # generate the ID column for predicted_data
        self.predicted_data["var_id"]=list1
        # generate the ID column for var_data
        self.var_data["var_id"] = list2
        self.merged_data = pd.merge(self.var_data, self.predicted_data,on='var_id',how='left')
        return self.merged_data

    def addAnnotation(self):
        annotation, ftmap_annotation = getAnnotationList(self.merged_data)
        coverage, offset = isCovered(self.merged_data)
        hotspotpatch = getHotSpotPatch(self.merged_data, offset = offset)
        self.merged_data["annotation"] = annotation
        self.merged_data["hotspotpatch"] = hotspotpatch
        self.merged_data["Covered by the structure"] = coverage
        myannotation=transformAnnotation(self.merged_data)
        self.merged_data["Putative drug interaction"] = ftmap_annotation
        self.merged_data["red flags"] = myannotation
        #self.merged_data["putative_drug_binder"] = 

def _getUNIPROT(genename):
    """
    fetches the uniprot code from fasta file generated by 
    the auxiliary workflow
    """
    file = FASTAHOME+"/"+genename+"_refseq.fa"
    with open(file, "r") as mapfile: 
        uniprot_code= (mapfile.readlines()[0].split("|")[0].replace(">", ""))
    return uniprot_code

def _getAnnotation(genename, res_position):
    """
    It generates a formatted text for each amino acid
    in the annotated list
    """
    uniprot = _getUNIPROT(genename)
    
    annotation_data = pd.read_csv(ANNOTATIONHOME+"/"+uniprot+"_annotation.csv")
    
    annot_list=[]
    if (annotation_data["residue_id"]== res_position ).any() : 
        subdf = annotation_data[annotation_data["residue_id"]== res_position]

        for i in range(0, len(subdf)):
            Annot_id = int(subdf.iloc[i].annotation_tag)
            try : # this will insure to add only the tag if the content of the note column in the annotation file is empty
                all_annotation= str(features[Annot_id])+": "+subdf.iloc[i].note
            except: 
                all_annotation= str(features[Annot_id])
            annot_list.append(all_annotation)    
    return '\n'.join(annot_list)


def isCovered(combineddataframe):
    inrange_list=[]
    for i in range(0,len(combineddataframe)):
        gene = combineddataframe.iloc[i]['gene_name_x']
        
        mapfile = UNIPROT2PDBHOME+"/"+gene+"*.tsv"
        path_to_map=glob.glob(mapfile)[0]    # because pandas does not read wildcards 
        map_data=pd.read_csv(path_to_map, sep="\t", comment='#')
        list_of_covered_residues = list(map_data.IDref)
        if combineddataframe.iloc[i]['AA_position'] in list_of_covered_residues: 
            inrange_list.append(True)
        else: 
            inrange_list.append(False)
    offset = map_data.IDref[0]- map_data.ID[0]
    return inrange_list, offset

def _getFtmap(genename, res_position): 
    """
    This function parses the outputs of ftmap 
    """
    uniprot = _getUNIPROT(genename)
    annotation_data = pd.read_csv(FTMAPHOME+"/"+genename+"_nb.csv")
    annot_list=[]
    if (annotation_data["resID"] == res_position ).any() :
        subdf = annotation_data[annotation_data["resID"]== res_position]
        probe_contacts = int(subdf["probe_contacts"])
        zscore = float(subdf["zscore"])
        percentile = float(subdf["percentile_score"])
        annot_list.append("Number of probe molecules = {0}, Zscore = {1}, Percentile Score = {2}%".format(probe_contacts, zscore, percentile ))
    return '\n'.join(annot_list)
   
def getAnnotationList(combineddataframe):
    """
    walks through the variants in the list 
    to generate the annotation of aa
    returns a list of strs 
    """
    annotation_list=[]
    ftmap_annotation = []
    for i in range(0,len(combineddataframe)) : 
        gene = combineddataframe.iloc[i]['gene_name_x'] 
        aa_position = combineddataframe.iloc[i]['AA_position']
        annotation_list.append(_getAnnotation(gene, aa_position) )
        ftmap_annotation.append(_getFtmap(gene, aa_position) )

    return annotation_list, ftmap_annotation

def getHotSpotPatch(combineddataframe, offset=0): 
    hs_list = []
    for index in range(0,len(combineddataframe)) :
        gene = combineddataframe.iloc[index]['gene_name_x']
        uniprot= _getUNIPROT(gene) 
        hotspot_file = HOTSPOTSPATCHES+"/"+uniprot+"_residueClus.csv"
        hs_data = pd.read_csv(hotspot_file)
        corrected_to_real_positions = list(hs_data["res_ID"]+offset )  # offset to convert to real coordinates in the hotspost patch files
        aa_position = combineddataframe.iloc[index]['AA_position'] 
        if aa_position in corrected_to_real_positions: 
            hs_list.append("1")
        else: 
            hs_list.append("0") 
    return hs_list

def transformAnnotation(combineddataframe):
    salt_bridge_diff = combineddataframe["sb_wt"]-combineddataframe["sb_mut"] 
    combineddataframe.insert(47, 'Salt bridge breakage', salt_bridge_diff)  
    annotations = combineddataframe.iloc[:,list(range(35, 48))+[50]]  # indexes indicaes the 3Dmissense features and the hotspotpatch   
    annotation_table_for_raws = []
    for i in range(0,len(annotations)):
        raw = annotations.iloc[i].to_dict()
        threed_missense = []
        for key in raw: 
            if raw[key] in [1, "1"] : 
                threed_missense.append("<span style='color: #c71100;'>&#9873;</span>"+key)  # this is where you can put the red flag
            elif raw[key] in [-1, "-1"] :
                threed_missense.append("<span style='color: #c71100;'>&#9873;</span>Salt bridge formation")  # -1 value is only in Salt bridge breakage column 
        if not threed_missense: 
            annotation_table_for_raws.append("")
        else: 
            annotators = ','.join(threed_missense)
            annotation_table_for_raws.append(annotators)
    return annotation_table_for_raws

######################## 
## A class for interactive plotting
#######################

def readData(): 
    return pd.read_csv(DATAHOME)

def getPdf(data):
    x=data.sort_values()
    #remove outlayers
    mu, std = norm.fit(x)
    final_list = [energy for energy in x if (energy > mu - 3 * std)]
    final_list = [energy for energy in final_list if (energy < mu + 3 * std)]
    mu, std = norm.fit(final_list)
    return np.array(final_list), norm.pdf(final_list, mu, std)
    
    
class Plot: 
    """
    This is a warapping class to plot interactive graphs 
    for SWAAT report.
    takes as argument the dataframe of the processed variants
    """
    def __init__(self, processsed_dataframe):
        self.processsed_dataframe = processsed_dataframe
        self.data = readData()
        self.dist_dG = getPdf(self.data['dG'])
        self.y_coor_dg = self.dist_dG[1].max()      
        self.dist_dS = getPdf(self.data['dS'])
        self.y_coor_ds = self.dist_dS[1].max()
        ## Wrapping 
        Plot.prepareDataSource(self)
        Plot.renderGraph(self)
        Plot.embedingCode(self)
         
    def prepareDataSource(self): 
        y_coors_dg = np.array( len( self.processsed_dataframe )*[self.y_coor_dg])-0.2
        y_coors_ds = np.array( len( self.processsed_dataframe )*[self.y_coor_ds])-0.7
        aa_variants = self.processsed_dataframe["Reference residue"] + self.processsed_dataframe["Residue position"].astype(str)+ self.processsed_dataframe["Residue variant"] 
        self.source = ColumnDataSource(data = {"dG": self.processsed_dataframe["dG (kcal/mol)"],
                                "dS": self.processsed_dataframe["dS (kcal/mol)"],
                                "Chromosome": self.processsed_dataframe.Chromosome, 
                                "Genome_position": self.processsed_dataframe.position, 
                                "y_coors_dg":y_coors_dg ,
                                "y_coors_ds":y_coors_ds ,
                                "aa_variants": aa_variants } )
    
    def PlotGraph(self, title, x_axis_label, y_axis_label, x, y, energy_tag, energy_column): 
        p = figure(title=title, x_axis_label=x_axis_label, 
                   y_axis_label=y_axis_label,
                  plot_width=400, plot_height=400)
        g1 = Cross(x=x, y=y)
        p.add_glyph(source_or_glyph=self.source, glyph=g1)
        g1_r = p.add_glyph(source_or_glyph=self.source, glyph=g1)
        g1_hover = HoverTool(renderers=[g1_r],
                         tooltips=[(energy_tag, energy_column),
                   ("Chromosome", "@Chromosome"),
                   ("Position", "@Genome_position"),
                   ("Variant", "@aa_variants")])
        p.add_tools(g1_hover)
        if energy_tag == "dG": 
            fit_data = self.dist_dG
            legend_label = "PDF(dG)"
        elif energy_tag == "dS": 
            fit_data = self.dist_dS
            legend_label = "PDF(dS)"        
        p.line(fit_data[0], fit_data[1], legend_label=legend_label, line_width=2, color="green")
        p.circle(x=x, y= y, source=self.source, size=5, color="red")
        p.vbar(x=x, bottom=0, top=y, source=self.source, width=0.01, color="red" )
        return p
    
    def renderGraph(self): 
        plot1= self.PlotGraph(title="", x_axis_label='dG(kcal/mol)', 
                          y_axis_label="Probability", x="dG", y='y_coors_dg', 
                          energy_tag="dG", energy_column="@dG")
        plot2= self.PlotGraph(title="", x_axis_label='dS(kcal/mol)', 
                          y_axis_label="Probability", x="dS", y='y_coors_ds', 
                          energy_tag="dS", energy_column="@dS")
        self.plot = row(plot1, plot2)
        #output_file("test.html")
        #save(self.plot)
    
    def embedingCode(self): 
        self.script, self.div = components(self.plot)
        

    
##################
# End of Plot class
##################


class formatHtML: 
    def __init__(self, dataframe):
        """self.dataframe is inhirited from the cleanData class"""
        self.dataframe = dataframe.sort_values(by=["position_x"])

    def list_of_genes(self) : 
        self.genes =  list( set(self.dataframe["gene_name_x"]) )

    def _generalStatistics(self, gene):
        non_processed_aa= pd.DataFrame(columns = self.dataframe.columns)
        processed_aa = pd.DataFrame(columns = self.dataframe.columns)
        gene_dataframe = self.dataframe[self.dataframe["gene_name_x"] == gene]
        total_number_of_coding_variants = len(gene_dataframe)
        number_of_indels = 0
        for i in range(0, len(gene_dataframe)) : 
            if "_" in  gene_dataframe.iloc[i].mutant_AA : 
                number_of_indels =+ 1

            if  pd.isnull(gene_dataframe.iloc[i].swaat_prediction) :
                non_processed_aa = non_processed_aa.append(gene_dataframe.iloc[i]) 
            else: 
                #print(gene_dataframe.iloc[i])
                processed_aa = processed_aa.append(gene_dataframe.iloc[i], ignore_index=True) 
        return total_number_of_coding_variants, number_of_indels, len(non_processed_aa), len(processed_aa), non_processed_aa, processed_aa

    def _cleanHtmlDf(self, dataframe, mode):
        print(dataframe)
        to_format_columns = ["gene_name_x", "gene_name_y", "wt_res","var_id", "mut_res", "position_y", 
        "subScore", "grantham", "sneath", "classWT",  "classMut", "sasa_mut", "sasa_wt", "hyrophob_WT",
         "hyrophob_Mut",  "volume_WT", "volume_Mut", "pssm_mut", "pssm_wt", "Covered by the structure", 
         'disulfide_breakage', 'buried_Pro_introduced', 'buried_glycine_replaced', 'buried_hydrophilic_introduced',  
         'buried_charge_introduced', 'buried_charge_switch', 'sec_struct_change', 'buried_charge_replaced',
           'buried_exposed_switch', 'exposed_hydrophilic_introduced', 'Buried_salt_bridge_breakage', "Large_helical_penality_in_alpha_helix", 
           'hotspotpatch', 'Salt bridge breakage']

        if mode == "processed":
            is_annotation_empty = list(dataframe["annotation"] == "")
            clean_dataframe = dataframe.drop(columns=to_format_columns)
            clean_dataframe.columns =["Chromosome", "position", "Reference Allele", "Aternative allele", "Reference residue", 
            "Residue position", "Residue variant",  "Chain", "dG (kcal/mol)", "Secondary structure", "dS (kcal/mol)", 
            "#hydrogen bonds ref", "#hydrogen bonds var", "#salt bridges var", "#salt bridges ref", "SASA ratio", "ML prediction", 
            "Annotation", "Putative drug binder", "Red Flags"  ]
            return clean_dataframe

        elif mode == "non-covered": 
            to_format_columns = to_format_columns+["chain", "dG","SecStruc","dS","hb_mut","hb_wt","sb_mut","sb_wt","sasa_ratio","swaat_prediction", "red flags"]
            col_names = ["Chromosome", "Position", "Reference Allele", "Aternative allele", "Reference residue", "Residue position", "Residue variant" ,"Putative drug binder", 
            "Annotation"] 
            is_annotation_empty = list(dataframe["annotation"] == "")  # erturns a list of booleans
            if all(is_annotation_empty) : 
                to_format_columns=to_format_columns + ["annotation"]
                col_names.remove('annotation')
            is_drug_empty = list(dataframe["Putative drug binder"] == "")
            if all(is_drug_empty) :
                to_format_columns=to_format_columns+ ["Putative drug binder"]
                col_names.remove('Putative drug binder')             
            clean_dataframe = dataframe.drop(columns=to_format_columns)
            clean_dataframe.columns = col_names
            return clean_dataframe


    def _cleanIndels(self, dataframe): 
        newdf = (dataframe[["gene_name_x" ,"chromosome" ,"position_x" ,"ref_allele" ,"alt_allele", "ref_AA" ,"AA_position", "annotation"]])
        newdf.columns= ["Gene name", "Chromosome", "Position", "Reference Allele", "Aternative allele", "Reference residue", "Residue position", "annotation" ]
        return newdf

    def _cleanNonProcessed(self, dataframe):
        newdf = (dataframe[["gene_name_x" ,"chromosome" ,"position_x" ,"ref_allele" ,"alt_allele", "annotation"]])
        newdf.columns= ["Gene name", "Chromosome", "Position", "Reference Allele", "Aternative allele", "annotation" ]
        return newdf


    def outputHtml(self): 
        """
        A arapping mehod to generate HTML report
        """
        with open("swaat2.html", "w") as file:
            file.write(template_header) 
            file.write("<ul>")
            for gene in self.genes : 
                file.write("<li><a href='#{0}'>{0}</a> </li>".format(gene))
            file.write("</ul>")

            for gene in self.genes: 
                # reporting counts of variants by category 
                n_variants, n_indels, n_non_processed, n_processed, df_nonprecessed, df_processed = self._generalStatistics(gene)
                file.write("<h3 id='{0}'>{0}</h3>".format(gene))
                file.write("<p style='margin-left: 40px' >Total number of coding variants: {0}  \
                    <ul> <li> Processed variants: {1} </li> \
                         <li> Non processed variants: {2} </li> \
                         <li> indels: {3} </li> </ul> </p>".format(n_variants, n_processed, n_non_processed, n_indels))
                # subset variants only for 'gene'
                vars_for_gene = self.dataframe[self.dataframe["gene_name_x"] == gene]

                # report non covered 
                file.write("<h4>Non processed variants summary</h4>")
                if vars_for_gene[vars_for_gene["Covered by the structure"] == False].empty : 
                    file.write("<p>All variants are covered by the PDB structure</p>")
                else: 
                    df_non_processed = self._cleanHtmlDf( vars_for_gene[vars_for_gene["Covered by the structure"] == False], mode="non-covered" )
                    html_non_processed = df_non_processed.to_html(index=False, escape=False)
                    file.write(html_non_processed)

                # report indels         
                indel_table = vars_for_gene[vars_for_gene["mutant_AA"] == "_"]
                if not indel_table.empty:
                    file.write("<h4>Indels summary</h4>")
                    indel_table = self._cleanIndels(indel_table)
                    file.write(indel_table.to_html(index=False, escape=False))

                # reporting details for processed variants
                self.df_processed = self._cleanHtmlDf(df_processed, mode="processed")
                self.df_processed.sort_values(by=["Residue position"],  inplace=True)
                #print(self.df_processed.columns)
                file.write("<h4>Processed variants summary</h4>")

                html_processed = self.df_processed.to_html(index=False, escape=False)
                file.write(html_processed)
                file.write("<hr>")

                # integrate interactive plot (needs bokeh library)
                try: 
                    if bokeh :
                        file.write("Hover to explore")
                    interactive_plot = Plot(self.df_processed)
                    file.write(interactive_plot.div)
                    file.write(interactive_plot.script)
                except: 
                    print("Install Bokeh to explore the result interactively")



parser = argparse.ArgumentParser(description=" This script formats the output of SWAAT to interactive html")

        
if __name__ == "__main__":
    #instance = cleanData("/home/houcemeddine/BILIM/SWAAT/main/work/f6/dc7d37010ef1cd778088353e8fd993//predicted_outcomes.csv",
    #"/home/houcemeddine/BILIM/SWAAT/main/work/06/5ea0bd667b002ae218083b19714295/allVariantsInOneFile.csv" )

    # add long and short argument
    parser.add_argument("--prediction", help="A file containing the calculated features and the prediction by the ML model")
    parser.add_argument("--variants", help="A list of variants")
    args = parser.parse_args()
    instance = cleanData(args.prediction, args.variants )
    """ exapmle 
    python formatOutput --prediction predicted_outcomes.csv --variants allVariantsInOneFile.csv
    """
    instance.readFiles()
    instance.cleanDf()
    instance.addAnnotation()
    html = formatHtML( instance.merged_data)
    html.list_of_genes()
    html.outputHtml()
